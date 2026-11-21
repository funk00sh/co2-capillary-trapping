"""Pore-geometry sensitivity of the Sgr-vs-angle trend (numerical experiment on networks A and B).

Knobs (one at a time): coord (random throat deletion), radius (throat radius scale), shape (throat shape-factor scale).
Per network variant: (1) calibrate the primary-drainage Pc that gives Sgi = TARGET, (2) run the four angle ranges x 3 seeds
at that Pc. Resumable: finished runs are appended to results.jsonl.
Usage: python run_geometry_sweep.py <work_dir>
"""
import json
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

import run_matrix as rm
from make_variant_network import build
from parse_filling_events import parse

SRC = {"A": rm.HERE / "network_pnextract", "B": rm.HERE / "network_B"}
KNOBS: Dict[str, List[float]] = {"coord": [0.1, 0.2, 0.3, 0.4], "radius": [0.7, 0.85, 1.2, 1.4], "shape": [0.5, 0.75, 1.5, 2.0]}
SEEDS = [11, 37, 50]
# second round: push the knobs towards where the first round suggested a crossing
EXTRA: Dict[str, Dict[str, List[float]]] = {"A": {"radius": [0.55, 0.45, 0.35], "shape": [0.25, 0.12]},
                                            "B": {"radius": [1.7, 2.0, 2.4], "shape": [2.5, 3.0]}}
TARGET = 0.60
PC_CAL_KPA = 40.0


def variants() -> List[dict]:
    """All network variants (baseline included once per network)."""
    out = []
    for n in "AB":
        out.append(dict(vid=f"{n}_base", net=n, knob="base", level=0.0, real=0))
        for knob, levels in KNOBS.items():
            for lv in levels + EXTRA.get(n, {}).get(knob, []):
                for r in (SEEDS if knob == "coord" else [0]):
                    out.append(dict(vid=f"{n}_{knob}_{lv:g}" + (f"_r{r}" if knob == "coord" else ""), net=n, knob=knob, level=lv, real=r))
    return out


def deck(pc_kpa: float, wett: str, seed: int) -> str:
    """PNflow input deck (identical to the production matrix)."""
    return f"""TITLE: run;
NETWORK: F Berea;
cycle1: 0.0 {pc_kpa * 1000:.1f} 0.005 T F;
cycle2: 1.0 -200000 0.005 T F;
cycle1_BC: T F T T DP 1.0 2.0;
cycle2_BC: T F T T DP 2.0 1.0;
CALC_BOX: 0.1 0.9;
INIT_CONT_ANG: {rm.RECEDING};
EQUIL_CON_ANG: {rm.WETTABILITY[wett]};
Water: {rm.WATER};
Oil: {rm.CO2};
WaterOilInterface: {rm.BASE_IFT};
DRAIN_SINGLETS: T;
RAND_SEED: {seed};
RES_FORMAT: Excel;
"""


def run_pnflow(netdir: Path, rundir: Path, pc_kpa: float, wett: str, seed: int) -> Optional[Tuple[list, list, str]]:
    """Run PNflow in rundir with the network copied in; return drainage table, imbibition table, log."""
    rundir.mkdir(parents=True, exist_ok=True)
    for f in netdir.glob("Berea_*.dat"):
        shutil.copy(f, rundir / f.name)
    (rundir / "input.dat").write_text(deck(pc_kpa, wett, seed))
    try:
        subprocess.run([str(rm.EXE), "input.dat"], cwd=rundir, capture_output=True, timeout=2400)
    except subprocess.TimeoutExpired:
        return None
    for f in rundir.glob("Berea_*.dat"):
        f.unlink()
    try:
        return (rm.read_curve(rundir / "run_cycle1_drain.csv"), rm.read_curve(rundir / "run_cycle2_imb.csv"),
                (rundir / "run_pnflow.prt").read_text(errors="replace"))
    except FileNotFoundError:
        return None


def calibrate(v: dict, work: Path) -> dict:
    """Build the network, run a drainage-only calibration and find Pc for Sgi = TARGET."""
    nd = work / "nets" / v["vid"]
    stats = build(SRC[v["net"]], nd, v["knob"] if v["knob"] != "base" else "coord", v["level"], v["real"])
    seed = v["real"] if v["knob"] == "coord" and v["real"] else 11
    res = run_pnflow(nd, work / "cal" / v["vid"], PC_CAL_KPA, "WW_20-40", seed)
    dr = np.array(res[0])[:, :2]
    sgi, pc = 1 - dr[:, 0], dr[:, 1] / 1000
    target = min(TARGET, float(sgi.max()) - 0.03)
    i = int(np.argmax(sgi >= target))
    pc_t = float(np.interp(target, sgi[max(i - 1, 0):i + 1], pc[max(i - 1, 0):i + 1])) if i > 0 else float(pc[0])
    return dict(v, pc_kpa=round(pc_t, 3), sgi_target=round(target, 3), sgi_max=float(sgi.max()), stats=stats)


def run_one(job: Tuple[dict, str, int, Path]) -> dict:
    """One production run; the result row."""
    v, wett, seed, work = job
    seed_eff = v["real"] if v["knob"] == "coord" and v["real"] else seed
    rd = work / "runs" / f"{v['vid']}_{wett}_s{seed_eff}"
    res = run_pnflow(work / "nets" / v["vid"], rd, v["pc_kpa"], wett, seed_eff)
    if res is None:
        return {}
    dr, im, prt = res
    row = dict(vid=v["vid"], net=v["net"], knob=v["knob"], level=v["level"], real=v["real"], wett=wett, seed=seed_eff,
               pc_kpa=v["pc_kpa"], Sgi=1 - dr[-1][0], Sgr=1 - im[-1][0], **parse(prt))
    for f in ("run_res",):
        shutil.rmtree(rd / f, ignore_errors=True)
    return row


def main() -> None:
    """Calibrate, run, append results."""
    work = Path(sys.argv[1])
    work.mkdir(parents=True, exist_ok=True)
    cal_file = work / "calibration.json"
    vs = variants()
    cal = json.loads(cal_file.read_text()) if cal_file.exists() else []
    done = {c["vid"] for c in cal}
    todo = [v for v in vs if v["vid"] not in done]
    print(f"{len(vs)} variants, {len(todo)} to calibrate", flush=True)
    with ThreadPoolExecutor(8) as pool:
        for c in pool.map(lambda v: calibrate(v, work), todo):
            cal.append(c)
            print("cal", c["vid"], c["pc_kpa"], round(c["sgi_max"], 3), flush=True)
    cal_file.write_text(json.dumps(cal, indent=1))
    res_file = work / "results.jsonl"
    seen = set()
    if res_file.exists():
        for l in res_file.read_text().splitlines():
            r = json.loads(l)
            seen.add((r["vid"], r["wett"], r["seed"]))
    jobs = []
    for c in cal:
        seeds = [c["real"]] if c["knob"] == "coord" and c["real"] else SEEDS
        for w in rm.WETTABILITY:
            for s in seeds:
                if (c["vid"], w, s) not in seen:
                    jobs.append((c, w, s, work))
    print(len(jobs), "runs to do", flush=True)
    with ThreadPoolExecutor(8) as pool, open(res_file, "a") as fh:
        for r in pool.map(run_one, jobs):
            if r:
                fh.write(json.dumps(r) + "\n")
                fh.flush()
    print("finished", flush=True)


if __name__ == "__main__":
    main()
