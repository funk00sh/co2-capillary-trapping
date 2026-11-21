"""Pc-controlled drainage/imbibition matrix on the Berea network re-extracted with pnextract.

Unmodified upstream PNflow. Drainage stops at a prescribed maximum capillary
pressure (Sgi is an output); water re-invasion then runs to Pc = -200 kPa.
Usage: python scripts/run_matrix.py <work_dir>
"""
import os
import csv
import logging
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
EXE = Path(os.environ.get("PNFLOW_EXE", "pnflow.exe" if os.name == "nt" else "pnflow"))
NET = Path(os.environ.get("PNFLOW_NETWORK_A", ROOT / "data/network_A")).resolve()

# 50 C, 10 MPa reference properties (retained CoolProp fluid card)
WATER = "0.00054853798402524831 1.2 992.30800190305467"
CO2 = "2.8368236573179568e-05 1000 384.3271516766676"
DELTA_RHO = 992.308 - 384.327
BASE_IFT = 0.03737
PC_LADDER_KPA = [4.4, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.8,
                 6.0, 6.5, 7.0, 8.0, 10.0, 15.0, 20.0, 30.0]
RECEDING = "1 0 10 -0.2 -3 rand 0"  # primary drainage, strongly water-wet
WETTABILITY: Dict[str, str] = {  # advancing angles through water (model 4)
    "WW_20-40": "4 20 40 -0.2 -3 rand 0",
    "MW_45-65": "4 45 65 -0.2 -3 rand 0",
    "IW_70-90": "4 70 90 -0.2 -3 rand 0",
    "NW_95-115": "4 95 115 -0.2 -3 rand 0",
}
SEEDS = [11, 23, 37, 50, 71]
IFT_CASES = [0.025, 0.050]
MAX_WORKERS = 8


@dataclass(frozen=True)
class Case:
    """One PNflow run."""

    name: str
    wett: str
    ift: float
    pc_max_kpa: float
    seed: int


def build_cases() -> List[Case]:
    """Wettability x Pc x seed block plus an IFT block at base wettability."""
    cases = [
        Case(f"{w}_ift37_pc{pc:g}_s{s}", w, BASE_IFT, pc, s)
        for w in WETTABILITY for pc in PC_LADDER_KPA for s in SEEDS
    ]
    cases += [
        Case(f"MW_45-65_ift{int(i * 1000)}_pc{pc:g}_s50", "MW_45-65", i, pc, 50)
        for i in IFT_CASES for pc in PC_LADDER_KPA
    ]
    return cases


def write_input(case: Case, run_dir: Path) -> None:
    """Write the PNflow deck and copy the network."""
    run_dir.mkdir(parents=True, exist_ok=True)
    for f in list(NET.glob("Berea_node*.dat")) + list(NET.glob("Berea_link*.dat")):
        shutil.copy(f, run_dir / f.name)
    deck = f"""TITLE: run;
NETWORK: F Berea;
cycle1: 0.0 {case.pc_max_kpa * 1000:.1f} 0.005 T F;
cycle2: 1.0 -200000 0.005 T F;
cycle1_BC: T F T T DP 1.0 2.0;
cycle2_BC: T F T T DP 2.0 1.0;
CALC_BOX: 0.1 0.9;
INIT_CONT_ANG: {RECEDING};
EQUIL_CON_ANG: {WETTABILITY[case.wett]};
Water: {WATER};
Oil: {CO2};
WaterOilInterface: {case.ift};
DRAIN_SINGLETS: T;
RAND_SEED: {case.seed};
RES_FORMAT: Excel;
"""
    (run_dir / "input.dat").write_text(deck)


def read_curve(path: Path) -> List[List[float]]:
    """Read a PNflow Excel-format curve table (Sw, Pc, krw, krg)."""
    rows = []
    for line in path.read_text().splitlines():
        if re.match(r"^\s*-?\d", line) and line.count(",") >= 4:
            rows.append([float(x) for x in line.strip().strip(",").split(",")])
    return rows


def filling_counts(prt: str) -> Dict[str, int]:
    """Throat snap-off / piston counts from the imbibition section."""
    sec = prt.split("Imbibition  =====")[-1].split("Throat Filling Process")[-1]
    out = {}
    for key, pat in (("snap", r"Snap off:\s+(\d+)"), ("piston", r"Piston type displacement:\s+(\d+)")):
        m = re.search(pat, sec)
        out[key] = int(m.group(1)) if m else -1
    return out


def run_case(case: Case, work: Path) -> Optional[Dict[str, object]]:
    """Run one case; keep curves, log and deck; drop the network copy."""
    run_dir = work / case.name
    write_input(case, run_dir)
    try:
        with open(run_dir / "stdout.log", "w") as so, open(run_dir / "stderr.log", "w") as se:
            rc = subprocess.run([str(EXE), "input.dat"], cwd=run_dir, stdout=so, stderr=se,
                                timeout=1200).returncode
    except subprocess.TimeoutExpired:
        logger.error("timeout: %s", case.name)
        return None
    for f in list(run_dir.glob("Berea_*.dat")):
        f.unlink()
    try:
        dr = read_curve(run_dir / "run_cycle1_drain.csv")
        im = read_curve(run_dir / "run_cycle2_imb.csv")
        fill = filling_counts((run_dir / "run_pnflow.prt").read_text())
    except FileNotFoundError:
        logger.error("no output (rc=%s): %s", rc, case.name)
        return None
    brk = next((r[0] for r in dr if r[3] > 0), float("nan"))
    return {
        "case": case.name, "wett": case.wett, "ift_mN_m": case.ift * 1000,
        "pc_max_kPa": case.pc_max_kpa, "column_m": case.pc_max_kpa * 1000 / (DELTA_RHO * 9.81),
        "seed": case.seed, "rc": rc, "Sgi": 1 - dr[-1][0], "Sgr": 1 - im[-1][0],
        "pc_reached_kPa": dr[-1][1] / 1000, "krg_at_Sgi": dr[-1][3], "krw_at_Sgi": dr[-1][2],
        "krw_at_Sgr": im[-1][2], "Sw_breakthrough": brk,
        "throat_snap": fill["snap"], "throat_piston": fill["piston"],
    }


def main() -> None:
    """Run the full matrix and write summary.csv."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    work = Path(sys.argv[1])
    cases = build_cases()
    logger.info("%d cases", len(cases))
    with ThreadPoolExecutor(MAX_WORKERS) as pool:
        results = [r for r in pool.map(lambda c: run_case(c, work), cases) if r]
    with open(work / "summary.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    logger.info("%d/%d completed", len(results), len(cases))


if __name__ == "__main__":
    main()
