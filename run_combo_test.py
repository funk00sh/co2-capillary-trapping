"""Can network A behave 'normally' at a realistic permeability? Combined B-like changes to A with K held near the reference.

Changes (on A): throat radius x f, throat G x g, throat volume share -> s (porosity conserved), and internal throat length
scaled so the absolute permeability stays at K_TARGET (B's value). The fill sequence depends only on radii, G, topology and
angles; throat length and volumes only affect conductance and saturation weights.
Usage: python run_combo_test.py <work_dir>
"""
import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

import run_matrix as rm
from make_variant_network import build_combo
from parse_filling_events import parse
from run_geometry_sweep import SEEDS, SRC, TARGET, deck, run_pnflow

K_TARGET = 1144.0  # mD, network B
B_SHARE = 0.112
VARIANTS = [  # name, radius factor, G factor, throat volume share (None = unchanged)
    ("K-matched control", 1.0, 1.0, None),
    ("volume share 11%", 1.0, 1.0, B_SHARE),
    ("radius x0.78 + share 11%", 0.78, 1.0, B_SHARE),
    ("radius x0.78 + G x0.69 + share 11% (B-like)", 0.78, 0.69, B_SHARE),
    ("radius x0.65 + share 11%", 0.65, 1.0, B_SHARE),
    ("radius x0.65 + G x0.69 + share 11%", 0.65, 0.69, B_SHARE),
    ("radius x0.55 + share 11%", 0.55, 1.0, B_SHARE),
]


def k_of(prt: str) -> float:
    """Absolute permeability (mD) from the log."""
    return float(re.search(r"Absolute permeability:.*?,\s*([\d.eE+-]+)\s*\(mD\)", prt).group(1))


def tune(v, work: Path) -> dict:
    """Find the throat-length factor that gives K_TARGET, then calibrate Pc for Sgi = TARGET."""
    name, f, g, s = v
    tag = re.sub(r"[^A-Za-z0-9]+", "_", name)
    c = 1.0
    for it in range(4):
        st = build_combo(SRC["A"], work / "nets" / tag, f, g, s, c)
        res = run_pnflow(work / "nets" / tag, work / "cal" / f"{tag}_k{it}", 40.0 if it == 3 else 0.6, "WW_20-40", 11)
        k = k_of(res[2])
        if abs(k / K_TARGET - 1) < 0.03:
            break
        c *= k / K_TARGET  # K is ~ 1 / throat length when throats dominate the resistance
    st = build_combo(SRC["A"], work / "nets" / tag, f, g, s, c)
    res = run_pnflow(work / "nets" / tag, work / "cal" / f"{tag}_final", 40.0, "WW_20-40", 11)
    k = k_of(res[2])
    dr = np.array(res[0])[:, :2]
    sgi, pc = 1 - dr[:, 0], dr[:, 1] / 1000
    target = min(TARGET, float(sgi.max()) - 0.03)
    i = int(np.argmax(sgi >= target))
    pc_t = float(np.interp(target, sgi[i - 1:i + 1], pc[i - 1:i + 1]))
    return dict(name=name, tag=tag, f=f, g=g, share=s, c=c, K=k, pc_kpa=round(pc_t, 3), sgi_target=target, stats=st)


def one(job) -> dict:
    """One production run."""
    v, wett, seed, work = job
    rd = work / "runs" / f"{v['tag']}_{wett}_s{seed}"
    res = run_pnflow(work / "nets" / v["tag"], rd, v["pc_kpa"], wett, seed)
    dr, im, prt = res
    return dict(name=v["name"], tag=v["tag"], wett=wett, seed=seed, K=v["K"], Sgi=1 - dr[-1][0], Sgr=1 - im[-1][0], **parse(prt))


def main() -> None:
    """Tune, run, save."""
    work = Path(sys.argv[1])
    work.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(7) as pool:
        tuned = list(pool.map(lambda v: tune(v, work), VARIANTS))
    for t in tuned:
        print(t["name"], "c=%.3f K=%.0f mD pc=%.2f kPa" % (t["c"], t["K"], t["pc_kpa"]), flush=True)
    (work / "combo_tuned.json").write_text(json.dumps(tuned, indent=1))
    jobs = [(t, w, s, work) for t in tuned for w in rm.WETTABILITY for s in SEEDS]
    with ThreadPoolExecutor(8) as pool:
        rows = list(pool.map(one, jobs))
    (work / "combo_results.json").write_text(json.dumps(rows, indent=1))
    print("finished", len(rows), flush=True)


if __name__ == "__main__":
    main()
