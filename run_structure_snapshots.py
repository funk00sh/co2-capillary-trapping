"""Re-run the 4 angle ranges on networks A and B at matched Sgi ~ 0.75 with per-element state output.

Usage: python run_structure_snapshots.py <work_dir>
A: Pc 5.8 kPa, B: Pc 8.0 kPa (both give Sgi ~ 0.75); seeds 11, 37, 50. Same deck as the production matrix.
"""
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import run_matrix as rm

NETS = {"A": (rm.HERE / "network_pnextract", 5.8),
        "B": (rm.HERE / "network_B", 8.0)}
SEEDS = [11, 37, 50]


def run_one(job) -> str:
    """Run one case with visualisation output."""
    net, wett, seed, work = job
    src, pc = NETS[net]
    d = work / f"{net}_{wett}_s{seed}"
    d.mkdir(parents=True, exist_ok=True)
    for f in list(src.glob("Berea_node*.dat")) + list(src.glob("Berea_link*.dat")):
        shutil.copy(f, d / f.name)
    (d / "input.dat").write_text(f"""TITLE: run;
NETWORK: F Berea;
cycle1: 0.0 {pc * 1000:.1f} 0.005 T F;
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
visuaLight: T T T F F;
""")
    with open(d / "stdout.log", "w") as so, open(d / "stderr.log", "w") as se:
        rc = subprocess.run([str(rm.EXE), "input.dat"], cwd=d, stdout=so, stderr=se, timeout=3000).returncode
    return f"{d.name} rc={rc} vtu={len(list((d / 'run_res').glob('*.vtu')))}"


def main() -> None:
    """Run all 24 cases."""
    work = Path(sys.argv[1])
    jobs = [(n, w, s, work) for n in NETS for w in rm.WETTABILITY for s in SEEDS]
    with ThreadPoolExecutor(8) as pool:
        for r in pool.map(run_one, jobs):
            print(r, flush=True)


if __name__ == "__main__":
    main()
