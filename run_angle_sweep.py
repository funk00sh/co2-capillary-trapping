"""Single-angle sweep: drainage to 5 kPa, then imbibition with a constant advancing angle.

Usage: python run_angle_sweep.py <work_dir> [pc_kPa]
Same network, fluids, receding angle (0-10 deg), IFT and seeds as the 462-run matrix.
"""
import csv
import logging
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from run_matrix import SEEDS, WETTABILITY, Case, run_case

logger = logging.getLogger(__name__)
ANGLES = [20, 40, 60, 80, 100, 120]
PC_KPA = 5.0
IFT = 0.03737


def main() -> None:
    """Run angles x seeds and write summary.csv."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    work = Path(sys.argv[1])
    work.mkdir(parents=True, exist_ok=True)
    pc = float(sys.argv[2]) if len(sys.argv) > 2 else PC_KPA
    for a in ANGLES:
        WETTABILITY[f"A{a}"] = f"4 {a} {a} -0.2 -3 rand 0"
    cases = [Case(f"A{a}_pc{pc:g}_s{s}", f"A{a}", IFT, pc, s) for a in ANGLES for s in SEEDS]
    with ThreadPoolExecutor(8) as pool:
        results = [r for r in pool.map(lambda c: run_case(c, work), cases) if r]
    with open(work / "summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0]))
        w.writeheader()
        w.writerows(results)
    logger.info("%d/%d completed", len(results), len(cases))
    for r in results[::5]:
        logger.info("%s Sgi %.3f Sgr %.3f rc %s", r["case"], r["Sgi"], r["Sgr"], r["rc"])


if __name__ == "__main__":
    main()
