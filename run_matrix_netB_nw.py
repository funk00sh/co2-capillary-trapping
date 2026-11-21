"""Add the 95-115 deg advancing-angle block to the Dong & Blunt network (B) matrix.

Same deck, Pc ladder and seeds as the earlier 306-run matrix on network B
(old_network_DongBlunt/), so the four angle ranges are directly comparable.
Usage: python run_matrix_netB_nw.py <work_dir>
"""
import csv
import logging
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import run_matrix as rm

logger = logging.getLogger(__name__)
NET_B = rm.HERE / "network_B"
LADDER_B = [5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.3, 7.6, 8.0, 8.5, 9.0, 10.0, 12.0, 15.0, 20.0, 30.0, 50.0]


def main() -> None:
    """Run the NW_95-115 cases on network B and write summary.csv."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    work = Path(sys.argv[1])
    rm.NET = NET_B
    cases = [rm.Case(f"NW_95-115_ift37_pc{pc:g}_s{s}", "NW_95-115", rm.BASE_IFT, pc, s)
             for pc in LADDER_B for s in rm.SEEDS]
    logger.info("%d cases", len(cases))
    with ThreadPoolExecutor(rm.MAX_WORKERS) as pool:
        results = [r for r in pool.map(lambda c: rm.run_case(c, work), cases) if r]
    with open(work / "summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0]))
        w.writeheader()
        w.writerows(results)
    logger.info("%d/%d completed", len(results), len(cases))


if __name__ == "__main__":
    main()
