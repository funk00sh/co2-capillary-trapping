"""Angle sweep on the older Dong & Blunt Berea network (6,298 pores) for comparison with the pnextract network.

Usage: python run_angle_sweep_oldnet.py <work_dir> [pc_kPa]   (default 8 kPa, which gives Sgi ~ 0.75 on this network)
"""
import sys

import run_angle_sweep
import run_matrix

OLD_NET = run_matrix.HERE / "network_B"

if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.argv.append("8.0")
    run_matrix.NET = OLD_NET
    run_angle_sweep.main()
