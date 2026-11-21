"""Check that rewritten networks reproduce the original runs and that modified variants load in PNflow."""
import subprocess
import sys
from pathlib import Path

import pandas as pd

import run_matrix as rm
from make_variant_network import build

W = Path(sys.argv[1])
SRC = {"A": rm.HERE / "network_pnextract", "B": rm.HERE / "network_B"}
PC = {"A": 5.8, "B": 8.0}


def go(tag, net, knob, level, seed=11, wett="WW_20-40"):
    d = W / tag
    stats = build(SRC[net], d, knob, level, seed)
    (d / "input.dat").write_text(f"""TITLE: run;
NETWORK: F Berea;
cycle1: 0.0 {PC[net] * 1000:.1f} 0.005 T F;
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
""")
    rc = subprocess.run([str(rm.EXE), "input.dat"], cwd=d, capture_output=True, timeout=1200).returncode
    dr, im = rm.read_curve(d / "run_cycle1_drain.csv"), rm.read_curve(d / "run_cycle2_imb.csv")
    print(tag, "rc", rc, "Sgi %.4f Sgr %.4f" % (1 - dr[-1][0], 1 - im[-1][0]), {k: round(v, 3) for k, v in stats.items()}, flush=True)


ref = pd.concat([pd.read_csv(rm.HERE / "summary_all_runs.csv"), pd.read_csv(rm.HERE / "old_network_DongBlunt/summary_all_runs.csv")])
for n, case in (("A", "WW_20-40_ift37_pc5.8_s11"), ("B", "WW_20-40_ift37_pc8_s11")):
    r = ref[ref.case == case].iloc[0]
    print("REFERENCE", n, "Sgi %.4f Sgr %.4f" % (r.Sgi, r.Sgr))
for n in "AB":
    go(f"t_{n}_base", n, "coord", 0)
    go(f"t_{n}_coord20", n, "coord", 0.2)
    go(f"t_{n}_rad14", n, "radius", 1.4)
    go(f"t_{n}_shape2", n, "shape", 2.0)
