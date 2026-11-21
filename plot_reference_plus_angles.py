"""Reference case (fig6: theta_a 45-65 deg, drainage to 5.8 kPa) with extra CO2 imbibition curves for other angles.

Usage: python plot_reference_plus_angles.py <matrix_dir> <sweep_dir>
The reference curves are drawn with exactly the fig6 fitting code; the added angle curves use the same method
(Corey fit to the 5-seed mean curve) on the single-angle sweep runs. No water imbibition curve is drawn.
"""
import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyse import BASE, BLUE, LOOP_PC, RED
from plot_fit_only import curve, fit_branches
from run_angle_sweep import ANGLES
from run_matrix import SEEDS

logger = logging.getLogger(__name__)
HERE = Path(__file__).resolve().parent


def sweep_frame(angle: int) -> pd.DataFrame:
    """Run table for one sweep angle in the format fit_branches expects."""
    return pd.DataFrame([dict(case=f"A{angle}_pc{LOOP_PC:g}_s{s}", wett=f"A{angle}", pc_max_kPa=LOOP_PC,
                              ift_mN_m=37.37, seed=s) for s in SEEDS])


def main() -> None:
    """Draw the reference figure plus the added angle curves."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    matrix, sweep = Path(sys.argv[1]), Path(sys.argv[2])
    ref = fit_branches(matrix, pd.read_csv(matrix / "summary.csv"), BASE, LOOP_PC)  # identical to fig6
    swi, sw_end = ref["krw_drainage"]["lo"], ref["krw_imbibition"]["hi"]

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    x, y = curve(ref["krw_drainage"])
    ax.plot(x, y, color=BLUE, ls="-", lw=2.4, label=f"Water, drainage  (n = {ref['krw_drainage']['n']:.2f})")
    x, y = curve(ref["krg_drainage"])
    ax.plot(x, y, color=RED, ls="-", lw=2.4, label=f"CO₂, drainage  (n = {ref['krg_drainage']['n']:.2f})")
    x, y = curve(ref["krg_imbibition"])
    ax.plot(x, y, color=RED, ls="--", lw=2.4,
            label=f"CO₂, imbibition, θa 45–65° (reference)  (n = {ref['krg_imbibition']['n']:.2f}, Sgr = {1 - sw_end:.2f})")
    ax.axvline(swi, color="grey", lw=0.6, ls=":")
    ax.axvline(sw_end, color="grey", lw=0.6, ls=":")
    ax.text(swi + 0.01, 0.62, f"Swi = {swi:.2f}", fontsize=8, color="grey")
    ax.text(sw_end + 0.01, 0.62, f"1 − Sgr = {sw_end:.2f}", fontsize=8, color="grey")

    cols = plt.cm.viridis(np.linspace(0.05, 0.92, len(ANGLES)))
    rows = []
    for a, c in zip(ANGLES, cols):
        p = fit_branches(sweep, sweep_frame(a), f"A{a}", LOOP_PC)["krg_imbibition"]
        x, y = curve(p)
        ax.plot(x, y, color=c, ls="--", lw=1.8, label=f"CO₂, imbibition, θa = {a}°  (n = {p['n']:.2f}, Sgr = {1 - p['hi']:.2f})")
        rows.append(dict(angle_deg=a, corey_n=p["n"], krg_start=p["k"], Sw_end=p["hi"], Sgr=1 - p["hi"]))
    pd.DataFrame(rows).to_csv(HERE / "reference_plus_angles_table.csv", index=False, float_format="%.4f")
    logger.info("reference: n_gas_imb %.3f Sgr %.3f\n%s", ref["krg_imbibition"]["n"], 1 - sw_end, pd.DataFrame(rows).round(3).to_string(index=False))

    ax.set(xlabel="Water saturation, Sw", ylabel="Relative permeability", xlim=(0, 1), ylim=(0, 1),
           title=f"Berea CO₂–water, fitted Corey curves · drainage to {LOOP_PC:g} kPa\nreference θa 45–65° with CO₂ imbibition at other angles")
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(HERE / "fig8_reference_plus_angle_curves.png", dpi=200)


if __name__ == "__main__":
    main()
