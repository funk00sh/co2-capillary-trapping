"""Base-case kr loop drawn from the Corey fits only (no raw curves, no dotted overlay).

Usage: python plot_fit_only_clean.py <work_dir> [out.png]
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from analyse import BASE, BLUE, LOOP_PC, RED
from plot_fit_only import curve, fit_branches

HERE = Path(__file__).resolve().parent


def main() -> None:
    """Draw the four Corey branches in the style of the base-case loop figure."""
    work = Path(sys.argv[1])
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "fig1b_kr_corey_only.png"
    df = pd.read_csv(work / "summary.csv")
    p = fit_branches(work, df, BASE, LOOP_PC)
    plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
    fig, ax = plt.subplots(figsize=(6.3, 3.9))
    style = [("krw_drainage", BLUE, "-", "Water, drainage"), ("krg_drainage", RED, "-", "CO₂, drainage"),
             ("krw_imbibition", BLUE, "--", "Water, imbibition"), ("krg_imbibition", RED, "--", "CO₂, imbibition")]
    for key, col, ls, lab in style:
        x, y = curve(p[key])
        ax.plot(x, y, color=col, ls=ls, lw=2.4, label=lab)
    ax.set(xlabel="Brine saturation, Sw", ylabel="Relative permeability", xlim=(0, 1), ylim=(0, 1))
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(out, dpi=185)


if __name__ == "__main__":
    main()
