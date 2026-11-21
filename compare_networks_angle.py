"""Residual CO2 versus advancing angle on two Berea networks (same solver, fluids, seeds, receding angle).

Usage: python compare_networks_angle.py <new_sweep_dir> <old_sweep_dir>
New network (pnextract, 4,670 pores): drainage to 5.8 kPa. Old network (Dong & Blunt, 6,298 pores): 8 kPa (matched Sgi ~ 0.75).
"""
import logging
import sys
from pathlib import Path
from typing import Dict, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyse import RED
from plot_fit_only import curve, fit_branches
from run_angle_sweep import ANGLES
from run_matrix import SEEDS, read_curve

logger = logging.getLogger(__name__)
HERE = Path(__file__).resolve().parent
NETS = {"new": ("pnextract network, 4,670 pores", 5.8, "#1f77b4"), "old": ("2007 Dong & Blunt network, 6,298 pores", 8.0, "#d95f02")}


def sweep_frame(angle: int, pc: float) -> pd.DataFrame:
    """Run table in the format fit_branches expects."""
    return pd.DataFrame([dict(case=f"A{angle}_pc{pc:g}_s{s}", wett=f"A{angle}", pc_max_kPa=pc, ift_mN_m=37.37, seed=s)
                         for s in SEEDS])


def sgr_stats(work: Path, pc: float) -> List[Dict[str, float]]:
    """Mean and sd of Sgi and Sgr over seeds for each angle."""
    out = []
    for a in ANGLES:
        sgi, sgr = [], []
        for s in SEEDS:
            d = work / f"A{a}_pc{pc:g}_s{s}"
            sgi.append(1 - read_curve(d / "run_cycle1_drain.csv")[-1][0])
            sgr.append(1 - read_curve(d / "run_cycle2_imb.csv")[-1][0])
        out.append(dict(angle=a, Sgi=np.mean(sgi), Sgr=np.mean(sgr), Sgr_sd=np.std(sgr, ddof=1)))
    return out


def fan(ax: plt.Axes, work: Path, pc: float, title: str) -> None:
    """Fitted CO2 drainage and imbibition curves for all angles (same method as fig8)."""
    cols = plt.cm.viridis(np.linspace(0.05, 0.92, len(ANGLES)))
    ref = None
    for a, c in zip(ANGLES, cols):
        p = fit_branches(work, sweep_frame(a, pc), f"A{a}", pc)
        ref = ref or p
        x, y = curve(p["krg_imbibition"])
        ax.plot(x, y, color=c, ls="--", lw=1.9, label=f"θa = {a}°  (Sgr = {1 - p['krg_imbibition']['hi']:.2f})")
    x, y = curve(ref["krg_drainage"])
    ax.plot(x, y, color=RED, lw=2.4, label="CO₂, drainage")
    ax.set(xlabel="Water saturation, Sw", ylabel="CO₂ relative permeability", xlim=(0, 1), ylim=(0, 1), title=title)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    ax.grid(alpha=0.3)


def main() -> None:
    """Draw the comparison and save the table."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    dirs = {"new": Path(sys.argv[1]), "old": Path(sys.argv[2])}
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.9))
    rows = []
    for k, (label, pc, col) in NETS.items():
        st = sgr_stats(dirs[k], pc)
        ax[0].errorbar([s["angle"] for s in st], [s["Sgr"] for s in st], yerr=[s["Sgr_sd"] for s in st], marker="o", ms=5,
                       color=col, capsize=3, label=f"{label}\ndrainage to {pc:g} kPa, Sgi = {st[0]['Sgi']:.2f}")
        rows += [dict(network=k, drainage_kPa=pc, **s) for s in st]
    ax[0].set(xlabel="Advancing contact angle θa, °", ylabel="Residual CO₂ saturation, Sgr", xlim=(10, 130), ylim=(0, 0.65),
              title="Sgr versus advancing angle (mean ± sd, 5 seeds)")
    ax[0].axvline(90, color="grey", lw=0.6, ls=":")
    ax[0].legend(frameon=False, fontsize=8, loc="lower left")
    ax[0].grid(alpha=0.3)
    fan(ax[1], dirs["new"], 5.8, "pnextract network (new): Sgi = 0.75")
    fan(ax[2], dirs["old"], 8.0, "Dong & Blunt network (old): Sgi = 0.75")
    fig.tight_layout()
    fig.savefig(HERE / "fig9_network_comparison_angle.png", dpi=200)
    pd.DataFrame(rows).to_csv(HERE / "network_comparison_angle_table.csv", index=False, float_format="%.4f")
    logger.info("\n%s", pd.DataFrame(rows).pivot(index="angle", columns="network", values=["Sgi", "Sgr", "Sgr_sd"]).round(3).to_string())


if __name__ == "__main__":
    main()
