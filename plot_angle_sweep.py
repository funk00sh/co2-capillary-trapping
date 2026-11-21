"""CO2 imbibition relative permeability versus advancing contact angle (drainage to 5 kPa), fitted curves only.

Usage: python plot_angle_sweep.py <sweep_dir> [pc_kPa]
Shows fitted CO2 imbibition curves for each angle plus the two (shared) drainage curves; no water imbibition curve.
"""
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from analyse import corey_fit
from run_angle_sweep import ANGLES
from run_matrix import SEEDS, read_curve

logger = logging.getLogger(__name__)
PC = 5.0
HERE = Path(__file__).resolve().parent
RED_D, BLUE = "#d62728", "#1f77b4"


def load(work: Path, angle: int, seed: int) -> Tuple[np.ndarray, np.ndarray]:
    """Drainage and imbibition tables (Sw, Pc, krw, krg) for one run."""
    d = work / f"A{angle}_pc{PC:g}_s{seed}"
    return np.array(read_curve(d / "run_cycle1_drain.csv")), np.array(read_curve(d / "run_cycle2_imb.csv"))


def resample(sw: np.ndarray, kr: np.ndarray, lo: float, hi: float, n: int = 150) -> Tuple[np.ndarray, np.ndarray]:
    """Uniform Sw grid for one branch (tables are densest near the endpoints)."""
    o = np.argsort(sw)
    x = np.linspace(lo, hi, n)
    return x, np.interp(x, sw[o], kr[o])


def fit_run(dr: np.ndarray, im: np.ndarray) -> Dict[str, float]:
    """Corey fits for CO2 drainage and imbibition and a shifted-Corey fit for water drainage, one run."""
    s0, s_end = dr[:, 0].min(), im[:, 0].max()
    gas = dr[dr[:, 3] > 0]
    s_bt = gas[:, 0].max()
    x, y = resample(gas[:, 0], gas[:, 3], s0, s_bt)
    ng_d, kg0, _ = corey_fit(x, y, s0, s_bt, False)
    x, y = resample(im[:, 0], im[:, 3], s0, s_end)
    ng_i, kg_i, _ = corey_fit(x, y, s0, s_end, False)
    x, y = resample(dr[:, 0], dr[:, 2], s0, 1.0)
    kw0 = float(y[0])
    se = (x - s0) / (1 - s0)
    (nw,), _ = curve_fit(lambda s, n: kw0 + (1 - kw0) * s ** n, se, y, p0=[2.0], bounds=(0.3, 12))
    return dict(s0=s0, s_end=s_end, s_bt=s_bt, ng_d=ng_d, kg0=kg0, ng_i=ng_i, kg_i=kg_i, nw=float(nw), kw0=kw0)


def mean_fit(work: Path, angle: int) -> Dict[str, float]:
    """Mean over seeds of the per-run fit parameters."""
    fits = [fit_run(*load(work, angle, s)) for s in SEEDS]
    return {k: float(np.mean([f[k] for f in fits])) for k in fits[0]}


def main() -> None:
    """Fit, tabulate and plot."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    global PC
    work = Path(sys.argv[1])
    PC = float(sys.argv[2]) if len(sys.argv) > 2 else PC
    # drainage must be identical across angles for the same seed (same receding angles and seed)
    diffs = []
    for s in SEEDS:
        ref = load(work, ANGLES[0], s)[0]
        for a in ANGLES[1:]:
            d = load(work, a, s)[0]
            diffs.append(np.abs(ref - d).max() if ref.shape == d.shape else np.inf)
    logger.info("max drainage-table difference across angles (same seed): %.3g", max(diffs))

    fits = {a: mean_fit(work, a) for a in ANGLES}
    table = []
    for a in ANGLES:
        runs = [fit_run(*load(work, a, s)) for s in SEEDS]
        sgr = [1 - r["s_end"] for r in runs]
        table.append(dict(advancing_angle_deg=a, Sgi=1 - fits[a]["s0"], Sgr_mean=np.mean(sgr), Sgr_sd=np.std(sgr, ddof=1),
                          krg_start=fits[a]["kg_i"], n_gas_imbibition=fits[a]["ng_i"], Sw_end=fits[a]["s_end"]))
    pd.DataFrame(table).to_csv(HERE / f"angle_sweep_table_pc{PC:g}.csv", index=False, float_format="%.4f")
    (HERE / f"angle_sweep_fits_pc{PC:g}.json").write_text(json.dumps(fits, indent=2))
    logger.info("\n%s", pd.DataFrame(table).round(3).to_string(index=False))

    # Figure in the style of fig6 (full 0-1 axes, red = CO2, blue = water, solid = drainage, dashed = imbibition)
    ref = fits[ANGLES[0]]
    s0 = ref["s0"]
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    x = np.linspace(s0, 1, 200)
    se = (x - s0) / (1 - s0)
    ax.plot(x, ref["kw0"] + (1 - ref["kw0"]) * se ** ref["nw"], color=BLUE, lw=2.4, label=f"Water, drainage  (n = {ref['nw']:.2f})")
    x = np.linspace(s0, ref["s_bt"], 200)
    ax.plot(x, ref["kg0"] * (1 - (x - s0) / (ref["s_bt"] - s0)) ** ref["ng_d"], color=RED_D, lw=2.4,
            label=f"CO₂, drainage  (n = {ref['ng_d']:.2f})")
    shades = plt.cm.Reds(np.linspace(0.3, 0.95, len(ANGLES)))
    for a, c, t in zip(ANGLES, shades, table):
        f = fits[a]
        x = np.linspace(s0, f["s_end"], 200)
        ax.plot(x, f["kg_i"] * (1 - (x - s0) / (f["s_end"] - s0)) ** f["ng_i"], color=c, lw=2.2, ls="--",
                label=f"CO₂, imbibition, θa = {a}°  (n = {f['ng_i']:.2f}, Sgr = {t['Sgr_mean']:.2f})")
    ax.axvline(s0, color="grey", lw=0.6, ls=":")
    ax.text(s0 + 0.01, 0.62, f"Swi = {s0:.2f}", fontsize=8, color="grey")
    ax.set(xlabel="Water saturation, Sw", ylabel="Relative permeability", xlim=(0, 1), ylim=(0, 1),
           title=f"Berea CO₂–water, fitted Corey curves · drainage to {PC:g} kPa (Sgi = {1 - s0:.2f})\nCO₂ imbibition for different advancing angles θa")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(HERE / f"fig7_co2_imbibition_vs_angle_pc{PC:g}kPa.png", dpi=200)


if __name__ == "__main__":
    main()
