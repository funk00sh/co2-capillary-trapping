"""Sgr vs Sgi with Land fits for network A (pnextract) and network B (Dong & Blunt), side by side.

Network B = earlier 306-run matrix (3 angle ranges) + the new 95-115 deg block.
Usage: python plot_trapping_A_vs_B.py <netB_nw_work_dir>
"""
import json
import sys
from pathlib import Path
from typing import Dict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from analyse import COLORS, LAND_SGI_MAX, WETT_LABEL, land

HERE = Path(__file__).resolve().parent


def aggregate(df: pd.DataFrame) -> pd.DataFrame:
    """Seed-mean Sgi and Sgr per angle range and maximum capillary pressure."""
    base = df[df.ift_mN_m.round(2) == 37.37]
    return base.groupby(["wett", "pc_max_kPa"]).agg(Sgi=("Sgi", "mean"), Sgi_sd=("Sgi", "std"),
                                                    Sgr=("Sgr", "mean"), Sgr_sd=("Sgr", "std"),
                                                    n=("Sgr", "size")).reset_index()


def fit_land(agg: pd.DataFrame) -> Dict[str, Dict[str, float]]:
    """Land coefficient per angle range, fitted on Sgi <= LAND_SGI_MAX."""
    out = {}
    for w, g in agg.groupby("wett"):
        g = g[g.Sgi <= LAND_SGI_MAX]
        (c,), cov = curve_fit(land, g.Sgi.values, g.Sgr.values, p0=[1.0])
        out[w] = {"land_C": float(c), "se": float(np.sqrt(cov[0, 0])),
                  "rmse": float(np.sqrt(np.mean((land(g.Sgi.values, c) - g.Sgr.values) ** 2)))}
    return out


def panel(ax, agg: pd.DataFrame, fits: Dict[str, Dict[str, float]], title: str) -> None:
    """One trapping-curve panel."""
    x = np.linspace(0, 1, 100)
    ax.plot(x, x, "k:", lw=0.8)
    ax.plot(x, land(x, 1.0), color="grey", lw=1.4, ls="-.", label="Land C = 1")
    for w in WETT_LABEL:
        g = agg[agg.wett == w]
        lo, hi = g[g.Sgi <= LAND_SGI_MAX], g[g.Sgi > LAND_SGI_MAX]
        ax.errorbar(lo.Sgi, lo.Sgr, yerr=lo.Sgr_sd, fmt="o", ms=4, color=COLORS[w], capsize=2,
                    label=f"{WETT_LABEL[w]}  C = {fits[w]['land_C']:.2f}")
        ax.errorbar(hi.Sgi, hi.Sgr, yerr=hi.Sgr_sd, fmt="o", ms=4, mfc="none", color=COLORS[w], capsize=2, alpha=0.6)
        xf = x[x <= LAND_SGI_MAX]
        ax.plot(xf, land(xf, fits[w]["land_C"]), color=COLORS[w], lw=1.3)
    ax.axvspan(LAND_SGI_MAX, 1, color="grey", alpha=0.08)
    ax.text(0.9, 0.025, "not fitted", ha="center", fontsize=9, color="grey")
    ax.set(xlabel="Initial CO₂ saturation, Sgi", xlim=(0, 1), ylim=(0, 0.7), title=title)
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    ax.grid(alpha=0.3)


def main() -> None:
    """Build the two-panel figure and the comparison table."""
    old = HERE / "old_network_DongBlunt"
    df_b = pd.concat([pd.read_csv(old / "summary_all_runs.csv"), pd.read_csv(Path(sys.argv[1]) / "summary.csv")],
                     ignore_index=True)
    agg_b = aggregate(df_b)
    fits_b = fit_land(agg_b)
    agg_a = pd.read_csv(HERE / "trapping_table.csv")
    fits_a = {w: {"land_C": v["land_C"]} for w, v in json.loads((HERE / "fits.json").read_text())["land"].items()}
    agg_b.to_csv(HERE / "trapping_table_networkB_4angles.csv", index=False, float_format="%.4f")
    (HERE / "fits_networkB_4angles.json").write_text(json.dumps(fits_b, indent=2))

    plt.rcParams.update({"font.size": 11, "axes.titlesize": 12.5})
    fig, axs = plt.subplots(1, 2, figsize=(12.4, 5.0), sharey=True)
    panel(axs[0], agg_a, fits_a, "Network A (pnextract re-extraction)")
    panel(axs[1], agg_b, fits_b, "Network B (Dong & Blunt)")
    axs[0].set_ylabel("Residual CO₂ saturation, Sgr")
    fig.tight_layout()
    fig.savefig(HERE / "fig10_trapping_networkA_vs_B.png", dpi=200)

    rows = []
    for name, agg, fits in (("A", agg_a, fits_a), ("B", agg_b, fits_b)):
        for w in WETT_LABEL:
            g = agg[agg.wett == w].sort_values("Sgi")
            sgr75 = float(np.interp(0.75, g.Sgi, g.Sgr))
            rows.append(dict(network=name, angle=WETT_LABEL[w], land_C=fits[w]["land_C"], Sgr_at_Sgi_0p75=sgr75,
                             Sgr_max=float(g.Sgr.max()), Sgi_at_Sgr_max=float(g.Sgi[g.Sgr.idxmax()])))
    t = pd.DataFrame(rows)
    t.to_csv(HERE / "trapping_networkA_vs_B_table.csv", index=False, float_format="%.3f")
    print(t.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
