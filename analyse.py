"""Analyse the Berea Pc-ladder matrix: trapping curves, Land and Corey fits, figures.

Usage: python analyse.py <work_dir>   (figures and tables are written next to this file)
"""
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from run_matrix import DELTA_RHO, read_curve

logger = logging.getLogger(__name__)
HERE = Path(__file__).resolve().parent
WETT_LABEL = {"WW_20-40": "θa 20–40°", "MW_45-65": "θa 45–65°", "IW_70-90": "θa 70–90°",
              "NW_95-115": "θa 95–115°"}
COLORS = {"WW_20-40": "#1f77b4", "MW_45-65": "#2ca02c", "IW_70-90": "#ff7f0e", "NW_95-115": "#d62728"}
BASE = "MW_45-65"
LOOP_PC = 5.8  # kPa, main hysteresis-loop case
SCAN_PCS = [4.8, 5.0, 5.2, 5.8, 10.0]
PC_TOP = 30
LAND_SGI_MAX = 0.80  # Land fit uses the rising part of the trapping curve only
BLUE, RED = "#1f77b4", "#d62728"


def land(sgi: np.ndarray, c: float) -> np.ndarray:
    """Land (1968) trapping model."""
    return sgi / (1 + c * sgi)


def curves(work: Path, case: str) -> Tuple[np.ndarray, np.ndarray]:
    """Drainage and imbibition tables (Sw, Pc, krw, krg) for one run."""
    d = work / case
    return np.array(read_curve(d / "run_cycle1_drain.csv")), np.array(read_curve(d / "run_cycle2_imb.csv"))


def seed_mean(tabs: List[np.ndarray], n: int = 120) -> np.ndarray:
    """Average several curves (Sw, krw, krg) on the Sw interval they all cover."""
    lo = max(t[:, 0].min() for t in tabs)
    hi = min(t[:, 0].max() for t in tabs)
    sw = np.linspace(lo, hi, n)
    cols = []
    for k in (2, 3):
        vals = []
        for t in tabs:
            o = np.argsort(t[:, 0])
            vals.append(np.interp(sw, t[o, 0], t[o, k]))
        cols.append(np.mean(vals, axis=0))
    return np.column_stack([sw, cols[0], cols[1]])


def corey_fit(sw: np.ndarray, kr: np.ndarray, s_lo: float, s_hi: float,
              rising: bool) -> Tuple[float, float, np.ndarray]:
    """Fit kr = k_end * Se^n between the immobile (s_lo) and endpoint (s_hi) saturations."""
    se = np.clip((sw - s_lo) / (s_hi - s_lo), 0, 1)
    if not rising:
        se = 1 - se
    k_end = float(kr[np.argmax(se)])
    (n,), _ = curve_fit(lambda x, n: k_end * x ** n, se, kr, p0=[2.0], bounds=(0.3, 12))
    return float(n), k_end, k_end * se ** n


def loop_data(work: Path, df: pd.DataFrame, wett: str, pc: float) -> Tuple[tuple, tuple, np.ndarray, np.ndarray]:
    """Raw seed curves and seed means for one wettability and Pc level."""
    sel = df[(df.wett == wett) & (df.pc_max_kPa == pc) & (df.ift_mN_m.round(2) == 37.37)]
    drs, ims = zip(*[curves(work, c) for c in sel.case])
    return drs, ims, seed_mean(list(drs)), seed_mean(list(ims))


def plot_loop(axes, work: Path, df: pd.DataFrame, wett: str, pc: float) -> Dict[str, float]:
    """kr hysteresis loop: seeds (thin), seed mean (thick), Corey fits (black dotted)."""
    drs, ims, md, mi = loop_data(work, df, wett, pc)
    swi, sw_end = md[0, 0], mi[-1, 0]
    gas = md[md[:, 2] > 0]
    sw_bt = gas[-1, 0]
    fits = {"Swi": float(swi), "Sgr": float(1 - sw_end), "Sgc_drainage": float(1 - sw_bt)}
    fit_lines = []
    n, _, y = corey_fit(md[:, 0], md[:, 1], swi, 1.0, True)
    fits.update(nw_drainage=n)
    fit_lines.append((md[:, 0], y))
    n, k, y = corey_fit(gas[:, 0], gas[:, 2], swi, sw_bt, False)
    fits.update(ng_drainage=n, krg_at_Swi=k)
    fit_lines.append((gas[:, 0], y))
    n, k, y = corey_fit(mi[:, 0], mi[:, 1], swi, sw_end, True)
    fits.update(nw_imbibition=n, krw_at_Sgr=k)
    fit_lines.append((mi[:, 0], y))
    n, _, y = corey_fit(mi[:, 0], mi[:, 2], swi, sw_end, False)
    fits.update(ng_imbibition=n)
    fit_lines.append((mi[:, 0], y))
    for ax in axes:
        for t in drs + ims:
            ax.plot(t[:, 0], t[:, 2], color=BLUE, lw=0.5, alpha=0.2)
            ax.plot(t[:, 0], t[:, 3], color=RED, lw=0.5, alpha=0.2)
        ax.plot(md[:, 0], md[:, 1], color=BLUE, lw=2.2, label="Water, drainage")
        ax.plot(md[:, 0], md[:, 2], color=RED, lw=2.2, label="CO₂, drainage")
        ax.plot(mi[:, 0], mi[:, 1], color=BLUE, lw=2.2, ls="--", label="Water, imbibition")
        ax.plot(mi[:, 0], mi[:, 2], color=RED, lw=2.2, ls="--", label="CO₂, imbibition")
        for i, (x, y) in enumerate(fit_lines):
            ax.plot(x, y, color="k", lw=0.9, ls=":", label="Corey fit" if i == 0 else None)
        ax.set(xlabel="Water saturation, Sw", ylabel="Relative permeability", xlim=(0, 1))
        ax.grid(alpha=0.3)
    axes[0].set_ylim(0, 1)
    axes[0].legend(frameon=False, fontsize=9)
    axes[1].set_yscale("log")
    axes[1].set_ylim(1e-4, 1)
    return fits


def main() -> None:
    """Build tables and figures."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    work = Path(sys.argv[1])
    df = pd.read_csv(work / "summary.csv")
    df.to_csv(HERE / "summary_all_runs.csv", index=False)
    base = df[df.ift_mN_m.round(2) == 37.37]
    agg = base.groupby(["wett", "pc_max_kPa"]).agg(
        column_m=("column_m", "first"), Sgi=("Sgi", "mean"), Sgi_sd=("Sgi", "std"),
        Sgr=("Sgr", "mean"), Sgr_sd=("Sgr", "std"), krg_at_Sgi=("krg_at_Sgi", "mean"),
        krw_at_Sgr=("krw_at_Sgr", "mean"), n=("Sgr", "size")).reset_index()
    agg.to_csv(HERE / "trapping_table.csv", index=False, float_format="%.4f")

    fits: Dict[str, Dict[str, float]] = {}
    for w, g in agg.groupby("wett"):
        g = g[g.Sgi <= LAND_SGI_MAX]
        (c,), cov = curve_fit(land, g.Sgi.values, g.Sgr.values, p0=[1.0])
        rmse = float(np.sqrt(np.mean((land(g.Sgi.values, c) - g.Sgr.values) ** 2)))
        fits[w] = {"land_C": float(c), "land_C_se": float(np.sqrt(cov[0, 0])), "rmse": rmse,
                   "Sgr_max": float(g.Sgr.max()), "Sgi_at_max": float(g.Sgi[g.Sgr.idxmax()])}

    # Fig 1: kr hysteresis loop with Corey fits
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    corey = plot_loop(ax, work, df, BASE, LOOP_PC)
    fig.suptitle(f"Berea CO₂–water relative permeability · {WETT_LABEL[BASE]} · drainage to {LOOP_PC:g} kPa · "
                 f"Sgi {1 - corey['Swi']:.2f} → Sgr {corey['Sgr']:.2f} · mean of 5 seeds", fontsize=10)
    fig.tight_layout()
    fig.savefig(HERE / "fig1_kr_hysteresis.png", dpi=200)
    (HERE / "fits.json").write_text(json.dumps({"land": fits, "corey_base_case": corey}, indent=2))
    logger.info(json.dumps({"land": fits, "corey_base_case": corey}, indent=2))

    # Fig 2: Pc hysteresis with scanning curves
    fig, ax = plt.subplots(figsize=(6.5, 4.6))
    cmap = plt.cm.viridis(np.linspace(0.1, 0.9, len(SCAN_PCS)))
    dr, _ = curves(work, f"{BASE}_ift37_pc{PC_TOP:g}_s50")
    ax.plot(dr[:, 0], dr[:, 1] / 1000, "k-", lw=2, label="primary drainage")
    for pc, col in zip(SCAN_PCS, cmap):
        _, im = curves(work, f"{BASE}_ift37_pc{pc:g}_s50")
        ax.plot(im[:, 0], im[:, 1] / 1000, color=col, ls="--", label=f"imbibition from {pc:g} kPa")
    ax.axhline(0, color="grey", lw=0.6)
    ax.set(xlabel="Water saturation, Sw", ylabel="Pc = P(CO₂) − P(water), kPa", xlim=(0, 1), ylim=(-4, 15),
           title=f"Capillary-pressure hysteresis · {WETT_LABEL[BASE]}")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(HERE / "fig2_pc_hysteresis.png", dpi=200)

    # Fig 3: trapping curves with Land fits
    fig, ax = plt.subplots(figsize=(6.8, 5))
    x = np.linspace(0, 1, 100)
    ax.plot(x, x, "k:", lw=0.8, label="100% trapped")
    ax.plot(x, land(x, 1.0), color="grey", lw=1.4, ls="-.", label="Berea experiments, Land C ≈ 1 (Krevor et al. 2012)")
    for w in WETT_LABEL:
        g = agg[agg.wett == w]
        lo, hi = g[g.Sgi <= LAND_SGI_MAX], g[g.Sgi > LAND_SGI_MAX]
        ax.errorbar(lo.Sgi, lo.Sgr, xerr=lo.Sgi_sd, yerr=lo.Sgr_sd, fmt="o", ms=4, color=COLORS[w], capsize=2,
                    label=f"{WETT_LABEL[w]}: Land C = {fits[w]['land_C']:.2f}")
        ax.errorbar(hi.Sgi, hi.Sgr, yerr=hi.Sgr_sd, fmt="o", ms=4, mfc="none", color=COLORS[w], capsize=2, alpha=0.6)
        xf = x[x <= LAND_SGI_MAX]
        ax.plot(xf, land(xf, fits[w]["land_C"]), color=COLORS[w], lw=1)
    ax.axvspan(LAND_SGI_MAX, 1, color="grey", alpha=0.08)
    ax.text(0.9, 0.03, "Sgi > 0.8: not fitted", ha="center", fontsize=8, color="grey")
    ax.set(xlabel="Initial CO₂ saturation, Sgi", ylabel="Residual CO₂ saturation, Sgr", xlim=(0, 1), ylim=(0, 0.6),
           title="Trapping curve: pore-network model (5 seeds) and Land fits")
    
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(HERE / "fig3_trapping_curve.png", dpi=200)

    # Fig 4: wettability comparison
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    for w in WETT_LABEL:
        _, _, _, mi = loop_data(work, df, w, LOOP_PC)
        ax[0].plot(mi[:, 0], mi[:, 2], color=COLORS[w], label=WETT_LABEL[w])
        ax[1].plot(mi[:, 0], mi[:, 1], color=COLORS[w], label=WETT_LABEL[w])
        g = agg[agg.wett == w]
        ax[2].plot(g.Sgi, g.Sgr / g.Sgi, "o-", ms=4, color=COLORS[w], label=WETT_LABEL[w])
    _, _, md, _ = loop_data(work, df, BASE, LOOP_PC)
    ax[0].plot(md[:, 0], md[:, 2], "k-", lw=1, label="drainage")
    ax[1].plot(md[:, 0], md[:, 1], "k-", lw=1, label="drainage")
    ax[0].set(xlabel="Sw", ylabel="CO₂ relative permeability", xlim=(0, 1), ylim=(0, 0.7),
              title=f"Imbibition krg after drainage to {LOOP_PC:g} kPa")
    ax[1].set(xlabel="Sw", ylabel="Water relative permeability", xlim=(0, 1), ylim=(0, 0.7),
              title="Imbibition krw")
    ax[2].set(xlabel="Sgi", ylabel="Trapping efficiency, Sgr / Sgi", xlim=(0, 1), ylim=(0, 1),
              title="Fraction of injected CO₂ trapped")
    for a in ax:
        a.grid(alpha=0.3)
        a.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / "fig4_wettability.png", dpi=200)

    # Fig 5: IFT
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    ift_df = df[(df.wett == BASE) & (df.seed == 50)]
    for (ift, g), col in zip(ift_df.groupby("ift_mN_m"), ["#9467bd", "#2ca02c", "#8c564b"]):
        g = g.sort_values("pc_max_kPa")
        lab = f"IFT {ift:.1f} mN/m"
        dr, _ = curves(work, g.case.iloc[-1])
        ax[0].plot(dr[:, 0], dr[:, 1] / 1000, color=col, label=lab)
        ax[1].plot(g.pc_max_kPa, g.Sgi, "o-", ms=4, color=col, label=lab)
        ax[2].plot(g.Sgi, g.Sgr, "o", ms=5, mfc="none", color=col, label=lab)
    ax[0].set(xlabel="Sw", ylabel="Drainage Pc, kPa", xlim=(0, 1), ylim=(0, 15), title="Pc scales with IFT")
    ax[1].set(xlabel="Applied maximum Pc, kPa", ylabel="Sgi reached", xscale="log", ylim=(0, 1))
    ax[1].set_xticks([5, 6, 8, 10, 15, 20, 30])
    ax[1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax[1].minorticks_off()
    k = 1000 / (DELTA_RHO * 9.81)
    ax2 = ax[1].secondary_xaxis("top", functions=(lambda p: p * k, lambda h: h / k))
    ax2.set_xticks([0.8, 1, 1.5, 2, 3, 5])
    ax2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax2.minorticks_off()
    ax2.set_xlabel("Same Pc → different Sgi · equivalent CO₂ column height, m", fontsize=9)
    ax[2].set(xlabel="Sgi", ylabel="Sgr", xlim=(0, 1), ylim=(0, 0.6), title="Sgr(Sgi) is independent of IFT")
    for a in ax:
        a.grid(alpha=0.3)
        a.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / "fig5_ift.png", dpi=200)


if __name__ == "__main__":
    main()
