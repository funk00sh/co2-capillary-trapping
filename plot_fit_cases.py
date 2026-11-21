"""One fitted-only relative permeability figure per run (seed-mean Corey fits).

Usage: python plot_fit_cases.py <work_dir>   (writes kr_fits_per_run/*.png and corey_parameters_all.csv)
"""
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

import numpy as np

from analyse import BLUE, RED, WETT_LABEL, corey_fit, curves
from plot_fit_only import curve

logger = logging.getLogger(__name__)
HERE = Path(__file__).resolve().parent
OUT = HERE / "kr_fits_per_run"
CASES: List[Tuple[str, float]] = [("MW_45-65", 5.8), ("MW_45-65", 7.0), ("MW_45-65", 10.0),
                                  ("WW_20-40", 5.8), ("IW_70-90", 5.8), ("NW_95-115", 5.8)]


def grid(sw: np.ndarray, kr: np.ndarray, lo: float, hi: float, n: int = 150) -> Tuple[np.ndarray, np.ndarray]:
    """Resample one kr branch on a uniform Sw grid (the tables are densest near the endpoints)."""
    o = np.argsort(sw)
    x = np.linspace(lo, hi, n)
    return x, np.interp(x, sw[o], kr[o])


def fit_seed(dr: np.ndarray, im: np.ndarray) -> Dict[str, Dict[str, float]]:
    """Corey parameters for the four branches of a single run."""
    swi, sw_end = dr[:, 0].min(), im[:, 0].max()
    gas = dr[dr[:, 3] > 0]
    sw_bt = gas[:, 0].max()
    out = {}
    for key, tab, k, lo, hi, rising in (("krw_drainage", dr, 2, swi, 1.0, True), ("krg_drainage", gas, 3, swi, sw_bt, False),
                                         ("krw_imbibition", im, 2, swi, sw_end, True), ("krg_imbibition", im, 3, swi, sw_end, False)):
        x, y = grid(tab[:, 0], tab[:, k], lo, hi)
        n, ke, _ = corey_fit(x, y, lo, hi, rising)
        out[key] = dict(n=n, k=ke, lo=float(lo), hi=float(hi), rising=rising)
    return out


def fit_case(work: Path, df: pd.DataFrame, wett: str, pc: float) -> Dict[str, Dict[str, float]]:
    """Mean of the per-seed Corey parameters."""
    sel = df[(df.wett == wett) & (df.pc_max_kPa == pc) & (df.ift_mN_m.round(2) == 37.37)]
    fits = [fit_seed(*curves(work, c)) for c in sel.case]
    return {b: {f: (fits[0][b][f] if f == "rising" else float(np.mean([x[b][f] for x in fits]))) for f in fits[0][b]}
            for b in fits[0]}


def main() -> None:
    """Fit and draw each case on identical axes."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    work = Path(sys.argv[1])
    df = pd.read_csv(work / "summary.csv")
    OUT.mkdir(exist_ok=True)
    rows = []
    for i, (wett, pc) in enumerate(CASES, 1):
        p = fit_case(work, df, wett, pc)
        swi, sw_end, sw_bt = p["krw_drainage"]["lo"], p["krw_imbibition"]["hi"], p["krg_drainage"]["hi"]
        sgi = float(df[(df.wett == wett) & (df.pc_max_kPa == pc) & (df.ift_mN_m.round(2) == 37.37)].Sgi.mean())
        fig, ax = plt.subplots(figsize=(7.2, 5.2))
        style = {"krw_drainage": (BLUE, "-", "Water, drainage"), "krg_drainage": (RED, "-", "CO₂, drainage"),
                 "krw_imbibition": (BLUE, "--", "Water, imbibition"), "krg_imbibition": (RED, "--", "CO₂, imbibition")}
        for key, (col, ls, lab) in style.items():
            x, y = curve(p[key])
            ax.plot(x, y, color=col, ls=ls, lw=2.4, label=f"{lab}  (n = {p[key]['n']:.2f})")
        ax.axvline(swi, color="grey", lw=0.6, ls=":")
        ax.axvline(sw_end, color="grey", lw=0.6, ls=":")
        ax.text(swi + 0.01, 0.62, f"Swi = {swi:.2f}", fontsize=8, color="grey")
        ax.text(sw_end + 0.01, 0.62, f"1 − Sgr = {sw_end:.2f}", fontsize=8, color="grey")
        ax.set(xlabel="Water saturation, Sw", ylabel="Relative permeability", xlim=(0, 1), ylim=(0, 1),
               title=f"Berea CO₂–water · fitted Corey curves · {WETT_LABEL[wett]}\n"
                     f"drainage to {pc:g} kPa · Sgi {sgi:.2f} → Sgr {1 - sw_end:.2f}")
        ax.legend(frameon=False, fontsize=9, loc="upper center")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        name = f"kr_fit_{i}_{wett}_pc{str(pc).replace('.', 'p')}.png"
        fig.savefig(OUT / name, dpi=200)
        plt.close(fig)
        for k, v in p.items():
            rows.append(dict(case=name[:-4], wettability=wett, pc_max_kPa=pc, Sgi=sgi, branch=k, corey_exponent=v["n"],
                             kr_endpoint=v["k"], Sw_min=v["lo"], Sw_max=v["hi"], Sgr=1 - sw_end, Sgc_drainage=1 - sw_bt))
        logger.info("%s  Sgi %.2f Sgr %.2f Sgc %.2f  n: %s", name, sgi, 1 - sw_end, 1 - sw_bt,
                    {k: round(v["n"], 2) for k, v in p.items()})
    pd.DataFrame(rows).to_csv(OUT / "corey_parameters_all.csv", index=False, float_format="%.4f")


if __name__ == "__main__":
    main()
