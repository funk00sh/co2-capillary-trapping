"""Plot only the fitted (Corey) relative permeability curves for the base case.

Usage: python plot_fit_only.py <work_dir>
"""
import logging
import sys
from pathlib import Path
from typing import Dict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyse import BASE, BLUE, LOOP_PC, RED, WETT_LABEL, corey_fit, loop_data

logger = logging.getLogger(__name__)
HERE = Path(__file__).resolve().parent


def fit_branches(work: Path, df: pd.DataFrame, wett: str, pc: float) -> Dict[str, Dict[str, float]]:
    """Corey parameters (n, kr endpoint, Sw range) for the four kr branches."""
    _, _, md, mi = loop_data(work, df, wett, pc)
    swi, sw_end = md[0, 0], mi[-1, 0]
    gas = md[md[:, 2] > 0]
    sw_bt = gas[-1, 0]
    out = {}
    n, k, _ = corey_fit(md[:, 0], md[:, 1], swi, 1.0, True)
    out["krw_drainage"] = dict(n=n, k=k, lo=swi, hi=1.0, rising=True)
    n, k, _ = corey_fit(gas[:, 0], gas[:, 2], swi, sw_bt, False)
    out["krg_drainage"] = dict(n=n, k=k, lo=swi, hi=sw_bt, rising=False)
    n, k, _ = corey_fit(mi[:, 0], mi[:, 1], swi, sw_end, True)
    out["krw_imbibition"] = dict(n=n, k=k, lo=swi, hi=sw_end, rising=True)
    n, k, _ = corey_fit(mi[:, 0], mi[:, 2], swi, sw_end, False)
    out["krg_imbibition"] = dict(n=n, k=k, lo=swi, hi=sw_end, rising=False)
    return out


def curve(p: Dict[str, float], npts: int = 200) -> tuple:
    """Smooth Corey curve over its saturation range."""
    sw = np.linspace(p["lo"], p["hi"], npts)
    se = (sw - p["lo"]) / (p["hi"] - p["lo"])
    return sw, p["k"] * (se if p["rising"] else 1 - se) ** p["n"]


def main() -> None:
    """Draw the fit-only figure and write the parameter table."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    work = Path(sys.argv[1])
    df = pd.read_csv(work / "summary.csv")
    p = fit_branches(work, df, BASE, LOOP_PC)
    swi, sw_end, sw_bt = p["krw_drainage"]["lo"], p["krw_imbibition"]["hi"], p["krg_drainage"]["hi"]

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    style = {"krw_drainage": (BLUE, "-", "Water, drainage"), "krg_drainage": (RED, "-", "CO₂, drainage"),
             "krw_imbibition": (BLUE, "--", "Water, imbibition"), "krg_imbibition": (RED, "--", "CO₂, imbibition")}
    for key, (col, ls, lab) in style.items():
        x, y = curve(p[key])
        ax.plot(x, y, color=col, ls=ls, lw=2.4, label=f"{lab}  (n = {p[key]['n']:.2f})")
    ax.set(xlabel="Water saturation, Sw", ylabel="Relative permeability", xlim=(0, 1), ylim=(0, 1),
           title=f"Berea CO₂–water, fitted Corey curves · {WETT_LABEL[BASE]} · drainage to {LOOP_PC:g} kPa")
    ax.axvline(swi, color="grey", lw=0.6, ls=":")
    ax.text(swi + 0.01, 0.62, f"Swi = {swi:.2f}", fontsize=8, color="grey")
    ax.axvline(sw_end, color="grey", lw=0.6, ls=":")
    ax.text(sw_end + 0.01, 0.62, f"1 − Sgr = {sw_end:.2f}", fontsize=8, color="grey")
    ax.legend(frameon=False, fontsize=9, loc="upper center")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(HERE / "fig6_kr_fit_only.png", dpi=200)

    rows = [{"branch": k, "corey_exponent": v["n"], "kr_endpoint": v["k"], "Sw_min": v["lo"], "Sw_max": v["hi"]}
            for k, v in p.items()]
    pd.DataFrame(rows).to_csv(HERE / "corey_parameters.csv", index=False, float_format="%.4f")
    logger.info(pd.DataFrame(rows).round(3).to_string(index=False))
    logger.info("Swi %.3f  Sgc(drainage) %.3f  Sgr %.3f", swi, 1 - sw_bt, 1 - sw_end)


if __name__ == "__main__":
    main()
