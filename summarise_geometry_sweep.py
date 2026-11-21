"""Summarise the pore-geometry sweep: Sgr/Sgi per angle vs structure, and where the 95-115 vs 20-40 deg trend flips.

Usage: python summarise_geometry_sweep.py <work_dir>
"""
import json
import sys
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ORDER = ["WW_20-40", "MW_45-65", "IW_70-90", "NW_95-115"]
LAB = ["20–40°", "45–65°", "70–90°", "95–115°"]
XVAR = {"coord": ("coord_mean", "Mean coordination number"), "radius": ("aspect_median", "Median pore/throat radius ratio"),
        "shape": ("G_median", "Median throat shape factor G")}
TITLE = {"coord": "Coordination number (throats deleted)", "radius": "Aspect ratio (throat radii scaled)", "shape": "Throat angularity (G scaled)"}


def flip_point(x: np.ndarray, d: np.ndarray) -> Optional[float]:
    """x where d crosses zero (nearest crossing to the baseline is not needed: report all via first)."""
    s = np.where(np.diff(np.sign(d)) != 0)[0]
    if not len(s):
        return None
    i = s[0]
    return float(x[i] + (0 - d[i]) * (x[i + 1] - x[i]) / (d[i + 1] - d[i]))


def main() -> None:
    """Build tables and the figure."""
    work = Path(sys.argv[1])
    res = pd.DataFrame([json.loads(l) for l in (work / "results.jsonl").read_text().splitlines()])
    cal = pd.DataFrame(json.loads((work / "calibration.json").read_text()))
    st = pd.concat([cal.drop(columns="stats"), pd.DataFrame(list(cal.stats))], axis=1)
    res["ratio"] = res.Sgr / res.Sgi
    base = res[res.knob == "base"]
    tab = []
    for n in "AB":
        for knob in XVAR:
            sub = pd.concat([base[base.net == n].assign(knob=knob, level=0.0 if knob == "coord" else 1.0), res[(res.net == n) & (res.knob == knob)]])
            sts = pd.concat([st[(st.net == n) & (st.knob == "base")].assign(knob=knob, level=0.0 if knob == "coord" else 1.0), st[(st.net == n) & (st.knob == knob)]])
            for lv, g in sub.groupby("level"):
                s = sts[sts.level == lv].mean(numeric_only=True)
                row = dict(net=n, knob=knob, level=lv, Sgi=g.Sgi.mean(), x=s[XVAR[knob][0]], coord_mean=s.coord_mean, aspect=s.aspect_median,
                           G=s.G_median, vol_share=s.throat_vol_share, pc_kpa=g.pc_kpa.mean())
                for w in ORDER:
                    gw = g[g.wett == w]
                    row[f"Sgr_{w}"], row[f"ratio_{w}"], row[f"Sgr_sd_{w}"] = gw.Sgr.mean(), gw.ratio.mean(), gw.ratio.std()
                tab.append(row)
    t = pd.DataFrame(tab).sort_values(["net", "knob", "level"])
    t["delta_ratio_NW_minus_WW"] = t["ratio_NW_95-115"] - t["ratio_WW_20-40"]
    t.to_csv(HERE / "geometry_sweep_table.csv", index=False, float_format="%.4g")
    pd.set_option("display.width", 250)
    print(t[["net", "knob", "level", "x", "Sgi", "ratio_WW_20-40", "ratio_MW_45-65", "ratio_IW_70-90", "ratio_NW_95-115", "delta_ratio_NW_minus_WW"]].round(3).to_string(index=False))
    fig, axs = plt.subplots(2, 3, figsize=(14, 7.2), sharey=True)
    flips = {}
    for i, n in enumerate("AB"):
        for j, knob in enumerate(XVAR):
            ax = axs[i, j]
            g = t[(t.net == n) & (t.knob == knob)].sort_values("x")
            for w, lab in zip(ORDER, LAB):
                ax.errorbar(g.x, g[f"ratio_{w}"], yerr=g[f"Sgr_sd_{w}"], marker="o", ms=4, lw=1.8, capsize=2, label=lab)
            fp = flip_point(g.x.values, g.delta_ratio_NW_minus_WW.values)
            flips[(n, knob)] = fp
            b = g[g.level == (0.0 if knob == "coord" else 1.0)].x.iloc[0]
            ax.axvline(b, color="k", ls=":", lw=1)
            if fp is not None:
                ax.axvline(fp, color="r", ls="--", lw=1)
                ax.text(fp, 0.97, f" flip {fp:.2f}", color="r", fontsize=9, va="top")
            ax.set(ylim=(0, 1), title=f"Network {n}: {TITLE[knob]}" if i == 0 else f"Network {n}")
            ax.set_xlabel(XVAR[knob][1])
            ax.grid(alpha=0.3)
        axs[i, 0].set_ylabel("Sgr / Sgi")
    axs[0, 0].legend(title="advancing angle", frameon=False, fontsize=9, loc="lower left")
    fig.tight_layout()
    fig.savefig(HERE / "fig13_geometry_sweep_flip.png", dpi=190)
    print({f"{k[0]}_{k[1]}": (None if v is None else round(v, 3)) for k, v in flips.items()})


if __name__ == "__main__":
    main()
