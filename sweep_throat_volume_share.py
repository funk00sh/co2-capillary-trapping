"""Continuous sweep of the throat volume share using the saved snapshot runs (exact: fill sequence is volume-independent).

Usage: python sweep_throat_volume_share.py <work_dir from run_structure_snapshots.py>
Prints/plots Sgr (and Sgr/Sgi) at each angle vs throat volume share, and the share where the 95-115 vs 20-40 trend flips.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyse_cause_tests import load
from analyse_trapping_structure import HERE

ORDER = ["WW_20-40", "MW_45-65", "IW_70-90", "NW_95-115"]
LAB = ["20–40°", "45–65°", "70–90°", "95–115°"]


def main() -> None:
    """Sweep and plot."""
    work = Path(sys.argv[1])
    rows = []
    for d in sorted(work.glob("[AB]_*_s*")):
        net_id, wett = d.name.split("_")[0], "_".join(d.name.split("_")[1:3])
        net, n, m, vol, mask, g0, g1, trapped = load(d)
        isp = np.r_[np.ones(n, bool), np.zeros(m, bool)]
        rows.append(dict(network=net_id, wett=wett, Vp=vol[mask & isp].sum(), Vt=vol[mask & ~isp].sum(),
                         G0p=g0[mask & isp].sum(), G0t=g0[mask & ~isp].sum(),
                         G1p=g1[mask & isp & trapped].sum(), G1t=g1[mask & ~isp & trapped].sum()))
    df = pd.DataFrame(rows).groupby(["network", "wett"]).mean(numeric_only=True).reset_index()
    shares = np.linspace(0.02, 0.9, 89)
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.3), sharey=True)
    res = {}
    for ax, n in zip(axs, "AB"):
        cur = {}
        for w, lab in zip(ORDER, LAB):
            r = df[(df.network == n) & (df.wett == w)].iloc[0]
            k = shares / (1 - shares) * r.Vp / r.Vt
            tot = r.Vp + k * r.Vt
            sgr = (r.G1p + k * r.G1t) / tot
            sgi = (r.G0p + k * r.G0t) / tot
            cur[w] = sgr / sgi
            ax.plot(shares * 100, sgr / sgi, lw=2, label=lab)
        d = cur["NW_95-115"] - cur["WW_20-40"]
        sign = np.where(np.diff(np.sign(d)) != 0)[0]
        flip = float(np.interp(0, [d[sign[0]], d[sign[0] + 1]], [shares[sign[0]], shares[sign[0] + 1]]) * 100) if len(sign) else None
        res[n] = (flip, float(df[df.network == n].iloc[0].Vt / (df[df.network == n].iloc[0].Vp + df[df.network == n].iloc[0].Vt) * 100))
        ax.axvline(res[n][1], color="k", ls=":", lw=1)
        ax.text(res[n][1] + 1, 0.05, f"actual {res[n][1]:.0f}%", fontsize=9)
        if flip:
            ax.axvline(flip, color="r", ls="--", lw=1)
            ax.text(flip + 1, 0.9, f"flip at {flip:.0f}%", color="r", fontsize=9)
        ax.set(title=f"Network {n} fill sequence", xlabel="Throat share of pore volume (%)", ylim=(0, 1))
        ax.grid(alpha=0.3)
    axs[0].set_ylabel("Sgr / Sgi (fraction of CO₂ trapped)")
    axs[0].legend(title="advancing angle", frameon=False, loc="lower left", bbox_to_anchor=(0.0, 0.1))
    fig.tight_layout()
    fig.savefig(HERE / "fig12_throat_volume_share_sweep.png", dpi=200)
    print(res)


if __name__ == "__main__":
    main()
