"""Trapping efficiency by element size class (radius terciles within each network), pores and throats separately.

Usage: python analyse_trapping_by_size.py <work_dir>
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from analyse_trapping_structure import HERE, components, gas_state
from build_viewer_cases import load_network, region_mask
from run_matrix import read_curve


def run_rows(d: Path) -> list:
    """Per-size-class initial gas, trapped gas and efficiency for one run."""
    net = load_network(d)
    n, m = len(net["rp"]), len(net["rt"])
    vol = np.r_[net["vp"], net["vt"]]
    files = sorted((d / "run_res").glob("*.vtu"), key=lambda p: (int(p.name[0]), int(re.search(r"(\d+)\.vtu$", p.name)[1])))
    pg0, tg0, sw0 = gas_state([f for f in files if f.name[0] == "1"][-1], n)
    pg1, tg1, sw1 = gas_state(files[-1], n)
    sw_native = read_curve(next(d.glob("*cycle2_*.csv")))[-1][0]
    best = None
    for flag in (True, False):
        pm, tm = region_mask(net, flag)
        mk = np.r_[pm, tm]
        s = float((sw1 * vol)[mk].sum() / vol[mk].sum())
        if best is None or abs(s - sw_native) < best[0]:
            best = (abs(s - sw_native), mk)
    mask = best[1]
    tot = vol[mask].sum()
    lab, conn = components(net, pg1, tg1)
    trapped = np.r_[pg1, tg1] & ~conn
    gas0, gas1 = (1 - sw0) * vol, (1 - sw1) * vol
    rows = []
    for kind, sl, r in (("pore", slice(0, n), net["rp"]), ("throat", slice(n, n + m), net["rt"])):
        mk = mask[sl]
        q = np.quantile(r[mk], [1 / 3, 2 / 3])
        cls = np.digitize(r, q)
        for c, name in enumerate(("small", "medium", "large")):
            sel = mk & (cls == c)
            g0, g1 = gas0[sl][sel].sum() / tot, gas1[sl][sel & trapped[sl]].sum() / tot
            rows.append(dict(kind=kind, size=name, vol_share_of_network=vol[sl][sel].sum() / tot, Sgi_class=g0, Sgr_class=g1,
                             eff=g1 / g0 if g0 > 0 else np.nan))
    return rows


def main() -> None:
    """Aggregate over seeds and print."""
    work = Path(sys.argv[1])
    out = []
    for d in sorted(work.glob("[AB]_*_s*")):
        net, wett = d.name.split("_")[0], "_".join(d.name.split("_")[1:3])
        out += [dict(network=net, wett=wett, **r) for r in run_rows(d)]
    df = pd.DataFrame(out)
    g = df.groupby(["network", "kind", "size", "wett"]).mean(numeric_only=True).reset_index()
    g.to_csv(HERE / "trapping_by_size_class.csv", index=False, float_format="%.4g")
    pd.set_option("display.width", 250)
    for k in ("pore", "throat"):
        print(k.upper())
        print(g[g.kind == k].pivot_table(index=["network", "size"], columns="wett", values=["Sgr_class", "eff"]).round(3).to_string())
        print(g[(g.kind == k) & (g.wett == "WW_20-40")][["network", "size", "vol_share_of_network", "Sgi_class"]].round(3).to_string(index=False))
        print()


if __name__ == "__main__":
    main()
