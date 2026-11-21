"""Two tests of why the angle trend in Sgr differs between networks A and B.

(1) Volume reweighting: the quasi-static fill sequence does not depend on element volume, so giving network A the
    pore/throat volume split of network B (and vice versa) is an exact change to the Sw bookkeeping.
(2) Pore trapping efficiency vs pore coordination number.

Usage: python analyse_cause_tests.py <work_dir>
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from analyse_trapping_structure import HERE, components, gas_state
from build_viewer_cases import load_network, region_mask
from run_matrix import read_curve


def load(d: Path):
    """Network, masks, and gas volumes (pores then throats) at end of drainage and end of imbibition."""
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
    _, conn = components(net, pg1, tg1)
    trapped = np.r_[pg1, tg1] & ~conn
    return net, n, m, vol, best[1], (1 - sw0) * vol, (1 - sw1) * vol, trapped


def main() -> None:
    """Run both tests and print the tables."""
    work = Path(sys.argv[1])
    res, eff_rows = [], []
    share = {}
    for d in sorted(work.glob("[AB]_*_s*")):
        net_id, wett = d.name.split("_")[0], "_".join(d.name.split("_")[1:3])
        net, n, m, vol, mask, g0, g1, trapped = load(d)
        isp = np.r_[np.ones(n, bool), np.zeros(m, bool)]
        Vp, Vt = vol[mask & isp].sum(), vol[mask & ~isp].sum()
        share[net_id] = (Vp, Vt)
        G0p, G0t = g0[mask & isp].sum(), g0[mask & ~isp].sum()
        G1p, G1t = g1[mask & isp & trapped].sum(), g1[mask & ~isp & trapped].sum()
        res.append(dict(network=net_id, wett=wett, Vp=Vp, Vt=Vt, G0p=G0p, G0t=G0t, G1p=G1p, G1t=G1t))
        deg = np.zeros(n)
        for a, b in net["pairs"]:
            for e in (a, b):
                if e > 0:
                    deg[e - 1] += 1
        z = np.digitize(deg, [3.5, 5.5, 7.5])  # 1-3, 4-5, 6-7, 8+
        for c, name in enumerate(("z<=3", "z4-5", "z6-7", "z>=8")):
            sel = mask[:n] & (z == c) & (g0[:n] > 0)
            eff_rows.append(dict(network=net_id, wett=wett, zclass=name, n_drained_pores=sel.sum(),
                                 gas0=g0[:n][sel].sum() / (Vp + Vt), eff=g1[:n][sel & trapped[:n]].sum() / g0[:n][sel].sum() if g0[:n][sel].sum() > 0 else np.nan))
    df = pd.DataFrame(res).groupby(["network", "wett"]).mean(numeric_only=True).reset_index()
    pd.set_option("display.width", 250)
    ts = {n: share[n][1] / sum(share[n]) for n in share}
    print("throat volume share of network (calc box):", {k: round(v, 3) for k, v in ts.items()})
    rows = []
    for _, r in df.iterrows():
        for label, k in (("own split", 1.0), ("A split on B", None)):
            pass
    # reweight: factor k on throat volumes; target throat share s -> k = s/(1-s) * Vp/Vt
    out = []
    for _, r in df.iterrows():
        own = r.network
        other = "B" if own == "A" else "A"
        s_t = ts[other]
        k = s_t / (1 - s_t) * r.Vp / r.Vt
        tot_o, tot_k = r.Vp + r.Vt, r.Vp + k * r.Vt
        out.append(dict(network=own, wett=r.wett, Sgi_own=(r.G0p + r.G0t) / tot_o, Sgr_own=(r.G1p + r.G1t) / tot_o,
                        Sgr_with_other_throat_share=(r.G1p + k * r.G1t) / tot_k, Sgi_with_other_throat_share=(r.G0p + k * r.G0t) / tot_k,
                        Sgr_pores_only=r.G1p / r.Vp / 1.0 * 0 + (r.G1p / r.G0p) * 1.0))
    o = pd.DataFrame(out)
    o["Sgr_over_Sgi_own"] = o.Sgr_own / o.Sgi_own
    o["Sgr_over_Sgi_reweighted"] = o.Sgr_with_other_throat_share / o.Sgi_with_other_throat_share
    o = o.rename(columns={"Sgr_pores_only": "pore_eff"})
    order = ["WW_20-40", "MW_45-65", "IW_70-90", "NW_95-115"]
    o["wett"] = pd.Categorical(o.wett, order)
    o = o.sort_values(["network", "wett"])
    o.to_csv(HERE / "cause_test_reweighting.csv", index=False, float_format="%.4g")
    print("\nTEST 1: reweight volumes to the other network's pore/throat split")
    print(o[["network", "wett", "Sgi_own", "Sgr_own", "Sgi_with_other_throat_share", "Sgr_with_other_throat_share", "Sgr_over_Sgi_own", "Sgr_over_Sgi_reweighted"]].round(3).to_string(index=False))
    e = pd.DataFrame(eff_rows).groupby(["network", "wett", "zclass"]).mean(numeric_only=True).reset_index()
    e.to_csv(HERE / "cause_test_pore_efficiency_by_coordination.csv", index=False, float_format="%.4g")
    print("\nTEST 2: pore trapping efficiency by coordination number (drained pores); share of initial gas in class")
    print(e.pivot_table(index=["network", "zclass"], columns="wett", values="eff")[order].round(3).to_string())
    print(e[e.wett == "WW_20-40"].pivot_table(index="network", columns="zclass", values="gas0").round(3).to_string())


if __name__ == "__main__":
    main()
