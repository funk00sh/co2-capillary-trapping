"""Where does the trapped CO2 sit? Pore/throat split, element size and cluster size on networks A and B.

Usage: python analyse_trapping_structure.py <work_dir from run_structure_snapshots.py>
Uses the end-of-drainage and end-of-imbibition PNflow states (VTU) of each run.
"""
import re
import sys
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from build_viewer_cases import load_network, read_arrays, region_mask
from run_matrix import read_curve

HERE = Path(__file__).resolve().parent


def gas_state(path: Path, npore: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Gas flags of pores and throats and local Sw of every element from one VTU snapshot."""
    a = read_arrays(path, {"PointData": ["Sw", "ffaz"], "CellData": ["Sw", "ffaz"]})
    pg = a["PointData"]["ffaz"][2:2 + npore] == 2
    tg = a["CellData"]["ffaz"] == 2
    sw = np.r_[a["PointData"]["Sw"][2:2 + npore], a["CellData"]["Sw"]]
    return pg, tg, sw


def components(net: Dict[str, np.ndarray], pg: np.ndarray, tg: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Cluster label of every element (pores then throats) and a flag for clusters touching an open face."""
    n, m = len(pg), len(tg)
    c = net["pairs"]
    rows, cols = [], []
    face_el = np.zeros(n + m, bool)
    for t in range(m):
        if not tg[t]:
            continue
        for e in c[t]:
            if e > 0 and pg[e - 1]:
                rows.append(e - 1), cols.append(n + t)
        if c[t][0] <= 0 or c[t][1] <= 0:
            face_el[n + t] = True
    _, lab = connected_components(coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n + m, n + m)), directed=False)
    face_lab = np.unique(lab[face_el])
    return lab, np.isin(lab, face_lab)


def analyse(d: Path) -> Dict[str, float]:
    """Metrics for one run folder."""
    net = load_network(d)
    n, m = len(net["rp"]), len(net["rt"])
    vol = np.r_[net["vp"], net["vt"]]
    files = sorted((d / "run_res").glob("*.vtu"), key=lambda p: (int(p.name[0]), int(re.search(r"(\d+)\.vtu$", p.name)[1])))
    last_drain = [f for f in files if f.name[0] == "1"][-1]
    sw_native = read_curve(next(d.glob("*cycle2_*.csv")))[-1][0]
    pg0, tg0, sw0 = gas_state(last_drain, n)
    pg1, tg1, sw1 = gas_state(files[-1], n)
    best = None
    for flag in (True, False):
        pm, tm = region_mask(net, flag)
        mk = np.r_[pm, tm]
        s = float((sw1 * vol)[mk].sum() / vol[mk].sum())
        if best is None or abs(s - sw_native) < best[0]:
            best = (abs(s - sw_native), mk)
    mask = best[1]
    tot = vol[mask].sum()
    gas0, gas1 = (1 - sw0) * vol, (1 - sw1) * vol
    lab, conn = components(net, pg1, tg1)
    g1 = np.r_[pg1, tg1]
    trapped = g1 & ~conn
    isp = np.r_[np.ones(n, bool), np.zeros(m, bool)]
    out = dict(Sgi=gas0[mask].sum() / tot, Sgr=gas1[mask & trapped].sum() / tot, face_connected=gas1[mask & conn].sum() / tot,
               mask_mismatch=best[0])
    for tag, sel in (("p", isp), ("t", ~isp)):
        out[f"Sgi_{tag}"] = gas0[mask & sel].sum() / tot
        out[f"Sgr_{tag}"] = gas1[mask & trapped & sel].sum() / tot
        out[f"eff_{tag}"] = out[f"Sgr_{tag}"] / out[f"Sgi_{tag}"]
    # element size: radius of trapped pores / throats against all gas-filled pores / throats at end of drainage
    rp, rt = net["rp"] * 1e6, net["rt"] * 1e6
    out["r_pore_all_drained"] = rp[(gas0[:n] > 0) & mask[:n]].mean()
    out["r_pore_trapped"] = rp[trapped[:n] & mask[:n]].mean()
    out["r_thr_all_drained"] = rt[(gas0[n:] > 0) & mask[n:]].mean()
    out["r_thr_trapped"] = rt[trapped[n:] & mask[n:]].mean()
    # connectivity of trapped pores
    deg = np.zeros(n)
    for a, b in net["pairs"]:
        for e in (a, b):
            if e > 0:
                deg[e - 1] += 1
    out["coord_trapped_pores"] = deg[trapped[:n] & mask[:n]].mean()
    out["coord_all_pores"] = deg[mask[:n]].mean()
    # clusters of trapped gas
    tl = lab[trapped & mask]
    ids, inv = np.unique(tl, return_inverse=True)
    cvol = np.bincount(inv, weights=gas1[trapped & mask])
    cpores = np.bincount(inv, weights=isp[trapped & mask].astype(float))
    out["n_clusters"] = len(ids)
    out["mean_cluster_vol_um3"] = float(cvol.mean() * 1e18)
    out["vol_share_in_clusters_ge5_pores"] = float(cvol[cpores >= 5].sum() / cvol.sum())
    out["vol_share_in_throat_only_clusters"] = float(cvol[cpores == 0].sum() / cvol.sum())
    out["vol_share_in_single_pore_clusters"] = float(cvol[cpores == 1].sum() / cvol.sum())
    out["largest_cluster_share"] = float(cvol.max() / cvol.sum())
    return out


def main() -> None:
    """Analyse every run folder and write the tables."""
    work = Path(sys.argv[1])
    rows = []
    for d in sorted(work.glob("[AB]_*_s*")):
        net, wett, seed = d.name.split("_")[0], "_".join(d.name.split("_")[1:3]), int(d.name.split("_s")[-1])
        rows.append(dict(network=net, wett=wett, seed=seed, **analyse(d)))
        print(d.name, round(rows[-1]["Sgi"], 3), round(rows[-1]["Sgr"], 3), flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(HERE / "trapping_structure_runs.csv", index=False, float_format="%.5g")
    g = df.groupby(["network", "wett"]).mean(numeric_only=True).drop(columns="seed")
    g.to_csv(HERE / "trapping_structure_mean.csv", float_format="%.5g")
    pd.set_option("display.width", 250, "display.max_columns", 40)
    print(g.round(3).T.to_string())


if __name__ == "__main__":
    main()
