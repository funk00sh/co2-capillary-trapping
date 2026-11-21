"""Side-by-side structural statistics of networks A and B (same CT image, different extraction)."""
from pathlib import Path

import numpy as np

SIM = Path("../..")
NETS = {"A": Path("network_pnextract"), "B": HERE / "network_B"}
L = 2138e-6


def load(d: Path):
    """Arrays: pore volume/radius/G, throat radius/G/length/volume/pore-end lengths, pairs."""
    n2 = np.array([l.split()[:5] for l in (d / "Berea_node2.dat").read_text().splitlines() if l.strip()], float)
    l1 = np.array([l.split()[:6] for l in (d / "Berea_link1.dat").read_text().splitlines()[1:] if l.strip()], float)
    l2 = np.array([l.split()[:8] for l in (d / "Berea_link2.dat").read_text().splitlines() if l.strip()], float)
    return n2[:, 1], n2[:, 2], n2[:, 3], l1[:, 1:3].astype(int), l1[:, 3], l1[:, 4], l1[:, 5], l2[:, 6], l2[:, 3], l2[:, 4]


def q(x, p=(10, 50, 90), s=1e6):
    return " / ".join(f"{np.percentile(x, k) * s:.1f}" for k in p)


for name, d in NETS.items():
    vp, rp, Gp, pr, rt, Gt, lt, vt, lp1, lp2 = load(d)
    internal = (pr[:, 0] > 0) & (pr[:, 1] > 0)
    a, b = pr[internal, 0] - 1, pr[internal, 1] - 1
    z = np.bincount(np.r_[a, b], minlength=len(rp))
    rmax = np.maximum(rp[a], rp[b])
    rmin = np.minimum(rp[a], rp[b])
    asp = rmax / rt[internal]
    At = rt[internal] ** 2 / (4 * Gt[internal])
    Ap_end = rp[a] ** 2 / (4 * Gp[a])  # pore-end cross-section (approx, pore radius & G of pore a)
    Ap_end2 = rp[b] ** 2 / (4 * Gp[b])
    g_t = At ** 2 * Gt[internal] / lt[internal]
    g_p1 = Ap_end ** 2 * Gp[a] / np.maximum(lp1[internal], 1e-9)
    g_p2 = Ap_end2 ** 2 * Gp[b] / np.maximum(lp2[internal], 1e-9)
    res_t = 1 / g_t
    res_p = 1 / g_p1 + 1 / g_p2
    print(f"=== Network {name}: {len(rp)} pores, {internal.sum()} internal throats")
    print(f" porosity pores/throats/total: {vp.sum() / L**3:.3f} / {vt.sum() / L**3:.3f} / {(vp.sum() + vt.sum()) / L**3:.3f}; throat vol share {vt.sum() / (vp.sum() + vt.sum()):.3f}")
    print(f" pore radius um (p10/50/90): {q(rp)}   throat radius um: {q(rt[internal])}   throat length um: {q(lt[internal])}")
    print(f" pore-end length lp1+lp2 um: {q((lp1 + lp2)[internal])}   throat length / pore radius (median): {np.median(lt[internal] / rmax):.1f}")
    print(f" aspect ratio r_pore_max/r_throat (p10/50/90): {np.percentile(asp, 10):.2f} / {np.median(asp):.2f} / {np.percentile(asp, 90):.2f}")
    print(f" coordination mean {z.mean():.2f}; share of pores z=0 {np.mean(z == 0):.3f}, z<=2 {np.mean(z <= 2):.3f}, z>=6 {np.mean(z >= 6):.3f}, z>=8 {np.mean(z >= 8):.3f}; max {z.max()}")
    big = rp > np.percentile(rp, 75)
    print(f" top-quartile pores (r>{np.percentile(rp, 75) * 1e6:.1f} um): hold {vp[big].sum() / vp.sum():.2f} of pore volume, mean z {z[big].mean():.1f}")
    print(f" throat G median {np.median(Gt[internal]):.3f}  (share angular G<0.048: {np.mean(Gt[internal] < 0.048):.2f}); pore G median {np.median(Gp):.3f}")
    print(f" pore volume / sphere volume (median) {np.median(vp / (4 / 3 * np.pi * rp**3)):.2f}; throat vol / (area*length) median {np.median(vt[internal] / (At * lt[internal])):.2f}")
    print(f" hydraulic resistance share of throat (median over throats) {np.median(res_t / (res_t + res_p)):.2f}")
    print(f" throat radius / min adjacent pore radius (median) {np.median(rt[internal] / rmin):.2f}")
