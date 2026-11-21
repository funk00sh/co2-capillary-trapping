"""One-factor-at-a-time variants of a PNflow network (numerical experiment, not a real rock).

Knobs
  coord : delete a random fraction of the internal throats (lower coordination number)
  radius: scale all throat radii (pore-to-throat aspect ratio); volumes are left unchanged to isolate the fill sequence
  shape : scale all throat shape factors G (angularity); volumes are left unchanged
Files keep the original token strings except the modified columns.
"""
from pathlib import Path
from typing import Dict, List

import numpy as np

G_MIN, G_MAX = 0.004, 0.0795  # circle = 1/(4 pi) = 0.0796


def _lines(path: Path) -> List[str]:
    return [l for l in path.read_text().splitlines() if l.strip()]


def build(src: Path, dst: Path, knob: str, level: float, seed: int = 0) -> Dict[str, float]:
    """Write the variant network into dst and return simple structure statistics."""
    dst.mkdir(parents=True, exist_ok=True)
    n1, n2 = _lines(src / "Berea_node1.dat"), _lines(src / "Berea_node2.dat")
    l1, l2 = _lines(src / "Berea_link1.dat"), _lines(src / "Berea_link2.dat")
    head1, body1 = n1[0], [l.split() for l in n1[1:]]
    nl = int(l1[0].split()[0])
    link1 = [l.split() for l in l1[1:]]
    link2 = [l.split() for l in l2]
    node2 = [l.split() for l in n2]
    assert len(link1) == nl and len(link2) == nl, (len(link1), len(link2), nl)
    rp = np.array([float(t[2]) for t in node2])
    internal = np.array([int(t[1]) > 0 and int(t[2]) > 0 for t in link1])

    if knob == "coord" and level > 0:
        rng = np.random.default_rng(seed)
        ids = np.flatnonzero(internal)
        drop = set(rng.choice(ids, size=int(round(level * len(ids))), replace=False).tolist())
        newid = {}
        for i in range(nl):
            if i not in drop:
                newid[i + 1] = len(newid) + 1
        link1 = [[str(newid[i + 1])] + t[1:] for i, t in enumerate(link1) if i not in drop]
        link2 = [[str(newid[i + 1])] + t[1:] for i, t in enumerate(link2) if i not in drop]
        new_body = []
        for t in body1:
            nc = int(t[4])
            assert len(t) == 5 + 2 * nc + 2, (t[:6], nc, len(t))
            nb, fl, lk = t[5:5 + nc], t[5 + nc:7 + nc], t[7 + nc:7 + 2 * nc]
            keep = [j for j in range(nc) if int(lk[j]) in newid]
            new_body.append(t[:4] + [str(len(keep))] + [nb[j] for j in keep] + fl + [str(newid[int(lk[j])]) for j in keep])
        body1 = new_body
        internal = np.array([int(t[1]) > 0 and int(t[2]) > 0 for t in link1])
    elif knob == "radius" and level != 1:
        for t in link1:
            p1, p2 = int(t[1]), int(t[2])
            r = float(t[3]) * level
            if p1 > 0 and p2 > 0:
                r = min(r, 0.99 * min(rp[p1 - 1], rp[p2 - 1]))
            t[3] = f"{r:.6E}"
    elif knob == "shape" and level != 1:
        for t in link1:
            t[4] = f"{min(max(float(t[4]) * level, G_MIN), G_MAX):.6E}"

    (dst / "Berea_node1.dat").write_text(head1 + "\n" + "\n".join("\t".join(t) for t in body1) + "\n")
    (dst / "Berea_node2.dat").write_text("\n".join(" ".join(t) for t in node2) + "\n")
    (dst / "Berea_link1.dat").write_text(f"{len(link1)}\n" + "\n".join("\t".join(t) for t in link1) + "\n")
    (dst / "Berea_link2.dat").write_text("\n".join("\t".join(t) for t in link2) + "\n")

    rt = np.array([float(t[3]) for t in link1])
    pr = np.array([[int(t[1]), int(t[2])] for t in link1])
    mx = np.array([max(rp[a - 1], rp[b - 1]) if a > 0 and b > 0 else np.nan for a, b in pr])
    vt = np.array([float(t[6]) for t in link2])
    vp = np.array([float(t[1]) for t in node2])
    # internal coordination number of pores
    z = np.zeros(len(rp))
    for a, b in pr[internal]:
        z[a - 1] += 1
        z[b - 1] += 1
    return dict(n_throats_internal=int(internal.sum()), coord_mean=float(2 * internal.sum() / len(rp)),
                throat_r_median_um=float(np.median(rt[internal]) * 1e6), pore_r_median_um=float(np.median(rp) * 1e6),
                aspect_median=float(np.median(mx[internal] / rt[internal])), G_median=float(np.median([float(t[4]) for t in link1])),
                throat_vol_share=float(vt.sum() / (vt.sum() + vp.sum())), frac_pores_z0=float((z == 0).mean()))


def build_combo(src: Path, dst: Path, f_radius: float = 1.0, g_scale: float = 1.0, share: float = None,
                c_len: float = 1.0) -> Dict[str, float]:
    """Combined variant: throat radius x f_radius, throat G x g_scale, throat volume share set to `share`
    (porosity conserved: volume moved between throats and pores), internal throat length x c_len (to hold K)."""
    dst.mkdir(parents=True, exist_ok=True)
    n1, n2 = _lines(src / "Berea_node1.dat"), _lines(src / "Berea_node2.dat")
    l1, l2 = _lines(src / "Berea_link1.dat"), _lines(src / "Berea_link2.dat")
    link1, link2, node2 = [l.split() for l in l1[1:]], [l.split() for l in l2], [l.split() for l in n2]
    rp = np.array([float(t[2]) for t in node2])
    vp = np.array([float(t[1]) for t in node2])
    vt = np.array([float(t[6]) for t in link2])
    internal = np.array([int(t[1]) > 0 and int(t[2]) > 0 for t in link1])
    for t, t2, inn in zip(link1, link2, internal):
        p1, p2 = int(t[1]), int(t[2])
        r = float(t[3]) * f_radius
        if inn:
            r = min(r, 0.99 * min(rp[p1 - 1], rp[p2 - 1]))
            t[5] = f"{float(t[5]) * c_len:.6E}"
            t2[5] = f"{float(t2[5]) * c_len:.6E}"
        t[3] = f"{r:.6E}"
        t[4] = f"{min(max(float(t[4]) * g_scale, G_MIN), G_MAX):.6E}"
    if share is not None:
        tot = vp.sum() + vt.sum()
        k, m = share * tot / vt.sum(), (1 - share) * tot / vp.sum()
        for t in link2:
            t[6] = f"{float(t[6]) * k:.6E}"
        for t in node2:
            t[1] = f"{float(t[1]) * m:.6E}"
    (dst / "Berea_node1.dat").write_text("\n".join(n1) + "\n")
    (dst / "Berea_node2.dat").write_text("\n".join(" ".join(t) for t in node2) + "\n")
    (dst / "Berea_link1.dat").write_text(f"{len(link1)}\n" + "\n".join("\t".join(t) for t in link1) + "\n")
    (dst / "Berea_link2.dat").write_text("\n".join("\t".join(t) for t in link2) + "\n")
    rt = np.array([float(t[3]) for t in link1])
    mx = np.array([max(rp[int(t[1]) - 1], rp[int(t[2]) - 1]) if inn else np.nan for t, inn in zip(link1, internal)])
    vt2 = np.array([float(t[6]) for t in link2])
    vp2 = np.array([float(t[1]) for t in node2])
    return dict(coord_mean=float(2 * internal.sum() / len(rp)), throat_r_median_um=float(np.median(rt[internal]) * 1e6),
                aspect_median=float(np.median(mx[internal] / rt[internal])),
                G_median=float(np.median([float(t[4]) for t in link1])), throat_vol_share=float(vt2.sum() / (vt2.sum() + vp2.sum())),
                length_scale=c_len)
