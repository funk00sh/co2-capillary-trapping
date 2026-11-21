"""Summary figure: where the initial CO2 sits and how each location traps it vs advancing angle (networks A and B)."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyse_trapping_structure import HERE

ORDER = ["WW_20-40", "MW_45-65", "IW_70-90", "NW_95-115"]
LAB = ["20–40°", "45–65°", "70–90°", "95–115°"]
m = pd.read_csv(HERE / "trapping_structure_mean.csv").set_index(["network", "wett"])
z = pd.read_csv(HERE / "cause_test_pore_efficiency_by_coordination.csv")

plt.rcParams.update({"font.size": 10.5, "axes.grid": True, "grid.alpha": 0.3, "axes.titlesize": 11.5, "axes.titleweight": "bold"})
fig, ax = plt.subplots(1, 3, figsize=(14.2, 4.5), gridspec_kw={"width_ratios": [0.8, 1, 1]})

# (a) where the initial CO2 sits (WW run, drainage is identical for all angles)
parts = {}
for n in "AB":
    zz = z[(z.network == n) & (z.wett == "WW_20-40")].set_index("zclass").gas0
    sgi = m.loc[(n, "WW_20-40"), "Sgi"]
    parts[n] = [m.loc[(n, "WW_20-40"), "Sgi_t"] / sgi, zz[["z<=3", "z4-5", "z6-7"]].sum() / sgi, zz["z>=8"] / sgi]
cols = ["#7f7f7f", "#6baed6", "#08519c"]
bottom = np.zeros(2)
for i, (name, c) in enumerate(zip(["throats", "pores, ≤ 7 throats", "pores, ≥ 8 throats"], cols)):
    vals = np.array([parts["A"][i], parts["B"][i]])
    ax[0].bar(["Network A", "Network B"], vals * 100, bottom=bottom * 100, color=c, label=name, width=0.6)
    for k, v in enumerate(vals):
        ax[0].text(k, (bottom[k] + v / 2) * 100, f"{v * 100:.0f}%", ha="center", va="center", color="white", fontweight="bold")
    bottom += vals
ax[0].set(ylabel="Share of initial CO₂ (%)", title="(a) Where the CO₂ sits at Sgi ≈ 0.75", ylim=(0, 100))
ax[0].legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), fontsize=9)
ax[0].grid(False)

# (b) trapping efficiency by location
def series(n, kind):
    if kind == "t":
        return [m.loc[(n, w), "eff_t"] for w in ORDER]
    zc = z[(z.network == n)].pivot_table(index="wett", columns="zclass", values="eff")
    return [zc.loc[w, kind] for w in ORDER]

for kind, col, name in (("t", "#7f7f7f", "throats"), ("z<=3", "#6baed6", "pores ≤ 3 throats"), ("z>=8", "#08519c", "pores ≥ 8 throats")):
    ax[1].plot(LAB, series("A", kind), "-o", color=col, lw=2, label=f"A · {name}")
    ax[1].plot(LAB, series("B", kind), "--s", color=col, lw=2, mfc="white", label=f"B · {name}")
ax[1].set(ylabel="Fraction of that CO₂ trapped", xlabel="Advancing contact angle", ylim=(0, 1), title="(b) Trapping efficiency by location")
ax[1].legend(frameon=False, fontsize=8.5, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.17))

# (c) net result
for n, ls, mk in (("A", "-", "o"), ("B", "--", "s")):
    ax[2].plot(LAB, [m.loc[(n, w), "Sgr"] for w in ORDER], ls, marker=mk, lw=2.2, color="#d62728" if n == "A" else "#1f77b4",
               mfc="white" if n == "B" else None, label=f"Network {n}")
ax[2].set(ylabel="Residual CO₂ saturation, Sgr", xlabel="Advancing contact angle", ylim=(0, 0.7), title="(c) Net result: Sgr at Sgi ≈ 0.75")
ax[2].legend(frameon=False, loc="upper center")
fig.tight_layout()
fig.savefig(HERE / "fig11_cause_trapping_A_vs_B.png", dpi=200)
