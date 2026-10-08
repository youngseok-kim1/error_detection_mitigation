import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
SURF, INK, GRID = "#fbfaf7", "#2a2a2a", "#e6e3dc"
C = {"ZNE exp (1,2)": "#2a78d6", "ZNE cum2 (1,2,3)": "#4a3aa7", "twirled + PEC(inv)": "#e34948", "repaired full twirl": "#0b0b0b"}
mk = {"ZNE exp (1,2)": "s", "ZNE cum2 (1,2,3)": "^", "twirled + PEC(inv)": "*", "repaired full twirl": "X"}
def style(ax):
    ax.set_facecolor(SURF); [ax.spines[s].set_visible(False) for s in ("top", "right")]; ax.grid(True, color=GRID, lw=.8); ax.tick_params(colors=INK, labelsize=9)
# 1% from the main sweep
base = json.load(open("demo_results_2_4_6_8_10_12_14_16.json")); rep = json.load(open("repaired_sweep.json"))["depths"]
def pt(bias, std, N=10000, eps=0.02): return N*std*std/(eps*eps-bias*bias) if abs(bias) < eps else np.inf
data = {0.01: {}}
for d, r in base["depths"].items():
    I = np.mean(r["ideal"][:5]); d = int(d); data[0.01][d] = {}
    for k in ("ZNE exp (1,2)", "ZNE cum2 (1,2,3)"): s = r["stage1"][k]; data[0.01][d][k] = (s["mx"]-I, s["mx_std"])
    s = r["stage3b"]["twirled stratified + PEC(inv)"]; data[0.01][d]["twirled + PEC(inv)"] = (s["mx"]-I, r["stage3"]["twirled stratified"]["mx_std"]*np.sqrt(s["cost_factor"]))
    s = rep[str(d)]["twirled stratified"]; data[0.01][d]["repaired full twirl"] = (s["mx"]-I, s["mx_std"]); data[0.01][d]["q0"] = rep[str(d)]["species"][0]
for p in (0.005, 0.002):
    J = json.load(open(f"threshold_p{p}.json")); data[p] = {}
    for d, r in J["depths"].items():
        d = int(d); I = r["ideal"]; data[p][d] = {k: (r[k]["mx"]-I, r[k]["mx_std"]) for k in C}; data[p][d]["q0"] = r["repaired full twirl"]["q0"]
fig, axs = plt.subplots(1, 4, figsize=(15, 4.2), facecolor=SURF)
for ax, p in zip(axs[:3], (0.01, 0.005, 0.002)):
    ds = sorted(data[p])
    for k in C:
        y = [pt(*data[p][d][k]) for d in ds]; yy = [v if np.isfinite(v) else np.nan for v in y]
        ax.plot(ds, yy, marker=mk[k], ms=9 if k == "twirled + PEC(inv)" else 7, lw=1.8, color=C[k], markeredgecolor=SURF, markeredgewidth=1, label=k)
        for d, v in zip(ds, y):
            if not np.isfinite(v): ax.scatter([d + (list(C).index(k)-1.5)*0.25], [2e7], marker="x", s=40, color=C[k])
    ax.set_yscale("log"); ax.set_ylim(8e2, 4e7); ax.set_xticks(ds); style(ax)
    ax.set_title(f"shots for RMSE ≤ 0.02,  p = q = {p*100:g}%", fontsize=10.5, loc="left", color=INK); ax.set_xlabel("Trotter steps", fontsize=10, color=INK)
axs[0].set_ylabel("shots (before post-selection)", fontsize=10, color=INK); axs[0].legend(frameon=False, fontsize=8, loc="upper left")
ax = axs[3]
for p, m, c in ((0.01, "o", "#898781"), (0.005, "s", "#eda100"), (0.002, "D", "#1baf7a")):
    ds = sorted(data[p]); q0 = [data[p][d]["q0"] for d in ds]; ratio = [data[p][d]["repaired full twirl"][1]/data[p][d]["ZNE exp (1,2)"][1] for d in ds]
    ax.plot(q0, ratio, marker=m, ms=7, lw=0, color=c, markeredgecolor=SURF, label=f"p = {p*100:g}%")
    for d, x, y in zip(ds, q0, ratio): ax.annotate(f"d={d}", (x, y), textcoords="offset points", xytext=(5, 3), fontsize=7, color=INK)
ax.axhline(1, color=INK, lw=.8, ls=(0, (4, 3))); ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.015, 0.8); ax.set_ylim(0.4, 300); style(ax)
import matplotlib.ticker as mt; ax.set_xticks([0.02,0.05,0.1,0.2,0.5]); ax.xaxis.set_major_formatter(mt.FormatStrFormatter("%g")); ax.xaxis.set_minor_formatter(mt.NullFormatter()); ax.set_yticks([0.5,1,2,5,10,100]); ax.yaxis.set_major_formatter(mt.FormatStrFormatter("%g")); ax.yaxis.set_minor_formatter(mt.NullFormatter())
ax.axvspan(0.25,0.4,color="#eda100",alpha=0.12,lw=0); ax.text(0.26,16,"break-even\n$q_0\\approx0.3$",fontsize=8,color=INK)
ax.set_xlabel("clean fraction $q_0$ of the repaired circuit", fontsize=10, color=INK); ax.set_ylabel("std(repaired twirl) / std(ZNE exp)", fontsize=10, color=INK)
ax.set_title("break-even against exponential ZNE", fontsize=10.5, loc="left", color=INK); ax.legend(frameon=False, fontsize=8, loc="upper right")
fig.tight_layout(); fig.savefig("fig_threshold.png", dpi=170, facecolor=SURF)
# table
rows = []
for p in (0.01, 0.005, 0.002):
    for d in sorted(data[p]):
        r = data[p][d]; rows.append(f"| {p*100:g}% | {d} | {r['q0']:.2f} | " + " | ".join(f"{r[k][0]:+.4f} / {r[k][1]:.3f} / " + ("∞" if not np.isfinite(pt(*r[k])) else f"{pt(*r[k])/1e3:.0f}k") for k in C) + " |")
open("threshold_table.md", "w").write("| p = q | d | $q_0$ (rep.) | ZNE exp: bias / std / shots | ZNE cum2 | twirled + PEC(inv) | repaired full twirl |\n|---|---|---|---|---|---|---|\n" + "\n".join(rows) + "\n")
print(open("threshold_table.md").read())
