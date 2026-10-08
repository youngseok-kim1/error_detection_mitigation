"""Figure: (a) <m_x> vs -ln(alpha) over all check subsets; (b) |bias| of m_x vs Trotter depth"""
import sys, json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
src, out = sys.argv[1], sys.argv[2]
recs = {r["tag"]: r for r in map(json.loads, open(src))}
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#4a3aa7"]; INK, MUTED, SURF, GRID = "#0b0b0b", "#898781", "#fcfcfb", "#e6e5e0"
fig, (a, b) = plt.subplots(1, 2, figsize=(11.2, 4.5), facecolor=SURF)
# ---- (a) subsets
r = recs["flag d=6 p=0.004"]; n = r["n"]; c = r["parity"]; k = n; I = r["ideal"][n]
xs = np.array([s[0] for s in c["subsets"]]); ys = np.array([s[1][n] for s in c["subsets"]]); mask = np.arange(len(xs))
par = (mask >> k) & 1; full = len(xs) - 1
a.axhline(I, color=INK, lw=1, ls=(0, (4, 3))); a.text(0.002, I + 0.004, "ideal", fontsize=9, color=INK, va="bottom")
rt = c["r_pool"][n]; est = c["est"]["ext, site-pooled training ratio"][n]; xstar = xs[full] * (1 + rt)
a.plot([0, xstar], [ys[0], est], color=MUTED, lw=1.2, zorder=1)
sel = (par == 0) & (mask != 0)
a.scatter(xs[sel], ys[sel], s=30, color=C[0], edgecolor=SURF, lw=.8, zorder=2, label="subsets of the 5 ancilla checks")
sel = (par == 1) & (mask != full)
a.scatter(xs[sel], ys[sel], s=34, marker="s", color=C[1], edgecolor=SURF, lw=.8, zorder=2, label="same subsets + parity check")
a.scatter([0, xs[full]], [ys[0], ys[full]], s=80, color=INK, edgecolor=SURF, lw=1.5, zorder=3, label="no post-selection / all 6 checks")
a.scatter([xstar], [est], s=95, marker="D", color=C[2], edgecolor=SURF, lw=1.5, zorder=4, label="extrapolated, ratio from Clifford training")
a.annotate("no post-selection", (0, ys[0]), xytext=(9, -3), textcoords="offset points", fontsize=9, color=INK, va="top")
a.annotate(f"all checks: bias {ys[full]-I:+.3f}", (xs[full], ys[full]), xytext=(10, -8), textcoords="offset points", fontsize=9, color=INK, va="top")
a.annotate(f"bias {est-I:+.4f}", (xstar, est), xytext=(10, -8), textcoords="offset points", fontsize=9, color=INK, va="top", ha="left")
a.set_xlabel(r"$-\ln\alpha$  (measured acceptance)", fontsize=10, color=INK); a.set_ylabel(r"mean magnetization $\langle m_x\rangle$", fontsize=10, color=INK)
a.set_title("Which check you post-select on matters\n5-qubit Ising ring, 6 Trotter steps, 120 CX, p = 0.004", fontsize=10.5, loc="left", color=INK)
a.set_xlim(-0.015, xstar * 1.32); a.set_ylim(ys[0] - 0.03, I + 0.03)
a.legend(frameon=False, fontsize=8.5, loc="lower right")
# ---- (b) bias vs depth
series = [("ED", "error detection only", "o"), ("ext, coverage ratio", "extrapolated, coverage ratio", "s"),
          ("ext, site-pooled training ratio", "extrapolated, ratio from Clifford training", "D"),
          ("PEC 1st order (exact model)", "1st-order PEC, exact noise model", "^"), ("ZNE raw, exp fit (G=1,2,3)", "ZNE, all noise amplified (no checks)", "v")]
ds = [2, 4, 6]
for (key, lab, mk), col in zip(series, [INK, C[1], C[2], C[0], C[4]]):
    y = [abs(recs[f"flag d={d} p=0.004"]["parity"]["est"][key][n] - recs[f"flag d={d} p=0.004"]["ideal"][n]) for d in ds]
    b.plot(ds, y, marker=mk, ms=7, lw=2, color=col, markeredgecolor=SURF, markeredgewidth=1.2, label=lab)
b.set_yscale("log"); b.set_xticks(ds); b.set_xlabel("Trotter steps", fontsize=10, color=INK); b.set_ylabel(r"|bias| of $\langle m_x\rangle$", fontsize=10, color=INK)
b.set_title("Bias vs depth, ancilla + parity checks\nsame ring, p = 0.004, infinite shots", fontsize=10.5, loc="left", color=INK)
b.legend(frameon=False, fontsize=8.5, loc="lower right")
for ax in (a, b):
    ax.set_facecolor(SURF); ax.grid(True, color=GRID, lw=.6, which="major"); ax.set_axisbelow(True)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=9, length=0, which="both")
fig.tight_layout(); fig.savefig(out, dpi=170, facecolor=SURF)
