import sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from edzne import *
out = sys.argv[1]
BLUE, ORANGE, INK, MUTED, SURF = "#2a78d6", "#eb6834", "#0b0b0b", "#898781", "#fcfcfb"
n, k, reps, p, wL = 10, 6, 5, 0.004, 2; M = 2 * reps
cases = [("Checks wrap the whole circuit", None, 3), ("Checks cover only the second half", [(M // 2, M)] * k, 3)]
fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.3), sharey=True, facecolor=SURF)
for ax, (title, win, seed) in zip(axs, cases):
    c, n2 = build(n, k, reps, p, p, wL, seed=seed, windows=win)
    b, info, cls, P = analyse(c, n2, k, p, p)
    xs, ys = [], []
    for mask in range(1 << k):
        a, o = ps_stats(P, k, mask); xs.append(-np.log(a)); ys.append(np.log(o))
    xs, ys = np.array(xs), np.array(ys); xf, yf, y0 = xs[-1], ys[-1], ys[0]
    xstar = xf * (1 + info["rO"]); ystar = yf + info["rO"] * (yf - y0)
    ax.set_facecolor(SURF)
    ax.axhline(0, color=INK, lw=1, ls=(0, (4, 3)))
    ax.plot([0, xstar], [y0, ystar], color=MUTED, lw=1.2, zorder=1)
    ax.scatter(xs[1:-1], ys[1:-1], s=26, color=BLUE, alpha=.55, edgecolor=SURF, lw=.8, zorder=2, label="post-select on a subset of checks")
    ax.scatter([0, xf], [y0, yf], s=70, color=BLUE, edgecolor=SURF, lw=1.5, zorder=3, label="no post-selection / all checks")
    ax.scatter([xstar], [ystar], s=90, marker="D", color=ORANGE, edgecolor=SURF, lw=1.5, zorder=4, label="extrapolated to all faults removed")
    ax.annotate("no post-selection", (0, y0), xytext=(8, -4), textcoords="offset points", fontsize=9, color=INK, va="top")
    ax.annotate(f"all {k} checks: bias {b['ED']*100:+.1f}%", (xf, yf), xytext=(10, -6), textcoords="offset points", fontsize=9, color=INK, va="top", ha="left")
    ax.annotate(f"extrapolated: bias {b['ED+ext(O-specific r)']*100:+.2f}%", (xstar, ystar), xytext=(10, -6), textcoords="offset points", fontsize=9, color=INK, va="top", ha="left")
    ax.text(0.0, 0.006, "ideal", fontsize=9, color=INK, va="bottom")
    ax.set_title(f"{title}\ncoverage {info['cov']*100:.0f}%, acceptance {info['alpha']*100:.0f}%", fontsize=10.5, loc="left", color=INK)
    ax.set_xlabel(r"detected-fault weight removed,  $-\ln\alpha$  (measured)", fontsize=10, color=INK)
    ax.grid(True, color="#e6e5e0", lw=.6); ax.set_axisbelow(True)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=9, length=0)
    ax.set_xlim(-0.02, xstar * 1.55); ax.set_ylim(y0 - 0.035, 0.035)
axs[0].set_ylabel(r"$\ln\langle O\rangle$   (ideal = 0)", fontsize=10, color=INK)
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.01))
fig.tight_layout(rect=(0, 0.06, 1, 1)); fig.savefig(out, dpi=170, facecolor=SURF)
