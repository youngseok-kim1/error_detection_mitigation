"""Figure: rms |bias| of <X_i> for every estimator, eight scenarios (small multiples)"""
import sys, json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from collections import defaultdict
out = sys.argv[1]
recs = [json.loads(l) for fn in ("flagship.jsonl", "generic_lo.jsonl", "generic_hi.jsonl") for l in open(fn)]
G = defaultdict(list)
for r in recs: G[r["tag"]].append(r)
C = ["#2a78d6", "#eb6834", "#1baf7a"]; INK, MUTED, SURF, GRID = "#0b0b0b", "#898781", "#fcfcfb", "#e6e5e0"
rows = [("raw", "unmitigated", 0), ("ED", "error detection only", 0), ("ZNE raw, exp fit (G=1,2,3)", "standard ZNE, no checks", 0),
        ("ext, coverage ratio", "extrapolated, coverage ratio", 1), ("ext, site-pooled training ratio", "extrapolated, Clifford-trained ratio", 1),
        ("CDR on ED (pooled slope)", "Clifford regression on ED data", 1),
        ("PEC 1st order (exact model)", "1st-order PEC, exact model", 2), ("PEC 1st order (model 30% high)", "1st-order PEC, model 30% high", 2),
        ("u-amp ZNE, exp fit (G=1,2,3)", "ZNE, undetected faults amplified", 2)]
groups = ["baseline", "post-selection as the knob (no noise model)", "needs a noise model"]; marks = ["o", "D", "s"]
cols = [("flag d=6 p=0.004", "paper angles, 6 steps"), ("generic p=0.004", "random angles + rates, 4 steps")]
def rms(tag, cfg, key):
    n = G[tag][0]["n"]; B = np.array([np.array(r[cfg]["est"][key])[:n] - np.array(r["ideal"])[:n] for r in G[tag]]); return np.sqrt((B ** 2).mean())
fig, axs = plt.subplots(2, 4, figsize=(13.5, 6.6), sharey=True, sharex=True, facecolor=SURF)
y = np.arange(len(rows))[::-1]
for i, hi in enumerate((False, True)):
    for j, ((tag, lab), cfg) in enumerate([(c, cfg) for c in cols for cfg in ("anc", "parity")]):
        ax = axs[i, j]; t = tag.replace("0.004", "0.012") if hi else tag
        lam = np.mean([r[cfg]["Ld"] + r[cfg]["Lu"] for r in G[t]]); al = np.mean([r[cfg]["alpha"] for r in G[t]])
        v = [rms(t, cfg, k) for k, _, _ in rows]
        ax.axvline(v[1], color=MUTED, lw=1, ls=(0, (4, 3)), zorder=1)
        for g in range(3):
            sel = [m for m, (_, _, gg) in enumerate(rows) if gg == g]
            ax.scatter([v[m] for m in sel], y[sel], s=62, marker=marks[g], color=C[g] if g else INK, edgecolor=SURF, lw=1.2, zorder=3, label=groups[g])
        for m in range(len(rows)): ax.plot([1e-4, v[m]], [y[m], y[m]], color=GRID, lw=1, zorder=0)
        ax.set_xscale("log"); ax.set_xlim(4e-4, 0.6); ax.set_facecolor(SURF)
        ax.set_title(f"{lab}\n{'ancilla checks' if cfg == 'anc' else 'ancilla + parity'},  $\\Lambda$={lam:.2f},  $\\alpha$={al:.2f}", fontsize=9.3, loc="left", color=INK)
        ax.grid(True, axis="x", color=GRID, lw=.6, which="major"); ax.set_axisbelow(True)
        for s in ax.spines.values(): s.set_visible(False)
        ax.tick_params(colors=MUTED, labelsize=9, length=0, which="both")
        if i == 1: ax.set_xlabel(r"rms |bias| of $\langle X_i\rangle$", fontsize=10, color=INK)
    axs[i, 0].set_yticks(y); axs[i, 0].set_yticklabels([l for _, l, _ in rows], fontsize=9.3, color=INK)
from matplotlib.lines import Line2D
h, l = axs[0, 0].get_legend_handles_labels()
h.append(Line2D([0], [0], color=MUTED, lw=1, ls=(0, (4, 3)))); l.append("level of error detection only")
fig.legend(h, l, loc="lower center", ncol=4, frameon=False, fontsize=9.5, bbox_to_anchor=(0.55, -0.005))
fig.tight_layout(rect=(0, 0.045, 1, 1)); fig.savefig(out, dpi=160, facecolor=SURF)
