"""Figure: |bias| vs physical error rate for the Clifford toy (one instance), with p and p^2 guides"""
import sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from edzne import *
out = sys.argv[1]
C = ["#2a78d6", "#eb6834", "#1baf7a"]; INK, MUTED, SURF, GRID = "#0b0b0b", "#898781", "#fcfcfb", "#e6e5e0"
n, k, reps, wL = 10, 6, 5, 2
ps = np.geomspace(5e-4, 1.6e-2, 11); rows = {"ED": [], "ED+ext(O-specific r)": [], "ED+ext(O-spec r)+2nd": []}; lam = []
for p in ps:
    c, n2 = build(n, k, reps, p, p, wL, seed=3); b, info, _, _ = analyse(c, n2, k, p, p); lam.append(info["Lam"])
    for key in rows: rows[key].append(abs(b[key]))
fig, ax = plt.subplots(figsize=(6.6, 4.6), facecolor=SURF); ax.set_facecolor(SURF)
labels = {"ED": "error detection only", "ED+ext(O-specific r)": "extrapolated (O-specific ratio)", "ED+ext(O-spec r)+2nd": "+ pair term from rejected shots"}
for (key, y), col, mk in zip(rows.items(), [INK, C[0], C[1]], ["o", "s", "D"]):
    ax.loglog(ps, y, marker=mk, ms=6.5, lw=2, color=col, markeredgecolor=SURF, markeredgewidth=1.1, label=labels[key])
for e, y0, txt in ((1, rows["ED"][0], r"$\propto p$"), (2, rows["ED+ext(O-specific r)"][0], r"$\propto p^2$")):
    ax.loglog(ps, y0 * 0.55 * (ps / ps[0]) ** e, color=MUTED, lw=1, ls=(0, (4, 3)))
    ax.text(ps[-1] * 1.05, y0 * 0.55 * (ps[-1] / ps[0]) ** e, txt, fontsize=9, color=MUTED, va="center")
ax.set_xlabel("two-qubit error rate p", fontsize=10, color=INK); ax.set_ylabel(r"|bias| of $\langle O\rangle$  (ideal = 1)", fontsize=10, color=INK)
ax.set_title(f"Clifford toy: bias vs noise strength\n10 data qubits, 6 checks, total fault rate {lam[0]:.2f} to {lam[-1]:.2f}", fontsize=10.5, loc="left", color=INK)
ax.set_xlim(ps[0] * 0.85, ps[-1] * 1.6); ax.legend(frameon=False, fontsize=9, loc="lower right")
ax.grid(True, color=GRID, lw=.6, which="major"); ax.set_axisbelow(True)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(colors=MUTED, labelsize=9, length=0, which="both")
fig.tight_layout(); fig.savefig(out, dpi=170, facecolor=SURF)
for p, a, b_, c_ in zip(ps, *rows.values()): print(f"p={p:.5f}  ED={a:.2e}  ext={b_:.2e}  +pair={c_:.2e}")
