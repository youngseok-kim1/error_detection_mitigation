import json, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from math import comb
out = sys.argv[1]
R = json.load(open("theory_results.json")); S = json.load(open("strat_results.json")); Nn = json.load(open("nested_results.json"))
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#4a3aa7"]; INK, MUTED, SURF, GRID = "#0b0b0b", "#898781", "#fcfcfb", "#e6e5e0"
def style(ax):
    ax.set_facecolor(SURF); ax.grid(True, color=GRID, lw=.6); ax.set_axisbelow(True)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=9, length=0, which="both")

# ---------- Figure 1: the post-selection family as a ZNE parameter
fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.6), facecolor=SURF)
rows = R["E1_M1000"]; Oall = rows[0][3]
a.plot([0, 1], [1, Oall], color=INK, lw=1.3, zorder=1, label=r"exact mixture law  $O_1 - u\,(O_1 - O_{all})$")
sub = np.array(R["E1_subsets"]); a.scatter(sub[:, 2], sub[:, 3], s=14, color=C[0], alpha=.6, edgecolor="none", zorder=2, label="twirled checks: all 64 check subsets")
a.scatter([r[2] for r in rows], [r[3] for r in rows], s=60, color=C[0], edgecolor=SURF, lw=1.2, zorder=4, label="twirled: first m checks (m = 0…6)")
thr = R["E1_thr_pm0.0"]; a.scatter([t[2] for t in thr], [t[3] for t in thr], s=60, marker="s", color=C[2], edgecolor=SURF, lw=1.2, zorder=4, label="twirled: accept ≤ t flags (t = 0…6)")
fx = np.array(R["E3_fixed_subsets"]); a.scatter(fx[:, 2], fx[:, 3], s=22, marker="D", color=C[1], alpha=.85, edgecolor="none", zorder=3, label="one fixed flavor: all 64 subsets")
a.scatter([0], [1], s=90, marker="*", color=INK, zorder=5); a.annotate("ideal (u = 0)", (0, 1), xytext=(14, 4), textcoords="offset points", fontsize=9, va="bottom")
a.annotate("no post-selection (u = 1)", (1, Oall), xytext=(-14, 2), textcoords="offset points", fontsize=9, ha="right", va="bottom")
a.set_xlabel(r"$u$ = (survival of a faulty shot) / (acceptance)   [measured]", fontsize=10, color=INK); a.set_ylabel(r"post-selected $\langle O\rangle$", fontsize=10, color=INK)
a.set_title("Post-selection as the extrapolation parameter\n10-qubit Clifford payload, 6 full-weight checks, p = 0.004", fontsize=10.5, loc="left", color=INK)
a.legend(frameon=False, fontsize=8.3, loc="lower left"); a.set_xlim(-0.03, 1.05); a.set_ylim(Oall - 0.012, 1.012); style(a)
# inset-like zoom: residual from the line
ins = a.inset_axes([0.54, 0.5, 0.44, 0.3]); ins.set_facecolor(SURF)
resid = lambda pts: np.array(pts)[:, 3] - (1 - np.array(pts)[:, 2] * (1 - Oall))
ins.axhline(0, color=INK, lw=.8); ins.scatter(sub[:, 2], resid(sub) * 1e3, s=10, color=C[0], alpha=.7, edgecolor="none"); ins.scatter(fx[:, 2], resid(fx) * 1e3, s=14, marker="D", color=C[1], alpha=.85, edgecolor="none")
ins.set_title("distance from the line, ×10⁻³ (twirled vs fixed)", fontsize=8, color=INK, loc="left"); ins.tick_params(labelsize=7, colors=MUTED, length=0); ins.grid(True, color=GRID, lw=.5)
for s_ in ins.spines.values(): s_.set_visible(False)
# ---------- (b) syndrome-weight histogram
h = R["E5"]; k = 6; w = np.arange(1, k + 1); binom = np.array([comb(k, i) / 2 ** k for i in range(k + 1)]); binom = binom[1:] / binom[1:].sum()
full = np.array(h["full"])[1:]; full /= full.sum(); loc = np.array(h["local"])[1:]; loc /= loc.sum()
b.bar(w - 0.22, binom, 0.2, color=INK, label="Binomial(6, ½), w ≥ 1"); b.bar(w, full, 0.2, color=C[0], label="twirled full-weight checks"); b.bar(w + 0.22, loc, 0.2, color=C[1], label="twirled weight-2 checks")
b.set_xlabel("number of flagged checks w, among flagged shots", fontsize=10, color=INK); b.set_ylabel("fraction", fontsize=10, color=INK)
b.set_title("Diagnostic: flag-count distribution tests uniform detection\nsame payload and noise", fontsize=10.5, loc="left", color=INK); b.legend(frameon=False, fontsize=8.5); style(b); b.grid(False, axis="x")
fig.tight_layout(); fig.savefig(f"{out}/fig_theory_mixture.png", dpi=170, facecolor=SURF)

# ---------- Figure 2: gadget noise and cost
fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.4), facecolor=SURF)
for pc, col in (("0.001", C[0]), ("0.004", C[1])):
    rows = R[f"E2_pc{pc}"]; ms = [r[0] for r in rows]
    a.plot(ms, [r[3] for r in rows], marker="o", ms=5, lw=1.2, ls=(0, (3, 2)), color=col, markeredgecolor=SURF, label=f"post-selection only, p_check = {pc}")
    a.plot(ms, [r[4] for r in rows], marker="s", ms=6, lw=2, color=col, markeredgecolor=SURF, label=f"single-species estimate, p_check = {pc}")
    mm = [2, 3, 4, 6]; a.plot(mm, [Nn[f"{pc}_{m}"]["nested"] for m in mm], marker="D", ms=7, lw=0, color=col, markeredgecolor=INK, markeredgewidth=.8, label=f"nested stratification, p_check = {pc}")
a.axhline(1, color=INK, lw=1, ls=(0, (4, 3))); a.text(6.05, 1.001, "ideal", fontsize=9, va="bottom", ha="right")
a.set_xlabel("number of twirled full-weight checks m", fontsize=10, color=INK); a.set_ylabel(r"$\langle O\rangle$ estimate", fontsize=10, color=INK); a.set_xticks([1, 2, 3, 4, 6])
a.set_title("Noisy check gadgets break uniform detection\npayload p = 0.004; gadgets add 15 two-qubit gates per check", fontsize=10.5, loc="left", color=INK); a.legend(frameon=False, fontsize=8, loc="lower right"); style(a)
# (b) cost
nv = R["E1_nvar"]; rows = R["E1_M1000"]; ms = range(1, 7)
plainED = [(1 - rows[m][3] ** 2) / rows[m][1] for m in ms]
b.plot(list(ms), nv, marker="o", ms=6, lw=2, color=C[0], markeredgecolor=SURF, label="leak-subtraction estimator (unbiased)")
b.plot(list(ms), plainED, marker="s", ms=6, lw=2, color=INK, markeredgecolor=SURF, label="plain post-selection (biased)")
b.set_yscale("log"); b.set_xlabel("number of twirled full-weight checks m", fontsize=10, color=INK); b.set_ylabel(r"$N\,\mathrm{Var}$ of the estimate", fontsize=10, color=INK)
b.set_title("Cost of the exact extrapolation\nnoiseless gadgets, p = 0.004; ideal value 1", fontsize=10.5, loc="left", color=INK); b.legend(frameon=False, fontsize=8.5); style(b)
fig.tight_layout(); fig.savefig(f"{out}/fig_theory_gadgets_cost.png", dpi=170, facecolor=SURF)
print("ok")
