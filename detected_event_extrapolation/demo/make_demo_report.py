"""Tables and figures from demo_results_*.json"""
import json, sys, glob, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
fn = sys.argv[1]; outdir = sys.argv[2]; rep_fn = sys.argv[3] if len(sys.argv) > 3 else None
R = json.load(open(fn)); D = {int(k): v for k, v in R["depths"].items()}; depths = sorted(D); N = R["N"]; n = R["n"]
REP = {int(k): v for k, v in json.load(open(rep_fn))["depths"].items()} if rep_fn else {}
C = {"raw": "#898781", "ZNE exp (1,2)": "#2a78d6", "ZNE cum2 (1,2,3)": "#4a3aa7", "fixed X PS": "#eda100", "ED+PEC opt": "#1baf7a", "twirled stratified": "#eb6834", "twirled + PEC(inv)": "#e34948", "repaired full twirl": "#0b0b0b"}
INK, MUTED, SURF, GRID = "#0b0b0b", "#898781", "#fcfcfb", "#e6e5e0"
def style(ax):
    ax.set_facecolor(SURF); ax.grid(True, color=GRID, lw=.6, which="major"); ax.set_axisbelow(True)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=9, length=0, which="both")

rows = {}   # protocol -> per depth dict(bias, std, nvar, rmse, Nreq02, Nreq03)
def add(name, d, bias, std, extra_factor=1.0):
    nv = N * std ** 2 * extra_factor; rm = np.sqrt(bias ** 2 + std ** 2 * extra_factor)
    def nreq(eps): return nv / (eps ** 2 - bias ** 2) if abs(bias) < eps else np.inf
    rows.setdefault(name, {})[d] = dict(bias=bias, std=std * np.sqrt(extra_factor), nvar=nv, rmse=rm, N02=nreq(0.02), N03=nreq(0.03), N05=nreq(0.05))
for d in depths:
    rec = D[d]; I = float(np.mean(rec["ideal"][:n]))
    for k in ("raw", "ZNE exp (1,2)", "ZNE cum2 (1,2,3)"): s = rec["stage1"][k]; add(k, d, s["mx"] - I, s["mx_std"])
    s = rec["stage2"]["fixed X PS"]; add("fixed X PS", d, s["mx"] - I, s["mx_std"])
    s = rec["stage2"]["fixed X ED+PEC optimistic"]; O = s["mx"]; g2a = s["cost"]; nv = g2a - O ** 2 / s["alpha"]   # (gamma^2 - O^2)/alpha with gamma^2/alpha = cost
    add("ED+PEC opt", d, O - I, np.sqrt(max(nv, 1e-9) / N))
    s = rec["stage2"]["fixed X ED+PEC payload-only"]; add("ED+PEC payload-only", d, s["mx"] - I, np.sqrt(max(s["cost"] - s["mx"] ** 2 / s["alpha"], 1e-9) / N))
    s = rec["stage3"]["twirled PS only (t=0)"]; add("twirled PS only", d, s["mx"] - I, s["mx_std"])
    s = rec["stage3"]["twirled stratified"]; add("twirled stratified", d, s["mx"] - I, s["mx_std"])
    sb = rec["stage3b"]["twirled stratified + PEC(inv)"]; add("twirled + PEC(inv)", d, sb["mx"] - I, s["mx_std"], extra_factor=sb["cost_factor"])
    if d in REP:
        sr = REP[d]["twirled stratified"]; add("repaired full twirl", d, sr["mx"] - I, sr["mx_std"])

# ---------- tables
L = ["## Results tables (mean magnetization $m_x$, $N = 10^4$ shots)\n"]
L.append("Ideal $m_x$ by depth: " + ", ".join(f"d={d}: {np.mean(D[d]['ideal'][:n]):.3f}" for d in depths) + "\n")
for key, title in (("bias", "Infinite-shot bias"), ("std", "Standard deviation at 10k shots"), ("rmse", "RMSE at 10k shots"), ("N02", "Shots needed for RMSE ≤ 0.02"), ("N03", "Shots needed for RMSE ≤ 0.03")):
    L.append(f"\n**{title}**\n\n| protocol | " + " | ".join(f"d={d}" for d in depths) + " |\n|---|" + "---|" * len(depths))
    for name in rows:
        vals = []
        for d in depths:
            v = rows[name][d][key]
            vals.append("∞" if np.isinf(v) else (f"{v:+.4f}" if key == "bias" else f"{v:.4f}" if key in ("std", "rmse") else f"{v/1e3:.0f}k" if v < 1e6 else f"{v/1e6:.1f}M"))
        L.append(f"| {name} | " + " | ".join(vals) + " |")
if REP: L.append("\n**Repaired full twirl** (any flavor, controlled-X repair; species f ≤ %d): " % json.load(open(rep_fn))["F"] + "; ".join(f"d={d}: gadget gates/step {REP[d]['gates_cp']:.1f}+{REP[d]['gates_cv']:.1f}, α(t=0)={REP[d]['alpha_t'][0]:.3f}, q₀={REP[d]['species'][0]:.3f}" for d in depths if d in REP))
L.append("\n**Acceptance and species** (twirled protocol): " + "; ".join(f"d={d}: α(t=0)={D[d]['stage3']['alpha_t'][0]:.3f}, clean fraction q₀={D[d]['stage3']['species'][0]:.3f}" for d in depths))
L.append("\n**Fixed-check acceptance**: " + ", ".join(f"d={d}: {D[d]['stage2']['fixed X PS']['alpha']:.3f}" for d in depths))
L.append("\n**ED+PEC optimistic reference cost** $e^{4\\Lambda_u}/\\alpha$: " + ", ".join(f"d={d}: {D[d]['stage2']['fixed X ED+PEC optimistic']['cost']:.1f}" for d in depths))
L.append("\n**Invisible-sector PEC cost factor** $e^{4\\Lambda_{inv}}$: " + ", ".join(f"d={d}: {D[d]['stage3b']['twirled stratified + PEC(inv)']['cost_factor']:.2f}" for d in depths))
open(f"{outdir}/demo_tables.md", "w").write("\n".join(L) + "\n")
json.dump({k: {str(d): v for d, v in vv.items()} for k, vv in rows.items()}, open(f"{outdir}/demo_summary.json", "w"), default=float)

# ---------- Figure 1: bias and RMSE vs depth
fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.5), facecolor=SURF)
show = ["raw", "ZNE exp (1,2)", "ZNE cum2 (1,2,3)", "fixed X PS", "ED+PEC opt", "twirled stratified", "twirled + PEC(inv)"] + (["repaired full twirl"] if REP else [])
mk = {"raw": "o", "ZNE exp (1,2)": "s", "ZNE cum2 (1,2,3)": "^", "fixed X PS": "v", "ED+PEC opt": "P", "twirled stratified": "D", "twirled + PEC(inv)": "*", "repaired full twirl": "X"}
for name in show:
    y = [abs(rows[name][d]["bias"]) for d in depths]; a.plot(depths, y, marker=mk[name], ms=7 if name != "twirled + PEC(inv)" else 10, lw=1.8, color=C[name], markeredgecolor=SURF, markeredgewidth=1, label=name)
    y = [rows[name][d]["rmse"] for d in depths]; b.plot(depths, y, marker=mk[name], ms=7 if name != "twirled + PEC(inv)" else 10, lw=1.8, color=C[name], markeredgecolor=SURF, markeredgewidth=1, label=name)
for ax, t in ((a, "Infinite-shot |bias| of $m_x$"), (b, f"RMSE of $m_x$ at {N//1000}k shots")):
    ax.set_yscale("log"); ax.set_xlabel("Trotter steps", fontsize=10, color=INK); ax.set_xticks(depths); ax.set_title(t + f"\n5-qubit TFIM ring, p = {R['p']}, readout {R['q']*100:.0f}%", fontsize=10.5, loc="left", color=INK); style(ax)
    ax.axhline(0.02, color=INK, lw=.8, ls=(0, (4, 3)))
a.text(depths[0], 0.022, "0.02", fontsize=8, color=INK); b.legend(frameon=False, fontsize=8, loc="lower right", ncol=2)
fig.tight_layout(); fig.savefig(f"{outdir}/fig_demo_bias_rmse.png", dpi=170, facecolor=SURF)

# ---------- Figure 2: shots needed for RMSE <= 0.02 / 0.03
fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.5), facecolor=SURF, sharey=True)
for ax, key, eps in ((a, "N02", 0.02), (b, "N03", 0.03)):
    for name in show:
        y = np.array([rows[name][d][key] for d in depths]); yy = np.where(np.isinf(y), np.nan, y)
        ax.plot(depths, yy, marker=mk[name], ms=7 if name != "twirled + PEC(inv)" else 10, lw=1.8, color=C[name], markeredgecolor=SURF, markeredgewidth=1, label=name)
        off = (show.index(name) - 3) * 0.12
        for d, v in zip(depths, y):
            if np.isinf(v): ax.scatter([d + off], [1.5e10], marker="x", s=40, color=C[name])
    ax.set_yscale("log"); ax.set_ylim(5e2, 3e10); ax.set_xlabel("Trotter steps", fontsize=10, color=INK); ax.set_xticks(depths)
    ax.set_title(f"Shots needed for RMSE ≤ {eps}\n(× at top: unreachable, bias ≥ {eps})", fontsize=10.5, loc="left", color=INK); style(ax)
a.set_ylabel("shots (before post-selection)", fontsize=10, color=INK); a.legend(frameon=False, fontsize=8.3, loc="upper left", bbox_to_anchor=(0.0, 0.91), ncol=2)
fig.tight_layout(); fig.savefig(f"{outdir}/fig_demo_shots.png", dpi=170, facecolor=SURF)
print("\n".join(L[:3]))
