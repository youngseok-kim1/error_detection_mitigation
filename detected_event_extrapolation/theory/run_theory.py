import numpy as np, json, sys, time
from edtheory import *
n, reps, k = 10, 5, 6; M2 = 2 * reps; FULL = (1 << k) - 1
res = {}
def rep(name, E, extra=""):
    print(f"\n== {name}  gates: payload {E.gates[0]:.0f}, checks {E.gates[1]:.1f} {extra}")
    O_all = E.O_all(); print(f"   O_all = {O_all:.5f}")
    rows = []
    for m in range(k + 1):
        T = (1 << m) - 1; al, O = E.stats(subset_rule(k, T)); B = np.prod(E.b[:m]); u = B / al
        est = E.single_species(T); rows.append((m, al, u, O, est)); print(f"   m={m}  alpha={al:.4f}  u={u:.4f}  O_PS={O:.5f}  single-species O1={est:.5f}")
    return rows

t0 = time.time()
# ---- E1: full-weight twirled checks, noiseless gadgets and readout: exact mixture law
for pm in (0.0, 0.03):
    E = Ensemble(n, reps, k, [list(range(n))] * k, [(0, M2)] * k, p=0.004, p_check=0.0, pm=pm, M=300, seed=1)
    rows = rep(f"E1 full-weight twirled, p=0.004, p_check=0, pm={pm}", E); res[f"E1_pm{pm}"] = rows
    # threshold family: accept iff at most t flags
    thr = []
    for t in range(k + 1):
        al, O = E.stats(threshold_rule(k, t)); beta = sum(np.math.comb(k, j) for j in range(t + 1)) / 2 ** k if False else None
        # exact faulty survival for threshold t: Bin(k, 1-sigma) <= t (sigma ~ 1/2) ; fault-free: Bin(k, pm) <= t
        from math import comb
        bfault = sum(comb(k, j) * (1 - E.b[0]) ** j * E.b[0] ** (k - j) for j in range(t + 1)); afree = sum(comb(k, j) * pm ** j * (1 - pm) ** (k - j) for j in range(t + 1))
        est = ((al * O - bfault * E.O_all()) / (afree - bfault)) / ((al - bfault) / (afree - bfault))
        thr.append((t, al, bfault / al, O, est)); print(f"   t<={t}: alpha={al:.4f} u={bfault/al:.4f} O={O:.5f} est={est:.5f}")
    res[f"E1_thr_pm{pm}"] = thr
# ---- E1b: 3x noise
E = Ensemble(n, reps, k, [list(range(n))] * k, [(0, M2)] * k, p=0.012, p_check=0.0, pm=0.0, M=300, seed=1); res["E1_p012"] = rep("E1 full-weight twirled, p=0.012", E)
# ---- E2: gadget noise.  Circuits with exactly m twirled checks, all gadgets noisy; reference = same
#          circuit with noise only on the outermost right gadget (the part no check can see).
for pc in (0.001, 0.004):
    rows = []
    print(f"\n== E2 gadget noise p_check={pc}: m checks, all gadgets noisy (payload p=0.004)")
    for m in (1, 2, 3, 4, 6):
        T = (1 << m) - 1
        E = Ensemble(n, reps, m, [list(range(n))] * m, [(0, M2)] * m, p=0.004, p_check=pc, pm=0.0, M=300, seed=1)
        R = Ensemble(n, reps, m, [list(range(n))] * m, [(0, M2)] * m, p=0.0, p_check=0.0, pm=0.0, M=300, seed=1, gadget_p=lambda j, side, m=m, pc=pc: pc if (j == m - 1 and side == 'R') else 0.0)
        al, O = E.stats(subset_rule(m, T)); est = E.single_species(T); ref = R.O_all()
        rows.append((m, E.gates[1], al, O, est, ref)); print(f"   m={m}: check gates {E.gates[1]:.1f}  alpha={al:.4f}  O_all={E.O_all():.5f}  O_PS={O:.5f}  estimate={est:.5f}  reference (outer right gadget only)={ref:.5f}")
    res[f"E2_pc{pc}"] = rows
# ---- E3: fixed (untwirled) full-weight checks: subsets not on a line, single-species biased
E = Ensemble(n, reps, k, [list(range(n))] * k, [(0, M2)] * k, p=0.004, p_check=0.0, pm=0.0, M=1, seed=5, fixed=True); res["E3_fixed"] = rep("E3 fixed full-weight checks", E)
# all 64 subsets for the figure (fixed) and for the twirled case
def all_subsets(E):
    out = []
    for T in range(1 << k):
        al, O = E.stats(subset_rule(k, T)); B = np.prod([E.b[j] for j in range(k) if T >> j & 1]); out.append((T, al, B / al, O))
    return out
res["E3_fixed_subsets"] = all_subsets(E)
E = Ensemble(n, reps, k, [list(range(n))] * k, [(0, M2)] * k, p=0.004, p_check=0.0, pm=0.0, M=1000, seed=1); res["E1_subsets"] = all_subsets(E)
res["E1_M1000"] = rep("E1 full-weight twirled, p=0.004, M=1000 flavor draws", E)
print("   N x Var of the estimator vs m:", [round(E.nvar((1 << m) - 1), 3) for m in range(1, k + 1)], "  plain ED cost (1-O^2)/alpha at m=6:", round((1 - E.stats(subset_rule(k, FULL))[1] ** 2) / E.stats(subset_rule(k, FULL))[0], 3))
res["E1_nvar"] = [E.nvar((1 << m) - 1) for m in range(1, k + 1)]
# ---- E4: local twirled checks (weight-2 supports), k=3: single-species vs stratified
k3 = 3; rng = np.random.default_rng(11)
for label, sup in (("same support", [[0, 1]] * k3), ("distinct supports", [[0, 1], [4, 5], [7, 8]])):
    E = Ensemble(n, reps, k3, sup, [(0, M2)] * k3, p=0.004, p_check=0.0, pm=0.0, M=600, seed=2)
    print(f"\n== E4 local twirled, {label}: gates {E.gates}")
    out = {}
    for T in range(1 << k3):
        al, O = E.stats(subset_rule(k3, T)); est = E.single_species(T) if T else O; out[T] = (al, O, est)
        print(f"   T={T:03b} alpha={al:.4f} O_PS={O:.5f} single-species={est:.5f}")
    O1, q, W, cond = E.stratified(); print(f"   stratified O_1 = {O1:.5f}   species weights q_V = {np.round(q, 4)}   cond(M) = {cond:.1f}")
    res[f"E4_{label}"] = dict(subsets={str(T): v for T, v in out.items()}, strat=O1, q=q.tolist(), cond=cond)
# ---- E5: detection-count histogram among flagged shots: twirled full vs local
def flag_hist(E):
    k_ = E.k; Ps = E.Pm[:1 << k_] + E.Pm[1 << k_:]; h = np.zeros(k_ + 1)
    for s in range(1 << k_): h[bin(s).count("1")] += Ps[s]
    return h
E = Ensemble(n, reps, k, [list(range(n))] * k, [(0, M2)] * k, p=0.004, p_check=0.0, pm=0.0, M=300, seed=1); hf = flag_hist(E)
El = Ensemble(n, reps, k, [[(2 * j) % n, (2 * j + 1) % n] for j in range(k)], [(0, M2)] * k, p=0.004, p_check=0.0, pm=0.0, M=300, seed=3); hl = flag_hist(El)
res["E5"] = dict(full=hf.tolist(), local=hl.tolist())
from math import comb
print("\n== E5 syndrome-weight distribution among flagged shots (w = 1..k), vs Binomial(k,1/2) conditioned on w>=1")
binom = np.array([comb(k, w) / 2 ** k for w in range(k + 1)]); binom[0] = 0; binom /= binom.sum()
print("   twirled full :", np.round(hf[1:] / hf[1:].sum(), 4)); print("   twirled local:", np.round(hl[1:] / hl[1:].sum(), 4)); print("   Bin(6,1/2)   :", np.round(binom[1:], 4))
print(f"\n{time.time()-t0:.0f}s")
json.dump(res, open("theory_results.json", "w"))
