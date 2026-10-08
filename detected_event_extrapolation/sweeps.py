import numpy as np, sys
from edzne import *

def gam_generic(info):
    """N x Var of the two-point estimator for a generic +-1 observable (per-shot variance -> 1)"""
    r, a = info["rO"], info["alpha"]; Dps, Dall = info["O_ps"], info["O_all"]
    return (1 + r) ** 2 / (a * Dps ** 2) + r ** 2 / Dall ** 2 - 2 * r * (1 + r) / (Dps * Dall)

def ensemble(label, n, k, reps, p, pm, wL, windows=None, seeds=range(40)):
    B, I = [], []
    for s in seeds:
        c, n2 = build(n, k, reps, p, pm, wL, seed=s, windows=windows)
        b, info, _, _ = analyse(c, n2, k, p, pm); info["G_ext_generic"] = gam_generic(info)
        B.append(b); I.append(info)
    print(f"\n== {label}: n={n} k={k} reps={reps} p={p}  ({len(B)} random instances)")
    for key in ["n2", "Lam", "cov", "alpha", "fd", "fu", "rO", "O_all", "O_ps", "G_pec", "G_edpec", "G_ext_generic"]:
        v = np.array([i[key] for i in I]); print(f"   {key:14s} median={np.median(v):9.4f}  [{v.min():.4f}, {v.max():.4f}]")
    for key in B[0]:
        v = np.array([b[key] for b in B]); print(f"   bias[{key:28s}] mean={v.mean():+.5f}  rms={np.sqrt((v**2).mean()):.5f}  max|.|={np.abs(v).max():.5f}")
    return B, I

n, k, reps, wL = 10, 6, 5, 2
M = 2 * reps
ensemble("A: checks wrap whole circuit, Lam~0.5", n, k, reps, 0.004, 0.004, wL)
ensemble("B: same, 3x noisier, Lam~1.4", n, k, reps, 0.012, 0.012, wL)
ensemble("C: only 2 checks (lower coverage)", n, 2, reps, 0.004, 0.004, wL)
ensemble("D: checks cover only 2nd half of circuit", n, k, reps, 0.004, 0.004, wL, windows=[(M // 2, M)] * k)
ensemble("E: checks cover only 1st half of circuit", n, k, reps, 0.004, 0.004, wL, windows=[(0, M // 2)] * k)

print("\n== scaling of bias with p (single instance, seed 3)")
for p in [0.001, 0.002, 0.004, 0.008, 0.016]:
    c, n2 = build(n, k, reps, p, p, wL, seed=3); b, info, _, _ = analyse(c, n2, k, p, p)
    print(f"  p={p:.3f} Lam={info['Lam']:.2f} alpha={info['alpha']:.3f} | ED {b['ED']:+.2e} | ext(O-spec) {b['ED+ext(O-specific r)']:+.2e} | +2nd {b['ED+ext(O-spec r)+2nd']:+.2e} | generic r {b['ED+ext(generic r)']:+.2e}")
