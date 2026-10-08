"""ZNE from below: per-fault heralding scales every rate by G<1 (acceptance alpha(G) = prod(1-(1-G)p)).
Compare extrapolation to G=0 from points above the native noise (standard ZNE) and from points below it,
with the same ansatze, on the paper-angle Ising ring (exact density matrix)."""
import numpy as np, json, time
from ising_ed import *
from study import obs_matrix
TH, PH = 3 * np.pi / 8, -np.pi / 4; n = 5; out = {}
def fit_exp(Gs, ys):                     # ln y linear in G -> intercept
    return np.sign(ys[0]) * np.exp(np.polyfit(Gs, np.log(np.abs(ys)), 1)[1])
def fit_cum2(Gs, ys):                    # ln y quadratic in G (2nd cumulant)
    return np.sign(ys[0]) * np.exp(np.polyfit(Gs, np.log(np.abs(ys)), 2)[2])
def fit_lin(Gs, ys): return np.polyfit(Gs, ys, 1)[1]
for d in (4, 6):
    for p in (0.004, 0.012):
        M = Ising(n, d, TH, PH); W, names = obs_matrix(n); O0 = W @ ideal(n, d, TH, PH)
        Gs = [0.125, 0.25, 0.5, 1.0, 2.0, 3.0]; C = {}
        t = time.time()
        for G in Gs: C[G] = W @ post(M.run(p, p, G, G), n, False)[2]          # unpostselected value at scaled rates
        nCX = len(M.cx); alpha = {G: ((1 - (1 - G) * p) ** (nCX + M.k) if G < 1 else 1.0) for G in Gs}   # heralding acceptance below native noise; amplification above costs nothing per shot
        sel = slice(0, n)                                                     # single-site <X_i>
        res = {}
        for label, pts, f in [("exp above (1,2)", [1, 2], fit_exp), ("exp below (1,1/2)", [1, .5], fit_exp), ("exp below (1,1/4)", [1, .25], fit_exp),
                              ("exp bracket (1/2,2)", [.5, 2], fit_exp), ("cum2 above (1,2,3)", [1, 2, 3], fit_cum2), ("cum2 below (1,1/2,1/4)", [1, .5, .25], fit_cum2),
                              ("cum2 bracket (1/4,1,3)", [.25, 1, 3], fit_cum2), ("linear below (1,1/8)", [1, .125], fit_lin)]:
            est = np.array([f(np.array(pts), np.array([C[G][i] for G in pts])) for i in range(2 * n + 2)])
            # cost proxy: Lagrange weights for intercept of a degree-(len-1) polynomial through the points, N Var = sum h^2 (1-C^2)/alpha
            P_ = np.array(pts); h = [np.prod([-P_[j] / (P_[i] - P_[j]) for j in range(len(P_)) if j != i]) for i in range(len(P_))]
            cost = sum(h[i] ** 2 * (1 - C[P_[i]][n // 2] ** 2) / alpha[P_[i]] for i in range(len(P_)))
            res[label] = dict(rms=float(np.sqrt(((est[sel] - O0[sel]) ** 2).mean())), cost=float(cost))
        res["raw"] = dict(rms=float(np.sqrt(((C[1.0][sel] - O0[sel]) ** 2).mean())), cost=float(1 - C[1.0][n // 2] ** 2))
        out[f"d={d} p={p}"] = dict(res=res, alpha={str(G): alpha[G] for G in Gs}, C_mid={str(G): float(C[G][n // 2]) for G in Gs}, ideal=float(O0[n // 2]))
        print(f"\n== d={d} p={p}  ({time.time()-t:.0f}s)  ideal <X_2>={O0[n//2]:.4f}   C(G) at mid site:", {G: round(float(C[G][n // 2]), 4) for G in Gs})
        print("   per-fault acceptance alpha(G):", {G: round(alpha[G], 3) for G in Gs})
        for k, v in res.items(): print(f"   {k:26s} rms bias {v['rms']:.5f}   N Var {v['cost']:.2f}")
json.dump(out, open("zne_below.json", "w"))
