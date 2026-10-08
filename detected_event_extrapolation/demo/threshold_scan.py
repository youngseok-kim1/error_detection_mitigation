"""Noise-level scan ("threshold" test): repaired full twirl vs constrained twirl (+PEC inv) vs ZNE at p = q in {1%, 0.5%, 0.2%}.
usage: python3 threshold_scan.py <p> <depths comma-separated> [M draws]"""
import numpy as np, json, time, sys
from tfim_demo import *; from protocols import *; from repaired import TFIMR, draw_full
n, TH, PH = 5, 3*np.pi/8, -np.pi/4
p = float(sys.argv[1]); q = p; depths = [int(x) for x in sys.argv[2].split(",")]; MD = int(sys.argv[3]) if len(sys.argv) > 3 else 200
N, R, F5, F7 = 10_000, 300, 5, 7
out = {"p": p, "q": q, "N": N, "R": R, "M": MD, "depths": {}}
mx = lambda v: float(np.mean(v[:n]))
def summ(inf, fin=None):
    d_ = dict(mx=mx(inf))
    if fin is not None: d_["mx_std"] = float(np.asarray(fin)[:, :n].mean(1).std())
    return d_
for d in depths:
    t0 = time.time(); rng = np.random.default_rng(1000 + d); I = mx(ideal_values(n, d, TH, PH)); rec = {"ideal": I}
    # ZNE
    dz = zne_dists(n, d, TH, PH, p); inf = zne_estimates({G: obs_from_dist(dz[G], n) for G in dz}); fin = zne_finite(dz, n, N, R, rng)
    for k in ("raw", "ZNE exp (1,2)", "ZNE cum2 (1,2,3)"): rec[k] = summ(inf[k], fin[k])
    # constrained twirl, stratified, + PEC(inv)
    pool = twirled_pool(n, d, TH, PH, p, q, MD, seed=d); pd = pooled(pool, d); inf3, al3, qv3 = twirled_estimates(pd, n, d, q, F5); fin3 = twirled_finite(pool, n, d, q, F5, N, R, rng)
    rec["twirled stratified"] = dict(summ(inf3["twirled stratified"], fin3["twirled stratified"]), alpha0=float(al3[0]), q0=float(qv3[0]))
    poolb = twirled_pool(n, d, TH, PH, p, q, MD, seed=d, inv_scale=0.0); infb, _, _ = twirled_estimates(pooled(poolb, d), n, d, q, F5)
    Lam_inv = 2*p/15*len(bonds(n))*d; cf = float(np.exp(4*Lam_inv))
    rec["twirled + PEC(inv)"] = dict(mx=mx(infb["twirled stratified"]), mx_std=rec["twirled stratified"]["mx_std"]*np.sqrt(cf), cost_factor=cf)
    # repaired full twirl
    M = TFIMR(n, d, TH, PH, p, q); poolr = []; rr = np.random.default_rng(d)
    for _ in range(MD):
        br = M.run_checked([draw_full(rr, n) for _ in range(d)]); poolr.append({c: xdist(v, n) for c, v in br.items()})
    pdr = pooled(poolr, d); infr, alr, qvr = twirled_estimates(pdr, n, d, q, F7); finr = twirled_finite(poolr, n, d, q, F7, N, R, rng)
    rec["repaired full twirl"] = dict(summ(infr["twirled stratified"], finr["twirled stratified"]), alpha0=float(alr[0]), q0=float(qvr[0]))
    out["depths"][d] = rec
    def line(k): r = rec[k]; b = r["mx"]-I; s = r["mx_std"]; nreq = N*s*s/(0.02**2-b*b) if abs(b) < 0.02 else float("inf"); return f"{k}: bias {b:+.4f} std {s:.4f} N(0.02) {nreq/1e3:.0f}k"
    print(f"p={p} d={d} ({time.time()-t0:.0f}s) ideal {I:.3f} | " + " | ".join(line(k) for k in ("ZNE exp (1,2)", "ZNE cum2 (1,2,3)", "twirled + PEC(inv)", "repaired full twirl")) + f" | q0 constr {rec['twirled stratified']['q0']:.3f} rep {rec['repaired full twirl']['q0']:.3f}", flush=True)
    json.dump(out, open(f"threshold_p{p}.json", "w"))
