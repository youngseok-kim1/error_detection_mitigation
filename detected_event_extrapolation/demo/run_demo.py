import numpy as np, json, time, sys
from tfim_demo import *; from protocols import *
n, TH, PH, p, q = 5, 3 * np.pi / 8, -np.pi / 4, 0.01, 0.01
depths = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [2, 4, 6, 8, 10, 12]
N, R, MDRAWS, F = 10_000, 300, 200, 5
out = {"n": n, "p": p, "q": q, "N": N, "R": R, "M": MDRAWS, "F": F, "depths": {}}
def summ(inf, fin=None):
    """inf: vector of 2n observables; fin: array (R, 2n) of finite-shot estimates -> site-averaged summaries"""
    d_ = dict(inf=np.asarray(inf).tolist(), mx=float(np.mean(inf[:n])), cxx=float(np.mean(inf[n:])))
    if fin is not None:
        fin = np.asarray(fin); mxs = fin[:, :n].mean(1); cxs = fin[:, n:].mean(1)
        d_.update(mean=fin.mean(0).tolist(), std=fin.std(0).tolist(), mx_mean=float(mxs.mean()), mx_std=float(mxs.std()), cxx_mean=float(cxs.mean()), cxx_std=float(cxs.std()))
    return d_
rng = np.random.default_rng(0)
for d in depths:
    t0 = time.time(); ideal = ideal_values(n, d, TH, PH); rec = {"ideal": ideal.tolist()}
    # ---- stage 1
    dz = zne_dists(n, d, TH, PH, p); inf = zne_estimates({G: obs_from_dist(dz[G], n) for G in dz}); fin = zne_finite(dz, n, N, R, rng)
    rec["stage1"] = {k: summ(inf[k], fin[k]) for k in inf}
    # ---- stage 2
    rec["stage2"] = {}
    for kind in ("X", "Z"):
        dd = fixed_dists(n, d, TH, PH, p, q, kind); al, Ops = ps_estimate(dd, n); fe, fa = fixed_finite(dd, n, N, R, rng)
        rec["stage2"][f"fixed {kind} PS"] = dict(summ(Ops, fe), alpha=float(al), alpha_finite=float(fa.mean()))
        # ED + first-order PEC reference: undetected payload faults set to zero (expectation of PEC), cost e^{4 Lam_u}/alpha
        PL = fixed_flavor(n, kind); und1, und2 = classify_fixed(PL, 0, 1)      # same for every bond (flavor is uniform over qubits)
        Mp = TFIM(n, d, TH, PH, p, q, zero1=und1, zero2=und2, gadget_p=0.0); brp = Mp.run_checked([PL] * d); P0 = xdist(brp[0], n)   # optimistic: gadget noise fully mitigated
        Lam_u = (len(und1) + len(und2)) * p / 15 * len(bonds(n)) * d
        rec["stage2"][f"fixed {kind} ED+PEC optimistic"] = dict(summ(obs_from_dist(P0 / P0.sum(), n)), alpha=float(P0.sum()), Lam_u=Lam_u, cost=float(np.exp(4 * Lam_u) / P0.sum()), n_undetected=[len(und1), len(und2)])
        # realistic variant: gadgets noisy and not mitigated, payload undetected faults cancelled (bias from gadget noise remains); cost as above
        Mr = TFIM(n, d, TH, PH, p, q, zero1=und1, zero2=und2); brr = Mr.run_checked([PL] * d); P0r = xdist(brr[0], n)
        rec["stage2"][f"fixed {kind} ED+PEC payload-only"] = dict(summ(obs_from_dist(P0r / P0r.sum(), n)), alpha=float(P0r.sum()), cost=float(np.exp(4 * Lam_u) / P0r.sum()))
    # ---- stage 3
    pool = twirled_pool(n, d, TH, PH, p, q, MDRAWS, seed=d); pd = pooled(pool, d); inf3, al3, qv = twirled_estimates(pd, n, d, q, F); fin3 = twirled_finite(pool, n, d, q, F, N, R, rng)
    rec["stage3"] = {k: summ(inf3[k], fin3[k]) for k in inf3}
    rec["stage3"]["alpha_t"] = al3.tolist(); rec["stage3"]["species"] = qv.tolist()
    # 3b: twirled + first-order PEC on the invisible sector (Z_k after CX#1, Z_jZ_k after CX#2): expectation = those rates set to 0
    poolb = twirled_pool(n, d, TH, PH, p, q, MDRAWS, seed=d, inv_scale=0.0); pdb = pooled(poolb, d); infb, _, _ = twirled_estimates(pdb, n, d, q, F)
    Lam_inv = 2 * p / 15 * len(bonds(n)) * d
    rec["stage3b"] = {"twirled stratified + PEC(inv)": dict(summ(infb["twirled stratified"]), Lam_inv=Lam_inv, cost_factor=float(np.exp(4 * Lam_inv)))}
    out["depths"][d] = rec
    mx = lambda v: float(np.mean(v[:n])); I = mx(ideal)
    print(f"d={d} ({time.time()-t0:.0f}s) ideal m_x={I:.4f} | raw {mx(inf['raw'])-I:+.4f} ZNEexp {mx(inf['ZNE exp (1,2)'])-I:+.4f} (10k std of m_x {rec['stage1']['ZNE exp (1,2)']['mx_std']:.4f}) cum2 {mx(inf['ZNE cum2 (1,2,3)'])-I:+.4f} (std {rec['stage1']['ZNE cum2 (1,2,3)']['mx_std']:.4f})"
          f" | fixedX PS {mx(rec['stage2']['fixed X PS']['inf'])-I:+.4f} a={rec['stage2']['fixed X PS']['alpha']:.3f} ED+PEC opt {rec['stage2']['fixed X ED+PEC optimistic']['mx']-I:+.4f} cost {rec['stage2']['fixed X ED+PEC optimistic']['cost']:.1f} payload-only {rec['stage2']['fixed X ED+PEC payload-only']['mx']-I:+.4f}"
          f" | twirled PS {mx(inf3['twirled PS only (t=0)'])-I:+.4f} a0={al3[0]:.3f} single {mx(inf3['twirled single-species'])-I:+.4f} strat {mx(inf3['twirled stratified'])-I:+.4f} (10k std of m_x {rec['stage3']['twirled stratified']['mx_std']:.4f}) +PECinv {mx(infb['twirled stratified'])-I:+.4f}", flush=True)
    json.dump(out, open(f"demo_results_{'_'.join(map(str, depths))}.json", "w"))
