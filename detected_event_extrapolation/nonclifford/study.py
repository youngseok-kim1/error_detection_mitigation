"""One instance (n, d, angles, noise) -> every estimator, for both check configurations."""
import numpy as np, json, sys, time
from ising_ed import *

EPS = 1e-3

def obs_matrix(n):
    """rows: X_i (n), mean X, X_i X_{i+1} (n-1), mean XX  -> weights over X-string masks"""
    W = np.zeros((2 * n + 2, 2 ** n)); names = []
    for i in range(n): W[i, 1 << (n - 1 - i)] = 1; names.append(f"X{i}")
    W[n] = W[:n].mean(0); names.append("mx")
    for i in range(n): W[n + 1 + i, (1 << (n - 1 - i)) | (1 << (n - 1 - (i + 1) % n))] = 1; names.append(f"XX{i}")
    W[2 * n + 1] = W[n + 1:2 * n + 1].mean(0); names.append("mxx")
    return W, names

def expfit(G, Y):
    """intercept at G=0 of a*exp(-bG) (log-linear least squares); linear fit if signs differ"""
    out = np.zeros(Y.shape[1])
    for j in range(Y.shape[1]):
        y = Y[:, j]
        if np.all(y > 0) or np.all(y < 0): out[j] = np.sign(y[0]) * np.exp(np.polyfit(G, np.log(np.abs(y)), 1)[1])
        else: out[j] = np.polyfit(G, y, 1)[1]
    return out

def study(n, d, theta, phi, p, pm, n_train=0, seed=0, noise_seed=None, verbose=False):
    W, names = obs_matrix(n); O0 = W @ ideal(n, d, theta, phi)
    base = Ising(n, d, theta, phi, noise_seed=noise_seed)
    R = {g: base.run(p, pm, g, g) for g in (1.0, 2.0, 3.0, EPS)}          # independent of the check configuration
    # Clifford training circuits: same structure, every angle a random multiple of pi/2
    rng = np.random.default_rng(seed); train = []
    for _ in range(n_train):
        th = rng.integers(0, 4, d * n) * np.pi / 2; ph = rng.integers(0, 4, d * base.k) * np.pi / 2
        train.append((W @ ideal(n, d, th, ph), Ising(n, d, th, ph, noise_seed=noise_seed).run(p, pm)))
    out = {"names": names, "ideal": O0.tolist(), "n": n, "d": d, "p": p, "pm": pm}
    for pc in (False, True):
        M = Ising(n, d, theta, phi, parity_check=pc, noise_seed=noise_seed); Ld, Lu = M.rates(p, pm)
        run = lambda gu, gd=1.0: post(M.run(p, pm, gu, gd), n, pc)
        alpha, ps, al = post(R[1.0], n, pc); Ops, Oall = W @ ps, W @ al
        est = {"raw": Oall, "ED": Ops}
        # --- two-point extrapolation with post-selection as the knob
        r_cov = Lu / Ld
        _, _, al_e = post(R[EPS], n, pc); _, _, al_u = post(M.run(p, pm, EPS, 0.0), n, pc)
        S_tot = (O0 - W @ al_e) / EPS; S_u = (O0 - W @ al_u) / EPS; r1 = S_u / (S_tot - S_u)
        lin = lambda r: Ops + r * (Ops - Oall)
        est["ext, coverage ratio"] = lin(r_cov)
        est["ext, r = 1"] = lin(1.0)
        est["ext, oracle 1st-order ratio"] = lin(r1)
        ratio = np.where(Ops / Oall > 0, Ops / Oall, 1.0)
        est["ext(exp form), oracle ratio"] = np.where(Ops / Oall > 0, Ops * ratio ** r1, lin(r1))
        P, E = table(R[1.0], n, pc); Eo = E @ W.T; ok = P > 1e-14; ok[0] = False
        Rs = (Eo[ok] / P[ok, None]) / Ops
        c2 = 0.5 * np.sum((P[ok, None] / P[0]) ** 2 * (1 - Rs ** 2), axis=0)
        est["ext, oracle ratio + pair term"] = lin(r1) * np.exp(c2)
        if train:
            A = np.array([(W @ post(T, n, pc)[1]) - (W @ post(T, n, pc)[2]) for _, T in train])      # O_ps - O_all
            B = np.array([o0 - (W @ post(T, n, pc)[1]) for o0, T in train])                           # O_0  - O_ps
            r_tr = (A * B).sum(0) / np.maximum((A * A).sum(0), 1e-30)
            est["ext, ratio from Clifford training"] = lin(r_tr)
            pool = lambda sl: (A[:, sl] * B[:, sl]).sum() / (A[:, sl] ** 2).sum()
            r_pool = np.concatenate([np.full(n + 1, pool(slice(0, n))), np.full(n + 1, pool(slice(n + 1, 2 * n + 1)))])
            est["ext, site-pooled training ratio"] = lin(r_pool)
            # plain Clifford data regression baselines on the same training set (one slope through the origin)
            T0 = np.array([o0 for o0, _ in train]); Tps = T0 - B; Tall = Tps - A
            def slope(Xf):
                sl = lambda s_: (Xf[:, s_] * T0[:, s_]).sum() / (Xf[:, s_] ** 2).sum()
                return np.concatenate([np.full(n + 1, sl(slice(0, n))), np.full(n + 1, sl(slice(n + 1, 2 * n + 1)))])
            est["CDR on raw (pooled slope)"] = slope(Tall) * Oall
            est["CDR on ED (pooled slope)"] = slope(Tps) * Ops
        else: r_tr = r_pool = np.full(len(O0), np.nan)
        # --- model-based references
        a0, ps0, _ = run(0.0); est["PEC 1st order (exact model)"] = W @ ps0
        _, psm, _ = run(-0.3); est["PEC 1st order (model 30% high)"] = W @ psm
        # --- amplify only the undetected sector, post-select, extrapolate
        y = {1.0: Ops}; alphas = {1.0: alpha}
        for g in (2.0, 3.0, 2.3, 3.6): a_g, ps_g, _ = run(g); y[g] = W @ ps_g; alphas[g] = a_g
        G = np.array([1., 2., 3.])
        est["u-amp ZNE, linear (G=1,2)"] = 2 * y[1.0] - y[2.0]
        est["u-amp ZNE, exp fit (G=1,2,3)"] = expfit(G, np.array([y[1.0], y[2.0], y[3.0]]))
        est["u-amp ZNE, linear (model 30% high)"] = 2 * y[1.0] - y[2.3]
        # --- ordinary amplification of all noise
        raws = np.array([W @ post(R[g], n, pc)[2] for g in (1.0, 2.0, 3.0)]); pss = np.array([W @ post(R[g], n, pc)[1] for g in (1.0, 2.0, 3.0)])
        est["ZNE raw, exp fit (G=1,2,3)"] = expfit(G, raws)
        est["ZNE on ED, exp fit (G=1,2,3)"] = expfit(G, pss)
        # --- sampling-cost proxies (N x variance) for a +-1 observable; central site
        c = n // 2; s1 = 1 - Ops[c] ** 2; r = r1[c]
        cost = {"PEC full": float(np.exp(4 * (Ld + Lu))), "PEC 1st order + ED": float(np.exp(4 * Lu) / alpha),
                "ED": float(s1 / alpha), "ext two-point": float((1 + r) ** 2 * s1 / alpha + r ** 2 * (1 - Oall[c] ** 2) - 2 * r * (1 + r) * s1),
                "u-amp linear": float((2 * np.sqrt(s1) + np.sqrt(1 - y[2.0][c] ** 2)) ** 2 / alpha)}
        wfit = np.array([4 / 3, 1 / 3, -2 / 3])
        def cost_exp(Y, A):
            yc = np.array([v[c] for v in Y]); return float(3 * np.sum(wfit ** 2 * (1 - yc ** 2) / (np.array(A) * yc ** 2)) * expfit(G, np.array(Y))[c] ** 2)
        cost["u-amp exp fit"] = cost_exp([y[1.0], y[2.0], y[3.0]], [alpha] * 3)
        cost["ZNE raw exp fit"] = cost_exp(list(raws), [1, 1, 1])
        cost["ZNE on ED exp fit"] = cost_exp(list(pss), [post(R[g], n, pc)[0] for g in (1.0, 2.0, 3.0)])
        subs = []
        k = len(P).bit_length() - 1
        for mask in range(len(P)):
            sel = (np.arange(len(P)) & mask) == 0
            subs.append((float(-np.log(P[sel].sum() / P.sum())), (Eo[sel].sum(0) / P[sel].sum()).tolist()))
        out["parity" if pc else "anc"] = dict(est={a: b.tolist() for a, b in est.items()}, cost=cost, alpha=float(alpha), Ld=float(Ld), Lu=float(Lu),
            coverage=float(Ld / (Ld + Lu)), r_cov=float(r_cov), r1=r1.tolist(), r_train=np.asarray(r_tr).tolist(), r_pool=np.asarray(r_pool).tolist(), c2=c2.tolist(),
            alpha_uamp={str(g): float(a) for g, a in alphas.items()}, alpha_gu0=float(a0), subsets=subs,
            alpha_G={str(g): float(post(R[g], n, pc)[0]) for g in (1.0, 2.0, 3.0)})
    return out

if __name__ == "__main__":
    t = time.time(); o = study(4, 3, 3 * np.pi / 8, -np.pi / 4, 0.004, 0.004, n_train=30, seed=1)
    print(f"{time.time()-t:.1f}s")
    for cfg in ("anc", "parity"):
        c = o[cfg]; print(f"\n[{cfg}] coverage={c['coverage']:.3f} alpha={c['alpha']:.3f} Lam={c['Ld']+c['Lu']:.3f} r_cov={c['r_cov']:.3f} r1(X_i)={np.round(c['r1'][:4],3)} r_train={np.round(c['r_train'][:4],3)}")
        print("   alpha under u-amp:", c["alpha_uamp"], " cost:", {a: round(b, 2) for a, b in c["cost"].items()})
        I = np.array(o["ideal"])
        for name, v in c["est"].items():
            b = np.array(v) - I; print(f"   {name:36s} X_i rms {np.sqrt((b[:4]**2).mean()):.5f}   m_x {b[4]:+.5f}   XX rms {np.sqrt((b[5:9]**2).mean()):.5f}")
    print("ideal X_i:", np.round(o["ideal"][:5], 4))
