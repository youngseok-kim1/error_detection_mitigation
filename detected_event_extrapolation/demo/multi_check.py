"""K twirled checks per step (K ancillas bracketing the same ZZ layer): survival per visible faulty step 2^-K."""
import numpy as np, sys, time
from math import comb
from tfim_demo import *; from protocols import *

class TFIMK(TFIM):
    def __init__(self, *a, K=2, **kw):
        super().__init__(*a, **kw); self.K = K; self.N = self.n + K; self.ancs = list(range(self.n, self.n + K))
    def gadget_k(self, rho, P, anc, N):
        rg = depol_rates(self.gp * self.G)
        for qb in range(self.n):
            if P[qb]:
                U = np.zeros((4, 4), complex); U[:2, :2] = I2; U[2:, 2:] = PAULI[P[qb]]
                rho = apply(rho, sup(U), [anc, qb], N)
                if self.gp > 0: rho = apply(rho, pauli2_sup(rg), [anc, qb], N)
        return rho
    def run_checked(self, flavors):
        n, N, K = self.n, self.N, self.K; rho0 = np.ones((2 ** n, 2 ** n)) / 2 ** n
        br = {0: rho0.reshape((2,) * (2 * n)).astype(complex)}; plusK = PLUS
        for _ in range(K - 1): plusK = np.kron(plusK, PLUS)
        for step in range(self.d):
            Ps = flavors[step]; new = {}
            for c, rd in br.items():
                rd = self.rx_layer(rd, n)
                r = np.kron(rd.reshape(2 ** n, 2 ** n), plusK).reshape((2,) * (2 * N))
                for a, P in zip(self.ancs, Ps): r = self.gadget_k(r, P, a, N)
                r = self.zz_layer(r, N)
                for a, P in zip(reversed(self.ancs), reversed(Ps)): r = self.gadget_k(r, P, a, N)
                # measure all K ancillas in X basis with readout flips; branch by number of 1s
                outs = {0: r}
                for a in self.ancs:
                    nxt = {}
                    for cc, rr in outs.items():
                        rp = apply(rr, PROJ_P, [a], N); rm = apply(rr, PROJ_M, [a], N)
                        nxt[cc] = nxt.get(cc, 0) + (1 - self.q) * rp + self.q * rm; nxt[cc + 1] = nxt.get(cc + 1, 0) + self.q * rp + (1 - self.q) * rm
                    outs = nxt
                for cc, rr in outs.items():
                    rd2 = np.trace(rr.reshape(2 ** n, 2 ** K, 2 ** n, 2 ** K), axis1=1, axis2=3).reshape((2,) * (2 * n))
                    new[c + cc] = new.get(c + cc, 0) + rd2
            br = new
        return br

def surv_matrix_K(d, q, F, K):
    M = np.zeros((K * d + 1, F + 1))
    for f in range(F + 1):
        pmf = np.zeros(K * d + 1)
        for a in range(K * f + 1):
            for c in range(K * (d - f) + 1): pmf[a + c] += comb(K * f, a) / 2 ** (K * f) * comb(K * (d - f), c) * q ** c * (1 - q) ** (K * (d - f) - c)
        M[:, f] = np.cumsum(pmf)
    return M

def run(n, d, th, ph, p, q, K, Mdraws, N, R, F=5, seed=0, inv_scale=1.0):
    rng = np.random.default_rng(seed); M = TFIMK(n, d, th, ph, p, q, K=K, inv_scale=inv_scale); pool = []
    for _ in range(Mdraws):
        br = M.run_checked([[draw_flavor(rng, n) for _ in range(K)] for _ in range(d)]); pool.append({c: xdist(v, n) for c, v in br.items()})
    nc = K * d + 1; L = 2 ** n
    def est(dc):
        al = np.zeros(nc); alO = np.zeros((nc, 2 * n)); acc = np.zeros(L)
        for t in range(nc):
            acc = acc + dc.get(t, 0); al[t] = acc.sum(); alO[t] = obs_from_dist(acc, n)
        Mx = surv_matrix_K(d, q, F, K); qv = np.linalg.lstsq(Mx, al, rcond=None)[0]; W = np.linalg.lstsq(Mx, alO, rcond=None)[0]
        return W[0] / qv[0], al[0], alO[0] / al[0]
    pd = {c: sum(br.get(c, 0) for br in pool) / Mdraws for c in range(nc)}; inf, a0, ps = est(pd)
    flat = np.array([np.concatenate([br.get(c, np.zeros(L)) for c in range(nc)]) for br in pool]); cdf = np.cumsum(flat, axis=1); cdf /= cdf[:, -1:]
    ests = []
    for _ in range(R):
        m = rng.integers(0, Mdraws, N); u = rng.random(N); idx = np.empty(N, int)
        for mm in np.unique(m):
            sel = m == mm; idx[sel] = np.searchsorted(cdf[mm], u[sel])
        idx = np.minimum(idx, nc * L - 1); cnt = np.bincount(idx, minlength=nc * L).reshape(nc, L) / N
        ests.append(est({c: cnt[c] for c in range(nc)})[0])
    ests = np.array(ests); return inf, a0, ps, ests.mean(0), ests.std(0)

if __name__ == "__main__":
    n, TH, PH, p, q = 5, 3 * np.pi / 8, -np.pi / 4, 0.01, 0.01
    for d in (8,):
        I = ideal_values(n, d, TH, PH); mx = lambda v: float(np.mean(v[:n]))
        for K in (1, 2, 3):
            t = time.time(); inf, a0, ps, mean, std = run(n, d, TH, PH, p, q, K, 120, 10_000, 200, seed=d)
            print(f"d={d} K={K}: ({time.time()-t:.0f}s) alpha0={a0:.3f}  PS-only bias {mx(ps)-mx(I):+.4f}  stratified bias {mx(inf)-mx(I):+.4f}  10k-shot mean bias {mx(mean)-mx(I):+.4f} std {mx(std):.4f}  N Var {1e4*mx(std)**2:.1f}", flush=True)
