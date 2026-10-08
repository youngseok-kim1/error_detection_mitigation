"""Variance of the stratified estimator: threshold-LS (current) vs count-basis weighted LS, same shots."""
import numpy as np, json, time
from math import comb
from tfim_demo import *; from protocols import *
n, TH, PH, p, q = 5, 3 * np.pi / 8, -np.pi / 4, 0.01, 0.01
def count_matrix(d, q, F):
    Q = np.zeros((d + 1, F + 1))
    for f in range(F + 1):
        for a in range(f + 1):
            for c in range(d - f + 1): Q[a + c, f] += comb(f, a) / 2 ** f * comb(d - f, c) * q ** c * (1 - q) ** (d - f - c)
    return Q
def est_threshold(dc, d, F):
    al, alO = threshold_stats(dc, n, d); return stratified(al, alO, d, q, F)[0]
def est_count_wls(dc, d, F, Ncounts):
    Q = count_matrix(d, q, F); nc = np.array([dc[c].sum() for c in range(d + 1)]); Sc = np.array([obs_from_dist(dc[c], n) for c in range(d + 1)])
    w = 1 / np.sqrt(np.maximum(Ncounts * nc, 1.0))                  # multinomial: var(n_c/N) ~ n_c/N^2
    qv = np.linalg.lstsq(Q * w[:, None], nc * w, rcond=None)[0]; W = np.linalg.lstsq(Q * w[:, None], Sc * w[:, None], rcond=None)[0]
    return W[0] / qv[0]
for d in (8, 12):
    I = ideal_values(n, d, TH, PH); mx = lambda v: float(np.mean(v[:n])); rng = np.random.default_rng(1)
    pool = twirled_pool(n, d, TH, PH, p, q, 100, seed=d); L = 2 ** n; N = 10_000
    flat = np.array([np.concatenate([br.get(c, np.zeros(L)) for c in range(d + 1)]) for br in pool]); cdf = np.cumsum(flat, axis=1); cdf /= cdf[:, -1:]
    res = {F: {"thr": [], "wls": []} for F in (3, 4, 5, 6)}
    for _ in range(200):
        m = rng.integers(0, 100, N); u = rng.random(N); idx = np.empty(N, int)
        for mm in np.unique(m):
            sel = m == mm; idx[sel] = np.searchsorted(cdf[mm], u[sel])
        idx = np.minimum(idx, (d + 1) * L - 1); cnt = np.bincount(idx, minlength=(d + 1) * L).reshape(d + 1, L) / N; dc = {c: cnt[c] for c in range(d + 1)}
        for F in res: res[F]["thr"].append(mx(est_threshold(dc, d, F))); res[F]["wls"].append(mx(est_count_wls(dc, d, F, N)))
    pd = pooled(pool, d)
    for F in res:
        inf_t = mx(est_threshold(pd, d, F)); inf_w = mx(est_count_wls(pd, d, F, 1e12))
        print(f"d={d} F={F}: threshold-LS bias {inf_t-mx(I):+.4f} 10k std {np.std(res[F]['thr']):.4f} | count-WLS bias {inf_w-mx(I):+.4f} 10k std {np.std(res[F]['wls']):.4f}", flush=True)
