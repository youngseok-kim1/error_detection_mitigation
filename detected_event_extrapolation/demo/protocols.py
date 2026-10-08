"""Estimators and finite-shot sampling for the three protocols."""
import numpy as np
from tfim_demo import *

def fit_exp(G, y):
    y = np.asarray(y); s = np.sign(y[0]) if y[0] != 0 else 1.0
    if np.all(np.sign(y) == s) and np.all(y != 0): return s * np.exp(np.polyfit(G, np.log(np.abs(y)), 1)[1])
    return np.polyfit(G, y, 1)[1]
def fit_cum2(G, y):
    y = np.asarray(y); s = np.sign(y[0]) if y[0] != 0 else 1.0
    if np.all(np.sign(y) == s) and np.all(y != 0): return s * np.exp(np.polyfit(G, np.log(np.abs(y)), 2)[2])
    return np.polyfit(G, y, 2)[2]

def sample_counts(P, N, rng):
    """multinomial counts over a probability vector"""
    return rng.multinomial(N, P / P.sum())

# ---------------- stage 1: ZNE without checks
def zne_dists(n, d, th, ph, p, Gs=(1.0, 2.0, 3.0)):
    return {G: xdist(TFIM(n, d, th, ph, p, 0.0, G=G).run_plain(), n) for G in Gs}
def zne_estimates(obs_by_G):
    G = np.array(sorted(obs_by_G)); Y = np.array([obs_by_G[g] for g in G])     # rows: gains; cols: observables
    out = {}
    out["raw"] = Y[0]
    out["ZNE exp (1,2)"] = np.array([fit_exp(G[:2], Y[:2, i]) for i in range(Y.shape[1])])
    out["ZNE cum2 (1,2,3)"] = np.array([fit_cum2(G, Y[:, i]) for i in range(Y.shape[1])])
    return out
def zne_finite(dists, n, N, R, rng):
    Gs = sorted(dists); per = N // len(Gs); ests = {}
    for _ in range(R):
        obs = {G: obs_from_dist(sample_counts(dists[G], per, rng) / per, n) for G in Gs}
        for k, v in zne_estimates(obs).items(): ests.setdefault(k, []).append(v)
    return {k: np.array(v) for k, v in ests.items()}

# ---------------- stage 2: fixed checks, full post-selection
def fixed_dists(n, d, th, ph, p, q, kind):
    M = TFIM(n, d, th, ph, p, q); fl = [fixed_flavor(n, kind)] * d; br = M.run_checked(fl)
    return {c: xdist(v, n) for c, v in br.items()}
def ps_estimate(dists_by_count, n):
    P0 = dists_by_count[0]; alpha = P0.sum(); return alpha, obs_from_dist(P0 / alpha, n)
def fixed_finite(dists_by_count, n, N, R, rng):
    counts = sorted(dists_by_count); P = np.concatenate([dists_by_count[c] for c in counts]); ests = []; alphas = []
    for _ in range(R):
        cnt = sample_counts(P, N, rng).reshape(len(counts), -1); c0 = cnt[0]; alphas.append(c0.sum() / N)
        ests.append(obs_from_dist(c0 / max(c0.sum(), 1), n))
    return np.array(ests), np.array(alphas)

# ---------------- stage 3: twirled checks, threshold family, stratified extrapolation
def twirled_pool(n, d, th, ph, p, q, Mdraws, seed, inv_scale=1.0):
    """list over flavor draws of {count: X-basis distribution (unnormalized)}"""
    rng = np.random.default_rng(seed); M = TFIM(n, d, th, ph, p, q, inv_scale=inv_scale); pool = []
    for _ in range(Mdraws):
        br = M.run_checked([draw_flavor(rng, n) for _ in range(d)]); pool.append({c: xdist(v, n) for c, v in br.items()})
    return pool
def pooled(pool, d):
    out = {c: np.zeros_like(pool[0][0]) for c in range(d + 1)}
    for br in pool:
        for c, v in br.items(): out[c] += v / len(pool)
    return out
def threshold_stats(dists_by_count, n, d):
    """alpha_t and alpha_t * O_t for t = 0..d (accept iff at most t flags)"""
    al = np.zeros(d + 1); alO = np.zeros((d + 1, 2 * n)); acc = np.zeros_like(dists_by_count[0])
    for t in range(d + 1):
        acc = acc + dists_by_count.get(t, 0); al[t] = acc.sum(); alO[t] = obs_from_dist(acc, n)    # obs_from_dist of unnormalized vector = alpha*O
    return al, alO
def stratified(al, alO, d, q, F):
    M = surv_matrix(d, q, F); qv = np.linalg.lstsq(M, al, rcond=None)[0]; W = np.linalg.lstsq(M, alO, rcond=None)[0]
    return W[0] / qv[0], qv
def single_species(al, alO):
    """assumes every faulty shot has exactly one visible faulty step (survives the t=0 rule with prob 1/2)"""
    Oall = alO[-1]; return (alO[0] - 0.5 * Oall) / (al[0] - 0.5)
def twirled_estimates(dists_by_count, n, d, q, F):
    al, alO = threshold_stats(dists_by_count, n, d); strat, qv = stratified(al, alO, d, q, F)
    return {"twirled PS only (t=0)": alO[0] / al[0], "twirled single-species": single_species(al, alO), "twirled stratified": strat}, al, qv
def twirled_finite(pool, n, d, q, F, N, R, rng):
    """each shot draws its own flavor from the pool, then (count, outcome) from that flavor's exact distribution"""
    Mdraws = len(pool); L = 2 ** n
    flat = np.array([np.concatenate([br.get(c, np.zeros(L)) for c in range(d + 1)]) for br in pool]); cdf = np.cumsum(flat, axis=1); cdf /= cdf[:, -1:]
    ests = {}
    for _ in range(R):
        m = rng.integers(0, Mdraws, N); u = rng.random(N); idx = np.empty(N, int)
        for mm in np.unique(m):
            sel = m == mm; idx[sel] = np.searchsorted(cdf[mm], u[sel])
        idx = np.minimum(idx, (d + 1) * L - 1); cnt = np.bincount(idx, minlength=(d + 1) * L).reshape(d + 1, L) / N
        dc = {c: cnt[c] for c in range(d + 1)}
        e, _, _ = twirled_estimates(dc, n, d, q, F)
        for k, v in e.items(): ests.setdefault(k, []).append(v)
    return {k: np.array(v) for k, v in ests.items()}
