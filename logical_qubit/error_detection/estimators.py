"""Detectors from shot records, the fault-count estimate k_hat, and the estimators compared in Stage 1.

All encoded estimators apply the blockwise terminal parity (absolute S_Z of the final data readout), as [R1] does
for every cut level. o is the +-1 value of the logical observable on each shot."""
import numpy as np
from circuits import detector_spec, parity_spec, OBS, skeleton


def unpack(rec, c, obs='Z0'):
    specs, names = detector_spec(c)
    V = np.stack([np.bitwise_xor.reduce([rec[t] for t in sp]) for sp in specs], 1) if specs else np.zeros((len(rec[('d', 0)]), 0), np.uint8)
    par = np.bitwise_xor.reduce([rec[t] for t in parity_spec()]) if c.n == 6 else np.zeros_like(rec[('d', 0)])
    o = 1 - 2 * np.bitwise_xor.reduce([rec[t] for t in OBS[obs]]).astype(np.int8) if c.n == 6 else None
    return V, par, o, names

def unencoded_obs(rec, obs='Z0'):
    b = {'Z0': rec[('d', 0)], 'Z1': rec[('d', 1)], 'ZZ': rec[('d', 0)] ^ rec[('d', 1)]}[obs]
    return 1 - 2 * b.astype(np.int8)


# ------------------------------------------------------------------ fault-count estimate
class KHat:
    """k_hat(s) = minimum number of detector-error-model mechanisms whose detector signatures XOR to s
    (exhaustive BFS over all syndromes, feasible for one block: <= 2^20 syndromes)."""
    def __init__(self, c, prm, kmax=8):
        dem = skeleton(c, prm)[0].detector_error_model()
        K = dem.num_detectors; sigs = set()
        for inst in dem.flattened():
            if inst.type != 'error': continue
            m = 0
            for t in inst.targets_copy():
                if t.is_relative_detector_id(): m ^= 1 << t.val
            if m: sigs.add(m)
        self.K = K; sig = np.array(sorted(sigs), np.int64)
        dist = np.full(1 << K, 255, np.uint8); dist[0] = 0; frontier = np.array([0], np.int64)
        for k in range(1, kmax + 1):
            nxt = np.unique((frontier[:, None] ^ sig[None, :]).ravel())
            nxt = nxt[dist[nxt] == 255]
            if not len(nxt): break
            dist[nxt] = k; frontier = nxt
        self.dist = dist; self.n_mech = len(sig)
    def __call__(self, V):
        w = (V.astype(np.int64) << np.arange(V.shape[1], dtype=np.int64)[None, :]).sum(1)
        return self.dist[w]


# ------------------------------------------------------------------ estimators (each returns a float)
def m(o, sel): return o[sel].mean() if sel.any() else np.nan

def est_raw(V, par, o): return o.mean()
def est_par(V, par, o): return m(o, par == 0)
def est_ps(V, par, o): return m(o, (par == 0) & ~V.any(1))

def orp_curve(V, par, o):
    """[R1] Methods A: rank detectors by |Delta_j| (on the given shots), cumulative cuts; returns order, O_k, sigma_k, f_k"""
    ok = par == 0; Vo, oo = V[ok], o[ok]
    delta = np.array([abs(m(oo, Vo[:, j] == 0) - m(oo, Vo[:, j] == 1)) if Vo[:, j].any() else 0.0 for j in range(V.shape[1])])
    order = np.argsort(-np.nan_to_num(delta))
    keep = np.ones(len(oo), bool); Ok, sk, fk = [], [], []
    Ok.append(oo.mean()); sk.append(oo.std() / np.sqrt(len(oo))); fk.append(1.0)
    for j in order:
        keep &= Vo[:, j] == 0; n = keep.sum()
        Ok.append(oo[keep].mean() if n else np.nan); sk.append(oo[keep].std() / np.sqrt(n) if n > 1 else np.nan); fk.append(n / len(oo))
    return order, np.array(Ok), np.array(sk), np.array(fk)

def orp_select(Ok, sk, fk, fmin=0.005, fref=0.05, z=1.0):
    """plateau finder, our reading of [R1] Methods A 2: reference = inverse-variance mean of the levels with
    fmin <= f_k <= fref (most heavily filtered usable levels); plateau = longest run of consecutive k, scanned from
    large k down, with |O_k - O_ref| <= z sqrt(sigma_k^2 + sigma_ref^2) and f_k >= fmin; k* = most permissive k of it;
    fallback k* = argmin (O_k - O_ref)^2 + sigma_k^2."""
    ks = np.arange(len(Ok)); use = (fk >= fmin) & np.isfinite(Ok) & np.isfinite(sk) & (sk > 0)
    ref = use & (fk <= fref)
    if not ref.any(): ref = use & (fk <= np.quantile(fk[use], 0.25)) if use.any() else use
    if not ref.any(): return 0
    w = 1 / sk[ref] ** 2; Oref = (w * Ok[ref]).sum() / w.sum(); sref = 1 / np.sqrt(w.sum())
    good = use & (np.abs(Ok - Oref) <= z * np.sqrt(sk ** 2 + sref ** 2))
    best, cur, best_run = None, [], []
    for k in ks[::-1]:
        if good[k]: cur.append(k)
        else:
            if len(cur) > len(best_run): best_run = cur
            cur = []
    if len(cur) > len(best_run): best_run = cur
    if best_run: return min(best_run)
    return int(np.nanargmin(np.where(use, (Ok - Oref) ** 2 + sk ** 2, np.inf)))

def est_orp(V, par, o, rng):
    """held-out procedure: choose k* on one half, report <O>_{k*} on the other half (same ranking from the first half)"""
    idx = rng.permutation(len(o)); h = len(o) // 2; A, B = idx[:h], idx[h:]
    order, Ok, sk, fk = orp_curve(V[A], par[A], o[A]); k = orp_select(Ok, sk, fk)
    ok = (par[B] == 0) & ~V[B][:, order[:k]].any(1)
    return m(o[B], ok)

def est_soft(V, par, o):
    """per-detector log-linear weights (theory/qed_extrapolation.md sec. 8): rho_j = E[o|v_j=1]/E[o|v_j=0],
    O_hat = sum w o / sum w^2, w = prod rho_j^{v_j}"""
    ok = par == 0; Vo, oo = V[ok].astype(bool), o[ok].astype(float)
    rho = np.ones(V.shape[1])
    for j in range(V.shape[1]):
        a, b = m(oo, ~Vo[:, j]), m(oo, Vo[:, j])
        if np.isfinite(a) and np.isfinite(b) and abs(a) > 1e-9: rho[j] = np.clip(b / a, -1, 1)
    w = np.prod(np.where(Vo, rho[None, :], 1.0), 1)
    return (w * oo).sum() / (w * w).sum()

def est_geo(V, par, o, khat, kfit=3):
    """geometric law in the fault-count estimate: weighted fit of ln E[o | k_hat] = ln O_0 + k ln rho, k = 0..kfit"""
    ok = par == 0; k = khat(V[ok]); oo = o[ok].astype(float)
    xs, ys, ws = [], [], []
    for kk in range(kfit + 1):
        s = k == kk; n = s.sum()
        if n < 30: continue
        mu = oo[s].mean()
        if mu <= 0.02: continue
        xs.append(kk); ys.append(np.log(mu)); ws.append(n * mu * mu / max(1 - mu * mu, 1e-3))
    if len(xs) < 2: return est_ps(V, par, o)
    A = np.vstack([np.ones(len(xs)), xs]).T; W = np.diag(ws)
    coef = np.linalg.solve(A.T @ W @ A, A.T @ W @ np.array(ys))
    return float(np.exp(coef[0]))

def est_ext(V, par, o, r):
    """two-point detection extrapolation O_PS (O_PS / O_par)^r"""
    a, b = est_ps(V, par, o), est_par(V, par, o)
    return a * (a / b) ** r if a * b > 0 else a


def fit_exp(G, y):
    y = np.asarray(y, float)
    if np.all(y > 0): return float(np.exp(np.polyfit(G, np.log(y), 1)[1]))
    if np.all(y < 0): return float(-np.exp(np.polyfit(G, np.log(-y), 1)[1]))
    return float(np.polyfit(G, y, 1)[1])

def fit_cum2(G, y):
    y = np.asarray(y, float)
    if np.all(y > 0): return float(np.exp(np.polyfit(G, np.log(y), 2)[2]))
    if np.all(y < 0): return float(-np.exp(np.polyfit(G, np.log(-y), 2)[2]))
    return float(np.polyfit(G, y, 2)[2])
