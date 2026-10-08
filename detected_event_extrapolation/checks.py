import numpy as np
from edzne import *
n, k, reps, p, wL = 10, 6, 5, 0.004, 2
c, n2 = build(n, k, reps, p, p, wL, seed=3)
b, info, cls, P = analyse(c, n2, k, p, p); F = 1 << k

# (1) Monte Carlo check of the exact DEM/Walsh-Hadamard evaluation
N = 4_000_000
det, obs = c.compile_detector_sampler(seed=5).sample(N, separate_observables=True)
acc = ~det.any(axis=1); o = 1 - 2 * obs[:, 0].astype(float)
print(f"sampled: alpha={acc.mean():.4f} (exact {info['alpha']:.4f})  <O>raw={o.mean():.4f} (exact {info['O_all']:.4f})  <O>ED={o[acc].mean():.4f} (exact {info['O_ps']:.4f})")
est = o[acc].mean() ** (1 + info['rO']) / o.mean() ** info['rO']
# batch std of the two-point estimator
B = 200; e = [o[i::B][acc[i::B]].mean() ** (1 + info['rO']) / o[i::B].mean() ** info['rO'] for i in range(B)]
print(f"two-point estimator on samples: {est:.4f}  (batch std x sqrt(N/B) -> N*Var = {np.var(e) * N / B:.3f})")

# (2) physical noise amplification on post-selected data: exponent is a polynomial in G
print("\nG   alpha(G)   -ln(alpha)/G   ln<O>_raw/G   ln<O>_ED/G")
Gs = np.array([0.5, 1, 1.5, 2, 3]); ya, yp, xa = [], [], []
for G in Gs:
    PG = joint(cls, k, G); a, op = ps_stats(PG, k, F - 1); _, oa = ps_stats(PG, k, 0)
    ya.append(np.log(oa)); yp.append(np.log(op)); xa.append(-np.log(a))
    print(f"{G:3.1f}  {a:7.4f}   {-np.log(a)/G:8.4f}     {np.log(oa)/G:8.4f}     {np.log(op)/G:8.4f}")
ya, yp = np.array(ya), np.array(yp)
sel = Gs >= 1
lin = np.polyfit(Gs[sel], yp[sel], 1)[1]; quad = np.polyfit(Gs[sel], yp[sel], 2)[2]
print(f"ED data, G in {Gs[sel]}: exponential fit intercept bias {np.exp(lin)-1:+.4f}; exponent quadratic in G: {np.exp(quad)-1:+.5f}")
print(f"first-order coefficient 2*Lam_O^u = {2*info['fu']*info['Lam']*(1-info['cov']):.4f}; ED slope at G->0 from fit = {-np.polyfit(Gs, yp, 3)[2]:.4f}")
print(f"relative N*Var at G=3: raw ~ exp(-2 y_raw) = {np.exp(-2*ya[-1]):.1f}; ED ~ exp(-2 y_ED)/alpha = {np.exp(-2*yp[-1]+xa[-1]):.1f}")
