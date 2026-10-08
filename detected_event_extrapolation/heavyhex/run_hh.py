"""Sweep: conventional ZNE vs Pauli-check detection-based protocols on the heavy-hex plaquette TFIM (6-ring).
usage: python3 run_hh.py <p> <depths comma-separated> [N shots] [R repetitions]

Every protocol runs on the same mediated circuit (mediators are forced by connectivity). Two check configurations:
  M  : mediator flags only (free: no extra gates; the anchor paper's checks, measured every step);
  MP : mediators + one peripheral check per data qubit per step (12 extra cZ per step, 6 extra readouts).
Estimators (m_x = mean_i <X_i>):
  raw / ZNE exp (G=1,2) / ZNE cum2 (G=1,2,3)          config M, flags ignored            [conventional]
  SV, SV+ZNE ...                                       same, post-selected on parity +1   [conventional + symmetry]
  PS[c], PS+SV[c]                                      post-selected on no flag (+ parity)
  PS+ZNE exp/cum2[c], PS+SV+ZNE exp/cum2[c]            ZNE applied to the post-selected data
  ext r[c], ext r+SV[c]                                two-point detection extrapolation O_PS (O_PS/O_all)^r with the
                                                       structural ratio r = (undetected harmful)/(detected) Pauli count
                                                       (no rates needed, assumes uniform depolarizing noise)
  PEC oracle[c], PEC+SV oracle[c]                      undetected-harmful rates set to 0 (first-order PEC with exact
                                                       rates, anchor-paper ED+PEC), variance x exp(4 Lambda_u)
Shots: N per protocol in total (split equally over gains), R repetitions; N_req = N Var / (eps^2 - bias^2).
"""
import sys, os, json, time, numpy as np
from hh_sim import *
from hh_sim import PARITY
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'demo'))
from protocols import fit_exp, fit_cum2

p = float(sys.argv[1]); depths = [int(x) for x in sys.argv[2].split(',')]
N = int(sys.argv[3]) if len(sys.argv) > 3 else 10_000; R = int(sys.argv[4]) if len(sys.argv) > 4 else 300
EPS = 0.02; KM = KMAX
CL = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'classify.json')))
LET = 'IXYZ'
def tup(s): return tuple(LET.index(c) for c in s)
def zero_set(cfg, par): return {loc: {tup(s) for s in v} for loc, v in CL[cfg]['undetected_parity' if par else 'undetected'].items()}
lam = -0.5 * np.log(1 - 2 * p / 15)                       # Pauli-Lindblad rate of one Pauli of a 2q depolarizing channel
def counts(cfg, par):
    """per Trotter step: number of (undetected harmful) and (detected) Pauli fault types"""
    t = CL[cfg]['table']; mult = {'cx1': 6, 'cx2': 6, 'cx3': 6, 'cx4': 6, 'cz1': 6, 'cz2': 6}
    nu = sum(mult[l] * sum(1 for (f, pa, h) in tt.values() if not f and h and not (par and pa)) for l, tt in t.items())
    nd = sum(mult[l] * sum(1 for (f, pa, h) in tt.values() if f or (par and pa)) for l, tt in t.items())
    return nu, nd

# ---- estimators on {G: {k: dist}} (dists may be probabilities or counts)
PAR = (PARITY > 0).astype(float)
def mval(D, ps, sv):
    v = D[0] if ps else sum(D[k] for k in range(KM + 1))
    if sv: v = v * PAR
    s = v.sum(); return mx_of(v) / s if s > 0 else np.nan
def make_estimators():
    E = {}
    E['raw'] = ('M', (1,), lambda Ds: mval(Ds[1], 0, 0))
    E['ZNE exp'] = ('M', (1, 2), lambda Ds: fit_exp([1, 2], [mval(Ds[G], 0, 0) for G in (1, 2)]))
    E['ZNE cum2'] = ('M', (1, 2, 3), lambda Ds: fit_cum2([1, 2, 3], [mval(Ds[G], 0, 0) for G in (1, 2, 3)]))
    E['SV'] = ('M', (1,), lambda Ds: mval(Ds[1], 0, 1))
    E['SV+ZNE exp'] = ('M', (1, 2), lambda Ds: fit_exp([1, 2], [mval(Ds[G], 0, 1) for G in (1, 2)]))
    E['SV+ZNE cum2'] = ('M', (1, 2, 3), lambda Ds: fit_cum2([1, 2, 3], [mval(Ds[G], 0, 1) for G in (1, 2, 3)]))
    for c in ('M', 'MP'):
        for sv in (0, 1):
            t = '+SV' if sv else ''
            E[f'PS{t}[{c}]'] = (c, (1,), lambda Ds, sv=sv: mval(Ds[1], 1, sv))
            E[f'PS{t}+ZNE exp[{c}]'] = (c, (1, 2), lambda Ds, sv=sv: fit_exp([1, 2], [mval(Ds[G], 1, sv) for G in (1, 2)]))
            E[f'PS{t}+ZNE cum2[{c}]'] = (c, (1, 2, 3), lambda Ds, sv=sv: fit_cum2([1, 2, 3], [mval(Ds[G], 1, sv) for G in (1, 2, 3)]))
            nu, nd = counts(c, sv); r = nu / nd
            E[f'ext r{t}[{c}]'] = (c, (1,), lambda Ds, sv=sv, r=r: mval(Ds[1], 1, sv) * (mval(Ds[1], 1, sv) / mval(Ds[1], 0, 0)) ** r)
            E[f'PEC{t} oracle[{c}]'] = (f'{c}-pec{t}', (1,), lambda Ds, sv=sv: mval(Ds[1], 1, sv))
    return E
EST = make_estimators()

def sample(Ds, gains, nper, rng):
    """R finite-shot copies of {G: {k: counts}}"""
    out = [dict() for _ in range(R)]
    for G in gains:
        P = np.concatenate([Ds[G][k] for k in range(KM + 1)]); P = np.clip(P, 0, None); P /= P.sum()
        C = rng.multinomial(nper, P, size=R).reshape(R, KM + 1, -1)
        for i in range(R): out[i][G] = {k: C[i, k].astype(float) for k in range(KM + 1)}
    return out

out = {"p": p, "q": p, "N": N, "R": R, "eps": EPS, "depths": {}}
fn = f"hh_p{p}.json"
if os.path.exists(fn): out = json.load(open(fn))
for d in depths:
    if str(d) in out["depths"]: continue
    t0 = time.time(); I = ideal_mx(d); rng = np.random.default_rng(1000 + d)
    runs = {}
    for c, per in (('M', False), ('MP', True)):
        runs[c] = {G: HeavyHexTFIM(d, p=p, G=G, periph=per).run() for G in (1, 2, 3)}
        for sv in (0, 1):
            t = '+SV' if sv else ''
            runs[f'{c}-pec{t}'] = {1: HeavyHexTFIM(d, p=p, periph=per, zero=zero_set(c, sv)).run()}
    rec = {"ideal": I, "alpha": {c: float(runs[c][1][0].sum()) for c in ('M', 'MP')},
           "alpha_sv": {c: float((runs[c][1][0] * PAR).sum()) for c in ('M', 'MP')}, "est": {}}
    for name, (c, gains, f) in EST.items():
        inf = f(runs[c]); fin = np.array([f(Ds) for Ds in sample(runs[c], gains, N // len(gains), rng)])
        std = float(np.nanstd(fin)); bias = float(inf - I)
        if 'oracle' in name:
            sv = '+SV' in name; nu, _ = counts(c.split('-')[0], sv); cf = float(np.exp(4 * lam * nu * d)); std *= np.sqrt(cf)
        nreq = N * std ** 2 / (EPS ** 2 - bias ** 2) if abs(bias) < EPS else float('inf')
        rec["est"][name] = dict(bias=bias, std=std, nreq=nreq, nan_frac=float(np.mean(np.isnan(fin))))
    out["depths"][str(d)] = rec
    json.dump(out, open(fn, "w"), indent=1)
    best = sorted(((v['nreq'], k) for k, v in rec["est"].items()))[:4]
    print(f"p={p} d={d} ({time.time()-t0:.0f}s) ideal {I:+.3f} alpha M {rec['alpha']['M']:.3f} MP {rec['alpha']['MP']:.3f} | "
          + " | ".join(f"{k}: b {rec['est'][k]['bias']:+.4f} s {rec['est'][k]['std']:.4f}" for k in ('raw', 'ZNE exp', 'ZNE cum2', 'PS+SV[M]', 'PS+SV+ZNE exp[M]', 'PEC+SV oracle[M]'))
          + " || cheapest: " + ", ".join(f"{k} {n/1e3:.0f}k" for n, k in best), flush=True)
