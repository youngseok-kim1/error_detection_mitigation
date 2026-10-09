"""Stage 1 sweep: one [[4,2,2]] block, two-site MFIM, encoded protocols vs unencoded ZNE.
usage: python3 stage1_block.py <p> [steps comma-separated] [schedules comma-separated: 1,2,end] [pool shots] [procs]
Output: stage1_p<p>.json (appended per configuration).

Per configuration (steps, schedule): encoded trajectories at gains G = 1, 2, 3 and with the undetected faults removed
(PEC reference), unencoded trajectories at G = 1, 2, 3. Infinite-shot values from the pool; finite-shot spread from
R = 200 subsamples of N = 10^4 shots per protocol (split equally over gains); N_req = N Var / (eps^2 - bias^2)."""
import sys, os, json, time, numpy as np
from multiprocessing import Pool
from circuits import *
from sim_mc import sample
from estimators import *
from dm_exact import values

p = float(sys.argv[1])
STEPS = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else [1, 2, 3, 4, 6, 8]
SCHED = [(int(x) if x.isdigit() else x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else [1, 2, 'end']
POOL = int(sys.argv[4]) if len(sys.argv) > 4 else 1_000_000
PROCS = int(sys.argv[5]) if len(sys.argv) > 5 else 12
N, R, EPS = 10_000, 200, 0.02
OBSN = ('Z0', 'ZZ')
GAINS = (1.0, 2.0, 3.0)


def rates(c, prm, und):
    """total Pauli-Lindblad weight of detected and undetected faults (structural, from the skeleton classification)"""
    l2 = -0.5 * np.log(1 - 2 * prm.p2 / 15); l1 = -0.5 * np.log(1 - 2 * prm.p1 / 3)
    Lu = Ld = 0.0
    for i, op in enumerate(c.ops):
        if op[0] not in ('N1', 'N2'): continue
        lam, n = (l1, 3) if op[0] == 'N1' else (l2, 15)
        nu = len(und.get(i, ())); Lu += nu * lam; Ld += (n - nu) * lam
    return Lu, Ld


def run_config(args):
    steps, sched = args; t0 = time.time(); seed = 1000 * steps + (0 if sched == 'end' else sched)
    c = encoded_mfim(steps, sched); u = unencoded_mfim(steps); I = ideal_logical(steps)
    und_all, und = classify_faults(c); Lu, Ld = rates(c, Params(p), und); Lu_all, _ = rates(c, Params(p), und_all); Ld = Ld - (Lu_all - Lu); r_struct = Lu / Ld
    kh = KHat(c, Params(p))
    enc = {G: sample(c, Params(p, G=G), POOL, seed=seed + int(10 * G)) for G in GAINS}
    pec = sample(c, Params(p), POOL, seed=seed + 7, zero=und)
    une = {G: sample(u, Params(p, G=G), POOL, seed=seed + 50 + int(10 * G)) for G in GAINS}
    q = Params(p).q; rem = {'Z0': 1 / (1 - 2 * q), 'Z1': 1 / (1 - 2 * q), 'ZZ': 1 / (1 - 2 * q) ** 2}
    out = dict(steps=steps, sched=str(sched), rounds=c.rounds, cx_enc=sum(1 for o in c.ops if o[0] == 'CX'), cx_une=2 * steps,
               Lambda_u=Lu, Lambda_u_all=Lu_all, Lambda_d=Ld, r_struct=r_struct, obs={})
    rng = np.random.default_rng(seed + 99)
    for name in OBSN:
        E = {G: unpack(enc[G], c, name)[:3] for G in GAINS}; Pp = unpack(pec, c, name)[:3]
        Uo = {G: unencoded_obs(une[G], name) * rem[name] for G in GAINS}
        def enc_est(f):                       # f(V, par, o) -> value, on each gain
            return [f(E[G][0], E[G][1], E[G][2]) for G in GAINS]
        protos = {
            'U-raw':        (('u',), lambda d: d['u'][1.0].mean()),
            'U-ZNE exp':    (('u',), lambda d: fit_exp([1, 2], [d['u'][G].mean() for G in (1.0, 2.0)])),
            'U-ZNE cum2':   (('u',), lambda d: fit_cum2([1, 2, 3], [d['u'][G].mean() for G in GAINS])),
            'E-raw':        (('e1',), lambda d: est_raw(*d['e'][1.0])),
            'E-par':        (('e1',), lambda d: est_par(*d['e'][1.0])),
            'E-PS':         (('e1',), lambda d: est_ps(*d['e'][1.0])),
            'E-ORP':        (('e1',), lambda d: est_orp(*d['e'][1.0], d['rng'])),
            'E-soft':       (('e1',), lambda d: est_soft(*d['e'][1.0])),
            'E-geo':        (('e1',), lambda d: est_geo(*d['e'][1.0], kh)),
            'E-ext':        (('e1',), lambda d: est_ext(*d['e'][1.0], r_struct)),
            'E-PS+ZNE exp': (('e',), lambda d: fit_exp([1, 2], [est_ps(*d['e'][G]) for G in (1.0, 2.0)])),
            'E-PS+ZNE cum2':(('e',), lambda d: fit_cum2([1, 2, 3], [est_ps(*d['e'][G]) for G in GAINS])),
            'E-PS+PEC':     (('pec',), lambda d: est_ps(*d['pec'])),
        }
        ngain = {'U-ZNE exp': 2, 'U-ZNE cum2': 3, 'E-PS+ZNE exp': 2, 'E-PS+ZNE cum2': 3}
        full = dict(u=Uo, e=E, pec=Pp, rng=rng)
        ex = {G: values(c, Params(p, G=G), (name,), 'ps')[name] for G in GAINS}
        ex_all = values(c, Params(p), (name,), 'all')[name]
        ex_pec = values(c, Params(p), (name,), 'ps', und)[name]
        exu = {G: values(u, Params(p, G=G), (name,), 'all')[name]['o_all'] * rem[name] for G in GAINS}
        ps1 = ex[1.0]['o_par']
        exact = {'U-raw': exu[1.0], 'U-ZNE exp': fit_exp([1, 2], [exu[1.0], exu[2.0]]), 'U-ZNE cum2': fit_cum2([1, 2, 3], [exu[G] for G in GAINS]),
                 'E-raw': ex_all['o_all'], 'E-par': ex_all['o_par'], 'E-PS': ps1,
                 'E-ext': ps1 * (ps1 / ex_all['o_par']) ** r_struct if ps1 * ex_all['o_par'] > 0 else ps1,
                 'E-PS+ZNE exp': fit_exp([1, 2], [ex[G]['o_par'] for G in (1.0, 2.0)]),
                 'E-PS+ZNE cum2': fit_cum2([1, 2, 3], [ex[G]['o_par'] for G in GAINS]), 'E-PS+PEC': ex_pec['o_par']}
        rec = {}
        for k, (_, f) in protos.items():
            inf = float(exact[k]) if k in exact else float(f(full)); gains_used = ngain.get(k, 1); nper = N // gains_used
            vals = []
            for _ in range(R):
                sub = {}
                sub['u'] = {G: Uo[G][rng.integers(0, POOL, nper)] for G in GAINS}
                sub['e'] = {}
                for G in GAINS:
                    ix = rng.integers(0, POOL, nper); sub['e'][G] = (E[G][0][ix], E[G][1][ix], E[G][2][ix])
                ix = rng.integers(0, POOL, nper); sub['pec'] = (Pp[0][ix], Pp[1][ix], Pp[2][ix]); sub['rng'] = rng
                vals.append(f(sub))
            vals = np.array(vals, float); std = float(np.nanstd(vals))
            if k == 'E-PS+PEC': std *= float(np.exp(2 * Lu))          # PEC sampling overhead gamma = exp(2 Lambda_u)
            bias = inf - I[name]
            nreq = N * std ** 2 / (EPS ** 2 - bias ** 2) if abs(bias) < EPS and np.isfinite(std) else float('inf')
            rec[k] = dict(val=inf, bias=bias, std=std, nreq=nreq, exact=k in exact, nan=float(np.mean(~np.isfinite(vals))))
        acc = float(ex[1.0]['alpha_par'])
        out['obs'][name] = dict(ideal=I[name], acc_ps=acc, est=rec)
    out['time'] = time.time() - t0
    return out


if __name__ == '__main__':
    fn = f'stage1_p{p}.json'
    res = json.load(open(fn)) if os.path.exists(fn) else {'p': p, 'configs': {}}
    todo = [(s, sc) for s in STEPS for sc in SCHED if f'{s}_{sc}' not in res['configs']]
    with Pool(PROCS) as pool:
        for out in pool.imap_unordered(run_config, todo):
            res['configs'][f"{out['steps']}_{out['sched']}"] = out
            json.dump(res, open(fn, 'w'), indent=1)
            e = out['obs']['Z0']['est']; f = lambda k: (f"{e[k]['nreq']/1e3:.0f}k" if np.isfinite(e[k]['nreq']) else 'inf') + f"({e[k]['bias']:+.3f})"
            print(f"p={p} steps={out['steps']} sched={out['sched']} ({out['time']:.0f}s) acc {out['obs']['Z0']['acc_ps']:.3f} r={out['r_struct']:.3f} | Z0: "
                  + " ".join(f"{k}:{f(k)}" for k in ('U-ZNE exp', 'U-ZNE cum2', 'E-PS', 'E-ORP', 'E-soft', 'E-geo', 'E-ext', 'E-PS+ZNE exp', 'E-PS+ZNE cum2', 'E-PS+PEC')), flush=True)
