"""Repaired full twirl (any flavor P_L, controlled-X repair around anticommuting bonds) over the full depth sweep,
same noise, shots and repetitions as run_demo.py; species f <= 7."""
import numpy as np, json, time, sys
from tfim_demo import *; from protocols import *; from repaired import TFIMR, draw_full
n, TH, PH, p, q = 5, 3 * np.pi / 8, -np.pi / 4, 0.01, 0.01
depths = [int(x) for x in sys.argv[1].split(",")]; N, R, MDRAWS, F = 10_000, 300, 200, 7; out = {"F": F, "M": MDRAWS, "depths": {}}
for d in depths:
    t0 = time.time(); rng = np.random.default_rng(d); M = TFIMR(n, d, TH, PH, p, q); pool = []; gcv = gcp = 0
    for _ in range(MDRAWS):
        br = M.run_checked([draw_full(rng, n) for _ in range(d)]); pool.append({c: xdist(v, n) for c, v in br.items()}); gcv += M.gates_cv; gcp += M.gates_cp
    pd = pooled(pool, d); inf, al, qv = twirled_estimates(pd, n, d, q, F); fin = twirled_finite(pool, n, d, q, F, N, R, np.random.default_rng(100 + d))
    I = ideal_values(n, d, TH, PH); mx = lambda v: float(np.mean(v[:n]))
    rec = {}
    for k in inf:
        f = fin[k]; rec[k] = dict(inf=inf[k].tolist(), mx=mx(inf[k]), cxx=float(np.mean(inf[k][n:])), mx_std=float(f[:, :n].mean(1).std()), cxx_std=float(f[:, n:].mean(1).std()))
    rec.update(ideal=I.tolist(), alpha_t=al.tolist(), species=qv.tolist(), gates_cv=gcv / MDRAWS / d, gates_cp=gcp / MDRAWS / d)
    out["depths"][d] = rec
    s = rec["twirled stratified"]; print(f"d={d} ({time.time()-t0:.0f}s) gadget/step cP {rec['gates_cp']:.1f}+cX {rec['gates_cv']:.1f} alpha0={al[0]:.3f} q0={qv[0]:.3f} | PS-only {rec['twirled PS only (t=0)']['mx']-mx(I):+.4f} | stratified bias {s['mx']-mx(I):+.4f} std {s['mx_std']:.4f}", flush=True)
    json.dump(out, open("repaired_sweep.json", "w"))
