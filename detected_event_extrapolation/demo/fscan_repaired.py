"""Species-count scan for the repaired full twirl at fixed depth: same pool, F = 3..7."""
import numpy as np, json, sys, time
from tfim_demo import *; from protocols import *; from repaired import TFIMR, draw_full
n, TH, PH, p, q = 5, 3*np.pi/8, -np.pi/4, 0.01, 0.01
N, R, MDRAWS = 10_000, 300, 200; out = {}
for d in [int(x) for x in sys.argv[1].split(",")]:
    rng = np.random.default_rng(d); M = TFIMR(n, d, TH, PH, p, q); pool = []
    for _ in range(MDRAWS):
        br = M.run_checked([draw_full(rng, n) for _ in range(d)]); pool.append({c: xdist(v, n) for c, v in br.items()})
    pd = pooled(pool, d); I = ideal_values(n, d, TH, PH); mI = float(np.mean(I[:n])); out[d] = {}
    for F in range(3, 8):
        inf, al, qv = twirled_estimates(pd, n, d, q, F); fin = twirled_finite(pool, n, d, q, F, N, R, np.random.default_rng(100+d))
        s = inf["twirled stratified"]; f = fin["twirled stratified"]
        b = float(np.mean(s[:n]) - mI); sd = float(f[:, :n].mean(1).std()); out[d][F] = (b, sd)
        print(f"d={d} F={F} bias {b:+.4f} std {sd:.4f} rmse {np.hypot(b,sd):.4f}", flush=True)
json.dump(out, open("fscan_repaired.json", "w"))
