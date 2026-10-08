"""Drivers.  usage: python run_all.py flagship|generic:0.004|generic:0.012|ensemble4  out.jsonl"""
import sys, json, numpy as np, time
from study import study
mode, out = sys.argv[1], sys.argv[2]
TH, PH = 3 * np.pi / 8, -np.pi / 4            # angles used in arXiv:2609.13108
jobs = []
if mode == "flagship":                         # ring of 5, paper angles, depth and noise scans
    jobs = [dict(tag=f"flag d={d} p={p}", n=5, d=d, theta=TH, phi=PH, p=p, pm=p, n_train=30, seed=d) for p in (0.004, 0.012) for d in (2, 4, 6)]
elif mode.startswith("generic"):               # ring of 5, depth 4, random angles AND non-uniform Pauli rates
    p = float(mode.split(":")[1]); rng = np.random.default_rng(2026)
    ang = [(float(rng.uniform(0.2, 1.4)), -float(rng.uniform(0.2, 1.4))) for _ in range(16)]
    jobs = [dict(tag=f"generic p={p}", n=5, d=4, theta=t, phi=f, p=p, pm=p, n_train=30, seed=100 + i, noise_seed=500 + i) for i, (t, f) in enumerate(ang)]
elif mode == "ensemble4":                      # ring of 4, depth 4, random angles (fast)
    rng = np.random.default_rng(7)
    ang = [(float(rng.uniform(0.2, 1.4)), -float(rng.uniform(0.2, 1.4))) for _ in range(40)]
    jobs = [dict(tag=f"ens4 p={p}", n=4, d=4, theta=t, phi=f, p=p, pm=p, n_train=30, seed=100 + i) for p in (0.004, 0.012) for i, (t, f) in enumerate(ang)]
with open(out, "w") as f:
    for j in jobs:
        t = time.time(); tag = j.pop("tag"); r = study(**j); r.update(tag=tag, theta=j["theta"], phi=j["phi"])
        f.write(json.dumps(r) + "\n"); f.flush(); print(tag, f"{time.time()-t:.0f}s", flush=True)
