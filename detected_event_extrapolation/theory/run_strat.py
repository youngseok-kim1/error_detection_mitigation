import numpy as np, json
from edtheory import *
n, reps = 10, 5; M2 = 2 * reps; out = {}
print("E2s: gadget noise, stratified inversion over all 2^m subsets vs single-species vs reference")
for pc in (0.001, 0.004):
    for m in (2, 3, 4):
        E = Ensemble(n, reps, m, [list(range(n))] * m, [(0, M2)] * m, p=0.004, p_check=pc, pm=0.0, M=300, seed=1)
        R = Ensemble(n, reps, m, [list(range(n))] * m, [(0, M2)] * m, p=0.0, p_check=0.0, pm=0.0, M=300, seed=1, gadget_p=lambda j, side, m=m, pc=pc: pc if (j == m - 1 and side == 'R') else 0.0)
        O1, q, W, cond = E.stratified(); ss = E.single_species((1 << m) - 1); ref = R.O_all()
        out[f"{pc}_{m}"] = dict(strat=O1, single=ss, ref=ref, q=q.tolist(), cond=cond)
        print(f"  p_check={pc} m={m}: single-species {ss:.5f}   stratified {O1:.5f}   reference {ref:.5f}   cond {cond:.0f}   q_V={np.round(q, 4)}")
json.dump(out, open("strat_results.json", "w"))
