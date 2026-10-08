"""Nested stratification: species = detection multiplicity, points = the m+1 nested subsets (outermost i checks)."""
import numpy as np, json
from edtheory import *
n, reps = 10, 5; M2 = 2 * reps; out = {}
def nested(E):
    m = E.k; b = E.b[0]
    # species s = 0..m by multiplicity (number of checks that can see the fault); survival under the outermost-i subset = b^min(i, s)
    Ts = [((1 << i) - 1) << (m - i) for i in range(m + 1)]                   # outermost i checks (check m-1 is outermost)
    Mx = np.array([[b ** min(i, s) for s in range(m + 1)] for i in range(m + 1)])
    al = np.array([accept_stats(E.Pm, m, subset_rule(m, T))[0] for T in Ts]); alO = np.array([accept_stats(E.Pm, m, subset_rule(m, T))[1] for T in Ts])
    q = np.linalg.solve(Mx, al); W = np.linalg.solve(Mx, alO); return W[0] / q[0], np.linalg.cond(Mx), q
print("nested stratification (m+1 points, m+1 species) under gadget noise")
for pc in (0.001, 0.004):
    for m in (2, 3, 4, 6):
        E = Ensemble(n, reps, m, [list(range(n))] * m, [(0, M2)] * m, p=0.004, p_check=pc, pm=0.0, M=300, seed=1)
        O1, cond, q = nested(E); ss = E.single_species((1 << m) - 1)
        out[f"{pc}_{m}"] = dict(nested=O1, single=ss, cond=cond, q=q.tolist()); print(f"  p_check={pc} m={m}: single-species {ss:.5f}  nested-stratified {O1:.5f}  cond {cond:.0f}  q_s={np.round(q, 4)}")
json.dump(out, open("nested_results.json", "w"))
