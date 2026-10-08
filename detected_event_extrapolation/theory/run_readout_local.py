"""Readout error on local twirled checks: stratified inversion with calibrated vs mis-calibrated readout error."""
import numpy as np, json
from run_pauli import EnsG, make_noise
n, reps = 10, 5; M2 = 2 * reps; k3 = 3; sup = [[0, 1], [4, 5], [7, 8]]; out = {}
for pm, pm_as in ((0.0, 0.0), (0.01, 0.01), (0.01, 0.0), (0.01, 0.02), (0.03, 0.03), (0.03, 0.0)):
    E = EnsG(n, reps, k3, sup, [(0, M2)] * k3, 0.004, 0.0, pm, 600, make_noise("depol"), seed=2, pm_assumed=pm_as)
    O1, q, W, cond = E.stratified(); ss = E.single_species(7)
    out[f"{pm}_{pm_as}"] = dict(strat=O1, single=ss, q0=q[0]); print(f"readout {pm*100:.0f}%, estimator assumes {pm_as*100:.0f}%:  stratified O_1 = {O1:.5f}   (single-species, all 3 checks: {ss:.4f};  species-0 weight {q[0]:.4f})")
json.dump(out, open("readout_local.json", "w"))
