"""Nested stratification with m twirled full-weight checks: depolarizing vs biased gadget noise, with 1% readout."""
import numpy as np, json
from run_pauli import EnsG, make_noise
from edtheory import accept_stats, subset_rule
n, reps = 10, 5; M2 = 2 * reps; out = {}
def nested(E):
    m = E.k; b = E.b[0]; a = E.a[0]
    Ts = [((1 << i) - 1) << (m - i) for i in range(m + 1)]
    Mx = np.array([[b ** min(i, s) * a ** (i - min(i, s)) for s in range(m + 1)] for i in range(m + 1)])
    al = np.array([accept_stats(E.Pm, m, subset_rule(m, T))[0] for T in Ts]); alO = np.array([accept_stats(E.Pm, m, subset_rule(m, T))[1] for T in Ts])
    q = np.linalg.solve(Mx, al); W = np.linalg.solve(Mx, alO); return W[0] / q[0]
for label, gn, pm in (("depolarizing gadgets, readout 0", make_noise("depol"), 0.0), ("depolarizing gadgets, readout 1%", make_noise("depol"), 0.01),
                      ("biased gadgets, readout 1%", make_noise("biased", 9), 0.01)):
    row = []
    for m in (1, 2, 3, 4):
        E = EnsG(n, reps, m, [list(range(n))] * m, [(0, M2)] * m, 0.004, 0.004, pm, 300, make_noise("depol"), seed=1, gadget_noise=gn)
        row.append((m, E.single_species((1 << m) - 1), nested(E)))
    print(f"{label:36s}" + "   ".join(f"m={m}: single {s:.4f} nested {nn:.4f}" for m, s, nn in row)); out[label] = row
json.dump(out, open("nested2.json", "w"))
