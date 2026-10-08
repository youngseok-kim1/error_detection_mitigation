"""Single-fault classification on the heavy-hex plaquette: insert one Pauli at one location of an otherwise noiseless
circuit and record (i) whether it raises a flag (mediator or peripheral), (ii) whether it flips the final parity
prod_i X_i, (iii) whether it changes the final X-basis distribution at all (harmful).
Output: classify.json with, per check configuration, the set of undetected-harmful Paulis per location type
(the sector that post-selection cannot remove and that PEC or extrapolation must handle)."""
import json, numpy as np
from hh_sim import *
from hh_sim import PARITY

d, STEP = 3, 1
def run(periph, ins):
    D = HeavyHexTFIM(d, p=0.0, q=0.0, periph=periph, insert=ins).run()
    return D
LET = 'IXYZ'
res = {}
for periph in (False, True):
    ref = run(periph, None)[0]
    tab = {}
    for loc in LOCS if periph else LOCS[:4]:
        for where in ((0, 3) if loc.startswith('cx') else (0,)):           # bond 0 (even sub-layer) and bond 3 (odd); data qubit 0
            for P in P2[1:]:
                D = run(periph, (STEP, where, loc, P))
                flag = 1 - D[0].sum()
                Pall = sum(D.values()); par = float(np.sum(Pall * PARITY))
                harm = float(np.abs(Pall - ref).max()) > 1e-12
                assert abs(flag) < 1e-12 or abs(flag - 1) < 1e-12, (loc, P, flag)
                assert abs(abs(par) - 1) < 1e-12, (loc, P, par)
                tab.setdefault(loc, {}).setdefault(''.join(LET[a] for a in P), []).append((round(flag), int(par < 0), int(harm)))
    # consistency between bond 0 and bond 3
    for loc in tab:
        for k, v in tab[loc].items(): assert len(set(v)) == 1, (loc, k, v)
    tab = {loc: {k: v[0] for k, v in t.items()} for loc, t in tab.items()}
    res['MP' if periph else 'M'] = tab

out = {}
for cfg, tab in res.items():
    print(f"\n=== checks: {'mediators + peripheral' if cfg == 'MP' else 'mediators only'} ===")
    print("   location : undetected & harmful (no parity) | ... also missed by parity | undetected but harmless")
    und, und_par = {}, {}
    for loc, t in tab.items():
        u = [k for k, (f, pa, h) in t.items() if not f and h]
        up = [k for k, (f, pa, h) in t.items() if not f and not pa and h]
        hl = [k for k, (f, pa, h) in t.items() if not f and not h]
        und[loc], und_par[loc] = u, up
        print(f"   {loc}: {u} | {up} | {hl}")
    nU = sum(len(v) for v in und.values()); nUp = sum(len(v) for v in und_par.values())
    nloc = len(tab); nD = sum(1 for t in tab.values() for (f, pa, h) in t.values() if f)
    print(f"   per gadget-location set: {nD} detected, {nU} undetected-harmful ({nUp} also missed by parity) of {15*nloc}")
    out[cfg] = dict(table={loc: {k: list(v) for k, v in t.items()} for loc, t in tab.items()}, undetected=und, undetected_parity=und_par)
json.dump(out, open('classify.json', 'w'), indent=1)
