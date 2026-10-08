"""Search for a 1-FT simultaneous S_X / S_Z extraction gadget for one [[4,2,2]] block with two ancillas
(a_z in |0>, measured in Z; a_x in |+>, measured in X) that flag each other, as in Froland et al. (arXiv:2607.24947),
whose exact circuit is only given as a figure.
Candidates: the four CX(d_i -> a_z) (Z part, fixed order d0..d3 by symmetry), the four CX(a_x -> d_i) (X part, any order),
any interleaving, and up to two CX(a_x -> a_z) flag couplings at any slot.
Criterion (distance 2, memory, one round, perfect preparation and terminal readout, in the Z and X bases):
  (a) both ancilla outcomes are deterministic (they are detectors), and
  (b) no single fault (2q depolarizing after every CX, measurement flips) flips a logical observable without firing a detector.
Output: the passing gadgets with the fewest CXs, printed as gate lists."""
import itertools, stim, sys
D, AZ, AX = [0, 1, 2, 3], 4, 5

def memory(gadget, basis, p=1e-3):
    c = stim.Circuit()
    c.append("R", range(6))
    c.append("H", [0]); c.append("CX", [0, 1, 0, 2, 0, 3])          # GHZ(4) = logical |00>  (noiseless here)
    if basis == "X": c.append("H", D)                               # H^4 GHZ = logical |++>
    c.append("H", [AX])
    for g in gadget:
        c.append("CX", g); c.append("DEPOLARIZE2", g, p)
    c.append("H", [AX])
    c.append("X_ERROR", [AZ, AX], p); c.append("M", [AZ, AX])
    c.append("DETECTOR", [stim.target_rec(-2)]); c.append("DETECTOR", [stim.target_rec(-1)])
    if basis == "X": c.append("H", D)
    c.append("M", D)
    r = lambda q: stim.target_rec(-4 + q)
    c.append("DETECTOR", [r(0), r(1), r(2), r(3)])                  # terminal stabilizer of the readout basis
    if basis == "Z":
        c.append("OBSERVABLE_INCLUDE", [r(0), r(1)], 0); c.append("OBSERVABLE_INCLUDE", [r(0), r(2)], 1)   # Zbar_0 = Z0Z1, Zbar_1 = Z0Z2
    else:
        c.append("OBSERVABLE_INCLUDE", [r(1), r(3)], 0); c.append("OBSERVABLE_INCLUDE", [r(2), r(3)], 1)   # Xbar_0 = X1X3, Xbar_1 = X2X3
    return c

def is_ft(gadget):
    for basis in ("Z", "X"):
        try: dem = memory(gadget, basis).detector_error_model()
        except ValueError: return False                                # non-deterministic detector
        for inst in dem.flattened():
            if inst.type == "error":
                t = inst.targets_copy()
                if any(x.is_logical_observable_id() for x in t) and not any(x.is_relative_detector_id() for x in t): return False
    return True

zpart = [[d, AZ] for d in D]
found = []
for xperm in itertools.permutations(D):
    xpart = [[AX, d] for d in xperm]
    for zpos in itertools.combinations(range(8), 4):
        seq, zi, xi = [], 0, 0
        for k in range(8):
            if k in zpos: seq.append(zpart[zi]); zi += 1
            else: seq.append(xpart[xi]); xi += 1
        for nf in (0, 1, 2):
            for fpos in itertools.combinations_with_replacement(range(9), nf):
                g = list(seq)
                for f in sorted(fpos, reverse=True): g.insert(f, [AX, AZ])
                if is_ft(g): found.append(g)
        if found and len(found[0]) == 8: break
    if len(found) > 50: break
print(f"{len(found)} passing gadgets found (search stopped early once enough were found)")
best = sorted(found, key=len)[:5]
for g in best: print(len(g), g)
