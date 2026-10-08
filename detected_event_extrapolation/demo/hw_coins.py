"""Flag probabilities of every fault location in one ZZ layer of the ring, for periphery-ancilla checks
(one ancilla per data qubit, controlled-Pauli on its own qubit only) with a per-qubit frame twirl
(physical rotation axis of qubit i drawn uniformly from {Z,X,Y}) and each check enabled with probability b.
Logical image of check i: Z_i before/after the layer, Z_{i-1}Z_i inside gadget (i-1,i) (i is the target), Z_i inside gadget (i,i+1)."""
import itertools, numpy as np
AX = ['X', 'Y', 'Z']
def anti(a, b):  # single-qubit Paulis as letters, 'I' commutes
    return a != 'I' and b != 'I' and a != b
def flag_prob(fault, images, b):
    """fault: dict qubit->letter (physical). images: list over checks of dict qubit->'Z' (logical image, Z-type on listed qubits).
    Frame f_q in AX maps logical Z_q to physical f_q. Average over frames (uniform) and enables (Bernoulli b) of P[odd number of flags]."""
    qs = sorted(set(fault) | {q for im in images for q in im}); tot = 0.0; cnt = 0
    for frames in itertools.product(AX, repeat=len(qs)):
        fr = dict(zip(qs, frames))
        for en in itertools.product([0, 1], repeat=len(images)):
            pe = np.prod([b if e else 1 - b for e in en]); par = 0
            for e, im in zip(en, images):
                if not e: continue
                par ^= sum(anti(fault.get(q, 'I'), fr[q]) for q in im) % 2
            tot += pe * par; cnt += 1
    return tot / 3 ** len(qs)
# locations in gadget (i,i+1): A = after CX1 (images: check i -> Z_i ; check i+1 -> Z_i Z_{i+1}), B = after CX2 (both single-site)
P2 = [a + b for a in 'IXYZ' for b in 'IXYZ' if a + b != 'II']
for b in (0.5, 0.75, 1.0):
    print(f"\nenable probability b = {b}")
    for loc, ims in (("after CX1 (between CXs)", [{0: 'Z'}, {0: 'Z', 1: 'Z'}]), ("after CX2", [{0: 'Z'}, {1: 'Z'}])):
        res = {}
        for pp in P2:
            f = {0: pp[0], 1: pp[1]}; f = {q: l for q, l in f.items() if l != 'I'}
            res[pp] = round(flag_prob(f, ims, b), 4)
        vals = sorted(set(res.values())); print(f"  {loc}: coins {vals};  " + ", ".join(f"{k}:{v}" for k, v in res.items() if v != 0.5))
    # gadget-gate fault on data qubit i after the controlled-Pauli of check i (check i certainly enabled), before the layer: only check i covers it
    print("  fault on q_i after its own cP (check i enabled): coin =", round(flag_prob({0: 'X'}, [{0: 'Z'}], 1.0), 4), "(same for Y, Z)")

print("\n=== compilation twirl: gadget A (rotation on target, check i+1 spreads) and gadget B (rotation on control, check i spreads), 50/50 ===")
P2 = [a + b for a in 'IXYZ' for b in 'IXYZ' if a + b != 'II']
for b in (0.5, 0.75):
    res = {}
    for pp in P2:
        f = {q: l for q, l in {0: pp[0], 1: pp[1]}.items() if l != 'I'}
        res[pp] = round(0.5 * flag_prob(f, [{0: 'Z'}, {0: 'Z', 1: 'Z'}], b) + 0.5 * flag_prob(f, [{0: 'Z', 1: 'Z'}, {1: 'Z'}], b), 4)
    print(f"b={b} loc A (between CXs): " + ", ".join(f"{k}:{v}" for k, v in res.items()), " | mean", round(np.mean(list(res.values())), 4))
    resB = {pp: round(flag_prob({q: l for q, l in {0: pp[0], 1: pp[1]}.items() if l != 'I'}, [{0: 'Z'}, {1: 'Z'}], b), 4) for pp in P2}
    print(f"b={b} loc B (after CX2):   coins {sorted(set(resB.values()))}  mean {np.mean(list(resB.values())):.4f}")
