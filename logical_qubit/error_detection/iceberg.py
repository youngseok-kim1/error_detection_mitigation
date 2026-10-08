"""One [[4,2,2]] Iceberg block in stim: preparation, syndrome extraction, detectors, terminal readout.

Code: S_X = X0X1X2X3, S_Z = Z0Z1Z2Z3; Xbar_i = X_{i+1} X_3, Zbar_i = Z_0 Z_{i+1} (i = 0, 1); logical |00> = GHZ(4),
logical |++> = H^4 GHZ(4). Qubits: data d0..d3 = 0..3, a_z = 4 (measures S_Z, prepared |0>, read in Z),
a_x = 5 (measures S_X, prepared |+>, read in X).

Syndrome extraction (EXTRACT): eight CX, Z part CX(d -> a_z) and X part CX(a_x -> d) interleaved so that each ancilla
flags the other's hook errors through the data (found by search_extraction.py; 1-FT for memory by exhaustive single-fault
check, test_iceberg.py). Froland et al. (arXiv:2607.24947) use the same two-ancilla, mutual-flag structure; their exact
gate order is only given as a figure and is not reproduced here. Connectivity (heavy-hex embedding, compute/syndrome
configuration switching) is not modelled in Stage 0.

Detectors. sigma_r = instantaneous stabilizer value reported in round r.
  with reset    : s_r = sigma_r                ; detector D_r = s_r xor s_{r-1}
  without reset : s_r = s_{r-1} xor sigma_r    ; detector D_r = s_r xor s_{r-2}   (two rounds apart, as in Froland et al.)
with s_0 = s_{-1} = 0 (both stabilizers are fixed to +1 by the GHZ preparation). The terminal data readout gives the
stabilizer of the readout basis, sigma_T, closing D_T = sigma_T xor sigma_R.

Noise (circuit-level, all Pauli): DEPOLARIZE2(p) after every CX, DEPOLARIZE1(p1) after every single-qubit gate and
reset, X_ERROR(q) before every measurement. Gain G multiplies p, p1 (not q). Defaults p1 = p/10, q = p.
"""
import numpy as np, stim

D, AZ, AX = [0, 1, 2, 3], 4, 5
EXTRACT = [[0, 4], [5, 0], [5, 1], [1, 4], [2, 4], [5, 2], [5, 3], [3, 4]]
GHZ_CHAIN = [[0, 1], [1, 2], [2, 3]]                 # non-FT GHZ(4): depth 3 chain
GHZ_VERIFY = [[0, 4], [3, 4]]                        # Z0Z3 parity onto a_z: catches the weight-2 spread of the chain


class Noise:
    def __init__(self, p=1e-3, p1=None, q=None, G=1.0):
        self.p, self.p1, self.q, self.G = p * G, (p / 10 if p1 is None else p1) * G, (p if q is None else q), G


def _cx(c, pairs, nz):
    for g in pairs:
        c.append("CX", g)
        if nz.p > 0: c.append("DEPOLARIZE2", g, nz.p)

def _1q(c, gate, qs, nz):
    c.append(gate, qs)
    if nz.p1 > 0: c.append("DEPOLARIZE1", qs, nz.p1)

def _meas(c, basis, qs, nz):
    if basis == "X": _1q(c, "H", qs, nz)
    if nz.q > 0: c.append("X_ERROR", qs, nz.q)
    c.append("M", qs)


def memory_circuit(basis="Z", rounds=1, noise=None, reset=False, prep="nft", obs=(0, 1), absolute=False):
    """Iceberg memory experiment.
    basis : 'Z' (logical |00>, terminal Z readout, observables Zbar_i) or 'X' (logical |++>, X readout, Xbar_i)
    prep  : 'perfect' (noiseless GHZ), 'nft' (noisy chain), 'ft' (noisy chain + Z0Z3 verification detector)
    Returns (circuit, meta) where meta lists detector labels (kind, round) in order."""
    nz = noise or Noise(0.0, 0.0, 0.0)
    c = stim.Circuit(); labels = []
    c.append("R", range(6))
    pz = nz if prep != "perfect" else Noise(0.0, 0.0, 0.0)
    if pz.p1 > 0: c.append("DEPOLARIZE1", range(6), pz.p1)
    _1q(c, "H", [0], pz); _cx(c, GHZ_CHAIN, pz)
    if prep == "ft":
        _cx(c, GHZ_VERIFY, pz); _meas(c, "Z", [AZ], pz)
        c.append("DETECTOR", [stim.target_rec(-1)], (AZ, -1)); labels.append(("prep", -1))
        c.append("R", [AZ])
        if pz.p1 > 0: c.append("DEPOLARIZE1", [AZ], pz.p1)
    if basis == "X": _1q(c, "H", D, pz)
    nmeas = 0; srec = {"Z": [], "X": []}               # absolute measurement indices of s_r per stabilizer
    m_count = 1 if prep == "ft" else 0
    for r in range(rounds):
        _1q(c, "H", [AX], nz)
        _cx(c, EXTRACT, nz)
        _1q(c, "H", [AX], nz)
        if nz.q > 0: c.append("X_ERROR", [AZ, AX], nz.q)
        c.append("M", [AZ, AX]); srec["Z"].append(m_count); srec["X"].append(m_count + 1); m_count += 2
        if reset:
            c.append("R", [AZ, AX])
            if nz.p1 > 0: c.append("DEPOLARIZE1", [AZ, AX], nz.p1)
        for g, col in (("Z", AZ), ("X", AX)):
            cur = [srec[g][r]]
            lag = 1 if reset else 2
            prev = [srec[g][r - lag]] if r - lag >= 0 else []
            tg = [stim.target_rec(i - m_count) for i in cur + prev]
            c.append("DETECTOR", tg, (col, r)); labels.append((g, r))
    _meas(c, basis, D, nz); m_count += 4
    data = [stim.target_rec(-4 + i) for i in range(4)]
    # terminal detector: sigma_T xor sigma_R, sigma_R = s_R (reset) or s_R xor s_{R-1} (no reset)
    g = basis
    if rounds == 0: prevs = []
    elif reset: prevs = [srec[g][-1]]
    else: prevs = srec[g][-2:]
    c.append("DETECTOR", data + [stim.target_rec(i - m_count) for i in prevs], (9, rounds)); labels.append(("T" + g, rounds))
    if absolute:
        # absolute stabilizer values against the preparation (+1): terminal parity sigma_T and last-round sigma_R.
        # These are the checks of [R1] App. B ("post-select on the final round only" + blockwise terminal parity):
        # any odd number of faults since preparation fires them, so only cancelling pairs survive.
        c.append("DETECTOR", data, (9, rounds + 1)); labels.append(("absT" + g, rounds))
        if rounds:
            cur = [srec[g][-1]] if reset else ([srec[g][-1]] + ([srec[g][-2]] if rounds >= 2 else []))
            c.append("DETECTOR", [stim.target_rec(i - m_count) for i in cur], (9, rounds + 2)); labels.append(("absR" + g, rounds - 1))
    LOG = {"Z": [[0, 1], [0, 2]], "X": [[1, 3], [2, 3]]}[basis]
    for k in obs: c.append("OBSERVABLE_INCLUDE", [data[i] for i in LOG[k]], k)
    return c, {"labels": labels, "basis": basis, "rounds": rounds, "reset": reset, "prep": prep}


def single_fault_report(circuit):
    """Exhaustive single-fault check from the detector error model: returns the list of error mechanisms that flip
    an observable without firing any detector (empty = no first-order undetected logical error)."""
    bad = []
    for inst in circuit.detector_error_model().flattened():
        if inst.type != "error": continue
        t = inst.targets_copy()
        if any(x.is_logical_observable_id() for x in t) and not any(x.is_relative_detector_id() for x in t):
            bad.append((inst.args_copy()[0], [x.val for x in t if x.is_logical_observable_id()]))
    return bad


def dem_mechanisms(circuit):
    """Detector error model as arrays: probabilities p_m, detector incidence H (M x K), observable incidence L (M x n_obs)."""
    dem = circuit.detector_error_model()
    K, nO = dem.num_detectors, dem.num_observables
    P, Hs, Ls = [], [], []
    for inst in dem.flattened():
        if inst.type != "error": continue
        h = np.zeros(K, np.uint8); l = np.zeros(nO, np.uint8)
        for x in inst.targets_copy():
            if x.is_relative_detector_id(): h[x.val] ^= 1
            elif x.is_logical_observable_id(): l[x.val] ^= 1
        P.append(inst.args_copy()[0]); Hs.append(h); Ls.append(l)
    return np.array(P), np.array(Hs).astype(np.int64), np.array(Ls).astype(np.int64)
