"""Exact density-matrix values for the estimators that need no per-detector resolution.

Post-selection on all detectors being 0 is equivalent to every recorded ancilla outcome being 0 (no-reset rule
D_r = s_r xor s_{r-2} with s_0 = s_{-1} = 0) together with an even terminal data parity. A measurement with readout
flip q recorded as 0 maps rho -> (1-q) P0 rho P0 + q P1 rho P1, so the accepted (unnormalized) state stays one density
matrix. Modes: 'ps' (record must be 0), 'all' (no conditioning). The final data readout includes the readout flips.
Returns, for the logical observables, alpha-weighted sums from which raw / parity / PS values follow exactly."""
import numpy as np, functools
from circuits import PAULI1, PAULI2, OBS, D4

I2 = np.eye(2); X = np.array([[0, 1], [1, 0.]]); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1., -1])
PM = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
def sup(U): return np.kron(U, U.conj())

def apply(rho, S, qs, n):
    m = len(qs); ax = list(qs) + [n + q for q in qs]
    rho = np.moveaxis(rho, ax, range(2 * m)); sh = rho.shape
    return np.moveaxis((S @ rho.reshape(4 ** m, -1)).reshape(sh), range(2 * m), ax)

@functools.lru_cache(maxsize=None)
def dep1(p, excl=()):
    S = np.zeros((4, 4), complex); kept = 0.0
    for P in PAULI1:
        if P in excl: continue
        S += p / 3 * sup(PM[P]); kept += p / 3
    return S + (1 - kept) * np.eye(4)                 # removed Paulis fall back to the identity, as in sim_mc

@functools.lru_cache(maxsize=None)
def dep2(p, excl=()):
    S = np.zeros((16, 16), complex); kept = 0.0
    for P in PAULI2:
        if P in excl: continue
        S = S + p / 15 * sup(np.kron(PM[P[0]], PM[P[1]])); kept += p / 15
    return S + (1 - kept) * np.eye(16)

CXM = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0.]]); HM = np.array([[1, 1], [1, -1]]) / np.sqrt(2)


def run_dm(c, prm, mode='ps', zero=None):
    """returns the joint distribution over the final data readout (after readout flips) of the accepted ensemble:
    a vector of length 2^(#data) (unnormalized; its sum is the acceptance before the terminal parity)"""
    n = c.n; zero = zero or {}
    rho = np.zeros((2,) * (2 * n), complex); rho[(0,) * (2 * n)] = 1
    for i, op in enumerate(c.ops):
        k = op[0]
        if k == 'R': pass                                                     # all qubits start in |0>; resets only at start / after prep
        elif k == 'H': rho = apply(rho, sup(HM), [op[1]], n)
        elif k == 'CX': rho = apply(rho, sup(CXM), [op[1], op[2]], n)
        elif k == 'RZ': rho = apply(rho, sup(np.diag([np.exp(-1j * op[2] / 2), np.exp(1j * op[2] / 2)])), [op[1]], n)
        elif k == 'RX':
            th = op[2]; rho = apply(rho, sup(np.array([[np.cos(th / 2), -1j * np.sin(th / 2)], [-1j * np.sin(th / 2), np.cos(th / 2)]])), [op[1]], n)
        elif k == 'N1' and prm.p1 > 0: rho = apply(rho, dep1(prm.p1, tuple(sorted(zero.get(i, ())))), [op[1]], n)
        elif k == 'N2' and prm.p2 > 0: rho = apply(rho, dep2(prm.p2, tuple(sorted(zero.get(i, ())))), [op[1], op[2]], n)
        elif k == 'M' and op[2][0] != 'd':
            q = op[1]; P0 = sup(np.diag([1., 0])); P1 = sup(np.diag([0., 1]))
            r0, r1 = apply(rho, P0, [q], n), apply(rho, P1, [q], n)
            if mode == 'ps': rho = (1 - prm.q) * r0 + prm.q * r1
            else: rho = r0 + r1
            if op[2][0] == 'prep':                                           # verification ancilla is reset after the check
                rho = apply(rho, sup(np.array([[1, 0], [0, 0.]])) + sup(np.array([[0, 1], [0, 0.]])), [q], n)   # reset to |0>
    # final data readout: diagonal over data qubits, trace out ancillas
    data = [q for q in range(n) if any(op[0] == 'M' and op[1] == q and op[2][0] == 'd' for op in c.ops)]
    R = rho.reshape(2 ** n, 2 ** n); diag = np.real(np.diagonal(R)).reshape((2,) * n)
    anc = [q for q in range(n) if q not in data]
    P = diag.sum(axis=tuple(anc)) if anc else diag
    # readout flips on each data qubit
    for j in range(len(data)):
        P = (1 - prm.q) * P + prm.q * np.flip(P, axis=j)
    return P.reshape(-1), len(data)


def values(c, prm, obs, mode='ps', zero=None):
    """(alpha, <o>) for: 'all' = no conditioning, 'par' = even terminal parity, for the ensemble of `mode`"""
    P, nd = run_dm(c, prm, mode, zero); idx = np.arange(2 ** nd)
    bits = [(idx >> (nd - 1 - j)) & 1 for j in range(nd)]
    out = {}
    if nd == 4:
        par = (bits[0] ^ bits[1] ^ bits[2] ^ bits[3]) == 0
        for name in obs:
            sgn = 1 - 2 * np.bitwise_xor.reduce([bits[t[1]] for t in OBS[name]])
            out[name] = dict(alpha_all=P.sum(), o_all=(P * sgn).sum() / P.sum(), alpha_par=P[par].sum(), o_par=(P * sgn)[par].sum() / P[par].sum())
    else:
        for name in obs:
            sgn = {'Z0': 1 - 2 * bits[0], 'Z1': 1 - 2 * bits[1], 'ZZ': 1 - 2 * (bits[0] ^ bits[1])}[name]
            out[name] = dict(alpha_all=P.sum(), o_all=(P * sgn).sum() / P.sum())
    return out
