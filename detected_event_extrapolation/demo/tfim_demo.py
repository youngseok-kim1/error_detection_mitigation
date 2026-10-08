"""
Three-stage demo: standard ZNE vs fixed Pauli checks vs twirled checks with partial post-selection.

Circuit : transverse-field Ising ring, n data qubits in |+>^n.  One Trotter step =
          Rx(theta) on every data qubit, then one ZZ layer: for every bond (j,k) the gadget
          CX(j,k) Rz_k(phi) CX(j,k) = exp(-i phi/2 Z_j Z_k).  Bonds in the order even, odd, closing.
Checks  : one check ancilla, reused.  Per Trotter step it brackets the ZZ layer:
          |+>, controlled-P_L, [ZZ layer], controlled-P_L (= P_R since P_L commutes with every
          Z_j Z_k), X-basis measurement with readout flip q, reset.
          Flavor family (all commute with Z_j Z_k): Z-strings on any subset (identity allowed) and
          strings with X or Y on every qubit.  Drawing a Z-string (uniform over 2^n) with prob 1/2
          and an XY-string (uniform over 2^n) with prob 1/2 flags every net fault that is not an
          even-weight Z-string with probability exactly 1/2, and never flags even-weight Z-strings.
Noise   : two-qubit depolarizing p after every CX and every controlled-Pauli of the gadgets,
          readout flip q on the check ancilla, perfect reset and single-qubit gates.
Output  : exact unnormalized data density matrices branched by the number of flagged checks,
          so every acceptance rule "at most t flags" is evaluated exactly, and the X-basis outcome
          distribution of each branch, so finite-shot experiments can be sampled exactly.
"""
import numpy as np, functools, itertools
from math import comb

I2 = np.eye(2); X = np.array([[0, 1], [1, 0.]]); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1., -1]); H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
PAULI = [I2, X, Y, Z]
CXM = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0.]])
P2 = [(a, b) for a in range(4) for b in range(4)]
def sup(U): return np.kron(U, U.conj())
def kron2(a, b): return np.kron(a, b)

def apply(rho, S, qs, N):
    m = len(qs); ax = list(qs) + [N + q for q in qs]
    rho = np.moveaxis(rho, ax, range(2 * m)); sh = rho.shape
    return np.moveaxis((S @ rho.reshape(4 ** m, -1)).reshape(sh), range(2 * m), ax)

@functools.lru_cache(maxsize=None)
def pauli2_sup(rates):
    """two-qubit Pauli channel; rates = tuple of 16 probabilities (index 0 unused -> identity remainder)"""
    S = (1 - sum(rates[1:])) * np.eye(16, dtype=complex)
    for i in range(1, 16): S = S + rates[i] * sup(kron2(PAULI[P2[i][0]], PAULI[P2[i][1]]))
    return S

def depol_rates(p, scale=None):
    r = [0.0] + [p / 15] * 15
    if scale:
        for i, s in scale.items(): r[i] *= s
    return tuple(r)

SCX = sup(CXM); SH = sup(H)
def RX(t): return sup(np.cos(t / 2) * I2 - 1j * np.sin(t / 2) * X)
def RZ(f): return sup(np.diag([np.exp(-1j * f / 2), np.exp(1j * f / 2)]))
PLUS = np.array([[0.5, 0.5], [0.5, 0.5]])
PROJ_P = sup(np.array([[0.5, 0.5], [0.5, 0.5]])); PROJ_M = sup(np.array([[0.5, -0.5], [-0.5, 0.5]]))

def bonds(n): return [(b, b + 1) for b in range(0, n - 1, 2)] + [(b, b + 1) for b in range(1, n - 1, 2)] + ([(n - 1, 0)] if n > 2 else [])

# fault index conventions for the two CX of a bond gadget: the fault invisible to every valid check is
#   after CX#1: Z on the target k  (index of (I,Z) = 3);   after CX#2: Z_j Z_k (index of (Z,Z) = 15)
INV1, INV2 = 3, 15

class TFIM:
    def __init__(self, n, d, theta, phi, p, q, G=1.0, inv_scale=1.0, gadget_p=None, zero1=(), zero2=()):
        self.n, self.d, self.theta, self.phi, self.p, self.q = n, d, theta, phi, p, q
        self.G, self.inv_scale = G, inv_scale; self.gp = p if gadget_p is None else gadget_p
        self.zero1, self.zero2 = tuple(zero1), tuple(zero2)        # Pauli indices whose rate is set to 0 after CX#1 / CX#2 (PEC expectation)
        self.bonds = bonds(n); self.N = n + 1; self.anc = n

    # ---- pieces
    def rx_layer(self, rho, N):
        S = RX(self.theta)
        for qb in range(self.n): rho = apply(rho, S, [qb], N)
        return rho
    def zz_layer(self, rho, N):
        sc1 = {INV1: self.inv_scale}; sc1.update({i: 0.0 for i in self.zero1}); sc2 = {INV2: self.inv_scale}; sc2.update({i: 0.0 for i in self.zero2})
        r1 = depol_rates(self.p * self.G, sc1); r2 = depol_rates(self.p * self.G, sc2)
        for (j, k) in self.bonds:
            rho = apply(rho, SCX, [j, k], N); rho = apply(rho, pauli2_sup(r1), [j, k], N)
            rho = apply(rho, RZ(self.phi), [k], N)
            rho = apply(rho, SCX, [j, k], N); rho = apply(rho, pauli2_sup(r2), [j, k], N)
        return rho
    def gadget(self, rho, P, N):
        """controlled-P from the ancilla, one controlled-Pauli per non-identity site, with gadget noise"""
        rg = depol_rates(self.gp * self.G)
        for qb in range(self.n):
            if P[qb]:
                U = np.zeros((4, 4), complex); U[:2, :2] = I2; U[2:, 2:] = PAULI[P[qb]]       # ancilla is the control (first factor)
                rho = apply(rho, sup(U), [self.anc, qb], N)
                if self.gp > 0: rho = apply(rho, pauli2_sup(rg), [self.anc, qb], N)
        return rho

    # ---- full runs
    def run_plain(self):
        """no checks: data density matrix (n qubits)"""
        n = self.n; rho = np.ones((2 ** n, 2 ** n)) / 2 ** n; rho = rho.reshape((2,) * (2 * n)).astype(complex)
        for _ in range(self.d):
            rho = self.rx_layer(rho, n); rho = self.zz_layer(rho, n)
        return rho
    def run_checked(self, flavors):
        """flavors: list of d Paulis (tuples of length n, 0..3).  Returns {flag count: unnormalized data DM}"""
        n, N = self.n, self.N; rho0 = np.ones((2 ** n, 2 ** n)) / 2 ** n
        br = {0: rho0.reshape((2,) * (2 * n)).astype(complex)}
        for step in range(self.d):
            P = flavors[step]; new = {}
            for c, rd in br.items():
                rd = self.rx_layer(rd, n)
                # attach ancilla |+>
                r = np.kron(rd.reshape(2 ** n, 2 ** n), PLUS).reshape((2,) * (2 * N))
                r = self.gadget(r, P, N); r = self.zz_layer(r, N); r = self.gadget(r, P, N)
                # X-basis measurement of the ancilla, readout flip q, trace out ancilla
                rp = apply(r, PROJ_P, [self.anc], N); rm = apply(r, PROJ_M, [self.anc], N)
                tr = lambda x: np.trace(x.reshape(2 ** N, 2 ** N).reshape(2 ** n, 2, 2 ** n, 2), axis1=1, axis2=3).reshape((2,) * (2 * n))
                b0, b1 = tr(rp), tr(rm); q = self.q
                out0 = (1 - q) * b0 + q * b1; out1 = q * b0 + (1 - q) * b1
                new[c] = new.get(c, 0) + out0; new[c + 1] = new.get(c + 1, 0) + out1
            br = new
        return br

# ---- flavors
def draw_flavor(rng, n):
    if rng.random() < 0.5: return tuple(int(v) * 3 for v in rng.integers(0, 2, n))        # Z-string, identity allowed
    return tuple(int(v) for v in rng.integers(1, 3, n))                                      # X or Y on every qubit
def fixed_flavor(n, kind): return {"X": (1,) * n, "Z": (3,) * n, "Y": (2,) * n}[kind]

# ---- observables and distributions
def xdist(rd, n):
    """X-basis outcome probabilities (length 2^n, unnormalized like rd) of a data DM"""
    r = rd
    for qb in range(n): r = apply(r, SH, [qb], n)
    return np.real(np.diagonal(r.reshape(2 ** n, 2 ** n)))
def obs_from_dist(P, n):
    """<X_i> for all i and nearest-neighbour <X_i X_{i+1}> on the ring, from an X-basis distribution (normalized)"""
    b = np.arange(2 ** n); bit = lambda qb: 1 - 2 * ((b >> (n - 1 - qb)) & 1)
    xi = np.array([np.sum(P * bit(qb)) for qb in range(n)]); xx = np.array([np.sum(P * bit(qb) * bit((qb + 1) % n)) for qb in range(n)])
    return np.concatenate([xi, xx])
def ideal_values(n, d, theta, phi):
    return obs_from_dist(xdist(TFIM(n, d, theta, phi, 0.0, 0.0).run_plain(), n), n)

# ---- survival matrix for the threshold family
def surv_matrix(d, q, F):
    """M[t, f] = P[ Bin(f,1/2) + Bin(d-f, q) <= t ]  for rules t = 0..d and species f = 0..F
    (f steps with a visible net fault, each flagged with prob 1/2; d-f clean steps, each flagged with prob q)"""
    M = np.zeros((d + 1, F + 1))
    for f in range(F + 1):
        pmf = np.zeros(d + 1)
        for a in range(f + 1):
            for c in range(d - f + 1): pmf[a + c] += comb(f, a) / 2 ** f * comb(d - f, c) * q ** c * (1 - q) ** (d - f - c)
        M[:, f] = np.cumsum(pmf)
    return M

def classify_fixed(PL, j, k):
    """for a fixed check flavor PL: indices (in P2 order) of the two-qubit Paulis after CX#1 and after CX#2 of
    bond (j,k) that the check does NOT flag (their Clifford-back-propagated image commutes with PL)"""
    def anti(a, b): return (a and b and a != b)
    und1, und2 = [], []
    for i, (a, b) in enumerate(P2[1:], start=1):
        E = kron2(PAULI[a], PAULI[b]); Eb = CXM @ E @ CXM                               # back-propagate through CX#1
        # identify Eb as a Pauli pair
        for i2, (a2, b2) in enumerate(P2):
            Q = kron2(PAULI[a2], PAULI[b2])
            if np.allclose(np.abs(np.trace(Q.conj().T @ Eb)), 4): break
        bit1 = (anti(PL[j], a2) + anti(PL[k], b2)) % 2; bit2 = (anti(PL[j], a) + anti(PL[k], b)) % 2
        if not bit1: und1.append(i)
        if not bit2: und2.append(i)
    return und1, und2
