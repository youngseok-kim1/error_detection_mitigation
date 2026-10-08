"""
Heavy-hex plaquette TFIM: exact density-matrix simulator with flag branching.

Layout (paper App. B, Fig. heavyhex): 6 data qubits q0..q5 on the vertices of one heavy-hex plaquette (TFIM ring of 6),
6 mediators m_i on the edges (bond (i,i+1) is implemented through m_i), 6 peripheral ancillas p_i on the outgoing edges.

One Trotter step:
  Rx(theta) on every data qubit (noiseless);
  [peripheral checks, optional] controlled-Z from p_i (in |+>) onto q_i, for every i;
  ZZ layer, bonds (0,1),(2,3),(4,5) then (1,2),(3,4),(5,0); bond (i,j) through its mediator m (in |0>):
      CX(i->m) CX(j->m) Rz_m(phi) CX(j->m) CX(i->m) = exp(-i phi/2 Z_i Z_j) on the data,
      then m is measured in Z (flag = outcome, readout flip q) and reset;
  [peripheral checks] controlled-Z from p_i onto q_i again, X-basis readout of p_i (flag, readout flip q).
Final X-basis readout of the data (noiseless). Observable m_x = (1/6) sum_i <X_i>; the parity prod_i X_i is a symmetry
(ideal value +1) and is available for free from the same readout.

Noise: two-qubit depolarizing of total probability p (p/15 per Pauli) after every CX (on (q, m)) and after every
controlled-Z of a peripheral check (on (p_i, q_i)); readout flip q on every mediator and peripheral measurement.
Gain G scales every two-qubit rate (readout flips are not scaled).

Simulation: the data (6 qubits) + one mediator slot are evolved exactly as a density matrix; the mediator is attached in
|0> for each bond and measured/reset after it, which is exact because the mediators of different bonds act on disjoint
qubits and are used once per step. The peripheral ancillas are not simulated as qubits: their checks are valid Pauli
checks (image Z_i on q_i throughout the ZZ layer, since q_i is only ever a CX control), so for Pauli noise the flag is
exactly determined by the image rule:
  * data fault with X or Y on q_i inside the window (after the first cZ, before the second): flips flag i;
  * fault (a on p_i, e on q_i) after the first cZ: data gets e, times Z_i if a in {X,Y} (the ancilla X swaps the two
    branches, i.e. applies P_L = Z_i to the data); flag ^= [e in {X,Y}] ^ [a in {Z,Y}];
  * fault after the second cZ: data gets e (outside the window, no flag); flag ^= [a in {Z,Y}];
  * readout flip q.
The state is branched by key (k, pbits): k = number of flags so far, capped at KMAX (k = KMAX means ">= KMAX"),
pbits = the running flag bits of the 6 peripheral checks within the current step.
"""
import os, sys, functools, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'demo'))
from tfim_demo import sup, PAULI, CXM, I2, X, Z, H

ND = 6                                      # data qubits
NQ = ND + 1                                 # + mediator slot (index 6)
MED = ND
DIM = 2 ** NQ
BONDS = [(0, 1), (2, 3), (4, 5), (1, 2), (3, 4), (5, 0)]
KMAX = 3
XY = (1, 2)                                 # Pauli indices with an X component
ZY = (3, 2)                                 # Pauli indices with a Z component
PMUL = np.array([[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]])   # Pauli product up to phase
LOCS = ('cx1', 'cx2', 'cx3', 'cx4', 'cz1', 'cz2')   # location types: four CX of a bond gadget, two cZ of a peripheral check
P2 = [(a, b) for a in range(4) for b in range(4)]  # (first qubit letter, second qubit letter); index 0 = II

# ---------------------------------------------------------------- batched superoperator application
def apply_b(R, S, qs):
    """R: (B,) + (2,)*2NQ batch of density matrices; S: superoperator (4^m x 4^m) acting on qubits qs (kron(U, U*) order)."""
    m = len(qs); ax = [1 + q for q in qs] + [1 + NQ + q for q in qs]
    R = np.moveaxis(R, ax, range(1, 2 * m + 1)); sh = R.shape
    R = (S @ R.reshape(sh[0], 4 ** m, -1)).reshape(sh)
    return np.moveaxis(R, range(1, 2 * m + 1), ax)

def rx(t): return sup(np.cos(t / 2) * I2 - 1j * np.sin(t / 2) * X)
def rz(f): return sup(np.diag([np.exp(-1j * f / 2), np.exp(1j * f / 2)]))
SCX = sup(CXM)

@functools.lru_cache(maxsize=None)
def pauli_sup(idx):
    """superoperator of a 1- or 2-qubit Pauli given as tuple of letters"""
    U = PAULI[idx[0]]
    for a in idx[1:]: U = np.kron(U, PAULI[a])
    return sup(U)

def split_channel(rates, groups):
    """rates: dict Pauli-tuple -> probability. groups: dict Pauli-tuple -> (data Pauli tuple, flag bit).
    Returns {flag: superoperator on the data qubits of the channel} (identity remainder goes to the flag of 'identity')."""
    out = {}
    tot = sum(rates.values())
    for P, r in rates.items():
        if r == 0: continue
        D, f = groups[P]
        out[f] = out.get(f, 0) + r * pauli_sup(D)
    ident = tuple(0 for _ in next(iter(groups.values()))[0])
    out[0] = out.get(0, 0) + (1 - tot) * pauli_sup(ident)
    return out

# ---------------------------------------------------------------- the simulator
class HeavyHexTFIM:
    def __init__(self, d, theta=3 * np.pi / 8, phi=-np.pi / 4, p=0.005, q=None, G=1.0, periph=False,
                 zero=None, kmax=KMAX, insert=None):
        """zero: dict loc -> set of Pauli tuples whose rate is set to 0 (PEC expectation with exact rates).
        insert: (step, bond index or data qubit, loc, Pauli tuple) to insert one deterministic fault in a noiseless run."""
        self.d, self.theta, self.phi, self.p = d, theta, phi, p
        self.q = p if q is None else q
        self.G, self.periph, self.kmax = G, periph, kmax
        self.zero = zero or {}; self.insert = insert
        self._build()

    def _rates(self, loc):
        r = {P: self.p * self.G / 15 for P in P2[1:]}
        for P in self.zero.get(loc, ()): r[P] = 0.0
        return r

    def _build(self):
        # CX noise on (q, m): data letter a, mediator letter b; peripheral flag of q flips iff a in XY (if periph checks on)
        self.cx_ch = {}
        for loc in ('cx1', 'cx2', 'cx3', 'cx4'):
            g = {P: (P, int(self.periph and P[0] in XY)) for P in P2[1:]}
            self.cx_ch[loc] = split_channel(self._rates(loc), g)
        # peripheral cZ noise on (p_i, q_i): ancilla letter a, data letter e -> effective single-qubit data channel
        g1 = {(a, e): ((int(PMUL[e, 3]) if a in XY else e,), int(e in XY) ^ int(a in ZY)) for (a, e) in P2[1:]}
        g2 = {(a, e): ((e,), int(a in ZY)) for (a, e) in P2[1:]}
        self.cz_ch = {'cz1': split_channel(self._rates('cz1'), g1), 'cz2': split_channel(self._rates('cz2'), g2)}
        self.RX, self.RZ = rx(self.theta), rz(self.phi)

    # ---- branch bookkeeping: keys are (k, pbits); states are 7-qubit DMs stacked along axis 0
    @staticmethod
    def _merge(keys, R):
        uk = sorted(set(keys)); idx = {k: i for i, k in enumerate(uk)}
        out = np.zeros((len(uk),) + R.shape[1:], R.dtype)
        np.add.at(out, np.array([idx[k] for k in keys]), R)
        live = [i for i in range(len(uk)) if np.any(out[i])]                 # drop branches of exactly zero weight
        return [uk[i] for i in live], out[live]

    def _noise(self, keys, R, ch, qs, bitq):
        """apply a flag-split channel on qubits qs; flag 1 toggles peripheral bit of data qubit bitq"""
        nk, parts = [], []
        for f, S in ch.items():
            parts.append(apply_b(R, S, qs))
            nk += [(k, pb ^ (f << bitq)) for (k, pb) in keys]
        return self._merge(nk, np.concatenate(parts))

    def _ins(self, keys, R, step, where, loc, qs, bitq, kind):
        """deterministic insertion of one Pauli (classification runs)"""
        if self.insert is None: return keys, R
        s, w, l, P = self.insert
        if (s, w, l) != (step, where, loc): return keys, R
        if kind == 'cx':
            f = int(self.periph and P[0] in XY); S = pauli_sup(P)
        else:
            a, e = P
            if loc == 'cz1': D, f = ((int(PMUL[e, 3]) if a in XY else e,), int(e in XY) ^ int(a in ZY))
            else: D, f = ((e,), int(a in ZY))
            S = pauli_sup(D)
        return [(k, pb ^ (f << bitq)) for (k, pb) in keys], apply_b(R, S, qs)

    def _bond(self, keys, R, step, b):
        i, j = BONDS[b]
        seq = (([i, MED], 'cx1', i), ([j, MED], 'cx2', j), None, ([j, MED], 'cx3', j), ([i, MED], 'cx4', i))
        for item in seq:
            if item is None: R = apply_b(R, self.RZ, [MED]); continue
            qs, loc, dq = item
            R = apply_b(R, SCX, qs)
            keys, R = self._noise(keys, R, self.cx_ch[loc], qs, dq)
            keys, R = self._ins(keys, R, step, b, loc, qs, dq, 'cx')
        # measure the mediator in Z with readout flip q, reset to |0>
        r = R.reshape(len(keys), 2 ** ND, 2, 2 ** ND, 2)
        r0, r1 = r[:, :, 0, :, 0], r[:, :, 1, :, 1]; q = self.q
        good, bad = (1 - q) * r0 + q * r1, q * r0 + (1 - q) * r1
        nk = keys + [(min(k + 1, self.kmax), pb) for (k, pb) in keys]
        return self._merge(nk, self._attach(np.concatenate([good, bad])))

    @staticmethod
    def _attach(rd):
        """(B, 64, 64) data DMs -> (B,) + (2,)*14 with the mediator in |0>"""
        B = rd.shape[0]; out = np.zeros((B, 2 ** ND, 2, 2 ** ND, 2), rd.dtype); out[:, :, 0, :, 0] = rd
        return out.reshape((B,) + (2,) * (2 * NQ))

    def _cz_layer(self, keys, R, step, loc):
        for i in range(ND):
            keys, R = self._noise(keys, R, self.cz_ch[loc], [i], i)
            keys, R = self._ins(keys, R, step, i, loc, [i], i, 'cz')
        return keys, R

    def run(self):
        """returns {k: X-basis outcome distribution of the data (unnormalized, length 64)} for k = 0..kmax"""
        rd = np.ones((1, 2 ** ND, 2 ** ND), complex) / 2 ** ND          # |+>^6
        keys, R = [(0, 0)], self._attach(rd)
        for s in range(self.d):
            for qb in range(ND): R = apply_b(R, self.RX, [qb])
            if self.periph: keys, R = self._cz_layer(keys, R, s, 'cz1')
            for b in range(len(BONDS)): keys, R = self._bond(keys, R, s, b)
            if self.periph:
                keys, R = self._cz_layer(keys, R, s, 'cz2')
                # read out the peripheral ancillas: flag bits with readout flip q, then fold into k
                keys, R = self._readout_periph(keys, R)
        return self._xdists(keys, R)

    def _readout_periph(self, keys, R):
        """fold pbits (with independent readout flips q on each of the 6 bits) into k; returns keys with pbits = 0"""
        q = self.q; out_keys, parts = [], []
        # distribution of the number of observed flags given the true bits: sum of independent flips
        for idx, (k, pb) in enumerate(keys):
            ones = bin(pb).count('1'); zeros = ND - ones
            # P(observed count = c) = sum_{a,b} Bin(ones, 1-q)[a] Bin(zeros, q)[b] with a + b = c
            pa = np.array([_binom(ones, a) * (1 - q) ** a * q ** (ones - a) for a in range(ones + 1)])
            pz = np.array([_binom(zeros, b) * q ** b * (1 - q) ** (zeros - b) for b in range(zeros + 1)])
            pc = np.convolve(pa, pz)
            for c, w in enumerate(pc):
                if w == 0: continue
                out_keys.append((min(k + c, self.kmax), 0)); parts.append(w * R[idx])
        return self._merge(out_keys, np.stack(parts))

    def _xdists(self, keys, R):
        H6 = functools.reduce(np.kron, [H] * ND)
        out = {k: np.zeros(2 ** ND) for k in range(self.kmax + 1)}
        for (k, pb), r in zip(keys, R):
            assert pb == 0
            rd = r.reshape(2 ** ND, 2, 2 ** ND, 2)[:, 0, :, 0]           # mediator is |0> between bonds
            out[k] += np.real(np.diagonal(H6 @ rd @ H6))
        return out

def _binom(n, k):
    from math import comb
    return comb(n, k)

# ---------------------------------------------------------------- observables from X-basis distributions
_b = np.arange(2 ** ND)
_BITS = np.array([1 - 2 * ((_b >> (ND - 1 - i)) & 1) for i in range(ND)])     # X_i eigenvalue of each outcome
PARITY = np.prod(_BITS, axis=0)                                               # prod_i X_i
def mx_of(P):
    """alpha * m_x of an unnormalized X-basis distribution (sum over outcomes of weight * mean_i X_i)"""
    return float(np.sum(P * _BITS.mean(0)))
def ideal_mx(d, theta=3 * np.pi / 8, phi=-np.pi / 4):
    P = HeavyHexTFIM(d, theta, phi, p=0.0, q=0.0).run()[0]
    return mx_of(P) / P.sum()
