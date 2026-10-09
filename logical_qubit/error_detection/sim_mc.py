"""Batched Pauli-fault trajectory simulator for circuits.OpList (n <= ~13 qubits).

Every shot carries its own statevector (B x 2^n). Noise is sampled per shot: N1 = single-qubit depolarizing p1
(X, Y, Z each p1/3), N2 = two-qubit depolarizing p2 (each of 15 Paulis p2/15), M = projective Z measurement with
classical readout flip q (the qubit is not reset and keeps its post-measurement state). Gain G is in Params.
`zero` (dict noise-op index -> set of Pauli strings) removes those Paulis from the channel: the infinite-shot
expectation of PEC that cancels them exactly.
Returns the measurement record as a dict tag -> (B,) uint8 array."""
import numpy as np
from circuits import PAULI1, PAULI2

LET = {'I': 0, 'X': 1, 'Y': 2, 'Z': 3}


class Sim:
    def __init__(self, n):
        self.n = n; N = 2 ** n; idx = np.arange(N)
        self.bit = [(idx >> (n - 1 - q)) & 1 for q in range(n)]
        self.flip = [idx ^ (1 << (n - 1 - q)) for q in range(n)]
        self.sign = [1 - 2 * b for b in self.bit]
        self.i0 = [idx[b == 0] for b in self.bit]; self.i1 = [idx[b == 1] for b in self.bit]

    def u1(self, S, q, U):
        a0, a1 = S[:, self.i0[q]], S[:, self.i1[q]]
        S[:, self.i0[q]], S[:, self.i1[q]] = U[0, 0] * a0 + U[0, 1] * a1, U[1, 0] * a0 + U[1, 1] * a1

    def cx(self, S, c, t):
        sel = self.bit[c] == 1; idx = np.arange(2 ** self.n); perm = np.where(sel, self.flip[t], idx)
        S[:] = S[:, perm]

    def pauli(self, S, mask, q, P):
        """apply Pauli letter P (1..3) on qubit q to the shots in mask"""
        if not mask.any(): return
        if P in (1, 2): S[mask] = S[mask][:, self.flip[q]]
        if P in (2, 3): S[mask] *= self.sign[q]

    def run(self, c, prm, B, rng, zero=None):
        n = self.n; S = np.zeros((B, 2 ** n), complex); S[:, 0] = 1
        rec = {}; zero = zero or {}
        H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
        for i, op in enumerate(c.ops):
            k = op[0]
            if k == 'R':
                q = op[1]                                            # reset: measure and flip back to |0> (only used at start)
                p1 = (np.abs(S[:, self.i1[q]]) ** 2).sum(1); one = rng.random(B) < p1
                self.pauli(S, one, q, 1)
            elif k == 'H': self.u1(S, op[1], H)
            elif k == 'CX': self.cx(S, op[1], op[2])
            elif k == 'RZ':
                th = op[2]; self.u1(S, op[1], np.diag([np.exp(-1j * th / 2), np.exp(1j * th / 2)]))
            elif k == 'RX':
                th = op[2]; self.u1(S, op[1], np.array([[np.cos(th / 2), -1j * np.sin(th / 2)], [-1j * np.sin(th / 2), np.cos(th / 2)]]))
            elif k == 'N1':
                if prm.p1 <= 0: continue
                hit = rng.random(B) < prm.p1; which = rng.integers(0, 3, B)
                for w in range(3):
                    if PAULI1[w] in zero.get(i, ()): continue
                    self.pauli(S, hit & (which == w), op[1], w + 1)
            elif k == 'N2':
                if prm.p2 <= 0: continue
                hit = rng.random(B) < prm.p2; which = rng.integers(0, 15, B)
                for w in range(15):
                    P = PAULI2[w]
                    if P in zero.get(i, ()): continue
                    m = hit & (which == w)
                    if not m.any(): continue
                    for L, q in zip(P, op[1:]):
                        if L != 'I': self.pauli(S, m, q, LET[L])
            elif k == 'M':
                q = op[1]; p1 = (np.abs(S[:, self.i1[q]]) ** 2).sum(1)
                out = (rng.random(B) < p1)
                S[out[:, None] & (self.bit[q] == 0)[None, :]] = 0; S[~out[:, None] & (self.bit[q] == 1)[None, :]] = 0
                S /= np.linalg.norm(S, axis=1, keepdims=True)
                flipped = out ^ (rng.random(B) < prm.q)
                rec[op[2]] = flipped.astype(np.uint8)
        return rec


def sample(c, prm, shots, seed=0, batch=50_000, zero=None):
    """shot records for `shots` shots: dict tag -> (shots,) uint8"""
    rng = np.random.default_rng(seed); sim = Sim(c.n); out = {}
    done = 0
    while done < shots:
        b = min(batch, shots - done); r = sim.run(c, prm, b, rng, zero)
        for t, v in r.items(): out.setdefault(t, []).append(v)
        done += b
    return {t: np.concatenate(v) for t, v in out.items()}
