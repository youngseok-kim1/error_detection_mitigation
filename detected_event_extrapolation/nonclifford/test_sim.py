"""Sanity checks of the simulator (run: python test_sim.py)"""
import numpy as np, time
from scipy.linalg import expm
from ising_ed import *
th, ph = 3 * np.pi / 8, -np.pi / 4
def kronall(ms):
    out = np.array([[1.]])
    for m in ms: out = np.kron(out, m)
    return out
# 1) ideal(): data-only statevector vs brute-force matrix exponentials (chain and ring)
for ring in (False, True):
    n, d = 4, 3; psi = np.ones(2 ** n, complex) / np.sqrt(2 ** n)
    for _ in range(d):
        for (j, k) in bond_list(n, ring): psi = expm(-1j * ph / 2 * kronall([Z if q in (j, k) else I2 for q in range(n)])) @ psi
        for q in range(n): psi = expm(-1j * th / 2 * kronall([X if r == q else I2 for r in range(n)])) @ psi
    ref = np.array([np.real(psi.conj() @ kronall([X if r == q else I2 for r in range(n)]) @ psi) for q in range(n)])
    idl = ideal(n, d, th, ph, ring=ring); mine = np.array([idl[1 << (n - 1 - q)] for q in range(n)])
    T0 = Ising(n, d, th, ph, ring=ring).run(0, 0)
    print(f"ring={ring}: |ideal - expm| = {np.abs(mine-ref).max():.1e};  noiseless circuit with ancillas: P(s=0) = {T0[0,0]:.12f}, |run - ideal| over all X strings = {np.abs(T0[0]-idl).max():.1e}")
# 2) fault classification: with only 'undetected' faults switched on, acceptance must be exactly 1
for n, d in ((4, 3), (5, 2)):
    for pc in (False, True):
        for ns in (None, 3):
            M = Ising(n, d, th, ph, parity_check=pc, noise_seed=ns); a = post(M.run(0.01, 0.01, 1, 0), n, pc)[0]
            Ld, Lu = M.rates(0.01, 0); print(f"n={n} d={d} parity={pc} nonuniform={ns is not None}: coverage {Ld/(Ld+Lu):.3f}, acceptance with only undetected faults = {a:.12f}")
# 3) with undetected faults removed (what 1st-order PEC does on average) the ED bias must be O(p^2)
n, d = 4, 3; idl = ideal(n, d, th, ph); M = Ising(n, d, th, ph); m0 = 1 << (n - 1)
for p in (0.002, 0.004, 0.008):
    print(f"p={p}: ED bias X_0 = {post(M.run(p, p), n, False)[1][m0]-idl[m0]:+.2e}   undetected removed: {post(M.run(p, p, 0, 1), n, False)[1][m0]-idl[m0]:+.2e}")
# 4) noise channel at full rate == two-qubit depolarizing with f = exp(-16 p / 15)
p = 0.03; S = noise_sup(tuple([False] + [True] * 15), tuple([0] + [p / 15] * 15), 1.0, 1.0); f = np.exp(-16 * p / 15)
rng = np.random.default_rng(0); A = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)); r = A @ A.conj().T; r /= np.trace(r)
print("depolarizing check:", np.abs((S @ r.reshape(-1)).reshape(4, 4) - (f * r + (1 - f) * np.eye(4) / 4)).max())
