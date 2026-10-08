"""Checks of hh_sim.py.
1. ideal m_x against exact matrix exponentials of the 6-site TFIM ring;
2. noiseless circuits never flag; probabilities sum to one with noise;
3. the image rule used for the peripheral checks against an explicit simulation of the ancilla
   (one bond, peripheral check on q_i, all qubits explicit, depolarizing noise everywhere)."""
import numpy as np, functools
from scipy.linalg import expm
from hh_sim import *
from hh_sim import _BITS
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'demo'))
from tfim_demo import apply, pauli2_sup, depol_rates

TH, PH = 3 * np.pi / 8, -np.pi / 4
def op(P, i):
    return functools.reduce(np.kron, [P if k == i else I2 for k in range(ND)])
def exact_mx(d):
    Hzz = sum(op(Z, i) @ op(Z, j) for (i, j) in BONDS); Hx = sum(op(X, i) for i in range(ND))
    Uzz = expm(-1j * PH / 2 * Hzz); Ux = expm(-1j * TH / 2 * Hx)
    psi = np.ones(2 ** ND) / 2 ** (ND / 2)
    for _ in range(d): psi = Uzz @ (Ux @ psi)
    return np.mean([np.real(psi.conj() @ op(X, i) @ psi) for i in range(ND)])

print("1. ideal m_x:")
for d in (1, 2, 3, 5):
    a, b = ideal_mx(d), exact_mx(d); print(f"   d={d}: sim {a:+.10f} exact {b:+.10f}  diff {abs(a-b):.1e}"); assert abs(a - b) < 1e-10

print("2. normalization / noiseless flags:")
for per in (False, True):
    D0 = HeavyHexTFIM(2, p=0.0, q=0.0, periph=per).run(); assert abs(D0[0].sum() - 1) < 1e-12
    D = HeavyHexTFIM(2, p=0.01, periph=per).run(); tot = sum(v.sum() for v in D.values())
    print(f"   periph={per}: noiseless P(k=0)={D0[0].sum():.12f}; noisy total prob {tot:.12f}, P(k=0)={D[0].sum():.4f}"); assert abs(tot - 1) < 1e-10

print("3. image rule for the peripheral check vs explicit ancilla:")
# explicit: qubits 0 = q_i, 1 = q_j, 2 = mediator, 3 = peripheral ancilla of q_i;  data q_i,q_j entangled with reference 4,5
p, q = 0.03, 0.0
D2 = pauli2_sup(depol_rates(p)); SCX = sup(CXM); SCZ = sup(np.diag([1, 1, 1, -1.]))
def explicit():
    N = 6; psi = np.zeros((2,) * 6, complex)
    for a in range(4): psi[a >> 1, a & 1, 0, 0, a >> 1, a & 1] = 0.5
    rho = np.einsum('abcdef,ghijkl->abcdefghijkl', psi, psi.conj())
    rho = apply(rho, sup(H), [3], N)                                         # ancilla |+>
    rho = apply(rho, SCZ, [3, 0], N); rho = apply(rho, D2, [3, 0], N)        # cZ(p -> q_i) + noise
    for qs in ([0, 2], [1, 2], 'rz', [1, 2], [0, 2]):
        if qs == 'rz': rho = apply(rho, rz(PH), [2], N); continue
        rho = apply(rho, SCX, qs, N); rho = apply(rho, D2, qs, N)
    rho = apply(rho, SCZ, [3, 0], N); rho = apply(rho, D2, [3, 0], N)
    rho = apply(rho, sup(H), [3], N)                                         # X-basis readout of the ancilla
    r = rho.reshape(2, 2, 2, 2, 4, 2, 2, 2, 2, 4)
    return {(mf, pf): r[:, :, mf, pf, :, :, :, mf, pf, :].reshape(16, 16) for mf in (0, 1) for pf in (0, 1)}
def image_rule():
    """same bond with hh_sim's channels: 2 data + mediator, peripheral bit tracked classically"""
    sim = HeavyHexTFIM(1, p=p, q=q, periph=True)
    N = 5; psi = np.zeros((2,) * 5, complex)
    for a in range(4): psi[a >> 1, a & 1, 0, a >> 1, a & 1] = 0.5
    br = {0: np.einsum('abcde,fghij->abcdefghij', psi, psi.conj())}
    def noise(br, ch, qs, flip_bit):
        out = {}
        for f, S in ch.items():
            for pb, r in br.items():
                k = pb ^ (f if flip_bit else 0); out[k] = out.get(k, 0) + apply(r, S, qs, N)
        return out
    br = noise(br, sim.cz_ch['cz1'], [0], True)
    for qs, loc in (([0, 2], 'cx1'), ([1, 2], 'cx2'), ('rz', None), ([1, 2], 'cx3'), ([0, 2], 'cx4')):
        if qs == 'rz': br = {k: apply(r, rz(PH), [2], N) for k, r in br.items()}; continue
        br = {k: apply(r, SCX, qs, N) for k, r in br.items()}
        # in the sim the peripheral bit of the CX's data qubit is toggled; here only q_i (qubit 0) carries a check
        ch = sim.cx_ch[loc] if qs[0] == 0 else split_channel(sim._rates(loc), {P: (P, 0) for P in P2[1:]})
        br = noise(br, ch, qs, True)
    br = noise(br, sim.cz_ch['cz2'], [0], True)
    out = {}
    for pf, r in br.items():
        rr = r.reshape(2, 2, 2, 4, 2, 2, 2, 4)
        for mf in (0, 1): out[(mf, pf)] = rr[:, :, mf, :, :, :, mf, :].reshape(16, 16)
    return out
E, I_ = explicit(), image_rule()
for key in sorted(E):
    print(f"   (mediator flag, peripheral flag) = {key}: prob {np.trace(E[key]).real:.6f} vs {np.trace(I_[key]).real:.6f}, max|dChoi| = {np.abs(E[key]-I_[key]).max():.1e}")
    assert np.abs(E[key] - I_[key]).max() < 1e-12
print("all tests passed")
