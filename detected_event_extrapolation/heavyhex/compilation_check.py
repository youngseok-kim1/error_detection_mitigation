"""Is the compilation twirl of the mediated bond gadget (paper App. B, Table hh-coins) a real twirl under depolarizing noise?
Compilation A: m in |0>, CX(i->m) CX(j->m) Rz_m(phi) CX(j->m) CX(i->m), measure m in Z.
Compilation B: data in the H frame, m in |+>, CX(m->i) CX(m->j) Rx_m(phi) CX(m->j) CX(m->i), measure m in X.
Both implement exp(-i phi/2 Z_i Z_j) on the data. Two-qubit depolarizing p after every CX, on (data, m).
We compare the flag-resolved data channels (Choi matrices of 'mediator reads 0' and 'reads 1') in the logical frame."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'demo'))
from tfim_demo import apply, sup, SCX, RZ, RX, SH, pauli2_sup, depol_rates, CXM
p, phi = 0.02, -np.pi / 4
D = pauli2_sup(depol_rates(p))
SCXr = sup(np.array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0.]]))   # CX with control = second qubit
def choi(comp):
    # qubits: 0,1 = data i,j ; 2 = mediator ; 3,4 = reference.  Max-entangle data with reference.
    N = 5; phi_ = np.zeros(16); 
    for a in range(4): phi_[a * 4 + a] = 0.5                  # |00>+|11>+... over (01)(34)
    v = np.einsum('ab,c->abc', phi_.reshape(4, 4), np.array([1., 0]))  # (data01, ref34, m) order -> rearrange below
    psi = np.zeros((2,) * 5, complex)
    for a in range(4):
        i, j = a >> 1, a & 1; psi[i, j, 0, i, j] = 0.5
    rho = np.einsum('abcde,fghij->abcdefghij', psi, psi.conj())
    if comp == 'B':
        rho = apply(rho, SH, [0], N); rho = apply(rho, SH, [1], N); rho = apply(rho, SH, [2], N)   # data -> H frame, m |0> -> |+>
        seq = [([2, 0], SCX), ([2, 1], SCX), ('rot', RX(phi)), ([2, 1], SCX), ([2, 0], SCX)]
    else:
        seq = [([0, 2], SCX), ([1, 2], SCX), ('rot', RZ(phi)), ([1, 2], SCX), ([0, 2], SCX)]
    for qs, S in seq:
        if qs == 'rot': rho = apply(rho, S, [2], N); continue
        rho = apply(rho, S, qs, N)
        dq = qs[0] if qs[0] != 2 else qs[1]
        rho = apply(rho, D, [dq, 2], N)
    if comp == 'B':
        rho = apply(rho, SH, [2], N); rho = apply(rho, SH, [0], N); rho = apply(rho, SH, [1], N)   # measure m in X, data back to Z frame
    r = rho.reshape(2, 2, 2, 4, 2, 2, 2, 4)
    out = []
    for o in (0, 1):
        out.append(r[:, :, o, :, :, :, o, :].reshape(16, 16))
    return out
A, B = choi('A'), choi('B')
for o in (0, 1):
    print(f"mediator flag {o}: trace A {np.trace(A[o]).real:.6f}  B {np.trace(B[o]).real:.6f}   max|Choi_A - Choi_B| = {np.abs(A[o] - B[o]).max():.2e}")
