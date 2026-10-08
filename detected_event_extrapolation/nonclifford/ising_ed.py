"""
Non-Clifford test: Trotterized transverse-field Ising chain with ancilla-mediated ZZ rotations
that double as error-detecting checks (small ring version of the circuit in arXiv:2609.13108, Fig. 3).

  step:  for each bond (j,k) with ancilla a:  CX(j,a) CX(k,a) Rz_a(phi) CX(k,a) CX(j,a)
         then Rx(theta) on every data qubit.          data start in |+>, ancillas in |0>.
  checks: every ancilla must read 0 at the end (terminal measurement only, no reset);
          optionally the parity prod_i X_i = +1 (a symmetry of the dynamics, free with X readout).
  noise : Pauli-Lindblad two-qubit depolarizing after every CX (15 Paulis, rate p/15 each),
          readout flip pm on each ancilla.  Rates of individually undetected / detected faults
          can be scaled separately by (gu, gd): gu=0 is what first-order spacetime PEC gives in
          expectation, gu=G>1 is amplification of the undetected sector, gu=gd=G is ordinary
          noise amplification.  Negative gu is allowed (linear map; PEC with a wrong model).

Everything is an exact density-matrix evolution (no shots).  Output of run():
  T[s, m] = Tr[ X^m  <s| rho |s> ]   (s = ancilla syndrome, m = bitmask of an X string on data)
"""
import numpy as np, itertools, functools

I2 = np.eye(2); X = np.array([[0, 1], [1, 0.]]); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1., -1])
PAULI = [I2, X, Y, Z]; XB = [0, 1, 1, 0]; ZB = [0, 0, 1, 1]          # x/z bits of I,X,Y,Z
CXM = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0.]])
P2 = [(a, b) for a in range(4) for b in range(4)]                      # (control Pauli, target Pauli)

def sup(U): return np.kron(U, U.conj())

def apply(rho, S, qs, N):
    m = len(qs); ax = list(qs) + [N + q for q in qs]
    rho = np.moveaxis(rho, ax, range(2 * m)); sh = rho.shape
    return np.moveaxis((S @ rho.reshape(4 ** m, -1)).reshape(sh), range(2 * m), ax)

def anticomm(a, b):      # two-qubit Paulis as (c,t) index pairs
    return (XB[a[0]] * ZB[b[0]] + ZB[a[0]] * XB[b[0]] + XB[a[1]] * ZB[b[1]] + ZB[a[1]] * XB[b[1]]) & 1

@functools.lru_cache(maxsize=None)
def noise_sup(detmask, lams, gu, gd):
    """superoperator of exp(sum_nu g_nu lam_nu (P_nu - 1)) over the 15 two-qubit Paulis (lams[0] unused)"""
    f = np.array([np.exp(-2 * sum(lams[i] * (gd if detmask[i] else gu) * anticomm(P2[i], Q) for i in range(1, 16))) for Q in P2])
    c = np.array([sum(f[q] * (-1) ** anticomm(P2[i], P2[q]) for q in range(16)) for i in range(16)]) / 16
    return sum(c[i] * sup(np.kron(PAULI[P2[i][0]], PAULI[P2[i][1]])) for i in range(16))

def bond_list(n, ring):
    """bonds (j,k) in application order: even bonds, then odd bonds; a ring adds the closing bond"""
    nb = n if ring else n - 1
    return [(b, (b + 1) % n) for b in list(range(0, nb, 2)) + list(range(1, nb, 2))]

class Ising:
    def __init__(self, n, d, theta, phi, parity_check=False, ring=True, noise_seed=None, sigma=1.0):
        self.n, self.d, self.ring = n, d, ring
        self.bonds = bond_list(n, ring)
        self.k = len(self.bonds); self.N = n + self.k; self.parity_check = parity_check
        self.thetas = np.broadcast_to(np.asarray(theta, float), (d * n,)) if np.ndim(theta) == 0 else np.asarray(theta, float)
        self.phis = np.broadcast_to(np.asarray(phi, float), (d * self.k,)) if np.ndim(phi) == 0 else np.asarray(phi, float)
        ops = []
        for _ in range(d):
            for a, (j, kq) in enumerate(self.bonds, start=n):
                ops += [("cx", j, a), ("cx", kq, a), ("rz", a), ("cx", kq, a), ("cx", j, a)]
            ops += [("rx", q) for q in range(n)]
        self.ops = ops
        self.cx = [(i, o[1], o[2]) for i, o in enumerate(ops) if o[0] == "cx"]
        # per-location Pauli weights (sum to 1 at each location): uniform = depolarizing,
        # noise_seed -> log-normal spread of the 15 rates, fixed per CX location
        nl = len(self.cx)
        if noise_seed is None: self.w = np.full((nl, 16), 1 / 15)
        else:
            w = np.random.default_rng(noise_seed).lognormal(0, sigma, (nl, 16)); w[:, 0] = 0
            self.w = w / w.sum(1, keepdims=True)
        self.w[:, 0] = 0
        self.classify()

    def classify(self):
        """syndrome of each of the 15 Paulis after each CX, by propagation through the CX skeleton"""
        n, N = self.n, self.N; self.synd = {}; self.detmask = {}
        for g, (pos, c, t) in enumerate(self.cx):
            masks = [False]
            for (pc, pt) in P2[1:]:
                x = np.zeros(N, int); z = np.zeros(N, int)
                x[c], z[c], x[t], z[t] = XB[pc], ZB[pc], XB[pt], ZB[pt]
                for (pos2, c2, t2) in self.cx[g + 1:]:
                    x[t2] ^= x[c2]; z[c2] ^= z[t2]
                s = tuple(x[n:]); par = int(z[:n].sum() & 1)
                self.synd[(g, pc, pt)] = (s, par)
                masks.append(any(s) or (self.parity_check and par == 1))
            self.detmask[g] = tuple(masks)
        self.n_det = sum(sum(m) for m in self.detmask.values()); self.n_tot = 15 * len(self.cx)
        self.w_det = float(sum(self.w[g][np.array(self.detmask[g])].sum() for g in self.detmask)); self.w_tot = float(len(self.cx))

    def run(self, p, pm, gu=1.0, gd=1.0):
        n, N, k = self.n, self.N, self.k
        psi = np.zeros(2 ** N, complex)
        plus = np.ones(2 ** n) / np.sqrt(2 ** n); zero = np.zeros(2 ** k); zero[0] = 1
        psi = np.kron(plus, zero)
        rho = np.outer(psi, psi.conj()).reshape((2,) * (2 * N))
        RZ = lambda f: sup(np.diag([np.exp(-1j * f / 2), np.exp(1j * f / 2)]))
        RX = lambda t: sup(np.cos(t / 2) * I2 - 1j * np.sin(t / 2) * X)
        SCX = sup(CXM); g = 0; iz = 0; ix = 0
        for o in self.ops:
            if o[0] == "cx":
                S = noise_sup(self.detmask[g], tuple(p * self.w[g]), gu, gd) @ SCX if p > 0 else SCX
                rho = apply(rho, S, [o[1], o[2]], N); g += 1
            elif o[0] == "rz": rho = apply(rho, RZ(self.phis[iz]), [o[1]], N); iz += 1
            else: rho = apply(rho, RX(self.thetas[ix]), [o[1]], N); ix += 1
        # ancilla-diagonal blocks -> T[s, m] = Tr[X^m rho_s]
        r = rho.reshape(2 ** n, 2 ** k, 2 ** n, 2 ** k)
        blocks = np.einsum("asbs->sab", r)
        idx = np.arange(2 ** n)
        T = np.array([[blocks[s][idx, idx ^ m].sum().real for m in range(2 ** n)] for s in range(2 ** k)])
        # ancilla readout flips (detected faults): convolve over syndrome bits
        if pm > 0:
            q = (1 - np.exp(2 * gd * 0.5 * np.log(1 - 2 * pm))) / 2
            T = T.reshape((2,) * k + (2 ** n,))
            for b in range(k): T = (1 - q) * T + q * np.flip(T, axis=b)
            T = T.reshape(2 ** k, 2 ** n)
        return T

    def rates(self, p, pm):
        """total generator weight of detected / undetected faults"""
        Ld = self.w_det * p + self.k * (-0.5 * np.log(1 - 2 * pm) if pm > 0 else 0)
        return Ld, (self.w_tot - self.w_det) * p

def table(T, n, parity_check):
    """-> P[sigma], E[sigma, m]: unnormalised weight and X-string value per extended syndrome"""
    full = 2 ** n - 1; m = np.arange(2 ** n)
    if not parity_check: return T[:, 0].copy(), T.copy()
    Pp, Pm = (T[:, 0] + T[:, full]) / 2, (T[:, 0] - T[:, full]) / 2
    Ep, Em = (T + T[:, m ^ full]) / 2, (T - T[:, m ^ full]) / 2
    return np.concatenate([Pp, Pm]), np.concatenate([Ep, Em])       # sigma = s + 2^k * parity_bit

def post(T, n, parity_check):
    """(alpha, <X^m> post-selected, <X^m> raw) for all X strings m"""
    P, E = table(T, n, parity_check)
    return P[0] / P.sum(), E[0] / P[0], E.sum(0) / P.sum()

def ideal(n, d, theta, phi, ring=True):
    """noiseless <X^m> for all X strings m, from a data-only statevector (independent of run())"""
    thetas = np.broadcast_to(np.asarray(theta, float), (d * n,)) if np.ndim(theta) == 0 else np.asarray(theta, float)
    bonds = bond_list(n, ring)
    phis = np.broadcast_to(np.asarray(phi, float), (d * len(bonds),)) if np.ndim(phi) == 0 else np.asarray(phi, float)
    b = np.arange(2 ** n); bit = lambda q: (b >> (n - 1 - q)) & 1
    psi = np.ones(2 ** n, complex) / np.sqrt(2 ** n); iz = ix = 0
    for _ in range(d):
        for (j, kq) in bonds:
            psi = psi * np.exp(-1j * phis[iz] / 2 * (1 - 2 * (bit(j) ^ bit(kq)))); iz += 1
        for q in range(n):
            t = thetas[ix]; ix += 1
            psi = np.cos(t / 2) * psi - 1j * np.sin(t / 2) * psi[b ^ (1 << (n - 1 - q))]
    return np.array([np.real(np.vdot(psi, psi[b ^ m])) for m in range(2 ** n)])
