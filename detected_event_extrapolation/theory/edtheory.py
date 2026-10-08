"""
Numerical companion to THEORY.md: post-selection on (fixed or flavor-randomized) coherent Pauli
checks as a one-parameter noise family.

Payload : random Clifford brickwork on n data qubits (ring).  Checks: two-sided coherent Pauli
checks, P_L before a window of layers and P_R = U_win P_L U_win^dag after it (nested windows).
Flavors : P_L uniform over the non-identity Paulis on a chosen support (a "twirled" check),
          or one fixed Pauli (a "fixed" check).
Noise   : two-qubit depolarizing p after each payload CZ, p_check after each check gate,
          readout flip pm on the check qubits.
Output  : exact joint distribution P[s, f] of syndrome s and observable flip f (stim DEM +
          Walsh-Hadamard), from which every acceptance rule is evaluated at infinite shots.
"""
import numpy as np, stim, itertools

C1 = ["I","X","Y","Z","H","H_XY","H_YZ","H_NXY","H_NXZ","H_NYZ","S","S_DAG","SQRT_X","SQRT_X_DAG",
      "SQRT_Y","SQRT_Y_DAG","C_XYZ","C_ZYX","C_NXYZ","C_XNYZ","C_XYNZ","C_NZYX","C_ZNYX","C_ZYNX"]

def uniform_pauli(rng, n, support):
    """uniformly random non-identity Pauli supported inside `support`"""
    while True:
        ps = stim.PauliString(n)
        for q in support: ps[int(q)] = int(rng.integers(0, 4))
        if any(ps[q] for q in support): return ps

def payload_layers(n, reps, seed):
    rng = np.random.default_rng(seed); layers = []
    for i in range(2 * reps):
        layers.append(([C1[j] for j in rng.integers(0, 24, n)], [(q, (q + 1) % n) for q in range(i % 2, n, 2)]))
    return layers

def seg(layers, n, a, b, noise):
    c = stim.Circuit(); c.append("I", [n - 1]); n2 = 0
    for oneq, pairs in layers[a:b]:
        for q, g in enumerate(oneq): c.append(g, [q])
        for x, y in pairs:
            c.append("CZ", [x, y]); n2 += 1
            if noise > 0: c.append("DEPOLARIZE2", [x, y], noise)
    return c, n2

def ctrl_pauli(c, anc, P, p):
    n2 = 0
    for q in range(len(P)):
        g = {1: "CX", 2: "CY", 3: "CZ"}.get(P[q])
        if g:
            c.append(g, [anc, q]); n2 += 1
            if p > 0: c.append("DEPOLARIZE2", [anc, q], p)
    return n2

def build(n, layers, PL, windows, p, p_check, pm, obs_seed=77, gadget_p=None):
    gp = gadget_p or (lambda j, side: p_check)
    """checks j = 0..k-1 with left Paulis PL[j] and windows (a_j, b_j); nested or disjoint"""
    k, M = len(PL), len(layers)
    tab = lambda a, b: stim.Tableau.from_circuit(seg(layers, n, a, b, 0)[0])
    PR = [tab(a, b)(P) for P, (a, b) in zip(PL, windows)]
    orng = np.random.default_rng(obs_seed); zs = stim.PauliString(n)
    for q in orng.choice(n, size=n // 2, replace=False): zs[int(q)] = 3
    O = tab(0, M)(zs); anc = list(range(n, n + k))
    c = stim.Circuit(); c.append("R", range(n + k)); c.append("H", anc); n_pay = n_chk = 0
    for t in range(M + 1):
        closing = sorted([j for j in range(k) if windows[j][1] == t], key=lambda j: (-windows[j][0], j))
        opening = sorted([j for j in range(k) if windows[j][0] == t], key=lambda j: (-windows[j][1], -j))
        for j in closing: n_chk += ctrl_pauli(c, anc[j], PR[j], gp(j, 'R'))
        for j in opening: n_chk += ctrl_pauli(c, anc[j], PL[j], gp(j, 'L'))
        if t < M:
            s, m = seg(layers, n, t, t + 1, p); c += s; n_pay += m
    c.append("H", anc)
    if pm > 0: c.append("X_ERROR", anc, pm)
    c.append("M", anc)
    for j in range(k): c.append("DETECTOR", [stim.target_rec(-k + j)])
    tg = []
    for q in range(n):
        if O[q]: tg += [{1: stim.target_x, 2: stim.target_y, 3: stim.target_z}[O[q]](q), stim.target_combiner()]
    c.append("MPP", tg[:-1]); c.append("OBSERVABLE_INCLUDE", [stim.target_rec(-1)], 0)
    return c, n_pay, n_chk

def joint(c, k):
    """exact P[s + 2^k f] from the detector error model"""
    cls = {}
    for inst in c.detector_error_model().flattened():
        if inst.type != "error": continue
        m = 0
        for t in inst.targets_copy():
            if t.is_relative_detector_id(): m ^= 1 << t.val
            elif t.is_logical_observable_id(): m ^= 1 << k
        pr = inst.args_copy()[0]; q = cls.get(m, 0.0); cls[m] = q * (1 - pr) + pr * (1 - q)
    N = 1 << (k + 1); z = np.arange(N); logc = np.zeros(N)
    for m, pr in cls.items(): logc += np.log(1 - 2 * pr) * np.array([bin(x & m).count("1") & 1 for x in z])
    a = np.exp(logc); h = 1
    while h < N:
        a = a.reshape(-1, 2, h); a = np.stack([a[:, 0] + a[:, 1], a[:, 0] - a[:, 1]], 1).reshape(-1); h *= 2
    return a / N

def accept_stats(P, k, accept):
    """accept: boolean array over syndromes -> (alpha, alpha*<O>)"""
    P0, P1 = P[:1 << k][accept].sum(), P[1 << k:][accept].sum()
    return P0 + P1, P0 - P1

def subset_rule(k, T):     # accept iff all checks in T read 0
    return (np.arange(1 << k) & T) == 0
def threshold_rule(k, t):  # accept iff at most t checks read 1
    return np.array([bin(s).count("1") <= t for s in range(1 << k)])

def sigma(w):
    """P[uniform non-identity Pauli on w qubits commutes with a fixed non-identity Pauli there]"""
    return (2 ** (2 * w - 1) - 1) / (4 ** w - 1)

class Ensemble:
    """pooled exact statistics over M flavor draws (or one fixed flavor)"""
    def __init__(self, n, reps, k, supports, windows, p, p_check, pm, M, seed=0, fixed=False, payload_seed=3, gadget_p=None):
        self.n, self.k, self.pm = n, k, pm
        layers = payload_layers(n, reps, payload_seed); rng = np.random.default_rng(seed)
        self.P = []; self.gates = None
        if fixed: PL = [uniform_pauli(rng, n, s) for s in supports]
        for _ in range(M):
            if not fixed: PL = [uniform_pauli(rng, n, s) for s in supports]
            c, n_pay, n_chk = build(n, layers, PL, windows, p, p_check, pm, gadget_p=gadget_p)
            self.gates = (n_pay, n_chk) if self.gates is None else tuple(np.add(self.gates, (n_pay, n_chk)))
            self.P.append(joint(c, k))
        self.gates = (self.gates[0] / M, self.gates[1] / M)
        self.Pm = np.mean(self.P, axis=0)                 # pooled joint distribution
        self.w = [len(s) for s in supports]
        self.b = [sigma(w) * (1 - pm) + (1 - sigma(w)) * pm for w in self.w]   # P[read 0 | visible fault]
        self.a = [1 - pm] * k                                                   # P[read 0 | no visible fault]
    def stats(self, accept):
        al, alO = accept_stats(self.Pm, self.k, accept); return al, alO / al
    def O_all(self): return self.stats(subset_rule(self.k, 0))[1]
    def single_species(self, T):
        """estimate of O_1 assuming every faulty shot is visible to every check in T"""
        al, alO = accept_stats(self.Pm, self.k, subset_rule(self.k, T)); B = np.prod([self.b[j] for j in range(self.k) if T >> j & 1]); A = np.prod([self.a[j] for j in range(self.k) if T >> j & 1])
        # fault-free shots survive with A, visible-faulty with B:  alO = A W0 + B (Oall - W0),  al = A q0 + B (1 - q0)
        return ((alO - B * self.O_all()) / (A - B)) / ((al - B) / (A - B))
    def stratified(self):
        """exact species inversion over all 2^k subsets -> O_1 = E[O | no visible fault]"""
        k = self.k; Ts = range(1 << k)
        Mx = np.array([[np.prod([(self.b[j] if V >> j & 1 else self.a[j]) if T >> j & 1 else 1.0 for j in range(k)]) for V in Ts] for T in Ts])
        al = np.array([accept_stats(self.Pm, k, subset_rule(k, T))[0] for T in Ts]); alO = np.array([accept_stats(self.Pm, k, subset_rule(k, T))[1] for T in Ts])
        q = np.linalg.solve(Mx, al); W = np.linalg.solve(Mx, alO)
        return W[0] / q[0], q, W, np.linalg.cond(Mx)
    def nvar(self, T, O1=None):
        """N x Var of the leak-subtraction estimator  sum_i (A_i - B) o_i / sum_i (A_i - B)  (delta method)"""
        al, alO = accept_stats(self.Pm, self.k, subset_rule(self.k, T)); B = np.prod([self.b[j] for j in range(self.k) if T >> j & 1])
        O1 = self.single_species(T) if O1 is None else O1; Oall = self.O_all(); Ops = alO / al
        Orej = (Oall - alO) / (1 - al) if al < 1 else 0.0
        m2acc = 1 - 2 * O1 * Ops + O1 ** 2; m2rej = 1 - 2 * O1 * Orej + O1 ** 2
        return ((1 - B) ** 2 * al * m2acc + B ** 2 * (1 - al) * m2rej) / (al - B) ** 2
    def u_of(self, accept, T=None):
        """ZNE parameter u = (faulty-shot survival)/alpha for an accept rule (threshold rules: full checks only)"""
        al, _ = accept_stats(self.Pm, self.k, accept); return al
