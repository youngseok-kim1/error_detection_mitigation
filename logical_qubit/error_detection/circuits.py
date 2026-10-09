"""Stage 1 circuits: the two-site mixed-field Ising model (MFIM) on the two logical qubits of one [[4,2,2]] block,
and the same Trotter circuit unencoded on two physical qubits.

H = -J Zbar_0 Zbar_1 - sum_i (g_x Xbar_i + g_z Zbar_i), first-order Trotter, step dt:
    U_step = exp(i J dt Zbar0 Zbar1) exp(i g_z dt Zbar0) exp(i g_z dt Zbar1) exp(i g_x dt Xbar0) exp(i g_x dt Xbar1).

Encoded logical rotations (non-FT, the [R1] fan-in style, one physical rotation on the parity-holding qubit):
  Z part, one fan-in for all three Z-type terms (Zbar0 = Z0Z1, Zbar1 = Z0Z2, Zbar0 Zbar1 = Z1Z2):
     CX(0->1) RZ_1(Zbar0) CX(0->2) RZ_2(Zbar1) CX(1->2) RZ_2(Zbar0Zbar1) CX(1->2) CX(0->2) CX(0->1)       (6 CX)
  X part: Xbar0 = X1X3 via CX(1->3) RX_1 CX(1->3); Xbar1 = X2X3 via CX(2->3) RX_2 CX(2->3)                 (4 CX)
  10 CX per Trotter step, against 2 CX (one RZZ gadget) unencoded.
Syndrome rounds: the 1-FT extraction of iceberg.py (8 CX), ancillas never reset, every `every` Trotter steps
('end' = a single round after the last step, [R1]'s R = 1). Preparation: FT GHZ (Z0Z3-verified) by default.

An op list drives both the trajectory simulator (sim_mc.py) and a stim "skeleton" in which every RZ/RX is dropped:
a logical rotation exp(-i theta P) commutes with the stabilizers, and a Pauli fault E inside a rotation gadget passes
the physical rotation as E R(theta) = R(+-theta) E, so every fault's detector pattern is that of the skeleton
(checked in test_stage1.py). The skeleton gives the detector error model, the fault-count estimator and the
detected/undetected classification of every fault.
"""
import numpy as np, stim

D4, AZ, AX = [0, 1, 2, 3], 4, 5
EXTRACT = [(0, 4), (5, 0), (5, 1), (1, 4), (2, 4), (5, 2), (5, 3), (3, 4)]
MFIM = dict(J=1.0, gx=0.75, gz=0.5, dt=0.5)          # [R1]: J = 1, g_x = 0.75, g_z = 0.5 (melting) / 1.5, dt = 0.5


class Params:
    def __init__(self, p=3e-3, p1=None, q=None, G=1.0):
        self.p2 = p * G; self.p1 = (p / 10 if p1 is None else p1) * G; self.q = p if q is None else q; self.G = G


def _angles(m):
    J, gx, gz, dt = m['J'], m['gx'], m['gz'], m['dt']
    return dict(zz=-2 * J * dt, z=-2 * gz * dt, x=-2 * gx * dt)      # RZ(theta) = exp(-i theta Z / 2)


class OpList:
    """ops: ('R', q) | ('H', q) | ('CX', c, t) | ('RZ', q, th) | ('RX', q, th) | ('N1', q) | ('N2', a, b) |
            ('M', q, tag)  (Z measurement, classical readout flip q, no reset)"""
    def __init__(self, n): self.n, self.ops = n, []
    def h(self, q, nz=True): self.ops.append(('H', q)); nz and self.ops.append(('N1', q))
    def cx(self, c, t): self.ops += [('CX', c, t), ('N2', c, t)]
    def rz(self, q, th): self.ops += [('RZ', q, th), ('N1', q)]
    def rx(self, q, th): self.ops += [('RX', q, th), ('N1', q)]
    def m(self, q, tag): self.ops.append(('M', q, tag))


def encoded_mfim(steps, every=1, model=MFIM, prep="ft"):
    """every: int k (a syndrome round after every k-th step) or 'end' (one round after the last step) or None (no rounds)"""
    a = _angles(model); c = OpList(6)
    for q in range(6): c.ops += [('R', q), ('N1', q)]
    c.h(0); c.cx(0, 1); c.cx(1, 2); c.cx(2, 3)
    if prep == "ft": c.cx(0, 4); c.cx(3, 4); c.m(4, ('prep', -1)); c.ops += [('R', 4), ('N1', 4)]
    rnd = 0
    for s in range(steps):
        c.cx(0, 1); c.rz(1, a['z']); c.cx(0, 2); c.rz(2, a['z']); c.cx(1, 2); c.rz(2, a['zz']); c.cx(1, 2); c.cx(0, 2); c.cx(0, 1)
        c.cx(1, 3); c.rx(1, a['x']); c.cx(1, 3)
        c.cx(2, 3); c.rx(2, a['x']); c.cx(2, 3)
        last = s == steps - 1
        if (every == 'end' and last) or (isinstance(every, int) and (s + 1) % every == 0):
            c.h(AX)
            for g in EXTRACT: c.cx(*g)
            c.h(AX); c.m(AZ, ('z', rnd)); c.m(AX, ('x', rnd)); rnd += 1
    for q in D4: c.m(q, ('d', q))
    c.rounds = rnd; c.prep = prep
    return c


def unencoded_mfim(steps, model=MFIM):
    a = _angles(model); c = OpList(2)
    for q in range(2): c.ops += [('R', q), ('N1', q)]
    for s in range(steps):
        c.cx(0, 1); c.rz(1, a['zz']); c.cx(0, 1)
        c.rz(0, a['z']); c.rz(1, a['z']); c.rx(0, a['x']); c.rx(1, a['x'])
    for q in range(2): c.m(q, ('d', q))
    c.rounds = 0; c.prep = None
    return c


def ideal_logical(steps, model=MFIM):
    """exact <Zbar0>, <Zbar1>, <Zbar0 Zbar1> of the logical Trotter circuit from |00>"""
    a = _angles(model)
    I = np.eye(2); Xm = np.array([[0, 1], [1, 0.]]); Zm = np.diag([1., -1])
    def rot(P, th): return np.cos(th / 2) * np.eye(4) - 1j * np.sin(th / 2) * P
    Z0, Z1, X0, X1 = np.kron(Zm, I), np.kron(I, Zm), np.kron(Xm, I), np.kron(I, Xm)
    U = rot(X1, a['x']) @ rot(X0, a['x']) @ rot(Z0 @ Z1, a['zz']) @ rot(Z1, a['z']) @ rot(Z0, a['z'])
    psi = np.zeros(4, complex); psi[0] = 1
    for _ in range(steps): psi = U @ psi
    ev = lambda O: float(np.real(psi.conj() @ O @ psi))
    return dict(Z0=ev(Z0), Z1=ev(Z1), ZZ=ev(Z0 @ Z1))


# ------------------------------------------------------------------ detectors from a measurement record
def detector_spec(c):
    """list of detectors, each a list of record tags to XOR. No-reset rule D_r = s_r xor s_{r-2}; terminal
    D_T = sigma_T xor sigma_R with sigma_R = s_R xor s_{R-1}; plus the absolute terminal parity ('par')."""
    specs, names = [], []
    if c.prep == "ft": specs.append([('prep', -1)]); names.append(('prep', -1))
    for r in range(c.rounds):
        for g in ('z', 'x'):
            specs.append([(g, r)] + ([(g, r - 2)] if r >= 2 else [])); names.append((g, r))
    term = [('d', q) for q in D4]
    if c.rounds >= 2: term += [('z', c.rounds - 1), ('z', c.rounds - 2)]
    elif c.rounds == 1: term += [('z', 0)]
    specs.append(term); names.append(('T', c.rounds))
    return specs, names

def parity_spec(): return [('d', q) for q in D4]                       # blockwise terminal parity (absolute S_Z)
OBS = {'Z0': [('d', 0), ('d', 1)], 'Z1': [('d', 0), ('d', 2)], 'ZZ': [('d', 1), ('d', 2)]}


# ------------------------------------------------------------------ stim skeleton
PAULI2 = [(a, b) for a in 'IXYZ' for b in 'IXYZ'][1:]
PAULI1 = ['X', 'Y', 'Z']

def skeleton(c, prm, insert=None):
    """stim circuit with the rotations dropped. Noise from prm (DEPOLARIZE after gates, M(q) readout flips) unless
    insert = (noise-op index, Pauli string) is given, in which case the circuit is noiseless apart from that one fault."""
    s = stim.Circuit(); tags = []
    for i, op in enumerate(c.ops):
        k = op[0]
        if k == 'R': s.append("R", [op[1]])
        elif k == 'H': s.append("H", [op[1]])
        elif k == 'CX': s.append("CX", [op[1], op[2]])
        elif k in ('RZ', 'RX'): pass
        elif k == 'N1':
            if insert is None and prm.p1 > 0: s.append("DEPOLARIZE1", [op[1]], prm.p1)
            elif insert is not None and insert[0] == i: s.append("CORRELATED_ERROR", [stim.target_pauli(op[1], insert[1])], 1.0)
        elif k == 'N2':
            if insert is None and prm.p2 > 0: s.append("DEPOLARIZE2", [op[1], op[2]], prm.p2)
            elif insert is not None and insert[0] == i:
                tg = [stim.target_pauli(q, P) for P, q in zip(insert[1], op[1:]) if P != 'I']
                s.append("CORRELATED_ERROR", tg, 1.0)
        elif k == 'M':
            s.append("M", [op[1]], prm.q if insert is None else 0.0); tags.append(op[2])
    pos = {t: j for j, t in enumerate(tags)}; nm = len(tags)
    specs, names = detector_spec(c)
    for sp in specs: s.append("DETECTOR", [stim.target_rec(pos[t] - nm) for t in sp])
    for j, (name, sp) in enumerate(OBS.items()): s.append("OBSERVABLE_INCLUDE", [stim.target_rec(pos[t] - nm) for t in sp], j)
    return s, tags


def _xrun(c, insert):
    """skeleton prepared in logical |++> (H on the data after preparation) and read out in X: observables Xbar0 = X1X3,
    Xbar1 = X2X3 are deterministic; no detectors (the Z-basis run supplies them)."""
    s = stim.Circuit(); tags = []; prepped = False
    first_round = next((i for i, op in enumerate(c.ops) if op[0] == 'CX' and op[1:] == (0, 1)), None)
    # end of preparation = the op before the first Trotter gadget (CX(0,1) after the GHZ chain and verification)
    cx01 = [i for i, op in enumerate(c.ops) if op[0] == 'CX' and op[1:] == (0, 1)]
    prep_end = cx01[1] if len(cx01) > 1 else len(c.ops)
    for i, op in enumerate(c.ops):
        if i == prep_end: s.append("H", D4); prepped = True
        k = op[0]
        if k == 'R': s.append("R", [op[1]])
        elif k == 'H': s.append("H", [op[1]])
        elif k == 'CX': s.append("CX", [op[1], op[2]])
        elif k in ('N1', 'N2') and insert[0] == i:
            P = insert[1]
            tg = [stim.target_pauli(q, L) for L, q in zip(P if k == 'N2' else (P,), op[1:]) if L != 'I']
            s.append("CORRELATED_ERROR", tg, 1.0)
        elif k == 'M':
            if op[2][0] == 'd' and op[2][1] == 0: s.append("H", D4)
            s.append("M", [op[1]]); tags.append(op[2])
    if not prepped: s.append("H", D4)
    pos = {t: j for j, t in enumerate(tags)}; nm = len(tags)
    for j, pair in enumerate(((1, 3), (2, 3))):
        s.append("OBSERVABLE_INCLUDE", [stim.target_rec(pos[('d', q)] - nm) for q in pair], j)
    return s


def classify_faults(c):
    """for every noise op and every Pauli: does the fault fire any detector or the terminal parity (skeleton), and
    if not, is it harmful (flips a Zbar in the |00> run or an Xbar in the |++> run, i.e. its propagated form is a
    non-trivial logical operator)? Returns (undetected, undetected_harmful): dicts noise-op index -> set of Paulis."""
    und, harm = {}, {}
    prm0 = Params(0.0, 0.0, 0.0)
    par_tags = parity_spec()
    for i, op in enumerate(c.ops):
        if op[0] not in ('N1', 'N2'): continue
        Ps = PAULI1 if op[0] == 'N1' else PAULI2
        for P in Ps:
            s, tags = skeleton(c, prm0, insert=(i, P))
            det, obs = s.compile_detector_sampler().sample(1, separate_observables=True)
            meas = s.compile_sampler().sample(1)[0]; pos = {t: j for j, t in enumerate(tags)}
            par = sum(int(meas[pos[t]]) for t in par_tags) % 2          # noiseless parity is 0 (GHZ: even weight)
            if det.any() or par: continue
            und.setdefault(i, set()).add(P)
            if c.n == 6:
                xo = _xrun(c, (i, P)).compile_detector_sampler().sample(1, separate_observables=True)[1]
                if obs.any() or xo.any(): harm.setdefault(i, set()).add(P)
    return und, harm
