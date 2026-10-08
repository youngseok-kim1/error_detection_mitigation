"""
Detected-event-driven extrapolation for error-detected circuits: exact toy study.

Circuit   : random Clifford brickwork payload on n data qubits (ring), protected by k two-sided
            coherent Pauli checks (left Pauli P_L before a window of layers, right Pauli
            P_R = U_win P_L U_win^dag after it; windows nested or disjoint so syndromes are
            deterministic).
Noise     : 2q depolarizing(p) after every 2q gate (payload and check gates), readout flip pm
            on each check qubit.
Observable: O = U Z_A U^dag, ideal value exactly +1, measured noiselessly.

stim's detector error model merges all faults with the same (syndrome s, O-flip o) into one
independent mechanism, so the DEM *is* the table of class rates lambda_(s,o).  The exact joint
distribution of (syndrome, O-flip) follows from a Walsh-Hadamard transform, so all estimators
are evaluated at infinite shots (pure bias), and variances come from the delta method.
"""
import numpy as np, stim

C1 = ["I","X","Y","Z","H","H_XY","H_YZ","H_NXY","H_NXZ","H_NYZ","S","S_DAG","SQRT_X","SQRT_X_DAG",
      "SQRT_Y","SQRT_Y_DAG","C_XYZ","C_ZYX","C_NXYZ","C_XNYZ","C_XYNZ","C_NZYX","C_ZNYX","C_ZYNX"]

def rand_pauli(rng, n, w):
    ps = stim.PauliString(n)
    for q in rng.choice(n, size=w, replace=False): ps[int(q)] = int(rng.integers(1, 4))
    return ps

def ctrl_pauli(c, anc, P, p):
    n2 = 0
    for q in range(len(P)):
        g = {1: "CX", 2: "CY", 3: "CZ"}.get(P[q])
        if g:
            c.append(g, [anc, q]); n2 += 1
            if p > 0: c.append("DEPOLARIZE2", [anc, q], p)
    return n2

def build(n, k, reps, p, pm, wL, seed, check_seed=None, windows=None, obs_w=None):
    rng = np.random.default_rng(seed)
    layers = []
    for _ in range(2 * reps):
        par = len(layers) % 2
        layers.append(([C1[i] for i in rng.integers(0, 24, n)], [(q, (q + 1) % n) for q in range(par, n, 2)]))
    M = len(layers)
    def seg(a, b, noise):
        c = stim.Circuit(); n2 = 0
        c.append("I", [n - 1])
        for oneq, pairs in layers[a:b]:
            for q, g in enumerate(oneq): c.append(g, [q])
            for x, y in pairs:
                c.append("CZ", [x, y]); n2 += 1
                if noise > 0: c.append("DEPOLARIZE2", [x, y], noise)
        return c, n2
    tab = lambda a, b: stim.Tableau.from_circuit(seg(a, b, 0)[0])
    windows = windows or [(0, M)] * k
    crng = np.random.default_rng(check_seed if check_seed is not None else seed + 10_000)
    PL = [rand_pauli(crng, n, wL) for _ in range(k)]
    PR = [tab(a, b)(P) for P, (a, b) in zip(PL, windows)]
    orng = np.random.default_rng(seed + 77)
    zs = stim.PauliString(n)
    for q in orng.choice(n, size=obs_w or n // 2, replace=False): zs[int(q)] = 3
    O = tab(0, M)(zs)
    anc = list(range(n, n + k))
    c = stim.Circuit(); n2 = 0
    c.append("R", range(n + k)); c.append("H", anc)
    for t in range(M + 1):
        closing = sorted([j for j in range(k) if windows[j][1] == t], key=lambda j: (-windows[j][0], j))
        opening = sorted([j for j in range(k) if windows[j][0] == t], key=lambda j: (-windows[j][1], -j))
        for j in closing: n2 += ctrl_pauli(c, anc[j], PR[j], p)
        for j in opening: n2 += ctrl_pauli(c, anc[j], PL[j], p)
        if t < M:
            s, m = seg(t, t + 1, p); c += s; n2 += m
    c.append("H", anc)
    if pm > 0: c.append("X_ERROR", anc, pm)
    c.append("M", anc)
    for j in range(k): c.append("DETECTOR", [stim.target_rec(-k + j)])
    tg = []
    for q in range(n):
        if O[q]: tg += [{1: stim.target_x, 2: stim.target_y, 3: stim.target_z}[O[q]](q), stim.target_combiner()]
    c.append("MPP", tg[:-1]); c.append("OBSERVABLE_INCLUDE", [stim.target_rec(-1)], 0)
    return c, n2

def classes(c, k):
    """dict mask -> probability; mask bits 0..k-1 = syndrome, bit k = O flip"""
    out = {}
    for inst in c.detector_error_model().flattened():
        if inst.type != "error": continue
        m = 0
        for t in inst.targets_copy():
            if t.is_relative_detector_id(): m ^= 1 << t.val
            elif t.is_logical_observable_id(): m ^= 1 << k
        pr = inst.args_copy()[0]; q = out.get(m, 0.0); out[m] = q * (1 - pr) + pr * (1 - q)
    return out

def wht(a):
    a = a.copy(); h = 1
    while h < len(a):
        a = a.reshape(-1, 2, h); a = np.stack([a[:, 0] + a[:, 1], a[:, 0] - a[:, 1]], 1).reshape(-1); h *= 2
    return a

PAR = {}
def joint(cls, k, G=1.0):
    """exact P[s + 2^k f]; all class generator rates scaled by G (G != 1: physical amplification)"""
    N = 1 << (k + 1)
    if N not in PAR:
        z = np.arange(N); PAR[N] = np.array([[bin(x & m).count("1") & 1 for x in z] for m in range(N)])
    logc = np.zeros(N)
    for m, pr in cls.items(): logc += np.log(1 - 2 * pr) * G * PAR[N][m]
    return wht(np.exp(logc)) / N

def ps_stats(P, k, mask):
    s = np.arange(1 << k); ok = (s & mask) == 0
    P0, P1 = P[: 1 << k][ok].sum(), P[1 << k:][ok].sum()
    return P0 + P1, (P0 - P1) / (P0 + P1)

lam_of = lambda p: -0.5 * np.log(1 - 2 * p)

def analyse(c, n2, k, p, pm, model_scale=1.3):
    cls = classes(c, k); F = 1 << k; full = F - 1
    Lam = n2 * p + k * pm
    Ld = sum(lam_of(v) for m, v in cls.items() if m & full)
    LdO = sum(lam_of(v) for m, v in cls.items() if (m & full) and (m & F))
    LuO = lam_of(cls.get(F, 0.0)); Lu = Lam - Ld
    P = joint(cls, k)
    _, O_all = ps_stats(P, k, 0); alpha, O_ps = ps_stats(P, k, full)
    r, rO = Lu / Ld, LuO / LdO
    y_all, y_ps = np.log(O_all), np.log(O_ps)
    s = np.arange(1, F); Ps0, Ps1 = P[s], P[s + F]
    Ls = (Ps0 + Ps1) / (P[0] + P[F]); Rs = (Ps0 - Ps1) / (Ps0 + Ps1) / O_ps
    c2 = 0.5 * np.sum(Ls ** 2 * (1 - Rs ** 2))          # in-situ estimate of 2nd-order attenuation
    est = {
        "raw": y_all,
        "ED": y_ps,
        "ED+ext(generic r)": y_ps + r * (y_ps - y_all),
        "ED+ext(O-specific r)": y_ps + rO * (y_ps - y_all),
        "ED+ext(O-spec r)+2nd": y_ps + rO * (y_ps - y_all) + c2,
        "ED+model rescale": y_ps + 2 * LuO,
        "ED+model rescale, scale off": y_ps + 2 * LuO * model_scale,
        "ED+ext(O-spec r), scale off": y_ps + rO * (y_ps - y_all),   # unchanged: only ratios enter
    }
    # delta-method variance x N of the two-point estimator (per-shot influence function)
    def gamma(rr):
        o = np.concatenate([np.ones(F), -np.ones(F)]); A = np.concatenate([np.arange(F) == 0] * 2).astype(float)
        psi = (1 + rr) * A * (o - O_ps) / (alpha * O_ps) - rr * (o - O_all) / O_all
        return float(np.sum(P * psi ** 2) - np.sum(P * psi) ** 2)
    info = dict(n2=n2, Lam=Lam, cov=Ld / Lam, alpha=alpha, fd=LdO / Ld, fu=LuO / Lu, r=r, rO=rO, O_all=O_all, O_ps=O_ps,
                G_pec=np.exp(4 * Lam), G_edpec=np.exp(4 * Lu) / alpha, G_ext=gamma(rO), G_ed=(1 - O_ps ** 2) / alpha / O_ps ** 2)
    return {a: np.exp(b) - 1 for a, b in est.items()}, info, cls, P

if __name__ == "__main__":
    n, k, reps, p, pm, wL = 10, 6, 5, 0.004, 0.004, 2
    c, n2 = build(n, k, reps, p, pm, wL, seed=1)
    b, info, _, _ = analyse(c, n2, k, p, pm)
    print({a: round(float(v), 4) for a, v in info.items()})
    for a, v in b.items(): print(f"  bias[{a:30s}] = {v:+.5f}")
