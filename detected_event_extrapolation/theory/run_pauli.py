"""(1) Is the mixture law exact beyond depolarizing noise?  Payload noise: biased two-qubit Pauli channel with
random per-location rates + correlated three-qubit Pauli faults.  (3) Gadget noise and 1% readout error."""
import numpy as np, stim, json
import edtheory as T
from edtheory import *

def make_noise(kind, seed=0):
    """returns noise(circuit, qubits, p) appending a two-qubit Pauli channel; and a 3-qubit correlated fault adder"""
    rng = np.random.default_rng(seed)
    if kind == "depol":
        return lambda c, qs, p: c.append("DEPOLARIZE2", qs, p)
    P2 = [(a, b) for a in range(4) for b in range(4)][1:]          # IX, IY, IZ, XI, ... (15 two-qubit Paulis)
    tgt = {1: stim.target_x, 2: stim.target_y, 3: stim.target_z}
    def biased(c, qs, p):
        """15 independent Pauli faults with random log-normal rates, Z-type faults 4x heavier (exact in the DEM)"""
        w = rng.lognormal(0, 1.0, 15)
        for i, (a, b) in enumerate(P2):
            if a in (0, 3) and b in (0, 3): w[i] *= 4
        w = w / w.sum() * p
        for i, (a, b) in enumerate(P2):
            t = ([tgt[a](qs[0])] if a else []) + ([tgt[b](qs[1])] if b else [])
            c.append("CORRELATED_ERROR", t, w[i])
    return biased

def build_general(n, layers, PL, windows, p, p_check, pm, noise, p3=0.0, seed3=0, gadget_noise=None, obs_seed=77):
    """as edtheory.build but with a general payload channel and optional correlated 3-qubit faults after each layer"""
    rng3 = np.random.default_rng(seed3); k, M = len(PL), len(layers)
    tab = lambda a, b: stim.Tableau.from_circuit(seg(layers, n, a, b, 0)[0])
    PR = [tab(a, b)(P) for P, (a, b) in zip(PL, windows)]
    orng = np.random.default_rng(obs_seed); zs = stim.PauliString(n)
    for q in orng.choice(n, size=n // 2, replace=False): zs[int(q)] = 3
    O = tab(0, M)(zs); anc = list(range(n, n + k)); gn = gadget_noise or noise
    def ctrl(c, a, P, pc):
        for q in range(len(P)):
            g = {1: "CX", 2: "CY", 3: "CZ"}.get(P[q])
            if g:
                c.append(g, [a, q])
                if pc > 0: gn(c, [a, q], pc)
    c = stim.Circuit(); c.append("R", range(n + k)); c.append("H", anc)
    for t in range(M + 1):
        closing = sorted([j for j in range(k) if windows[j][1] == t], key=lambda j: (-windows[j][0], j))
        opening = sorted([j for j in range(k) if windows[j][0] == t], key=lambda j: (-windows[j][1], -j))
        for j in closing: ctrl(c, anc[j], PR[j], p_check)
        for j in opening: ctrl(c, anc[j], PL[j], p_check)
        if t < M:
            oneq, pairs = layers[t]
            for q, g in enumerate(oneq): c.append(g, [q])
            for x, y in pairs:
                c.append("CZ", [x, y]); noise(c, [x, y], p)
            if p3 > 0:                                                   # one random weight-3 correlated Pauli fault per layer
                qs = rng3.choice(n, 3, replace=False); ps = rng3.integers(1, 4, 3)
                c.append("CORRELATED_ERROR", [{1: stim.target_x, 2: stim.target_y, 3: stim.target_z}[int(pp)](int(q)) for q, pp in zip(qs, ps)], p3)
    c.append("H", anc)
    if pm > 0: c.append("X_ERROR", anc, pm)
    c.append("M", anc)
    for j in range(k): c.append("DETECTOR", [stim.target_rec(-k + j)])
    tg = []
    for q in range(n):
        if O[q]: tg += [{1: stim.target_x, 2: stim.target_y, 3: stim.target_z}[O[q]](q), stim.target_combiner()]
    c.append("MPP", tg[:-1]); c.append("OBSERVABLE_INCLUDE", [stim.target_rec(-1)], 0)
    return c

class EnsG(Ensemble):
    def __init__(self, n, reps, k, supports, windows, p, p_check, pm, M, noise, p3=0.0, seed=0, payload_seed=3, gadget_noise=None, pm_assumed=None):
        self.n, self.k, self.pm = n, k, pm
        layers = payload_layers(n, reps, payload_seed); rng = np.random.default_rng(seed); self.P = []
        for i in range(M):
            PL = [uniform_pauli(rng, n, s) for s in supports]
            self.P.append(joint(build_general(n, layers, PL, windows, p, p_check, pm, noise, p3, seed3=7, gadget_noise=gadget_noise), k))
        self.Pm = np.mean(self.P, axis=0); self.w = [len(s) for s in supports]
        q = pm if pm_assumed is None else pm_assumed                    # readout error assumed by the estimator
        self.b = [sigma(w) * (1 - q) + (1 - sigma(w)) * q for w in self.w]; self.a = [1 - q] * k; self.gates = (0, 0)

if __name__ == "__main__":
  n, reps, k = 10, 5, 6; M2 = 2 * reps; FULL = (1 << k) - 1; out = {}
  def line(tag, E):
      O_all = E.O_all(); r = []
      for m in (1, 2, 3, 6):
          T_ = (1 << m) - 1; al, O = E.stats(subset_rule(k, T_)); r.append((m, al, O, E.single_species(T_)))
      print(f"{tag:70s} O_all={O_all:.4f} | " + "  ".join(f"m={m}: PS {O:.4f} est {e:.5f}" for m, al, O, e in r)); out[tag] = r

  print("(1) full-weight twirled checks, noiseless gadgets, M=400 flavor draws, ideal = 1")
  line("depolarizing p=0.004 (reference)", EnsG(n, reps, k, [list(range(n))] * k, [(0, M2)] * k, 0.004, 0, 0, 400, make_noise("depol"), seed=1))
  line("biased random Pauli rates p=0.004", EnsG(n, reps, k, [list(range(n))] * k, [(0, M2)] * k, 0.004, 0, 0, 400, make_noise("biased", 5), seed=1))
  line("biased + correlated 3-qubit faults p3=0.003/layer", EnsG(n, reps, k, [list(range(n))] * k, [(0, M2)] * k, 0.004, 0, 0, 400, make_noise("biased", 5), p3=0.003, seed=1))
  line("biased + correlated, 3x noise (p=0.012, p3=0.009)", EnsG(n, reps, k, [list(range(n))] * k, [(0, M2)] * k, 0.012, 0, 0, 400, make_noise("biased", 5), p3=0.009, seed=1))
  print("\n(3) gadget noise p_check = p = 0.004 and readout error 1%  (full-weight twirled, M=400)")
  line("depol gadgets, readout 0", EnsG(n, reps, k, [list(range(n))] * k, [(0, M2)] * k, 0.004, 0.004, 0.0, 400, make_noise("depol"), seed=1))
  line("depol gadgets, readout 1%", EnsG(n, reps, k, [list(range(n))] * k, [(0, M2)] * k, 0.004, 0.004, 0.01, 400, make_noise("depol"), seed=1))
  line("depol gadgets, readout 1% but estimator assumes 0", EnsG(n, reps, k, [list(range(n))] * k, [(0, M2)] * k, 0.004, 0.004, 0.01, 400, make_noise("depol"), seed=1, pm_assumed=0.0))
  line("biased payload + biased gadgets, readout 1%", EnsG(n, reps, k, [list(range(n))] * k, [(0, M2)] * k, 0.004, 0.004, 0.01, 400, make_noise("biased", 5), seed=1, gadget_noise=make_noise("biased", 9)))
  json.dump(out, open("pauli_results.json", "w"))
