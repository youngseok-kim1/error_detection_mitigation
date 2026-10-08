"""Option B: full Pauli twirl of the check flavor with the controlled-X repair around every R_z whose bond
anticommutes with the flavor.  Every payload fault, including the rotation-axis ones, is then a fair coin."""
import numpy as np, json, sys, time
from tfim_demo import *; from protocols import *

def anticomm(P, j, k):
    """does Pauli string P (tuple 0..3) anticommute with Z_j Z_k ?  (X or Y on exactly one of j,k)"""
    return ((P[j] in (1, 2)) + (P[k] in (1, 2))) % 2 == 1

class TFIMR(TFIM):
    def zz_layer_repaired(self, rho, P, N):
        r1 = depol_rates(self.p * self.G, {INV1: self.inv_scale}); r2 = depol_rates(self.p * self.G, {INV2: self.inv_scale}); rg = depol_rates(self.gp * self.G)
        CXa = np.zeros((4, 4), complex); CXa[:2, :2] = I2; CXa[2:, 2:] = X; SCXa = sup(CXa)       # controlled-X from the ancilla
        self.n_cv = 0
        for (j, k) in self.bonds:
            rho = apply(rho, SCX, [j, k], N); rho = apply(rho, pauli2_sup(r1), [j, k], N)
            fix = anticomm(P, j, k)
            if fix:
                rho = apply(rho, SCXa, [self.anc, k], N)
                if self.gp > 0: rho = apply(rho, pauli2_sup(rg), [self.anc, k], N)
            rho = apply(rho, RZ(self.phi), [k], N)
            if fix:
                rho = apply(rho, SCXa, [self.anc, k], N)
                if self.gp > 0: rho = apply(rho, pauli2_sup(rg), [self.anc, k], N)
                self.n_cv += 2
            rho = apply(rho, SCX, [j, k], N); rho = apply(rho, pauli2_sup(r2), [j, k], N)
        return rho
    def run_checked(self, flavors):
        n, N = self.n, self.N; rho0 = np.ones((2 ** n, 2 ** n)) / 2 ** n
        br = {0: rho0.reshape((2,) * (2 * n)).astype(complex)}; self.gates_cv = 0; self.gates_cp = 0
        for step in range(self.d):
            P = flavors[step]; new = {}; self.gates_cp += 2 * sum(1 for v in P if v)
            for c, rd in br.items():
                rd = self.rx_layer(rd, n)
                r = np.kron(rd.reshape(2 ** n, 2 ** n), PLUS).reshape((2,) * (2 * N))
                r = self.gadget(r, P, N); r = self.zz_layer_repaired(r, P, N); r = self.gadget(r, P, N)
                rp = apply(r, PROJ_P, [self.anc], N); rm = apply(r, PROJ_M, [self.anc], N)
                tr = lambda x: np.trace(x.reshape(2 ** n, 2, 2 ** n, 2), axis1=1, axis2=3).reshape((2,) * (2 * n))
                b0, b1 = tr(rp), tr(rm); q = self.q
                new[c] = new.get(c, 0) + (1 - q) * b0 + q * b1; new[c + 1] = new.get(c + 1, 0) + q * b0 + (1 - q) * b1
            self.gates_cv += self.n_cv; br = new
        return br

def draw_full(rng, n): return tuple(int(v) for v in rng.integers(0, 4, n))          # uniform over all 4^n Paulis (identity = no check)

def coin_check(n=5, M=4000):
    """flag probability of single faults after CX#1 of bond (0,1) and after CX#2, over the full twirl with repair"""
    rng = np.random.default_rng(0); j, k = 0, 1
    def flag_after_cx1(E):   # image at that location is CX P CX -> flag iff CX E CX anticommutes with P ; repair does not change it
        Q = P2_backprop(E)
        return lambda P: sum(1 for qq in (j, k) if P[qq] and Q[qq] and P[qq] != Q[qq]) % 2
    res = {}
    for name, E in (("Z_k after CX#1", (0, 3)), ("Z_jZ_k after CX#2", (3, 3)), ("X_j after CX#1", (1, 0))):
        hits = 0
        for _ in range(M):
            P = draw_full(rng, n)
            if name.endswith("CX#2"): hits += sum(1 for qq in (j, k) if P[qq] and E[qq - j] and P[qq] != E[qq - j]) % 2
            else: hits += flag_after_cx1(E)(P)
        res[name] = hits / M
    return res
def P2_backprop(E):
    """Clifford back-propagation of a two-qubit Pauli (a,b) on (j,k) through CX(j,k) -> Pauli pair"""
    M = CXM @ kron2(PAULI[E[0]], PAULI[E[1]]) @ CXM
    for a in range(4):
        for b in range(4):
            if abs(abs(np.trace(kron2(PAULI[a], PAULI[b]).conj().T @ M)) - 4) < 1e-9: return (a, b)

if __name__ == "__main__":
    print("coin check (analytic flag rule, full twirl):", coin_check())
    n, TH, PH, p, q = 5, 3 * np.pi / 8, -np.pi / 4, 0.01, 0.01
    depths = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [4, 8, 12]
    N, R, MDRAWS, F = 10_000, 300, int(sys.argv[2]) if len(sys.argv) > 2 else 200, 5; out = {}
    for d in depths:
        t0 = time.time(); rng = np.random.default_rng(d); M = TFIMR(n, d, TH, PH, p, q); pool = []; gcv = gcp = 0
        for _ in range(MDRAWS):
            br = M.run_checked([draw_full(rng, n) for _ in range(d)]); pool.append({c: xdist(v, n) for c, v in br.items()}); gcv += M.gates_cv; gcp += M.gates_cp
        pd = pooled(pool, d); inf, al, qv = twirled_estimates(pd, n, d, q, F); fin = twirled_finite(pool, n, d, q, F, N, R, np.random.default_rng(100 + d))
        I = ideal_values(n, d, TH, PH); mx = lambda v: float(np.mean(v[:n]))
        strat = inf["twirled stratified"]; fs = fin["twirled stratified"]; mxs = fs[:, :n].mean(1)
        out[d] = dict(bias=mx(strat) - mx(I), std=float(mxs.std()), alpha0=float(al[0]), q0=float(qv[0]), gates_cv=gcv / MDRAWS / d, gates_cp=gcp / MDRAWS / d, ps_bias=mx(inf["twirled PS only (t=0)"]) - mx(I))
        print(f"d={d} ({time.time()-t0:.0f}s): gadget gates/step: cP {gcp/MDRAWS/d:.1f} + cX repair {gcv/MDRAWS/d:.1f}  alpha0={al[0]:.3f} q0={qv[0]:.3f} | PS-only bias {out[d]['ps_bias']:+.4f} | stratified bias {out[d]['bias']:+.4f}  10k std {out[d]['std']:.4f}  N Var {1e4*out[d]['std']**2:.1f}", flush=True)
    json.dump(out, open(f"repaired_results_M{MDRAWS}.json", "w"))
