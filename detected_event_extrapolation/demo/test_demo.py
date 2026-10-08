import numpy as np, time
from tfim_demo import *
n, d, th, ph = 5, 3, 3 * np.pi / 8, -np.pi / 4
# 1) ideal vs independent statevector (reuse nonclifford/ising_ed.ideal with ring=True, same bond order? ising_ed uses even,odd + closing bond (b,(b+1)%n))
import sys; sys.path.insert(0, "../nc"); from ising_ed import ideal as ideal_sv
iv = ideal_values(n, d, th, ph); sv = ideal_sv(n, d, th, ph, ring=True); svx = np.array([sv[1 << (n - 1 - i)] for i in range(n)])
print("ideal <X_i> demo:", np.round(iv[:n], 6), " statevector:", np.round(svx, 6))
# 2) noiseless checked run: every flag count other than 0 must be empty, branch 0 == plain
M = TFIM(n, d, th, ph, 0.0, 0.0); rng = np.random.default_rng(0); fl = [draw_flavor(rng, n) for _ in range(d)]
br = M.run_checked(fl); print("noiseless: counts present", {c: round(float(np.trace(v.reshape(32, 32)).real), 6) for c, v in br.items()}, " max|branch0-plain|", np.abs(br[0] - M.run_plain()).max())
# 3) detection coin: single fault inserted by hand -> flagged with prob 1/2 over flavors (visible), 0 for Z_jZ_k
def single_fault_flag_prob(pauli_on, M_draw=2000):
    """flag probability of one fault at the start of the ZZ layer of step 0 (noiseless otherwise), over random flavors"""
    rng = np.random.default_rng(1); hits = 0
    for _ in range(M_draw):
        P = draw_flavor(rng, n)
        # bit = anticommutation of P with the fault
        ac = sum(1 for qb in range(n) if P[qb] and pauli_on[qb] and P[qb] != pauli_on[qb]) % 2
        hits += ac
    return hits / M_draw
print("flag prob X0:", single_fault_flag_prob((1, 0, 0, 0, 0)), " Z0:", single_fault_flag_prob((3, 0, 0, 0, 0)), " Y0Y1:", single_fault_flag_prob((2, 2, 0, 0, 0)), " Z0Z1 (invisible):", single_fault_flag_prob((3, 3, 0, 0, 0)), " Z0Z1Z2:", single_fault_flag_prob((3, 3, 3, 0, 0)))
# 4) timing
t = time.time(); M = TFIM(n, 8, th, ph, 0.01, 0.01); br = M.run_checked([draw_flavor(rng, n) for _ in range(8)]); print(f"d=8 checked run: {time.time()-t:.2f}s, total trace {sum(float(np.trace(v.reshape(32,32)).real) for v in br.values()):.6f}")
t = time.time(); M.run_plain(); print(f"d=8 plain run: {time.time()-t:.2f}s")
print("surv matrix d=8 q=0.01 F=4 cond:", np.linalg.cond(surv_matrix(8, 0.01, 4)))
