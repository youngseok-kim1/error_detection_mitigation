# Detected-event extrapolation: exact toy study

**The compiled write-up is `REPORT.md`** (theory, both tests, figures, full code listing). It is rebuilt with `python build_report.py` from `report_template.md`. This README only describes the scripts. The theory of detection-driven extrapolation (fixed vs twirled checks, the mixture law, stratification) is in `theory/THEORY.md` with its own scripts in `theory/`; `theory/reactivity_connection.md` relates it to arXiv:2609.31934 (reactivity functions). The three-stage numerical test (ZNE vs fixed checks vs twirled checks on a deep TFIM circuit) is `demo/DEMO.md`. The hardware-native follow-up (6-ring TFIM on one heavy-hex plaquette with mediator and peripheral checks, p = 0.5–0.3%, detection-based protocols vs ZNE) is `heavyhex/HEAVYHEX.md`.

Exploration for the question: can detected events drive a ZNE-type mitigation of the logical
error left after post-selection (setting of arXiv:2609.13108)?

## Setup
Random Clifford brickwork payload (10 data qubits, 50 CZ) with k two-sided coherent Pauli
checks, 2q depolarizing noise on every 2q gate (payload and checks), readout flips on check
qubits. Observable O = U Z_A U^dag (ideal value 1). stim's detector error model merges faults
with equal (syndrome s, O-flip o) into one mechanism, so it is the table of class rates
Lambda_(s,o); the joint distribution of (syndrome, O-flip) is then exact via a Walsh-Hadamard
transform. All biases are infinite-shot values; `checks.py` cross-checks against sampling.

## Relations used (independent Pauli faults, Clifford circuit, Pauli observable)
    <O>_all = O_0 exp[-2 (L_O^u + L_O^d)]
    <O>_PS  = O_0 exp[-2 L_O^u] (1 - 2 sum_{s!=0} Lambda_(s,0) Lambda_(s,1) + O(lambda^3))
L_O^u, L_O^d: total rate of O-flipping faults that are undetected / detected.
From (syndrome, outcome) data alone, every Lambda_(s,o) with s != 0 is identifiable without
knowing O_0; the single unidentifiable combination is ln O_0 - 2 L_O^u.

Estimators
    two-point extrapolation   O_hat = <O>_PS^(1+r) / <O>_all^r ,   r = L_O^u / L_O^d
    in-situ pair correction   O_hat *= exp[ 1/2 sum_{s!=0} (N_s/N_0)^2 (1 - R_s^2) ],  R_s = <O>_s / <O>_PS
"generic r" uses the coverage ratio (all undetected rate)/(all detected rate) instead of the
O-specific one. "model rescale" multiplies <O>_PS by exp(2 L_O^u) from the noise model;
"scale off" evaluates it with all model rates 30% too large (ratios unchanged).

## Files
- `edzne.py`       circuit builder, class rates, exact joint distribution, estimators
- `sweeps.py`      ensembles over 40 random instances (5 scenarios) + scaling with p
- `checks.py`      Monte Carlo cross-check; physical noise amplification on post-selected data
- `fig_subsets.py` figure: ln<O> vs -ln(alpha) over all 64 check subsets
- `results.txt`    output of sweeps.py and checks.py
- `fig_subset_extrapolation.png`

Requires `stim`, `numpy`, `matplotlib`. Run e.g. `python sweeps.py`.

## Caveats
Clifford circuit with a stabilizer observable is the best case: every fault multiplies <O> by
+-1, so the noise response is a single exponential. Non-Clifford circuits are not tested.
The ratio r is taken from the true model. The in-situ pair term uses the full 2^k syndrome
histogram, which is only practical for small k.

---

# Non-Clifford test (folder `nonclifford/`)

Same footing as above (exact, infinite-shot bias), on a small version of the paper's own
non-Clifford circuit.

## Setup
Trotterized transverse-field Ising ring, 5 data qubits + 5 ancillas (10 qubits), data in |+>.
One step: for each bond, `CX(j,a) CX(k,a) Rz_a(phi) CX(k,a) CX(j,a)`, then `Rx(theta)` on all
data. Checks: each ancilla must read 0 at the end (terminal measurement, no reset). Optional
extra check: parity prod_i X_i = +1, a symmetry of the dynamics, free with X-basis readout.
Noise: Pauli-Lindblad on the two qubits of every CX (total rate p), ancilla readout flip p.
Evaluation: exact density matrix; `T[s, m] = Tr[X^m <s|rho|s>]` gives every syndrome-resolved
X-string expectation in one run. A ring is used because on an open chain a flipped ancilla
(sign flip of one ZZ coupling) is a gauge transformation for X observables and so harmless.

Scenarios
- `flagship.jsonl`: paper angles (theta=3pi/8, phi=-pi/4), 2/4/6 steps, uniform (depolarizing)
  rates, p = 0.004 and 0.012.
- `generic_lo/hi.jsonl`: 16 instances with random angles and log-normal (sigma=1) Pauli rates
  per CX location, 4 steps, p = 0.004 / 0.012.

## Estimators
- `ext, <ratio>`: two-point extrapolation `O_PS + r (O_PS - O_all)` with r from
  coverage (undetected rate / detected rate), from the exact first-order response of the ideal
  circuit (oracle, not available in practice), or fitted on 30 Clifford training circuits
  (same circuit, every angle a random multiple of pi/2; one ratio pooled over sites).
- `CDR on ED/raw`: plain Clifford data regression (one slope) on the same training set.
- `PEC 1st order`: undetected-fault rates set to 0 (exact model) or to -0.3x (model 30% high);
  this is the expectation value of first-order spacetime PEC.
- `u-amp ZNE`: undetected-fault rates scaled by G = 1,2,3, post-selected, extrapolated.
  Inserted faults have zero syndrome, so acceptance is unchanged (checked to 1e-15).
- `ZNE raw` / `ZNE on ED`: all rates scaled by G, without / with post-selection.
- `pair term`: the Clifford formula exp[1/2 sum_s (N_s/N_0)^2 (1 - R_s^2)] applied as is.

## Files
- `ising_ed.py`  simulator, fault classification, post-selection tables
- `study.py`     all estimators for one instance; `run_all.py` drivers; `report.py`, `summary.py` tables
- `test_sim.py`  checks: ideal circuit vs matrix exponentials, classification, O(p^2) scaling, channel
- `fig_nc.py`, `fig_nonclifford.png`, `results_nonclifford.txt`

## Caveats
10 qubits only. Noise sits on CX gates and does not depend on rotation angles, which is what
makes Clifford training circuits transferable. Training values are infinite-shot, and their
shots are not counted in the cost estimates. The pooled ratio assumes sites behave alike.
