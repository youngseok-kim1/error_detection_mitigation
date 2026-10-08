# error_detection — detection-conditioned extrapolation for encoded (logical-qubit) simulations

Project folder: `logical_qubit/error_detection/`. Started 8 October 2026. Companion to `../../detected_event_extrapolation/` (physical-qubit Pauli-check study), whose theory notes this document builds on.

## 1. Where this comes from

The physical-qubit study (`detected_event_extrapolation/`) asked whether error-detection events can drive a ZNE-type extrapolation. Its results, in short:

- Twirling the flavor of a coherent Pauli check turns detection into a fair coin and makes partial post-selection an exact mixture law; the fault-free value follows by linear extrapolation in the accepted fraction, model-free, to all orders for arbitrary stochastic Pauli noise [N1].
- On a deep non-Clifford TFIM circuit this beats exponential ZNE in accuracy at every depth, but only beats it in cost while the checked circuit carries about one expected fault, $(g_p+g)\,p\,d\lesssim1.3$; the cost exponents are all proportional to $pd$, so there is no threshold in the gate error alone, and the depth-independent comparison favours checks only if the gadget overhead ratio $r$ and the observable sensitivity $\kappa$ satisfy $1+r<G_{\max}\kappa$ [N2, N3].
- On heavy-hex hardware with ancilla-mediated bonds the flavor twirl is not available at all (every valid check image at a CZ noise location is $Z$-type), $Z$-type data faults are intrinsically invisible, and the check gadgets cost as much as they remove [N4, N5]. **Conclusion: for ZNE-style mitigation of physical circuits the approach brings no immediate advantage.**

The natural home of the same logic is a quantum error-detecting (QED) code: the detectors exist anyway, the circuit is already paying for them, and the problem is the opposite one — full post-selection on trivial syndromes leaves no shots once the circuit is large, and the detector record is thrown away together with the shots. The reference experiment is Froland et al. [R1]: 21 blocks of the $[[4,2,2]]$ Iceberg code on ibm_boston (42 logical qubits, up to 136 physical), Trotterized mixed-field Ising dynamics in 1+1D and 2+1D, non-fault-tolerant logical rotations and fault-tolerant syndrome extraction, and a selective-filtering method, Observable-Ranked Postselection (ORP), to avoid the exponential shot loss.

## 2. The reference paper: Froland et al., arXiv:2607.24947 (27 July 2026)

*Realizing Error Suppression in Partially Fault-Tolerant Quantum Simulations with IBM Quantum Computers*, IQuS / University of Washington. 42 pages; summary from the v1 PDF.

### 2.1 Claim in one paragraph

Iceberg-code error detection improves local-observable estimates of Ising dynamics on ibm_boston, with no noise-learning mitigation anywhere:
- **Scale:** 21 blocks of $[[4,2,2]]$, i.e. 42 logical qubits on up to 136 physical qubits.
- **Gadget mix:** syndrome extraction is fault-tolerant; logical rotations are not ("partially FT", $p_L=Ap+Bp^2$ with $A\ll B\ll1$ by design).
- **Shot loss:** full syndrome post-selection loses shots exponentially. Observable-Ranked Postselection (ORP) post-selects only on the detectors that correlate with the observable.
- **1+1D (42-site chain):** the encoded circuits are 2–6× deeper. They still win by 2–6% in the signal-survival factor α at intermediate times, and by up to 36% in median absolute error for $t\le4$. The melting regime loses beyond $t\approx4.5$.
- **2+1D (logical square grid):** the encoding removes the heavy-hex embedding overhead, and the unencoded circuit is 1.5× deeper. The improvement in α grows roughly linearly in time (≈23% per unit $t$) to over 200% at late times. Median absolute error improves by 57%.
- **Fairness caveat:** encoded and unencoded runs use equal shot counts. The authors note that the unencoded "extra" shots could be spent on noise-learning mitigation (ZNE, PEC, …), and they leave that comparison open. **That comparison is the question this project asks.**

### 2.2 Setup

| item | detail |
|---|---|
| code | $[[4,2,2]]$: $S_X=X^{\otimes4}$, $S_Z=Z^{\otimes4}$; $\bar X_i=X_{i+1}X_3$, $\bar Z_i=Z_0Z_{i+1}$; logical $\lvert00\rangle$ = GHZ(4); undetectable errors contain two $X$ or two $Z$ |
| block on heavy-hex | 4 data $d_0..d_3$ + ancillas $a_x,a_z$ (one degree-3 vertex + two degree-2 vertices); "bridge" qubits between some blocks |
| configurations | "compute" (logical gates) and "syndrome" (extraction); switching costs depth ≈15 each way |
| syndrome extraction | measures $S_X,S_Z$ simultaneously, $a_x,a_z$ act as mutual flags against hook errors, verified 1-FT by simulation |
| state prep | non-FT GHZ (CNOT depth 4); the FT version (depth 7 + 2 measurements) is worse at their operating point |
| logical rotations | non-FT: CNOT fan-in of the support parity onto one central qubit, one physical $R_z$, uncompute. Single-logical-qubit rotations and intra-block $R_{ZZ}$ are folded into the inter-block fan-ins at no extra two-qubit depth. A frame bit tracks SWAP$(d_1,d_2)$; odd steps are reversed to restore frames |
| model | MFIM $H=-J\sum\bar Z_i\bar Z_j-\sum(g_x\bar X_i+g_z\bar Z_i)$, $J=1$, $g_x=0.75$, $g_z=0.5$ (bubble melting) / $1.5$ (localized); false-vacuum bubble quench; first-order Trotter, $\delta t=0.5$, up to 16 steps ($t=8$) |
| observable | $\langle\bar Z_n\rangle$ per logical qubit; figure of merit α = slope of measured vs MPS values through the origin, plus residual σ |
| schedule | $R$ = number of syndrome rounds ($\le1$ per Trotter step); $R=1$ means a single round at the end |
| reset | ancillas **not** reset: outcomes are running parities, detectors $v_{j+1}=s_{j-1}\oplus s_{j+1}$ (two rounds apart), first detector padded with 0, terminal $Z$ readout of the data closes the $S_Z$ detectors |

### 2.3 ORP and its theory (Methods A, App. A)

- $\Delta_j=\langle O\rangle_{v_j=0}-\langle O\rangle_{v_j=1}$. Rank detectors by $|\Delta_j|$. Cut $k$ keeps the shots in which none of the top-$k$ detectors fired, and blockwise terminal parity is always enforced.
- Plateau finder on a held-out half of the data ($f_{\min}=0.005$, $f_{\rm ref}=0.05$, $z=1$) picks $k^\ast$, which is applied to the other half. If there is no plateau, $k^\ast=\arg\min[(\langle O\rangle_k-O_{\rm ref})^2+\sigma_k^2]$. Example: $\langle Z\rangle=0.85(2)$ against an ideal 0.88.
- Mechanism model: independent mechanisms $m$ with rates $p_m$ and flip probabilities $q_m$, each flipping a set of detectors (single-detector and pair processes). Exact solution for $E[Q\mid v_j]$. At weak noise $\Delta_j\approx2\langle O\rangle B_j\sum_{m\ni j}p_mq_m/\sum_{m\ni j}p_m$, the flip probability of the processes that fire $v_j$, where $B_j$ is the attenuation from everything else.
- Alternative metric $\beta_j$ (linear regression of $Q$ on all detectors). It differs from $\Delta_j$ only through pair processes and gives comparable results on hardware.
- The detectors ORP ranks highest coincide with the backward light cone of $O$ (checked by Pauli propagation, $\gamma_j$, at $p=0.003$ depolarizing).

### 2.4 Fault-tolerance benchmarks and ablations (App. B, C)

- **Leading-order logical error.** Non-FT components give $Ap$ and FT components $Bp^2$, with a crossover at $p^\ast\approx A/B$.
- **Memory experiment on hardware.** 60 rounds (761 CNOT depth). Dynamical decoupling and no-reset lower the TVD by more than 10×, and acceptance decays exponentially. Using only the last round amplifies the error, and the undetected fraction grows quadratically in the number of rounds (rounds 8–30): their evidence that the gadget suppresses errors at $O(p^2)$.
- **Ancilla-assisted rotations** (extra ancilla prepared in $\lvert0\rangle$ and post-selected) have a 2–3× smaller $A$ for intra-block $R_{ZZ}$ and are no better for inter-block ones. The standard gadget was kept.
- **Frequency of syndrome extraction.** More rounds help early and hurt late. Without resets, a fired ancilla stays in $\lvert1\rangle$, and every later rotation routed through it runs at $-\theta$, so $J_{\rm eff}=J(1-2p_{\rm flip})$. Adding resets does not help because they cost coherence time.
- **Where the gain comes from.** Without post-selection the 2D gain mostly disappears (α drops by up to 0.4), so the improvement comes from post-selection, not from the shallower circuit.
- **Hardware sensitivity.** A device with somewhat worse two-qubit error ($O(10^{-3})$) and coherence ($O(100\,\mu s)$) shows no encoding advantage at all.

### 2.5 How it maps onto the physical-qubit study

| physical-qubit study (`detected_event_extrapolation/`) | this paper |
|---|---|
| coherent Pauli check / mediator flag | stabilizer detector $v_j$ |
| post-selection = projection onto "no visible fault" | full syndrome post-selection |
| invisible sector: faults along the rotation axis (logical operators of the check) | the $Ap$ term: a fault on the parity qubit between fan-in and uncompute is a weight-2 logical error |
| heavy-hex: logical-$Z$ data faults invisible to every valid local check | same structure one level up: errors equivalent to a logical operator are never detected |
| mediator without per-step reset flips later rotation signs | no-reset ancilla gives $J_{\rm eff}=J(1-2p_{\rm flip})$ |
| detection extrapolation with a structural ratio $r$ fails when detected and undetected faults affect $O$ differently | ORP's plateau is $O_0$ (fully post-selected), not $O_{\rm ideal}$; same missing second step |
| post-selection + ZNE on the remainder: unbiased at every point, 3–9× cheaper than robust ZNE | not tried in the paper |

### 2.6 How the reference experiment evolves the Ising model with logical qubits

- Code: $S_X=X_0X_1X_2X_3$, $S_Z=Z_0Z_1Z_2Z_3$, logical $|00\rangle=$ GHZ(4), logical operators $\bar X_i=X_{i+1}X_3$, $\bar Z_i=Z_0Z_{i+1}$ ($i=0,1$). Every logical Pauli is a weight-2 physical Pauli; $\bar Z_0\bar Z_1=Z_1Z_2$; inter-block $\bar Z\bar Z$ is weight 4.
- Hamiltonian $H=-J\sum\bar Z_i\bar Z_j-\sum(g_x\bar X_i+g_z\bar Z_i)$, first-order Trotter; $J=1$, $g_x=0.75$, $g_z=0.5$ (bubble melting) or $1.5$ (localized).
- Non-Clifford logical gates are physical multi-qubit Pauli rotations done non-fault-tolerantly: a CNOT fan-in collects the parity of the support onto one central qubit (a data or bridge qubit, routed through the idle ancilla when needed), one physical $R_z$ is applied there, and the fan-in is uncomputed (depth 4–12 CNOTs per inter-block $R_{ZZ}$). The single-logical-qubit rotations and the intra-block $R_{ZZ}$ are folded into the same fan-in chains at the cost of single-qubit gates only. A per-block frame bit tracks the SWAP$(d_1,d_2)$ needed for $\bar Z_1,\bar X_0$. One fault on the parity qubit between fan-in and uncompute becomes a weight-2, undetectable logical error: the $Ap$ term of $p_L=Ap+Bp^2$. In our language, a fault along the logical rotation axis is a logical operator and hence invisible to the stabilizers — the rotation-axis sector of [N1] one level up.
- Only syndrome extraction is fault-tolerant (1-FT circuit, the two ancillas $a_x,a_z$ as mutual flags), run in a separate "syndrome" qubit configuration; switching costs depth $\approx15$ each way, so syndromes are extracted every $R$ Trotter steps. Ancillas are not reset; detectors are parities of syndrome bits two rounds apart, and the terminal data readout closes the last detector.
- Encoded circuits are 2–6× deeper than unencoded; local observables improve by 2–6% (1+1D) and over 200% at late times (2+1D), with no noise-learning mitigation.

## 3. ORP and what we add

ORP [R1, Methods A, App. A]: for each detector $v_j$ compute from data $\Delta_j=E[Q\mid v_j{=}0]-E[Q\mid v_j{=}1]$, rank detectors by $\Delta_j$, keep the shots in which none of the top-$k$ detectors fired, choose $k^\ast$ by a plateau-finding heuristic on held-out shots; blockwise terminal parity is always enforced. At weak noise $\Delta_j\approx2\langle O\rangle\sum_{m\ni j}p_mq_m/\sum_{m\ni j}p_m$, the flip probability of the processes firing $v_j$; the ranked detectors coincide with the backward light cone of $O$. An additive regression metric $\beta_j$ is defined and used only as an alternative ranking.

**Geometric law** (theory in [N6]; paper Sec. 10). In a Clifford code circuit with independent Pauli faults, let $k$ be the number of detected faults in a shot and $o=\pm1$ the logical outcome. Then, exactly,

$$E[o\mid k]=O_0\,\bar\rho^{\,k},\qquad\bar\rho=1-2\kappa_d ,$$

with $O_0$ the fully post-selected value and $\kappa_d$ the probability that a detected fault flips $O$ (conditioned on $k$, Poisson faults are $k$ iid draws from the detected-fault distribution, independent of the undetected ones, and each contributes a sign). Detector by detector, in the sparse regime ORP also assumes,

$$E[Q\mid v]=O_0\prod_j\rho_j^{\,v_j},\qquad\rho_j=1-\frac{\Delta_j}{E[Q\mid v_j{=}0]}=1-2\kappa_j ,$$

whose first-order expansion is ORP's additive model. Consequences:

1. **Soft ORP / regression estimator.** Fit $\rho_j$ from the detector-wise conditional means (which ORP already computes), then $\hat O_0=\sum_iw_iQ_i/\sum_iw_i^2$ with $w_i=\prod_j\rho_j^{v_{ij}}$ over all shots. Detectors outside the light cone have $\rho_j=1$ and cost nothing; light-cone detectors down-weight a shot by $\rho_j^2$ instead of discarding it; no cut level, no plateau search.
2. **Effective shots** $N\exp[-\sum_jp_j(1-\rho_j^2)]$ against $N\exp[-\sum_{j\le k^\ast}p_j]$ for the cut — a gain of $\exp[\sum_jp_j\rho_j^2]$. With $\bar\rho$ known, $N_{\rm eff}=N\alpha_0^{1-\bar\rho^2}$ against $N\alpha_0$ for full post-selection (e.g. $N=10^6$, $\alpha_0=10^{-3}$, $\bar\rho=0.9$: $2.7\times10^5$ against $10^3$).
3. **Built-in model check.** $\ln E[Q\mid |v|{=}k]$ over light-cone detectors must be linear in $k$; curvature flags overlapping faults or the breakdown of the product structure under non-Clifford rotations.
4. **What stays.** The plateau "set by undetectable errors" in ORP is $O_0$; neither method reaches $O_{\rm ideal}$ without a second step: (a) ZNE from below with the detector as the noise meter ($\Lambda_d\to G\Lambda_d$ measured from the detection rate, $\Lambda_u\propto G^2$ for a distance-2 code, one-parameter extrapolation in $G^2$); (b) model-assisted $\Lambda_u=c\Lambda_d^2$ with $c$ a circuit property; (c) rates learned from the detector record and first-order PEC on the undetected generator, as in the spacetime-mitigation paper [R2].
5. **QEC extension** ($d\ge3$): the decoder's soft output $g$ (complementary gap) gives $P_L(g)$ per shot and $E[o\mid g]=O_{\rm ideal}(1-2P_L(g))$ — a regression over all shots, no post-selection, effective count $N\,E[(1-2P_L)^2]$.

Assumptions, stated: independent Pauli faults (twirled gates); Clifford circuit for exactness, first order in light-cone overlap for the non-Clifford rotations of [R1]; a detector count that resolves the fault count (a single terminal round of two stabilizers hides $k$ and contaminates the $s=0$ bin by even-$k$ cancellations $\approx\Lambda_d^2/2$ per stabilizer — periodic checks are needed for this reason beyond detection itself); a separate treatment of the undetected sector.

## 4. Plan

### 4.1 Question

At equal shot budget, counted as shots executed before any post-selection, does an encoded $[[4,2,2]]$ simulation plus detection-based post-processing estimate a logical observable more accurately than the **unencoded circuit plus conventional ZNE**? [R1] compared only against raw unencoded results. The physical-qubit study answered the analogous question on heavy-hex: post-selection followed by ZNE on the remainder won, and detection extrapolation with a structural ratio failed where detected and undetected faults affect $O$ differently. The logical setting changes two things:
- **The detectors are free.** The code needs them anyway, and $g_{\rm check}$ is the syndrome-round overhead, not an added gadget.
- **The undetected sector has a known structure.** It is the $Ap$ term (faults on the parity qubit of non-FT rotations) plus $Bp^2$ (pairs).

### 4.2 Protocols compared (same set at every stage)

| tag | circuit | post-processing | needs |
|---|---|---|---|
| U-raw, U-ZNE exp / cum2 | unencoded logical circuit on physical qubits (same Trotter order, CX–$R_z$–CX for $R_{ZZ}$) | none / exponential (G = 1, 2) / 2nd cumulant (G = 1, 2, 3) | amplified copies |
| E-raw | encoded | none (terminal blockwise parity only) | – |
| E-PS | encoded | full post-selection on every detector | – |
| E-ORP | encoded | [R1]'s ranking + plateau finder, reimplemented exactly (held-out half, $f_{\min}$, $f_{\rm ref}$, $z$) | – |
| E-soft | encoded | geometric weights $w_i=\prod_j\rho_j^{v_{ij}}$ (§3, item 1) | – |
| E-ext | encoded | two-point detection extrapolation $O_{\rm PS}(O_{\rm PS}/O_{\rm all})^r$, with $r$ from structure or fitted on Clifford training circuits | structure / training |
| E-PS+ZNE, E-ORP+ZNE | encoded, amplified | ZNE on the post-selected data (G = 1, 2, 3); fit form below | amplified copies |
| E-PS+PEC (reference) | encoded | undetected-sector rates set to zero (exact model), variance × $e^{4\Lambda_u}$ | full model |

**Metric.** Infinite-shot bias, standard deviation at $N=10^4$, and $N_{\rm req}$ for RMSE ≤ ε (ε = 0.02, also 0.05), as in `../../detected_event_extrapolation/heavyhex/`.

**ZNE fit form for the encoded remainder.** $\Lambda_u(G)=aG+bG^2$: the $a$ term comes from the non-FT rotations, the $b$ term from detected-pair cancellations and FT-gadget pairs. A one-parameter $G^2$ extrapolation (§3, item 4(a)) is biased unless $a\approx0$, so the fits to use are exponential (G = 1, 2), 2nd cumulant (G = 1, 2, 3), and the two-term $\ln O=c_0-2(aG+bG^2)$. The detection rate measured at each gain serves as the noise meter: $\Lambda_d(G)/\Lambda_d(1)$ gives the realized gain, which matters when folding is inexact on hardware.

### 4.3 Stages, simplest first

**Stage 0: one block, memory, Clifford. Goal: infrastructure and detector bookkeeping.** Tool: stim. **Status: done, see `ED_stage0.md`.**
- Build the block circuits:
  - non-FT GHZ preparation;
  - the flagged simultaneous $S_X/S_Z$ extraction with $a_x,a_z$ as mutual flags (reconstruct from [R1] Fig. 1c and Ref. [75] of [R1]);
  - compute ↔ syndrome switching (SWAPs);
  - terminal $Z$ readout.
- Build the detectors for no-reset (two rounds apart, 0-padded) and reset (consecutive) ancillas.
- Checks:
  - (i) every single fault is detected or harmless (exhaustive enumeration, the 1-FT condition);
  - (ii) the post-selected logical error is $O(p^2)$, reproducing App. B Fig. 11d–e;
  - (iii) "last round only" amplification gives the quadratic growth of App. B Fig. 12c;
  - (iv) the geometric law $E[o\mid k]=O_0\bar\rho^{\,k}$ holds exactly from the detector error model.
- Exit criterion: all four reproduced.

**Stage 1: one block, two logical qubits, non-Clifford. Goal: the first real comparison.** Tool: exact density matrix, 4 data + 2 ancillas = 6 qubits, branched by detector pattern (the full $2^{2R+2}$ record is exact for $R\le4$).
- **Logical circuit:** two-site MFIM on $(\bar 0,\bar1)$. $R_{ZZ}$ on $\bar Z_0\bar Z_1=Z_1Z_2$ by CX fan-in, $R_z$, uncompute. $\bar X_i$ and $\bar Z_i$ rotations as weight-2 fan-ins, with the frame bit for SWAP$(d_1,d_2)$ as in [R1]. Syndrome round every $R$ Trotter steps, ancillas without reset. The unencoded baseline is the same Trotter circuit on 2 physical qubits.
- **Observables:** $\langle\bar Z_0\rangle$, $\langle\bar X_0\rangle$, $\langle\bar Z_0\bar Z_1\rangle$, so that the rotation-axis errors are sometimes harmless (Z-type observables) and sometimes harmful.
- **Noise:** two-qubit depolarizing $p\in\{0.1,0.2,0.3,0.5\}\%$, single-qubit $p/10$, readout $q=p$, perfect reset; idle noise in a second pass.
- **Depth:** 1–12 Trotter steps; $R\in\{1,2,\text{every step}\}$.
- **Questions:**
  - Q1: Is the undetected sector exactly the parity-qubit rotation-axis faults plus pairs? Use single-fault classification as in `../../detected_event_extrapolation/heavyhex/classify.py`.
  - Q2: Does $\ln E[Q\mid k]$ stay linear in $k$ at non-Clifford angles, and how large is the curvature?
  - Q3: Which protocol has the lowest $N_{\rm req}$, and where does the encoded route overtake U-ZNE (break-even in $p$ and depth)?
  - Q4: Does E-ext fail the way it did on heavy-hex, i.e. is $\kappa_u\ne\kappa_d$ for these observables?
- **Exit criterion:** a table like `../../detected_event_extrapolation/heavyhex/hh_tables.md` and a go/no-go on whether any encoded protocol beats U-ZNE at matched shots.

**Stage 2: two blocks, four logical qubits (open chain). Goal: inter-block gates, light cones, no-reset contamination.** Tool: Pauli fault-path Monte Carlo with a statevector. 8 data + 4 ancillas + 1 bridge = 13 qubits, $2^{13}$ amplitudes, batched over trajectories. This produces shot records $(v_i,Q_i)$ in the hardware format, so finite-shot analysis and ORP's held-out procedure run unchanged.
- **Additions over Stage 1:**
  - inter-block $R_{ZZ}$ via fan-in through the bridge qubit or ancilla;
  - the stale-ancilla mechanism: a fired, unreset ancilla flips the sign of later rotations routed through it ($J_{\rm eff}=J(1-2p_{\rm flip})$, [R1] App. C), as a non-Pauli, detector-correlated error;
  - reset vs no-reset;
  - detector ranking by $\Delta_j$ vs Pauli-propagation light cone.
- **Questions:** Q1–Q4 again, plus Q5: does stale-ancilla contamination break the product structure that E-soft and E-ext assume, and does E-PS+ZNE absorb it?

**Stage 3: many blocks, Clifford stand-in. Goal: shot-loss scaling at realistic acceptance.** Tool: stim, 4 to 21 blocks, rotation angles at multiples of π/2 (training-circuit style), [R1]'s detector layout and $R$.
- Measure $N_{\rm eff}$ of E-soft against ORP's $N\alpha_{k^\ast}$ and full post-selection at acceptance $\alpha_0$ down to $10^{-3}$.
- Measure the ORP plateau statistics and the cost of the second step (ZNE or PEC on the undetected sector) at scale.
- The noise per block matches Stage 1–2, so the non-Clifford bias measured there can be transferred.

**Stage 4: data reanalysis (parallel track).** Item 1 of the earlier plan: request the shot-level detector records of [R1] and apply E-soft and E-ORP+ZNE. This is not possible with E-PS+ZNE unless amplified runs exist.

### 4.4 Lessons carried over from the physical-qubit study

- **Count every gate.** Classify every single fault exhaustively before trusting any coin or visibility claim. The heavy-hex "compilation twirl" claim failed this test (`../../detected_event_extrapolation/heavyhex/HEAVYHEX.md` §3).
- **Treat structural-$r$ extrapolation as a hypothesis.** It worked where detected and undetected faults were equally harmful and failed by 0.03–0.05 elsewhere; report it per observable.
- **Give ZNE exact gain scaling in the baseline**, the best case for ZNE. Otherwise the comparison flatters the encoded route.
- **Report $N_{\rm req}$, not just bias.** Post-selection is cheap in bias and expensive in shots, and the break-even sits near a clean fraction of about 0.3.

### 4.5 Files (planned)

```
logical_qubit/error_detection/
  iceberg.py        # block layout, gadgets (GHZ prep, fan-in rotations, flagged extraction, switching), detectors
  sim_dm.py         # Stage 1 exact DM with detector branching
  sim_mc.py         # Stage 2 fault-path Monte Carlo, shot records
  estimators.py     # ORP (+ plateau finder), soft/geometric, two-point ext, ZNE fits, PEC reference
  stage0_memory.py, stage1_block.py, stage2_chain.py, stage3_scale.py
  test_iceberg.py   # 1-FT enumeration, ideal logical evolution vs 2-qubit exact, detector determinism
```

### 4.6 First concrete steps

1. `iceberg.py` with the single-block gadgets in stim, plus `test_iceberg.py`: the noiseless circuit gives deterministic detectors and the right logical evolution at Clifford angles, and single-fault enumeration shows which faults are undetected.
2. Stage 0 checks (i)–(iv).
3. `sim_dm.py` for one block at non-Clifford angles: verify the ideal logical evolution against the two-qubit Trotter circuit, then run the Stage 1 sweep at $p=0.3\%$ before the full noise scan.

## 5. Files

- Theory: `../../detected_event_extrapolation/theory/qed_extrapolation.md` (geometric law, estimator, variance, ORP relation, QEC extension) [N6]; `../../detected_event_extrapolation/theory/breakeven_theory.md` [N3]; `../../detected_event_extrapolation/theory/THEORY.md` [N1]; `../../detected_event_extrapolation/theory/hardware_native_checks.md`, `../../detected_event_extrapolation/theory/check_scheduling.md` [N4, N5].
- Paper: `../../detected_event_extrapolation/paper/main.pdf` (source `main.tex`), Sec. 10 "Detection-conditioned extrapolation in error-detecting codes" (`\label{sec:qed}`).
- Reference PDFs (not stored in this repository): [arXiv:2607.24947](https://arxiv.org/pdf/2607.24947), [arXiv:2609.13108](https://arxiv.org/pdf/2609.13108), [arXiv:2504.15725](https://arxiv.org/pdf/2504.15725), [arXiv:2212.03937](https://arxiv.org/pdf/2212.03937).

## References

**External**

- [R1] H. Froland, D. M. Grabowska, S. Grieninger, J. Hartse, A. L. Lashbrook, Z. Li, Z. Li, S. J. M. Powell, M. J. Savage, X. Yao, N. A. Zemlevskiy, *Realizing Error Suppression in Partially Fault-Tolerant Quantum Simulations with IBM Quantum Computers*, arXiv:2607.24947 (2026). — $[[4,2,2]]$ Iceberg blocks on heavy-hex, nFT logical rotations, FT syndrome extraction, Observable-Ranked Postselection.
- [R2] Fischer, Javadi-Abhari, Martiel, Seif, *Spacetime mitigation of logical errors*, arXiv:2609.13108 (2026). — Error detection with Pauli checks plus first-order PEC on the undetected generator; heavy-hex ancilla-mediated TFIM.
- [R3] S. Martiel, A. Javadi-Abhari, *Low-overhead error detection with spacetime codes*, arXiv:2504.15725 (2025). — Spacetime checks, periphery ancillas, controlled-$V$ repair for non-Clifford gates, greedy check selection.
- [R4] E. van den Berg, S. Bravyi, J. M. Gambetta, P. Jurcevic, D. Maslov, K. Temme, *Single-shot error mitigation by coherent Pauli checks*, arXiv:2212.03937, Phys. Rev. Research 5, 033193 (2023). — Gadget-noise model, critical payload error rate, one- vs two-sided checks.
- [R5] A. Gonzales, R. Shaydulin, Z. H. Saleem, M. Suchara, *Quantum error mitigation by Pauli check sandwiching*, arXiv:2206.00215, Sci. Rep. 13, 2122 (2023).
- [R6] D. M. Debroy, K. R. Brown, *Extended flag gadgets for low-overhead circuit verification*, arXiv:2009.07752, Phys. Rev. A 102, 052409 (2020).
- [R7] *Pauli Check Extrapolation*, arXiv:2406.14759 (2024).
- [R8] *Pauli-path reactivity functions*, arXiv:2609.31934 (2026). — Shared language with the physical-qubit study (post-selection filter, heralding as ZNE from below).

**Internal notes (this repository)**

- [N1] `detected_event_extrapolation/theory/THEORY.md` — detection twirl, mixture law, species and stratification, non-Clifford constraints.
- [N2] `detected_event_extrapolation/demo/THRESHOLD.md` — noise-level scan, break-even $q_0\approx0.3$.
- [N3] `detected_event_extrapolation/theory/breakeven_theory.md` — cost exponents, no threshold in $p$ alone, scaling condition $1+r<G_{\max}\kappa$.
- [N4] `detected_event_extrapolation/theory/hardware_native_checks.md` — heavy-hex fit (6-ring, mediators, periphery ancillas), why mediated bonds admit no compilation twirl.
- [N5] `detected_event_extrapolation/theory/check_scheduling.md` — check scheduling variants, enable-probability rule, terminal-vs-per-step measurement.
- [N6] `detected_event_extrapolation/theory/qed_extrapolation.md` — geometric law, soft ORP, second step, QEC extension.
