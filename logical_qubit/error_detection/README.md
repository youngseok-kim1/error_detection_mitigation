# error_detection — detection-conditioned extrapolation for encoded (logical-qubit) simulations

Project folder: `logical_qubit/error_detection/`. Started 8 October 2026. Companion to `../../detected_event_extrapolation/` (physical-qubit Pauli-check study), whose theory notes this document builds on.

## 1. Where this comes from

The physical-qubit study (`detected_event_extrapolation/`) asked whether error-detection events can drive a ZNE-type extrapolation. Its results, in short:

- Twirling the flavor of a coherent Pauli check turns detection into a fair coin and makes partial post-selection an exact mixture law; the fault-free value follows by linear extrapolation in the accepted fraction, model-free, to all orders for arbitrary stochastic Pauli noise [N1].
- On a deep non-Clifford TFIM circuit this beats exponential ZNE in accuracy at every depth, but only beats it in cost while the checked circuit carries about one expected fault, $(g_p+g)\,p\,d\lesssim1.3$; the cost exponents are all proportional to $pd$, so there is no threshold in the gate error alone, and the depth-independent comparison favours checks only if the gadget overhead ratio $r$ and the observable sensitivity $\kappa$ satisfy $1+r<G_{\max}\kappa$ [N2, N3].
- On heavy-hex hardware with ancilla-mediated bonds the flavor twirl is not available at all (every valid check image at a CZ noise location is $Z$-type), $Z$-type data faults are intrinsically invisible, and the check gadgets cost as much as they remove [N4, N5]. **Conclusion: for ZNE-style mitigation of physical circuits the approach brings no immediate advantage.**

The natural home of the same logic is a quantum error-detecting (QED) code: the detectors exist anyway, the circuit is already paying for them, and the problem is the opposite one — full post-selection on trivial syndromes leaves no shots once the circuit is large, and the detector record is thrown away together with the shots. The reference experiment is Froland et al. [R1]: 21 blocks of the $[[4,2,2]]$ Iceberg code on ibm_boston (42 logical qubits, up to 136 physical), Trotterized mixed-field Ising dynamics in 1+1D and 2+1D, non-fault-tolerant logical rotations and fault-tolerant syndrome extraction, and a selective-filtering method, Observable-Ranked Postselection (ORP), to avoid the exponential shot loss.

## 2. How the reference experiment evolves the Ising model with logical qubits

- Code: $S_X=X_0X_1X_2X_3$, $S_Z=Z_0Z_1Z_2Z_3$, logical $|00\rangle=$ GHZ(4), logical operators $\bar X_i=X_{i+1}X_3$, $\bar Z_i=Z_0Z_{i+1}$ ($i=0,1$). Every logical Pauli is a weight-2 physical Pauli; $\bar Z_0\bar Z_1=Z_1Z_2$; inter-block $\bar Z\bar Z$ is weight 4.
- Hamiltonian $H=-J\sum\bar Z_i\bar Z_j-\sum(g_x\bar X_i+g_z\bar Z_i)$, first-order Trotter; $J=1$, $g_x=0.75$, $g_z=0.5$ (bubble melting) or $1.5$ (localized).
- Non-Clifford logical gates are physical multi-qubit Pauli rotations done non-fault-tolerantly: a CNOT fan-in collects the parity of the support onto one central qubit (a data or bridge qubit, routed through the idle ancilla when needed), one physical $R_z$ is applied there, and the fan-in is uncomputed (depth 4–12 CNOTs per inter-block $R_{ZZ}$). The single-logical-qubit rotations and the intra-block $R_{ZZ}$ are folded into the same fan-in chains at the cost of single-qubit gates only. A per-block frame bit tracks the SWAP$(d_1,d_2)$ needed for $\bar Z_1,\bar X_0$. One fault on the parity qubit between fan-in and uncompute becomes a weight-2, undetectable logical error: the $Ap$ term of $p_L=Ap+Bp^2$. In our language, a fault along the logical rotation axis is a logical operator and hence invisible to the stabilizers — the rotation-axis sector of [N1] one level up.
- Only syndrome extraction is fault-tolerant (1-FT circuit, the two ancillas $a_x,a_z$ as mutual flags), run in a separate "syndrome" qubit configuration; switching costs depth $\approx15$ each way, so syndromes are extracted every $R$ Trotter steps. Ancillas are not reset; detectors are parities of syndrome bits two rounds apart, and the terminal data readout closes the last detector.
- Encoded circuits are 2–6× deeper than unencoded; local observables improve by 2–6% (1+1D) and over 200% at late times (2+1D), with no noise-learning mitigation.

## 3. ORP and what we add

ORP [R1, Methods A, App. A]: for each detector $v_j$ compute from data $\Delta_j=E[Q\mid v_j{=}0]-E[Q\mid v_j{=}1]$, rank detectors by $\Delta_j$, keep the shots in which none of the top-$k$ detectors fired, choose $k^\ast$ by a plateau-finding heuristic on held-out shots; blockwise terminal parity is always enforced. At weak noise $\Delta_j\approx2\langle O\rangle\sum_{m\ni j}p_mq_m/\sum_{m\ni j}p_m$, the flip probability of the processes firing $v_j$; the ranked detectors coincide with the backward light cone of $O$. An additive regression metric $\beta_j$ is defined and used only as an alternative ranking.

**Geometric law** (theory in [N6]; paper Sec. 9). In a Clifford code circuit with independent Pauli faults, let $k$ be the number of detected faults in a shot and $o=\pm1$ the logical outcome. Then, exactly,

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

1. **Reanalysis of existing data (cheapest, highest value).** The estimator needs only the shot-level $(v_i,Q_i)$ records ORP uses. Ask the authors of [R1] for the detector-level records of the ibm_boston runs; compare $\hat O_0\pm$ bootstrap with $\langle O\rangle_{k^\ast}\pm\sigma_{k^\ast}$ for the observables of their Figs. 2–4 and run the linearity diagnostic. Expected: same central values, error bars reduced by $\sqrt{\exp[\sum_{j\le k^\ast}p_j\rho_j^2]}$ — a factor of several at their late-time acceptance fractions.
2. **Controlled simulation of the assumptions.** One and two $[[4,2,2]]$ blocks (4 data + $a_x,a_z$ each) modelled as in [R1]: GHZ preparation, fan-in / $R_z$ / uncompute logical rotations, flagged extraction every $R$ steps, no-reset detectors, two-qubit depolarizing noise at $0.1$–$1\%$, readout flips. Exact density matrix branched by detector pattern (6–12 qubits), or stim for Clifford variants at scale. Measure: linearity of $\ln E[Q\mid k]$ and $\rho_j$ against the fault enumeration; bias and variance of soft ORP vs full post-selection vs ORP cuts; the Clifford vs non-Clifford (rotation angle) breakdown; the $G^2$ extrapolation with detector-measured gain.
3. **Scale test with stim.** 21-block Clifford stand-in (rotation angles $\pi/2$) with the real detector layout and the acceptance fractions of [R1], to confirm the $N_{\rm eff}$ formula at realistic $\alpha_0$.
4. **Write-up.** Fold the results into Sec. 9 of `detected_event_extrapolation/paper/main.tex` or split into a standalone note, depending on what the reanalysis shows.

## 5. Files

- Theory: `../../detected_event_extrapolation/theory/qed_extrapolation.md` (geometric law, estimator, variance, ORP relation, QEC extension) [N6]; `breakeven_theory.md` [N3]; `THEORY.md` [N1]; `hardware_native_checks.md`, `check_scheduling.md` [N4, N5].
- Paper: `../../detected_event_extrapolation/paper/main.pdf`, Sec. 9 "Detection-conditioned extrapolation in error-detecting codes".
- Reference PDFs: `.../space_time_check/references/2607.24947v1.pdf`, `2609.13108v1.pdf`, `2504.15725v1.pdf`, `2212.03937v1.pdf`.

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
