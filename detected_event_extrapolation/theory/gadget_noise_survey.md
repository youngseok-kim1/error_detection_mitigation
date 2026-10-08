# How the spacetime-code / Pauli-check literature deals with gadget noise and break-even

Survey compiled 7 October 2026 from the papers listed at the end (arXiv versions). Motivation: in our TFIM demo the check gadgets (controlled-Paulis, controlled-X repair gates, ancilla readout) add as much noise as the checks remove at 1% two-qubit error, and the repaired full twirl is the most accurate protocol but also the most expensive one at depth. The question is whether the literature has a threshold or break-even notion for this, and how it keeps the gadget overhead down. Items marked [unverified] could not be confirmed from the fetched text.

## 1. Is the gadget noise modelled?

| paper | gadget noise | what is included |
|---|---|---|
| van den Berg, Bravyi, Gambetta, Jurcevic, Maslov, Temme, *Single-shot error mitigation by coherent Pauli checks* (PRR 2023, 2212.03937) | yes, analytically | Two-qubit depolarizing ε on every CNOT of the check; per-check clean probability $t_{ok}=(1-\varepsilon)^k$ with $k$ the CNOT count of the check (≈ $3n/2$ two-sided all-to-all, $3n/4$ left-only, $9n/2$ / $9n/4$ on a line). Detectable check-gate fault probability $8\varepsilon/15$ per CNOT; $t_d=\tfrac12(1-(1-2p)^k)$, $t_u=1-t_d-t_{ok}$. Readout as symmetric flip $m$. Idle/T1/T2 on the check qubits in the "realistic" model. Flag qubits added because X/Y faults on the control of a controlled-P propagate to the data. |
| Gonzales, Shaydulin, Saleem, Suchara, *Pauli check sandwiching* (Sci. Rep. 2023, 2206.00215) | theorem: no; numerics: yes | Exact recovery theorem assumes noiseless checks. Simulations put 1q/2q depolarizing after every gate including check gates (2q = 10× 1q); measurements noiseless. Checks searched in increasing weight to keep check-induced noise low. |
| Debroy & Brown, *Extended flag gadgets* (PRA 2020, 2009.07752) | yes, as fault count | Each two-qubit gate adds 6 possible faults (3 control, 3 target). Score $q=N_{\rm detected}-6w(P)-6w(P')$. |
| Martiel & Javadi-Abhari, *Low-overhead error detection with spacetime codes* (2504.15725) | yes, in check scoring | Checks scored on the circuit including their own gates and wires, with ε = 8×10⁻⁴ after 2q gates plus idle; greedy rule "if a check's score decreases compared to the prior check, abort check picking for this path". Hardware (Heron r2, CZ 0.23–0.27%, readout 0.56–0.61%). |
| Fischer, Javadi-Abhari, Martiel, Seif, *Spacetime mitigation of logical errors* (2609.13108, anchor) | yes, fully | Learned sparse Pauli–Lindblad model covers data and check qubits ("captures both the additional check-circuit noise and its suppression in the post-selected results"); check SPAM learned separately; TREX; leakage filter on data and ancillas. ibm_aachen: CZ 1.6×10⁻³, readout 0.49%. |
| Martiel et al., *Sampling hard circuits with verifiably high fidelity* (2607.25941) | yes [partly inferred] | Local Pauli channels after every 2q gate; 0.1% CZ, 0.3% readout in simulation. 70 data + 27 ancillas, 454 of 2869 CZ are syndrome extraction (19% overhead); acceptance 5.9×10⁻⁴, 29× fidelity gain at 860× sampling cost. |
| *Pauli Check Extrapolation* (2406.14759) | depolarizing on all gates incl. checks | "Each additional layer of checks introduces more noise into the circuit"; the analytic Markov model nevertheless assumes near-perfect checks. |
| Liao et al. (PRX Quantum 2025, 2411.14638) | only via SWAP comparison | Flags needing one extra SWAP are net negative. |

## 2. Break-even / threshold statements

Nobody states a two-qubit-error threshold of the form "checks help iff $p_2 < x$" for a general payload. What exists:

- **van den Berg et al.** is the only explicit construction. They define a *critical payload error rate* $P_{\rm critical}$: checks improve the logical error only when the payload error $P_{\min}(k)$ exceeds it; otherwise applying checks increases the logical error (Sec. 4.4; the value is read off a figure, no closed form). In the many-check limit the logical error is $E=t_u/(\tfrac12-t_d)$ **only if** $t_{ok}>\tfrac12$; if a check's own gates fail with probability ≥ ½ the logical error goes to 1. For small ε this gives a floor $E\approx\tfrac{14}{15}k\varepsilon\approx\tfrac{7}{5}n\varepsilon$ for two-sided all-to-all checks: set entirely by gadget noise and linear in payload width. One-sided performance saturates at about ten checks; two-sided checks on hardware got *worse* with more checks (idle between L and R).
- **Gonzales et al.**: case statements only. For small circuits "the fidelity gain is negative because the checks introduce more errors than they eliminate"; gain is concave in the number of layers.
- **Debroy & Brown**: qualitative, "determined by the number of errors detected by the flag relative to the number of errors added by the flag"; many random single flags were worse than no flag, flag pairs almost always helped.
- **Martiel & Javadi-Abhari**: "the success of this approach hinges on our ability to design checks that detect more errors than they introduce"; enforced greedily per check, not as a closed form.
- **Fischer et al.**: "ED+PEC is advantageous when this reduction outweighs the cost of syndrome rejection and the additional noise introduced by the checks, hence the need for gate-efficient error detection." Leading order $\Gamma_{\rm ED+PEC}\approx\alpha^{-1}\exp[4\sum_{\rm undetected}\lambda]$ vs $\Gamma_{\rm PEC}=\exp[4\sum_{\rm all}\lambda]$; advantage increases with volume (3.7×, 15.9×, 63× at $d=2,4,6$). No numeric crossover.

**Relation to our observation.** Our break-even is the first-order bookkeeping "removed detectable fault weight > added gadget fault weight + the $\ln(1/\alpha)$ sampling cost", which is exactly van den Berg's $P_{\rm critical}$ logic. Their gadget floor $\tfrac75 n\varepsilon$ is, for $n=5$ and $\varepsilon=0.01$, 7% per check round, comparable to the 10% per-step payload error, which is the regime in which our fixed checks are useless and the twirled protocols are cost-limited. The threshold scan in `demo/threshold_scan.py` makes this quantitative for our estimator.

## 3. Overhead-reduction techniques in the literature

- **Single-sided checks** (van den Berg; Martiel & Javadi-Abhari): the right check is evaluated classically on the measured bits; halves the gadget gate count and removes the L–R idle time. Two-sided checks "effectively double the gate and depth overhead".
- **Local low-weight checks on a dedicated ancilla** (Martiel & Javadi-Abhari; Martiel et al.; Fischer): each ancilla supports Paulis only on its nearest data neighbours; greedy search over syndrome-decoding solutions with the "stop when the score worsens" rule.
- **Dual-purpose ancillas** (Fischer): the ancillas that mediate ZZ on heavy-hex *are* the check qubits (27 checks for 22 data qubits), measured once at the end, no mid-circuit measurement or reset.
- **Half-SWAP routing** (Martiel SI, Figs S11–S12): CZs between two ancillas commute with controlled-Paulis in pairs; SWAP+CZ reduces to a cost-2 circuit; ancillas detect each other's errors.
- **Ancilla pairs per data qubit** on heavy-hex so the pair cross-checks itself (Martiel).
- **Flag qubits on the check qubits** (van den Berg; Debroy & Brown): lower logical error at small check number but lower acceptance; harmful under high idle noise.
- **LNN gate merging** (van den Berg Fig. 4): controlled-X/Y/Z followed by SWAP in 2 CNOTs instead of 4.
- **Check-count cap or extrapolation**: ~10 one-sided checks (van den Berg); the γ-curve stopping rule (qiskit-addon-paulice); Pauli Check Extrapolation fits 2–4 check layers and extrapolates to the all-checks limit.
- **Weighted ensembles instead of hard post-selection** (Langfitt et al., 2408.05565).
- **Readout-aware qubit choice**: fold parities onto the best-readout qubits (Martiel); layout chosen to maximize the product of measurement fidelities (Fischer).

## 4. Non-Clifford gates

- **Commutation constraint** (Debroy & Brown; Gonzales; Martiel; Fischer App. B.3): the back-propagated check must commute with every rotation it meets. Gonzales: with 15 $R_z$ gates no six-layer check set exists for ≥ 20-CNOT circuits. Martiel: "the number of valid checks decreases exponentially with the number of rotations" (Fig. S8); each rotation adds one linear mod-2 constraint.
- **Controlled-V repair** (Martiel SI §V.B): for a check anticommuting with $R_P$ on wire $w$, insert controlled-$V$ ($V$ anticommuting with $P$) just before and after; the added gates are Clifford, two per repaired rotation. Caveat stated there: a check compatible with a rotation of axis $P$ on $w$ cannot detect a $P$ error on $w$, and after twirling "any error due to the T gates remains undetected by valid checks". This is the repair we test as Stage 3c.
- **Doped placement** (Martiel et al. 2607.25941): T gates placed only where they commute with the check back-cumulants.
- **Trotter-native compatibility** (Fischer): ancilla-assisted ZZ with $R_z(\phi)$ on the ancilla keeps the propagated checks Z-type and the data $R_x(\theta)$ commute with checks on the check register. No repair gates needed — the compilation is chosen so that the rotation axis is never on a data wire the check must see through.
- **Flagging rotation products** (Debroy & Brown): flag the product of consecutive multi-qubit Pauli rotations to improve the detected/added ratio.

## 5. What we can borrow for the non-Clifford Trotter circuit

1. Compile so that checks commute with the rotations by construction (Fischer: ZZ via an ancilla carrying the $R_z$), which removes the need for the controlled-X repair and the invisible sector moves onto the ancilla's own rotation.
2. If repair is unavoidable, keep the number of anticommuting bonds per step small by biasing the flavor distribution; the repair cost is two gates per anticommuting bond.
3. Single-sided (terminal-measurement) checks halve the gadget count; our two-sided coherent checks pay for the right-hand half.
4. Dual-purpose ancillas measured once at the end avoid mid-circuit reset and readout on every step.
5. The van den Berg $P_{\min}/P_{\rm critical}$ bookkeeping with the measured ε decides a priori whether a window is worth checking; the greedy "stop when the score worsens" rule of Martiel does the same per check.
6. Cap the number of checks (saturation near ten one-sided checks) or extrapolate in the number of check layers rather than stacking them.

## Sources

- Fischer, Javadi-Abhari, Martiel, Seif, *Spacetime mitigation of logical errors*, arXiv:2609.13108
- van den Berg et al., *Single-shot error mitigation by coherent Pauli checks*, arXiv:2212.03937 (PRR 2023)
- Gonzales, Shaydulin, Saleem, Suchara, *Quantum error mitigation by Pauli check sandwiching*, arXiv:2206.00215 (Sci. Rep. 2023)
- Debroy & Brown, *Extended flag gadgets for low-overhead circuit verification*, arXiv:2009.07752 (PRA 2020)
- Martiel & Javadi-Abhari, *Low-overhead error detection with spacetime codes*, arXiv:2504.15725
- Martiel et al., *Sampling hard circuits with verifiably high fidelity*, arXiv:2607.25941
- Delfosse & Paetznick, *Spacetime codes of Clifford circuits*, arXiv:2304.05943
- Bacon, Flammia, Harrow, Shi, *Sparse quantum codes from quantum circuits*, arXiv:1411.3334
- *Pauli Check Extrapolation*, arXiv:2406.14759
- Langfitt et al., *Pauli Check Sandwiching for quantum characterization*, arXiv:2408.05565
- Liao et al., PRX Quantum 6, 020331 (2025), arXiv:2411.14638
- Aharonov et al., *Syndrome aware mitigation of logical errors*, arXiv:2512.23810
- Zhong et al., *Combining error detection and mitigation*, arXiv:2510.01181
- Froland et al., *The utility of sparse error detection*, arXiv:2608.02944
- IBM qiskit-addon-paulice guide, *Low-overhead error detection using spacetime codes*
