# How to schedule the checks in time on the heavy-hex 6-ring

7 October 2026. Design memo before building the simulator. Geometry fixed by `hardware_native_checks.md` §6: data $q_0..q_5$ on the plaquette vertices, mediators $m_i$ between $q_i$ and $q_{i+1}$, peripheral ancillas $p_i$ on $q_i$.

## 1. The time structure we have to respect

One Trotter step, transpiled:

| layer | what runs | CZ-depth | who is busy |
|---|---|---|---|
| R | $R_x(\theta)$ on all data; frame changes merged | 0 (single-qubit) | data only |
| E (even bonds $m_0,m_2,m_4$) | 3 gadgets in parallel, 4 CZ slots each: slot 1 CX($q_{\rm first}\!\to\! m$), slot 2 CX($q_{\rm second}\!\to\! m$), $R_z(m)$, slot 3 CX($q_{\rm second}\!\to\! m$), slot 4 CX($q_{\rm first}\!\to\! m$) | 4 | each data qubit busy in 2 of the 4 slots, idle in the other 2 |
| O (odd bonds $m_1,m_3,m_5$) | same | 4 | same, roles shifted by one |

Facts that drive the design:

- **Every data qubit idles for half of every sub-layer.** The "first" qubit of a gadget is busy in slots 1 and 4, the "second" in slots 2 and 3. Which qubit is first is our choice, per gadget and per step.
- **Mediators idle for the whole opposite sub-layer** (and for layer R); peripherals idle unless we use them.
- **A mid-circuit measurement plus reset costs roughly 10 CZ durations** on current IBM hardware (readout $\sim1\,\mu$s against CZ $\sim0.1\,\mu$s). If an ancilla is measured every step and the data wait for it, the idling alone adds the equivalent of a whole sub-layer of decoherence per step; this, more than gate count, is why the anchor paper measures only at the end. An ancilla can be measured *while* the data continue only if nobody needs it for about one step.
- **Check windows cannot cross layer R.** A single-wire image that is valid through a ZZ layer is the data-frame type ($Z$-type logically), which anticommutes with $R_x$; the $X$-type image that survives $R_x$ dies in the ZZ layers. So every check window lives inside one ZZ layer (one or both sub-layers), and its opening and closing controlled-Paulis sit at that window's boundaries.

## 2. Variants

**V0 — anchor.** Mediators only, fixed $Z$-type flavor, measured once at the end. Zero overhead. Undetected: $Z$-type data faults (7/15 per CX), data faults in layer R, and the data-only faults after a mediator's uncompute CXs. Cross-step parity cancellation on each mediator.

**V1 — V0 + compilation twirl.** Per data qubit and step draw the frame type $t_i\in\{X,Y,Z\}$ (and the mediator axis); compile each gadget to collect the parity of the two current-frame operators. Zero gate overhead (the conjugations merge into layer R). Every mediator-visible fault now has coin $\tfrac23$ (single-site, mediator-involving) or $\tfrac49$ (data$\otimes$mediator); data-only faults after the uncompute CXs (slots 3 and 4 of the *other* qubit, 3/15 at two of four locations, 10% of CZ faults) stay invisible. Still terminal measurement.

**V2 — idle-slot peripheral checks (zero depth overhead).** Let $p_i$ check $q_i$ during the gadget in which $q_i$ is the *second* qubit: open with a controlled-$t_i$ in slot 1 (when $q_i$ is idle), close with controlled-$t_i$ in slot 4 (idle again). The window covers slots 2 and 3, i.e. $q_i$'s two CX locations in that gadget, plus the mediator rotation time. No extra CZ depth: the two controlled-Paulis sit in slots the data qubit would otherwise spend waiting. Choose roles so that every qubit is second in exactly one of its two gadgets per step (even qubits second in E, odd in O, or the reverse) — then **every data qubit is checked once per step by its peripheral at zero depth cost**, 2 CZ per qubit per step (12 per step for the ring). Coverage gap: the locations where $q_i$ is the first qubit (slots 1 and 4 of its other gadget) are seen only by that mediator, and its slot-4 data-only faults are the invisible 10% of V1. Fix by **alternating roles across steps** (even qubits second in E on even steps, second in O on odd steps): every location is then peripheral-covered in half of the steps, coin $\tfrac12\cdot\tfrac23=\tfrac13$ at those locations, nothing invisible.

**V3 — idle mediators as coherent checks ("use every ancilla every layer").** $m_i$ is idle during the sub-layer in which its two neighbours sit in gadgets with $m_{i-1}$ and $m_{i+1}$. After its own gadget it is ideally in $|0\rangle$; a Hadamard turns it into a coherent-check ancilla, and controlled-$t_i$ on $q_i$ and controlled-$t_{i+1}$ on $q_{i+1}$ at the start and end of the opposite sub-layer (4 CZ, placed in those qubits' idle slots) make it a weight-2 check on both neighbours; a second Hadamard returns it to $|0\rangle$ for its next gadget, and the terminal $Z$ measurement reads the parity of its mediator flags and its coherent-check flags over the whole run. This is the "check further-away data synchronously" idea: the mediator reaches two data qubits the peripheral cannot, at the cost of merging all its detections into one parity bit. Zero depth overhead if the slots fit (they do: a second qubit's idle slots 1 and 4 and a first qubit's idle slots 2 and 3; a weight-2 check needs one of each at the window boundary, which the two-sub-layer structure provides only for the pair (first, second) — the scheduling constraint is real and is the first thing the simulator should check).

**V4 — time-sparse checks (enable twirl).** Any of V2/V3 with each check enabled with probability $b<1$ per step, independently; shots carry the schedule. This is the knob that trades gadget noise against detector strength (§3). Deterministic sparsity (check every $k$-th step) is *not* equivalent: faults in unchecked steps would be invisible and bias the estimate; random sparsity keeps every fault's coin nonzero.

**V5 — measurement pipelining.** With MCM, measure each ancilla every other step while the data continue (an ancilla busy with readout for $\sim$1 step), using the two halves of the ancillas alternately; or measure every $k$ steps with parity accumulation in between. Interpolates between per-step measurement (no cross-step cancellation, heavy idling) and terminal measurement (no idling, cancellation).

## 3. Two quantitative design rules

**Enable probability.** With a visible-fault weight $\Lambda_{\rm pay}$ in the payload, a gadget weight $\Lambda_g$ at $b=1$, and a per-fault detection probability $c$ when the covering check is on, the leak-subtraction variance is
$$N\,\mathrm{Var}\;\propto\;\frac{e^{2\Lambda_{\rm pay}}\,e^{2b\Lambda_g}}{(b\,c)^2},\qquad\text{minimized at}\quad b^\ast=\min\!\Big(1,\ \frac{1}{\Lambda_g}\Big).$$
Check every step only while the whole gadget budget is below one expected fault; beyond that, enable each check with probability $1/\Lambda_g$. For the 6-ring with peripherals, $\Lambda_g=2\cdot6\cdot\tfrac{14}{15}pd\approx11\,pd$: $b^\ast=1$ up to $pd\approx0.09$, $b^\ast\approx0.45$ at $pd=0.2$. The mediator checks are free and always on.

**Terminal versus per-step measurement.** A mediator accumulating over $d$ steps has visible weight $\lambda_m d$ with $\lambda_m\approx4\cdot\tfrac{14}{15}p\cdot\tfrac23\approx2.5p$ per step; two visible faults cancel with probability $\approx(\lambda_m d)^2/2$ per mediator, so the fraction of faulty shots that pass all mediator checks by cancellation is $\approx6\cdot(\lambda_md)^2/2\,/\,\Lambda_{\rm tot}$, i.e. about $8\%$ of faulty shots at $pd=0.1$ and $17\%$ at $pd=0.2$ (with $\Lambda_{\rm tot}\approx22\,pd$ for the mediated payload). That is a second-order invisible species comparable to the first-order residuals we have been removing, so at depth either the peripherals must be measured every $k\lesssim4$ steps (pipelined, V5) or the cancellation species must be included in the inversion (its survival factor is known). Per-step MCM for all twelve ancillas would cost $\sim10$ CZ durations of data idling per step and is not acceptable.

## 4. Recommendation and what to simulate

Baseline: **V1** (mediators, compilation twirl, terminal measurement) — zero overhead, already removes the model dependence of the anchor protocol except for the 10% invisible data-only locations.

Protocol: **V1 + V2 with role alternation + V4** — peripheral checks in the idle slots (12 CZ per step, zero depth), roles alternating so nothing is invisible, enable probability $b^\ast=\min(1,1/\Lambda_g)$; peripherals measured every other step in two pipelined halves (V5), mediators at the end. Gate budget per step: 24 CZ payload (connectivity-forced) + $12b$ CZ checks, no added depth, 3 ancilla readouts per step in flight.

Optional: **V3** on top, if the slot analysis confirms the weight-2 windows fit; it adds coverage of the first-qubit locations without peripherals but merges detections into one parity bit per mediator.

The simulator (fault-path Monte Carlo, `hardware_native_checks.md` §6) should take the schedule as data — for every check: ancilla, (qubit, slot) of each controlled-Pauli, enable probability, measurement step — so that V1–V5 are configurations, not code paths; the flag of every ancilla is the parity, over its windows, of the anticommutation of the sampled faults with the images; the two-dimensional threshold rules (mediator flags, peripheral flags) feed the species inversion. Compare against ZNE on the same mediated circuit at $p=q\in\{1\%,0.5\%,0.2\%\}$, $d\le12$.
