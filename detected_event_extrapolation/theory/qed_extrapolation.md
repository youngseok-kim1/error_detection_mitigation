# Detection-conditioned extrapolation in error-detecting and error-correcting codes

8 October 2026. Theory note. Motivation: in a distance-2 code (Iceberg $[[k{+}2,k,2]]$; the reference experiment is Froland et al., arXiv:2607.24947, 21 blocks of $[[4,2,2]]$ on ibm_boston, 42 logical qubits, Ising dynamics with non-fault-tolerant logical rotations and periodic fault-tolerant syndrome extraction) full post-selection on trivial syndromes leaves no shots once the circuit is large; the detector record is then thrown away together with the shots. The question is whether the detector events can instead be used the way we used Pauli-check flags: not to say *which* error happened, but *how much* noise a shot carried and, under stated assumptions, what that does to the observable, so that every shot contributes and the fault-free value is reached by extrapolation.

The answer has three layers: (1) an exact **geometric law** for the observable conditioned on the number of detected faults, which turns full post-selection into a regression over all shots and multiplies the effective shot count by $e^{\Lambda_d\bar\rho^2}$; (2) a **second step** for the undetected sector, which the detection count cannot see but can calibrate; (3) the **QEC version**, where the decoder's soft output replaces the detection count as the per-shot noise meter.

## 1. Setting

A Clifford code circuit (encoding, logical Clifford layers, periodic syndrome extraction, terminal measurement) with $K$ detectors $D_1,\dots,D_K$ (parities of measurement outcomes that are deterministic in the absence of faults) and a logical observable $O$ with outcome $o\in\{\pm1\}$. A shot returns $(s,o)$ with $s\in\{0,1\}^K$ the detector vector. Noise: independent stochastic Pauli faults $\nu$ at circuit locations with rates $\lambda_\nu$ (Pauli–Lindblad), as in the anchor paper. For each fault $\nu$:

- its **syndrome** $s_\nu\in\{0,1\}^K$ (which detectors it fires); $\nu$ is *detected* if $s_\nu\neq0$, *undetected* otherwise;
- its **logical effect** $\epsilon_\nu\in\{\pm1\}$: by the Pauli-frame push-through, a Pauli fault in a Clifford circuit either flips the outcome of $O$ or not.

Write $\Lambda_d=\sum_{\nu\ {\rm det}}\lambda_\nu$ and $\Lambda_u=\sum_{\nu\ {\rm undet}}\lambda_\nu$. For a distance-2 code every single-location fault of weight one is detected; the undetected set is weight-two combinations with trivial syndrome and the logical-operator-like faults, so $\Lambda_u=O(\Lambda^2)$ at the level of *fault configurations* (two independent detected faults whose syndromes cancel). To keep the bookkeeping first-order we treat syndrome cancellation as a second-order correction in Sec.~4.

Full post-selection keeps $s=0$: acceptance $\alpha_0=e^{-\Lambda_d}$ to first order, and the post-selected observable is
$$O_0\equiv E[o\mid s=0]=O_{\rm ideal}\prod_{\nu\ {\rm undet}}\bigl(1-2p_\nu\mathbf 1[\epsilon_\nu=-1]\bigr)\approx O_{\rm ideal}\,e^{-2\Lambda_u^{(O)}},$$
with $\Lambda_u^{(O)}$ the undetected weight that flips $O$. The "no shots" problem is $N\alpha_0\to0$.

## 2. The geometric law

Let $k(s)$ be the **number of detected faults** in the shot (Sec.~4 says how to get it from $s$). Poisson faults have the property that, conditioned on $k$ detected faults having occurred, the faults are $k$ independent draws from the detected-fault distribution $\pi_d(\nu)=\lambda_\nu/\Lambda_d$, independent of the undetected faults.

**Theorem 1 (geometric law).** Under independent Pauli faults in a Clifford circuit,
$$E[o\mid k]\;=\;O_0\,\bar\rho^{\,k},\qquad \bar\rho\;\equiv\;E_{\pi_d}[\epsilon_\nu]\;=\;1-2\,\Pr[\text{a detected fault flips }O].$$

*Proof.* $o=o_{\rm undet}\prod_{i=1}^k\epsilon_{\nu_i}$ with $o_{\rm undet}$ the outcome in the absence of the detected faults. Conditioning on $k$, the $\nu_i$ are iid $\pi_d$ and independent of the undetected configuration, so $E[o\mid k]=E[o_{\rm undet}]\prod_iE[\epsilon_{\nu_i}]=O_0\bar\rho^k$. $\square$

Remarks. (i) The law is exact to all orders in the detected faults; it does not expand in $\Lambda$. (ii) $\bar\rho$ is a property of the observable and the fault distribution, not of the circuit size: for a local logical observable in a large code, most detected faults are outside its light cone and have $\epsilon=+1$, so $\bar\rho$ is close to one. Writing $\kappa_d=\Pr[\text{detected fault flips }O]$, $\bar\rho=1-2\kappa_d$. (iii) $\ln E[o\mid k]$ is linear in $k$: the detection count is a *noise meter* and the observable follows a single exponential in it. This is the exact statement of "detection tells how much noise and what it does on average".

**Estimator.** Fit $\ln\bar o_k=\ln O_0+k\ln\bar\rho$ over the observed $k$ (weighted least squares with the binomial variances), or equivalently the two-parameter maximum-likelihood fit of $o_i\sim\pm1$ with mean $O_0\bar\rho^{k_i}$. The intercept is the fully post-selected value, now estimated from every shot. If $\bar\rho$ is known from a model it can be fixed; it need not be.

**Variance.** With $\bar\rho$ known, the Fisher information about $O_0$ is $\sum_iN_k\bar\rho^{2k}/\sigma^2$, i.e. the effective number of shots is
$$N_{\rm eff}=N\sum_kP(k)\,\bar\rho^{2k}=N\,e^{-\Lambda_d(1-\bar\rho^2)}=N\,\alpha_0^{\,1-\bar\rho^2},$$
against $N\alpha_0$ for full post-selection: a gain of $\alpha_0^{-\bar\rho^2}=e^{\Lambda_d\bar\rho^2}$. A shot with $k$ detections is worth $\bar\rho^{2k}$ of a clean shot. Example: $N=10^6$, $\alpha_0=10^{-3}$ (full post-selection leaves $10^3$ shots), $\kappa_d=0.05$ ($\bar\rho=0.9$): $N_{\rm eff}=10^6\cdot10^{-0.57}=2.7\times10^5$, 270 times more. With $\bar\rho$ fitted rather than known the gain is smaller by the usual factor for a two-parameter exponential fit (about 3–5 in variance when the $k$ range is a few units), still two orders of magnitude in this example.

**Relation to the mixture law.** In the language of THEORY.md, a code detector is a fixed check: coin 1 for detected faults, 0 for undetected. Post-selection is therefore a projection and cannot be rescaled to the ideal value; what the geometric law adds is that, with species labelled by the detected-fault count and iid damage per detected fault, the species contributions are a geometric sequence, so the $k=0$ species is recoverable from all the others. Theorem 1 is the stratified inversion of Sec.~4 of THEORY.md with a one-parameter species model.

## 3. What the geometric law does not do, and the second step

$O_0$ still carries the undetected damage $e^{-2\Lambda_u^{(O)}}$. The detection count cannot see it: conditioning on $k$ leaves the undetected configuration untouched (that is the content of the proof). Three ways to get from $O_0$ to $O_{\rm ideal}$:

**(a) Zero-noise extrapolation from below, with the detector as the noise meter.** Scale the physical noise by $G$ (gate folding, probabilistic error amplification). Then $\Lambda_d\to G\Lambda_d$, which is *measured* directly from the mean detection count $\bar k(G)=G\Lambda_d$, and $\Lambda_u^{(O)}\to G^2\Lambda_u^{(O)}$ for a distance-2 code (pair events) so that $\ln O_0(G)=\ln O_{\rm ideal}-2G^2\Lambda_u^{(O)}$. Fit at two or three gains and extrapolate to $G=0$. The detector removes the usual weakness of ZNE, the uncertainty in the achieved amplification: $G$ is read off the detection rate, not assumed. And because the extrapolation is in $G^2$ with the $G^1$ term already removed by the geometric law, it is a one-parameter extrapolation of a quantity that is quadratic in the noise, far better conditioned than extrapolating the raw observable.

**(b) Model-assisted correction.** $\Lambda_u^{(O)}=c\,\Lambda_d^2$ with $c$ a property of the code circuit (the fraction of detected-fault pairs with trivial syndrome and logical effect on $O$), computable from the circuit structure under an assumed noise *shape* (e.g. uniform depolarizing per gate), while $\Lambda_d$ is measured. One number from the model, the magnitude from the data.

**(c) Learning the rates from the detector record.** The detector statistics determine the per-location detected rates (the syndrome-based noise learning of Chen et al.); under the Pauli–Lindblad model the undetected rates at the same locations follow, and $O_0\to O_{\rm ideal}$ by first-order PEC on the undetected generator, as in the anchor paper but with the model learned from the same shots.

In all three, the *first* step has already done the expensive part: it replaced $N\alpha_0$ usable shots by $N\alpha_0^{1-\bar\rho^2}$.

## 4. From detector bits to the fault count, and the second-order corrections

The law is in the number of detected faults $k$, not in the number of fired detectors $|s|$. Three cases:

- *One syndrome round at the end (Iceberg-style, two stabilizers).* Then $|s|\le2$ whatever $k$ is: the detectors report only the parity of the number of detected faults in each stabilizer, and $k$ is hidden. The geometric law still holds for the hidden $k$, but what the experiment sees is its mixture, $E[o\mid s=0]=\sum_kP(k\mid s=0)\,O_0\bar\rho^{\,k}$, in which the even-$k$ cancellations contaminate the "clean" bin with probability $\approx\Lambda_d^2/2$ per stabilizer. Here syndrome cancellation is not a correction but the whole story, and the clean extrapolation needs more detectors. This is why periodic syndrome checks matter beyond their error-detection role: they make $k$ observable.
- *Periodic checks, one detector pair per check.* A fault fires the detectors of the window it falls in. If the windows are short enough that two faults rarely share one ($\Lambda_d/\#{\rm windows}\ll1$), $k\approx$ number of fired windows and the law holds up to the cancellation probability $\approx\Lambda_d^2/(2\,\#{\rm windows})$, which is a known second-order species (zero flags, two faults) with survival factor one; include it as in Sec.~4 of THEORY.md or keep $\Lambda_d/\#{\rm windows}$ small by design.
- *Surface-code-like detector lattices.* A data fault fires two adjacent detectors, a measurement fault two in time; $k$ is the number of *clusters* (the matching), not of detectors; a minimum-weight matching gives $k$ directly and the same law applies to its event count, with $\bar\rho$ now the mean logical effect per matched event.

A fourth subtlety: faults of different types have different $\epsilon$; the law needs only their *mixture* to be iid per detected fault, which is true for Poisson faults. If the noise is strongly location-dependent (one hot gate), the law still holds; $\bar\rho$ is then dominated by that gate.

## 5. Error-correcting codes: the decoder's soft output as the noise meter

For $d\ge3$ with a decoder, the detector record is used twice: to choose the correction and, through the decoder's **soft output** (the log-likelihood ratio $g$ between the chosen logical class and its complement, the "complementary gap"), to estimate the probability $P_L(g)=1/(1+e^{g})$ that the corrected shot carries a logical error. Then
$$E[o\mid g]=O_{\rm ideal}\,\bigl(1-2P_L(g)\bigr),$$
and a regression of $o$ on the known function $1-2P_L(g)$ over all shots gives $O_{\rm ideal}$ with effective shot count $N\,E[(1-2P_L)^2]$ — no post-selection, every shot weighted by the square of its decoder confidence. Post-selecting on $g>g_0$ is the hard-threshold version; the regression is the extrapolation version, and it is unbiased to the extent that the decoder's likelihoods are calibrated (which the detection statistics themselves can check: the fraction of shots with a given $g$ whose logical outcome disagrees with the majority must equal $P_L(g)$). This is the exact analogue of Theorem 1 with $\bar\rho^k$ replaced by the per-shot factor $1-2P_L(g)$.

## 6. Assumptions, stated

1. Independent Pauli faults (Pauli–Lindblad); achieved by Pauli twirling the gates.
2. Clifford circuit, so that each fault's logical effect is a sign. For non-Clifford circuits the effect of a fault on $\langle O\rangle$ is a factor in $[-1,1]$ and the product structure of the proof holds only when the light cones of distinct detected faults do not overlap before the measurement; then $E[o\mid k]\approx O_0\bar\rho^k$ with $\bar\rho$ the mean damage factor, and the deviation is second order in the overlap probability. The Iceberg demonstrations with arbitrary-angle logical rotations are in this regime only if faults are sparse.
3. The detector count resolves the fault count (Sec.~4), or the cancellation species is included.
4. The undetected sector is treated separately (Sec.~3); the geometric law alone recovers the fully post-selected value, not the ideal one.

## 7. What to test first

A stim simulation (exact detector statistics, fast) of a $[[4,2,2]]$ or $[[6,4,2]]$ Iceberg circuit with $m$ periodic syndrome checks, a logical Clifford payload of 50–200 two-qubit gates, depolarizing noise at $0.1$–$1\%$:
1. verify $\ln E[o\mid k]$ linear in $k$ and read off $\bar\rho$; compare with $1-2\kappa_d$ computed from the fault enumeration;
2. compare the variance of $\hat O_0$ from the fit with that of full post-selection at equal $N$, against the prediction $N\alpha_0^{1-\bar\rho^2}$;
3. add the noise-scaled runs and the $G^2$ extrapolation of Sec.~3(a), and check that the detector-measured $G$ removes the amplification uncertainty;
4. repeat with a non-Clifford logical rotation to measure the breakdown of the product structure.

## 8. Relation to Observable-Ranked Postselection (Froland et al., arXiv:2607.24947)

The reference experiment faces exactly the problem above and answers it with **ORP**: for each detector $v_j$ compute $\Delta_j=E[Q\mid v_j{=}0]-E[Q\mid v_j{=}1]$ from the data (their Eq. 5, 15), rank detectors by $\Delta_j$, keep only shots in which none of the top-$k$ detectors fired, and choose $k^\ast$ by a plateau-finding heuristic on a held-out half of the shots (Methods A 2); blockwise terminal parity post-selection is always applied. Their Appendix A 4 shows that at weak noise $\Delta_j\approx2\langle O\rangle\sum_{m\ni j}p_mq_m/\sum_{m\ni j}p_m$: the firing rate of the processes that fire $v_j$, weighted by their probability $q_m$ of flipping $O$ — in our notation $\Delta_j=2\langle O\rangle\kappa_j$. They also define an additive linear-regression metric $\beta_j$ (their Eq. 16, $Q_i=\beta_0+\sum_j\beta_jX_{ij}+\epsilon_i$) and use it only as an alternative ranking. Their Appendix A 5 shows that the detectors with large $\Delta_j$ coincide with the backward light cone of $O$.

The geometric law is the exact form of what ORP measures, and it suggests replacing the ranking-plus-cut by a regression:

**Per-detector log-linear law.** Let $v\in\{0,1\}^{N_d}$ be the detector vector. If every fired detector is attributable to a distinct fault (the sparse regime ORP also assumes), Theorem 1 applied detector by detector gives
$$E[Q\mid v]\;=\;O_0\prod_{j}\rho_j^{\,v_j},\qquad \rho_j=\frac{E[Q\mid v_j{=}1]}{E[Q\mid v_j{=}0]}=1-\frac{\Delta_j}{E[Q\mid v_j{=}0]}=1-2\kappa_j ,$$
i.e. $\ln E[Q\mid v]=\ln O_0+\sum_jv_j\ln\rho_j$. The additive model of their Eq. (16) is the first-order expansion of this log-linear law ($\beta_j\approx O_0(\rho_j-1)=-\Delta_j$ at leading order), which is why $\beta_j$ and $\Delta_j$ rank alike. Faults that fire a detector pair (a measurement fault fires $v_j$ and $v_{j+2}$ in their no-reset construction) are accommodated because the pair's log-damage can be split between its two detectors; what the model cannot represent is a detector fired by two fault classes with different damage, which is a second-order effect at the firing rates of the experiment.

**Estimator ("soft ORP").** Fit $(\ln O_0,\{\ln\rho_j\})$ by maximum likelihood on all shots ($Q_i=\pm1$ with mean $O_0\prod_j\rho_j^{v_{ij}}$), or in two stages: $\rho_j$ from the detector-wise conditional means (which ORP already computes), then $O_0$ from the weighted mean $\hat O_0=\sum_iw_iQ_i/\sum_iw_i\prod_j\rho_j^{v_{ij}}$ with $w_i=\prod_j\rho_j^{v_{ij}}$, which is the Fisher-optimal weighting at known $\rho$. Detectors outside the light cone have $\rho_j=1$ and their firing costs nothing; detectors inside are down-weighted by $\rho_j^2$ rather than cut. There is no cut level to choose and no plateau to find; the shots ORP discards at level $k^\ast$ contribute with weight $\rho^2$ instead of 0, and the shots it keeps are used identically. The effective shot count is $N\,E\!\left[\prod_j\rho_j^{2v_j}\right]=N\exp[-\sum_j p_j(1-\rho_j^2)]$ with $p_j$ the firing probabilities, against $N\exp[-\sum_{j\le k^\ast}p_j]$ for the cut; the two agree only when every selected detector has $\rho_j\approx0$ (a fired detector means a flipped observable), and the regression wins by $\exp[\sum_jp_j\rho_j^2]$ otherwise. For the light-cone detectors of a local observable $\rho_j$ is expected between about $0.5$ and $0.9$ (a fired detector in the light cone flips a local observable with probability $\kappa_j$ well below one, since most of the processes firing it act on other qubits of the block or commute with $O$), so most of the light-cone shots are retained at substantial weight; the actual values are read off the data as $1-\Delta_j/E[Q\mid v_j{=}0]$.

**Consistency check that comes for free.** Under the law, $\ln E[Q\mid |v|=k]$ for the light-cone detectors is linear in $k$ with slope $\overline{\ln\rho}$; curvature signals either overlapping faults (two faults in one window) or breakdown of the product structure for the non-Clifford rotations (Sec. 6, assumption 2). This replaces the plateau criterion by a model check.

**What stays the same.** The plateau level "set by the undetectable errors" in ORP is $O_0$ here; neither method reaches $O_{\rm ideal}$ without the second step of Sec. 3. The blockwise terminal parity post-selection is a detector like any other and enters the product.

**Proposed reanalysis.** The estimator needs only the shot-level $(v_i,Q_i)$ records that ORP already uses, so it can be run on the existing ibm_boston data set without new experiments: compare $\hat O_0$ and its bootstrap error from the regression with $\langle O\rangle_{k^\ast}$ and $\sigma_{k^\ast}$ for the observables of their Figs. 2–4, and check the linearity diagnostic. The expected outcome is the same central values with error bars reduced by the factor $\sqrt{\exp[\sum_{j\le k^\ast}p_j\rho_j^2]}$, which at their late-time acceptance fractions of a few per cent is a factor of several.
