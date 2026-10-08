# What to take from "Quantum error mitigation from information dynamics" (arXiv:2609.31934)

Reading note, 7 October 2026. Suchsland, Schuster, Rudolph, Noll, O'Brien, Minev (Google Quantum AI). Companion to `THEORY.md`; the numerical check is `../nonclifford/zne_below.py`, output in `zne_below.txt`.

## 1. What the paper does

- **Reactivity function.** Expand the expectation value in Pauli paths and group paths by their noise weight $w$: $C(\gamma)=\sum_w e^{-\gamma w}R(w)$. The noisy signal is the Laplace transform of $R(w)$, so the whole noise response of a circuit and observable is encoded in one function. For general Pauli noise $w$ becomes a real "exposure" $v$ (sum of Lindblad eigenvalues along the path) and the picture is unchanged.
- **Effective reactivity.** Finite shot precision makes the inverse Laplace transform ill-posed, so only features of $R(w)$ coarser than $1/\gamma$ are resolvable. In the regime where mitigation is affordable ($\gamma w^*\lesssim5$) that leaves $O(1)$ features, which is why fits with 2–8 parameters work.
- **Pauli-path ZNE.** Fit a physically motivated ansatz for $R(w)$ to data at amplified noise: a cumulant expansion (polynomial in $\gamma$ for $\ln C$) or a sum of a few exponentials. Single-exponential ZNE is the special case $R\approx\delta(w-W)$, which holds for Clifford circuits and little else; its bias comes from the width of $R$. Model order is chosen by convergence, with a telescoping-sum or AIC combination of orders.
- **Discrete insertions beat amplification.** Inserting $k$ noise events damps weight $w$ by $(1-4w/3V)^k$ instead of $e^{-\gamma w}$; amplification data are a Poisson mixture of insertion data (their Eq. 14), so insertions carry strictly more information and never cost more.
- **Filter functions and tunable error cancellation (TEC).** Any linear mitigation estimator $\hat C=\sum_r h_rC(\gamma_0 r)$ or $\sum_k h_kC(k)$ acts as a filter $h(w)$ on the reactivity; the bias is $\sum_w[1-e^{-\gamma_0w}h(w)]R(w)$ and the cost is $(\sum|h|)^2$. PEC is the filter $h=e^{\gamma_0w}$ at cost $e^{4K}$ ($K$ = mean number of fault events); Richardson ZNE matches $e^{\gamma_0w}$ to order $m$ in $\gamma_0w$. TEC is the $m\to\infty$ limit with tilted-Chebyshev spacing: it cancels noise exactly up to a threshold weight $w^*$ at cost $e^{4\gamma_0w^*}$, and interpolates between RZNE and PEC as $w^*$ grows. Their Fig. 14 is the useful picture: PEC is a signed combination of $k$-insertion experiments that cancels every non-zero net-error count and isolates the zero-error sector.

Across OTOC, Stark-MBL and 2D-TFIM circuits they remove ZNE biases up to 30% and reduce cost relative to PEC by up to $10^3$ at $10^{-2}$ RMSE.

## 2. Dictionary: our objects in their framework

**Twirled checks (THEORY.md, Section 4).** With acceptance rule $\mathcal A$ and $u=\beta_{\mathcal A}/\alpha_{\mathcal A}$,

$$
C_{\mathcal A}=\sum_w\big[(1-u)+u\,e^{-\gamma_0w}\big]R(w)\quad\Longleftrightarrow\quad h_{\mathcal A}(w)=(1-u)\,e^{\gamma_0w}+u .
$$

The leak subtraction $(\alpha C_{PS}-\beta C_{all})/(\alpha-\beta)$ is the linear combination with $h(w)\equiv e^{\gamma_0w}$: the PEC filter, exact at every weight, at cost $1/\alpha\approx e^{K}$ instead of $e^{4K}$. In the language of their Fig. 14, PEC isolates the zero-error sector by signed interference of insertion experiments; a twirled check tags every non-zero net error with a fair coin and the subtraction removes the tagged fraction by counting. Same filter, exponent reduced by four, paid in qubits and gadget gates instead of shots.

**Fixed checks (THEORY.md, Section 3; REPORT.md, Sections 4–7).** Post-selecting on a fixed check set $T$ removes the damping of each path by the faults $T$ detects. Writing a path's exposure to detected and undetected faults as $\omega_d(P)$ and $\omega_u(P)$,

$$
C_0=\sum_PA_P,\qquad C_{all}=\sum_PA_Pe^{-\omega_d(P)-\omega_u(P)},\qquad C_T=\sum_PA_Pe^{-\omega_u(P)}+O(\text{pairs}),
$$

a two-variable reactivity $R(\omega_d,\omega_u)$. The two-point extrapolation $C_T^{1+r}/C_{all}^{\,r}$ assumes $\omega_u=r\,\omega_d$ path by path. In the cumulant language of their Section II.B, $\ln C_T-\ln C_0=-\langle\omega_u\rangle+\tfrac12\mathrm{Var}(\omega_u)-\dots$ and $\ln C_{all}-\ln C_T=-\langle\omega_d\rangle+\tfrac12\mathrm{Var}(\omega_d)+\mathrm{Cov}(\omega_d,\omega_u)-\dots$, averages taken over the reactivity-weighted (signed) path distribution. The first-order ratio $r=\langle\omega_u\rangle/\langle\omega_d\rangle$ drops the second cumulants, which is why it failed in the non-Clifford test of `REPORT.md` Section 7; the ratio fitted on Clifford training circuits absorbs them, exactly as their fitted $\hat\kappa_m$ absorb the unresolved higher orders. The Clifford/Pauli-observable case is their $\delta$-function reactivity, where every ratio is exact.

**Two detection families.** A per-fault herald (each fault event flagged independently with probability $1-\sigma$, as in erasure conversion or a twirled check window around every single location) rescales every rate, $\lambda_\nu\to\sigma\lambda_\nu$: the ZNE family with gain $G=\sigma^k<1$. A single twirled check whose window covers everything gives the mixture family, linear in $u$. A general set of twirled windows interpolates between the two (the species of THEORY.md Section 9). So "error-detection ZNE" is one of two different things depending on the window structure: an exponential family entered *from below*, or a mixture family that is exact with two points.

## 3. What to take

### 3.1 ZNE from below

The bias of an order-$m$ fit extrapolated to $G=0$ scales with the product of the distances of the data points from zero. Heralded detection supplies points at $G<1$, so the extrapolation becomes short. Checked on the paper-angle Ising ring (`REPORT.md` Section 7) with exact rate scaling $G$ on every fault, the per-fault acceptance $\alpha(G)=\prod_\nu(1-(1-G)p_\nu)$ for $G<1$ and no overhead for $G>1$. Rms bias of $\langle X_i\rangle$, 5 sites; cost is $N\,\mathrm{Var}$ of the central site.

| estimator (points in $G$) | d=4, p=0.004 | d=4, p=0.012 | d=6, p=0.004 | d=6, p=0.012 | cost (d=6, p=0.012) |
|---|---|---|---|---|---|
| raw | 0.0381 | 0.1127 | 0.1680 | 0.4013 | 0.8 |
| exponential, above (1, 2) | 0.0033 | 0.0204 | 0.0055 | 0.0374 | 4.2 |
| exponential, below (1, ½) | 0.0009 | 0.0067 | 0.0015 | 0.0115 | 6.3 |
| exponential, below (1, ¼) | 0.0005 | 0.0035 | 0.0007 | 0.0060 | 2.9 |
| exponential, bracket (½, 2) | 0.0017 | 0.0112 | 0.0028 | 0.0203 | 2.5 |
| 2nd cumulant, above (1, 2, 3) | 0.0006 | 0.0088 | 0.0006 | 0.0136 | 16.8 |
| 2nd cumulant, below (1, ½, ¼) | 0.00002 | 0.0004 | 0.00001 | 0.0004 | 16.6 |
| 2nd cumulant, bracket (¼, 1, 3) | 0.0001 | 0.0015 | 0.0001 | 0.0019 | 3.5 |
| linear in $C$, below (1, ⅛) | 0.00004 | 0.0002 | 0.0028 | 0.0205 | 2.0 |

Same ansatz, same number of points: moving the points below the native noise reduces the bias 4–40×. The bracketed designs (one heralded point, the native point, one amplified point) keep most of the gain at a quarter of the cost, because the amplified point stabilises the slope. The plain linear fit in $C$ works at $d=4$ and fails at $d=6$, where the curvature of $C(G)$ near zero is large; this is the reactivity width showing, and the cumulant ansatz handles it.

Implications for check design: many short twirled windows (which is what mid-circuit check measurements give) are not a defect relative to one long window. They give the exponential-from-below family, whose extrapolation needs an ansatz but is well conditioned; one long window gives the mixture family, exact with two points but requiring full-weight checks. The acceptance $\alpha(G)$ measures $G$ in both cases.

### 3.2 Insert faults, do not amplify them

Our "undetected faults amplified" estimator in `REPORT.md` scales the undetected rates by $G=2,3$. The paper's Section I.E says to insert a known number $k$ of faults instead. For the undetected generator $\mathcal L_u=\sum_{\nu\in V_u}\lambda_\nu(\mathcal P_\nu-1)$ the relation is exact:

$$
e^{(G-1)\mathcal L_u}=\sum_{k\ge0}\frac{e^{-K}K^k}{k!}\,\mathcal D^k,\qquad K=(G-1)\Lambda_u,\quad \mathcal D=\frac1{\Lambda_u}\sum_{\nu\in V_u}\lambda_\nu\mathcal P_\nu ,
$$

so $C(k)$ (native noise plus $k$ Paulis drawn from $V_u$ with probabilities $\lambda_\nu/\Lambda_u$) determines the amplification data and not conversely. Insertions of undetected Paulis have zero syndrome, so acceptance is unchanged, as already verified for amplification. The fitting forms are their Eqs. 16 and 22 with the uniform $(1-4w/3V)$ replaced by the path's share of the insertion distribution. Not tested here; it is the natural replacement for the amplified variant, and the paper reports multi-exponential fits in $k$ outperforming fits in $\gamma$.

### 3.3 ED + TEC: a tunable version of ED + PEC

The anchor paper's first-order spacetime PEC inverts $\mathcal L_u$ exactly at cost $e^{4\Lambda_u}/\alpha$. In filter language that is the full PEC filter applied to the undetected sector. TEC restricted to $V_u$ cancels undetected noise up to a threshold weight $w^*$ at cost $e^{4\gamma_0w^*}/\alpha$, with $w^*$ as the convergence knob, and reduces to ED + Richardson ZNE at small $w^*$ and to ED + PEC at $w^*\gtrsim V_u$. The filter coefficients $h_k$ over $k$-insertion experiments can be taken from their analytic form or optimised numerically (a convex program). This is the most directly usable item for the grant's "syndrome-resolved mitigation" aim: it keeps the paper's syndrome classification and adds a tunable bias–cost trade-off on the part the checks cannot see. It also composes with the twirled-check mixture law: remove the visible sector exactly by leak subtraction, then apply TEC to the invisible sector.

### 3.4 An ansatz hierarchy for fixed-check extrapolation

The fixed-check subset family (`REPORT.md` Figure 1 and Section 7) currently has no convergence knob. The two-variable cumulant expansion of Section 2 supplies one: fit $\ln C_T$ over subsets $T$ as a polynomial in the detected exposure $-\ln\alpha_T$, order by order, and select the order by their telescoping rule (Eq. 19) or AIC (Eq. 20). The Clifford-trained ratio is the order-1 version with the slope fitted rather than computed; higher orders would use the training circuits to fit curvature as well. Their observation that only $O(1)$ features of the reactivity are resolvable explains why one pooled ratio sufficed in the non-Clifford test.

### 3.5 Reactivity as a check-placement criterion

Coverage in `REPORT.md` was counted in fault weight. The quantity that sets the residual bias is the reactivity-weighted exposure to invisible faults, $\langle\omega_u\rangle_R$, which is concentrated where the important Pauli paths live (the light cone of the observable, and for localised dynamics a small fraction of it, their Section III). Checks should be placed to cover those locations first; the visible fraction $f_O^d=\langle\omega_d\rangle/\langle\omega\rangle$ is measurable as the slope of $\ln C_T$ against $-\ln\alpha_T$ without any model. For twirled checks the question is moot inside the visible region and reduces to where to put the windows.

## 4. What does not transfer

- Their methods need amplification or insertion with a characterised Pauli model. The mixture law needs neither, but needs twirled checks and uniform visibility.
- The reactivity is signed (they allow complex cumulants), so the cumulant moments in Section 2 are not those of a probability distribution; convergence of the hierarchy is empirical, as in their work.
- Their non-unital and coherent extensions (Appendix A.2) do not cover leakage, and neither does anything here.

## 5. Suggested next steps, in order

1. Replace amplification of the undetected sector by $k$-insertion data and refit with their Eq. 22 (cheap in the existing simulator via the Poisson relation).
2. Implement ED + TEC on the Ising ring and map bias against cost as $w^*$ varies, next to ED + PEC and the trained-ratio extrapolation.
3. Add the order-2 term to the fixed-check subset extrapolation with order selection, and test on the random-angle ensemble.
4. For the twirled-check protocol, compare one long window against several short windows at equal acceptance, using the bracket design of Section 3.1 for the short-window case.
