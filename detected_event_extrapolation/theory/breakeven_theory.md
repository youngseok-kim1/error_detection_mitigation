# Break-even theory for detection-driven extrapolation

7 October 2026. Companion to `demo/THRESHOLD.md` (data) and `theory/gadget_noise_survey.md` (literature). Numerical check: `demo/breakeven_check.py` → `demo/breakeven_check.md`.

The question raised by the noise scan: is there a threshold, in the sense of quantum error correction, below which check-based extrapolation works at any depth? Short answer: **no threshold in the gate error alone exists, and none can in this class of methods; what exists is (i) a depth- and error-independent scaling condition on the gadget overhead and the observable, which decides whether checks or amplification win at large depth, and (ii) a finite-depth break-even at about one expected visible fault in the checked circuit, $\Lambda_{\rm tot}\approx1.2$–$1.5$, which is where the 1% results fail and the 0.2% results succeed.** The derivation below says why, and what would have to change to get a genuine threshold.

## 1. Setting and the three fault weights

A layered circuit of depth $d$ with $g_p$ noisy two-qubit payload gates per layer at error $p$, checked by gadgets with $g$ noisy two-qubit gates per layer, readout flip $q$ on the ancilla. Three fault weights (Pauli–Lindblad sums $\Lambda=\sum_\nu\lambda_\nu$) matter:

- $\Lambda_{\rm pay}=\tfrac{14}{15}\,g_p\,p\,d$ — the payload fault weight, what PEC has to cancel.
- $\Lambda_O=-\ln(O_{\rm raw}/O_{\rm ideal})\equiv\kappa\,\Lambda_{\rm pay}$ — the fault weight the observable actually feels; $\kappa\le1$ is the **observable sensitivity**. Faults that commute with the back-propagated observable, or land outside its light cone, do not count. For $m_x$ on the TFIM ring, $\kappa\approx0.45$ (measured: $\Lambda_O/\Lambda_{\rm tot}\approx0.15$–$0.26$ with $\Lambda_{\rm tot}=2.25\Lambda_{\rm pay}$).
- $\Lambda_{\rm tot}=\Lambda_{\rm vis,pay}+\Lambda_{\rm gad}=(1+r)\,\Lambda_{\rm pay}$ for the repaired twirl (every fault visible), with $r=g/g_p$ the **gadget overhead ratio**: $r=1.25$ for the repaired gadget, $0.75$ for the constrained family (where $\Lambda_{\rm vis,pay}=\tfrac{13}{15}$ of the payload weight). The clean fraction is $q_0=e^{-\Lambda_{\rm tot}}(1-q)^d$.

The key observation is that **post-selection is observable-agnostic**: it pays for every visible fault whether or not that fault would have changed $O$, while ZNE and PEC pay only for $\Lambda_O$ and $\Lambda_{\rm pay}$. This single fact fixes the whole break-even structure.

## 2. Cost exponents

Per-shot variance $\sigma^2$ of the raw outcome ($\approx0.2$ for the 5-site mean of $\pm1$ outcomes). All statements are leading order in the exponentials; $N$ is the total shot count.

**ZNE, two-point exponential (gains 1, 2).** $\hat O=O_1^2/O_2$ with $O_G=O_0e^{-G\Lambda_O}$ and $N/2$ shots per gain. Delta method:
$$N\,\mathrm{Var}(\hat O)\;\simeq\;2\sigma^2\bigl(4e^{2\Lambda_O}+e^{4\Lambda_O}\bigr).$$
Richardson with gains up to $G_{\max}$: $\propto e^{2G_{\max}\Lambda_O}$ (second cumulant, gains 1,2,3: $e^{6\Lambda_O}$ with a larger prefactor).

**PEC.** $N\,\mathrm{Var}\simeq\sigma^2e^{4\Lambda_{\rm pay}}$.

**Fixed checks + first-order PEC on the undetected faults.** $N\,\mathrm{Var}\simeq\sigma^2e^{4\Lambda_u}/\alpha$ (anchor paper).

**Twirled checks, leak subtraction / species inversion.** From the variance formula of THEORY.md with $\beta=\tfrac12$ and $\alpha-\beta=q_0(1-\beta)$: the numerator is $\simeq\sigma^2/4$ and the denominator $q_0^2/4$, so
$$N\,\mathrm{Var}(\hat O_\varnothing)\;\simeq\;\sigma^2\,A\,e^{2\Lambda_{\rm tot}},$$
with $A\ge1$ the species-inversion amplification ($A\approx1$ for $\Lambda_{\rm tot}\lesssim1$, $\approx2$–$3$ in variance at $\Lambda_{\rm tot}\approx2$, and unbounded at fixed $N$ once the inferred $\hat q_0$ is within a few standard errors of zero, $q_0\lesssim|c|\sqrt{\alpha/N}$).

Numerical check (`breakeven_check.md`, 13 points at 1%, 0.5%, 0.2%): with no fitted parameter the formulas reproduce the ZNE standard deviation within $1.6\times$ at every point (the formula ignores the finite-$N$ fit nonlinearity), the repaired-twirl standard deviation within $2\times$ for $\Lambda_{\rm tot}\le2.3$ (the excess is the species amplification $A\approx1.3$–$4$ in variance), and the ratio of the two within $1.3\times$ at eight of the eleven points with $q_0>0.09$ and within $1.7\times$ at the other three; beyond $q_0\approx0.09$ the heavy tail of $W_0/\hat q_0$ takes over, as the data show.

## 3. No threshold in $p$ alone

**Proposition 1.** Under layer-wise noise every cost exponent above is proportional to $p\,d$ (and $q\,d$): $\Lambda_O=\kappa\tfrac{14}{15}g_p\,pd$, $\Lambda_{\rm pay}=\tfrac{14}{15}g_p\,pd$, $\Lambda_{\rm tot}=(1+r)\tfrac{14}{15}g_p\,pd$. Hence for a fixed circuit family, observable and gadget, the ratio of the cost of any two methods is a function of $pd$ only, and the shot count needed by every method grows exponentially in $pd$.

This is the collapse of the scan onto $q_0$. Its consequences:

- At fixed shot budget the reachable depth of **every** method, ZNE and PEC included, scales as $d_{\max}\propto1/p$. Lowering the gate error buys depth linearly. That is the "scalable" sense in which the 0.2% results work where the 1% ones fail: nothing changed qualitatively, the break-even depth moved from 6 to 29.
- A QEC-type threshold is a statement that the *per-layer* logical error is suppressed, $p\to cp^2$, so that the exponent itself shrinks with $p$. Detection plus extrapolation does not suppress anything: the accepted set still contains faulty shots (a fraction $1-q_0/\alpha$, which is 84% at $1\%$, $d=16$), and the fault-free value is recovered by re-weighting. The per-layer visible error is $(1+r)\tfrac{14}{15}g_pp$ at every $p$, so the exponent cannot shrink. Rejection-based detection (QED) has the same property: its cost is $1/q_0=e^{\Lambda_{\rm tot}}$, exponential in $pd$, with no threshold.
- The one condition that is in $p$ alone is trivial here: the coin must stay informative, which needs the gadget's own fault probability per check below $\tfrac12$ (van den Berg's $t_{ok}>\tfrac12$), i.e. $gp\lesssim0.5$, $p\lesssim4\%$ for $g=12.5$.

## 4. The scaling condition: a threshold on the gadget, not on $p$

What *is* depth- and error-independent is the comparison of exponents.

**Proposition 2.** The twirled-check cost grows more slowly with depth than that of ZNE with maximum gain $G_{\max}$ if and only if
$$2(1+r)\,\Lambda_{\rm pay}\;<\;2G_{\max}\,\kappa\,\Lambda_{\rm pay}\quad\Longleftrightarrow\quad 1+r\;<\;G_{\max}\,\kappa .$$
Against PEC the condition is $2(1+r)<4$, i.e. $r<1$ (always true for the constrained family, false for the repaired gadget).

For the TFIM demo, $\kappa\approx0.45$, so $G_{\max}\kappa\approx0.9$ (two-point exponential) or $1.35$ (second cumulant), against $1+r=2.25$ (repaired) or $1.75$ (constrained). **ZNE wins asymptotically for every gadget we have, and would still win with a free gadget ($r=0$), because $\kappa<1/G_{\max}$: the observable feels fewer than half of the faults the checks reject.** The check-based protocols can therefore only win at finite depth, through the prefactor, which is what the scan shows: the ratio $\mathrm{std}_{\rm rep}/\mathrm{std}_{\rm ZNE}$ rises monotonically with $q_0^{-1}$ and crosses 1 at $q_0\approx0.3$. Against second-cumulant ZNE the slope is smaller and the prefactor larger, so the crossing sits at $q_0\approx0.1$, as observed.

## 5. The finite-depth break-even

Setting the two costs equal with $\Lambda_O=\tfrac{\kappa}{1+r}\Lambda_{\rm tot}\approx0.2\,\Lambda_{\rm tot}$:
$$A\,e^{2\Lambda_{\rm tot}}=2\bigl(4e^{0.4\Lambda_{\rm tot}}+e^{0.8\Lambda_{\rm tot}}\bigr)\;\Rightarrow\;\Lambda_{\rm tot}^\ast\approx1.5\ (A=1),\ \ \approx1.2\ (A=2),$$
i.e. $q_0^\ast\approx0.22$–$0.30$. The scan gives $q_0^\ast\approx0.3$, $\ln(1/q_0^\ast)\approx1.2$.

**Symbols.** $d$ = number of checked layers (Trotter steps); $g_p$ = noisy two-qubit payload gates per layer (10 here); $g$ = noisy two-qubit gadget gates per layer, averaged over the flavor distribution (controlled-Paulis of both check sides, plus controlled-$X$ repair gates for the repaired gadget: $g=7.5$ constrained, $12.5$ repaired); $p$ = two-qubit depolarizing probability per gate; $q$ = ancilla readout flip probability, one readout per layer.

**Repaired gadget** (every non-identity Pauli visible, probability $\tfrac{14}{15}p$ per gate): $q_0=(1-\tfrac{14}{15}p)^{(g_p+g)d}(1-q)^d$, so $q_0\ge q_0^\ast$ is
$$\tfrac{14}{15}(g_p+g)\,p\,d+q\,d\;\lesssim\;\ln(1/q_0^\ast)\approx1.2,$$
the expected number of flag-producing events (visible faults + readout flips) in the whole checked circuit. The readout enters through $(1-q)^d$: a readout flip rejects a clean shot exactly like a fault. With $q=p$ it is 1 unit against 21 for the gates (5%). The rounder form $(g_p+g)pd\lesssim1.3$ quoted in THRESHOLD.md drops the readout term and absorbs the $15/14$ into the constant. Numerically $pd\lesssim0.055$: $d\lesssim5$–6 at 1%, 11 at 0.5%, 27 at 0.2%, 55 at 0.1%.

**Constrained family + PEC(inv)**: $g=7.5$; only $\tfrac{13}{15}$ of payload faults are visible; the PEC factor $e^{4\Lambda_{\rm inv}}$, $\Lambda_{\rm inv}=\tfrac{2}{15}g_ppd$, shifts $\ln(1/q_0^\ast)$ by $-2\Lambda_{\rm inv}$ (5%). Condition: $(\tfrac{13}{15}g_p+\tfrac{14}{15}g)pd+qd\approx16.7\,pd\lesssim1.2$, $pd\lesssim0.07$: $d\lesssim7$ at 1% (measured crossing vs exponential ZNE between $d=8$ and 10), 14 at 0.5%, 36 at 0.2% — about 30% more depth than the repaired gadget at the same error, for one learned rate per gate.

Below the break-even the check-based estimate is both cheaper and unbiased. Above it ZNE is cheaper but carries its model bias, which is $O(\Lambda_O^2)$ for the exponential fit and not a function of $pd$ alone (it depends on how non-exponential the decay is: $-0.028$ at $1\%$, $d=12$; $+0.020$ at $0.5\%$, $d=16$; harmless at $0.2\%$). The check-based estimate has no such term, the mixture law being exact, so the regime above break-even is "ZNE if its bias is tolerable, checks otherwise".

## 6. What would give a genuine threshold-like advantage

Proposition 2 says the lever is not $p$ but $\kappa$ and $r$. Two routes:

1. **Observable-targeted checks: make $\Lambda_{\rm vis}\approx\Lambda_O$.** Place checks only on the backward light cone of $O$ and, within it, choose flavors that see the faults $O$ is sensitive to (for $m_x$, faults that anticommute with the propagated $X$-strings). Then the exponent becomes $2(1+r')\Lambda_O$ with $r'$ the gadget-to-protected-payload ratio, and the condition is $1+r'<G_{\max}$: checks beat two-point exponential ZNE asymptotically whenever $r'<1$. Single-sided checks ($r'\approx0.4$) and dual-purpose ancillas measured once ($r'\approx0$, Fischer et al.) satisfy it; the two-sided repaired gadget ($r'=1.25$) does not. This is a prediction to test, and the reason the anchor paper's gate-efficient, ancilla-local checks are the right design.
2. **Reduce $r$ at fixed coverage.** A flavor distribution that keeps the number of repaired bonds small, Z-string-heavy flavors (2.5 gates instead of 10), and shared ancillas all lower $r$ without changing $\kappa$. For the constrained family $r=0.75$ already satisfies the PEC condition $r<1$; the comparison to PEC is then won asymptotically with no noise model.

Neither route produces a threshold in $p$; both move the exponent, which is the only thing that matters at depth. The honest summary for the paper: detection-driven extrapolation is unbiased at every noise level; it is the cheapest method up to about one expected visible fault, $(g_p+g)pd\lesssim1.3$, a depth that scales as $1/p$; and it scales better than ZNE at large depth only if the checks are targeted so that the gadget overhead per protected fault satisfies $1+r<G_{\max}\kappa$.
