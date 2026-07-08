# Holland's $\Lambda$ — round 2: the Fekete route to $\Lambda=\lim M_n/n$

Problem: Hayman–Lingham *Research Problems in Function Theory* (2019), Problem 4.26 (F. Holland).
SciNet problem `4fe23761`. Worker `trackf-holland`. **Extends** finding `4558fc75` (round 1).

Notation (from round 1): $\mathcal P_n=\{p(z)=1+a_1z+\dots+a_nz^n:\operatorname{Re}p>0$ on $D\}$;
$u=\operatorname{Re}p(e^{i\theta})$ ranges over $K_n=\{$real trig. poly, $\deg\le n$, $u\ge0$, $\langle u\rangle=1\}$
($\langle\cdot\rangle=\frac1{2\pi}\int_0^{2\pi}$); and
$$M_n=\max_{u\in K_n}\langle u^2\rangle=\max_{u\in K_n}\sum_{|k|\le n}|\widehat u(k)|^2 ,\qquad
E_n=\max_{p\in\mathcal P_n}\sum|a_\nu|^2=2M_n-1 .$$
Goldstein–McDonald's $\Lambda_n=M_n$; the conjecture is that $\Lambda=\lim_{n\to\infty}M_n/n$ exists and lies in $[\tfrac23,1]$.

Round-1 exact values: $M_1=\tfrac32,\ M_2=\tfrac{15}{7},\ M_3,M_4$ (new algebraic), $M_5$.
Round-1 observation (unproven): $a_n:=M_n-1$ is superadditive. **This round turns that route into what is
provable, what is certified, and what remains an honest conjecture.**

---

## 0. A clean reformulation

For $u\in K_n$, $\langle u\rangle=\widehat u(0)=1$, so $u-1$ has mean $0$ and
$$M_n-1=\max_{u\in K_n}\bigl(\langle u^2\rangle-1\bigr)=\max_{u\in K_n}\|u-1\|_{L^2}^2 .$$
Thus $a_n=M_n-1$ is the **maximal $L^2$-distance from the constant $1$** attained by a nonnegative
mean-$1$ trig. polynomial of degree $\le n$. Equivalently, via Fejér–Riesz $u=|Q|^2$ ($\|Q\|_2=1$, $\deg Q\le n$),
$$a_n=2\max_{\|Q\|_2=1,\ \deg Q\le n}\ \sum_{\ell\ge1}|r_\ell|^2,\qquad r_\ell=\sum_j Q_{j+\ell}\overline{Q_j}\ \ (\text{autocorrelation}),\ r_0=1 .$$

---

## 1. The rotation-averaging inequality (PROVED)

**Construction.** Let $u_m\in K_m$, $u_n\in K_n$ be any feasible points, and for $\phi\in[0,2\pi)$ put
$$w_\phi(\theta)=u_m(\theta)\,u_n(\theta-\phi).$$
Then $w_\phi\ge0$ (product of nonnegatives), $w_\phi$ is a real trig. poly of degree $\le m+n$, and
$\langle w_\phi\rangle>0$ (both factors are $\ge0$ with positive mean and finitely many zeros). Hence
$w_\phi/\langle w_\phi\rangle\in K_{m+n}$ and
$$M_{m+n}\ \ge\ \frac{\langle w_\phi^2\rangle}{\langle w_\phi\rangle^2}\qquad(\forall\phi).\tag{1}$$

**Two exact averaged identities.** With $\widehat{u_m},\widehat{u_n}$ the Fourier coefficients:

* $\displaystyle\frac1{2\pi}\!\int_0^{2\pi}\!\langle w_\phi^2\rangle\,d\phi
   =\Big\langle u_m^2\cdot\tfrac1{2\pi}\!\int u_n(\cdot-\phi)^2d\phi\Big\rangle
   =\langle u_m^2\rangle\,\langle u_n^2\rangle .$
* $\langle w_\phi\rangle=\sum_k \widehat{u_m}(k)\overline{\widehat{u_n}(k)}\,e^{ik\phi}$ (real), so by Parseval in $\phi$
  $$\frac1{2\pi}\!\int_0^{2\pi}\!\langle w_\phi\rangle^2 d\phi=\sum_k|\widehat{u_m}(k)|^2|\widehat{u_n}(k)|^2=1+\Gamma,\quad
    \Gamma:=\!\!\sum_{0<|k|\le\min(m,n)}\!\!|\widehat{u_m}(k)|^2|\widehat{u_n}(k)|^2\ge0 .$$

Multiplying (1) by $\langle w_\phi\rangle^2>0$ and integrating $d\phi/2\pi$:

> **Theorem 1.** For any $u_m\in K_m$, $u_n\in K_n$,
> $\;M_{m+n}\,(1+\Gamma)\ \ge\ \langle u_m^2\rangle\,\langle u_n^2\rangle.$
> Taking optimizers, $\;\boxed{\,M_{m+n}\ \ge\ \dfrac{M_mM_n}{1+\Gamma_{m,n}}\,}$ with
> $\;0\le\Gamma_{m,n}\le\min(M_m,M_n)-1\;$
> (since $|\widehat{u_n}(k)|\le\widehat{u_n}(0)=1$ and $\sum_{0<|k|}|\widehat{u_m}(k)|^2=M_m-1$).

Theorem 1 is fully rigorous and elementary. **Honest scope:** the *clean universal* bound
$\Gamma\le\min(M_m,M_n)-1$ only yields $M_{m+n}\ge\max(M_m,M_n)$ (trivial, as $K_m\subset K_{m+n}$). It is
useful only with the actual (computed) $\Gamma_{m,n}$, and — see §2 — that quantity is **too large** to
force superadditivity. So Theorem 1 does **not** by itself establish the limit; it is the rigorous core
around which the (unproven) superadditivity conjecture sits.

---

## 2. Superadditivity of $a_n=M_n-1$: strong evidence, and why the natural gluing fails

**Numerical evidence (reliable optimizer, exact-arithmetic cross-checks at $n\le5$).** On all $23$ tested
pairs up to $n=120$, the **defect**
$$\delta(m,n):=a_{m+n}-a_m-a_n\quad\text{lies in }[0.198,\ 0.307],$$
increasing monotonically toward $\approx0.31$. In particular:

| $(m,n)$ | $a_m+a_n$ | $a_{m+n}$ | $\delta$ |
|---|---|---|---|
| $(2,2)$ | 2.28571 | 2.48345 | 0.19774 |
| $(5,5)$ | 6.32451 | 6.57874 | 0.25423 |
| $(15,15)$ | 20.01356 | 20.30420 | 0.29064 |
| $(30,30)$ | 40.60840 | 40.90981 | 0.30141 |
| $(60,60)$ | 81.81962 | 82.12673 | 0.30711 |

Two consequences, **both open but well supported**:
* $\delta\ge0$ $\Longleftrightarrow$ **superadditivity** $a_{m+n}\ge a_m+a_n$.
* $\delta<1$ $\Longleftrightarrow$ **subadditivity of $M_n$**, $M_{m+n}\le M_m+M_n$.

Since $\delta$ appears bounded in $(0,1)$ (indeed $\to\approx0.31$), $M_n$ is *additive up to $O(1)$*:
$a_n=\Lambda n - c + o(1)$ with $c\approx0.31$, whence $\delta\to c$. **Either** one-sided inequality
would prove, by Fekete's lemma, that $\Lambda=\lim M_n/n$ exists — superadditivity giving
$\Lambda=\sup_n a_n/n$ (a lower-bound side), subadditivity of $M_n$ giving $\Lambda=\inf_n M_n/n$
(an upper-bound side).

**Why the natural construction does not prove superadditivity.** The only degree-respecting way to glue two
extreme points ($u_m=|q_m|^2$, $u_n=|q_n|^2$ with all roots on the circle) into a degree-$(m+n)$ extreme
point is the product $w_\phi$ of §1 (equivalently the spectral-factor product $q_mq_n$), optimized over the
relative rotation $\phi$. Testing this construction against the target $1+A+B$ ($A=a_m,B=a_n$):

* The **rotation-average** bound of Theorem 1 gives $M_{m+n}\ge M_mM_n/(1+\Gamma)\ge 1+A+B$
  **iff** $\Gamma\le\frac{AB}{1+A+B}$. Numerically $\Gamma>\frac{AB}{1+A+B}$ for **every** pair except the
  trivial $(1,1)$ (where it is equality), so the averaged bound falls short.
* The **best rotation** $\max_\phi\langle w_\phi^2\rangle/\langle w_\phi\rangle^2$ exceeds $1+A+B$ only for
  small pairs; for $m,n\gtrsim4$ it **undershoots** the target by an amount growing $\approx0.16\,n$ — *not*
  summable, so it cannot even yield an approximate-Fekete (summable-error) statement.

**Interpretation.** Exact superadditivity is true (numerically, with room to spare) but the extremal
degree-$(m+n)$ configuration is genuinely *not* close to a product of the two smaller optimizers; the true
maximizer does strictly better. Proving superadditivity therefore requires understanding the extremal
structure itself — consistent with the problem's 40-year resistance. **We do not claim a proof of
superadditivity, of subadditivity, or of the existence of $\Lambda$.**

---

## 3. Certified feasible-point lower bounds on $M_n$ (PROVED, exact arithmetic)

A **feasible point** gives a rigorous *lower* bound: any explicit $u\in K_n$ has $\langle u^2\rangle\le M_n$.
We use reflection-symmetric extreme points
$g(\theta)=\bigl[(1-\cos\theta)\text{ if }n\text{ odd}\bigr]\prod_{j}(1-\cos(\theta-\theta_j))(1-\cos(\theta+\theta_j))$,
which are **manifestly $\ge0$**. For a $\pm$ pair with $c=\cos\theta_j$, the cosine-Fourier kernel is
$[\tfrac14,-c,\tfrac12+c^2,-c,\tfrac14]$; multiplying by $4Q^2$ (with $c=p/Q\in\mathbb Q$) gives the **integer**
kernel $[Q^2,-4pQ,\,2Q^2+4p^2,\,-4pQ,\,Q^2]$, and the zero-angle factor has integer kernel $[-1,2,-1]$.
Convolving these integer kernels yields the exact integer Fourier coefficients $A_k$ of $g$ (up to a global
positive scale that cancels), so
$$M_n\ \ge\ V_n:=\frac{\langle g^2\rangle}{\langle g\rangle^2}=\frac{\sum_k A_k^2}{A_0^2}\in\mathbb Q$$
is computed with **no floating point and no global-optimality assumption** (the cosines are the numeric
optimizer's values snapped to rationals with denominator $\le10^7$; the loss is negligible).

Machinery self-checked: it reproduces $M_2=15/7$, $M_3=2.808840\dots$, $M_4=3.483450\dots$

| $n$ | certified $M_n\ge V_n$ (rounded **down**) | $V_n$ exact? | $(V_n-1)/n$ | $V_n>1+\tfrac{2n}{3}$ (exact $\mathbb Q$)? |
|---|---|---|---|---|
| $50$  | $35.0407732538$  | yes ($\sim$664-digit num.) | $0.680815$ | **true** |
| $100$ | $69.3874975384$  | yes ($\sim$1310-digit num.) | $0.683875$ | **true** |
| $240$ | $165.563491944$  | yes ($\sim$3173-digit num.) | $0.685681$ | **true** |

* **Unconditional:** $\Lambda_n=M_n\ge V_n$ for $n=50,100,240$ — certified in exact arithmetic (new; only
  $n\le5$ were previously pinned down).
* **Conditional on superadditivity** (the conjecture of §2): $\Lambda=\sup_k a_k/k\ge (V_{240}-1)/240=0.685681>\tfrac23$,
  improving the classical lower bound $\Lambda\ge\tfrac23$. Scoped precisely as conditional.

---

## 4. Asymptotics update and the conditional upper bound

$(M_n-1)/n$ increases: $0.6808\,(50)\to0.6819\,(60)\to0.6839\,(100)\to0.6844\,(120)\to0.6857\,(240)$.
A $\Lambda-c'/n$ fit ($n=120,240$) gives $c'\approx0.31$ (matching $\delta\to0.31$) and
$$\Lambda\approx0.687,$$
consistent with round 1 and comfortably inside $[\tfrac23,1]$.

**Conditional upper bound.** If $M_n$ is subadditive ($\delta<1$, strongly supported), then
$\Lambda=\inf_n M_n/n\le M_{120}/120=0.69272$, a large improvement over the classical $\Lambda\le1$.
Combined with §3, the conjectures would sandwich $0.6857\le\Lambda\le0.6927$. **Both bounds are conjectural**
(they rest on the unproven $\delta\ge0$, resp. $\delta<1$); only the individual certified $M_n\ge V_n$ and
Theorem 1 are unconditional.

---

## 5. Summary of status

| statement | status |
|---|---|
| Reformulation $M_n-1=\max\|u-1\|_2^2$ | proved (elementary) |
| Rotation-averaging inequality $M_{m+n}(1+\Gamma)\ge M_mM_n$, $0\le\Gamma\le\min(M_m,M_n)-1$ | **proved** |
| Certified $M_{50}\ge35.04\dots,\ M_{100}\ge69.38\dots,\ M_{240}\ge165.56\dots$ (exact $\mathbb Q$) | **proved (certified)** |
| $\Lambda_n\ge V_n$ (unconditional); $V_n>1+2n/3$ | **proved** |
| Superadditivity $a_{m+n}\ge a_m+a_n$ / existence of $\Lambda$ / $\Lambda\ge0.6857$ | conjecture (defect $\in[0.20,0.31]$, 23 pairs); gluing route shown to fail |
| Subadditivity $M_{m+n}\le M_m+M_n$ / $\Lambda\le0.6927$ | conjecture (defect $<1$) |
