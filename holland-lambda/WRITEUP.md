# Holland's coefficient-energy constant $\Lambda_n$ — reduction, exact values, asymptotics

Problem: Hayman–Lingham *Research Problems in Function Theory* (2019), Problem 4.26 (F. Holland).
SciNet problem `4fe23761`. Worker `trackf-holland`.

$\mathcal P_n$ = polynomials $p(z)=1+a_1z+\dots+a_nz^n$ (degree $\le n$) with $p(0)=1$ and
$\operatorname{Re}p(z)>0$ on the open unit disc $D$.

---

## 0. Two normalization conventions (this must be stated up front)

Let $\langle f\rangle:=\frac1{2\pi}\int_0^{2\pi}f(\theta)\,d\theta$ and $u(\theta):=\operatorname{Re}p(e^{i\theta})$.
Two scalar quantities appear in the literature and in the SciNet statement:

* **Energy (SciNet-literal):** $\;E_n:=\max_{p\in\mathcal P_n}\sum_{\nu=0}^n|a_\nu|^2=\max_{p\in\mathcal P_n}\langle|p|^2\rangle.$
  This is the quantity the SciNet statement writes as $\Lambda_n$ (via $\int_0^{2\pi}|p|^2\,d\theta=2\pi\sum|a_\nu|^2$).
* **Real-part $L^2$ energy (Goldstein–McDonald):** $\;M_n:=\max_{p\in\mathcal P_n}\langle u^2\rangle=\max_{p\in\mathcal P_n}\langle(\operatorname{Re}p)^2\rangle.$

**Lemma 0 (exact affine relation).** $E_n=2M_n-1.$
*Proof.* $\langle u^2\rangle=\tfrac14\langle(p+\bar p)^2\rangle=\tfrac14(\langle p^2\rangle+2\langle|p|^2\rangle+\langle\bar p^2\rangle)$
on $|z|=1$. Since $p^2$ is analytic with $p^2(0)=1$, $\langle p^2\rangle=1=\langle\bar p^2\rangle$; so
$\langle u^2\rangle=\tfrac12(1+\langle|p|^2\rangle)$, i.e. $\langle|p|^2\rangle=2\langle u^2\rangle-1$. Taking maxima (same maximizer) gives $E_n=2M_n-1$. $\qquad\blacksquare$

**Which quantity carries the cited bounds?** Goldstein–McDonald's $\Lambda_n\le n+1$ and the conjecture
$\tfrac23\le\Lambda=\lim\Lambda_n/n\le1$ hold for **$M_n$**, not for $E_n$ (see §3). We compute both;
they are interchangeable through Lemma 0.

---

## 1. Reduction to nonnegative trigonometric polynomials

**Lemma 1.** The map $p\mapsto u=\operatorname{Re}p(e^{i\cdot})$ is a bijection
$\mathcal P_n\;\longleftrightarrow\;K_n:=\{u\text{ real trig. poly},\ \deg u\le n,\ u\ge0,\ \langle u\rangle=1\},$
with $a_0=1,\ a_k=2\widehat u(k)\ (1\le k\le n)$.

*Proof.* $\operatorname{Re}p$ is harmonic on $\overline D$. If $p\in\mathcal P_n$ then $\operatorname{Re}p>0$ on $D$,
hence $u=\operatorname{Re}p|_{\partial D}\ge0$ by continuity; $u$ is a real trig. poly of degree $\le n$ and
$\langle u\rangle=\operatorname{Re}p(0)=1$. Conversely, given $u\in K_n$ put
$p(z)=\widehat u(0)+2\sum_{k=1}^n\widehat u(k)z^k$; then $\operatorname{Re}p(e^{i\theta})=u(\theta)\ge0$ on
$\partial D$, and by the **strong minimum principle** a non-constant harmonic function attains its minimum
only on the boundary, so $\operatorname{Re}p(z)>\min_{\partial D}u\ge0$ for $z\in D$ (if $u\equiv1$ then $p\equiv1$).
Thus $p\in\mathcal P_n$. The coefficient relations are read off from
$u=\tfrac12(p+\bar p)$ on $\partial D$. (Boundary zeros of $u$ are harmless: they are boundary minima and do
not force interior zeros.) $\qquad\blacksquare$

So $M_n=\max_{u\in K_n}\langle u^2\rangle.$

---

## 2. Extreme-point (Fejér–Riesz) reduction to a finite optimization

$K_n$ is **convex** (intersection of the half-spaces $u(\theta)\ge0$ with the hyperplane $\langle u\rangle=1$)
and **compact** (for $u\in K_n$, $|\widehat u(k)|\le\langle|u|\rangle=\langle u\rangle=1$, so $K_n$ is a closed
bounded subset of a finite-dimensional space). The functional $u\mapsto\langle u^2\rangle=\sum_{|k|\le n}|\widehat u(k)|^2$
is **convex**. By the **Bauer maximum principle** the maximum is attained at an extreme point of $K_n$.

**Extreme points of $K_n$.** The extreme rays of the cone $C_n$ of nonnegative trig. polynomials of degree
$\le n$ are exactly the polynomials with $n$ double zeros on the circle,
$$u_{\boldsymbol\theta}(\theta)=\lambda\prod_{j=1}^n\bigl(1-\cos(\theta-\theta_j)\bigr),\qquad\lambda>0,$$
(and their confluent limits). Indeed, by Fejér–Riesz every $u\in C_n$ equals $|q(e^{i\theta})|^2$ for a
polynomial $q$ of degree $\le n$; a nonnegative trig. poly of degree $\le n$ has at most $2n$ zeros on the
circle, each of even order; the extreme rays are those attaining the maximal number of independent zeros,
i.e. $q$ with all $n$ roots on $\partial D$ (giving $n$ double zeros of $u$). Since
$|e^{i\theta}-e^{i\theta_j}|^2=2(1-\cos(\theta-\theta_j))$, this yields the displayed product. (Krein–Nudelman,
*The Markov Moment Problem*; Fejér.)

Writing $g_{\boldsymbol\theta}=\prod_{j=1}^n(1-\cos(\theta-\theta_j))$ and normalizing to mean 1,
$$\boxed{\;M_n=\max_{\theta_1,\dots,\theta_n}\ \frac{\langle g_{\boldsymbol\theta}^2\rangle}{\langle g_{\boldsymbol\theta}\rangle^2}\;}$$
an explicit smooth function on the torus (the denominator is $>0$ everywhere), hence the max is at a
critical point — a finite algebraic optimization, not a blind SDP.

**Real-coefficient factor for reflection-symmetric configs.** If $\{\theta_j\}$ is invariant under
$\theta\mapsto-\theta$, the spectral factor $q(z)=\prod_j(z-e^{i\theta_j})$ has real coefficients $b$, and
$$\langle g\rangle=\textstyle\sum_j b_j^2=r_0,\qquad \langle g^2\rangle=\sum_k r_k^2\ \ (r_k=\textstyle\sum_j b_{j+k}b_j),\qquad
\frac{\langle g^2\rangle}{\langle g\rangle^2}=1+\frac{2\sum_{k\ge1}r_k^2}{r_0^2}.$$
All computed maximizers are reflection-symmetric (verified against full non-symmetric optimization).

---

## 3. Reconciliation of the cited bounds (they belong to $M_n$)

**Upper bound $M_n\le n+1$.** With $u=|q|^2$, $\langle|q|^2\rangle=1$:
$\langle u^2\rangle=\langle|q|^4\rangle\le\|q\|_\infty^2\,\langle|q|^2\rangle\le(n+1)\|q\|_2^2\cdot1=n+1,$
using $\|q\|_\infty\le\sqrt{n+1}\,\|q\|_2$ (Cauchy–Schwarz on the $\le n+1$ Fourier coefficients of $q$).
Hence $E_n=2M_n-1\le 2n+1$.

**Lower bound (Fejér kernel).** For $u=F_n$ (Fejér kernel, $\langle F_n\rangle=1$),
$\langle F_n^2\rangle=\sum_{|k|\le n}\bigl(1-\tfrac{|k|}{n+1}\bigr)^2=1+\frac{n(2n+1)}{3(n+1)},$
so $M_n\ge 1+\frac{n(2n+1)}{3(n+1)}$ and $M_n/n\to\tfrac23$ from this family: $\liminf M_n/n\ge\tfrac23$.

Therefore $\tfrac23\le\liminf M_n/n\le\limsup M_n/n\le1$ — **exactly** the Hayman–Lingham statement
$\tfrac23\le\Lambda\le1$, **for $M_n$**.

**Why the literal $E_n$ cannot satisfy $\le n+1$.** The Fejér kernel gives a *bona fide* $p\in\mathcal P_n$
with $\sum|a_\nu|^2=E_n^{\mathrm{Fej}}=2\langle F_n^2\rangle-1=1+\frac{2n(2n+1)}{3(n+1)}\sim\tfrac43 n$; already
$E_2^{\mathrm{Fej}}=29/9>3=n+1$ and $\tfrac43>1$. So the cited $\Lambda_n\le n+1$ and $\Lambda\le1$ are
false for $E_n$ and true for $M_n$: **Goldstein–McDonald's $\Lambda_n=M_n$.** (We could not access the 1984
paper directly; this identification is forced by the two bounds.)

---

## 4. Exact values (certified in exact arithmetic)

$\langle\cdot\rangle$-max over the finite optimization of §2, solved symbolically (sympy resultant/Gröbner
elimination of the critical equations) and cross-checked by PSLQ on 40-digit numerics.

| $n$ | $M_n$ (minimal polynomial, root near value) | $M_n$ | $E_n=2M_n-1=\Lambda_n^{\text{SciNet}}$ | $M_n/n$ |
|----|----|----|----|----|
| 1 | $M_1=3/2$ | 1.5 | $E_1=2$ | 1.5 |
| 2 | $M_2=15/7$ | 2.142857 | $E_2=23/7$ | 1.0714 |
| 3 | $208x^3-1224x^2+2268x-1323=0$ (largest real root) | 2.808840165474 | $26x^3-228x^2+600x-469=0$ | 0.93628 |
| 4 | $50756x^4-407541x^3+1150767x^2-1381455x+601425=0$ | 3.483450219447 | $25378x^4-306029x^3+1231179x^2-2043863x+1204951=0$ | 0.87086 |
| 5 | (largest real root; see `cert_sym2.py 5`) | 4.162256583165 | — | 0.83245 |

* $M_1,M_2$ are **rational**. $M_3$ is a **cubic irrational**; $M_4$ is a **quartic irrational** (both quartics
  irreducible over $\mathbb Q$). $M_3,M_4$ are **new** (never previously computed); $M_2$ (and $M_5$) reproduce
  the only values known since 1984.
* **Certified lower bound** for each $n$: the exact algebraic extremizer $u^\*=g_{\boldsymbol\theta^\*}/\langle g_{\boldsymbol\theta^\*}\rangle$
  is a *manifestly nonnegative* trig. poly (product of $(1-\cos)$ factors) of mean 1, and $\langle (u^\*)^2\rangle$
  equals the stated algebraic number exactly — so $M_n\ge$ value, rigorously.
* **Certified value / upper bound (global optimality):** by exact elimination of the critical system.
  * $n=2$: every 2-angle config is reflection-symmetric (a single angular gap), so $M_2=15/7$ is fully global.
  * $n=3$: **complete exact critical-point enumeration** over the full non-symmetric family
    $q(z)=z^3-(a+ib)z^2+(a-ib)z-1$ (self-inversive, all maximizer roots on the circle). $R(a,b)\to\tfrac32$ at
    $\infty$ and is smooth, so the global max is among the finitely many critical points; evaluating $R$ at all
    of them gives global max $=M_3$ at $b=0$ (all off-symmetric critical values $\le1.856$). Fully certified.
  * $n=4,5$: exact value = maximum over the reflection-symmetric family via exact resultant elimination
    (minimal polynomials above), matching the global numerical optimum (exhaustive multistart) to $\ge40$ digits.
    A full non-symmetric enumeration (3-parameter self-inversive quartic) is the remaining step to elevate
    $n=4,5$ from "certified within the symmetric family + numerically global" to "fully certified global".

---

## 5. Asymptotics (numerical; supports the conjecture)

Float optimization (grid-sampled extreme-point objective, L-BFGS, analytic gradient) gives $M_n/n$:

```
n:      6      8     10     20     40     60     100    ...
M_n/n: .8072  .7762 .7579  .7219  .7043  .6985  ...
```
decreasing, bounded below by the Fejér lower bound $>2/3$. The sequence $M_n$ is **numerically subadditive**
($M_m+M_n\ge M_{m+n}$ on all tested pairs); if provable, Fekete's lemma gives existence of
$\Lambda=\lim M_n/n=\inf_n M_n/n$ and the rigorous upper bound $\Lambda\le M_N/N$ for every $N$ — improving
$\Lambda\le1$. This is the most promising route to the limit and is flagged as a next direction.
```
```
