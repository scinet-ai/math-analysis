# Rippon 7.54: structural theorems, a verified certificate, and proof-route analysis

**Problem (Hayman–Lingham, Research Problems in Function Theory, 50th anniv. ed., Problem
7.54, posed by P. J. Rippon).** Let $\varphi_t(z)=e^{tz}-1$ and
$\varphi_t^{1}=\varphi_t,\ \varphi_t^{n+1}=\varphi_t\circ\varphi_t^{n}$. Each
$\varphi_t^{n}(-1)$ is a formal power series in $t$, and by the substitution identity

$$\varphi_t^{n+1}(-1)=e^{t\,\varphi_t^{n}(-1)}-1 .$$

Write $f_n(t):=\varphi_t^{n}(-1)=\sum_{k\ge 0}c^{(n)}_k t^k$ and $c^{(n)}_k:=[t^k]f_n(t)$.
**Question.** Is $\bigl|c^{(n)}_k\bigr|\le 1$ for every $n\ge1$ and $k\ge0$?

The book's Update reads, in full: *"No progress on this problem has been reported to us."*
This note records the first structural progress: three exact lemmas, a reduction of the
conjecture to two explicit sub-statements, a dual-method (exact + interval) verified
certificate over a large window, a precise extremal-structure map, and an analysis of
three proof routes.

---

## 0. Preliminaries: valuations

$f_1=e^{-t}-1$ has valuation $1$. If $f_n$ has valuation $n$ then $t f_n$ has valuation
$n+1$ and $f_{n+1}=e^{tf_n}-1=\sum_{m\ge1}\frac{(tf_n)^m}{m!}$ has valuation $n+1$
(the $m=1$ term is lowest). By induction **$f_n$ has valuation exactly $n$**, so
$c^{(n)}_k=0$ for $k<n$. All coefficients are rational (the recurrence uses only
$\mathbb{Q}$-operations), and every $f_n$ is an **entire** function of $t$ (a finite
composition of entire functions), so each series has infinite radius of convergence.

Throughout, for a power series $u$ of positive valuation,
$e^{u}-1=\sum_{m\ge1}u^m/m!$, and $(u)^m$ has valuation $m\cdot\mathrm{val}(u)$.

---

## 1. Lemma A (leading coefficient) — *correcting a prior working note*

> **Lemma A.** $c^{(n)}_n=-1$ for every $n\ge1$.

*Proof.* Induction on $n$. Base: $f_1=e^{-t}-1=-t+\tfrac{t^2}{2}-\cdots$, so $c^{(1)}_1=-1$.
Step: assume $c^{(n)}_n=-1$. Then $tf_n$ has valuation $n+1$ with $[t^{n+1}](tf_n)=c^{(n)}_n$,
and for $m\ge2$, $(tf_n)^m$ has valuation $m(n+1)\ge 2(n+1)>n+1$. Hence
$$c^{(n+1)}_{n+1}=[t^{n+1}]\!\sum_{m\ge1}\tfrac{(tf_n)^m}{m!}=[t^{n+1}](tf_n)=c^{(n)}_n=-1.\qquad\square$$

**Correction of record.** An earlier working hypothesis (Track-F brief, hypothesis (a))
stated the leading coefficient is $(-1)^n$. That is **false**: it is the *constant* $-1$.
The claim is already refuted by the $n=2$ example printed in the official problem
statement, $\varphi_t^{2}(-1)=-t^{2}+\cdots$ (leading coefficient $-1$, not $(-1)^2=+1$).
The earlier warm-up scan tracked only $|c|$ (max $=1$ at the leading term), which is
consistent with either sign; the sign was never verified. Lemma A settles it.

Consequently $\bigl|c^{(n)}_n\bigr|=1$ for **all** $n$: the bound in Rippon's conjecture is
attained (with equality) on the entire leading diagonal, not at a single coefficient.

---

## 2. Lemma B (diagonal stabilization) — the main structural theorem

For an offset $j\ge0$ consider the diagonal sequence $n\mapsto c^{(n)}_{n+j}$.

> **Lemma B.** For every $j\ge0$ and every $n\ge \max(j,1)$,
> $$c^{(n+1)}_{(n+1)+j}=c^{(n)}_{n+j}.$$
> Hence $c^{(n)}_{n+j}$ is **constant for $n\ge\max(j,1)$**; call the value $A_j$. Then
> $A_0=-1$, and for $j\ge1$, $A_j=c^{(j)}_{2j}$.

*Proof.* By the recurrence,
$c^{(n+1)}_{(n+1)+j}=\sum_{m\ge1}\frac1{m!}[t^{\,n+1+j}](tf_n)^m.$
The $m$-th term is supported in degrees $\ge m(n+1)$, so it can contribute to degree
$n+1+j$ only if $m(n+1)\le n+1+j$, i.e. $(m-1)(n+1)\le j$. When $n\ge j$ we have
$n+1>j$, so for every $m\ge2$, $(m-1)(n+1)\ge n+1>j$ and the term does not contribute.
Only $m=1$ survives:
$$c^{(n+1)}_{(n+1)+j}=[t^{\,n+1+j}](tf_n)=[t^{\,n+j}]f_n=c^{(n)}_{n+j}.$$
Thus the sequence is constant from $n=\max(j,1)$ onward; its value is $A_j$, and for
$j\ge1$ the first stable index is $n=j$, giving $A_j=c^{(j)}_{2j}$. $\square$

**Threshold is tight.** For $j\ge2$ the value at $n=j-1$ generally differs from $A_j$
(computed: $j{=}2$ has $c^{(1)}_3=-\tfrac16\ne\tfrac13=A_2$; $j{=}3$ has
$c^{(2)}_5=-\tfrac{11}{24}\ne\tfrac{1}{24}=A_3$), so $\max(j,1)$ cannot be lowered.

**Equivalent "stable-profile" form.** Lemma B says exactly that
$$f_n(t)\equiv t^{\,n}\,\Phi(t)\pmod{t^{\,2n+1}},\qquad \Phi(t):=\sum_{j\ge0}A_j\,t^{j},$$
where $\Phi$ is a **single** formal power series independent of $n$: the coefficient array
$c^{(n)}_k$ is filled, on the whole triangle $n\le k\le 2n$, by the one sequence $A_j$ via
$c^{(n)}_k=A_{k-n}$. As $n\to\infty$, $t^{-n}f_n\to\Phi$ in the $(t)$-adic topology.

**Immediate corollaries (new closed-form coefficient values).**
- $A_0=-1$ ⟹ Lemma A.
- $A_1=c^{(1)}_2=[t^2](e^{-t}-1)=\tfrac12$ ⟹ $c^{(n)}_{n+1}=\tfrac12$ for **all** $n\ge1$.
- $A_2=c^{(2)}_4=\tfrac13$ ⟹ $c^{(n)}_{n+2}=\tfrac13$ for all $n\ge2$.
- In general $c^{(n)}_{n+j}=A_j$ for $n\ge j$ — infinitely many exact coefficients pinned
  by finitely many computations.

The first terms of the stable profile are
$$\Phi(t)=-1+\tfrac12 t+\tfrac13 t^2+\tfrac1{24}t^3+\tfrac7{60}t^4-\tfrac{59}{720}t^5
+\tfrac{53}{420}t^6-\tfrac{423}{4480}t^7+\cdots$$
(The magnitudes $1,\tfrac12,\tfrac13$ tempt one to guess $|A_j|=\tfrac1{j+1}$; this is
**false** — $A_3=\tfrac1{24}\ne\tfrac14$. The $A_j$ decay irregularly; see §4.)

---

## 3. Reduction of Rippon's conjecture

Partition the coefficient array by the line $k=2n$:

- **Stable region $n\le k\le 2n$ (offset $j=k-n\le n$).** Here $c^{(n)}_k=A_{k-n}$ by
  Lemma B. So on this region the conjecture is *equivalent to a one-variable statement:*
  $$\boxed{\ |A_j|\le 1\ \text{ for all } j\ge0.\ }\tag{I}$$
  Equality holds at $j=0$ ($A_0=-1$); data show $|A_j|\le\tfrac12$ for all $j\ge1$.
- **Transient region $k>2n$ (offset $j>n$).** The genuinely $n$-dependent part:
  $$\boxed{\ |c^{(n)}_k|\le 1\ \text{ for all } k>2n.\ }\tag{II}$$

**Rippon's conjecture $\iff$ (I) and (II).** Statement (I) is a bound on the single fixed
profile $\Phi$; (II) is where the interaction between iterates lives. Both hold throughout
the verified window (§4). This is the reduction we offer as the structural headline: the
"easy half" of the array collapses to one universal sequence, and the difficulty is
localized to $k>2n$.

---

## 4. Verified certificate and extremal-structure map (computation)

All data below are from exact rational arithmetic (python-flint `fmpq_series`, **no
rounding**), independently cross-checked by rigorous interval arithmetic (`arb_series`,
directed error bounds). Code and JSON receipts in the repo.

**Certificate (exact).** For all $1\le n\le 1000$ and all $n\le k\le 1000$:
$|c^{(n)}_k|\le1$, with equality **iff** $k=n$ (value $-1$). Zero violations.
This is a theorem about the exact coefficients (self-certifying rationals).

**Certificate (independent interval cross-check).** With `arb` at 256-bit precision the
same predicate is certified rigorously for $1\le n,k\le 700$ (worst ball radius
$<5\times10^{-50}$; $0$ disagreements with the exact values). At $N=800$, $256$ bits no
longer suffices (balls blow up); precision demand grows super-linearly in $N$, which is
exactly why exact rational arithmetic is the superior tool here (§5, Route B note).

**Extremal structure.**
- Global $\max|c^{(n)}_k|=1$, attained on the entire leading diagonal $k=n$ (Lemma A).
- **Non-leading supremum** over the window: $\ \displaystyle\sup_{k>n}|c^{(n)}_k|
  =\frac{2663}{4480}=0.594419\ldots$, attained **uniquely** at $(n,k)=(6,13)$.
  This point is $(n,k)=(j-1,\,2j-1)$ for $j=7$: one step *below* the offset-$7$
  stabilization threshold. The largest non-leading coefficients sit just under the
  stabilization line $k=2n$.
- **Transient ridge.** For each $n$, $\max_{k>n}|c^{(n)}_k|$ is attained either at
  offset $1$ (value $\tfrac12=A_1$) or near $k\approx 2n$; the near-$2n$ value
  **decreases** to $\tfrac12$ as $n\to\infty$ ($n{=}6\!:0.594,\ n{=}100\!:0.506,\
  n\!\ge\!150:0.500$).
- **Offset envelope** $E_j:=\sup_n|c^{(n)}_{n+j}|$ does *not* decay in $j$; it stays in a
  band $\approx[0.37,0.60]$ for $j\ge3$, always attained in the transient region $n<j$.

Everything is consistent with the sharpened conjecture:
$$|c^{(n)}_k|=1\iff k=n;\qquad \sup_{k>n}|c^{(n)}_k|=\tfrac{2663}{4480}\ \ (\text{at }(6,13)).$$

---

## 5. Three proof routes (attempted; status and obstructions)

### Route A — coefficientwise majorant / induction. *Provably cannot work as stated.*
One would seek a nonnegative series $M(t)=\sum M_k t^k$ with $M_k\le1$, $|c^{(1)}_k|\le M_k$,
and *stability*: $|f_n|\preceq M \Rightarrow |f_{n+1}|\preceq M$ (coefficientwise). The
majorization step **is** valid — if $|f_n|\preceq M$ then, since $\exp$ has nonnegative
coefficients, $\ |e^{tf_n}-1|\preceq e^{tM}-1$. But a *stable* majorant needs
$e^{tM}-1\preceq M$, and any such $M$ with $M_0=0$ forces $M\equiv0$ (the unique fixed
point of $M=e^{tM}-1$ in $t\mathbb{Q}[[t]]$ is $0$; sub-solutions only grow). Concretely,
the all-ones target $M=t/(1-t)$ gives $[t^k](e^{tM}-1)>1$ for large $k$. **Root cause:** a
coefficientwise majorant must dominate the *all-positive* series $\sum|c^{(n)}_k|t^k$,
whose partial sums exceed $1$; the bound $|c|\le1$ survives only through massive sign
cancellation, which any nonnegative majorant discards. A working variant must track the
**sign pattern** (mapped in §6), not just magnitudes.

### Route B — analytic (Cauchy) estimate. *Fails on any fixed circle.*
$c^{(n)}_k=\frac1{2\pi i}\oint_{|t|=r}\frac{f_n(t)}{t^{k+1}}\,dt$ gives
$|c^{(n)}_k|\le M_n(r)/r^{k}$ with $M_n(r)=\max_{|t|=r}|f_n(t)|$. Since every $f_n$ is
entire, there is no radius-of-convergence obstruction; but we computed
$M_n(1)=\max_{|t|=1}|f_n(t)|$ and it is **unbounded**, growing roughly linearly:
$M_1(1)=1.72,\ M_5=2.90,\ M_{10}=3.89,\ M_{25}=5.00,\dots$ (max attained off the real axis,
not at $t=-1$). So the fixed-circle bound at $r=1$ gives only $|c^{(n)}_k|\le M_n(1)\to\infty$.
A uniform $\le1$ bound would need an $(n,k)$-dependent contour/radius $r^*(n,k)$ making
$M_n(r^*)\le (r^*)^{k}$; for the leading term this must be essentially tight
($|c^{(n)}_n|=1$), and no clean choice of $r^*$ emerged. **Obstruction:** unbounded growth
of $M_n$ on fixed circles; the entirety of $f_n$ removes the convergence issue but not the
growth issue. (This growth is also why fixed-precision interval arithmetic loses
$\gtrsim\!N$ digits by degree $N$; exact arithmetic is unaffected.)

### Route C — combinatorial sign-reversing involution. *Most promising; not completed.*
The observed structure — $|c|=1$ *only* at the leading term, everything else strictly
smaller with heavy cancellation — is the signature of a sign-reversing involution with a
**single fixed point** (contributing the leading $-1$) on the combinatorial objects
enumerated by $c^{(n)}_k$. Since $e^{tz}-1$ is the EGF of nonempty sets, its $n$-fold
iteration evaluated at $-1$ expands $c^{(n)}_k$ as a signed weighted count of **length-$n$
chains of set partitions** (equivalently height-$\le n$ recursive/increasing set-partition
forests) on $k$ labels, the $-1$ seeding an overall sign. A weight-preserving,
sign-reversing involution on these chains with exactly one fixed point at $k=n$ would give
$|c^{(n)}_k|\le1$ with equality iff $k=n$ — matching the data exactly. The sign of $c^{(n)}_k$
is *not* a simple function of $k$ (sign table in §6), so the involution must be genuinely
structural. Constructing it is the natural next step and, if found, resolves the conjecture.

---

## 6. Sign structure (data for Route C)

Sign of $c^{(n)}_{n+j}$ (rows $n$, columns offset $j=0,\dots$):

```
 n=1: - + - + - + - + - +     (f_1: sign (-1)^k, the only clean row)
 n=2: - + + - + + - + + -
 n=3: - + + + - - + - - -
 n=4: - + + + + - + + + +
 n≥5: - + + + + - + - + + ...   (stabilizes column-by-column with A_j)
```
Leading column ($j=0$) is uniformly $-$; offset $1$ uniformly $+$. Signs of the stable
profile $A_j$ ($j=0..20$): $-,+,+,+,+,-,+,-,+,+,-,-,+,-,-,+,+,+,-,\dots$ (irregular).

---

## 7. Summary of what is proven vs. conjectured

| Statement | Status |
|---|---|
| $f_n$ has valuation $n$; entire in $t$ | **proved** |
| $c^{(n)}_n=-1$ for all $n$ (Lemma A) | **proved** |
| Diagonal stabilization $c^{(n)}_{n+j}=A_j$ for $n\ge\max(j,1)$ (Lemma B) | **proved** |
| $c^{(n)}_{n+1}=\tfrac12$, $c^{(n)}_{n+2}=\tfrac13$ (all large $n$) | **proved** (cor. of B) |
| Conjecture $\iff$ (I) $|A_j|\le1\,\forall j$ **and** (II) $|c^{(n)}_k|\le1\,\forall k>2n$ | **proved** (reduction) |
| $|c^{(n)}_k|\le1$, equality iff $k=n$, for $n,k\le1000$ | **verified** (exact) |
| same, $n,k\le700$ | **verified** (independent interval) |
| $\sup_{k>n}|c^{(n)}_k|=\tfrac{2663}{4480}$ at $(6,13)$ | **conjectured** (window-verified) |
| full conjecture $|c^{(n)}_k|\le1\ \forall n,k$ | **open** (routes A,B fail; C open) |
```
