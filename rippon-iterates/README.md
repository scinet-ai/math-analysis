# rippon-iterates — Rippon's iterated exponential (Hayman–Lingham Problem 7.54)

First recorded structural progress on **Problem 7.54** (P. J. Rippon), Hayman & Lingham,
*Research Problems in Function Theory* (50th anniversary ed., Springer 2019; arXiv:1809.07200).
Book update: *"No progress on this problem has been reported to us."*

**Problem.** $\varphi_t(z)=e^{tz}-1$; $\varphi_t^{n+1}=\varphi_t\circ\varphi_t^{n}$;
$f_n(t):=\varphi_t^n(-1)=\sum_k c^{(n)}_k t^k$ (formal power series in $t$, valuation $n$,
via $f_{n+1}=e^{t f_n}-1$). **Is $|c^{(n)}_k|\le 1$ for all $n\ge1,\,k\ge0$?**

## ROUND 2 results (see `proofs/proofs2.md`; extends finding cda0ff02)

Round 2 replaces the round-1 reduction with a complete **band decomposition** of the
coefficient array via a *formal Koenigs linearization at multiplier $t$*:

- **Theorem 1 (exact factorization).** $f_n=-t^n\prod_{k<n}h(tf_k)$ with $h(x)=(e^x-1)/x$;
  hence $\Phi=-\prod_{k\ge0}h(tf_k)$ ($t$-adic) — new proof of stabilization. *Proved*;
  product verified exactly to 120 profile terms (`rippon/product_form.py`).
- **Theorem 2 (formal Koenigs).** Unique $P(w)=w+\sum_{m\ge2}p_m(t)w^m$, $p_m\in\mathbb Q[[t]]$,
  with $P(tw)=e^{tP(w)}-1$; $\operatorname{val}p_m\ge m-1$; $p_2=\frac{t}{2(t-1)}$,
  $p_3=\frac{t^2(t+2)}{6(t-1)^2(t+1)}$ (poles only at roots of unity). *Proved.*
- **Theorem 3 (band decomposition).** $f_n=\sum_m p_m t^{mn}\Phi^m$, so every coefficient
  is a **finite** sum $c^{(n)}_k=\sum_{m\le(k+1)/(n+1)}[t^{k-mn}]\Psi_m$ of coefficients of
  the universal band profiles $\Psi_m=p_m\Phi^m$. **Rippon's conjecture is equivalent to a
  family of inequalities about these fixed universal series** — the iterate index $n$ only
  selects which coefficients are summed. *Proved.*
- **Theorem 4 (transient lines).** $c^{(n)}_{2n+i}=A_{n+i}+\tau_i$ ($i\le n+1$) and
  $c^{(n)}_{3n+i}=A_{2n+i}+\tau_{n+i}+\sigma_i$ ($i\le n+2$), $\tau=\Psi_2$-, $\sigma=\Psi_3$-
  coefficients, $\tau_i=-\frac12\sum_{l<i}[t^l]\Phi^2$. Corollaries: $c^{(n)}_{2n+1}=A_{n+1}-\frac12$,
  $c^{(n)}_{2n+2}=A_{n+2}$, ridge $\to\frac12$, and the round-1 extremal
  $\frac{2663}{4480}$ at $(6,13)$ **explained** ($=|A_7-\frac12|$). *Proved* + exact check
  (3960 cells, 0 failures, `rippon/bands.py`).
- **Theorem 5 (analytic profile).** (a) radius of convergence $R\ge\ln2$; (b)
  $-1<\Phi(t)<0$ on the real segment $(0,1)$ — the profile bound holds there; (c)
  **certified** (arb balls + proved tail bound, `rippon/phi_cert.py`):
  $\Phi(-\frac12)<-1.159<-1$, so $\max_{|t|=r}|\Phi|\ge1.159$ for all $r\in[\frac12,R)$ —
  **sup-norm/Cauchy routes provably cannot prove $|A_j|\le1$**, even in the $r\uparrow R$ limit.
- **Profile data (exact, 600 terms).** $|A_j|\le1$ for $j\le600$; $\max_{j\ge1}|A_j|=\frac12$
  only at $j=1$; envelope decay $\sim j^{-0.88}$, $|A_j|^{1/j}\to0.99$ ($R=1$ conjectured);
  $\Phi$ empirically **not P-recursive** (no recurrence, order $\le8$); radial limits at
  root-of-unity directions $\to0$ (parabolic boundary zeros; natural-boundary picture);
  $\Phi(t)\sim2(t-1)$ at $t=1^-$; $|\tau_i|\le\frac12$, $|\sigma_i|\le\frac13$ ($i\le600$);
  $L^1$ circle means of $|\Phi|$ are $>1$ (H^1 route obstructed too, by only ~13%).
- **Lean 4 (round 2).** `lean/Rippon2.lean` machine-checks (`native_decide`): factorization
  window ($n\le30$, deg $\le60$), band-2 identity ($n\le20$), corollaries ($n\le25$).

## Round-1 results (see `proofs/proofs.md`)

- **Lemma A (leading coefficient).** $c^{(n)}_n=-1$ for all $n$ — the *constant* $-1$.
  This **corrects** a prior working note that guessed $(-1)^n$ (already refuted by the
  official statement's $\varphi_t^2(-1)=-t^2+\cdots$). *Proved.*
- **Lemma B (diagonal stabilization).** For each offset $j$, $c^{(n)}_{n+j}=A_j$ is
  constant for all $n\ge\max(j,1)$; equivalently $f_n\equiv t^n\Phi(t)\pmod{t^{2n+1}}$ for a
  single fixed profile $\Phi=\sum_j A_j t^j$. Hence $c^{(n)}_{n+1}=\tfrac12$ and
  $c^{(n)}_{n+2}=\tfrac13$ for all large $n$ — infinitely many exact coefficients pinned.
  *Proved.*
- **Reduction.** Rippon's conjecture $\iff$ (I) $|A_j|\le1\ \forall j$ (a one-variable
  statement about $\Phi$, covering the whole region $k\le2n$) **and** (II)
  $|c^{(n)}_k|\le1\ \forall k>2n$. *Proved.*
- **Verified certificate (exact).** $|c^{(n)}_k|\le1$ with equality **iff** $k=n$ (value
  $-1$), for **all $1\le n,k\le 1000$**. Exact rational arithmetic (no rounding).
- **Independent interval cross-check.** Same predicate certified rigorously by `arb`
  ball arithmetic (256-bit, directed rounding) for $n,k\le 700$; $0$ disagreements with
  the exact values.
- **Extremal structure.** Global $\max|c|=1$ on the entire leading diagonal; the
  **non-leading supremum is $\tfrac{2663}{4480}=0.59442\ldots$, attained uniquely at
  $(n,k)=(6,13)$** (one step below the offset-7 stabilization threshold). The transient
  ridge near $k\approx 2n$ decreases to $\tfrac12$ as $n\to\infty$.
- **Proof routes.** Coefficientwise majorant — *provably cannot work* (kills the sign
  cancellation the bound relies on). Analytic Cauchy on $|t|=1$ — *fails* ($\max_{|t|=1}|f_n|$
  grows unboundedly, $\to 5$ by $n{=}25$). Combinatorial sign-reversing involution — *open,
  most promising* (matches "equality only at leading term").
- **Lean 4.** `lean/Rippon.lean` machine-checks (`native_decide`) the no-violation +
  leading$=-1$ + equality-only-at-leading claims for the window $n,k\le 40$. Trusted base:
  Lean kernel + native compiler (**`native_decide` — disclosed, not pure-kernel**).

## Layout
```
rippon/engine.py       exact (fmpq_series) iterate engine
rippon/certificate.py  exact rational certificate + A_j + envelope   (headline)
rippon/arb_scan.py     independent rigorous arb-ball cross-check
rippon/structure.py    full extremal-structure map (stabilization audit, A_j, E_j)
rippon/extremal.py     non-leading supremum trend tracker
proofs/proofs.md       theorems, reduction, proof-route analysis, sign data
lean/Rippon.lean       native_decide finite certificate (window 40)
data/*.json            saved certificate / structure receipts
reproduce.sh           zero-download smoke (exact + arb + structure)
```

## Reproduce (needs [`uv`](https://docs.astral.sh/uv/))
```bash
./reproduce.sh          # smoke: exact + arb cross-check to N=250   (seconds)
./reproduce.sh full     # headline: exact N=1000 + arb N=700         (~30 s)
```
Pinned env in `requirements.txt` (python-flint 0.9.0, numpy 2.4.6).

Direct invocations:
```bash
uv run --no-project --with python-flint python3 -m rippon.certificate 1000 out.json
uv run --no-project --with python-flint python3 -m rippon.arb_scan 700 256 700
uv run --no-project --with python-flint python3 -m rippon.structure 250
# Lean:
cd lean && elan run leanprover/lean4:v4.32.0-rc1 lean Rippon.lean
```

## Honesty notes
- The full conjecture is **open**. What is *proved* are Lemmas A, B and the reduction;
  the $\le1$ bound itself is *verified* over a finite window, not proved for all $(n,k)$.
- "Verified" = exact rational or rigorous interval arithmetic over a stated window;
  floating point is used only for the (clearly-labeled) analytic circle-map heuristic.
- The Lean certificate uses `native_decide` (compiler-trust), disclosed in the finding's
  `formal` block; it is not a pure-kernel proof.
