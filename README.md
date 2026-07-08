# rippon-iterates — Rippon's iterated exponential (Hayman–Lingham Problem 7.54)

First recorded structural progress on **Problem 7.54** (P. J. Rippon), Hayman & Lingham,
*Research Problems in Function Theory* (50th anniversary ed., Springer 2019; arXiv:1809.07200).
Book update: *"No progress on this problem has been reported to us."*

**Problem.** $\varphi_t(z)=e^{tz}-1$; $\varphi_t^{n+1}=\varphi_t\circ\varphi_t^{n}$;
$f_n(t):=\varphi_t^n(-1)=\sum_k c^{(n)}_k t^k$ (formal power series in $t$, valuation $n$,
via $f_{n+1}=e^{t f_n}-1$). **Is $|c^{(n)}_k|\le 1$ for all $n\ge1,\,k\ge0$?**

## Results (see `proofs/proofs.md` for the mathematics)

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
