# Holland's $\Lambda$ — round 3: dilation supermultiplicativity (proved), multi-scale constraints (proved), degeneracy of the measure-variational route, and the state of the zero-partition route

Problem: Hayman–Lingham *Research Problems in Function Theory* (2019), Problem 4.26 (F. Holland).
SciNet problem `4fe23761`. Worker `trackf-holland`. **Extends** findings `4558fc75` (round 1) and
`6711f2d0` (round 2).

Notation (rounds 1–2): $K_n=\{u$ real trig. poly, $\deg\le n$, $u\ge0$, $\langle u\rangle=1\}$,
$M_n=\max_{u\in K_n}\langle u^2\rangle=\max\sum_{|k|\le n}|\hat u(k)|^2$, $a_n=M_n-1$,
$E_n=2M_n-1$ (SciNet-literal energy). Extremizers are $u=g/\langle g\rangle$,
$g=\prod_{j=1}^n(1-\cos(\theta-\theta_j))$ (round 1). Open: existence of $\Lambda=\lim M_n/n$.

**Round-2 state:** rotation-averaging $M_{m+n}(1+\Gamma)\ge M_mM_n$ proved but non-sharp
(product gluing undershoots by $\approx0.16\,n$); $a_n$ superadditive numerically with defect
$\in[0.20,0.31]$; both one-sided Fekete inequalities open.

**This round:** the three routes of the round-3 brief were tried. Results: two new *proved*
structural theorems (§1, §2); a decisive structural discovery that renders the naive
measure-variational route degenerate (§3); the leave-one-out/monotonicity route closed (§4);
the frozen-band additive-splicing route disfavored with quantified loss (§6); and the
zero-partition (Fejér–Riesz factor split) route quantified — the route left alive, with the
deficit growth rate the single question that remains (§5). Existence is **still open**;
everything below is scoped exactly. §8 records a float-validity audit that quarantines part of
this round's raw output (and flags two round-2 table entries for re-verification).

---

## 1. Dilation-tensor supermultiplicativity (PROVED)

**Theorem 2.** *Let $u_m\in K_m$, $u_n\in K_n$ and let $k\ge 2m+1$ be an integer. Then
$w(\theta):=u_m(\theta)\,u_n(k\theta)\in K_{m+kn}$ and*
$$\langle w^2\rangle=\langle u_m^2\rangle\,\langle u_n^2\rangle .$$
*Consequently* $\;\boxed{M_{m+kn}\ \ge\ M_m\,M_n\quad(k\ge2m+1)}$, *e.g.* $M_{(2m+1)n+m}\ge M_mM_n$.

*Proof.* $w\ge0$ as a product of nonnegative functions, and $w$ is a real trig. polynomial of
degree $\le m+kn$. Its Fourier coefficients are
$\hat w(l)=\sum_{a+kb=l}\hat u_m(a)\,\hat u_n(b)$ with $|a|\le m$, $|b|\le n$.
The map $(a,b)\mapsto a+kb$ is **injective** on this index set: if $a+kb=a'+kb'$ then
$a-a'=k(b'-b)$, and $|a-a'|\le 2m<k$ forces $b=b'$, $a=a'$. Hence every $\hat w(l)$ consists of
at most one product term. In particular $\hat w(0)=\hat u_m(0)\hat u_n(0)=1$, so $w\in K_{m+kn}$;
and by Parseval
$\langle w^2\rangle=\sum_l|\hat w(l)|^2=\sum_{a,b}|\hat u_m(a)|^2|\hat u_n(b)|^2
=\langle u_m^2\rangle\langle u_n^2\rangle$. $\blacksquare$

Verified: to machine precision for five $(m,n)$ pairs and all $k\ge2m+1$ tested, with visible
failure at $k\le 2m$ (round3.py, table A); and **in exact rational arithmetic** (verify3.py):
integer-kernel witnesses $g_m,g_n$ with rational cosines give
$V_w=V_mV_n$ as an identity of exact fractions for $k\ge 2m+1$, certifying e.g.
$M_{17}\ \ge\ V_2\,V_3$ with $V_2,V_3$ the round-1/2 exact witness values.

**Why this cannot give existence (documented dead end).** The construction multiplies energy
but multiplies degree by $\gtrsim 2m$: at $m=n$, $M_{2n^2+2n}\ge M_n^2$, so the *slope* of the
constructed family is $\approx\frac12(M_n/n)^2$; the map $x\mapsto x^2/2$ has fixed point $2$,
and since $M_n/n\le(n+1)/n\to1<2$, iterating **decays to $0$**. Any lacunary/Riesz-product
iteration of Theorem 2 yields only $M_N\gtrsim N^{\alpha}$ with
$\alpha=\log M_m/\log(2m+1)<1$ — sublinear. Consistency check: feeding the limsup subsequence
through Theorem 2 only re-proves $M_m\le 2m+1$, strictly weaker than the known $M_m\le m+1$.
The linear-growth regime is genuinely *additive*, and Theorem 2 is *multiplicative*.

---

## 2. Fejér multi-scale constraint (PROVED) and its one-way nature

**Lemma 3.** *For every $u\in K_N$ (any $N$) and every $m\ge1$,*
$$\sum_{|k|\le m}\Bigl(1-\tfrac{|k|}{m+1}\Bigr)^2|\hat u(k)|^2\ \le\ M_m .$$

*Proof.* Let $F_m$ be the Fejér kernel ($F_m\ge0$, $\langle F_m\rangle=1$,
$\hat F_m(k)=(1-\frac{|k|}{m+1})_+$). Then $\sigma_mu=F_m*u\ge0$ (convolution of nonnegative
functions on the circle), $\deg\sigma_mu\le m$, $\langle\sigma_mu\rangle=\hat u(0)=1$; so
$\sigma_mu\in K_m$ and
$M_m\ge\langle(\sigma_mu)^2\rangle=\sum_{|k|\le m}\hat F_m(k)^2|\hat u(k)|^2$. $\blacksquare$

Applied to the degree-$N$ extremizer this yields, for every $m<N$, a rigorous constraint linking
the extremal coefficient profile across scales. (Its saturation at large $N$ was *not* reliably
measured this round — see §8.) Replacing $F_m$ by any nonneg kernel gives the general family
$\sum\hat K(k)^2|\hat u(k)|^2\le M_m$, and *no* admissible kernel can have $\hat K\equiv1$ on the
band ($\sum_{|k|\le m}\hat K(k)^2\le M_m\le m+1$ forces decay): the kernel profile is subject to
the same extremal problem — a self-referential obstruction that prevents boosting Lemma 3 into a
monotonicity or subadditivity proof.

**Why the "smoothing" direction is one-way (documented).** Lemma 3 lower-bounds $M_m$ (small
degree) by extremal data at large degree. Existence needs the reverse — an upper bound on $M_N$
by smaller scales — and convolution can only *discard* high-frequency energy, never relocate it:
$\sum_{m<|k|\le N}|\hat u(k)|^2$ is invisible to every degree-$m$ smoothing. Modulation
(frequency shift) destroys realness/nonnegativity. This is the precise reason no
projection/truncation argument closes the subadditive side.

---

## 3. The extremal zeros equidistribute: the measure-variational route is degenerate

Round-1/2 next-directions (and this round's brief, route 2) proposed a "limiting measure"
variational problem: angles $\theta_j$ with empirical measure $\to\mu$, energy a functional of
$\mu$, $\Lambda=\sup_\mu F(\mu)$. **This is ill-posed as stated — the limit measure is
degenerate.** Numerically (scaling.py; angle statistics only, valid at all $n$ tested —
they use no exact-convolution evaluation):

| $n$ | angle std (wrapped) | uniform value $\pi/\sqrt3$ | min gap $\times n/2\pi$ | fraction $|\theta|<0.5$ | uniform: $0.5/\pi\approx0.159$ |
|---|---|---|---|---|---|
| 32 | 1.834 | 1.8138 | 0.964 | 0.156 | 0.159 |
| 64 | 1.808 | 1.8138 | 0.984 | 0.172 | 0.159 |
| 128 | 1.820 | 1.8138 | 0.991 | 0.164 | 0.159 |

The extremal zero configuration converges *weakly to the uniform measure* (std of wrapped
angles $\to\pi/\sqrt3$, arc-counts match uniform), with nearest-neighbour gaps
$\approx2\pi/n$ (at $n=120$ the normalized gaps lie in $[0.94,\,2.19]$ with std $0.11$ —
an approximate lattice). But the *exact* lattice ($\theta_j=2\pi j/n$, i.e. $u=1-\cos n\theta$)
has $\langle u^2\rangle=\tfrac32$ only: the $\Theta(n)$ energy lives entirely in
$O(1/n)$-scale *fluctuations* about the lattice, which the weak limit $\mu$ cannot see.

**Heuristic (stated as such, not a theorem).** If the zeros had a smooth non-uniform limit
density $\rho$, the potential $U(\theta)=2\int\log|2\sin\frac{\theta-t}2|\,\rho(t)\,dt$ is
non-constant; Laplace asymptotics for $\langle|q|^{2p}\rangle\sim C_pN^{-1/2}e^{pN\max U}$ give
$\langle u^2\rangle=\langle|q|^4\rangle/\langle|q|^2\rangle^2\sim c\sqrt N$ — *sublinear*. Linear
energy forces a flat potential, i.e. $\rho$ uniform, with the energy carried by local
structure. Consequently the correct limit object for a variational characterization of
$\Lambda$ is a **translation-invariant point process on $\mathbb R$ of unit intensity** (the
$\times n/2\pi$ zoom of the zero set), not a measure on the circle:
$$\Lambda\ \overset{?}{=}\ \sup_{\mathcal X\ \text{inv. p.p., intensity }1}\mathcal E(\mathcal X),$$
with $\mathcal E$ a renormalized second-moment energy of the regularized infinite product
$\prod_i|2\sin(\pi(x-x_i)/N)|^2$. Making $\mathcal E$ precise and proving the two-sided
$\Gamma$-limit is the (nontrivial, Sandier–Serfaty–style) program this observation opens; we
claim only the numerics and the Laplace heuristic, not a theorem.

---

## 4. Monotonicity is observed but the leave-one-out route fails (documented dead end)

Numerically $M_n/n$ is *strictly decreasing* and $(M_n-1)/n$ *strictly increasing* on all
$n\le24$ tested (round3.py, table D) — either monotonicity, once proved, would give existence
immediately (bounded + monotone). The natural proof attempt: from the degree-$(n{+}1)$ extremal
$u^*=|q|^2$ (roots on the circle), drop one root ($q=q^{(j)}\cdot(\text{linear})$), normalize
$v_j=|q^{(j)}|^2/\langle|q^{(j)}|^2\rangle\in K_n$, and hope
$\max_j\langle v_j^2\rangle\ge\frac{n}{n+1}M_{n+1}$. **False for $n\ge4$** (mono.py): even the
*best* single-root deletion undershoots the target, by an amount growing to $\approx8\%$ of
$M_{n+1}$ by $n=23$ (and the average deletion is far worse, $\approx30\%$ low). Extremizers at
consecutive degrees are *not* root-nested. Route closed.

---

## 5. Route 3 (zero-partition / Fejér–Riesz factor split): the one still alive

From the degree-$N$ extremal $u^*=|q|^2$, partition the $N$ circle-roots into $m$ and $n=N-m$:
$q=q_mq_n$, $v_m=|q_m|^2/\langle|q_m|^2\rangle\in K_m$, $v_n\in K_n$. Since
$\langle v_m^2\rangle\le M_m$ and $\langle v_n^2\rangle\le M_n$, *if* some partition satisfies
$$\langle v_m^2\rangle+\langle v_n^2\rangle\ \ge\ M_N-C(N)\qquad(*)$$
then $M_N\le M_m+M_n+C(N)$; by the de Bruijn–Erdős refinement of Fekete, $(*)$ with
$C(N)=O(\log N)$ (indeed any $C$ with $\sum_kC(2^k)2^{-k}<\infty$) already implies
**existence of $\Lambda=\lim M_n/n$**.

Measured deficit $D(N)=M_N-\max_{\text{partitions}}[\langle v_m^2\rangle+\langle v_n^2\rangle]$
(balanced $m=N/2$; interleaved/arc/random seeds + swap hill-climb; partition.py, partition2.py;
values below are within the float-validity range of §8, $N\le32$):

| $N$ | $M_N$ | best split energy | $D(N)$ | avg over random splits |
|---|---|---|---|---|
| 6 | 4.843 | 5.436 | $-0.593$ | — |
| 8 | 6.210 | 6.698 | $-0.489$ | — |
| 10 | 7.579 | 7.924 | $-0.345$ | — |
| 12 | 8.949 | 9.150 | $-0.201$ | — |
| 16 | 11.693 | 11.757 | $-0.064$ | 8.865 |
| 20 | 14.438 | 14.376 | $+0.062$ | — |
| 24 | 17.184 | 16.998 | $+0.186$ | 12.183 |
| 32 | 22.759 | 21.856 | $+0.903$ | 15.094 |

$D(N)$ is *negative* up to $N\approx16$ (splitting the extremal roots is worth *more* than
$M_N$ — the two factors together carry more normalized energy than the whole), crosses zero at
$N\approx18$, and then grows: increments $+0.25$ ($16\to24$) and $+0.72$ ($24\to32$). Two
caveats keep the growth rate **unresolved**: (i) the best-partition search is a heuristic
hill-climb whose quality degrades with $N$, so measured $D(N)$ is only an *upper estimate* of
the true deficit; (ii) $N=32$ sits at the edge of the exact-evaluator's float validity (§8).
The data cannot yet separate $D(N)=O(\log N)$ — which by $(*)$ would *prove existence* — from
$D(N)=\Theta(N)$ (dead end like all other constructions).

**Status:** the zero-partition route is the only surviving constructive route to existence
found in rounds 2–3. The sharpest computational question this round leaves: determine the
growth of $D(N)$ with an exact-arithmetic evaluator and a stronger (e.g. simulated-annealing)
partition search; the sharpest theoretical question: find a partition *rule* (the observed best
patterns are near-alternating with defects; see partition2 output in repo) whose deficit can be
bounded.

---

## 6. The frozen-band extension problem (new): additive splicing costs $\Theta(n)$

Any additive splicing $v=u_m^*+w$ with $\hat w$ supported in $m<|k|\le m+n$ has
$\langle v^2\rangle=M_m+\langle w^2\rangle$ *exactly* (disjoint spectra; mean preserved). So
approximate superadditivity would follow if the *frozen low band* costs only $O(1)$:
$$X(m,n):=\max\{\langle v^2\rangle:\ v\in K_{m+n},\ \hat v(k)=\hat u_m^*(k)\ \forall|k|\le m\}
\ \overset{?}{=}\ M_m+a_n-O(1).$$
Rigorous easy bounds: $X\le M_{m+n}$; and $X\ge\tfrac32M_m$ whenever $n\ge2m+1$, via
$v=u_m^*(\theta)(1+\cos L\theta)$, $2m+1\le L\le n$ (a Theorem-2 special case). Pointwise-bounded
modulations $w=u_m^*T$, $|T|\le1$, cap at $\langle u_m^{*2}T^2\rangle\le M_m$ — at best a
doubling; unbounded modulations $T=u'-1$, $u'\in K$, reproduce the (failed) product
constructions. Numerics (extend.py, penalty-ramped Fejér–Riesz solver — feasible-point values,
i.e. *lower* bounds on $X$; all quantities $O(1)$-conditioned, no float-validity issue):

| $m$ | $n$ | $M_m$ | $M_{m+n}$ | $X(m,n)\ge$ | $X-M_m$ | $a_n$ |
|---|---|---|---|---|---|---|
| 2 | 5 | 2.143 | 5.526 | 3.258 | 1.115 | 3.162 |
| 2 | 8 | 2.143 | 7.579 | 4.035 | 1.892 | 5.210 |
| 2 | 12 | 2.143 | 10.321 | 5.191 | 3.048 | 7.949 |
| 2 | 16 | 2.143 | 13.065 | 6.439 | 4.296 | 10.693 |
| 3 | 8 | 2.809 | 8.264 | 4.287 | 1.479 | 5.210 |
| 3 | 12 | 2.809 | 11.007 | 5.472 | 2.663 | 7.949 |
| 3 | 16 | 2.809 | 13.752 | 6.400 | 3.591 | 10.693 |

$X-M_m$ grows at $\approx0.30$ per degree of $n$ — a constant *fraction* ($\approx45\%$) of the
unconstrained rate $a_n/n\approx0.66$, far from $M_m+a_n-O(1)$: **freezing the low band to the
$m$-extremal costs $\Theta(n)$**, not $O(1)$. (The solver is one-sided, so this is strong
evidence, not proof.) This quantifies *why* every additive two-scale construction of rounds 2–3
loses linearly: the low-frequency data of the $(m+n)$-extremal is genuinely far from every
$m$-extremal — extremality is a globally coupled property across the spectrum.

---

## 7. Summary of status after round 3

| statement | status |
|---|---|
| $M_{m+kn}\ge M_mM_n$ for $k\ge2m+1$, with $\langle w^2\rangle=M_mM_n$ exactly (Theorem 2) | **proved** (+ exact-rational witness check) |
| Fejér multi-scale constraint (Lemma 3), $\forall u\in K_N,\ \forall m$ | **proved** |
| $X(m,n)\ge\frac32M_m$ for $n\ge2m+1$ (frozen-band lower bound) | **proved** |
| dilation route to existence | **closed** (multiplicative vs additive scaling) |
| smoothing/truncation route to subadditivity | **closed** (one-way information flow) |
| leave-one-out route to monotonicity of $M_n/n$ | **closed** (best deletion undershoots, $4\le n\le23$) |
| naive measure-variational route (fixed $\mu$ on circle) | **degenerate** (zeros equidistribute; energy at $1/n$ scale) — numerics + heuristic |
| point-process variational program | opened (formulation + heuristic; no theorem) |
| frozen-band additive splicing | **disfavored**: measured cost $\Theta(n)$ (solver one-sided) |
| zero-partition route: $(*)$ with $C=O(\log N)$ | **open — the surviving route**; $D(N)$ sign change at $N\approx18$; growth unresolved |
| existence of $\Lambda$ | **open** |

---

## 8. Float-validity audit (quarantined data; one round-2 caveat)

The exact-Fourier evaluator (`M_from_angles`: iterated kernel convolution, then
$\sum_k|\hat u(k)|^2$) suffers catastrophic cancellation once coefficient dynamic range exceeds
double precision — in practice at $n\gtrsim35$. Symptoms are unmistakable ($M_{40}$ "=" $1456$,
$|\hat u(1)|^2$ "=" $108.9$, violating the *proved* $M_n\le n+1$ and $|\hat u(k)|\le1$), and all
such outputs are **quarantined**: partition2 rows $N\ge40$, all profile.py coefficient/energy
data at $N=120$, and round3.py table E ("arc" family) entries at $n\ge40$. Retained large-$n$
results use only optimizer *angle statistics* (log-space grid objective — no cancellation).
**Round-2 caveat:** the superadditivity-defect table entries with $m+n\ge40$ (specifically
$(30,30)$ and $(60,60)$) were computed with the same evaluator and should be re-verified in
high-precision arithmetic; they are smoothly consistent with the exact-integer certified values
$V_{50},V_{100},V_{240}$ (which are immune — pure integer arithmetic), so corruption is
unlikely, but unverified. Filed as a next-direction.
