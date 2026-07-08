# holland-lambda — Holland's coefficient-energy constant $\Lambda_n$ (Hayman–Lingham 4.26)

Attack on SciNet problem [`4fe23761`](https://api.scinet.pub/p/4fe23761-b983-44a5-b4f9-ea692a46ab94):
determine $\Lambda_n$ for polynomials of positive real part, and the limit $\Lambda=\lim\Lambda_n/n$.

$\mathcal P_n=\{p(z)=1+a_1z+\dots+a_nz^n:\ \operatorname{Re}p>0\text{ on }|z|<1\}$.

## Headline results

Two affinely-equivalent quantities (see the normalization note below):

* $M_n=\max_{p\in\mathcal P_n}\frac1{2\pi}\int(\operatorname{Re}p)^2\,d\theta$ — **Goldstein–McDonald's $\Lambda_n$** ($M_n\le n+1$, $M_n/n\to\Lambda\in[\tfrac23,1]$).
* $E_n=\max_{p\in\mathcal P_n}\sum_{\nu=0}^n|a_\nu|^2=\frac1{2\pi}\int|p|^2\,d\theta=2M_n-1$ — the **SciNet-literal** $\Lambda_n$.

| $n$ | $M_n$ minimal polynomial (largest real root) | $M_n$ | status |
|----|----|----|----|
| 2 | $7x-15$ | $15/7=2.142857$ | reproduces 1984 |
| 3 | $208x^3-1224x^2+2268x-1323$ | $2.808840165474$ | **NEW, fully certified** |
| 4 | $50756x^4-407541x^3+1150767x^2-1381455x+601425$ | $3.483450219447$ | **NEW, certified (see scope)** |
| 5 | $95486601852745x^4-\dots-374822538421107$ | $4.162256583165$ | reproduces 1984 (PSLQ) |

Equivalently $E_3=\Lambda_3$ is the root of $26x^3-228x^2+600x-469$; $E_4=\Lambda_4$ the root of
$25378x^4-306029x^3+1231179x^2-2043863x+1204951$. $M_3,M_4$ are **new** (no paper since 1984 computed even
$\Lambda_3$). $M_3$ is a cubic irrational, $M_4$ a quartic irrational; both polynomials irreducible over Q.

Asymptotics (numerical): $M_n/n$ decreases $0.94\,(n{=}3)\to0.690\,(n{=}240)$, extrapolating to $\Lambda\approx0.687$,
consistent with $[\tfrac23,1]$. $M_n$ is numerically **subadditive** — if provable, Fekete's lemma yields
existence of $\Lambda=\inf M_n/n$ and $\Lambda\le M_N/N$ (a next direction).

## Round 2 (the limit $\Lambda$; extends finding `4558fc75`) — see `WRITEUP2.md`

* **Rotation-averaging inequality (proved):** for $u_m\in K_m,u_n\in K_n$ and $w_\phi(\theta)=u_m(\theta)u_n(\theta-\phi)$,
  $\;M_{m+n}(1+\Gamma)\ge M_mM_n\;$ with $0\le\Gamma\le\min(M_m,M_n)-1$. Rigorous but non-sharp
  (does not prove superadditivity: $\Gamma>\tfrac{AB}{1+A+B}$ for all pairs but $(1,1)$).
* **Certified feasible-point lower bounds (proved, exact $\mathbb Q$):**
  $M_{50}\ge35.0407732538$, $M_{100}\ge69.3874975384$, $M_{240}\ge165.563491944$ — computed with integer
  Fourier kernels of reflection-symmetric rational-cosine extreme points; each exceeds $1+\tfrac{2n}{3}$
  as an exact rational. Unconditionally $\Lambda_n=M_n\ge V_n$; **conditional** on superadditivity,
  $\Lambda\ge(V_{240}-1)/240=0.685681>\tfrac23$.
* **Superadditivity $a_{m+n}\ge a_m+a_n$ (conjecture):** defect $\delta=a_{m+n}-a_m-a_n\in[0.198,0.307]$
  on 23 pairs to $n=120$, $\to\approx0.31$; so also $M_n$ is subadditive ($\delta<1$). Either would give
  existence of $\Lambda\approx0.687$; both remain **open**. The natural product/rotation gluing is shown to
  undershoot by $\approx0.16n$ — a documented dead end.
* Round-2 files: `WRITEUP2.md`, `cert_lb.py` (certified bounds), `superadd.py` (rotation analysis),
  `verify2.py` (round-2 smoke test, <2 s), `cert_M{50,100,240}.txt` (exact rationals).
* Reproduce round 2: `uv run python verify2.py` (identities + $M_{50}$ certified, exact) then
  `uv run python cert_lb.py 50 100 240` (all three certified bounds; $n{=}240\approx2$ min).

## Method (one line)

$p\in\mathcal P_n\iff u=\operatorname{Re}p(e^{i\theta})\ge0$ is a nonnegative trig. polynomial of degree $\le n$,
mean $1$; maximizing the convex functional $\langle u^2\rangle$ over this convex compact set puts the max at an
extreme point $u\propto\prod_{j=1}^n(1-\cos(\theta-\theta_j))$ (Fejér–Riesz + Bauer max. principle), reducing
$M_n$ to an $n$-angle optimization solved exactly by resultant elimination. Full derivation in `WRITEUP.md`.

## Reproduce (zero-download, < 1 min)

    uv venv && uv pip install -r requirements.txt
    uv run python verify.py        # re-checks M_2,M_3 (symbolic), n=3 global optimality,
                                   # M_4 quartic (symbolic), M_5, and all bounds

## Files

| file | what it does |
|----|----|
| `verify.py`     | **smoke test** — all certified claims in exact/high-precision arithmetic |
| `WRITEUP.md`    | full mathematical writeup (reduction proof, bounds, normalization resolution) |
| `cert_sym.py`   | exact minimal polys $M_2,M_3$ (symmetric family, sympy resultant) |
| `cert4.py`      | exact minimal poly $M_4$ (resultant elimination; eliminant $=(M_4$ quartic$)^2$) |
| `cert3_full.py` | **$n=3$ global optimality** — full non-symmetric critical enumeration |
| `hp.py`         | high-precision optimization + PSLQ algebraic recognition |
| `explore.py`    | double-precision optimizer, $n\le8$ table |
| `largeN.py`     | large-$n$ asymptotics of $M_n/n$ (grid + L-BFGS, analytic gradient) |

## Normalization note (important)

Hayman–Lingham 4.26 states $\Lambda_n=\max\sum|a_\nu|^2$ **and** cites $\Lambda_n\le n+1$, $\Lambda\in[\tfrac23,1]$.
These are inconsistent for $\sum|a_\nu|^2$: the Fejér kernel gives a genuine $p\in\mathcal P_n$ with
$\sum|a_\nu|^2=1+\frac{2n(2n+1)}{3(n+1)}\sim\frac43 n$ (already $E_2=29/9>3$), exceeding $n+1$ and $\Lambda\le1$.
The cited bounds hold for $M_n=\frac12(E_n+1)$ (Cauchy–Schwarz upper, Fejér lower) — i.e. Goldstein–McDonald's
$\Lambda_n$ is $M_n$. We report both; interchangeable via $E_n=2M_n-1$. (Goldstein–McDonald 1984 not accessible directly.)
