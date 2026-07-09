"""Band-decomposition verification for Rippon 7.54 (round 2).

THEOREM (band decomposition, proved in proofs/proofs2.md): with Phi the stable
profile and p_m(t) the formal Koenigs/Poincare coefficients (p_1=1,
p_2 = t/(2(t-1)), p_3 = t^2(t+2)/(6(t-1)^2(t+1)), ...; val p_m >= m-1),

    f_n(t) = sum_{m>=1} p_m(t) t^{mn} Phi(t)^m        (t-adically convergent),

so with Psi_m := p_m Phi^m (universal series, independent of n):

    c^{(n)}_k = sum_{m>=1, k-mn >= val Psi_m} [t^{k-mn}] Psi_m.

Specializations verified here in EXACT rational arithmetic:
  band<=1 (k <= 2n):        c^{(n)}_k = A_{k-n}                       [Lemma B]
  band<=2 (2n < k <= 3n+1): c^{(n)}_{2n+i} = A_{n+i} + tau_i,  tau_i=[t^i]Psi_2
  band<=3 (3n+1 < k <= 4n+2): c^{(n)}_{3n+i} = A_{2n+i} + tau_{n+i} + sigma_i,
                              sigma_i = [t^i]Psi_3   (i <= n+2)
Also checks THRESHOLD TIGHTNESS: the band-2 formula must FAIL at i = n+2 by
exactly sigma_{...} etc., i.e. the first omitted band term is the exact defect.

Usage: python -m rippon.bands <NMAX> [out.json]
"""
from __future__ import annotations
import sys, json, time
from fractions import Fraction
from flint import fmpq_series, fmpq, ctx
from rippon.engine import iterates_exact
from rippon.profile import to_frac


def phi_series(J: int) -> fmpq_series:
    """Stable profile Phi as fmpq_series to order J, via Phi = t^{-J} f_J mod t^{J+1}.
    (c^{(J)}_{J+j} = A_j for j <= J by Lemma B.)"""
    ctx.cap = 2 * J + 1
    t = fmpq_series([0, 1])
    f = (fmpq_series([0, -1])).exp() - 1
    for n in range(1, J):
        f = (t * f).exp() - 1
    # f = f_J to order 2J; extract A_j = [t^{J+j}] f_J
    A = [f[J + j] for j in range(0, J + 1)]
    ctx.cap = J + 1
    return fmpq_series(A), A


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    t0 = time.time()
    # coefficient array to degree N = 4*NMAX + 6 to cover band-3 checks
    N = 4 * NMAX + 6
    C = {}
    for n, f in iterates_exact(N, NMAX):
        C[n] = [to_frac(f[k]) for k in range(N + 1)]

    # profile to order J = 3*NMAX + 8 (need A up to 2n+i <= 3NMAX+8)
    J = 3 * NMAX + 8
    Phi, Araw = phi_series(J)
    A = [to_frac(a) for a in Araw]

    ctx.cap = J + 1
    tser = fmpq_series([0, 1])
    one = fmpq_series([1])
    # p_2 = t / (2(t-1)) as series: -(1/2) t * 1/(1-t)
    inv_1mt = one / (one - tser)
    p2 = -(tser * inv_1mt) / 2
    # p_3 = t^2 (t+2) / (6 (t-1)^2 (t+1)):  (t-1)^2 = (1-t)^2 so sign ok
    p3 = (tser ** 2) * (tser + 2) * inv_1mt * inv_1mt * (one / (one + tser)) / 6
    Phi2 = Phi * Phi
    Phi3 = Phi2 * Phi
    Psi2 = p2 * Phi2
    Psi3 = p3 * Phi3
    tau = [to_frac(Psi2[i]) for i in range(0, J + 1)]
    sigma = [to_frac(Psi3[i]) for i in range(0, J + 1)]

    # independent tau check: tau_i = -(1/2) sum_{l<i} [t^l] Phi^2
    B = [to_frac(Phi2[i]) for i in range(0, J + 1)]
    tau_alt_ok = all(tau[i] == -Fraction(1, 2) * sum(B[:i]) for i in range(0, min(J, 400) + 1))

    results = {'NMAX': NMAX, 'N': N, 'J': J, 'tau_closed_form_ok': tau_alt_ok}

    # --- band <= 2 check: c^{(n)}_{2n+i} == A_{n+i} + tau_i for 0 <= i <= n+1 ---
    bad2 = []
    cnt2 = 0
    for n in range(1, NMAX + 1):
        for i in range(0, n + 2):
            k = 2 * n + i
            if k > N: break
            cnt2 += 1
            if C[n][k] != A[n + i] + tau[i]:
                bad2.append((n, i, str(C[n][k]), str(A[n + i] + tau[i])))
    results['band2_checked'] = cnt2
    results['band2_failures'] = bad2

    # --- threshold tightness: at i = n+2, band-2 formula defect == sigma-onset? ---
    # k = 2n + (n+2) = 3n+2, i.e. band-3 line i'=2: defect should be sigma_2 = [t^2]Psi_3
    tight = []
    for n in range(1, NMAX + 1):
        k = 3 * n + 2
        if k > N: break
        defect = C[n][k] - (A[2 * n + 2] + tau[n + 2])
        tight.append((n, str(defect), str(sigma[2]), defect == sigma[2]))
    results['band2_first_defect_equals_sigma2'] = all(x[3] for x in tight)
    results['tightness_rows'] = tight[:6]

    # --- band <= 3 check: c^{(n)}_{3n+i} == A_{2n+i} + tau_{n+i} + sigma_i, 2 <= i <= n+2
    #     (also i in {0,1}: sigma_0=sigma_1=0 since val Psi_3 >= 2; include them) ---
    bad3 = []
    cnt3 = 0
    for n in range(1, NMAX + 1):
        for i in range(0, n + 3):
            k = 3 * n + i
            if k > N: break
            cnt3 += 1
            rhs = A[2 * n + i] + tau[n + i] + sigma[i]
            if C[n][k] != rhs:
                bad3.append((n, i, str(C[n][k]), str(rhs)))
    results['band3_checked'] = cnt3
    results['band3_failures'] = bad3
    results['sigma_val_ok'] = (sigma[0] == 0 and sigma[1] == 0)

    # --- tau/sigma magnitude data ---
    abst = [abs(x) for x in tau]
    abss = [abs(x) for x in sigma]
    results['max_abs_tau'] = str(max(abst)); results['argmax_tau'] = abst.index(max(abst))
    results['max_abs_sigma'] = str(max(abss)); results['argmax_sigma'] = abss.index(max(abss))
    results['tau_leq_half_all'] = all(x <= Fraction(1, 2) for x in abst)
    results['sigma_leq_1_all'] = all(x <= 1 for x in abss)
    results['tau_tail_float'] = [float(tau[i]) for i in range(J - 5, J + 1)]
    results['sigma_tail_float'] = [float(sigma[i]) for i in range(J - 5, J + 1)]
    # band-sum bound check on lines: max |A_{n+i}+tau_i| over checked band-2 cells
    mx = Fraction(0); arg = None
    for n in range(1, NMAX + 1):
        for i in range(0, n + 2):
            if n + i <= J:
                v = abs(A[n + i] + tau[i])
                if v > mx: mx = v; arg = (n, i)
    results['max_abs_band2_value'] = str(mx); results['argmax_band2'] = arg
    results['wall_sec'] = time.time() - t0

    print(json.dumps({k: v for k, v in results.items() if k not in ('tightness_rows',)}, indent=1, default=str))
    print("tightness rows (n, defect, sigma_2, equal):", results['tightness_rows'])
    if len(sys.argv) > 2:
        json.dump(results, open(sys.argv[2], 'w'), default=str)
        print("wrote", sys.argv[2])


if __name__ == "__main__":
    main()
