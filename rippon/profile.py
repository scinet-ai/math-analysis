"""Stable-profile coefficients A_j for Rippon 7.54.

Round 1 (finding cda0ff02) proved diagonal stabilization: for every offset
j >= 0 and every n >= max(j,1),  c^{(n)}_{n+j} = A_j, a value independent of n.
Equivalently f_n(t) = t^n * Phi(t) mod t^{2n+1} with Phi(t) = sum_j A_j t^j the
universal profile.  Rippon's conjecture restricted to k <= 2n is exactly
(I): |A_j| <= 1 for all j.

This module computes A_0..A_J EXACTLY (python-flint fmpq_series, no rounding),
by iterating f_n to degree 2J and reading A_j = c^{(J)}_{J+j} (n=J >= j is in the
stable range for all j <= J).  It also RE-VERIFIES stabilization independently:
for each j it checks c^{(n)}_{n+j} is constant over max(j,1) <= n <= J.

Usage:  python -m rippon.profile <J> [out.json]
"""
from __future__ import annotations
import sys, json, time
from fractions import Fraction
from flint import fmpq_series, fmpq, ctx


def to_frac(q) -> Fraction:
    return Fraction(int(q.numer()), int(q.denom()))


def profile(J: int, verify_stab: bool = True):
    """Return (A, stab_ok, worst_stab) where A[j] is the exact Fraction A_j, j=0..J."""
    N = 2 * J
    ctx.cap = N + 1
    t = fmpq_series([0, 1])
    f = (fmpq_series([0, -1])).exp() - 1   # f_1 = e^{-t} - 1
    # Store the full coefficient array only if verifying; else we can stream.
    # We need c^{(n)}_{n+j} for all n>=max(j,1). Keep rows.
    rows = {}   # n -> dict{offset j: fmpq} restricted to j with n+j<=N
    for n in range(1, J + 1):
        row = {}
        for j in range(0, J + 1):
            k = n + j
            if k <= N:
                row[j] = f[k]
        rows[n] = row
        if n < J:
            f = (t * f).exp() - 1

    # A_j from n=J (n>=j for all j<=J)
    A = [to_frac(rows[J][j]) for j in range(0, J + 1)]

    stab_ok = True
    worst = None  # (j, n) where stabilization first-failed, if any
    if verify_stab:
        for j in range(0, J + 1):
            target = rows[J][j]
            n0 = max(j, 1)
            for n in range(n0, J + 1):
                if j in rows[n] and rows[n][j] != target:
                    stab_ok = False
                    worst = (j, n)
                    break
            if not stab_ok:
                break
    return A, stab_ok, worst


def analyze(J: int):
    t0 = time.time()
    A, stab_ok, worst = profile(J)
    wall = time.time() - t0
    # magnitudes
    absA = [abs(a) for a in A]
    # violation of |A_j|<=1?
    violations = [(j, str(A[j])) for j in range(len(A)) if absA[j] > 1]
    # max non-leading (j>=1)
    max_nl = Fraction(0); arg_nl = None
    for j in range(1, len(A)):
        if absA[j] > max_nl:
            max_nl = absA[j]; arg_nl = j
    # signs
    signs = ''.join('+' if a > 0 else ('-' if a < 0 else '0') for a in A)
    # growth: |A_j|^{1/j} for tail
    growth = []
    for j in range(1, len(A)):
        if absA[j] > 0:
            growth.append((j, float(absA[j]) ** (1.0 / j)))
    return {
        'J': J, 'N': 2 * J, 'wall_sec': wall,
        'stabilization_reverified': stab_ok, 'stab_worst': worst,
        'A_leq_1_all': len(violations) == 0,
        'violations_absA_gt_1': violations,
        'max_nonleading_absA': str(max_nl), 'max_nonleading_absA_float': float(max_nl),
        'argmax_nonleading_j': arg_nl,
        'A': [str(a) for a in A],
        'signs': signs,
        'growth_tail': growth[-10:],
    }


if __name__ == "__main__":
    J = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    res = analyze(J)
    print(f"J={res['J']}  N={res['N']}  wall={res['wall_sec']:.1f}s")
    print(f"stabilization re-verified (c^(n)_(n+j) const for n>=max(j,1), n<=J): {res['stabilization_reverified']}  worst={res['stab_worst']}")
    print(f"|A_j|<=1 for ALL j=0..{J}: {res['A_leq_1_all']}   violations: {res['violations_absA_gt_1'][:5]}")
    print(f"max non-leading |A_j| = {res['max_nonleading_absA']} ~ {res['max_nonleading_absA_float']:.8f} at j={res['argmax_nonleading_j']}")
    print("signs A_0..A_{min(60,J)}:")
    print("  " + res['signs'][:61])
    print("first 12 A_j:")
    for j in range(0, min(12, J + 1)):
        print(f"  A_{j:<3d} = {res['A'][j]:>26s}  ~ {float(Fraction(res['A'][j])):+.10f}")
    print("growth |A_j|^(1/j) tail:", [(j, round(g, 5)) for j, g in res['growth_tail']])
    if len(sys.argv) > 2:
        with open(sys.argv[2], "w") as fh:
            json.dump(res, fh)
        print("wrote", sys.argv[2])
