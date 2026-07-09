"""Verify the Koenigs product representation of the stable profile Phi.

CLAIM (derived from Koenigs linearizer K(z)=lim phi_t^n(z)/t^n, telescoped):

    Phi(t) = - prod_{k>=0} h(t f_k(t)),   h(x) = (e^x - 1)/x = sum_{m>=0} x^m/(m+1)!,

where f_0 = -1 and f_{k+1} = e^{t f_k} - 1 (so f_k = phi_t^k(-1)).

Since t f_k has valuation k+1, h(t f_k) = 1 + O(t^{k+1}), so the infinite product
converges t-adically and only k = 0..J-1 affect coefficients through degree J.

This module computes the RHS product to order J with exact rationals (fmpq_series)
and checks it equals Phi (the A_j from rippon.profile) exactly.  A match to high J is
strong verification of the closed-form derivation.

Usage: python -m rippon.product_form <J>
"""
from __future__ import annotations
import sys
from fractions import Fraction
from flint import fmpq_series, fmpq, ctx
from rippon.profile import profile, to_frac


def h_of(u: fmpq_series, N: int) -> fmpq_series:
    """h(u) = sum_{m>=0} u^m/(m+1)!, truncated to order N (u must have positive valuation)."""
    out = fmpq_series([1])           # m=0 term
    term = fmpq_series([1])          # u^0
    fact = 1
    m = 1
    while m <= N:
        fact *= (m + 1)              # (m+1)!
        term = term * u
        # if term is entirely beyond truncation it is 0; flint handles cap
        out = out + term / fact
        m += 1
        # early stop: valuation of u is v>=1, u^m has valuation m*v; once m*val(u) > N, done.
    return out


def phi_product(J: int):
    """Return list of exact Fractions [b_0..b_J] where -prod_{k} h(t f_k) = sum b_j t^j."""
    N = J
    ctx.cap = N + 1
    t = fmpq_series([0, 1])
    # f_0 = -1
    f = fmpq_series([-1])
    prod = fmpq_series([1])
    for k in range(0, J):          # k with k+1 <= J can contribute
        u = t * f                  # t f_k, valuation k+1
        prod = prod * h_of(u, N)
        # advance f -> f_{k+1} = e^{t f_k} - 1  (note: uses same u = t f_k)
        f = u.exp() - 1
    prod = -prod
    return [to_frac(prod[j]) for j in range(0, J + 1)]


if __name__ == "__main__":
    J = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    b = phi_product(J)
    A, _, _ = profile(J)
    mism = [j for j in range(J + 1) if b[j] != A[j]]
    print(f"J={J}: product representation vs exact profile A_j")
    print(f"  first 8 from product: {[str(b[j]) for j in range(min(8,J+1))]}")
    print(f"  mismatches: {len(mism)}  {mism[:10]}")
    print("  MATCH: product formula reproduces Phi exactly" if not mism else "  MISMATCH!")
