"""Coefficient engine for Rippon 7.54 (Hayman-Lingham).

phi_t(z) = e^{tz} - 1.  Iterates phi_t^{n+1}(-1) = e^{t * phi_t^n(-1)} - 1,
with phi_t^1(-1) = e^{-t} - 1.  Each phi_t^n(-1) is a formal power series in t
of valuation n.  We write c^{(n)}_k = [t^k] phi_t^n(-1).

Exact-rational iterates via python-flint fmpq_series (arbitrary-precision Q).
Rigorous ball iterates via python-flint arb_series (real ball arithmetic with
directed rounding); the imaginary part is identically 0 here since all
coefficients are rational, so real arb balls suffice and each ball rigorously
encloses the exact coefficient.
"""
from __future__ import annotations
from flint import fmpq_series, fmpq, ctx


def iterates_exact(N: int, nmax: int):
    """Yield (n, f_n) for n = 1..nmax, f_n an fmpq_series truncated at order N.

    f_n[k] is the EXACT rational c^{(n)}_k.  Uses global ctx.cap = N+1.
    """
    ctx.cap = N + 1
    t = fmpq_series([0, 1])               # the series t
    f = (fmpq_series([0, -1])).exp() - 1  # f_1 = e^{-t} - 1
    for n in range(1, nmax + 1):
        yield n, f
        if n < nmax:
            f = (t * f).exp() - 1         # f_{n+1} = e^{t f_n} - 1


def coeffs_exact(f, N: int):
    """Return list [f[0], ..., f[N]] of exact fmpq coefficients."""
    return [f[k] for k in range(N + 1)]
