"""Region (II) structural checks: express the first transient lines k=2n+i via the
profile A_j.  Derived identity 1 (proved from the recurrence):

    c^{(n)}_{2n+1} = A_{n+1} - 1/2      (all n>=1).

We verify it exactly, and search for the analogous closed form on the next line
k=2n+2 by fitting c^{(n)}_{2n+2} against {A_m, products} in exact arithmetic.
"""
import sys
from fractions import Fraction
from flint import fmpq_series, fmpq, ctx
from rippon.engine import iterates_exact
from rippon.profile import profile, to_frac

def main():
    NMAX = 30
    N = 3*NMAX + 6
    C = {}
    for n, f in iterates_exact(N, NMAX):
        C[n] = [to_frac(f[k]) for k in range(N+1)]
    A, _, _ = profile(2*NMAX+6)   # A_0..A_(2NMAX+6)

    print("=== identity 1: c^{(n)}_{2n+1} == A_{n+1} - 1/2 ===")
    ok1 = True
    for n in range(1, NMAX+1):
        lhs = C[n][2*n+1]
        rhs = A[n+1] - Fraction(1,2)
        mark = "" if lhs==rhs else "  MISMATCH"
        if lhs!=rhs: ok1=False
        if n<=8 or lhs!=rhs:
            print(f"  n={n:2d}: c^(n)_(2n+1)={str(lhs):>14s}  A_(n+1)-1/2={str(rhs):>14s}{mark}")
    print(f"  identity 1 holds for all n<=%d: {ok1}" % NMAX)

    print("\n=== line k=2n+2: values c^{(n)}_{2n+2} and A_{n+2} ===")
    for n in range(1, 12):
        c = C[n][2*n+2]
        print(f"  n={n:2d}: c^(n)_(2n+2)={str(c):>18s} ~{float(c):+.6f}   A_(n+2)={str(A[n+2]):>14s}   c-A_(n+2)={str(c-A[n+2]):>16s}")

    # Try: c^{(n)}_{2n+2} = A_{n+2} + alpha*A_1 + ... ; empirically inspect c - A_{n+2}
    print("\n  (difference d_n := c^(n)_(2n+2) - A_(n+2); look for pattern in n)")
    for n in range(1, 12):
        d = C[n][2*n+2]-A[n+2]
        print(f"    n={n:2d}: d_n={str(d):>18s} ~{float(d):+.6f}")

if __name__=="__main__":
    main()
