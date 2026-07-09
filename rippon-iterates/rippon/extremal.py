"""Track the NON-LEADING supremum B = sup_{n>=1, k>n} |c^{(n)}_k| and its trend.

The conjecture |c^{(n)}_k| <= 1 splits into:
  (leading)     |c^{(n)}_n| = 1     -- proven (leading lemma, value -1)
  (non-leading) |c^{(n)}_k| <= 1 for k > n
The whole question reduces to: is B <= 1, and does B creep toward 1 as the
window grows?  This script measures the per-row non-leading max M_n = max_{k>n}
|c^{(n)}_k| and the running non-leading sup, to see the trend in n.
"""
from __future__ import annotations
import sys
from fractions import Fraction
from flint import fmpq
from rippon.engine import iterates_exact


def run(N: int, nmax: int | None = None):
    if nmax is None:
        nmax = N
    B = fmpq(0); Barg = None
    rows = []
    for n, f in iterates_exact(N, nmax):
        Mn = fmpq(0); argk = None
        for k in range(n + 1, N + 1):
            a = abs(f[k])
            if a > Mn:
                Mn = a; argk = k
        if argk is not None:
            rows.append((n, Mn, argk, argk - n))
            if Mn > B:
                B = Mn; Barg = (n, argk)
    return B, Barg, rows


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 140
    nmax = N
    B, Barg, rows = run(N)
    fB = float(Fraction(int(B.numer()), int(B.denom())))
    print(f"N={N}  non-leading sup B = {fB:.8f} = {B}  at (n,k)={Barg} (offset {Barg[1]-Barg[0]})")
    print("per-row non-leading max M_n (k>n): n, M_n, argmax k, offset")
    # print a sampling: first 20, then every 10th, to see the trend
    for n, Mn, argk, off in rows:
        if n <= 20 or n % 10 == 0 or n >= nmax - 3:
            fMn = float(Fraction(int(Mn.numer()), int(Mn.denom())))
            print(f"  n={n:4d}  M_n={fMn:.6f}  at k={argk:4d}  offset={off:4d}  (k/n={argk/n:.3f})")
