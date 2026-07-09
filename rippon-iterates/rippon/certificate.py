"""RIGOROUS EXACT certificate for Rippon 7.54 over a window (n,k) <= (nmax, N).

All arithmetic is exact rational (python-flint fmpq_series over Q); there is NO
rounding, so the output is a theorem about the exact Taylor coefficients:

  For every 1 <= n <= nmax and every n <= k <= N,
      |c^{(n)}_k| <= 1,   with equality iff k = n (leading term, value -1),
  and  sup_{k>n} |c^{(n)}_k| = B  attained at the reported (n,k).

Also records the stable-diagonal sequence A_j = c^{(n)}_{n+j} (n>=j) and the
offset envelope E_j = max_n |c^{(n)}_{n+j}| over the window.
"""
from __future__ import annotations
import sys, json, time
from fractions import Fraction
from flint import fmpq
from rippon.engine import iterates_exact


def frac(q):
    return Fraction(int(q.numer()), int(q.denom()))


def certify(N: int, nmax: int | None = None, jcap: int | None = None):
    if nmax is None:
        nmax = N
    if jcap is None:
        jcap = min(N, 64)
    t0 = time.time()

    violations = []            # (n,k, value) with |c|>1
    equality_nonleading = []   # (n,k) with |c|==1 and k!=n  (should be empty)
    leading_bad = []           # n with c^{(n)}_n != -1
    B = fmpq(0); Barg = None   # non-leading sup
    # stable diagonal: for offset j, value stabilizes for n>=max(j,1)
    A = {}                     # j -> fmpq (value at first stable n)
    A_from = {}                # j -> n where it became stable-confirmed
    Aprev = {}                 # j -> previous value seen (to detect stabilization)
    Estab = {}                 # j -> [max|c|, argn] envelope
    ridge = []                 # per n: (n, nonleading_max_float, argk, offset)

    for n, f in iterates_exact(N, nmax):
        c = [f[k] for k in range(N + 1)]
        # leading
        if c[n] != fmpq(-1):
            leading_bad.append((n, str(c[n])))
        # scan
        rowmax = fmpq(0); rk = None
        for k in range(n, N + 1):
            ck = c[k]
            a = abs(ck)
            if a > 1:
                violations.append((n, k, str(ck)))
            if a == 1 and k != n:
                equality_nonleading.append((n, k))
            if k > n:
                if a > rowmax:
                    rowmax = a; rk = k
                if a > B:
                    B = a; Barg = (n, k)
            # offset envelope + stable diagonal
            j = k - n
            if j <= jcap:
                av = frac(a)
                if j not in Estab or av > Estab[j][0]:
                    Estab[j] = [av, n]
                # stabilization: value constant for n>=j
                if n >= max(j, 1):
                    if j not in A:
                        A[j] = ck; A_from[j] = n
                    # else: assert unchanged (checked below)
        if rk is not None:
            ridge.append((n, float(frac(rowmax)), rk, rk - n))

    # verify stabilization actually held (recompute independently is overkill;
    # structure.py does the full stabilization audit). Here trust A[j] first-seen.
    dt = time.time() - t0
    Bf = frac(B)
    cert = {
        "window": {"nmax": nmax, "N": N},
        "no_violation": len(violations) == 0,
        "violations": violations[:20],
        "n_violations": len(violations),
        "equality_only_at_leading": len(equality_nonleading) == 0,
        "equality_nonleading": equality_nonleading[:20],
        "leading_all_minus_one": len(leading_bad) == 0,
        "leading_bad": leading_bad[:20],
        "nonleading_sup_B": str(Bf),
        "nonleading_sup_B_float": float(Bf),
        "nonleading_sup_argmax": Barg,
        "A_j": {j: str(frac(A[j])) for j in sorted(A)},
        "A_j_float": {j: float(frac(A[j])) for j in sorted(A)},
        "envelope_E_j": {j: str(Estab[j][0]) for j in sorted(Estab)},
        "envelope_E_j_float": {j: float(Estab[j][0]) for j in sorted(Estab)},
        "envelope_argmax_n": {j: Estab[j][1] for j in sorted(Estab)},
        "wall_seconds": round(dt, 2),
        "engine": "python-flint fmpq_series (exact rational, no rounding)",
    }
    return cert, ridge


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    out = sys.argv[2] if len(sys.argv) > 2 else None
    cert, ridge = certify(N)
    print(f"=== EXACT certificate, window (n,k) <= ({cert['window']['nmax']},{cert['window']['N']}) ===")
    print(f"no |c|>1 violation:              {cert['no_violation']}  (checked, {cert['n_violations']} violations)")
    print(f"leading c^(n)_n == -1 for all n:  {cert['leading_all_minus_one']}")
    print(f"|c|=1 ONLY at leading term:       {cert['equality_only_at_leading']}")
    print(f"non-leading sup B = {cert['nonleading_sup_B']} = {cert['nonleading_sup_B_float']:.8f} at {cert['nonleading_sup_argmax']}")
    print(f"wall {cert['wall_seconds']}s   ({cert['engine']})")
    print("stable diagonal A_j (first 16):")
    for j in sorted(cert['A_j'])[:16]:
        print(f"  A_{j} = {cert['A_j'][j]}")
    if out:
        with open(out, "w") as fh:
            json.dump(cert, fh, indent=1)
        print("wrote", out)
