"""Extremal-structure map for Rippon 7.54.

Computes, from the exact-rational iterate array c^{(n)}_k = [t^k] phi_t^n(-1):

  * stabilization: for each offset j = k - n, the sequence n |-> c^{(n)}_{n+j};
    verify it becomes constant for n >= j and record the stable value A_j.
  * offset envelope E_j = max_n |c^{(n)}_{n+j}| (sup over the scanned window).
  * transient max: max_{k > 2n} |c^{(n)}_k| per n and overall.
  * global max |c| and its argmax.
"""
from __future__ import annotations
import sys, json
from fractions import Fraction
from flint import fmpq
from rippon.engine import iterates_exact


def to_frac(q: fmpq) -> Fraction:
    return Fraction(int(q.numer()), int(q.denom()))


def analyze(N: int, nmax: int | None = None):
    if nmax is None:
        nmax = N
    C = {}  # C[n] = list of fmpq coeffs 0..N
    for n, f in iterates_exact(N, nmax):
        C[n] = [f[k] for k in range(N + 1)]

    # --- stabilization + A_j ---
    # For offset j, gather (n, c^{(n)}_{n+j}) for all n with n+j <= N.
    stab = {}   # j -> {'A': Fraction or None, 'stable_from': n0, 'seq': [(n, Fraction)]}
    max_j = N - 1
    for j in range(0, max_j + 1):
        seq = []
        for n in range(1, nmax + 1):
            k = n + j
            if k > N:
                break
            seq.append((n, C[n][k]))
        if not seq:
            continue
        # theory: constant for n >= j. verify and record A_j = value at n=max(j,1).
        A = None
        stable_from = None
        # find the tail-constant value: take last value, check how far back it holds
        tailval = seq[-1][1]
        # walk back while equal
        idx = len(seq) - 1
        while idx - 1 >= 0 and seq[idx - 1][1] == tailval:
            idx -= 1
        stable_from = seq[idx][0]
        A = tailval
        stab[j] = {
            'A': to_frac(A),
            'stable_from_n': stable_from,
            'theory_threshold_n': max(j, 1),
            'n_samples': len(seq),
        }

    # verify stabilization theory: stable_from_n <= theory threshold max(j,1)
    # (i.e., value is already stable by n = max(j,1); may stabilize earlier)
    stab_ok = all(v['stable_from_n'] <= v['theory_threshold_n'] for v in stab.values())

    # --- offset envelope E_j = max_n |c^{(n)}_{n+j}| ---
    env = {}
    for j in range(0, max_j + 1):
        best = fmpq(0); argn = None
        for n in range(1, nmax + 1):
            k = n + j
            if k > N:
                break
            a = abs(C[n][k])
            if a > best:
                best = a; argn = n
        env[j] = {'E': to_frac(best), 'argmax_n': argn}

    # --- transient max (k > 2n) and global max ---
    gmax = fmpq(0); garg = None
    tmax = fmpq(0); targ = None
    per_n_max = {}
    violations = []
    for n in range(1, nmax + 1):
        row_max = fmpq(0); row_arg = None
        for k in range(n, N + 1):
            a = abs(C[n][k])
            if a > row_max:
                row_max = a; row_arg = k
            if a > gmax:
                gmax = a; garg = (n, k)
            if a > 1:
                violations.append((n, k, str(C[n][k])))
            if k > 2 * n and a > tmax:
                tmax = a; targ = (n, k)
        per_n_max[n] = (to_frac(row_max), row_arg)

    return {
        'N': N, 'nmax': nmax,
        'global_max': str(to_frac(gmax)), 'global_argmax': garg,
        'global_max_is_leading': (garg is not None and garg[0] == garg[1]),
        'violations': violations,
        'transient_max_over_2n': str(to_frac(tmax)), 'transient_argmax': targ,
        'stabilization_theory_holds': stab_ok,
        'A': {j: {'A': str(v['A']), 'A_float': float(v['A']),
                  'stable_from_n': v['stable_from_n'],
                  'theory_threshold_n': v['theory_threshold_n']}
              for j, v in stab.items()},
        'envelope': {j: {'E': str(v['E']), 'E_float': float(v['E']), 'argmax_n': v['argmax_n']}
                     for j, v in env.items()},
        'per_n_max': {n: (str(v[0]), v[1]) for n, v in per_n_max.items()},
    }


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    res = analyze(N)
    # print a compact summary
    print(f"N={res['N']} nmax={res['nmax']}")
    print(f"global max|c| = {res['global_max']} at {res['global_argmax']} (leading={res['global_max_is_leading']})")
    print(f"violations |c|>1: {len(res['violations'])} {res['violations'][:3]}")
    print(f"transient max (k>2n) = {res['transient_max_over_2n']} ~ {float(Fraction(res['transient_max_over_2n'])):.6f} at {res['transient_argmax']}")
    print(f"stabilization theory (stable by n=max(j,1)) holds: {res['stabilization_theory_holds']}")
    print("stable diagonal A_j = c^{(j)}_{2j}:")
    for j in range(0, min(24, len(res['A']))):
        a = res['A'][j]
        e = res['envelope'][j]
        print(f"  j={j:3d}  A_j={a['A']:>22s} |A_j|={abs(a['A_float']):.6e}   E_j(envelope)={e['E']:>22s} |E_j|={e['E_float']:.6e}")
    if len(sys.argv) > 2:
        with open(sys.argv[2], "w") as fh:
            json.dump(res, fh, indent=0)
        print("wrote", sys.argv[2])
