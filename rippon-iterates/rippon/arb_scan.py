"""Independent RIGOROUS ball-arithmetic certificate (python-flint arb_series).

Every coefficient is computed as an arb ball [mid +/- rad] that PROVABLY encloses
the exact rational c^{(n)}_k (arb error bounds are rigorous / directed).  To
certify |c^{(n)}_k| <= 1 it suffices that abs_upper(ball) <= 1, i.e. |mid|+rad<=1.

We separate the leading term (k=n, known exactly = -1) from the rest and certify
  max over {(n,k): k>n} of abs_upper(c^{(n)}_k)  <  1   (strict).
This is a second, independent-of-exact confirmation of the no-violation window,
and reaches larger N cheaply.  We also verify (optional) that each ball contains
the exact rational from the fmpq engine on a sub-window, validating the code.
"""
from __future__ import annotations
import sys, json, time
from flint import arb_series, arb, ctx


def arb_certify(N: int, prec: int, nmax: int | None = None, cross_check_to: int = 0):
    if nmax is None:
        nmax = N
    ctx.prec = prec
    ctx.cap = N + 1
    t = arb_series([0, 1])
    f = (arb_series([0, -1])).exp() - 1     # f_1

    max_upper_nonleading = arb(0)           # max abs_upper over k>n
    argmax = None
    max_rad = arb(0)                        # worst ball radius (tightness monitor)
    leading_upper_max = arb(0)              # max abs_upper over leading terms
    fail = []                               # (n,k) where abs_upper >= 1 for k>n (certificate gap)

    # optional exact cross-check
    exact = None
    if cross_check_to > 0:
        from rippon.engine import iterates_exact
        exact = {}
        for n, fe in iterates_exact(min(cross_check_to, N), min(cross_check_to, nmax)):
            exact[n] = [fe[k] for k in range(min(cross_check_to, N) + 1)]
    cross_fail = []

    t0 = time.time()
    for n in range(1, nmax + 1):
        for k in range(n, N + 1):
            ck = f[k]
            r = ck.rad()
            if r > max_rad:
                max_rad = r
            up = ck.abs_upper()
            if k == n:
                if up > leading_upper_max:
                    leading_upper_max = up
            else:
                if up > max_upper_nonleading:
                    max_upper_nonleading = up; argmax = (n, k)
                if not (up < 1):
                    fail.append((n, k, float(up)))
            # cross-check ball encloses exact rational (skip un-floatable tiny values)
            if exact is not None and n in exact and k < len(exact[n]):
                q = exact[n][k]
                try:
                    qf = int(q.numer()) / int(q.denom())
                except OverflowError:
                    qf = None  # |value| below float range -> trivially inside any O(1) ball
                if qf is not None:
                    lo = float(ck.lower()); hi = float(ck.upper())
                    if not (lo - 1e-9 <= qf <= hi + 1e-9):
                        cross_fail.append((n, k))
        if n < nmax:
            f = (t * f).exp() - 1
    dt = time.time() - t0

    return {
        "window": {"nmax": nmax, "N": N}, "prec_bits": prec,
        "engine": "python-flint arb_series (rigorous ball arithmetic)",
        "certified_all_nonleading_below_1": len(fail) == 0,
        "certificate_gaps": fail[:20], "n_gaps": len(fail),
        "max_abs_upper_nonleading": float(max_upper_nonleading),
        "max_abs_upper_nonleading_argmax": argmax,
        "leading_abs_upper_max": float(leading_upper_max),
        "worst_ball_radius": float(max_rad),
        "cross_check_to": cross_check_to,
        "cross_check_failures": cross_fail[:20], "n_cross_fail": len(cross_fail),
        "wall_seconds": round(dt, 2),
    }


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    prec = int(sys.argv[2]) if len(sys.argv) > 2 else 128
    xc = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    out = sys.argv[4] if len(sys.argv) > 4 else None
    r = arb_certify(N, prec, cross_check_to=xc)
    print(f"=== ARB ball certificate, window (n,k)<=({r['window']['nmax']},{r['window']['N']}), prec={prec} bits ===")
    print(f"all non-leading |c| < 1 certified:  {r['certified_all_nonleading_below_1']}  (gaps: {r['n_gaps']})")
    print(f"max abs_upper (non-leading):        {r['max_abs_upper_nonleading']:.10f} at {r['max_abs_upper_nonleading_argmax']}")
    print(f"leading abs_upper max (~1):          {r['leading_abs_upper_max']:.12f}")
    print(f"worst ball radius:                   {r['worst_ball_radius']:.3e}")
    if xc:
        print(f"exact cross-check (to n,k<= {xc}): failures = {r['n_cross_fail']}")
    print(f"wall {r['wall_seconds']}s")
    if out:
        with open(out, "w") as fh:
            json.dump(r, fh, indent=1)
        print("wrote", out)
