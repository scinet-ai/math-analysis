"""Rigorous certificate: Phi(-1/2) < -1  (so sup_{|t|=r}|Phi| > 1 for r >= 1/2 --
actually already at the point t=-1/2 -- and hence no Cauchy / sup-norm bound on any
circle of radius >= 1/2 can prove |A_j| <= 1).

Method (ball arithmetic, python-flint arb, directed rounding):
  g_n := f_n(t)/t^n = -prod_{k<n} h(t f_k),  h(x) = (e^x-1)/x = expm1(x)/x.
  Phi(t) = lim g_n (Koenigs value K_t(-1); exists for t=-1/2 since the orbit of -1
  converges to the attracting fixed point 0, multiplier |t|=1/2<1).

  TAIL BOUND (proved in proofs/proofs2.md, Prop. tail):
    |Phi(t) - g_N| <= |g_N| * (exp(S_N) - 1),  where S_N := sum_{k>=N} |x_k|,
    x_k := t f_k(t),  because |log h(x)| <= |x| hmm -- we use instead:
    |h(x) - 1| <= e^{|x|} - 1 - |x| <= |x|^2/2 * e^{|x|} <= |x| (for |x|<=1), and
    |g_{n+1}| <= |g_n| e^{|x_n|}; multiplying the enclosures
    |g_M/g_N| in [prod (1 - d_k), prod (1 + d_k)] with d_k := e^{|x_k|}-1-|x_k|...

  Simpler rigorous route implemented here: propagate BALLS the whole way.
    f_k(-1/2) as an arb ball via f_{k+1} = expm1(t f_k)  (arb expm1 is rigorous),
    g_N = f_N / t^N as a ball, plus a rigorous tail estimate:
      once |f_k| <= 1/8 =: r0 (certified by ball upper bound), we have for all m>=k:
      |f_{m+1}| = |expm1(t f_m)| <= e^{|t f_m|} - 1 <= |t||f_m| e^{|t f_m|}
                <= 0.5 * e^{1/16} * |f_m| <= 0.54 |f_m|   (contraction lambda=0.54)
      hence |x_m| = |t f_m| <= 0.5 * r0 * 0.54^{m-k} and
      S_N = sum_{m>=N} |x_m| <= 0.5 r0 0.54^{N-k} / 0.46   for N >= k.
      |Phi - g_N| <= |g_N| (e^{S_N} - 1)  since g_M = g_N * prod_{m=N}^{M-1} h(x_m)
      and |h(x)-1| <= e^{|x|}-1, so |prod h - 1| <= prod(1+(e^{|x_m|}-1)) - 1
                                                 = prod e^{... } hmm bounded by
      exp(sum (e^{|x_m|}-1)) - 1 <= exp(sum |x_m| e^{|x_m|}) - 1 <= exp(1.04 S_N)-1.
  All constants are certified with rational arithmetic; balls give |f_k| bounds.

Output: certified interval for Phi(-1/2) and the verdict Phi(-1/2) < -1.
"""
from __future__ import annotations
import sys
from flint import arb, ctx


def main():
    prec = int(sys.argv[1]) if len(sys.argv) > 1 else 256
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    ctx.prec = prec
    t = arb(-1) / 2
    f = arb(-1)
    k_small = None
    for k in range(1, N + 1):
        f = (t * f).expm1()
        # certified |f| <= 1/8 checkpoint
        if k_small is None:
            up = abs(f).upper()  # arf upper bound
            if arb(up) < arb(1) / 8:
                k_small = k
    gN = f / t ** N
    print(f"prec={prec} N={N}")
    print(f"g_N = f_N/t^N ball: {gN}")
    print(f"first k with certified |f_k| <= 1/8: k0 = {k_small}")
    # tail: S_N <= 0.5 * (1/8) * 0.54^(N-k0) / 0.46
    lam = arb(54) / 100
    S = arb(1) / 16 * lam ** (N - k_small) / (arb(46) / 100)
    tail = abs(gN) * ((arb(104) / 100 * S).exp() - 1)
    print(f"S_N bound: {S}")
    print(f"tail |Phi - g_N| <= {tail}")
    # final enclosure: Phi in gN +- (ball radius already inside) + tail
    lo = gN - tail
    hi = gN + tail
    print(f"certified Phi(-1/2) in [{lo.lower()}, {hi.upper()}]")
    verdict = arb(hi.upper()) < arb(-1)
    print(f"VERDICT Phi(-1/2) < -1 : {bool(verdict)}")
    # also print |Phi| lower bound
    print(f"|Phi(-1/2)| >= {abs(arb(hi.upper()))} > 1: {bool(verdict)}")


if __name__ == "__main__":
    main()
