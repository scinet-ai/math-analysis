"""Numeric probe of the profile Phi(t) = lim_n phi_t^n(-1)/t^n (Koenigs value K(-1)).

Two independent evaluators:
  (S) power-series sum  sum_j A_j t^j  (needs |t| < R = radius of conv.),
  (O) orbit limit       z_0=-1, z_{k+1}=e^{t z_k}-1, Phi = lim z_k / t^k
      (converges when the orbit of -1 tends to 0, i.e. t in the basin of the
       attracting/parabolic fixed point 0 of phi_t).

Cross-checking (S) vs (O) validates both.  Then we map |Phi| on circles |t|=r and
on the real axis, locate the dominant singularity, and estimate its exponent.
"""
import sys, json, cmath, math
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 80

def phi_orbit(t, kmax=20000, tol=mp.mpf('1e-40')):
    """Phi(t) via z_k/t^k limit.  Returns (value, converged, k_used).
    Uses expm1 to avoid catastrophic cancellation in e^{tz}-1 for tiny tz."""
    t = mp.mpc(t)
    if t == 0:
        return mp.mpc(-1), True, 0
    z = mp.mpc(-1)
    tk = mp.mpc(1)      # t^k
    prev = None
    for k in range(0, kmax):
        val = z / tk
        if prev is not None and abs(val - prev) < tol * (1 + abs(val)):
            return val, True, k
        prev = val
        z = mp.expm1(t*z)          # e^{tz} - 1, accurate for small tz
        tk = tk * t
        if abs(z) > mp.mpf('1e12'):
            return prev, False, k
    return prev, False, kmax

def phi_series(t, A):
    s = mp.mpc(0); tp = mp.mpc(1)
    for a in A:
        s += mp.mpf(a.numerator)/mp.mpf(a.denominator) * tp
        tp *= t
    return s

def main():
    d = json.load(open(sys.argv[1] if len(sys.argv)>1 else "data/profile_J300.json"))
    A = [Fraction(s) for s in d['A']]
    J = d['J']

    print("=== cross-check series (S) vs orbit (O) at sample points (|t|<1) ===")
    for t in [0.3, 0.6, 0.9, -0.5, complex(0.5,0.5), complex(0,0.8), complex(-0.7,0.3)]:
        vs = phi_series(t, A)
        vo, conv, k = phi_orbit(t)
        diff = abs(vs - vo) if (vo is not None) else float('nan')
        print(f"  t={str(t):>18s}: |Phi_series|={abs(vs):.6f}  orbit_conv={conv}(k={k})  |S-O|={mp.nstr(diff,3)}")

    print("\n=== |Phi| on real axis ===")
    for x in [0.0,0.2,0.4,0.6,0.8,0.9,0.95,0.99,0.999]:
        v,conv,k = phi_orbit(x)
        vv = v.real if hasattr(v,'real') else v
        print(f"  t=+{x:.3f}: Phi={mp.nstr(vv,8)}  (conv={conv}, k={k})")
    for x in [-0.2,-0.4,-0.6,-0.8,-0.9,-0.95,-0.99]:
        v,conv,k = phi_orbit(x)
        tag = "" if abs(v)<=1 else "   <-- |Phi|>1"
        print(f"  t={x:.3f}: Phi={mp.nstr(v.real,8)}  |Phi|={mp.nstr(abs(v),8)} (conv={conv}){tag}")

    print("\n=== max |Phi| on circles |t|=r (grid 720 pts, orbit evaluator) ===")
    for r in [0.5,0.7,0.85,0.95,0.99]:
        mx=0; arg=None; nconv=0; ntot=0
        for i in range(720):
            th = 2*math.pi*i/720
            t = mp.mpf(r)*mp.e**(mp.mpc(0,1)*th)
            v,conv,k = phi_orbit(t, kmax=6000)
            ntot+=1
            if conv:
                nconv+=1
                if abs(v)>mx: mx=abs(v); arg=th
        print(f"  r={r:.2f}: max|Phi|={mp.nstr(mx,6)} at theta={arg:.3f} ({arg/math.pi:.3f} pi); converged {nconv}/{ntot}")

    print("\n=== behavior approaching t=1 (real) ===")
    for x in [0.9,0.99,0.999,0.9999]:
        v,conv,k = phi_orbit(x, kmax=200000, tol=mp.mpf('1e-20'))
        print(f"  t={x}: Phi~{mp.nstr(v.real,8)} conv={conv} k={k}")

    print("\n=== singularity exponent fit: A_j ~ C j^(-a) (assuming R=1, real sing at t=1) ===")
    absA=[abs(a) for a in A]
    import numpy as np
    js=[j for j in range(150,J+1) if absA[j]>0]
    # local envelope: use running max over window to reduce oscillation
    y=np.array([math.log(float(absA[j])) for j in js]); x=np.log(np.array(js,float))
    b,a0=np.polyfit(x,y,1)
    print(f"  raw fit log|A_j| = {a0:.3f} + ({b:.3f}) log j  => |A_j| ~ j^({b:.3f})")

if __name__=="__main__":
    main()
