import sys, json, math
from fractions import Fraction
import mpmath as mp
mp.mp.dps = 50

d = json.load(open("data/profile_J600.json"))
A = [Fraction(s) for s in d['A']]
J = d['J']

def phi_series(t):
    t = mp.mpc(t); s = mp.mpc(0); tp = mp.mpc(1)
    for a in A:
        s += mp.mpf(a.numerator)/mp.mpf(a.denominator) * tp
        tp *= t
    return s

def phi_orbit(x, kmax=60000):
    t = mp.mpf(x); z = mp.mpf(-1); tk = mp.mpf(1); prev=None
    for k in range(kmax):
        val = z/tk
        if prev is not None and abs(val-prev) < mp.mpf('1e-30')*(1+abs(val)):
            return val, k
        prev=val; z = mp.expm1(t*z); tk*=t
    return prev, kmax

print("=== max |Phi| on circles |t|=r (series, 360 pts) — obstruction check ===", flush=True)
for r in [0.3,0.5,0.7,0.85,0.9,0.95]:
    mx=0; arg=None
    for i in range(360):
        th=2*math.pi*i/360
        v=phi_series(mp.mpf(r)*mp.e**(mp.mpc(0,1)*th))
        if abs(v)>mx: mx=abs(v); arg=th
    print(f"  r={r:.2f}: max|Phi|={mp.nstr(mx,6)} at theta={arg/math.pi:.3f}pi   (>1 => Cauchy/max-modulus cannot bound A_j)", flush=True)

print("\n=== Phi on positive real axis (orbit; product PROVES Phi in (-1,0)) ===", flush=True)
for x in [0.2,0.5,0.8,0.9,0.95,0.99,0.995,0.999]:
    v,k = phi_orbit(x)
    print(f"  t={x:.3f}: Phi={mp.nstr(v,10)}  (k={k})", flush=True)

print("\n=== Phi on negative real axis ===", flush=True)
for x in [-0.2,-0.5,-0.8,-0.9,-0.95,-0.99]:
    v,k = phi_orbit(x)
    tag = "  <-- |Phi|>1" if abs(v)>1 else ""
    print(f"  t={x:.3f}: Phi={mp.nstr(v,10)}  |Phi|={mp.nstr(abs(v),8)}{tag}", flush=True)

print("\n=== singularity exponent A_j ~ C j^(-a): local-slope estimate on envelope ===", flush=True)
absA=[abs(a) for a in A]
# smooth via max over blocks to suppress oscillation, then log-log slope
import numpy as np
blocks=[]
w=20
for lo in range(40,J-w,w):
    seg=[float(absA[j]) for j in range(lo,lo+w) if absA[j]>0]
    if seg: blocks.append((lo+w/2, max(seg)))
xs=np.log([b[0] for b in blocks]); ys=np.log([b[1] for b in blocks])
sl,ic=np.polyfit(xs,ys,1)
print(f"  envelope max|A_j| ~ j^({sl:.3f}) (block-max log-log fit)  => algebraic sing exponent ~ {-sl-1:.2f}", flush=True)
# also Phi(1) via Abel: partial sums
ps=0.0;
for j in range(len(A)): ps+=float(A[j])
print(f"  sum_j A_j (partial, Abel-suggestive of Phi(1)) = {ps:.6f}", flush=True)
