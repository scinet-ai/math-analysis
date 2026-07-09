import json, sys
from fractions import Fraction
import math

f = sys.argv[1] if len(sys.argv) > 1 else "data/profile_J300.json"
d = json.load(open(f))
A = [Fraction(s) for s in d['A']]
J = d['J']
absA = [abs(a) for a in A]

# max |A_j| over windows
print("=== magnitude trend ===")
for lo, hi in [(1,10),(10,30),(30,60),(60,100),(100,150),(150,200),(200,250),(250,J)]:
    hi = min(hi, J)
    seg = [(j, absA[j]) for j in range(lo, hi+1)]
    mx = max(seg, key=lambda x: x[1])
    print(f"  j in [{lo:3d},{hi:3d}]: max|A_j| = {float(mx[1]):.6e} at j={mx[0]}")

# overall max non-leading and its location, plus the running max
run = Fraction(0); runj = None
crossings = []
for j in range(1, J+1):
    if absA[j] > run:
        run = absA[j]; runj = j
        crossings.append((j, float(run)))
print("running max of |A_j| (j>=1) updated at:", crossings)

# ratio |A_{j+1}/A_j| tail and sign pattern period search
print("\n=== |A_j| tail values (last 20) ===")
for j in range(J-19, J+1):
    print(f"  A_{j} ~ {float(A[j]):+.8e}")

# Estimate decay exponent: fit log|A_j| ~ a + b*j + c*log j  on tail (j in [150,300], nonzero)
import numpy as np
js = [j for j in range(100, J+1) if absA[j] > 0]
y = np.array([math.log(float(absA[j])) for j in js])
X = np.column_stack([np.ones(len(js)), np.array(js, float), np.log(np.array(js,float))])
coef, *_ = np.linalg.lstsq(X, y, rcond=None)
print(f"\nfit log|A_j| ~ {coef[0]:.3f} + {coef[1]:.5f}*j + {coef[2]:.3f}*log j")
print(f"  => |A_j| ~ C * exp({coef[1]:.5f} j) * j^({coef[2]:.3f});  implied radius R=exp({-coef[1]:.5f})={math.exp(-coef[1]):.4f}")

# Also: is there an oscillation? examine sign run lengths
signs = d['signs']
print("\nsign string (first 80):", signs[:80])
# period of sign pattern via autocorrelation-ish
