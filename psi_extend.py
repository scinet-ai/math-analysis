"""Extend band profiles Psi_2, Psi_3 (tau_i, sigma_i) to i<=600 exactly from the
J=600 profile, plus partial-sum / zero-at-1 diagnostics."""
import json
from fractions import Fraction
from flint import fmpq_series, fmpq, ctx

d = json.load(open("data/profile_J600.json"))
A = [Fraction(s) for s in d['A']]
J = d['J']
ctx.cap = J + 1
Phi = fmpq_series([fmpq(a.numerator, a.denominator) for a in A])
t = fmpq_series([0, 1])
one = fmpq_series([1])
inv1mt = one / (one - t)
p2 = -(t * inv1mt) / 2
p3 = (t ** 2) * (t + 2) * inv1mt * inv1mt * (one / (one + t)) / 6
Phi2 = Phi * Phi
Psi2 = p2 * Phi2
Psi3 = p3 * (Phi2 * Phi)

def fr(q): return Fraction(int(q.numer()), int(q.denom()))

tau = [fr(Psi2[i]) for i in range(J + 1)]
sig = [fr(Psi3[i]) for i in range(J + 1)]
at = [abs(x) for x in tau]; asg = [abs(x) for x in sig]
print(f"tau_i, i<=600: max|tau|={max(at)} at i={at.index(max(at))}; all<=1/2: {all(x<=Fraction(1,2) for x in at)}")
print(f"  |tau| second largest: {sorted(at)[-2]} ; tail floats {[float(tau[i]) for i in (598,599,600)]}")
print(f"sigma_i, i<=600: max|sigma|={max(asg)} at i={asg.index(max(asg))}; all<=1/2: {all(x<=Fraction(1,2) for x in asg)}")
print(f"  tail floats {[float(sig[i]) for i in (598,599,600)]}")

# partial sums of A (coeffs of Phi/(1-t)): sup over J?
ps = Fraction(0); mx = Fraction(0); arg = None
for j in range(J + 1):
    ps += A[j]
    if abs(ps) > mx: mx = abs(ps); arg = j
print(f"partial sums S_J=sum_(j<=J) A_j: max|S_J| = {mx} ~ {float(mx):.6f} at J={arg}; S_600 ~ {float(ps):.6f}")

# Phi near 1: check Phi(r) ~ -2(1-r): ratio Phi(r)/(1-r)
for rr in [Fraction(9,10), Fraction(95,100), Fraction(99,100)]:
    val = sum(A[j] * rr ** j for j in range(J + 1))
    print(f"  Phi({float(rr)}) / (1-r) = {float(val/(1-rr)):.5f}   (simple-zero slope -> -2?)")

# save
json.dump({'tau': [str(x) for x in tau], 'sigma': [str(x) for x in sig]},
          open('data/psi_profiles_600.json', 'w'))
print("wrote data/psi_profiles_600.json")
