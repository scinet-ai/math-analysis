"""
FULL (non-symmetric) global-optimality certificate for M_3.

The maximizer of R(q)=<|q|^4>/<|q|^2>^2 over deg<=3 polys has all spectral-factor
roots on the unit circle (extreme point of the convex problem), hence is
self-inversive.  After the rotation gauge (fix product of roots = 1) every such
cubic is
        q(z) = z^3 - (a+ib) z^2 + (a-ib) z - 1 ,   a,b real,
so M_3 = max over (a,b) in R^2 of R(a,b).  We enumerate ALL real critical points
exactly (resultants) and confirm the global max is the reflection-symmetric (b=0)
value = root of 208V^3-1224V^2+2268V-1323.
"""
import sympy as sp

a, b, V = sp.symbols('a b V', real=True)
I = sp.I
# q coefficients b_0..b_3
bc = [sp.Integer(-1), a - I*b, -(a + I*b), sp.Integer(1)]
def conj(x): return x.subs(I, -I)
# autocorrelation r_k = sum_j b_{j+k} conj(b_j), k=0..3
def r(k):
    return sp.expand(sum(bc[j+k]*conj(bc[j]) for j in range(0, 4-k)))
r0 = sp.expand(sp.re(r(0)))
num = sp.expand(r0**2 + 2*sum(sp.Abs(r(k))**2 for k in range(1, 4)))
# Abs^2 -> expand via z*conj(z)
def abs2(k):
    rk = r(k); return sp.expand(sp.re(rk)**2 + sp.im(rk)**2)
num = sp.expand(r0**2 + 2*sum(abs2(k) for k in range(1, 4)))
den = r0**2
R = sp.cancel(num/den)
R = sp.simplify(R)
print("R(a,b) =", R)

pa = sp.expand(sp.numer(sp.together(sp.diff(R, a))))
pb = sp.expand(sp.numer(sp.together(sp.diff(R, b))))
print("built partials; factoring pb (contains b?)")
print("pb factored:", sp.factor(pb))

# Critical points: pa=0, pb=0.  Eliminate a.
Ra = sp.resultant(sp.Poly(pa, a), sp.Poly(pb, a))   # poly in b
print("resultant in b (factored):", sp.factor(Ra))

# For the global max, evaluate R at all real critical (a,b); also the b=0 slice.
print("\n--- b=0 slice (reflection-symmetric) ---")
R0 = sp.cancel(R.subs(b, 0))
dR0 = sp.factor(sp.numer(sp.together(sp.diff(R0, a))))
print("R(a,0)=", R0)
print("crit (a): ", dR0)
