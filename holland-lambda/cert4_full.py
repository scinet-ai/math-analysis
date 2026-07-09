"""
FULL (non-symmetric) global-optimality attempt for M_4.
Self-inversive quartic (all maximizer roots on circle, product=1, rotation gauge):
   q(z) = z^4 - (a+ib) z^3 + c z^2 - (a-ib) z + 1 ,  a,b,c real.
Maximize R(a,b,c)=<|q|^4>/<|q|^2>^2 over R^3.  Reflection symmetry = {b=0}.
Enumerate real critical points exactly; confirm global max = M_4 at b=0.
"""
import sympy as sp
import mpmath as mp

a,b,c = sp.symbols('a b c', real=True); I=sp.I
bc=[sp.Integer(1), -(a-I*b), c, -(a+I*b), sp.Integer(1)]  # b_0..b_4 of q (z^0..z^4)
def conj(x): return x.subs(I,-I)
def r(k): return sp.expand(sum(bc[j+k]*conj(bc[j]) for j in range(0,5-k)))
r0=sp.expand(sp.re(r(0)))
def abs2(k):
    rk=r(k); return sp.expand(sp.re(rk)**2+sp.im(rk)**2)
num=sp.expand(r0**2+2*sum(abs2(k) for k in range(1,5))); den=r0**2
R=sp.cancel(num/den)
print("R built; simplifying partials...", flush=True)
pa=sp.expand(sp.numer(sp.together(sp.diff(R,a))))
pb=sp.expand(sp.numer(sp.together(sp.diff(R,b))))
pc=sp.expand(sp.numer(sp.together(sp.diff(R,c))))
print("pb factored:", sp.factor(pb), flush=True)
# pb should factor as (stuff)*b*(branch); off-symmetric needs branch=0
# Enumerate: eliminate to get critical points. Try resultant elimination a,c on b=0 slice first.
print("\n--- b=0 slice ---", flush=True)
R0=sp.cancel(R.subs(b,0))
pa0=sp.numer(sp.together(sp.diff(R0,a))); pc0=sp.numer(sp.together(sp.diff(R0,c)))
# eliminate c then a
r_ac=sp.resultant(sp.Poly(pa0,c),sp.Poly(pc0,c))
print("b=0: resultant in a (factored):", sp.factor(r_ac), flush=True)

# full: pb = b * G(a,b,c) * (positive).  off-branch: G=0 with pa=pc=0.
Rf=sp.lambdify((a,b,c),R,'mpmath')
mp.mp.dps=30
# numeric global check via many restarts on R^3 (bounded, ->3/2 at infinity)
from scipy.optimize import minimize
import numpy as np
Rn=sp.lambdify((a,b,c),R,'numpy')
best=-1;ba=None
rng=np.random.default_rng(0)
for _ in range(400):
    x0=rng.uniform(-4,4,3)
    res=minimize(lambda x:-float(Rn(*x)),x0,method='Nelder-Mead',
                 options={'xatol':1e-10,'fatol':1e-12,'maxiter':5000})
    if -res.fun>best: best=-res.fun; ba=res.x
print("\nnumeric global max over R^3:", best, " at (a,b,c)=",ba)
print("M_4 target 3.483450219447 ; |b| at opt =", abs(ba[1]))
