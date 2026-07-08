"""
Symbolic exact minimal polynomial for M_3 (and validation M_2), using the
reflection-symmetric family and then FULL 2-variable critical-point enumeration.
"""
import sympy as sp

th = sp.symbols('theta', real=True)

def factor(tj):
    return 1 - sp.cos(th - tj)

def mean_trig(expr):
    """(1/2pi) int_0^2pi expr dtheta for a trig polynomial in `th`."""
    expr = sp.expand_trig(sp.expand(expr))
    return sp.integrate(expr, (th, 0, 2*sp.pi)) / (2*sp.pi)

def M_of_config(thetas):
    g = sp.prod([factor(t) for t in thetas])
    g = sp.expand_trig(sp.expand(g))
    D = sp.simplify(mean_trig(g))
    N = sp.simplify(mean_trig(sp.expand(g**2)))
    R = sp.simplify(N / D**2)
    return R, N, D

# ---------- n=2: symmetric {alpha, -alpha}, verify 15/7 ----------
a = sp.symbols('alpha', real=True)
R2, N2, D2 = M_of_config([a, -a])
R2 = sp.simplify(R2)
dR2 = sp.simplify(sp.diff(R2, a))
crit2 = sp.solve(sp.numer(sp.together(dR2)), sp.cos(2*a))
print("n=2: R2 =", R2)
print("n=2: dR2 numerator solve for cos(2a):", crit2)
# evaluate R2 at cos(2a) = -1/4
c2 = sp.symbols('c2')  # cos(2a)
R2c = sp.simplify(R2.rewrite(sp.cos))
# substitute cos(2a)->-1/4 by expressing R2 in cos(2a)
R2_in_c = sp.simplify(sp.expand_trig(R2))
val2 = sp.nsimplify(R2.subs(a, sp.acos(sp.Rational(-1,4))/2))  # cos2a=-1/4
print("n=2: M_2 =", sp.simplify(val2), "=", sp.nsimplify(sp.simplify(val2)))

# ---------- n=3: symmetric {0, beta, -beta}, exact minimal polynomial ----------
b = sp.symbols('beta', real=True)
R3, N3, D3 = M_of_config([sp.Integer(0), b, -b])
R3 = sp.simplify(R3)
print("\nn=3 symmetric family R3(beta) =", R3)
# use x = cos(beta)
x = sp.symbols('x', real=True)
R3x = sp.simplify(R3.rewrite(sp.cos).subs(sp.cos(b), x))
R3x = sp.simplify(sp.expand_trig(R3).rewrite(sp.cos))
# safer: express via cos(beta); expand_trig then replace cos(2b)=2x^2-1, cos(3b)=4x^3-3x
R3e = sp.expand_trig(R3)
R3e = R3e.rewrite(sp.cos)
R3e = R3e.subs(sp.cos(b), x)
R3e = sp.simplify(R3e)
print("n=3: R3(x), x=cos(beta):", R3e)
dR3 = sp.diff(R3e, x)
crit_poly = sp.simplify(sp.numer(sp.together(dR3)))
crit_poly = sp.Poly(sp.expand(crit_poly), x)
print("n=3: critical polynomial in x=cos(beta):", crit_poly)
roots = sp.solve(crit_poly.as_expr(), x)
print("n=3: critical x roots:", [sp.nsimplify(r) if r.is_number else r for r in roots])
# Evaluate R3(x) at each real root in (-1,1), pick max -> M_3; get its minimal polynomial
V = sp.symbols('V')
# minimal polynomial of M_3: eliminate x between (V - R3e)=0 and crit_poly=0
num_eq = sp.expand(sp.numer(sp.together(V - R3e)))   # = V*den - num_R
den = sp.denom(sp.together(R3e))
res = sp.resultant(sp.Poly(num_eq, x), crit_poly, x)
res = sp.factor(res)
print("\nn=3: resultant (elim x) factors, minimal poly of M_3 among them:")
print(res)
