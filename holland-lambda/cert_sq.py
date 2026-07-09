"""
Exact M_4, M_5 minimal polynomials using elementary-symmetric coords s=c1+c2, q=c1 c2
(removes the c1<->c2 swap symmetry that collapses naive resultants).
"""
import sympy as sp
from cert_sym import ratio_expr

def minpoly_sq(n, root_at_1, name, near):
    c1, c2, s, q, V = sp.symbols('c1 c2 s q V')
    R, r0, num = ratio_expr([c1, c2], root_at_1)
    R = sp.cancel(R)
    numR, denR = sp.fraction(sp.together(R))
    # rewrite symmetric poly P(c1,c2) in s=c1+c2, q=c1 c2 via c_{1,2}=(s±d)/2, d^2=s^2-4q
    d = sp.symbols('d')
    def sq(expr):
        e = sp.expand(expr.subs({c1: (s + d)/2, c2: (s - d)/2}))
        e = sp.Poly(e, d)
        out = 0
        for (k,), coeff in e.terms():
            assert k % 2 == 0, f"non-symmetric: odd power of d ({k})"
            out += coeff * (s**2 - 4*q)**(k//2)
        return sp.expand(out)
    numRs = sq(numR); denRs = sq(denR)
    Rs = numRs/denRs
    p_s = sp.expand(sp.numer(sp.together(sp.diff(Rs, s))))
    p_q = sp.expand(sp.numer(sp.together(sp.diff(Rs, q))))
    E0 = sp.expand(V*denRs - numRs)
    print(f"\n=== n={n} ({name}) via (s,q) resultants ===", flush=True)
    R1 = sp.resultant(sp.Poly(E0, s), sp.Poly(p_s, s))      # in q,V
    R2 = sp.resultant(sp.Poly(p_q, s), sp.Poly(p_s, s))     # in q
    print("elim s done; R2 zero?", R2 == 0, flush=True)
    Rf = sp.resultant(sp.Poly(R1, q), sp.Poly(R2, q))       # in V
    Rf = sp.factor(Rf)
    print("minimal-poly candidates in V:")
    for f in sp.Mul.make_args(Rf):
        if not f.has(V):
            continue
        pf = sp.Poly(f, V)
        roots = [r for r in sp.nroots(pf.as_expr()) if abs(sp.im(r)) < 1e-9]
        hit = any(abs(float(sp.re(r)) - near) < 1e-6 for r in roots)
        print(f"   deg {pf.degree()}: {sp.expand(f)}"
              f"{'   <== M_%d MINIMAL POLY' % n if hit else ''}")
    return Rf

if __name__ == "__main__":
    import sys
    w = sys.argv[1] if len(sys.argv) > 1 else "4"
    if w == "4":
        minpoly_sq(4, False, "{+-b1,+-b2}", near=3.483450219447)
    else:
        minpoly_sq(5, True, "{0,+-b1,+-b2}", near=4.162256583165)
