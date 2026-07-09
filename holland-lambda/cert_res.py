"""
Faster exact minimal polynomial for M_4, M_5 via RESULTANT chaining (2 params).
"""
import sympy as sp
from cert_sym import ratio_expr

def minpoly_res(n, root_at_1, name, near):
    c1, c2, V = sp.symbols('c1 c2 V')
    R, r0, num = ratio_expr([c1, c2], root_at_1)
    R = sp.cancel(R)
    numR, denR = sp.fraction(sp.together(R))
    p1 = sp.expand(sp.numer(sp.together(sp.diff(R, c1))))
    p2 = sp.expand(sp.numer(sp.together(sp.diff(R, c2))))
    E0 = sp.expand(V*denR - numR)
    print(f"\n=== n={n} ({name}) via resultants ===", flush=True)
    # eliminate c1
    R1 = sp.resultant(sp.Poly(E0, c1), sp.Poly(p1, c1))   # in c2,V
    R2 = sp.resultant(sp.Poly(p2, c1), sp.Poly(p1, c1))   # in c2
    print("eliminated c1", flush=True)
    # eliminate c2
    Rf = sp.resultant(sp.Poly(R1, c2), sp.Poly(R2, c2))   # in V
    print("eliminated c2", flush=True)
    Rf = sp.factor(Rf)
    print("resultant in V (factored):")
    # find the factor with a root near `near`
    facs = sp.Mul.make_args(Rf)
    for f in facs:
        pf = sp.Poly(f, V) if f.has(V) else None
        if pf is None:
            continue
        roots = [r for r in sp.nroots(pf.as_expr()) if abs(sp.im(r))<1e-9]
        hit = any(abs(float(sp.re(r))-near)<1e-6 for r in roots)
        tag = "   <-- MINIMAL POLY of M_%d" % n if hit else ""
        print(f"   deg {pf.degree()}: {sp.expand(f)}{tag}")
    return Rf

if __name__ == "__main__":
    import sys
    w = sys.argv[1] if len(sys.argv)>1 else "4"
    if w=="4":
        minpoly_res(4, False, "{+-b1,+-b2}", near=3.483450219447)
    else:
        minpoly_res(5, True, "{0,+-b1,+-b2}", near=4.162256583165)
