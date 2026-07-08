"""
Exact minimal polynomials for M_4, M_5 via reflection-symmetric families (2 params).
n=4: q = (z^2-2 c1 z+1)(z^2-2 c2 z+1)          [pairs +-b1, +-b2]
n=5: q = (z-1)(z^2-2 c1 z+1)(z^2-2 c2 z+1)      [0, +-b1, +-b2]
Maximize R(c1,c2); eliminate c1,c2 by Groebner to get minimal polynomial of M_n.
"""
import sympy as sp
from cert_sym import ratio_expr

def minpoly_2param(n, root_at_1, name, expected=None):
    c1, c2 = sp.symbols('c1 c2', real=True)
    V = sp.symbols('V')
    R, r0, num = ratio_expr([c1, c2], root_at_1)
    R = sp.cancel(R)
    numR, denR = sp.fraction(sp.together(R))
    # critical equations (numerators of partials)
    p1 = sp.Poly(sp.numer(sp.together(sp.diff(R, c1))), c1, c2)
    p2 = sp.Poly(sp.numer(sp.together(sp.diff(R, c2))), c1, c2)
    E0 = sp.Poly(sp.expand(V*denR - numR), c1, c2)
    print(f"\n=== n={n} ({name}) ===")
    # Eliminate c1, c2 via Groebner (lex, c1>c2>V) then take poly in V only
    G = sp.groebner([E0.as_expr(), p1.as_expr(), p2.as_expr()], c1, c2, V,
                    order='lex')
    vpolys = []
    for g in G.exprs:
        fs = g.free_symbols
        if fs <= {V}:
            vpolys.append(sp.factor(g))
    print("V-only Groebner elements (factored):")
    for vp in vpolys:
        print("   ", vp)
    if expected is not None:
        print("expected PSLQ minpoly:", expected)
    return vpolys

if __name__ == "__main__":
    import sys
    which = sys.argv[1] if len(sys.argv) > 1 else "4"
    if which == "4":
        minpoly_2param(4, False, "{+-b1,+-b2}",
            expected="50756 V^4 -407541 V^3 +1150767 V^2 -1381455 V +601425")
    elif which == "5":
        minpoly_2param(5, True, "{0,+-b1,+-b2}")
