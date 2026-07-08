"""
Exact minimal polynomial of M_n for the reflection-symmetric 2-pair family,
n=4 (root_at_1=False) and n=5 (root_at_1=True), via (s,q)=(c1+c2, c1 c2) resultants
with automatic stripping of the (always-positive) denominator common factor.
"""
import sympy as sp
from cert_sym import ratio_expr

def minpoly(n, root_at_1, near):
    c1,c2,s,q,V,d = sp.symbols('c1 c2 s q V d')
    R,_,_ = ratio_expr([c1,c2], root_at_1); R=sp.cancel(R)
    numR,denR = sp.fraction(sp.together(R))
    def sq(expr):
        e=sp.expand(expr.subs({c1:(s+d)/2,c2:(s-d)/2})); e=sp.Poly(e,d); out=0
        for (k,),coeff in e.terms(): out+=coeff*(s**2-4*q)**(k//2)
        return sp.expand(out)
    numRs=sq(numR); denRs=sq(denR)
    ps=sp.numer(sp.together(sp.diff(numRs/denRs,s)))
    pq=sp.numer(sp.together(sp.diff(numRs/denRs,q)))
    common=sp.gcd(ps,pq)
    Ps=sp.cancel(ps/common); Pq=sp.cancel(pq/common)
    # Ps typically has a factor s; keep full (s=0 branch handled by resultant)
    A=sp.Poly(sp.expand(sp.cancel(Ps/sp.gcd(Ps,s))) if Ps.has(s) else Ps, s)
    # safer: just use Ps, Pq directly
    A=sp.Poly(sp.expand(Ps), s); B=sp.Poly(sp.expand(Pq), s)
    E0=sp.Poly(sp.expand(V*denRs-numRs), s)
    rV=sp.resultant(E0, A)           # in q,V
    rq=sp.resultant(A, B)            # in q
    minV=sp.resultant(sp.Poly(rV,q), sp.Poly(rq,q))
    print(f"\n=== n={n} ===", flush=True)
    for f,m in sp.factor_list(minV)[1]:
        if not f.has(V): continue
        pf=sp.Poly(f,V); roots=[complex(r) for r in sp.nroots(pf.as_expr())]
        realr=[r.real for r in roots if abs(r.imag)<1e-9]
        hit=any(abs(rr-near)<1e-6 for rr in realr)
        if hit:
            prim=sp.primitive(sp.Poly(f,V))[1]
            print(f"  M_{n} minimal poly (deg {pf.degree()}): {sp.expand(prim)}")
            print(f"     real roots: {sorted(realr)}")
            print(f"     M_{n} = {max(r for r in realr if abs(r-near)<1e-6):.15f}")

if __name__=="__main__":
    import sys
    if "4" in sys.argv: minpoly(4, False, 3.483450219447)
    if "5" in sys.argv: minpoly(5, True, 4.162256583165)
