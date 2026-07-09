import sympy as sp
from cert_sym import ratio_expr
c1,c2,s,q,V,d=sp.symbols('c1 c2 s q V d')
R,_,_=ratio_expr([c1,c2],True); R=sp.cancel(R)
numR,denR=sp.fraction(sp.together(R))
def sq(expr):
    e=sp.expand(expr.subs({c1:(s+d)/2,c2:(s-d)/2})); e=sp.Poly(e,d); out=0
    for (k,),coeff in e.terms(): out+=coeff*(s**2-4*q)**(k//2)
    return sp.expand(out)
numRs=sq(numR); denRs=sq(denR); common=sp.gcd(sp.numer(sp.together(sp.diff(numRs/denRs,s))),sp.numer(sp.together(sp.diff(numRs/denRs,q))))
ps=sp.numer(sp.together(sp.diff(numRs/denRs,s))); pq=sp.numer(sp.together(sp.diff(numRs/denRs,q)))
A=sp.Poly(sp.expand(sp.cancel(ps/common)),s); B=sp.Poly(sp.expand(sp.cancel(pq/common)),s)
E0=sp.Poly(sp.expand(V*denRs-numRs),s)
print('deg A in s',A.degree(),'deg B in s',B.degree(),flush=True)
rV=sp.resultant(E0,A); print('rV done',flush=True)
rq=sp.resultant(A,B); print('rq done',flush=True)
minV=sp.resultant(sp.Poly(rV,q),sp.Poly(rq,q)); print('minV done',flush=True)
target=4.162256583165
for f,m in sp.factor_list(minV)[1]:
    if not f.has(V): continue
    pf=sp.Poly(f,V); roots=[complex(r) for r in sp.nroots(pf.as_expr())]
    realr=[r.real for r in roots if abs(r.imag)<1e-9]
    if any(abs(rr-target)<1e-6 for rr in realr):
        print('M_5 minimal poly (deg %d):'%pf.degree(), sp.expand(sp.primitive(sp.Poly(f,V))[1]))
        print('  real roots:', sorted(realr))
