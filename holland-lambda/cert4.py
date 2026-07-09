import sympy as sp
from cert_sym import ratio_expr
c1,c2,s,q,V,d=sp.symbols('c1 c2 s q V d')
R,r0,num=ratio_expr([c1,c2],False); R=sp.cancel(R)
numR,denR=sp.fraction(sp.together(R))
def sq(expr):
    e=sp.expand(expr.subs({c1:(s+d)/2,c2:(s-d)/2})); e=sp.Poly(e,d); out=0
    for (k,),coeff in e.terms(): out+=coeff*(s**2-4*q)**(k//2)
    return sp.expand(out)
numRs=sq(numR); denRs=sq(denR); common=8*q**2+8*q+4*s**2+3
ps=sp.numer(sp.together(sp.diff(numRs/denRs,s))); pq=sp.numer(sp.together(sp.diff(numRs/denRs,q)))
A=sp.expand(sp.cancel(sp.cancel(ps/common)/(8*s)))
B=sp.expand(sp.cancel(pq/common))
E0=sp.expand(V*denRs-numRs)

# Step 1: eliminate s from {E0=0, A=0}
rV=sp.resultant(sp.Poly(E0,s),sp.Poly(A,s))   # in q,V
print('rV has V:',rV.has(V),' deg_q:',sp.Poly(rV,q).degree(),' deg_V:',sp.Poly(rV,V).degree())
# Step 2: eliminate s from {A=0,B=0}
rq=sp.resultant(sp.Poly(A,s),sp.Poly(B,s))     # in q
rq_facs=[f for f,_ in sp.factor_list(rq)[1]]
print('rq factors (in q):',[sp.expand(f) for f in rq_facs])
# Step 3: for each q-factor, eliminate q from rV and that factor -> candidate minpoly(V)
target=3.483450219447
for f in rq_facs:
    if not f.has(q): continue
    mV=sp.resultant(sp.Poly(rV,q),sp.Poly(f,q))
    mV=sp.factor(mV)
    for g in sp.Mul.make_args(mV):
        if not g.has(V): continue
        pg=sp.Poly(g,V); roots=[r for r in sp.nroots(pg.as_expr()) if abs(sp.im(r))<1e-9]
        hit=any(abs(float(sp.re(r))-target)<1e-6 for r in roots)
        if hit:
            print('M_4 MINIMAL POLY (deg %d):'%pg.degree(), sp.expand(sp.primitive(g)[1]))
            # verify against PSLQ quartic
            pslq=50756*V**4-407541*V**3+1150767*V**2-1381455*V+601425
            print('  matches PSLQ quartic:', sp.simplify(sp.gcd(sp.Poly(g,V),sp.Poly(pslq,V)).as_expr())!=1)
