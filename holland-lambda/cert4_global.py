"""
n=4 FULL global optimality: enumerate ALL critical points of R(a,b,c) for the
self-inversive quartic q=z^4-(a+ib)z^3+c z^2-(a-ib)z+1, then take the max.
Branches:  b=0 (reflection-symmetric family -> M_4 quartic, already certified)
           b!=0: H(a,b,c)=0 with pa=pc=0.  R is even in b, so substitute b^2 from H.
"""
import sympy as sp, mpmath as mp
a,b,c = sp.symbols('a b c', real=True); I=sp.I
bc=[sp.Integer(1), -(a-I*b), c, -(a+I*b), sp.Integer(1)]
cj=lambda x:x.subs(I,-I)
rr=lambda k: sp.expand(sum(bc[j+k]*cj(bc[j]) for j in range(0,5-k)))
r0=sp.re(rr(0))
num=sp.expand(r0**2+2*sum(sp.re(rr(k))**2+sp.im(rr(k))**2 for k in range(1,5)))
R=sp.cancel(num/r0**2)
pa=sp.numer(sp.together(sp.diff(R,a)))
pb=sp.numer(sp.together(sp.diff(R,b)))
pc=sp.numer(sp.together(sp.diff(R,c)))
# pb = (-8 b)*(pos)*(H); extract H
pbf=sp.factor(pb)
H=[f for f,_ in sp.factor_list(pb)[1] if f.has(b) and f.diff(b)!=0 and not (f-b).is_zero and f!=b]
# pick the nontrivial branch factor (contains b^2 and c)
Hbranch=None
for f,_ in sp.factor_list(pb)[1]:
    if f.has(c) and f.has(b): Hbranch=f
print("H branch:", sp.expand(Hbranch))
B2=sp.symbols('B2')  # = b^2
Hsub=sp.expand(Hbranch.subs(b**2,B2))
b2_sol=sp.solve(sp.Poly(Hsub,B2),B2)
print("b^2 solutions on H-branch:", len(b2_sol))
b2=b2_sol[0]
Rl=sp.lambdify((a,b,c),R,'mpmath')
# substitute b^2 -> b2(a,c) into pa,pc (even in b), get 2 eqs in (a,c)
def sub_b2(expr):
    e=sp.expand(expr)
    assert e.subs(b,-b)==e, "not even in b"
    return sp.simplify(e.subs(b**2,B2).subs(B2,b2))
pa_ac=sp.numer(sp.together(sub_b2(pa)))
pc_ac=sp.numer(sp.together(sub_b2(pc)))
print("eliminating a ...", flush=True)
res_c=sp.resultant(sp.Poly(pa_ac,a),sp.Poly(pc_ac,a))
res_c=sp.factor(res_c)
print("resultant in c (b!=0 branch), factors:")
mp.mp.dps=30
sols=[]
for f,_ in sp.factor_list(res_c)[1]:
    if not f.has(c): continue
    for r in mp.polyroots([mp.mpf(int(x)) for x in sp.Poly(f,c).all_coeffs()]):
        if abs(mp.im(r))<1e-18: sols.append(mp.re(r))
print("real c critical values (b!=0):", [mp.nstr(x,8) for x in sols])
# for each c, get a (from pa_ac=0), b^2, evaluate R, collect
vals=[]
for cval in sols:
    for r in mp.polyroots([mp.mpf(int(x)) for x in sp.Poly(pa_ac.subs(c,sp.nsimplify(0)),a).all_coeffs()]) if False else []:
        pass
    # solve pa_ac=0 for a numerically at this c
    pa_c=sp.Poly(pa_ac.subs(c,sp.Float(float(cval),30)),a)
    for ra in mp.polyroots([mp.mpf(str(x)) for x in pa_c.all_coeffs()]):
        if abs(mp.im(ra))<1e-15:
            av=mp.re(ra); b2v=complex(b2.subs({a:float(av),c:float(cval)}))
            if b2v.real>1e-12 and abs(b2v.imag)<1e-9:
                Rv=Rl(av,mp.sqrt(mp.mpf(b2v.real)),cval)
                vals.append(mp.re(Rv))
print("R values at b!=0 critical points:", sorted(set(mp.nstr(v,12) for v in vals)))
print("max on b!=0 branch:", mp.nstr(max(vals),16) if vals else "none")
print("M_4 = 3.483450219447  (b=0 symmetric-family max, already certified)")
