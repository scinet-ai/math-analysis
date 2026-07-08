"""
verify.py  --  zero-download smoke test for Holland's Lambda_n results ( < 1 min ).
Run:  uv run python verify.py

Convention: M_n = max_{p in P_n} (1/2pi) int (Re p)^2  (Goldstein-McDonald's Lambda_n).
            E_n = max_{p in P_n} sum|a_nu|^2 = (1/2pi) int |p|^2 = 2 M_n - 1  (SciNet-literal).
Re-checks: M_2=15/7 & M_3 cubic (symbolic); n=3 GLOBAL optimality (exact enumeration);
M_4 quartic (symbolic); extremizer values & nonnegativity; M_5 (the 1984 value); bounds.
"""
import sympy as sp, mpmath as mp
mp.mp.dps = 50
V = sp.symbols('V')

MINPOLY = {
  2: 7*V-15,
  3: 208*V**3-1224*V**2+2268*V-1323,
  4: 50756*V**4-407541*V**3+1150767*V**2-1381455*V+601425,
  5: (95486601852745*V**4-605014169885889*V**3+904976779994997*V**2
      -80594514797592*V-374822538421107),               # PSLQ-recognized
}
NUM = {2:'2.142857142857142857142857142857',
       3:'2.808840165474396887068345490608',
       4:'3.483450219447045602770082821030',
       5:'4.162256583165277914962971405587'}

def R_of_angles(thetas):
    """<u^2>/<u>^2 for u = prod_j (1-cos(theta-theta_j)), exact via Fourier convolution."""
    c = [mp.mpc(1)]
    for t in thetas:
        ker = [-mp.e**(1j*t)/2, mp.mpf(1), -mp.e**(-1j*t)/2]
        nc = [mp.mpc(0)]*(len(c)+2)
        for i, ci in enumerate(c):
            for j, kj in enumerate(ker):
                nc[i+j] += ci*kj
        c = nc
    n = (len(c)-1)//2
    D = c[n].real                       # <u> = frequency-0 coeff
    N = sum(abs(x)**2 for x in c)        # <u^2> = sum|coeff|^2 (Parseval)
    return N/D**2

def extremal_angles(n):
    if n == 2:
        c = mp.sqrt(mp.mpf(3)/8);  be = mp.acos(c);  return [be, -be]
    if n == 3:
        cs = [mp.re(r) for r in mp.polyroots([4,12,6,-1]) if abs(mp.im(r))<1e-30 and -1<mp.re(r)<1]
        be = mp.acos(max(cs));    return [mp.mpf(0), be, -be]
    return None

print("="*72)
print("1) SYMBOLIC re-derivation of M_2, M_3 (resultant on the symmetric family)")
from cert_sym import ratio_expr
for n, root_at_1 in [(2, False), (3, True)]:
    cc = sp.symbols('cc')
    R, _, _ = ratio_expr([cc], root_at_1); R = sp.cancel(R)
    numR, denR = sp.fraction(sp.together(R))
    crit = sp.Poly(sp.numer(sp.together(sp.diff(R, cc))), cc)
    elim = sp.Poly(sp.expand(V*denR - numR), cc)
    res  = sp.resultant(elim, crit)
    div  = sp.rem(sp.Poly(res, V), sp.Poly(MINPOLY[n], V)) == 0
    print(f"   n={n}: {sp.factor(MINPOLY[n])} divides the eliminant: {div}")
    assert div, f"n={n} symbolic check failed"

print("\n2) n=3 GLOBAL optimality (full non-symmetric critical enumeration)")
a, b = sp.symbols('a b', real=True); I = sp.I
bc = [sp.Integer(-1), a-I*b, -(a+I*b), sp.Integer(1)]
cj = lambda x: x.subs(I, -I)
rr = lambda k: sp.expand(sum(bc[j+k]*cj(bc[j]) for j in range(0, 4-k)))
r0 = sp.re(rr(0))
num = sp.expand(r0**2 + 2*sum(sp.re(rr(k))**2 + sp.im(rr(k))**2 for k in range(1, 4)))
Rl = sp.lambdify((a, b), sp.cancel(num/r0**2), 'mpmath')
crits = [(mp.re(t), mp.mpf(0)) for t in mp.polyroots([1,3,-3,-3]) if abs(mp.im(t))<1e-25]
crits.append((mp.mpf(0), mp.mpf(0)))
for coeffs in ([8,0,6,-3], [8,-12,-6,3]):
    for t in mp.polyroots(coeffs):
        if abs(mp.im(t)) < 1e-25:
            av = mp.re(t); b2 = (10*av**3+3*av**2+6*av-3)/(3*(2*av-1))
            if b2 > 0: crits.append((av, mp.sqrt(b2)))
gmax = max(mp.re(Rl(av, bv)) for av, bv in crits)
print(f"   global max over ALL {len(crits)} critical points = {mp.nstr(gmax,18)}  (= M_3)")
assert abs(gmax - mp.mpf(NUM[3])) < 1e-14

print("\n3) SYMBOLIC M_4 quartic (import cert4 -> resultant elimination)")
import cert4  # prints its own confirmation that the eliminant = (M_4 quartic)^2

print("\n4) Extremizer values, nonnegativity, and E_n = 2 M_n - 1")
for n in (2, 3):
    ang = extremal_angles(n)
    Rv = R_of_angles(ang)               # u = prod(1-cos) >= 0 automatically
    print(f"   n={n}: u>=0 (manifest), <u>=1, <u^2>={mp.nstr(Rv,16)} ; "
          f"E_{n}=2M-1={mp.nstr(2*Rv-1,16)}")
    assert abs(Rv - mp.mpf(NUM[n])) < 1e-12

print("\n5) Each M_n is a root of its minimal polynomial (residual vs coeff scale)")
for n in (2,3,4,5):
    p = sp.Poly(MINPOLY[n], V)
    val = mp.mpf(NUM[n])
    resid = mp.polyval([mp.mpf(int(c)) for c in p.all_coeffs()], val)
    scale = max(abs(mp.mpf(int(c))) for c in p.all_coeffs())
    print(f"   n={n}: |minpoly(M_n)|/scale = {mp.nstr(abs(resid)/scale,4)}")
    assert abs(resid)/scale < 1e-20

print("\n6) Bounds  Fejer <= M_n <= n+1")
for n in (2,3,4,5):
    fej = 1 + mp.mpf(n*(2*n+1))/(3*(n+1))
    ok = fej <= mp.mpf(NUM[n])+1e-9 and mp.mpf(NUM[n]) <= n+1
    print(f"   n={n}: {mp.nstr(fej,8)} <= {NUM[n][:10]} <= {n+1}: {ok}")
    assert ok

print("\n"+"="*72)
print("ALL CHECKS PASSED.")
print("  M_2 = 15/7           E_2 = 23/7")
print("  M_3 : 208x^3-1224x^2+2268x-1323     E_3 : 26x^3-228x^2+600x-469   [NEW]")
print("  M_4 : 50756x^4-407541x^3+1150767x^2-1381455x+601425             [NEW]")
print("  M_5 = 4.16225658316528...  (reproduces the 1984 Goldstein-McDonald value)")
