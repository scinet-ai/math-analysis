"""
Exact minimal polynomials for M_n via the reflection-symmetric extremal family.

For a reflection-symmetric angle config, the spectral factor q(z)=prod(z-e^{i theta_j})
has REAL coefficients b.  Then
    mean(u) = <|q|^2> = sum_j b_j^2 = r_0
    <u^2>   = <|q|^4> = sum_k r_k^2       (r_k = real autocorrelation of b)
    R = <u^2>/mean^2 = r0^2 + 2*sum_{k>=1} r_k^2  over r0^2.
We build q as a product of quadratic factors (z^2 - 2 c_i z + 1) [conjugate pair at
angle +-beta_i, c_i=cos beta_i] times optional (z-1) [root at angle 0] for odd n.
Maximize R over the c_i in [-1,1]; get exact critical eqs and minimal polynomial of M_n.
"""
import sympy as sp

def build_q_coeffs(cos_syms, root_at_1):
    """Return list of real polynomial coefficients b_0..b_deg of
       q(z) = [ (z-1) if root_at_1 ] * prod_i (z^2 - 2 c_i z + 1)."""
    z = sp.symbols('z')
    q = sp.Integer(1)
    if root_at_1:
        q *= (z - 1)
    for c in cos_syms:
        q *= (z**2 - 2*c*z + 1)
    q = sp.expand(q)
    p = sp.Poly(q, z)
    deg = p.degree()
    b = [p.coeff_monomial(z**k) for k in range(deg+1)]
    return b

def autocorr(b):
    d = len(b) - 1
    r = {}
    for k in range(0, d+1):
        r[k] = sp.expand(sum(b[j+k]*b[j] for j in range(0, d+1-k)))
    return r  # r[0..d]

def ratio_expr(cos_syms, root_at_1):
    b = build_q_coeffs(cos_syms, root_at_1)
    r = autocorr(b)
    d = len(b) - 1
    r0 = r[0]
    num = r0**2 + 2*sum(r[k]**2 for k in range(1, d+1))
    R = sp.simplify(num / r0**2)
    return sp.simplify(R), r0, num

def minpoly_1param(n, root_at_1, name):
    c = sp.symbols('c', real=True)   # cos(beta)
    R, r0, num = ratio_expr([c], root_at_1)
    V = sp.symbols('V')
    R = sp.cancel(R)
    dR = sp.cancel(sp.diff(R, c))
    crit = sp.Poly(sp.numer(sp.together(dR)), c)
    print(f"\n=== n={n} ({name}) ===")
    print("R(c) =", R)
    print("critical poly in c=cos(beta):", sp.factor(crit.as_expr()))
    # minimal polynomial of M_n = R at the relevant critical c
    numR, denR = sp.fraction(sp.together(R))
    elim = sp.Poly(sp.expand(V*denR - numR), c)
    res = sp.resultant(elim, crit, c)
    res = sp.factor(res)
    print("resultant (V):", res)
    # numeric check
    sols = sp.nroots(crit.as_expr())
    best = None
    for s in sols:
        if abs(sp.im(s)) < 1e-12 and -1 <= sp.re(s) <= 1:
            val = float(R.subs(c, sp.re(s)))
            if best is None or val > best[0]:
                best = (val, float(sp.re(s)))
    print("max over real crit c in [-1,1]:  M =", best[0], " at c=", best[1],
          " beta(deg)=", sp.deg(sp.acos(best[1])) if False else round(float(sp.acos(best[1]))*180/3.141592653589793,4))
    return R, crit, res, best

if __name__ == "__main__":
    # n=2: one conjugate pair, no root at 1
    minpoly_1param(2, False, "pair {+-beta}")
    # n=3: root at 1 + one pair
    minpoly_1param(3, True,  "{0, +-beta}")
