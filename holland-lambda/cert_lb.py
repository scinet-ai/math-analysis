"""
cert_lb.py -- CERTIFIED feasible-point lower bounds on M_n at large n.

M_n = max over nonneg trig polys u (deg<=n, mean 1) of <u^2>.
A feasible point gives a rigorous LOWER bound: any explicit nonneg mean-1 trig
poly u of degree <=n has <u^2> <= M_n.

We use reflection-symmetric extreme-point configs u ~ prod_j (1-cos(theta-theta_j))
with the angles in +-pairs (and optionally one at 0). For such a config the
Fourier (cosine) coefficients of g are RATIONAL in the cosines c_j=cos(theta_j).
Snapping each c_j to a rational p_j/Q_j and clearing denominators, g has INTEGER
Fourier coefficients (up to a global scale that cancels in the scale-invariant
ratio <g^2>/<g>^2). Hence

        M_n  >=  V_n := <g^2>/<g>^2  =  (sum_k A_k^2) / A_0^2        (exact rational)

where A_k are the integer Fourier coefficients of g. This is fully certified in
exact integer arithmetic -- no global optimality, no floating point in the bound.

Per +-pair (cos = c): the trig-degree-2 cosine kernel is
    [1/4, -c, 1/2 + c^2, -c, 1/4]     (frequencies -2..2)
times 4Q^2 (c=p/Q) it becomes the INTEGER kernel
    [Q^2, -4 p Q, 2Q^2 + 4 p^2, -4 p Q, Q^2].
A zero-angle factor (1-cos theta) has integer kernel [-1, 2, -1] (freqs -1,0,1).
"""
from fractions import Fraction
import numpy as np
from scipy.optimize import minimize

# ---------------- numeric reflection-symmetric optimizer (to locate good angles) ------------
def R_and_grad_sym(halfangles, has_zero, Mgrid):
    """u = [(1-cos(th)) if has_zero] * prod_j (1-cos(th-a_j))(1-cos(th+a_j)); maximize <u^2>/<u>^2."""
    a = halfangles
    th = np.linspace(0, 2*np.pi, Mgrid, endpoint=False)
    logF = np.zeros(Mgrid)
    if has_zero:
        logF += np.log(np.clip(1-np.cos(th), 1e-300, None))
    diffp = th[None,:]-a[:,None]; diffm = th[None,:]+a[:,None]
    Fp = 1-np.cos(diffp); Fm = 1-np.cos(diffm)
    logF += np.log(np.clip(Fp,1e-300,None)).sum(0) + np.log(np.clip(Fm,1e-300,None)).sum(0)
    u = np.exp(logF)
    U1 = u.mean(); U2 = (u*u).mean(); R = U2/U1**2
    grad = np.zeros(len(a))
    for j in range(len(a)):
        # d/d a_j of logu = [-sin(th-a_j)*(-1)]/Fp_j *(-1)... compute directly:
        # d/da_j log Fp_j = d/da_j log(1-cos(th-a_j)) = -sin(th-a_j)/(1-cos(th-a_j))
        # d/da_j log Fm_j = +sin(th+a_j)/(1-cos(th+a_j))
        dlog = -np.sin(diffp[j])/np.clip(Fp[j],1e-300,None) + np.sin(diffm[j])/np.clip(Fm[j],1e-300,None)
        du = u*dlog
        dU1 = du.mean(); dU2 = (2*u*du).mean()
        grad[j] = (dU2*U1**2 - U2*2*U1*dU1)/U1**4
    return R, grad

def optimize_sym(n, restarts=10, seed=3):
    has_zero = (n % 2 == 1)
    k = n//2  # number of +- pairs
    if k == 0:                      # n<=1: only the zero-angle factor (or empty)
        M = 1.5 if n == 1 else 1.0  # M_1=<(1-cos)^2>=3/2 ; M_0=1
        return M, np.array([]), has_zero
    Mgrid = max(20*n, 1024)
    rng = np.random.default_rng(seed)
    best=-1; ba=None
    for r in range(restarts):
        if r==0:
            a0 = np.linspace(0.15, np.pi-0.05, k)      # spread
        elif r==1:
            a0 = np.linspace(0.02, 1.2, k)             # clustered near 0 (extremizers concentrate)
        else:
            a0 = np.sort(rng.uniform(0.02, np.pi-0.02, k))
        def obj(a):
            R,g = R_and_grad_sym(a, has_zero, Mgrid)
            return -R, -g
        res = minimize(obj, a0, jac=True, method='L-BFGS-B',
                       bounds=[(1e-3, np.pi-1e-3)]*k,
                       options={'maxiter':4000,'ftol':1e-15,'gtol':1e-11})
        R = -res.fun
        if R>best: best=R; ba=res.x
    return best, ba, has_zero

# ---------------- EXACT certified lower bound from a rational-cosine config ------------------
def int_conv(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i,ai in enumerate(a):
        if ai==0: continue
        for j,bj in enumerate(b):
            out[i+j]+=ai*bj
    return out

def certified_ratio(cos_fracs, has_zero):
    """cos_fracs: list of Fraction c_j=cos(theta_j) for the +-pairs.
       Returns exact Fraction V = <g^2>/<g>^2  (certified M_n >= V)."""
    coeffs = [1]  # integer Fourier coeff vector of g (freq 0), builds up; global scale irrelevant
    if has_zero:
        coeffs = int_conv(coeffs, [-1, 2, -1])
    for c in cos_fracs:
        p, Q = c.numerator, c.denominator
        ker = [Q*Q, -4*p*Q, 2*Q*Q + 4*p*p, -4*p*Q, Q*Q]  # = 4Q^2 * [1/4,-c,1/2+c^2,-c,1/4]
        coeffs = int_conv(coeffs, ker)
    L = len(coeffs); mid = L//2
    A0 = coeffs[mid]                       # <g> (times global scale)
    S = sum(a*a for a in coeffs)           # <g^2> (times scale^2)
    return Fraction(S, A0*A0)              # scale cancels; exact rational

def snap(cos_vals, Q):
    """Snap cosines to rationals with denominator <= Q (nearest)."""
    out=[]
    for x in cos_vals:
        # nearest fraction with bounded denominator
        f = Fraction(x).limit_denominator(Q)
        out.append(f)
    return out

def certify_n(n, Qsnap=10**7, restarts=12, seed=3):
    M_num, halfangles, has_zero = optimize_sym(n, restarts=restarts, seed=seed)
    cosv = np.cos(halfangles)
    cf = snap(cosv, Qsnap)                 # high-denominator rational cosines: loss negligible
    V = certified_ratio(cf, has_zero)      # EXACT rational lower bound on M_n
    Vf = float(V)
    return dict(n=n, M_num=M_num, V=V, Vf=Vf, ratio_num=M_num/n,
                cert_ratio=(Vf-1)/n, has_zero=has_zero, Qsnap=Qsnap)

if __name__ == "__main__":
    import sys, time
    ns = [int(x) for x in sys.argv[1:]] or [50, 100, 240]
    print(f"{'n':>4} {'M_num':>12} {'M_num/n':>9} {'V(cert)':>14} {'V/n':>9} {'(V-1)/n':>9} {'>2/3?':>6}")
    for n in ns:
        t0=time.time()
        d = certify_n(n)
        gt = d['cert_ratio'] > 2/3
        print(f"{n:>4} {d['M_num']:>12.6f} {d['ratio_num']:>9.6f} "
              f"{d['Vf']:>14.6f} {d['Vf']/n:>9.6f} {d['cert_ratio']:>9.6f} {str(gt):>6}"
              f"   [{time.time()-t0:.1f}s]", flush=True)
        # print the exact rational's size (digits) for the record
        print(f"      certified: M_{n} >= {d['Vf']:.10f}  (exact rational, "
              f"num~{len(str(d['V'].numerator))}digits/den~{len(str(d['V'].denominator))}digits), "
              f"Q_snap={d['Qsnap']}", flush=True)
