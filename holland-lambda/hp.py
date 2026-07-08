"""
High-precision extremal M_n = max <g^2>/<g>^2, g = prod_{j=1}^n (1 - cos(theta - theta_j)).
Polish double-precision optima with mpmath Newton on the stationarity equations
(rotation gauge fixed: theta_1 = 0), then algebraically recognize M_n.
"""
import mpmath as mp
import numpy as np
from scipy.optimize import minimize

mp.mp.dps = 60

def factor_kernel(tj):
    # Fourier coeffs (freqs -1,0,1) of 1 - cos(theta - tj)
    return [mp.mpf(-1)/2 * mp.e**(1j*tj), mp.mpf(1), mp.mpf(-1)/2 * mp.e**(-1j*tj)]

def dfactor_kernel(tj):
    # d/d tj of the above
    return [mp.mpf(-1)/2 * 1j * mp.e**(1j*tj), mp.mpf(0), mp.mpf(1)/2 * 1j * mp.e**(-1j*tj)]

def convolve(a, b):
    n = len(a) + len(b) - 1
    out = [mp.mpc(0) for _ in range(n)]
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i+j] += ai*bj
    return out

def g_coeffs_and_grad(thetas):
    """Return (c, dc) where c = Fourier coeffs of g (freqs -n..n),
    dc[m] = d c / d theta_m (list over m)."""
    n = len(thetas)
    kernels = [factor_kernel(t) for t in thetas]
    # prefix/suffix products of kernels
    # base coeffs
    c = [mp.mpc(1)]
    for k in kernels:
        c = convolve(c, k)
    # gradient: replace kernel m by its derivative
    dc = []
    for m in range(n):
        acc = [mp.mpc(1)]
        for idx, k in enumerate(kernels):
            kk = dfactor_kernel(thetas[m]) if idx == m else k
            acc = convolve(acc, kk)
        dc.append(acc)
    return c, dc

def ratio_and_grad(free, gauge0=True):
    """free = theta_2..theta_n (theta_1=0). Return (R, grad wrt free)."""
    thetas = [mp.mpf(0)] + list(free)
    n = len(thetas)
    c, dc = g_coeffs_and_grad(thetas)
    N = sum(abs(ck)**2 for ck in c)          # <g^2> = sum |c_k|^2
    D = c[n].real                             # <g> = freq-0 coeff (real)
    R = N / D**2
    grad = []
    for m in range(1, n):  # derivative wrt free vars theta_2..theta_n = index 1..n-1
        dN = sum(2*(ck.conjugate()*dc[m][k]).real for k, ck in enumerate(c))
        dD = dc[m][n].real
        dR = (dN*D - 2*N*dD)/D**3
        grad.append(dR)
    return R, grad

def polish(free0):
    free0 = [mp.mpf(str(x)) for x in free0]
    def F(*args):
        _, g = ratio_and_grad(list(args))
        return g
    sol = mp.findroot(F, free0, tol=mp.mpf(10)**-50)
    free = [sol[i] for i in range(len(free0))] if hasattr(sol, '__len__') else [sol]
    R, _ = ratio_and_grad(free)
    return R, free

# ---- double precision seed ----
def npg_fourier(thetas):
    coeff = np.array([1.0+0j])
    for tj in thetas:
        coeff = np.convolve(coeff, np.array([-0.5*np.exp(1j*tj),1.0+0j,-0.5*np.exp(-1j*tj)]))
    return coeff
def npratio(free):
    thetas = np.concatenate([[0.0], free])
    c = npg_fourier(thetas); n=(len(c)-1)//2
    return np.sum(np.abs(c)**2)/c[n].real**2
def seed(n, restarts, seed_):
    rng = np.random.default_rng(seed_)
    best=-1; bf=None
    for _ in range(restarts):
        x0=rng.uniform(0,2*np.pi,size=n-1)
        r=minimize(lambda x:-npratio(x),x0,method='Nelder-Mead',
                   options={'xatol':1e-11,'fatol':1e-13,'maxiter':40000,'maxfev':40000})
        if -r.fun>best: best=-r.fun; bf=r.x
    return best, bf

def recognize(x, name):
    print(f"  {name} = {mp.nstr(x, 40)}")
    # try rational
    for prec in (25, 35):
        try:
            fr = mp.pslq([x, 1], maxcoeff=10**8, maxsteps=10**5)  # placeholder
        except Exception:
            pass
    # minimal polynomial search: pslq on [1, x, x^2, ..., x^d]
    for d in range(1, 7):
        vec = [x**k for k in range(d+1)]
        rel = mp.pslq(vec, maxcoeff=10**12, maxsteps=10**6)
        if rel:
            # rel are integer coeffs c_0 + c_1 x + ... = 0
            poly = " + ".join(f"{rel[k]}*x^{k}" for k in range(d+1) if rel[k]!=0)
            print(f"    minimal poly (deg {d}): {poly} = 0   coeffs={rel}")
            return rel
    print("    no low-degree algebraic relation found (maxcoeff 1e12)")
    return None

if __name__ == "__main__":
    import sys
    ns = [int(a) for a in sys.argv[1:]] or [2,3,4,5]
    for n in ns:
        rs = {2:200,3:300,4:500,5:800,6:1200,7:1500}.get(n,2000)
        val, bf = seed(n, rs, 100+n)
        R, free = polish(bf)
        E = 2*R - 1
        print(f"n={n}: M_n = {mp.nstr(R,42)}")
        print(f"       E_n=2M-1 = {mp.nstr(E,42)}   M/n={mp.nstr(R/n,12)}")
        recognize(R, f"M_{n}")
        recognize(E, f"E_{n}")
        # print angles
        angs = sorted([float(mp.nstr(mp.degrees(t),10)) for t in ([mp.mpf(0)]+list(free))])
        print(f"       angles(deg): {[round(a%360,4) for a in angs]}")
        sys.stdout.flush()
