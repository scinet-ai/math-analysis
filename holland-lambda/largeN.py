"""
Large-n asymptotics of M_n/n via grid-sampled extreme-point optimization with
analytic gradient (L-BFGS).  u(theta)=prod_j (1-cos(theta-theta_j)); maximize
R=<u^2>/<u>^2.  Also a symmetric-ansatz cross-check for n=5.
"""
import numpy as np
from scipy.optimize import minimize

def RfromThetas_and_grad(thetas, Mgrid=None):
    n = len(thetas)
    if Mgrid is None:
        Mgrid = max(8*n, 256)
    th = np.linspace(0, 2*np.pi, Mgrid, endpoint=False)
    # factors F[j] = 1-cos(th-theta_j); u = prod_j F[j]
    diff = th[None,:] - thetas[:,None]            # n x M
    F = 1 - np.cos(diff)                          # n x M
    logu = np.log(np.clip(F, 1e-300, None)).sum(axis=0)
    u = np.exp(logu)
    mean = lambda x: x.mean()                     # quadrature ~ (1/2pi)int
    U1 = mean(u); U2 = mean(u*u)
    R = U2/U1**2
    # d u/d theta_j = u * [ -sin(th-theta_j)/F_j ]  (F_j = 1-cos)
    #   = -sin(th-theta_j) * prod_{k!=j} F_k
    grad = np.zeros(n)
    for j in range(n):
        duj = -np.sin(diff[j]) * np.exp(logu - np.log(np.clip(F[j],1e-300,None)))
        dU1 = mean(duj); dU2 = mean(2*u*duj)
        grad[j] = (dU2*U1**2 - U2*2*U1*dU1)/U1**4
    return R, grad

def optimize(n, restarts=8, seed=0, Mgrid=None):
    rng = np.random.default_rng(seed)
    best=-1; bth=None
    for r in range(restarts):
        if r==0:
            x0 = np.linspace(0,2*np.pi,n,endpoint=False) + rng.normal(0,0.05,n)
        else:
            x0 = np.sort(rng.uniform(0,2*np.pi,n))
        x0[0]=0.0
        def obj(x):
            xx=np.concatenate([[0.0],x])
            R,g=RfromThetas_and_grad(xx,Mgrid)
            return -R, -g[1:]
        res=minimize(obj,x0[1:],jac=True,method='L-BFGS-B',
                     options={'maxiter':2000,'ftol':1e-14,'gtol':1e-10})
        R=-res.fun
        if R>best: best=R; bth=np.concatenate([[0.0],res.x])
    return best,bth

def sym5_value():
    # symmetric ansatz {0, +-b1, +-b2}: cross-check global M_5
    from scipy.optimize import minimize as mz
    def R(bb):
        b1,b2=bb
        th=np.array([0.0,b1,-b1,b2,-b2])
        r,_=RfromThetas_and_grad(th,Mgrid=512)
        return -r
    best=1e9;ba=None
    rng=np.random.default_rng(1)
    for _ in range(60):
        x0=rng.uniform(0.3,3.0,2)
        res=mz(R,x0,method='Nelder-Mead',options={'xatol':1e-10,'fatol':1e-12})
        if res.fun<best: best=res.fun; ba=res.x
    return -best, ba

if __name__ == "__main__":
    import sys
    if sys.argv[1:] == ["sym5"]:
        v,a = sym5_value()
        print(f"n=5 symmetric ansatz {{0,+-b1,+-b2}}:  M_5 = {v:.12f}  (global numeric 4.162256583165)")
        print(f"   b1,b2 (deg) = {np.degrees(a)}")
        sys.exit()
    ns = [int(a) for a in sys.argv[1:]] or [10,20,40,60,80,100,150,200]
    print(f"{'n':>4} {'M_n':>16} {'M_n/n':>12} {'(n+1)/n':>9}")
    for n in ns:
        R,th = optimize(n, restarts=6, seed=7, Mgrid=max(16*n,512))
        print(f"{n:>4} {R:>16.8f} {R/n:>12.8f} {(n+1)/n:>9.5f}", flush=True)
