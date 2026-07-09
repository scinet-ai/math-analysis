"""
superadd.py -- numerical investigation of superadditivity of a_n = M_n - 1
and the rotation-averaging construction, for Holland's Lambda_n.

M_n = max over nonneg trig polys u (deg<=n, mean 1) of <u^2>.
Optimizer u = g/<g>, g = prod_j (1-cos(theta-theta_j)).

We test, for pairs (m,n):
  (A) exact superadditivity a_{m+n} >= a_m + a_n  (a_k=M_k-1), using optimized M.
  (B) Gamma = sum_{0<|k|<=min(m,n)} |uhat_m(k)|^2 |uhat_n(k)|^2, and whether the
      rigorous rotation-average bound M_mM_n/(1+Gamma) >= 1+A+B  (<=> Gamma <= AB/(1+A+B)).
  (C) whether the PRODUCT construction alone reaches it: max_phi <w_phi^2>/<w_phi>^2.
  (D) value at phi_0 (cross term c(phi)=0).
All Fourier coeffs computed EXACTLY from the angle config (kernel convolution).
"""
import numpy as np
from scipy.optimize import minimize

# ---------- exact Fourier coefficients of u = prod(1-cos(theta-theta_j)) / mean ----------
def kernel_coeffs(thetas):
    """Return complex Fourier coeffs c[k], k=-n..n, of g=prod(1-cos(theta-theta_j))."""
    c = np.array([1.0+0j])
    for t in thetas:
        ker = np.array([-np.exp(1j*t)/2, 1.0+0j, -np.exp(-1j*t)/2])  # freqs -1,0,+1
        c = np.convolve(c, ker)
    return c  # length 2n+1, index 0 -> freq -n

def u_hat(thetas):
    """Fourier coeffs uhat[k], k=-n..n of the NORMALIZED optimizer u=g/<g>."""
    c = kernel_coeffs(thetas)
    n = (len(c)-1)//2
    g0 = c[n].real  # mean of g = <g>
    return c/g0, n   # uhat[k] with index n -> freq 0

def M_from_thetas(thetas):
    uh, n = u_hat(thetas)
    return float(np.sum(np.abs(uh)**2))  # <u^2> = sum |uhat(k)|^2

# ---------- optimizer for M_n (grid objective, then exact eval) ----------
def RfromThetas_grid(thetas, Mgrid):
    th = np.linspace(0, 2*np.pi, Mgrid, endpoint=False)
    diff = th[None,:]-thetas[:,None]
    F = 1-np.cos(diff)
    logu = np.log(np.clip(F,1e-300,None)).sum(axis=0)
    u = np.exp(logu)
    U1=u.mean(); U2=(u*u).mean()
    R=U2/U1**2
    grad=np.zeros(len(thetas))
    for j in range(len(thetas)):
        duj=-np.sin(diff[j])*np.exp(logu-np.log(np.clip(F[j],1e-300,None)))
        dU1=duj.mean(); dU2=(2*u*duj).mean()
        grad[j]=(dU2*U1**2-U2*2*U1*dU1)/U1**4
    return R,grad

def optimize_M(n, restarts=8, seed=7):
    Mgrid=max(16*n,512)
    rng=np.random.default_rng(seed)
    best=-1; bth=None
    for r in range(restarts):
        if r==0:
            x0=np.linspace(0,2*np.pi,n,endpoint=False)+rng.normal(0,0.05,n)
        else:
            x0=np.sort(rng.uniform(0,2*np.pi,n))
        x0[0]=0.0
        def obj(x):
            xx=np.concatenate([[0.0],x])
            R,g=RfromThetas_grid(xx,Mgrid)
            return -R,-g[1:]
        res=minimize(obj,x0[1:],jac=True,method='L-BFGS-B',
                     options={'maxiter':3000,'ftol':1e-15,'gtol':1e-11})
        xx=np.concatenate([[0.0],res.x])
        R=M_from_thetas(xx)  # EXACT eval from angles
        if R>best: best=R; bth=xx
    return best,bth

# ---------- rotation construction quantities ----------
def rotation_analysis(thm, thn):
    """Given angle configs for u_m,u_n, compute Gamma, avg-bound, max_phi ratio, phi0 value."""
    am,_=u_hat(thm); bn,_=u_hat(thn)
    m=(len(am)-1)//2; n=(len(bn)-1)//2
    A=float(np.sum(np.abs(am)**2))-1.0
    B=float(np.sum(np.abs(bn)**2))-1.0
    # index arrays: freq k -> am[k+m], bn[k+n]
    kmax=min(m,n)
    Gamma=0.0
    for k in range(-kmax,kmax+1):
        if k==0: continue
        Gamma+=(abs(am[k+m])**2)*(abs(bn[k+n])**2)
    # rotation construction: w_phi = u_m(t) u_n(t-phi).
    # uhat_{u_n(.-phi)}(l) = bn[l+n] e^{-i l phi}.
    # <w_phi> = sum_k am[k+m]*conj(bn[k+n]) e^{ik phi}   (real)
    # <w_phi^2> = sum_l |What(l)|^2, What(l)=sum_k am[k+m] bn[l-k+n] e^{-i(l-k)phi}
    def ratio_at(phi):
        # mean
        mean=0j
        for k in range(-min(m,n),min(m,n)+1):
            mean+=am[k+m]*np.conj(bn[k+n])*np.exp(1j*k*phi)
        mean=mean.real
        # coeffs of w_phi, freq l in [-(m+n),(m+n)]
        L=m+n
        What=np.zeros(2*L+1,dtype=complex)
        for l in range(-L,L+1):
            s=0j
            for k in range(max(-m,l-n),min(m,l+n)+1):
                s+=am[k+m]*bn[(l-k)+n]*np.exp(-1j*(l-k)*phi)
            What[l+L]=s
        P=float(np.sum(np.abs(What)**2))
        return P/mean**2, mean, P
    # scan phi
    phis=np.linspace(0,2*np.pi,721,endpoint=False)
    vals=[]; means=[]
    for ph in phis:
        r,mn,P=ratio_at(ph); vals.append(r); means.append(mn)
    vals=np.array(vals); means=np.array(means)
    # refine max
    jmax=int(np.argmax(vals))
    from scipy.optimize import minimize_scalar
    lo=phis[(jmax-1)%len(phis)]; hi=phis[(jmax+1)%len(phis)]
    if hi<lo: hi+=2*np.pi
    rr=minimize_scalar(lambda p:-ratio_at(p)[0],bounds=(lo,hi),method='bounded')
    maxratio=-rr.fun
    # phi0 where mean crosses 1 (c(phi)=0)
    c=means-1.0
    phi0val=None
    for i in range(len(phis)):
        if c[i]==0 or c[i]*c[(i+1)%len(phis)]<0:
            # bisect
            a0=phis[i]; b0=phis[(i+1)%len(phis)]
            if b0<a0: b0+=2*np.pi
            for _ in range(60):
                mid=(a0+b0)/2
                cm=ratio_at(mid)[1]-1.0
                if (ratio_at(a0)[1]-1.0)*cm<=0: b0=mid
                else: a0=mid
            phi0=(a0+b0)/2
            v=ratio_at(phi0)[0]
            phi0val=v if phi0val is None else max(phi0val,v)
    return dict(A=A,B=B,Gamma=Gamma,thr=A*B/(1+A+B),
                avg_bound=(1+A)*(1+B)/(1+Gamma),
                maxratio=maxratio, phi0val=phi0val, target=1+A+B)

if __name__=="__main__":
    import sys
    # cache optimizers
    NS=sorted(set([1,2,3,4,5,6,7,8,9,10,12,15,16,20,24,30,40]))
    Mval={}; Th={}
    print("Optimizing M_n ...")
    for n in NS:
        M,th=optimize_M(n)
        Mval[n]=M; Th[n]=th
        print(f"  n={n:3d}  M_n={M:.8f}  M_n/n={M/n:.6f}  (M_n-1)/n={(M-1)/n:.6f}")
    print("\n(A) EXACT superadditivity a_{m+n} >= a_m + a_n  [a_k=M_k-1]")
    print(f"{'m':>3}{'n':>3}{'m+n':>5}{'a_m+a_n':>12}{'a_{m+n}':>12}{'defect':>12}{'ok':>4}")
    pairs=[(1,1),(1,2),(2,2),(1,3),(2,3),(3,3),(2,4),(4,4),(3,5),(5,5),(4,6),(5,7),(8,8),(10,10),(12,12),(15,15),(20,20)]
    for (m,n) in pairs:
        if m not in Mval or n not in Mval or (m+n) not in Mval: continue
        am=Mval[m]-1; an=Mval[n]-1; amn=Mval[m+n]-1
        defect=amn-(am+an)
        print(f"{m:>3}{n:>3}{m+n:>5}{am+an:>12.6f}{amn:>12.6f}{defect:>12.6f}{'Y' if defect>=-1e-6 else 'N':>4}")
    print("\n(B,C,D) rotation construction  [want Gamma<=thr, and maxratio/phi0val >= target]")
    print(f"{'m':>3}{'n':>3}{'Gamma':>10}{'thr=AB/(1+A+B)':>16}{'G<=thr':>7}{'avgbnd':>10}{'target':>10}{'maxratio':>10}{'phi0val':>10}")
    for (m,n) in pairs:
        if m not in Th or n not in Th: continue
        d=rotation_analysis(Th[m],Th[n])
        print(f"{m:>3}{n:>3}{d['Gamma']:>10.5f}{d['thr']:>16.5f}"
              f"{'Y' if d['Gamma']<=d['thr']+1e-9 else 'N':>7}"
              f"{d['avg_bound']:>10.5f}{d['target']:>10.5f}{d['maxratio']:>10.5f}"
              f"{(d['phi0val'] if d['phi0val'] else float('nan')):>10.5f}")
