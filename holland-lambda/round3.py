"""
round3.py -- Holland's Lambda, round 3.  Numerical tests of candidate rigorous
lemmas and the structural obstruction to existence of Lambda = lim M_n/n.

M_n = max over nonneg trig polys u (deg<=n, mean 1) of <u^2>.  Optimizer
u = g/<g>, g = prod_j (1-cos(theta-theta_j)); equivalently u=|Q|^2, deg Q<=n,
||Q||_2=1, and <u^2> = sum_ell |r_ell|^2, r_ell = autocorrelation of Q.

Tests:
 (A) DILATION-TENSOR: w(theta)=u_m(theta) u_n(k theta), k>=2m+1  =>  <w^2>/<w>^2 = M_m M_n
     (exact spectral disjointness). Verify equality for k>=2m+1, and undershoot for k<2m+1.
 (B) SQRT-SUBADDITIVITY: sqrt(M_{m+n}) <= sqrt(M_m)+sqrt(M_n).
 (C) OBSTRUCTION IDENTITY: split optimizer Q(deg N=m+n) coeffs into A=(Q_0..Q_m),
     S=(Q_{m+1}..Q_N).  Then M_N = <|A|^4> + <|S|^4> + 4 kappa + 2<f h>,
     f=|A|^2+|S|^2, h=2Re(conj(A)S), kappa=<|A|^2|S|^2>.  Show clean part
     <|A|^4>+<|S|^4> <= max(M_m,M_n) while cross = M_N - clean is Theta(n).
 (D) MONOTONICITY: is M_n/n decreasing and (M_n-1)/n increasing?
 (E) EXPLICIT FAMILIES: best provable slope lim <u^2>/n over closed-form nonneg
     families (Fejer, normalized Fejer^2, arc-concentrated) vs classical 2/3.
"""
import numpy as np
from numpy.polynomial import polynomial as P
from scipy.optimize import minimize

# ---------------- optimizer for M_n (angles), exact Fourier eval -------------------
def uhat_from_angles(thetas):
    """Complex Fourier coeffs of NORMALIZED u=g/<g>, g=prod(1-cos(theta-theta_j)).
       Returns array uh (freq -n..n, index n->freq0) and n."""
    c = np.array([1.0+0j])
    for t in thetas:
        c = np.convolve(c, np.array([-np.exp(1j*t)/2, 1.0+0j, -np.exp(-1j*t)/2]))
    n = (len(c)-1)//2
    return c/c[n].real, n

def M_from_angles(thetas):
    uh,_ = uhat_from_angles(thetas)
    return float(np.sum(np.abs(uh)**2))

def optimize_M(n, restarts=10, seed=7):
    if n==0: return 1.0, np.array([])
    Mgrid=max(16*n,512); rng=np.random.default_rng(seed)
    best=-1; bth=None
    for r in range(restarts):
        if r==0:  x0=np.linspace(0,2*np.pi,n,endpoint=False)+rng.normal(0,0.05,n)
        elif r==1:x0=np.concatenate([[0.0],np.sort(rng.uniform(0.02,1.5,n-1))]) if n>1 else np.array([0.0])
        else:     x0=np.sort(rng.uniform(0,2*np.pi,n))
        x0[0]=0.0
        th=np.linspace(0,2*np.pi,Mgrid,endpoint=False)
        def obj(x):
            xx=np.concatenate([[0.0],x])
            diff=th[None,:]-xx[:,None]; Fr=1-np.cos(diff)
            logu=np.log(np.clip(Fr,1e-300,None)).sum(0); u=np.exp(logu)
            U1=u.mean();U2=(u*u).mean();R=U2/U1**2
            g=np.zeros(len(xx))
            for j in range(len(xx)):
                duj=-np.sin(diff[j])*np.exp(logu-np.log(np.clip(Fr[j],1e-300,None)))
                g[j]=((2*u*duj).mean()*U1**2-U2*2*U1*duj.mean())/U1**4
            return -R,-g[1:]
        res=minimize(obj,x0[1:] if n>1 else np.array([]),jac=(n>1),method='L-BFGS-B',
                     options={'maxiter':3000,'ftol':1e-15,'gtol':1e-11}) if n>1 else None
        xx=np.concatenate([[0.0],res.x]) if n>1 else np.array([0.0])
        R=M_from_angles(xx)
        if R>best: best=R; bth=xx
    return best,bth

# ---------------- (A) dilation-tensor ----------------
def dilation_ratio(thm, thn, k):
    """<w^2>/<w>^2 for w=u_m(theta) u_n(k theta), exact via Fourier convolution."""
    am,m = uhat_from_angles(thm)         # freq -m..m
    bn,n = uhat_from_angles(thn)         # freq -n..n
    # u_n(k theta): freq k*j has coeff bn[j+n]; embed on grid -kn..kn
    Bfull = np.zeros(2*k*n+1,dtype=complex)
    for j in range(-n,n+1): Bfull[k*j + k*n] = bn[j+n]
    # w = u_m * (dilated u_n): convolve am(freq -m..m) with Bfull(freq -kn..kn)
    W = np.convolve(am, Bfull)           # freq -(m+kn)..(m+kn)
    mid=(len(W)-1)//2
    mean=W[mid].real
    w2=float(np.sum(np.abs(W)**2))
    return w2/mean**2, mean

# ---------------- (C) obstruction decomposition ----------------
def Q_from_angles(thetas):
    """Spectral factor Q (deg n) with |Q|^2 = g/<g>, ||Q||_2=1.
       g=prod(1-cos(th-th_j)) = |prod(e^{ith}-e^{ith_j})|^2 * (1/?) ; use q(z)=prod(z-e^{i th_j})."""
    q=np.array([1.0+0j])
    for t in thetas:
        q=np.convolve(q,[1.0,-np.exp(1j*t)])   # (z - e^{i t})
    # |q(e^{ith})|^2 = prod|e^{ith}-e^{ith_j}|^2 = prod 2(1-cos(th-th_j)) = 2^n g
    # so u = |q|^2 / <|q|^2>, and Q = q/||q||_2 gives |Q|^2 = u (mean 1).
    q=q/np.sqrt(np.sum(np.abs(q)**2))
    return q  # coeffs Q_0..Q_n

def poly_L2(coeffsA, coeffsB=None):
    """<A conj(B)> = sum_j A_j conj(B_j) (coeff dot). If B None, ||A||^2."""
    if coeffsB is None: coeffsB=coeffsA
    L=max(len(coeffsA),len(coeffsB))
    a=np.zeros(L,dtype=complex); a[:len(coeffsA)]=coeffsA
    b=np.zeros(L,dtype=complex); b[:len(coeffsB)]=coeffsB
    return np.vdot(b,a)   # sum a_j conj(b_j)

def mean_absprod(coeffsA, coeffsB):
    """<|A|^2 |B|^2> for A,B polynomials given by coeff vectors (freq via |.|^2)."""
    # |A|^2 has Fourier coeff = autocorrelation of A; product in freq = convolution.
    def autocorr(c):
        n=len(c)-1
        r=np.zeros(2*n+1,dtype=complex)
        for l in range(-n,n+1):
            s=0j
            for j in range(len(c)):
                if 0<=j-l<len(c): s+=c[j]*np.conj(c[j-l])
            r[l+n]=s
        return r,n
    rA,nA=autocorr(coeffsA); rB,nB=autocorr(coeffsB)
    # <|A|^2|B|^2> = sum_l rA[l] conj(rB[l])  (Parseval; both real functions -> real)
    L=max(nA,nB)
    RA=np.zeros(2*L+1,dtype=complex); RA[L-nA:L+nA+1]=rA
    RB=np.zeros(2*L+1,dtype=complex); RB[L-nB:L+nB+1]=rB
    return np.real(np.vdot(RB,RA))

def mean_abs4(coeffs):
    return mean_absprod(coeffs,coeffs)

def obstruction(thetas_N, m):
    """Split optimizer Q(deg N) coeffs into A(0..m),S(m+1..N)."""
    Q=Q_from_angles(thetas_N); N=len(Q)-1; n=N-m
    A=Q[:m+1]; Shigh=Q[m+1:]            # S = z^{m+1} * Shigh(z)
    # <|A|^4>, <|S|^4>=<|Shigh|^4>, kappa=<|A|^2|Shigh|^2>
    A4=mean_abs4(A); S4=mean_abs4(Shigh); kappa=mean_absprod(A,Shigh)
    MN=float(np.sum(np.abs(Q)**2)**0)  # placeholder
    MN=mean_abs4(Q)                     # = <|Q|^4> = M_N (since ||Q||=1)
    alpha=np.real(poly_L2(A))
    clean=A4+S4
    cross=MN-clean
    return dict(N=N,m=m,n=n,MN=MN,A4=A4,S4=S4,kappa=kappa,clean=clean,cross=cross,
                alpha=alpha,fourkappa=4*kappa)

# ---------------- (E) explicit families ----------------
def fejer_slope(n):
    """Fejer kernel F_n: uhat(k)=1-|k|/(n+1). <F^2>=sum uhat^2."""
    k=np.arange(-n,n+1); uh=1-np.abs(k)/(n+1)
    return float(np.sum(uh**2))

def normfejer2_slope(n):
    """u = c * F_{n/2}^2 (Fejer squared, degree n), normalized mean 1. (Jackson-type)"""
    h=n//2
    k=np.arange(-h,h+1); fh=1-np.abs(k)/(h+1)   # Fejer coeffs deg h
    # F^2 coeffs = convolution; then normalize so mean(uhat0)=1
    c2=np.convolve(fh,fh)                          # deg 2h
    c2=c2/c2[len(c2)//2]                            # normalize constant term to 1 (mean 1)
    return float(np.sum(c2**2))

def arc_concentrated_slope(n, w):
    """angles equally spaced on arc [0,w]; compute M via exact Fourier."""
    th=np.linspace(0.0, w, n, endpoint=True) if n>1 else np.array([0.0])
    return M_from_angles(th)

if __name__=="__main__":
    import sys
    print("Optimizing M_n (n up to 24)...")
    NS=[1,2,3,4,5,6,7,8,9,10,11,12,14,16,18,20,22,24]
    M={}; TH={}
    for n in NS:
        M[n],TH[n]=optimize_M(n)
    for n in NS:
        print(f"  n={n:3d}  M={M[n]:.6f}  M/n={M[n]/n:.6f}  (M-1)/n={(M[n]-1)/n:.6f}")

    print("\n(A) DILATION-TENSOR  <w^2>/<w>^2 vs M_m*M_n  (want EQ for k>=2m+1)")
    print(f"{'m':>2}{'n':>3}{'k':>3}{'2m+1':>5}{'ratio':>12}{'M_m*M_n':>12}{'rel.err':>11}")
    for (m,n) in [(2,3),(3,3),(2,4),(3,5),(4,4)]:
        for k in [2, 2*m, 2*m+1, 2*m+3, 4*m+1]:
            r,mean=dilation_ratio(TH[m],TH[n],k)
            prod=M[m]*M[n]
            print(f"{m:>2}{n:>3}{k:>3}{2*m+1:>5}{r:>12.6f}{prod:>12.6f}{abs(r-prod)/prod:>11.2e}")

    print("\n(B) SQRT-SUBADDITIVITY  sqrt(M_{m+n}) <= sqrt(M_m)+sqrt(M_n)")
    print(f"{'m':>3}{'n':>3}{'sqrt(M_mn)':>12}{'sqrtMm+sqrtMn':>15}{'slack':>10}{'ok':>4}")
    for (m,n) in [(1,1),(2,2),(3,3),(2,4),(4,4),(5,5),(6,6),(3,9),(2,10)]:
        if m+n not in M: continue
        lhs=np.sqrt(M[m+n]); rhs=np.sqrt(M[m])+np.sqrt(M[n])
        print(f"{m:>3}{n:>3}{lhs:>12.6f}{rhs:>15.6f}{rhs-lhs:>10.6f}{'Y' if lhs<=rhs+1e-6 else 'N':>4}")

    print("\n(C) OBSTRUCTION DECOMPOSITION  M_N = <|A|^4>+<|S|^4> + cross  (split at m)")
    print(f"{'N':>3}{'m':>3}{'M_N':>10}{'clean=A4+S4':>12}{'max(Mm,Mn)':>12}{'cross':>10}{'4kappa':>10}{'cross/N':>9}")
    for (m,n) in [(2,2),(3,3),(4,4),(5,5),(6,6),(4,8),(6,6)]:
        N=m+n
        if N not in TH: continue
        d=obstruction(TH[N],m)
        mm=max(M.get(m,0),M.get(n,0))
        print(f"{N:>3}{m:>3}{d['MN']:>10.5f}{d['clean']:>12.5f}{mm:>12.5f}{d['cross']:>10.5f}{d['fourkappa']:>10.5f}{d['cross']/N:>9.5f}")

    print("\n(D) MONOTONICITY   M_n/n decreasing?   (M_n-1)/n increasing?")
    dec=inc=True
    for i in range(1,len(NS)):
        a,b=NS[i-1],NS[i]
        if M[b]/b > M[a]/a+1e-9: dec=False
        if (M[b]-1)/b < (M[a]-1)/a-1e-9: inc=False
    print(f"   M_n/n monotonically decreasing on tested n: {dec}")
    print(f"   (M_n-1)/n monotonically increasing on tested n: {inc}")

    print("\n(E) EXPLICIT FAMILIES  slope <u^2>/n (want > 2/3=0.6667 for a rigorous LB improvement)")
    print(f"{'n':>4}{'Fejer/n':>10}{'normFej2/n':>12}{'arc(best w)/n':>14}{'opt M/n':>10}")
    for n in [20,40,80,160]:
        fj=fejer_slope(n)/n
        nf2=normfejer2_slope(n)/n
        # scan arc width
        bestarc=0; bw=0
        for w in np.linspace(0.5,2*np.pi,60):
            v=arc_concentrated_slope(n,w)/n
            if v>bestarc: bestarc=v; bw=w
        optn = M[max(k for k in M if k<=n)] / max(k for k in M if k<=n) if n<=24 else float('nan')
        print(f"{n:>4}{fj:>10.5f}{nf2:>12.5f}{bestarc:>14.5f}(w={bw:.2f}){'':2}")
