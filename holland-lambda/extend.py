"""
extend.py -- the FROZEN LOW-FREQUENCY EXTENSION problem:

    X(m,n) := max { <v^2> : v in K_{m+n},  vhat(k) = uhat_m*(k) for all |k| <= m }

where u_m* is the degree-m extremal. Any "additive splicing" proof of
approximate superadditivity M_{m+n} >= M_m + a_n - C would need the low band
frozen to the m-extremal to cost only O(1):  X ~ M_m + ~Lambda*n.
If instead X - M_m stalls (bounded in n), the additive route is dead.

  - X(m,n) <= M_{m+n} trivially.
  - X(m,n) >= (3/2)M_m via u_m*(1+cos Lθ), 2m+1 <= L <= n [proved, T-family].

Parametrization: v = |Q|^2, Q complex poly deg m+n (v>=0 automatic), Q
normalized (r_0 = ||Q||^2 = 1 = mean). Penalty-ramped L-BFGS on
   f(Q) = sum|r_l|^2 - lam * sum_{k=1..m} |r_k - target_k|^2 ,
lam ramped 10 -> 1e6; feasibility checked at the end (|r_k - t_k| < 3e-5).
Reported X is the exact <v^2> of the FEASIBLE-PROJECTED candidate: we take the
final Q, and report sum|r_l|^2 with the low band REPLACED by the targets and
the residual counted only if violation tiny (else discarded).
"""
import numpy as np
from scipy.optimize import minimize
from round3 import optimize_M, uhat_from_angles, Q_from_angles

def autocorr(Q):
    N=len(Q)-1
    r=np.zeros(N+1,dtype=complex)
    for l in range(N+1):
        r[l]=np.sum(Q[l:]*np.conj(Q[:N+1-l]))
    return r

def solve_extension(m, n, restarts=6, seed=0):
    N=m+n
    Mm,thm=optimize_M(m,restarts=10)
    uh,_=uhat_from_angles(thm)
    target=uh[m:2*m+1].copy()          # k=0..m, target[0]=1
    MN,thN=optimize_M(N,restarts=10)
    rng=np.random.default_rng(seed)
    qm=Q_from_angles(thm)
    best=-np.inf
    for rr in range(restarts):
        if rr==0:
            Q0=Q_from_angles(thN)
        elif rr==1:
            L=min(2*m+1,n)
            qL=np.zeros(L+1,dtype=complex); qL[0]=1/np.sqrt(2); qL[L]=1/np.sqrt(2)
            Q0=np.convolve(qm,qL); Q0=np.pad(Q0,(0,N+1-len(Q0)))
        elif rr==2:
            # u_m* times near-flat high tail
            tail=rng.normal(size=n+1)*0.1; tail[0]=1
            Q0=np.convolve(qm,tail/np.linalg.norm(tail)); Q0=Q0[:N+1]
        else:
            Q0=rng.normal(size=N+1)+1j*rng.normal(size=N+1)
        Q0=Q0/np.linalg.norm(Q0)
        x=np.concatenate([Q0.real,Q0.imag])
        for lam in [10.,100.,1e3,1e4,1e5,1e6]:
            def negf(x):
                Q=x[:N+1]+1j*x[N+1:]
                nrm=np.linalg.norm(Q); Q=Q/nrm
                r=autocorr(Q)
                obj=np.abs(r[0])**2+2*np.sum(np.abs(r[1:])**2)
                pen=np.sum(np.abs(r[1:m+1]-target[1:])**2)
                return -(obj-lam*pen)
            res=minimize(negf,x,method='L-BFGS-B',
                         options={'maxiter':2000,'ftol':1e-16,'gtol':1e-12})
            x=res.x
        Q=x[:N+1]+1j*x[N+1:]; Q=Q/np.linalg.norm(Q)
        r=autocorr(Q)
        viol=np.max(np.abs(r[1:m+1]-target[1:]))
        if viol<3e-5:
            val=float(np.abs(r[0])**2+2*np.sum(np.abs(r[1:])**2))
            if val>best: best=val
    return dict(m=m,n=n,Mm=Mm,MN=MN,X=best)

if __name__=="__main__":
    print("FROZEN LOW-FREQ EXTENSION  X(m,n): freeze vhat|_{|k|<=m} = m-extremal data")
    print(f"{'m':>3}{'n':>3}{'M_m':>9}{'M_(m+n)':>10}{'X(m,n)':>10}{'X-M_m':>9}{'M_mn-X':>9}{'a_n':>8}")
    for (m,n) in [(2,3),(2,5),(3,5),(2,8),(3,8),(4,8),(2,12),(3,12),(2,16),(3,16)]:
        d=solve_extension(m,n)
        Mn,_=optimize_M(n,restarts=8)
        print(f"{m:>3}{n:>3}{d['Mm']:>9.4f}{d['MN']:>10.4f}{d['X']:>10.4f}"
              f"{d['X']-d['Mm']:>9.4f}{d['MN']-d['X']:>9.4f}{Mn-1:>8.4f}",flush=True)
    print("\nX-M_m ~ 0.66n => additive splicing viable (frozen band costs O(1));")
    print("X-M_m bounded => additive route dead.")
