"""
profile.py -- structure of the extremizer at large n:
  (1) coefficient-energy profile:  phi_N(s) ~ |uhat(k)|^2 at s=k/N
      (M_N = 1 + 2 sum_{k>=1} |uhat(k)|^2 ~ N * 2*int_0^1 phi)
  (2) zero-gap statistics: normalized gaps  g_i = N*(theta_{(i+1)}-theta_{(i)})/(2pi)
      (exact lattice: all gaps = 1; Poisson: exponential; what does extremal do?)
  (3) Fejer-smoothing saturation: for m=tN, ratio
        [sum_{|k|<=m} (1-|k|/(m+1))^2 |uhat_N(k)|^2] / M_m   (provably <= 1)
"""
import numpy as np
from round3 import optimize_M, uhat_from_angles

if __name__ == "__main__":
    N = 120
    print(f"optimizing N={N} ...", flush=True)
    MN, TH = optimize_M(N, restarts=12, seed=13)
    uh, n = uhat_from_angles(TH)
    p2 = np.abs(uh[n:])**2          # |uhat(k)|^2, k=0..N
    print(f"M_{N} = {MN:.6f}   M/N = {MN/N:.6f}   (M-1)/N = {(MN-1)/N:.6f}")
    print("\n(1) coefficient-energy profile |uhat(k)|^2 vs s=k/N (binned):")
    print(f"{'s-bin':>12}{'mean |uhat|^2':>15}")
    bins = np.linspace(0,1,21)
    ks = np.arange(1,N+1)
    for i in range(20):
        sel = (ks/N > bins[i]) & (ks/N <= bins[i+1])
        if sel.sum():
            print(f"[{bins[i]:.2f},{bins[i+1]:.2f}]{np.mean(p2[1:][sel]):>15.5f}")
    print(f"   |uhat(1)|^2={p2[1]:.5f}  |uhat(2)|^2={p2[2]:.5f}  ... |uhat(N)|^2={p2[N]:.6f}")
    print(f"   check: 1+2*sum = {1+2*p2[1:].sum():.6f} = M_N")

    print("\n(2) zero-gap statistics (normalized so lattice=1):")
    th = np.sort(np.mod(TH, 2*np.pi))
    gaps = np.diff(np.concatenate([th, [th[0]+2*np.pi]]))*N/(2*np.pi)
    print(f"   min={gaps.min():.4f}  max={gaps.max():.4f}  mean={gaps.mean():.4f}  std={gaps.std():.4f}")
    hist,edges = np.histogram(gaps, bins=np.linspace(0,3,16))
    for h,e0,e1 in zip(hist,edges[:-1],edges[1:]):
        print(f"   gap in [{e0:.1f},{e1:.1f}): {h}")

    print("\n(3) Fejer-smoothing saturation  (PROVED: ratio <= 1; how tight?)")
    print(f"{'m':>5}{'t=m/N':>7}{'smoothed sum':>14}{'M_m':>10}{'ratio':>8}")
    for m in [10,20,30,40,60,80,100]:
        w = np.zeros(N+1); k=np.arange(0,m+1)
        Mm,_ = optimize_M(m, restarts=8)
        sm = (1-k/(m+1))**2
        val = p2[0]*sm[0] + 2*np.sum(sm[1:]*p2[1:m+1])
        print(f"{m:>5}{m/N:>7.3f}{val:>14.5f}{Mm:>10.5f}{val/Mm:>8.4f}")
