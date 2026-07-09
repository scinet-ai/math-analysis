"""
scaling.py -- how do the extremal angles concentrate as n grows?
Determines whether a FIXED limiting measure mu (Route 2 Gamma-limit) is well-posed,
or whether the concentration scale is coupled to n (making a fixed-mu limit fail).

For the optimizer u* of M_n (angles theta_1..theta_n, one pinned at 0), report:
  - angular spread: max gap-free arc / std of angles about their circular mean
  - the fraction of angles within arc [−w,w] for various w
  - the smallest angle scale (nearest neighbor near 0)
"""
import numpy as np
from round3 import optimize_M

def wrap(a):  # to (-pi,pi]
    return (a+np.pi)%(2*np.pi)-np.pi

if __name__=="__main__":
    print(f"{'n':>4}{'ang_std(rad)':>14}{'std*sqrt(n)':>12}{'std*n':>10}{'frac|θ|<0.5':>12}{'min|θ| (excl 0)':>16}{'minθ*n':>9}")
    for n in [8,12,16,24,32,48,64,96,128]:
        M,th = optimize_M(n, restarts=14, seed=11)
        w = np.sort(np.abs(wrap(th)))       # distances from 0
        std = np.std(wrap(th))
        frac = np.mean(np.abs(wrap(th))<0.5)
        nz = w[w>1e-6]
        mn = nz.min() if len(nz) else np.nan
        print(f"{n:>4}{std:>14.5f}{std*np.sqrt(n):>12.4f}{std*n:>10.3f}{frac:>12.3f}{mn:>16.5f}{mn*n:>9.3f}")
    print("\nIf ang_std ~ const  -> angles fill an O(1) arc: fixed-measure limit plausible (Route 2 well-posed).")
    print("If ang_std -> 0 (e.g. ~1/sqrt(n) or 1/n) -> concentration scale coupled to n: no fixed-mu Gamma-limit.")
