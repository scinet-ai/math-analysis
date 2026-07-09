"""
mono.py -- test whether M_n/n is monotonically DECREASING via a leave-one-out
construction, which (with bounded-below) would PROVE existence of Lambda.

Target: M_{n+1} <= (1+1/n) M_n, i.e. M_n >= (n/(n+1)) M_{n+1}.

Leave-one-out: optimal u* in K_{n+1}, u*=|Q|^2, Q=prod(z-zeta_j), deg n+1.
Drop root j -> Q^{(j)} deg n, v_j = |Q^{(j)}|^2 / <|Q^{(j)}|^2> in K_n, <v_j^2> <= M_n.
Test:
  (i)   max_j <v_j^2>            >= (n/(n+1)) M_{n+1}   ?   (suffices)
  (ii)  avg_j <v_j^2>            >= (n/(n+1)) M_{n+1}   ?   (stronger, => (i))
  (iii) also report best <v_j^2>/M_n  (is a leave-one-out already >= M_n? can't exceed)

Also test the DIRECT monotonicity gap (n+1)M_n - n M_{n+1} and its trend.
And test a CONVOLUTION restriction:  v = u* (*) g, g in K_n chosen adaptively.
"""
import numpy as np
from round3 import optimize_M, M_from_angles, uhat_from_angles

def leave_one_out(angles_np1, Mn):
    """angles for optimal u* in K_{n+1}. Return stats vs M_n."""
    N1 = len(angles_np1)           # = n+1
    n = N1 - 1
    Mnp1 = M_from_angles(angles_np1)
    vals = []
    for j in range(N1):
        sub = np.delete(angles_np1, j)
        vals.append(M_from_angles(sub))
    vals = np.array(vals)
    target = (n/(n+1)) * Mnp1      # want some v_j >= this
    return dict(n=n, Mnp1=Mnp1, Mn=Mn, maxv=vals.max(), avgv=vals.mean(),
                target=target, max_ok=vals.max()>=target-1e-9,
                avg_ok=vals.mean()>=target-1e-9,
                maxv_over_Mn=vals.max()/Mn, gap=(n+1)*Mn - n*Mnp1)

if __name__ == "__main__":
    NS = list(range(1,25))
    M={}; TH={}
    print("optimizing...")
    for n in NS: M[n],TH[n]=optimize_M(n, restarts=12)

    print("\nLEAVE-ONE-OUT test:  does dropping a root from u*_{n+1} reach (n/(n+1)) M_{n+1}?")
    print(f"{'n':>3}{'M_{n+1}':>10}{'M_n':>10}{'target':>10}{'max v_j':>10}{'avg v_j':>10}{'max>=tgt':>9}{'avg>=tgt':>9}{'gap(n+1)Mn-nMn1':>16}")
    for n in range(1,24):
        d = leave_one_out(TH[n+1], M[n])
        print(f"{n:>3}{d['Mnp1']:>10.5f}{d['Mn']:>10.5f}{d['target']:>10.5f}{d['maxv']:>10.5f}"
              f"{d['avgv']:>10.5f}{str(d['max_ok']):>9}{str(d['avg_ok']):>9}{d['gap']:>16.6f}")

    print("\nInterpretation:")
    print(" max_ok True on all n  => monotonicity provable via 'best leave-one-out' (need the mechanism)")
    print(" avg_ok True on all n  => provable via the AVERAGE identity sum_j <v_j^2> >= n M_{n+1} (cleanest)")
