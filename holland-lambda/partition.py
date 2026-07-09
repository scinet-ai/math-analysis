"""
partition.py -- ROUTE 3 test: Fejer-Riesz zero-partition subadditivity.

From the degree-N extremal u* = |q|^2 (all zeros on the circle), partition the
N zeros into sets of size m and n=N-m: q = q_m * q_n, and form the normalized
feasible factors  v_m = |q_m|^2/<|q_m|^2> in K_m,  v_n = |q_n|^2/<|q_n|^2> in K_n.

Since <v_m^2> <= M_m and <v_n^2> <= M_n, IF for some partition
      <v_m^2> + <v_n^2>  >=  M_N - C                                   (*)
with provable C = O(1), then M_N <= M_m + M_n + C: approximate subadditivity
=> existence of Lambda by Fekete-with-error.

We measure, over partition strategies (balanced m=n=N/2 and m=1 leave-one-out):
   D(N) = M_N - max_partition [ <v_m^2> + <v_n^2> ]
Strategies: interleaved (alternate zeros), contiguous (arc split), random (many),
plus a local-swap hill climb from the best seed.
If D(N) grows linearly in N, route 3 is a dead end (documented).

Also the MULTIPLICATIVE diagnostic per partition: M_N vs <v_m^2><v_n^2>/(1+Gamma)
was round 2; here we record the additive deficit only.
"""
import numpy as np
from round3 import optimize_M, M_from_angles

def split_energy(angles, idx_m):
    """Given all N zero angles and index subset for the m-part, return
       <v_m^2> + <v_n^2> (normalized factor energies)."""
    mask = np.zeros(len(angles), bool); mask[list(idx_m)] = True
    Em = M_from_angles(angles[mask])
    En = M_from_angles(angles[~mask])
    return Em + En, Em, En

def best_partition(angles, m, tries=200, seed=0, climb=True):
    N = len(angles); rng = np.random.default_rng(seed)
    order = np.argsort(angles)
    cands = []
    # interleaved: every other zero (for balanced), every k-th
    if m == N//2:
        cands.append(set(order[::2][:m]))
    # contiguous arc
    for s in range(0, N, max(1,N//16)):
        cands.append(set(order[(s+np.arange(m)) % N]))
    # random
    for _ in range(tries):
        cands.append(set(rng.choice(N, m, replace=False)))
    best = -np.inf; bidx = None
    for c in cands:
        if len(c) != m: continue
        E,_,_ = split_energy(angles, c)
        if E > best: best = E; bidx = c
    if climb:
        improved = True
        cur = set(bidx)
        while improved:
            improved = False
            outside = [j for j in range(N) if j not in cur]
            for i in list(cur):
                for j in outside:
                    c2 = (cur - {i}) | {j}
                    E,_,_ = split_energy(angles, c2)
                    if E > best + 1e-12:
                        best = E; cur = c2; improved = True
                        break
                if improved: break
        bidx = cur
    return best, bidx

if __name__ == "__main__":
    print("ROUTE 3: zero-partition subadditivity deficit  D(N)=M_N - max_part[E_m+E_n]")
    print(f"{'N':>4}{'m':>4}{'M_N':>10}{'bestE':>10}{'D(N)':>9}{'M_m+M_n':>10}{'subadd_slack':>13}")
    rows=[]
    for N in [6,8,10,12,16,20,24]:
        MN, TH = optimize_M(N, restarts=12)
        m = N//2
        Mm,_ = optimize_M(m, restarts=10)
        bestE, idx = best_partition(TH, m, tries=300)
        D = MN - bestE
        rows.append((N,m,MN,bestE,D))
        print(f"{N:>4}{m:>4}{MN:>10.5f}{bestE:>10.5f}{D:>9.5f}{2*Mm:>10.5f}{2*Mm-MN:>13.5f}", flush=True)
    print("\nGrowth of D(N):")
    for i in range(1,len(rows)):
        N0,_,_,_,D0 = rows[i-1]; N1,_,_,_,D1 = rows[i]
        print(f"   D({N1})-D({N0}) = {D1-D0:+.5f}   slope ~ {(D1-D0)/(N1-N0):+.5f}/deg")
    print("\nIf D(N) ~ c*N with c>0: NO partition achieves (*) with C=O(1); route 3 dead end.")
