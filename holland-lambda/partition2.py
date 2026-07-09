"""
partition2.py -- larger-N zero-partition deficit D(N) = M_N - max_part[E_m + E_n]
to distinguish D(N) = O(log N)  (=> approximate-Fekete existence still viable)
from D(N) ~ c N (dead end).  Balanced m = N/2. Stronger search:
seeds = interleaved, arc, random(400); hill climb with first-improvement passes.
Also records: avg over random partitions, and the best partition's structure
(interleaving index pattern) for proof guidance.
"""
import numpy as np
from round3 import optimize_M, M_from_angles

def split_E(angles, mask):
    return M_from_angles(angles[mask]) + M_from_angles(angles[~mask])

def best_partition(angles, m, tries=400, seed=0, max_passes=6):
    N=len(angles); rng=np.random.default_rng(seed)
    order=np.argsort(np.mod(angles,2*np.pi))
    def mask_from(idx):
        mk=np.zeros(N,bool); mk[list(idx)]=True; return mk
    cands=[set(order[::2][:m])]
    for s in range(0,N,max(1,N//8)):
        cands.append(set(order[(s+np.arange(m))%N]))
    rand_vals=[]
    for _ in range(tries):
        c=set(rng.choice(N,m,replace=False))
        cands.append(c)
    best=-np.inf; bset=None
    for c in cands:
        if len(c)!=m: continue
        E=split_E(angles,mask_from(c))
        rand_vals.append(E)
        if E>best: best=E; bset=set(c)
    # hill climb (first improvement)
    for _ in range(max_passes):
        improved=False
        cur=list(bset); out=[j for j in range(N) if j not in bset]
        rng.shuffle(cur); rng.shuffle(out)
        for i in cur:
            for j in out:
                c2=(bset-{i})|{j}
                E=split_E(angles,mask_from(c2))
                if E>best+1e-10:
                    best=E; bset=c2; improved=True; break
            if improved: break
        if not improved: break
    # structure: sorted-order membership pattern
    patt=''.join('1' if j in bset else '0' for j in order)
    return best, np.mean(rand_vals), patt

if __name__=="__main__":
    import sys, time
    NS=[int(x) for x in sys.argv[1:]] or [16,24,32,40,48,56,64]
    print(f"{'N':>4}{'M_N':>10}{'bestE':>10}{'D(N)':>9}{'avgE(rand)':>11}{'D/logN':>9}{'D/N':>9}")
    rows=[]
    for N in NS:
        t0=time.time()
        MN,TH=optimize_M(N,restarts=12,seed=17)
        m=N//2
        bestE,avgE,patt=best_partition(TH,m)
        D=MN-bestE
        rows.append((N,D))
        print(f"{N:>4}{MN:>10.5f}{bestE:>10.5f}{D:>9.5f}{avgE:>11.5f}"
              f"{D/np.log(N):>9.5f}{D/N:>9.5f}   [{time.time()-t0:.0f}s]",flush=True)
        print(f"      pattern: {patt}",flush=True)
    print("\nfit: successive D increments vs 0.72*ln(N1/N0) [log model] vs c*(N1-N0) [linear]")
    for i in range(1,len(rows)):
        N0,D0=rows[i-1]; N1,D1=rows[i]
        print(f"  N {N0}->{N1}: dD={D1-D0:+.4f}   log-model {0.72*np.log(N1/N0):+.4f}")
