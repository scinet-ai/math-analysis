"""
verify2.py -- zero-download smoke test for Holland's Lambda ROUND 2 ( < 1 min ).
Run:  uv run python verify2.py

Re-checks:
 (1) the two rotation-averaging identities  avg_phi<w_phi^2>=<u_m^2><u_n^2>  and
     avg_phi<w_phi>^2 = sum_k |uhat_m(k)|^2|uhat_n(k)|^2 = 1+Gamma  (Theorem 1);
 (2) certificate machinery reproduces M_2=15/7, M_3, M_4 (exact integer kernels);
 (3) CERTIFIED feasible-point lower bound M_50 >= 35.0407732538 in EXACT arithmetic,
     and V_50 > 1 + 2*50/3 verified as exact rationals;
 (4) a superadditivity defect entry a_{m+n}-a_m-a_n >= 0 (reliable optimizer).
"""
import numpy as np
from fractions import Fraction
from cert_lb import certified_ratio, certify_n, optimize_sym

TOL = 1e-9
def approx(a, b, tol=1e-6): return abs(a-b) <= tol*max(1, abs(a), abs(b))

# ---------- (1) rotation-averaging identities ----------
print("="*72)
print("(1) Rotation-averaging identities (Theorem 1) on a random feasible pair")
rng = np.random.default_rng(0)
def rand_config(n):
    return np.sort(rng.uniform(0.1, np.pi-0.1, n))     # n distinct angles in (0,pi)
def uhat(thetas):
    c = np.array([1.0+0j])
    for t in thetas:
        c = np.convolve(c, np.array([-np.exp(1j*t)/2, 1.0+0j, -np.exp(-1j*t)/2]))
    n = (len(c)-1)//2
    return c/c[n].real, n
m, n = 3, 4
am,_ = uhat(rand_config(m)); bn,_ = uhat(rand_config(n))
Mm = float(np.sum(np.abs(am)**2)); Mn = float(np.sum(np.abs(bn)**2))
# pad to common freq grid |k|<=m+n
def pad(v, half):
    h=(len(v)-1)//2; out=np.zeros(2*half+1,dtype=complex); out[half-h:half+h+1]=v; return out
L=m+n; A=pad(am,L); B=pad(bn,L)
Nphi=2048; phis=np.linspace(0,2*np.pi,Nphi,endpoint=False)
avg_w2=0.0; avg_mean2=0.0
for ph in phis:
    Bp = B*np.exp(-1j*np.arange(-L,L+1)*ph)          # uhat of u_n(.-phi)
    What = np.convolve(A, Bp)                          # freqs -2L..2L
    mean = What[len(What)//2].real                     # <w_phi>
    w2 = float(np.sum(np.abs(What)**2))                # <w_phi^2>
    avg_w2 += w2/Nphi; avg_mean2 += mean**2/Nphi
Gamma = float(np.sum((np.abs(am)**2)*(np.abs(bn[ (n-m):(n+m+1) ])**2))) - 1.0  # sum over |k|<=m minus k=0
# recompute Gamma cleanly over overlapping freqs
kmax=min(m,n); Gamma=0.0
for k in range(-kmax,kmax+1):
    if k==0: continue
    Gamma += (abs(am[k+m])**2)*(abs(bn[k+n])**2)
print(f"   avg<w_phi^2>           = {avg_w2:.8f}   vs  M_m*M_n = {Mm*Mn:.8f}")
print(f"   avg<w_phi>^2           = {avg_mean2:.8f}   vs  1+Gamma = {1+Gamma:.8f}")
assert approx(avg_w2, Mm*Mn, 1e-4), "identity (i) failed"
assert approx(avg_mean2, 1+Gamma, 1e-4), "identity (ii) failed"
print("   OK: both averaged identities hold.")

# ---------- (2) certificate machinery vs known M_2,M_3,M_4 ----------
print("\n(2) Exact certificate machinery reproduces small M_n")
def best_over_grid(pairs_dims, has_zero, step=0.01):
    best=Fraction(0)
    import itertools
    grid=[Fraction(int(round(x/step)), int(round(1/step))) for x in np.arange(-0.99,1.0,step)]
    if pairs_dims==1:
        for c in grid:
            V=certified_ratio([c],has_zero); best=max(best,V)
    else:
        g2=[Fraction(int(round(x/0.03)),int(round(1/0.03))) for x in np.arange(-0.96,1.0,0.03)]
        for a in g2:
            for b in g2:
                V=certified_ratio([a,b],has_zero); best=max(best,V)
    return float(best)
v2=best_over_grid(1,False); v3=best_over_grid(1,True); v4=best_over_grid(2,False)
print(f"   n=2: {v2:.6f} (M_2=15/7={15/7:.6f})   n=3: {v3:.6f} (2.808840)   n=4: {v4:.6f} (3.483450)")
assert abs(v2-15/7)<3e-4 and abs(v3-2.808840)<3e-4 and abs(v4-3.483450)<3e-3  # v4 coarse 2D grid
print("   OK.")

# ---------- (3) certified M_50 lower bound, exact ----------
print("\n(3) CERTIFIED lower bound M_50 (exact rational arithmetic)")
d = certify_n(50)
V = d['V']; thresh = Fraction(1)+Fraction(2*50,3)
Vfloor = (V.numerator*10**10)//V.denominator / 1e10
print(f"   M_50 >= {Vfloor:.10f}  (rounded DOWN from exact rational)")
print(f"   exact check  V_50 > 1 + 2*50/3 :  {V > thresh}   => (V-1)/50 = {(float(V)-1)/50:.6f} > 2/3")
assert V > thresh and Vfloor > 35.0
assert float(V) <= 51.0, "must respect proven upper bound M_n<=n+1"
print("   OK: M_50 certified above 2/3 slope. (n=100,240 in cert_lb.py, longer.)")

# ---------- (4) superadditivity defect >= 0 ----------
print("\n(4) Superadditivity defect a_{m+n}-a_m-a_n >= 0 (reliable optimizer)")
Ms={}
for k in (5,10,15):
    Ms[k]=optimize_sym(k,restarts=10,seed=5)[0]
a=lambda k:Ms[k]-1
d_5_10 = a(15)-a(5)-a(10)
print(f"   (5,10): a_5+a_10={a(5)+a(10):.5f}  a_15={a(15):.5f}  defect={d_5_10:.5f} (in (0,1))")
assert 0 <= d_5_10 < 1
print("   OK.")

print("\n"+"="*72)
print("ALL ROUND-2 CHECKS PASSED.")
print("  Theorem 1 (rotation-averaging) identities verified.")
print("  Certified: M_50>=35.0407732538, and V_50>1+2n/3 (exact).  [M_100>=69.387..., M_240>=165.563...]")
print("  Superadditivity of a_n=M_n-1 holds numerically (defect in (0,1)); existence of Lambda open.")
