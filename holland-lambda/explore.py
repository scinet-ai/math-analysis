"""
Holland Lambda_n exploration.

Reduction (see writeup):
  p in P_n  <=>  u(theta)=Re p(e^{i theta}) >= 0, trig poly deg <= n, mean(u)=1.
  Energy    E_n = max sum_{nu} |a_nu|^2 = max (1/2pi) int |p|^2  = 2 M_n - 1.
  GM value  M_n = max (1/2pi) int (Re p)^2 = max <u^2>  over u in K_n.
Extreme points of K_n (nonneg trig poly deg<=n, mean 1) are
  u_theta ∝ g(theta) = prod_{j=1}^n (1 - cos(theta - theta_j)),  n double zeros on circle.
So  M_n = max over (theta_1..theta_n) of  <g^2> / <g>^2 .

We compute <g^2>/<g>^2 exactly from Fourier coeffs of g (convolution of n
length-3 kernels [-1/2 e^{i tj}, 1, -1/2 e^{-i tj}]).
"""
import numpy as np
from scipy.optimize import minimize

def g_fourier(thetas):
    """Return Fourier coeffs of g = prod (1 - cos(theta - theta_j)) as complex array,
    index 0..2n corresponds to frequencies -n..n."""
    coeff = np.array([1.0 + 0j])  # frequency 0 only
    for tj in thetas:
        kernel = np.array([-0.5*np.exp(1j*tj), 1.0+0j, -0.5*np.exp(-1j*tj)])
        coeff = np.convolve(coeff, kernel)
    return coeff  # length 2n+1, freqs -n..n

def ratio(thetas):
    """<g^2>/<g>^2 where <.> is mean over circle."""
    c = g_fourier(thetas)
    n = (len(c)-1)//2
    mean_g = c[n].real  # frequency-0 coeff
    mean_g2 = np.sum(np.abs(c)**2)  # Parseval: <g^2> = sum |c_k|^2
    return mean_g2 / mean_g**2

def neg_ratio_reduced(x, n):
    # fix theta_1 = 0 (rotation invariance); optimize theta_2..theta_n
    thetas = np.concatenate([[0.0], x])
    return -ratio(thetas)

def optimize_n(n, restarts=200, seed=0):
    if n == 1:
        return ratio([0.0]), np.array([0.0])
    rng = np.random.default_rng(seed)
    best = -np.inf; best_th = None
    for _ in range(restarts):
        x0 = rng.uniform(0, 2*np.pi, size=n-1)
        res = minimize(neg_ratio_reduced, x0, args=(n,), method='Nelder-Mead',
                       options={'xatol':1e-10,'fatol':1e-12,'maxiter':20000,'maxfev':20000})
        val = -res.fun
        if val > best:
            best = val; best_th = np.concatenate([[0.0], res.x])
    return best, best_th

if __name__ == "__main__":
    print(f"{'n':>2} {'M_n=<u^2>max':>14} {'E_n=2M-1':>12} {'M_n/n':>10} {'E_n/n':>10} {'n+1':>5}")
    for n in range(1, 9):
        M, th = optimize_n(n, restarts=120 if n<=5 else 300, seed=42+n)
        E = 2*M - 1
        th_sorted = np.sort(np.mod(th, 2*np.pi))
        print(f"{n:>2} {M:>14.9f} {E:>12.9f} {M/n:>10.6f} {E/n:>10.6f} {n+1:>5}")
        print(f"     extremal angles (deg): {np.round(np.degrees(th_sorted),3)}")
