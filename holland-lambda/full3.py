"""
FULL (non-symmetric) critical-point enumeration for n=3, to certify that the
reflection-symmetric extremizer is the GLOBAL maximizer.

Extreme points: u ∝ prod_{j=1}^3 (1-cos(theta-theta_j)); gauge theta_1=0.
Free: theta_2, theta_3.  Variables c2=cos th2, s2=sin th2, c3, s3 with c^2+s^2=1.
R = <u^2>/<u>^2 (rational in c2,s2,c3,s3).  Solve grad_{theta}=0 exactly (Groebner),
list all real critical values, confirm max = M_3 (root of 208V^3-1224V^2+2268V-1323).
"""
import sympy as sp

th = sp.symbols('theta', real=True)
c2,s2,c3,s3 = sp.symbols('c2 s2 c3 s3', real=True)

def mean_trig(expr):
    expr = sp.expand(sp.expand_trig(expr))
    return sp.integrate(expr, (th, 0, 2*sp.pi))/(2*sp.pi)

# g = (1-cos th)(1-cos(th-th2))(1-cos(th-th3)), with cos(th-th2)=cos th*c2+sin th*s2
f1 = 1 - sp.cos(th)
f2 = 1 - (sp.cos(th)*c2 + sp.sin(th)*s2)
f3 = 1 - (sp.cos(th)*c3 + sp.sin(th)*s3)
g = sp.expand(f1*f2*f3)

D = sp.expand(mean_trig(g))            # <g>, poly in c2,s2,c3,s3
N = sp.expand(mean_trig(sp.expand(g**2)))  # <g^2>
print("built D,N")

# R = N/D^2 ; d/dtheta2 = -s2 dR/dc2 + c2 dR/ds2 ; numerator of stationarity:
# use  D^3 * dR/dtheta = (dN * D - 2 N dD) with the theta-derivative operators
def dtheta(expr, cc, ss):
    return -ss*sp.diff(expr, cc) + cc*sp.diff(expr, ss)

# stationarity numerators (times D^3>0):
St2 = sp.expand(dtheta(N,c2,s2)*D - 2*N*dtheta(D,c2,s2))
St3 = sp.expand(dtheta(N,c3,s3)*D - 2*N*dtheta(D,c3,s3))
print("built stationarity")

circ2 = c2**2+s2**2-1
circ3 = c3**2+s3**2-1

V = sp.symbols('V')
E0 = sp.expand(V*D**2 - N)

# Eliminate c2,s2,c3,s3 -> minimal poly in V
G = sp.groebner([E0, St2, St3, circ2, circ3], c2,s2,c3,s3,V, order='lex')
print("Groebner done; V-only elements:")
for gexpr in G.exprs:
    if gexpr.free_symbols <= {V}:
        print("   ", sp.factor(gexpr))
