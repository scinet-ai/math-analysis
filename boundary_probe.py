"""Boundary behavior of Phi (numeric, J=600 exact coefficients).

(1) L1 means M1(r) = (1/2pi) int |Phi(r e^{i theta})| dtheta vs r -> 1:
    if sup_r M1(r) <= 1 then |A_j| <= M1 < 1 would follow from the H^1 bound
    |A_j| <= ||Phi_r||_{L^1} -- is the H^1 route numerically viable?
(2) Radial limits at root-of-unity directions: Phi(r w) as r->1 for w=e^{2pi i p/q}:
    prediction (parabolic parameter): -> 0.  Supports dense boundary zeros +
    natural boundary picture (consistent with non-holonomicity).
(3) Radial behavior at an irrational angle (generic |t|=1): does |Phi| stay bounded?
"""
import json, math, cmath
from fractions import Fraction

d = json.load(open("data/profile_J600.json"))
A = [Fraction(s) for s in d['A']]
Af = [float(a) for a in A]
J = d['J']

def phi(t: complex) -> complex:
    s = 0j; tp = 1.0 + 0j
    for a in Af:
        s += a * tp
        tp *= t
    return s

print("=== (1) L1 and L2 means on |t|=r (trapezoid, 2048 pts; series J=600) ===")
print("    r      mean|Phi|   sqrt(mean|Phi|^2)   max|Phi|   tail_est")
for r in [0.5, 0.7, 0.8, 0.9, 0.95, 0.97, 0.985]:
    n = 2048
    s1 = s2 = mx = 0.0
    for i in range(n):
        v = abs(phi(r * cmath.exp(2j * math.pi * i / n)))
        s1 += v; s2 += v * v; mx = max(mx, v)
    # truncation error bound using |A_j|<=1/2 for j>=1 (verified to 600; heuristic beyond)
    tail = 0.5 * r ** (J + 1) / (1 - r)
    print(f"  {r:.3f}   {s1/n:.6f}     {math.sqrt(s2/n):.6f}      {mx:.4f}    {tail:.2e}")

print("\n=== (2) radial limits at root-of-unity directions e^(2 pi i p/q) ===")
for (p, q) in [(0,1),(1,2),(1,3),(1,4),(1,5),(2,5),(1,6)]:
    w = cmath.exp(2j * math.pi * p / q)
    vals = []
    for r in [0.9, 0.95, 0.98, 0.99]:
        vals.append(abs(phi(r * w)))
    print(f"  p/q={p}/{q}: |Phi(r w)| at r=0.9,0.95,0.98,0.99 -> " + ", ".join(f"{v:.4f}" for v in vals))

print("\n=== (3) radial behavior at 'generic' angles (irrational multiples of pi) ===")
for th in [1.0, 2.0, 0.5, 2.399963]:  # radians; 2.399963 ~ golden angle
    w = cmath.exp(1j * th)
    vals = []
    for r in [0.9, 0.95, 0.98, 0.99]:
        vals.append(abs(phi(r * w)))
    print(f"  theta={th:.4f}: |Phi| at r=0.9..0.99 -> " + ", ".join(f"{v:.4f}" for v in vals))

print("\n=== (4) partial sums sum_{j<=J} A_j (Phi(1^-) Abel check) ===")
ps = 0.0
marks = {}
for j in range(J + 1):
    ps += Af[j]
    if j in (50, 100, 200, 300, 400, 500, 600):
        marks[j] = ps
print("  partial sums:", {k: round(v, 6) for k, v in marks.items()})
