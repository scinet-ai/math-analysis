"""Guess a P-recurrence (holonomic/linear-ODE) for the profile A_j.

Seek  sum_{i=0}^{d} P_i(j) A_{j+i} = 0  with deg P_i <= r.  Detection is done in
F_p (prime) for speed: build the linear system over the unknown coefficients
c_{i,e} (of j^e in P_i), find its null space mod p, and require the relation to
hold on ALL computed A_j (including rows not used to solve) as a validity test.
Reports the smallest (d,r) admitting a consistent relation, or that none is found
up to the search bounds (evidence Phi is NOT holonomic of that complexity).
"""
import sys, json
from fractions import Fraction

P = (1 << 61) - 1  # Mersenne prime

def modinv(a, p=P): return pow(a % p, p - 2, p)
def to_mod(fr, p=P): return (fr.numerator % p) * modinv(fr.denominator, p) % p

def nullspace_modp(rows, ncols, p=P):
    """rows: list of length-ncols lists (mod p). Return list of null vectors (basis)."""
    M = [r[:] for r in rows]
    nrows = len(M)
    pivcol = []
    r = 0
    for c in range(ncols):
        piv = None
        for i in range(r, nrows):
            if M[i][c] % p != 0:
                piv = i; break
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        inv = modinv(M[r][c], p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(nrows):
            if i != r and M[i][c] % p != 0:
                f = M[i][c]
                M[i] = [(M[i][j] - f * M[r][j]) % p for j in range(ncols)]
        pivcol.append(c); r += 1
        if r == nrows: break
    free = [c for c in range(ncols) if c not in pivcol]
    basis = []
    for fc in free:
        v = [0] * ncols; v[fc] = 1
        for ri, pc in enumerate(pivcol):
            v[pc] = (-M[ri][fc]) % p
        basis.append(v)
    return basis, pivcol

def try_recurrence(A, d, r, p=P):
    """A: list of Fractions. unknowns c[i][e], i=0..d, e=0..r. col index = i*(r+1)+e.
       Row for shift j: sum_i sum_e c[i][e]*j^e*A[j+i]. Need j+d <= len-1."""
    ncols = (d + 1) * (r + 1)
    Amod = [to_mod(a, p) for a in A]
    Jmax = len(A) - 1 - d
    rows = []
    for j in range(0, Jmax + 1):
        row = [0] * ncols
        for i in range(d + 1):
            aij = Amod[j + i]
            je = 1
            for e in range(r + 1):
                row[i * (r + 1) + e] = (je * aij) % p
                je = (je * j) % p
        rows.append(row)
    if len(rows) < ncols + 4:   # want overdetermination for a real check
        return None
    # use first ncols+? rows to find nullspace, validate on the rest
    solve_rows = rows[: ncols + 2]
    basis, _ = nullspace_modp(solve_rows, ncols, p)
    valid = []
    for v in basis:
        ok = all(sum(v[c] * rows[jj][c] for c in range(ncols)) % p == 0 for jj in range(len(rows)))
        if ok:
            valid.append(v)
    return valid, ncols, len(rows)

def main():
    f = sys.argv[1] if len(sys.argv) > 1 else "data/profile_J600.json"
    A = [Fraction(s) for s in json.load(open(f))['A']]
    print(f"loaded {len(A)} profile coefficients A_0..A_{len(A)-1}")
    found = False
    for total in range(2, 14):        # search by increasing complexity d+r
        for d in range(1, min(total, 8) + 1):
            r = total - d
            if r < 0: continue
            res = try_recurrence(A, d, r)
            if not res: continue
            valid, ncols, nrows = res
            if valid:
                print(f"*** candidate holonomic recurrence: order d={d}, poly deg r={r} "
                      f"({ncols} unknowns, validated on {nrows} rows); nullity={len(valid)}")
                found = True
        if found:
            break
    if not found:
        print("No P-recurrence found with order<=8 and poly degree<= (12-order).")
        print("=> Phi is (empirically) NOT holonomic of that complexity; A_j not P-recursive within search box.")

if __name__ == "__main__":
    main()
