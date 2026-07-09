/-
  Rippon 7.54, ROUND 2: machine-checked finite certificates for the new
  structural identities (finding extending cda0ff02).

  f_0 = -1, f_{n+1} = e^{t f_n} - 1;  c^{(n)}_k = [t^k] f_n;  A_j = c^{(j)}_{2j}
  (stable profile, Lemma B).  h(x) = (e^x - 1)/x = sum_{m>=0} x^m/(m+1)!.

  Certified here (windows; exact rational arithmetic; native_decide):
    (F)  factorization      f_n = -t^n * prod_{k<n} h(t f_k)
    (B2) band-2 theorem     c^{(n)}_{2n+i} = A_{n+i} + tau_i   (0 <= i <= n+1)
         with tau_i = -(1/2) * sum_{l<i} [t^l] Phi^2
    (C1) corollary          c^{(n)}_{2n+1} = A_{n+1} - 1/2
    (C2) corollary          c^{(n)}_{2n+2} = A_{n+2}

  TRUST NOTE: native_decide trusts the Lean native compiler (ofReduceBool);
  trusted base = {Lean kernel, native compiler}, NOT pure-kernel.  Disclosed.
  Pure Lean 4 core - no Mathlib.
-/

abbrev PS := Array Rat   -- truncated series, PS[k] = coefficient of t^k

def mulTrunc (a b : PS) (N : Nat) : PS := Id.run do
  let mut out : PS := Array.replicate (N+1) (0 : Rat)
  for i in [0:N+1] do
    let ai := a[i]!
    if ai ≠ 0 then
      for j in [0:N+1-i] do
        let bj := b[j]!
        if bj ≠ 0 then
          out := out.set! (i+j) (out[i+j]! + ai * bj)
  return out

/-- e^{u} - 1, truncated to degree N (requires u[0] = 0). -/
def expMinus1 (u : PS) (N : Nat) : PS := Id.run do
  let mut out : PS := Array.replicate (N+1) (0 : Rat)
  let mut term : PS := (Array.replicate (N+1) (0 : Rat)).set! 0 (1 : Rat)
  let mut fact : Nat := 1
  for m in [1:N+1] do
    fact := fact * m
    term := mulTrunc term u N
    let inv : Rat := (1 : Rat) / (fact : Rat)
    for k in [0:N+1] do
      let tk := term[k]!
      if tk ≠ 0 then
        out := out.set! k (out[k]! + tk * inv)
  return out

/-- h(u) = (e^u - 1)/u = sum_{m>=0} u^m/(m+1)!, truncated (requires u[0] = 0). -/
def hOf (u : PS) (N : Nat) : PS := Id.run do
  let mut out : PS := (Array.replicate (N+1) (0 : Rat)).set! 0 (1 : Rat)
  let mut term : PS := (Array.replicate (N+1) (0 : Rat)).set! 0 (1 : Rat)  -- u^m
  let mut fact : Nat := 1   -- (m+1)!
  for m in [1:N+1] do
    fact := fact * (m + 1)
    term := mulTrunc term u N
    let inv : Rat := (1 : Rat) / (fact : Rat)
    for k in [0:N+1] do
      let tk := term[k]!
      if tk ≠ 0 then
        out := out.set! k (out[k]! + tk * inv)
  return out

/-- t * f  (shift up one degree). -/
def shiftT (f : PS) (N : Nat) : PS := Id.run do
  let mut u : PS := Array.replicate (N+1) (0 : Rat)
  for k in [0:N] do
    u := u.set! (k+1) (f[k]!)
  return u

/-- Iterates f_1 .. f_nmax (index n-1), truncated to degree N. -/
def iterates (N nmax : Nat) : Array PS := Id.run do
  let mut fs : Array PS := #[]
  let u1 : PS := (Array.replicate (N+1) (0 : Rat)).set! 1 (-1 : Rat)  -- t * f_0 = -t
  let mut f : PS := expMinus1 u1 N                                    -- f_1
  for _n in [1:nmax+1] do
    fs := fs.push f
    f := expMinus1 (shiftT f N) N
  return fs

/-- (F): f_n = - t^n * prod_{k=0}^{n-1} h(t f_k), checked to degree N for n <= nmax. -/
def factorizationOK (N nmax : Nat) : Bool := Id.run do
  let fs := iterates N nmax
  -- running product P_n = prod_{k<n} h(t f_k); f_0 = -1 so t f_0 = -t
  let mut prodP : PS := hOf ((Array.replicate (N+1) (0 : Rat)).set! 1 (-1 : Rat)) N
  let mut ok := true
  for n in [1:nmax+1] do
    -- expected: f_n[k] = - prodP[k - n]  (coefficient of t^{k-n} in prodP), 0 for k<n
    let f := fs[n-1]!
    for k in [0:N+1] do
      let expected : Rat := if k < n then 0 else -(prodP[k-n]!)
      if f[k]! ≠ expected then ok := false
    -- extend product with h(t f_n)
    if n < nmax then
      prodP := mulTrunc prodP (hOf (shiftT f N) N) N
  return ok

/-- Stable profile A_0..A_J read off f_J: A_j = c^{(J)}_{J+j} (Lemma B, n=J>=j).
    Needs truncation N >= 2J.  Returns array of length J+1. -/
def profileA (J : Nat) : PS := Id.run do
  let N := 2*J
  let fs := iterates N J
  let fJ := fs[J-1]!
  let mut A : PS := Array.replicate (J+1) (0 : Rat)
  for j in [0:J+1] do
    A := A.set! j (fJ[J+j]!)
  return A

/-- tau_i = -(1/2) sum_{l<i} [t^l] Phi^2, from the profile array A (needs i <= |A|-1). -/
def tauOf (A : PS) (I : Nat) : PS := Id.run do
  -- B_l = sum_{a+b=l} A_a A_b
  let mut tau : PS := Array.replicate (I+1) (0 : Rat)
  let mut acc : Rat := 0
  for i in [1:I+1] do
    -- add B_{i-1}
    let l := i - 1
    let mut bl : Rat := 0
    for a in [0:l+1] do
      bl := bl + A[a]! * A[l-a]!
    acc := acc + bl
    tau := tau.set! i (-acc / 2)
  return tau

/-- (B2): c^{(n)}_{2n+i} = A_{n+i} + tau_i for all 1 <= n <= nmax, 0 <= i <= n+1.
    Uses J = nmax+2 profile terms and truncation N = 3*nmax+3. -/
def band2OK (nmax : Nat) : Bool := Id.run do
  let N := 3*nmax + 3
  let J := 2*nmax + 3            -- need A up to n+i <= 2*nmax+3... wait i<=n+1 so n+i <= 2n+1 <= 2*nmax+1
  let A := profileA J            -- profileA computes its own iterates to degree 2J
  let tau := tauOf A (nmax+2)
  let fs := iterates N nmax
  let mut ok := true
  for n in [1:nmax+1] do
    let f := fs[n-1]!
    for i in [0:n+2] do
      let k := 2*n + i
      if k ≤ N then
        if f[k]! ≠ A[n+i]! + tau[i]! then ok := false
  return ok

/-- (C1)+(C2): c^{(n)}_{2n+1} = A_{n+1} - 1/2  and  c^{(n)}_{2n+2} = A_{n+2}. -/
def corollariesOK (nmax : Nat) : Bool := Id.run do
  let N := 2*nmax + 2
  let J := nmax + 2
  let A := profileA J
  let fs := iterates N nmax
  let mut ok := true
  for n in [1:nmax+1] do
    let f := fs[n-1]!
    if 2*n+1 ≤ N then
      if f[2*n+1]! ≠ A[n+1]! - (1:Rat)/2 then ok := false
    if 2*n+2 ≤ N then
      if f[2*n+2]! ≠ A[n+2]! then ok := false
  return ok

-- Sanity evals (tiny windows)
#eval factorizationOK 10 4
#eval band2OK 4
#eval corollariesOK 4

/-- MACHINE-CHECKED (F): the exact factorization f_n = -t^n prod_{k<n} h(t f_k)
    holds coefficientwise for all 1 <= n <= 30, degrees <= 60. -/
theorem rippon_factorization_30 : factorizationOK 60 30 = true := by native_decide

/-- MACHINE-CHECKED (B2): the band-2 (transient-line) identity
    c^{(n)}_{2n+i} = A_{n+i} + tau_i, tau_i = -(1/2) sum_{l<i} [t^l] Phi^2,
    holds for all 1 <= n <= 20, 0 <= i <= n+1. -/
theorem rippon_band2_20 : band2OK 20 = true := by native_decide

/-- MACHINE-CHECKED (C1,C2): c^{(n)}_{2n+1} = A_{n+1} - 1/2 and
    c^{(n)}_{2n+2} = A_{n+2} for all 1 <= n <= 25. -/
theorem rippon_corollaries_25 : corollariesOK 25 = true := by native_decide
