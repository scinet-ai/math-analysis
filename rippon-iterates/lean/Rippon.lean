/-
  Rippon 7.54 (Hayman–Lingham): machine-checked FINITE certificate.

  φ_t(z) = e^{tz} - 1;  f_n := φ_t^n(-1),  f_1 = e^{-t}-1,  f_{n+1} = e^{t f_n} - 1.
  c^{(n)}_k := [t^k] f_n(t) ∈ ℚ.

  We implement the exact rational coefficient recurrence on truncated power series
  (`Array Rat`, index = degree) and machine-check, by `native_decide`, that over a
  window all coefficients satisfy |c^{(n)}_k| ≤ 1 with equality only at the leading
  term (k=n, value -1).  Pure Lean 4 core — no Mathlib.

  TRUST NOTE: `native_decide` compiles the Bool computation to native code and trusts
  its result via the `ofReduceBool`/compiler axioms.  This certificate's trusted base is
  {Lean kernel, native compiler}, i.e. NOT pure-kernel.  Disclosed accordingly.
-/

abbrev PS := Array Rat   -- truncated series, PS[k] = coefficient of t^k, length N+1

/-- Truncated product of two series to degree N. -/
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

/-- e^{u} - 1 = ∑_{m≥1} u^m / m!, truncated to degree N (requires u[0] = 0). -/
def expMinus1 (u : PS) (N : Nat) : PS := Id.run do
  let mut out : PS := Array.replicate (N+1) (0 : Rat)
  let mut term : PS := (Array.replicate (N+1) (0 : Rat)).set! 0 (1 : Rat)  -- u^0
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

/-- The iterates f_1, …, f_{nmax} as an array (index n-1), truncated to degree N. -/
def iterates (N nmax : Nat) : Array PS := Id.run do
  let mut fs : Array PS := #[]
  let u1 : PS := (Array.replicate (N+1) (0 : Rat)).set! 1 (-1 : Rat)  -- -t
  let mut f : PS := expMinus1 u1 N                                    -- f_1
  for _n in [1:nmax+1] do
    fs := fs.push f
    -- u = t * f  (shift up by one degree)
    let mut u : PS := Array.replicate (N+1) (0 : Rat)
    for k in [0:N] do
      u := u.set! (k+1) (f[k]!)
    f := expMinus1 u N
  return fs

/-- |c| ≤ 1 for a rational c ⟺ |num| ≤ den. -/
@[inline] def leq1 (c : Rat) : Bool := c.num.natAbs ≤ c.den

/-- Certificate: over the window 1≤n≤nmax, all coefficients (deg ≤ N) satisfy |c|≤1. -/
def noViolation (N nmax : Nat) : Bool := Id.run do
  let fs := iterates N nmax
  let mut ok := true
  for n in [1:nmax+1] do
    let f := fs[n-1]!
    for k in [0:N+1] do
      if ¬ leq1 (f[k]!) then ok := false
  return ok

/-- Leading term c^{(n)}_n = -1 for all n in the window. -/
def leadingMinusOne (N nmax : Nat) : Bool := Id.run do
  let fs := iterates N nmax
  let mut ok := true
  for n in [1:nmax+1] do
    if (fs[n-1]!)[n]! ≠ (-1 : Rat) then ok := false
  return ok

/-- |c^{(n)}_k| = 1 occurs ONLY at the leading term k=n (window). -/
def equalityOnlyLeading (N nmax : Nat) : Bool := Id.run do
  let fs := iterates N nmax
  let mut ok := true
  for n in [1:nmax+1] do
    let f := fs[n-1]!
    for k in [0:N+1] do
      let c := f[k]!
      if c.num.natAbs = c.den ∧ k ≠ n then ok := false   -- |c|=1 off the diagonal
  return ok

-- Sanity evals (small window)
#eval iterates 6 3      -- inspect f_1,f_2,f_3 to degree 6
#eval noViolation 12 12
#eval leadingMinusOne 12 12
#eval equalityOnlyLeading 12 12

/-- MACHINE-CHECKED: no coefficient of modulus > 1 for 1 ≤ n,k ≤ 40. -/
theorem rippon_no_violation_40 : noViolation 40 40 = true := by native_decide

/-- MACHINE-CHECKED: every leading coefficient c^{(n)}_n = -1 for 1 ≤ n ≤ 40
    (refuting the (-1)^n guess; the value is the constant -1). -/
theorem rippon_leading_40 : leadingMinusOne 40 40 = true := by native_decide

/-- MACHINE-CHECKED: |c^{(n)}_k| = 1 holds only at the leading term k=n, for n,k ≤ 40. -/
theorem rippon_equality_only_leading_40 : equalityOnlyLeading 40 40 = true := by native_decide
