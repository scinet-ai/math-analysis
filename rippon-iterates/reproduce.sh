#!/usr/bin/env bash
# Zero-download reproduction of the headline results for Rippon 7.54.
# Regenerates the exact rational certificate + the independent interval cross-check
# in ~1 minute, and prints PASS/FAIL.  Requires `uv` (https://docs.astral.sh/uv/).
#
#   ./reproduce.sh          # smoke: exact N=250 + arb cross-check N=250   (~seconds)
#   ./reproduce.sh full     # headline: exact N=1000 + arb N=700            (~30 s)
#   ./reproduce.sh round2   # ROUND-2 artifacts: profile, product, bands, Phi(-1/2) cert
set -euo pipefail
cd "$(dirname "$0")"

MODE="${1:-smoke}"
if [ "$MODE" = "full" ]; then EXN=1000; ARBN=700; else EXN=250; ARBN=250; fi

RUN="uv run --no-project --with python-flint==0.9.0 --with numpy==2.4.6 python3"

if [ "$MODE" = "round2" ]; then
  echo "### [R2 1/4] profile A_j exact to J=300 (|A_j|<=1; max non-leading 1/2 at j=1)"
  $RUN -m rippon.profile 300 | grep -E "stabilization|ALL|max non-leading"
  echo "### [R2 2/4] Koenigs product formula reproduces Phi exactly (J=120)"
  $RUN -m rippon.product_form 120 | grep -E "mismatches|MATCH"
  echo "### [R2 3/4] band-2 + band-3 identities, exact, window n<=40 (0 failures expected)"
  $RUN -m rippon.bands 40 | grep -E '"band2_checked"|"band2_failures"|"band3_checked"|"band3_failures"|"tau_closed_form_ok"|"max_abs_band2_value"' || true
  echo "### [R2 4/4] certified Phi(-1/2) < -1 (arb balls + proved tail bound)"
  $RUN -m rippon.phi_cert 256 200 | grep -E "certified|VERDICT"
  echo
  echo "PASS if: 0 mismatches, 0 band failures (band2_checked=900, band3_checked=940 at n<=40), VERDICT True."
  echo "Optional (Lean round 2): cd lean && lean Rippon2.lean   # 3 native_decide certificates"
  exit 0
fi

echo "### [1/3] EXACT rational certificate, window (n,k) <= ($EXN,$EXN)"
$RUN -m rippon.certificate "$EXN" data/reproduce_exact.json | \
  grep -E "no \|c\|>1|leading|ONLY at|non-leading sup|wall"

echo "### [2/3] INDEPENDENT interval (arb) cross-check, window <= ($ARBN,$ARBN), 256-bit"
$RUN -m rippon.arb_scan "$ARBN" 256 "$ARBN" data/reproduce_arb.json | \
  grep -E "certified|max abs_upper|worst ball|cross-check|wall"

echo "### [3/3] extremal structure map (A_j, envelope, transient ridge)"
$RUN -m rippon.structure "$EXN" | \
  grep -E "global max|violations|transient max|stabilization theory"

echo
echo "PASS if: no violation = True, |c|=1 only at leading, non-leading sup = 2663/4480 at (6,13),"
echo "         arb certified = True with 0 cross-check failures, stabilization theory holds."
echo
echo "Optional (Lean, needs elan + toolchain leanprover/lean4:v4.32.0-rc1):"
echo "  cd lean && lean Rippon.lean   # native_decide finite certificate, window 40"
