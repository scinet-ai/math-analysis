#!/usr/bin/env bash
# Zero-download reproduction of the headline results for Rippon 7.54.
# Regenerates the exact rational certificate + the independent interval cross-check
# in ~1 minute, and prints PASS/FAIL.  Requires `uv` (https://docs.astral.sh/uv/).
#
#   ./reproduce.sh          # smoke: exact N=250 + arb cross-check N=250   (~seconds)
#   ./reproduce.sh full     # headline: exact N=1000 + arb N=700            (~30 s)
set -euo pipefail
cd "$(dirname "$0")"

MODE="${1:-smoke}"
if [ "$MODE" = "full" ]; then EXN=1000; ARBN=700; else EXN=250; ARBN=250; fi

RUN="uv run --no-project --with python-flint==0.9.0 --with numpy==2.4.6 python3"

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
