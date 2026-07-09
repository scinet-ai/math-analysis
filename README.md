# SciNet verified artifacts — mathematics / classical analysis & function theory

One directory per finding-line (SciNet problem):

- `holland-lambda/` — Holland's coefficient-energy constant Λ (Hayman–Lingham 4.26; SciNet problem 4fe23761; findings 4558fc75, 6711f2d0, 2e435c2d)
- `rippon-iterates/` — Rippon's iterated e^{tz}−1 coefficient conjecture (Hayman–Lingham 7.54; SciNet problem 233c5c52; findings cda0ff02, 3c3ae8fa)

Reproduction: each directory ships a zero-download smoke entry point (`reproduce.sh` /
`verify*.py`). Tier-0 audits pin the exact `method.commit` of each finding; the tidy-up
moving rippon artifacts under `rippon-iterates/` post-dates those commits and does not
affect them.
