# Baseline rerun record

Date: 2026-07-14

Frozen proof commit: 58dd3d7295960f3e53843695a9b73db57c06fbdd.

The archived round-one command initially resolved to Apple Python 3.9.6 and stopped before the mathematical checks because NumPy was unavailable. This is an environment failure, not a theorem failure.

Rerunning with the recorded Conda Python 3.13.13 environment reproduced NumPy 2.4.6 and SymPy 1.14.0. All twelve frozen proof checks and all seven round-one audit checks passed. The suite again detected the old rank-four unsliced-majorant defect, checked 250 graph slices, enumerated 7,692 small layered graphs and 65,536 two-chain placements, and verified the $1/24$ exponent ledger.

Interpretation: the frozen regressions are reproducible once the interpreter is pinned. They do not cover the arbitrary-map collision in RT-001 or establish the global reverse invariant in RT-002.
