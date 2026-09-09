# Finite parallel-probe certificate at N=4096

Read [the proof PDF](output/pdf/forr4_n4096_advantage.pdf), built from the consolidated [main.tex](main.tex).

Three single photons, each traversing two masks, solve the promise at hard dose six with error `81/256`. The recorded lower-bound certificate excludes a **single parallel probe**, block diagonal in total signal photon number, at the same hard dose. Its output TV is at most `0.260969224792207925`, giving equal-prior error at least `0.369515387603896037`.

The previous extension to outcome-dependent fresh probes has an invalid norm inference. It is now an open proof obligation, documented in [AUDIT.md](../AUDIT.md). The normalized strategy induction alone does not close it. Cross-number coherence, mean dose, quantum memory and multipass optimality remain outside the result.

The short proof states the exact interface to the 210-state matrix. [COEFFICIENTS.md](code/COEFFICIENTS.md) maps the 888 balanced high-sector entries to their supporting derivations. Numerical replay checks the recorded coefficient formulas and outward ledger; it is not an independent proof of all analytic coefficient families.

Run `make check` and `make build` at the repository root. The frozen source snapshot preserves historical claims and tests; only the checks selected by the current wrapper support the current result.
