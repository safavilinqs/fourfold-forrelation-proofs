# Provenance

The v3 package was regenerated from the previously produced v2 clean package and the user-supplied review `Pasted text(8).txt`.

No repository history, parent-directory proof notes, prior hidden scratch files, or public search for this problem was used. The only mathematical source material used in this revision was the v2 package, the supplied review, and calculations produced during this revision.

The supplied review was not accepted wholesale. Its mandatory corrections were checked against the source. Its proposed removal of the middle contraction was rejected after an explicit one-batch counterexample, and its proposed constant `c0=1/2` was not claimed because the review did not provide a uniform all-level estimate. The package instead proves the conservative improved constant `c0=2/15` using the existing uniform bounds and exact arithmetic.

## Version 3.1 repository audit (2026-08-15)

For repository release, every line of `main.tex`, `README.md`, `PROVENANCE.md`, and `verify_constants.py` was reviewed against the theorem interface stated in the main paper. The exact-arithmetic checker and a clean two-pass LaTeX build were rerun. This audit corrected the hard-dose supremum and integer-charge wording, the count of possible last-occurrence records, the repeated-coordinate wording, the zero-dose edge case in finite approximation, and the derivative convention used in the conditional variance calculation. It also expanded the constant checker and made the package's provisional, LLM-assisted status explicit.

This repository audit is not independent subject-matter verification. In particular, the fresh-pivot inequality and the recursive fresh-cut contraction remain the central claims requiring expert review. No computational check in this package proves those analytic statements.
