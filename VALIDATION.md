# Revision validation

Review baseline: `95c4c006d7566c7a01f68f0ee2ed465714dadf59`.
See [AUDIT.md](AUDIT.md) for the mathematical findings and limits.

## Length and presentation

| Main proof | Original PDF | Revised PDF | Approximate LaTeX source words |
|---|---:|---:|---:|
| Asymptotic | 17 pages | 5 pages | 6,832 → 2,123 (69% reduction) |
| Finite certificate | 14 pages | 5 pages | 8,329 → 2,091 (75% reduction) |

The PDFs use readable 11-point type and one-inch margins. Page reductions are factors of 3.4 and 2.8. Source-word counts use the same alphabetic-token rule before and after; they include LaTeX commands and are only an approximate measure of prose length. The finite baseline includes its former section and appendix files.

This is not a shorter proof of the former adaptive theorems: those extensions are no longer asserted as proved. The reduction also removes duplicate exposition and consolidates the supported arguments. The finite coefficient derivations remain external dependencies, mapped in [COEFFICIENTS.md](finite_n4096_certificate/code/COEFFICIENTS.md); their length is not included in the table. All 469 files in the frozen source snapshot retain their original Git blob hashes.

## Checks performed

- **Parallel Fourier record identity:** exact exhaustive comparison on 15,625 pairs of three-slot words over four nonzero labels and vacuum. This exercises repeated labels, cancellations, empty marked sides and unique last occurrences.
- **Invalid adaptive norm inference:** exact rational counterexamples at sizes 2, 3 and 16. These refute the stated abstract implication, not the physical adaptive conjecture.
- **Asymptotic constants:** exact rational verification of the concentration allowance, geometric-series bound and positive final gap.
- **Active protocol:** exact overlap/resource calculation and majority error `81/256` at dose six.
- **Finite replay:** the supported `make check` workflow passed, including all 888 accepted balanced high-sector coefficients, exclusion of 272 unbalanced incidences, the directed-arithmetic 210-state matrix, the committed positive Collatz candidate and the rational promise bound. [Completed replay](https://github.com/safavilinqs/fourfold-forrelation-proofs/actions/runs/34306207371) at `0b057e76ab39f4510adb3bd8f8930c394c415c3b`. Later edits expand the exact record check and refine wording/layout; they do not change the numerical backend.
- **Proof contract:** displayed decimal bounds are checked with exact fractions against the ledger; citations, coefficient counts and the current probe scope are checked as well.
- **Documents:** both PDFs rebuilt with LaTeX and inspected page by page. No overfull boxes or unresolved references. The comparison figure was regenerated with the corrected parallel-probe scope.
- **Repository integrity:** whitespace checks passed; the frozen snapshot is unchanged. No changes were made to the external sensing-paper repository.

The replay gives the outward upper bound

`TV <= 0.260969224792207925`,

and therefore the outward lower bound

`Bayes error >= 0.369515387603896037 > 1/3`.

The high-precision components are added before the displayed total is rounded. Adding the two separately displayed component upper bounds would be slightly looser.

## What these checks do not establish

The numerical replay does not independently prove every analytic coefficient bound. In particular, rounding a binary64 source value upward does not by itself certify the formula from which that value was computed. The retained family derivations remain necessary and have not all been independently rederived in this revision.

Neither the arithmetic nor the new structural regression supplies the missing adaptive contraction. The finite result also retains its restriction against coherence between total signal-number sectors. These are explicit limits of the revised claims, not issues that passing CI resolves.
