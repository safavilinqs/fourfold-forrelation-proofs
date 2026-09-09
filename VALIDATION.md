# Revision validation

Original baseline: `95c4c006d7566c7a01f68f0ee2ed465714dadf59`. Independent-review checkpoint: [`2b7a3cb`](https://github.com/safavilinqs/fourfold-forrelation-proofs/commit/2b7a3cb72e4f54b4989d66e213d4436333ad4554), committed and pushed before incorporation.

See [AUDIT.md](AUDIT.md) for mathematical findings and [REVIEW.md](REVIEW.md) for adopted reviewer inputs. The finite lower bound is conditional on unresolved kernel inequalities. Successful tests do not change that status.

## Presentation

The first revision reduced the main PDFs from 17 and 14 pages to five pages each, using 11-point type and one-inch margins. The independent review keeps the notes concise while adding the missing parity definitions and residual-support proof. Approximate current LaTeX source-word counts are 2,128 asymptotic and 1,914 finite, compared with 6,832 and 8,329 originally (about 69% and 77% reductions). Counts use the same alphabetic-token rule and include LaTeX commands.

The finite note now proves a conditional implication and states that its coefficient hypotheses are unresolved. The former stronger statements are not established by this shortening. Supporting coefficient derivations remain in the frozen source snapshot; their length is not included in the comparison.

## Current checks

- The unique last-occurrence record partition passes an exact exhaustive check on 15,625 word pairs, including repeated labels and vacuum.
- The new finite regression evaluates signed-permutation moments independently, without importing the frozen coefficient routines. It reproduces the physical-number/parity mismatch, the q=8 four-cubic counterexample, the incorrect Walsh phase identity, and the completed cubic amplitude's Gram eigenspaces at q=4. The audit supplies the general-q analytic derivation.
- Exact rational arithmetic verifies the simplified asymptotic ending: promise failure below 1/10, parallel upper gap 1/12, direct required gap above 1/5. The unproved adaptive Fourier target would give upper gap 1/6.
- The active majority error is exactly 81/256 at hard dose six. Manuscript decimals are checked outward against the unchanged ledger; the contract requires conditional status and parity-sector scope.
- The comparison figure uses conditional labels and rounds displayed upper and lower bounds in the appropriate directions.
- Whitespace checks pass. All 469 frozen source files are preserved unchanged. No changes were made to the external sensing manuscript.

The previous [completed ledger replay](https://github.com/safavilinqs/fourfold-forrelation-proofs/actions/runs/34306540292) reconstructs the same numerical backend and positive Collatz candidate. Its old success labels certify reproduction of the stored calculation, not the now-disproved intermediate proof premises. The new wrapper explicitly says that analytic coefficient hypotheses remain unresolved.

The recorded matrix and conditioning arithmetic would imply `TV <= 0.260969224792207925` and `Bayes error >= 0.369515387603896037` under the current theorem's hypotheses. These are conditional bounds, not an independently established physical separation.

## Review coverage and limits

Six focused reviews checked the Fourier, Gaussian, finite analytic, inherited-coefficient, advanced-coefficient and presentation arguments. The [coefficient guide](finite_n4096_certificate/code/COEFFICIENTS.md) separates sound norm reductions from false or unsupported premises and from scalar enumerations not exhaustively re-proved. Numerical upward rounding cannot repair a false analytic formula. The broader finite and adaptive lower bounds remain open.
