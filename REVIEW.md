# Standalone exposition review — September 16, 2026

A new round of five independent reviews covered narrative structure, the Fourier argument, Gaussian tools and concentration, finite reduction, and adversarial correctness. All five also reread the integrated `fourfold_forrelation.tex` draft.

Accepted changes:

- One progressive document: experiment, lower-bound strategy, parity example, Fourier collection, hard inputs, concentration, finite question, internal analytic appendix.
- Center the effect before the Fourier expansion; remove the separate zero-diagonal argument and the unused factorization-norm notation.
- Prove the Gaussian sign identity, Hermite contraction, Gaussian Poincare inequality, and variance tensorization inside the PDF.
- Replace references to a stored coefficient ledger with exact definitions of complete weighted kernel norms. Retain all support splits, proving a sufficient 210-state spectral criterion for arbitrary hard-dose-six parallel probes, including number coherence. This is broader than the earlier balanced conditional reduction, but no dose-six spectral estimate is asserted.
- Add a worked, proved finite dose-two obstruction: output total variation is below `1/100` at `N=4096`.
- Explain the physical-number/parity distinction and a false cubic entry bound within the PDF. Keep the larger historical audit separate from the reading path.
- Clarify that jointly prepared batches are included; a predetermined schedule alone does not make an adaptive protocol parallel. Derive adaptive polynomial degree from unnormalized amplitudes while leaving the required norm estimate open.

The pre-incorporation checkpoint was committed and pushed as [`2687e90`](https://github.com/safavilinqs/fourfold-forrelation-proofs/commit/2687e90f3e6540f17f9283e29298085513eacf27). The reviewers found no remaining mathematical gap in the established results. The finite dose-six spectral estimate and adaptive Fourier estimate remain unresolved.

---

# Independent review and adopted changes

Pre-review checkpoint: [`2b7a3cb`](https://github.com/safavilinqs/fourfold-forrelation-proofs/commit/2b7a3cb72e4f54b4989d66e213d4436333ad4554), already committed and pushed before these changes. The reviews used six independent agent assignments. Their conclusions were assessed against the mathematics and incorporated as follows.

| Review | Principal finding | Adopted change |
|---|---|---|
| Parallel Fourier argument | Record factorization is valid; nonconstant record symbols have zero diagonal | Center `E` at `I/2`, gaining a factor of two; clarify slots and idler vectors |
| Gaussian moments and concentration | No flaw found in the assigned steps; only even Fourier levels survive and direct promise-mass accounting is sharper | Explain the Hermite/translation steps and close with `1/12 < 1/5`, retaining `2/15` |
| Finite analytic reduction | Physical photon totals were incorrectly substituted for parity-support counts | Define the odd-occupation operator, correct the scope, and prove the residual-support lift explicitly |
| Inherited and local coefficients | Completed local-Walsh contraction is false; affected inherited and quintic entries remain in the ledger | Record the exact all-q norm counterexample and classify coefficient dependencies |
| Advanced coefficients and arithmetic | Four-cubic maximum is false; four residual entries lack a claimed normalization; outward arithmetic remains valid conditionally | Record an exact q=8 regression and separate the matrix certificate from a physical theorem |
| Presentation and simplicity | Conditional status and parity definitions must precede the finite argument; history obscures the main proof | Reorganize the finite note, show the actual input-independent interferometer preparation, and move detailed failures to the audit |

The proposed larger asymptotic constant `1/6` is valid with tighter arithmetic, but was not adopted: retaining `2/15` yields a cleaner ending and preserves the original theorem statement. The centered Fourier lemma is strengthened independently of that choice.

The reviews do not repair the adaptive extensions. A simple replacement for the false four-cubic maximum makes the affected incidence coefficients too large to recover the recorded certificate. Rather than reinstate an unsupported conclusion, the finite result is now explicitly conditional on its complete-kernel inequalities. A new joint bound or a different certificate is needed.

Coverage is recorded family by family in [COEFFICIENTS.md](finite_n4096_certificate/code/COEFFICIENTS.md). Sound structural reductions are distinguished from exhaustive revalidation of inherited scalar enumerations. In particular, 438 inherited high-sector entries and the inherited lower-degree package were not independently re-proved in this review. These reviews are substantive mathematical checks, not a formal proof certification.

The new exact regressions independently average signed-permutation moments and verify the failed premises. They are intended to prevent those premises being silently reinstated; they do not disprove the desired final lower bounds.
