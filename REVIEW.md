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
