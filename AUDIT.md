# Proof audit and scope correction

Base: `95c4c006d7566c7a01f68f0ee2ed465714dadf59`.

The two presentations have been consolidated around their actual proof dependencies. This revision is a mathematical scope correction as well as an exposition change. It does not certify the former adaptive claims.

## Findings

| Argument | Assessment | Revision |
|---|---|---|
| Folded two-pass protocol | Overlap and majority error are exact. | One overlap equation and one binomial calculation. |
| Signed-permutation plant | Exact; the original telescoping explanation hid the useful identity. | Prove `H L = R` explicitly, then telescope three links. |
| Finite promise concentration | The proxy has a short direct proof. | Reveal blocks successively; the four squared coefficient budgets sum to `(1+beta²)(1+beta⁴)/N`. |
| Asymptotic one-link moments and concentration | The Gaussian calculations give the stated bounds. | Preserve the Hermite, conditional-variance and translation-symmetry arguments, removing duplicated explanations. |
| Parallel Fourier growth | A direct Schur factorization proves the bound, including repeated labels and vacuum coherence. | Replace the network-cut argument by the explicit record symbol in equation (4). |
| Asymptotic adaptive cut | The claimed common contraction is not established by the supplied construction. | State the parallel theorem and isolate the exact extra Fourier inequality needed for adaptation. |
| Finite normalized strategy induction | The squared-mass induction is valid. | Retain the invariant as an observation, not an adaptive theorem. |
| Finite adaptive interface | A trace-norm estimate is incorrectly used to bound an unrestricted entry sum. | Remove the adaptive claim and give an explicit counterexample to that inference. |
| Numerical certificate | Reproduction and outward rounding are distinct from proof of analytic coefficient bounds. | Retain the frozen coefficients, candidate and arithmetic; make that boundary explicit. |

## The finite adaptive error

Let `K=J_n`, the all-ones matrix, and let every scalar feature be `U_i=V_i=1/sqrt(n)`. Both feature families have squared mass one. The kernel has Schur factorization norm one. With `p_i=w_i=1/n`, the weighted matrix is `J_n/n`, whose trace norm is one. But the entry sum is

`sum_ij K_ij <V_j,U_i> = n`.

Thus normalized feature masses and a contractive Schur symbol do **not** imply the needed scalar bound. Pairing a matrix with the all-ones observable can cost its operator norm, which is `n`. The old proof controls the matrix in its equation (3.4), then passes to the entry sum without controlling that pairing. The old regression checks the induction and individual Schur/Perron facts, but never tests this missing implication.

This is a counterexample to the abstract inference, not a counterexample to the desired physical adaptive lower bound. A repair must use additional causal/measurement structure that the invariant has discarded.

## The asymptotic adaptive gap

The original fresh-cut lemma introduces input/cut spaces and asserts a contraction after a zig-zag rearrangement. It does not give compatible explicit maps that establish the rearranged operator norm. In particular, an instrument isometry on its original input-output split is not automatically contractive across a new tensor split. Outcome-dependent preparations must be accounted for in that same norm.

The new proof avoids this issue for a parallel probe. Its input and effect are fixed, and a record kernel is explicitly the product of a pulled-back small coefficient matrix, a residual-parity equality kernel and two one-sided filters. This proves the same Fourier bound for that class. No counterexample to the broader adaptive Fourier inequality is asserted; it is recorded as an open obligation rather than an established lemma.

## What the finite arithmetic verifies

The retained numerical backend reconstructs the accepted 888 balanced high-sector rows, rejects unbalanced placeholders, checks exact squared coefficients for the residual and dual-endpoint families, reconstructs the 210-state matrix with directed arithmetic, evaluates the committed positive Collatz candidate, and compares the full regenerated ledger with its committed artifact. It also checks the exact rational promise relaxation.

Rounding a binary64 value upward certifies an upper bound on that floating value, not automatically on an analytic expression. The family proofs remain necessary. This revision does not claim an independent rederivation of all 888 coefficient cases. The separate coefficient guide maps each accepted family to its source.

## Scope of the revised results

- Asymptotic: one parallel batch, arbitrary number-sector coherence, idlers and measurements; hard dose at least `(2/15)N^(1/8)` for powers of two `N >= 2^30`.
- Finite: one parallel batch at `N=4096`, block diagonal in total signal number, hard dose six excluded by the recorded coefficient certificate.
- Both also cover predetermined collections of batches that can be prepared jointly within the same hard cap. The finite input must remain block diagonal in total signal number.
- Outcome-dependent fresh probes require an additional argument in both packages. Mean-dose bounds and coherent quantum memory are not established. The active construction proves an upper bound of six, not optimality.

The former presentations remain available in Git history. Frozen source snapshots are historical: their assertions such as `CERTIFIED` for the adaptive multiplier have been superseded by this audit and are not supported theorem verdicts. The external sensing manuscript was not modified.
