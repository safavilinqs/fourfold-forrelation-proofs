# Fourfold forrelation: proofs and finite certificate audit

| Status | Result | Note |
|---|---|---|
| Proved | Parallel hard dose `D >= (2/15) N^(1/8)` for powers of two `N >= 2^30`; arbitrary idlers, repeated modes and number coherence | [Asymptotic PDF](asymptotic_single_pass_floor/main.pdf) · [source](asymptotic_single_pass_floor/main.tex) |
| Proved | Three two-pass flags use hard dose six and have error `81/256` | [Finite PDF](finite_n4096_certificate/output/pdf/forr4_n4096_advantage.pdf) · [source](finite_n4096_certificate/main.tex) |
| Conditional | Proposed finite parallel obstruction at `N=4096`, assuming complete-kernel bounds and diagonality in **parity-support size** | Same finite note; the analytic coefficient hypotheses remain unresolved |
| Open | The intended broader finite and adaptive lower bounds | [Audit](AUDIT.md) |

Independent reviews found a photon-number/parity identification error and false intermediate coefficient bounds in the finite package. Its numerical replay succeeds, but **does not establish an unconditional finite separation**. The adaptive extensions also remain unproved. See [REVIEW.md](REVIEW.md) for the reviewers' findings and the changes adopted.

The asymptotic argument is a parallel Fourier bound, a correlated Gaussian sign law, and concentration onto the promise. Centering the acceptance effect halves the Fourier estimate; the proof closes with `1/12 < 1/5`.

The finite note proves the interferometer construction and the exact signed-permutation plant. It then gives a complete conditional reduction from parity supports to a 210-state matrix. [COEFFICIENTS.md](finite_n4096_certificate/code/COEFFICIENTS.md) records which coefficient arguments survive review and which require repair. Frozen source files retain historical verdicts; those verdicts are superseded by the current audit.

## Verify and build

```sh
make check
make build
```

Use Python with [requirements-check.txt](requirements-check.txt); building also requires LaTeX and latexmk. The checks verify structural identities, exact counterexamples, asymptotic constants, manuscript arithmetic and the recorded finite ledger. Passing them does not prove unresolved kernel inequalities. CI also builds both PDFs. [VALIDATION.md](VALIDATION.md) records the checks and presentation changes.

The notes were developed with substantial language-model assistance. Neither mean-dose bounds nor optimality of multipass dose six is established. Changes are confined to this repository; the external sensing manuscript is unchanged.
