# Fourfold forrelation: a self-contained proof

Read **[the complete PDF](fourfold_forrelation.pdf)**. Its [single LaTeX source](fourfold_forrelation.tex) contains the entire exposition, including the Gaussian tools, finite reduction and diagnostic counterexamples. No other repository files are needed to follow its mathematical arguments.

The document starts with the two-pass interference protocol, explains the lower-bound strategy, develops the Fourier bound through a worked parity example, and then constructs and analyzes the hard inputs.

| Status | Result |
|---|---|
| Proved | Hard dose six and error at most `81/256` using three two-pass photons |
| Proved | Parallel hard dose `D >= (2/15) N^(1/8)` for powers of two `N >= 2^30`, allowing arbitrary idlers, repeated modes and number coherence |
| Proved | At `N=4096`, every parallel probe of dose at most two has output TV below `1/100` under the constructed promise laws |
| Proved reduction | At `N=4096`, an explicit 210-state matrix bounds distinguishability for every hard-dose-six parallel probe; all coefficients are defined in the PDF |
| Open | The spectral estimate needed for a finite dose-six lower bound, and the adaptive extensions |

The finite section no longer asks readers to trust an external coefficient ledger. It defines exact complete-kernel norms and proves the reduction to them. It does **not** claim the remaining spectral bound has been established. Fixed physical photon number alone does not permit the earlier balanced-parity restriction.

## Build and verify

```sh
make build
make check
```

Building the standalone PDF requires LaTeX and latexmk. Verification uses Python with [requirements-check.txt](requirements-check.txt). The structural and arithmetic checks supplement the proofs; successful numerical replay does not prove unresolved coefficient inequalities.

## Review history

[REVIEW.md](REVIEW.md) records the independent mathematical and presentation reviews. [VALIDATION.md](VALIDATION.md) records builds and checks. [AUDIT.md](AUDIT.md) and the [coefficient guide](finite_n4096_certificate/code/COEFFICIENTS.md) preserve the detailed audit of the earlier finite certificate. These are development records, not prerequisites for reading the PDF.

The earlier separate notes remain in `asymptotic_single_pass_floor/` and `finite_n4096_certificate/` for comparison and numerical reproduction; `make build-legacy` rebuilds them. Their presentations are superseded by the standalone document. The 469-file frozen source snapshot is unchanged.

The proofs were developed with substantial language-model assistance. Neither mean-dose bounds nor optimality of dose six is established. Only this repository is modified; the external sensing manuscript is unchanged.
