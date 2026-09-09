# Fourfold forrelation: proofs and certificate

Start with the two short proof notes. Both count **hard photon-mask dose** and use `N` for modes per mask.

| Result | Supported probe class | Proof |
|---|---|---|
| `D >= (2/15) N^(1/8)`, powers of two `N >= 2^30` | One parallel quantum interrogation; arbitrary idlers, repeated modes and number coherence | [Asymptotic proof](asymptotic_single_pass_floor/main.pdf) · [source](asymptotic_single_pass_floor/main.tex) |
| At `N=4096`, dose six is insufficient | One parallel probe block diagonal in total signal photon number | [Finite certificate](finite_n4096_certificate/output/pdf/forr4_n4096_advantage.pdf) · [source](finite_n4096_certificate/main.tex) |
| Dose six suffices, with error `81/256` | Three single-photon two-pass flags | [Finite proof, Section 2](finite_n4096_certificate/main.tex) |

The proof audit found gaps in the previously claimed extensions to **outcome-dependent fresh batches**. Those extensions are now open proof obligations. See [AUDIT.md](AUDIT.md) for the exact missing steps and a counterexample to the finite argument's norm inference. The gaps do not disprove the desired adaptive bounds.

The presentations were developed with substantial language-model assistance. The finite coefficient derivations remain subject to independent review; a successful numerical replay is not a verification of every analytic lemma.

## Read the proofs

The asymptotic proof is: parallel Fourier bound → small Gaussian sign moments → concentration onto the promise. The finite proof is: folded interferometer → exact signed-permutation plant → occurrence matrix → outward spectral bound. [The coefficient guide](finite_n4096_certificate/code/COEFFICIENTS.md) identifies the finite proof's supporting family derivations.

For the finite parallel certificate,

`TV <= 0.260969224792207925`, hence `Bayes error >= 0.369515387603896037 > 1/3`.

The displayed bounds are rounded outward. The two packages do not establish a mean-dose lower bound, an adaptive optimum at `N=4096`, or optimality of multipass dose six.

## Verify and build

With Python and the packages in [requirements-check.txt](requirements-check.txt):

```sh
make check
make build
```

`make check` checks the parallel record identity, the counterexample guarding the adaptive norm interface, the asymptotic constants, manuscript bounds against the ledger, and the frozen finite arithmetic. It does not treat the old adaptive multiplier-one artifact as a theorem. `make build` additionally needs LaTeX, latexmk and the packages listed in the source preambles; CI builds both PDFs and publishes them as a workflow artifact.

The Round-2 through Round-4 source snapshot is retained unchanged for reproducibility. Its older adaptive and broader number-coherence claims are historical. The optional `finite_n4096_certificate/code/run_all.sh --full` is an archival regression, not evidence for those claims.
