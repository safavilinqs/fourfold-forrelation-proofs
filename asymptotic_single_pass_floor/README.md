# Fourfold passive-dose proof package, v3.1

This package contains a standalone provisional derivation of the following theorem:

\[
D \ge \frac{2}{15}N^{1/8}
\]

for every power of two \(N\ge 2^{30}\), for the full fresh passive model and exact promises \(\mathrm{Forr}_4\ge 1/4\) versus \(\mathrm{Forr}_4\le -1/4\).

The exponent \(1/8\) is an unconditional positive exponent. It does not meet the original stronger target \(\alpha>1/4\).

## Verification status

The derivation was produced with substantial large-language-model assistance. Its definitions, numerical ledger, and build have been internally reviewed, but it has not received independent end-to-end verification by a subject-matter expert. It must therefore be cited as provisional, as it is in the main paper. The new, load-bearing point most in need of independent review is the fresh-pivot inequality and especially the recursive fresh-cut contraction in Lemma 4.2. `verify_constants.py` checks arithmetic only; it does not verify that analytic argument.

Here \(N\) is the number of modes on each of the four masks. The proof temporarily writes \(M=4N\) for the total number of Boolean sign coordinates. Thus the symbol called \(M\) in the main paper corresponds to this package's \(N\), not to its temporary \(M=4N\).

The claim concerns branchwise hard photon-pass dose. It does not establish a mean-dose lower bound, does not cover protocols carrying quantum memory between batches, and assumes one direct-sum mask interaction per fresh batch with no coherent inter-mask transform inside that batch.

## Files

- `main.tex`: standalone LaTeX source.
- `main.pdf`: compiled proof.
- `verify_constants.py`: exact-rational checks for the final constant ledger and rational consequences of the stated classical bounds on \(e\) and \(\pi\).
- `PROVENANCE.md`: clean-room provenance and scope.

## Build

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
python3 verify_constants.py
```

The build requires Python 3 and a TeX installation providing `latexmk` and `pdflatex`. A successful check prints the exact margin, conditioning error, average-case upper bound, exact-promise lower bound, and their positive difference. LaTeX should finish with no unresolved references or overfull boxes.
