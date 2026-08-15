# Supporting proof packages for fourfold forrelation

This repository contains two complementary proof packages for the fourfold-forrelation results discussed in the associated manuscript. Both packages study fresh single-pass access, use branchwise hard photon-pass dose, and use \(N\) for the number of modes in each of the four masks (the manuscript uses \(M\)).

> [!IMPORTANT]
> Both derivations were produced with substantial large-language-model assistance. They passed the analytic, computational, and reproducibility checks recorded here, but have not received independent end-to-end verification by a subject-matter expert. The manuscript therefore states the results provisionally.

## Results and entry points

| package | result | rendered proof | verification |
|---|---|---|---|
| [`asymptotic_single_pass_floor/`](asymptotic_single_pass_floor/) | \(D \ge (2/15)N^{1/8}\) for every power of two \(N\ge 2^{30}\) in the full fresh single-pass model | [`main.pdf`](asymptotic_single_pass_floor/main.pdf) | [`verify_constants.py`](asymptotic_single_pass_floor/verify_constants.py) |
| [`finite_n4096_certificate/`](finite_n4096_certificate/) | at \(N=4096\), hard dose six is excluded for classically adaptive fresh single-pass batches that are block-diagonal in total signal photon number | [`forr4_n4096_advantage.pdf`](finite_n4096_certificate/output/pdf/forr4_n4096_advantage.pdf) | [`run_all.sh`](finite_n4096_certificate/code/run_all.sh) |

The asymptotic result assumes the exact promise \(\operatorname{Forr}_4\ge 1/4\) versus \(\operatorname{Forr}_4\le-1/4\). The finite certificate proves

$$
\operatorname{TV}(T_+,T_-)
\le 0.260969224792207924,
$$

and hence equal-prior Bayes error at least

$$
0.369515387603896038>\frac13.
$$

See each package README for its precise theorem, model, and internal file map:

- [Asymptotic package guide](asymptotic_single_pass_floor/README.md)
- [Finite certificate guide](finite_n4096_certificate/README.md)

## Quick verification

The portable focused checks require Python 3.13 and the packages in [`requirements-check.txt`](requirements-check.txt). They inspect committed artifacts and do not rewrite them.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-check.txt
make check
```

`make check` runs the asymptotic exact-rational checks and the finite package's focused paper, figure, and theorem checks. The same command runs in GitHub Actions.

The complete Round-2 through Round-4 replay remains available as an archival regression command:

```sh
cd finite_n4096_certificate/code
./run_all.sh --full
```

That replay is substantially longer and includes historical adversarial tests, so CI runs the theorem-critical focused checks. The recorded environment and results are documented in [`RUN_OF_RECORD.md`](finite_n4096_certificate/RUN_OF_RECORD.md), and the source inventory and artifact hashes are recorded in [`SOURCE_SNAPSHOT.md`](finite_n4096_certificate/code/SOURCE_SNAPSHOT.md).

## Building the PDFs

The document build additionally requires `latexmk` and a TeX distribution containing REVTeX 4.2, `amsmath`, `amssymb`, `booktabs`, `bm`, `graphicx`, `microtype`, `hyperref`, and `xcolor`.

```sh
make build
```

This rebuilds the asymptotic PDF, regenerates the finite package's vector figures, and refreshes the finite rendered proof at `finite_n4096_certificate/output/pdf/forr4_n4096_advantage.pdf`. Regenerable LaTeX auxiliary files and local Python environments are excluded by the repository `.gitignore`; the rendered proofs and theorem figures are intentionally committed.

## Scope

The asymptotic package treats the full fresh single-pass model. The finite \(N=4096\) certificate allows arbitrary states within each fixed total-signal-photon-number sector, idlers, repeated modes, collective measurements, and classical adaptation between fresh batches, but assumes that every fresh batch is block-diagonal in total signal photon number. It therefore does not cover coherent superpositions between different total signal-number sectors, including vacuum--nonvacuum coherence; extending the finite certificate to that unrestricted class is open. Both results concern branchwise hard photon-pass dose and exclude coherent quantum memory carried between fresh batches.

The packages do **not** establish a single-pass mean-dose lower bound, determine the unrestricted single-pass optimum at \(N=4096\), or prove that six is the optimal multipass dose.

The finite package includes its full frozen Round-2 through Round-4 source snapshot because the numerical certificate depends on that recorded dependency chain. Historical files in that snapshot preserve original environment paths and diagnostic language as provenance; the supported public entry point is `make check` at the repository root.

## Provenance, citation, and license

- The asymptotic package records its generation history in [`PROVENANCE.md`](asymptotic_single_pass_floor/PROVENANCE.md).
- The finite package records the accepted execution in [`RUN_OF_RECORD.md`](finite_n4096_certificate/RUN_OF_RECORD.md).
- Repository citation metadata is provided in [`CITATION.cff`](CITATION.cff).

No reuse license has been selected for this release. Publication on GitHub does not itself grant permission to copy, modify, or redistribute the contents; select an explicit license before inviting third-party reuse or contributions.
