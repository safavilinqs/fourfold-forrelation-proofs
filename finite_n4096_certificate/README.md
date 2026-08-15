# Finite-size active--passive advantage for four-fold forrelation at $N=4096$

This folder contains the manuscript, technical appendices, bibliography records, audited source snapshot, and verification commands for the finite-size theorem

$$
D_{\mathsf P}^{\rm hard}(4096)>6,
\qquad
D_{\mathsf A}^{\rm hard}(4096)\le6.
$$

The certified passive class consists of classically adaptive trees whose fresh signal--idler state at every node is block diagonal in total signal photon number. It permits arbitrary mixtures of number sectors, idlers, repeated modes, within-sector entanglement and coherence, collective measurements, and classical feed-forward, with no coherent quantum memory carried between batches. The cost is branchwise hard photon-pass dose. The active protocol uses three single-photon interferometric flags, two charged sign-mask traversals per flag, and majority decoding.

The broader claim for fresh passive probes with coherence between different signal-number sectors is open. The audit found 272 unbalanced high-sector split/state incidences, forming 136 occupation edges, that the original package incorrectly called irrelevant and evaluated with unproved routing placeholders. The corrected verifier excludes every unbalanced placeholder and fails if one enters the theorem ledger.

## Build

```bash
python3 figures/generate_figures.py
latexmk -pdf main.tex
```

The paper uses REVTeX 4.2, `amsmath`, `amssymb`, `booktabs`, `bm`, `graphicx`, `microtype`, `hyperref`, and `xcolor`. The figure generator uses NumPy and Matplotlib and writes vector PDF and SVG outputs.

## Verify

Focused paper and theorem checks:

```bash
cd code
./run_all.sh
```

Historical inherited Round-2 through Round-4 replay:

```bash
cd code
./run_all.sh --full
```

The full command is intentionally longer because it reruns the original structural, adversarial, and certificate regressions. It passed on 2026-07-19 but contains superseded broad-scope assertions, so it is retained for provenance and is not evidence for the corrected theorem. The supported theorem check is the focused command above, whose invariants reject all unbalanced placeholders.

## Folder map

| path | contents |
|---|---|
| `main.tex` | paper master file |
| `macros.tex` | notation, theorem environments, and claim-grade labels |
| `sections/` | abstract and main-text sections |
| `appendix/` | hard-pair details, ledger certificate, adaptive proof, notation, and verification record |
| `figures/` | deterministic generator, documentation, and committed vector figures |
| `output/pdf/forr4_n4096_advantage.pdf` | verified 14-page rendered paper |
| `references.bib` | machine-readable bibliography |
| `dev/PAPER_PLAN.md` | section plan and completion gates |
| `dev/STYLE_CONTRACT.md` | binding prose and citation conventions |
| `dev/BIBLIOGRAPHY_NOTES.md` | metadata and claim scope for every conventional bibliography citation |
| `dev/SOURCE_MAP.md` | mapping from paper claims to certified source files and artifacts |
| `code/SOURCE_SNAPSHOT.md` | source provenance, audit history, inventory, and theorem-artifact hashes |
| `code/source_snapshot/` | Round-2 and Round-3 history plus the audited Round-4 certificate source |
| `code/run_all.sh` | focused and full verification entry point |
| `code/check_paper_contract.py` | manuscript structure, citation, and numerical-constant checks |

## Claim grades

The manuscript uses three grades.

- `Theorem`: proved in the source snapshot and covered by the outward certificate.
- `Derived specification`: exact arithmetic derived from the theorem or active protocol, without a device-specific noise model.
- `Experimental evidence`: a dated comparison with published demonstrations, used only to assess implementation credibility.

The number-sector-incoherent mathematical separation is theorem-level. Its extension to passive probes with cross-number-sector coherence is open. The $0.904294855157$ active contrast threshold is a derived specification. The conclusion that the $N=4096$ implementation is not yet experimentally credible is an evidence-grade assessment and is not used in the proof.
