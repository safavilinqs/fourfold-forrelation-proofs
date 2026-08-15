# Style contract

## Prose

The manuscript uses sober, declarative scientific prose. Sentences may be long when their clauses carry distinct assumptions, mechanisms, or consequences, but each sentence must preserve a clear subject and a single logical direction. Headers are plain nominal descriptions such as `Definition`, `The model`, `Terms of the comparison`, `The hard pair`, `The one-batch certificate`, and `Limitations`.

Rhetorical dressing, management vocabulary, slogan-like compression, flourish vocabulary, dramatic fragments, and binary rhetorical contrast templates are excluded. A contrast is stated once, together with the physical or mathematical reason for the distinction.

Every formal object receives a physical reading near its definition. Jargon is introduced before use or omitted. Concrete values carry the interpretation: $N=4096$, $M=16{,}384$, active dose six, active error $81/256$, number-sector-incoherent passive transcript upper $0.260969224792208$, and active zero-bias contrast threshold $0.904294855157$.

## Claim scope

Each material statement is identified as theorem-level, derived specification, diagnostic history, or experimental evidence. Scope qualifiers remain attached to the claim, including the fresh-batch boundary, the absence of coherent inter-batch quantum memory, the branchwise hard-dose meter, and the distinction between an ideal resource specification and a demonstrated device.

The manuscript does not describe $N=4096$ as experimentally feasible. It states the measured requirements and the dated evidence assessment.

## Citations

Every external source is cited through a conventional `\cite{BibTeXKey}` command and appears in the numbered REVTeX bibliography. Every cited key is resolved in `dev/BIBLIOGRAPHY_NOTES.md`, which records the full citation, primary URL or DOI, and the exact sentence-level claim supported by the source. Forced inclusion through `\nocite{*}` and literal author--year placeholders are excluded.

## Figures

Scientific figures are generated from repository code and committed in vector PDF and SVG formats. Binary phase patterns are displayed on two-dimensional coordinate grids, with phases zero and $\pi$ identified explicitly; a square grid denotes the logical factorization $4096=64^2$ and is not described as a required physical layout. Figure captions distinguish exact theorem data, deterministic examples, schematic components, and evidence-grade comparisons.

## Mathematics

The theorem statement and numerical constants are copied from the certified artifacts. Rounded theorem values are labeled outward or lower as appropriate. Engineering calculations are separated from theorem dependencies. Proof summaries may compress algebra but may not change quantifiers or omit the physical access boundary.
