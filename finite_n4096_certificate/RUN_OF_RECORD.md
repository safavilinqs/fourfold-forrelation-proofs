# Run of record

Status: **PASS for the corrected, number-sector-incoherent theorem.**

Audit date: 2026-08-15, America/Los_Angeles.

## Scope correction

The original 2026-07-19 run certified all 888 balanced high-sector entries but also routed 5,128 unbalanced profile-splits through placeholders that were not coefficient theorems. Those splits induce 272 directed split/state incidences, or 136 undirected edges, between different total-number sectors. They therefore cannot support the originally stated theorem for fresh probes with coherence between total signal-photon-number sectors.

The corrected theorem assumes that every fresh signal--idler batch is block diagonal in total signal photon number. It still permits mixtures across number sectors, arbitrary within-sector coherence and entanglement, arbitrary idlers, repeated modes, collective measurements, and unlimited classical adaptation. `AUDIT_CORRECTION.md` gives the complete claim boundary and supersedes the broad statements preserved in the historical Round-2, Round-3, and pre-audit Round-4 source snapshot.

## Source state and environment

The source snapshot began at repository commit `6aff0c9` (`Close Round 4 theorem and feasibility decision`) and preserves the complete Round-2, Round-3, and Round-4 history. The 2026-08-15 review corrected the active theorem path, artifact schema, verification entry points, manuscript, and outward-facing documentation without rewriting the historical record.

Corrected verification environment:

- macOS 26.5.1, build 25F80;
- Python 3.13.13;
- NumPy 2.4.6;
- SciPy 1.17.1;
- SymPy 1.14.0;
- Matplotlib 3.10.9;
- Latexmk 4.88.

## Corrected focused verification

Command:

```bash
cd code
./run_all.sh
```

Result: **PASS**. The corrected run refuses every unbalanced placeholder and checks that the theorem high-sector map is exactly the 888 balanced entries. It reported

- Collatz--Wielandt Perron upper bound `0.25874409638557722267923079152925562403311987598701005323368022563713944348400092`;
- promise-conditioning loss upper bound `0.00222512840663070225495016586456093734557722903264557427546649971895038635139684472693787993`;
- total-variation upper bound `0.26096922479220792493418095739381656137869710501965562750914672535608982983539778`;
- reserve lower bound `0.071364108541125408399152375939516771954636228313677705824186607977243503497935553`;
- equal-prior Bayes-error lower bound `0.36951538760389603753290952130309171931065144749017218624542663732195508508230111`;
- adaptive frontier multiplier one; and
- active dose six with majority error `81/256`.

The repository-level `make check` command reran this focused suite and the asymptotic exact-constant checker successfully.

## Historical full suite

The inherited command `./run_all.sh --full` passed on 2026-07-19, but that historical suite contains superseded broad-scope assertions and is not evidence for the corrected theorem. It is retained for provenance and adversarial-test history only. The corrected publication claim relies on the focused verification above, whose invariants exclude the unproved unbalanced routes.

## Document build and inspection

Command:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

Result: **PASS**. The corrected document contains 14 letter-size pages. The log has no LaTeX errors, undefined controls, undefined references, undefined citations, overfull boxes, or unprocessed floats. Remaining underfull-box notices are confined to long code-like labels.

The corrected `output/pdf/forr4_n4096_advantage.pdf` has SHA-256 hash `e849b7c1f9853fc592774b58c62ab313a7c1f8c0fe3aacea732504cc320442f6`. The committed figure PDFs have hashes `92835475834b71d40aee4149eee722325fa5cb6839e3f74ab3eb3aba546354a0` for the active protocol, `8f0c0823a3cd5b71d2362d4814b575d0f6a409e6f9a666acd71ad621eaf7078d` for the corrected passive certificate, and `916867757b253bb6800f4c429d80c9174d319962efe8bbe9c6592489403483d0` for the signed-permutation phase grids.
