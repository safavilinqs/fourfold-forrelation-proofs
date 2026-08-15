# Paper plan

Date: 2026-07-19

## Product

A self-contained RevTeX paper presenting the certified finite-size active--passive separation for four-fold forrelation at $N=4096$, together with the exact hard distribution, the complete numerical certificate, the adaptive factorization, the active protocol, reproducible vector figures, a physical resource assessment, and the complete verification source snapshot.

## Sources of truth

1. `code/source_snapshot/open_problem_forr4_passive_floor_consolidation_round_4/notes/PAPER_CLOSEOUT_THEOREM_PACKAGE.md` for the final task, hard pair, theorem, corollary, and limitations.
2. `code/source_snapshot/open_problem_forr4_passive_floor_consolidation_round_4/notes/Q64_COMPLETE_OUTWARD_LEDGER.md` and its JSON artifact for the one-batch numerical certificate.
3. `code/source_snapshot/open_problem_forr4_passive_floor_consolidation_round_4/notes/Q64_ADAPTIVE_TREE_FRONTIER_THEOREM.md` and its JSON artifact for the adaptive lift.
4. `code/source_snapshot/open_problem_forr4_passive_floor_consolidation_round_4/notes/ACTIVE_SIX_DOSE_RESOURCE_ROW.md` and `ACTIVE_SIX_ROBUSTNESS_GATE.md` for the active protocol and derived contrast requirement.
5. `code/source_snapshot/open_problem_forr4_passive_floor_consolidation_round_4/notes/EXPERIMENTAL_FEASIBILITY_DECISION.md` for the dated implementation assessment.
6. The preceding asymptotic paper `scratch/20260718_forr4_floor_paper/` for model vocabulary and broader proof context; its numerical citation style is not inherited.

## Main text

1. `S0_abstract.tex`: finite-size result, exact classes, numerical certificate, active protocol, and implementation qualification.
2. `S1_introduction.tex`: physical question, relationship to forrelation and the asymptotic theorem, finite-size motivation, and contribution summary.
3. `S2_model.tex`: sign masks, Sylvester transform, promise, hard-dose meter, passive class, active class, and memory boundary.
4. `S3_statement.tex`: principal theorem, numerical certificate, common-scale error interpretation, and proof dependency map.
5. `S4_active.tex`: folded-state protocol, overlap identity, majority error, dose ledger, and resource count.
6. `S5_hard_pair.tex`: signed-permutation plant, attenuation, conditioning, and promise loss.
7. `S6_one_batch.tex`: occupation reduction, physical coefficient registry, outward Perron certificate, and one-batch total.
8. `S7_adaptive.tex`: normalized direct-sum strategy factorization, temporal occupation pullback, two-law Perron interface, and conditioning.
9. `S9_implementation.tex`: exact contrast threshold, mode and transform requirements, primary platform evidence, and acceptance checklist.
10. `S10_discussion.tex`: interpretation, relation to the asymptotic theorem, limitations, and next mathematical and experimental questions.
11. `S8_verification.tex`: artifacts, directed rounding, independent tests, figure contract, and reproduction commands, placed last so verification does not interrupt the scientific narrative.

## Appendices

- `A_hard_pair.tex`: exact signed-permutation construction and concentration bound.
- `B_ledger.tex`: registry structure, coefficient rounding, 210-state matrix, and Collatz--Wielandt certificate.
- `C_adaptive.tex`: complete normalized feature-factorization proof.
- `D_notation.tex`: notation and claim-grade table.
- `E_verification.tex`: code map and complete command record.

## Completion gates

1. Every theorem constant agrees with the audited JSON artifacts.
2. Every conventional citation is resolved in `references.bib` and `dev/BIBLIOGRAPHY_NOTES.md` and is supported by a primary source.
3. `latexmk -pdf main.tex` completes without undefined references or citations.
4. `cd code && ./run_all.sh` passes.
5. `cd code && ./run_all.sh --full` passes.
6. Correctness, scope, and prose reviews are complete under `dev/reviews/`.
7. No manuscript sentence describes the $N=4096$ device as experimentally demonstrated or currently credible.
8. Every committed figure is regenerated from code, phase patterns use two-dimensional grids, and the final PDF passes page-by-page visual inspection.
