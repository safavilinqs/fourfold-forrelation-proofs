# Round 2 plan

Status: closed on 2026-07-14.  Phases A--C succeeded after four substantive
repairs to the round-one proof.  Phase D produced useful finite-size
contractions and counterexamples but did not certify passive hard dose
greater than six at \(N=1024\).  See ROUND_2_SUMMARY.md for the outcome and
the round-three folder for continuing work.

## Phase A — Freeze and isolate the claim

- Pin the exact round-1 theorem, proof files, and commit used as input.
- Extract only the reverse-tree contraction and its prerequisites into a dependency map.
- Write the contraction as a finite tensor-network statement, including types, normalization, masks, base fibers, histories, and terminal functional.
- List every place a dimension, outcome count, time integral, or adaptive branch could introduce a hidden factor.

Deliverable: notes/reverse_tree_spec.md.

## Phase B — Try to break the reverse tree

- Re-derive the contraction independently from the leaves upward and from the root downward.
- Attack the PSD majorant, vectorization, cross-entry-fiber treatment, frame insertion, and vector-valued bilateral lemma separately.
- Enumerate all minimal diagrams and small dimensions where the claimed invariant is meaningful.
- Generate adversarial tensors and masks, concentrating on equality and near-equality cases rather than generic random samples.
- Compare the claimed bound with direct contraction and with independently formulated convex or semidefinite relaxations when available.
- Check adaptive histories and dose partitions for missing branch, simplex-volume, factorial, or depth multiplicities.

Deliverables: tests/, certificates/, and notes/counterexample_log.md.

## Phase C — Make a confidence decision

Classify every dependency as:

- rigorously proved;
- independently reproduced;
- exhaustively checked in a finite regime;
- numerically stress-tested only; or
- unresolved.

The $\Omega(N^{1/24})$ result is accepted only if every theorem-critical item is rigorously proved and the independent checks agree. Otherwise record the smallest failing instance or the exact missing lemma and stop propagating the claim.

Deliverable: CONFIDENCE_REPORT.md.

## Phase D — Target realistic sizes

The present asymptotic form must not be treated as useful at $N\sim10^3$ merely because it grows with $N$. First compute explicit constants and the exact transcript-distance threshold needed to rule out dose 6 at $N=1024$.

Then pursue, in order:

1. a stronger contraction with substantially better $N$ suppression and explicit constants;
2. removal of avoidable dose/depth losses in the current reverse-tree argument;
3. a modified or different hard instance whose low-dose indistinguishability is strong at finite $N$;
4. rigorous computer-assisted certificates for the finite-size regime, if they can be made uniform over all passive protocols.

The milestone is concrete: a proved transcript-distance bound below the discrimination threshold for every passive protocol of dose at most 6 at a stated realistic size near $10^3$.

Deliverable: REALISTIC_SIZE_PROGRAM.md and, if successful, a machine-checkable certificate.

## Explicit non-goals

- Manuscript-style polishing.
- Reorganizing prose that does not affect correctness.
- Large parameter sweeps without a stated falsification target.
- Reporting heuristic numerics as a proof.
