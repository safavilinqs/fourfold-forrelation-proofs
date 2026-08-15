# Line-item audit of the repaired reverse-tree contraction

Date: 2026-07-14

Verdict: **PASS for the finite-skeleton contraction after replacing the round-one mixed-component proof by the global dichotomy and expanding the distinct-label mask before tensor multiplication.** This verdict concerns the reverse-tree interface, not the usefulness of its constants at realistic $N$.

## Audit table

| Item | Check | Result |
|---|---|---|
| Dichotomy is exhaustive | Either every component-entry count is at most one or some maximum is at least two | Pass |
| Case I party structure | At most one vertex per component per physical entry makes vertical grouping valid | Pass |
| Distinct Fourier labels | Walsh-expand every same-layer inequality mask; coefficient mass is at most $2^{P_G}$ | Pass after RT-006 repair |
| Case I graph norm | After the mask expansion, component injective norms are bounded by natural-cut norms and vertical norms multiply termwise | Pass |
| Case I frame norm | Reverse frame expansion is projective over grouped physical entries from start to finish | Pass |
| Case I suppression | A four-layer component exists and contributes $N^{-1/2}$ | Pass |
| Case II assignment | Every component can be assigned to a maximum-occupancy entry with $k_C\ge1$ | Pass |
| Weak assignment | $v_C-1-e_C\le0$, so $k_C=1$ creates no dimension loss | Pass |
| Strong assignment | $v_C-k_C-e_C\le-1$ for $k_C\ge2$ | Pass |
| Case II suppression | At least one strong assignment exists by the case hypothesis | Pass |
| PSD majorant | Unit Hilbert preimage reshapes to an operator of norm at most one | Pass |
| Open frontier wires | Carry them in an operator-valued boundary; fix a unit input before the Hilbert reshape | Pass after operator-frontier repair |
| Cross-entry fibers | Fixed before the majorant; exact diagonal improvement cancels their multiplicity | Pass |
| Zero/one/two assigned entries at a node | Identity, unilateral, and bilateral complete-frame packing cover the three cases | Pass |
| Adaptive outcomes | Outcome-selected descendants are permitted by the Hilbert-valued lemma; unmarked nodes are stochastic | Pass |
| Mixed norm | Never used in Case II; never leaves projective norm in Case I | Pass; RT-003 avoided |
| Collision handling | Within-base fiber factor is factorial in at most twelve marks | Pass |
| Dose ledger | Exact base/ordered-mark identity keeps the factorial collision fiber inside one falling-factorial square mass | Pass |
| Adaptive marked-time sum | Tree-level binomial induction gives $(2D)^v$ with max/complete-frame control over children | Pass |
| Interpolation weights | After conditioning, terminal Stein weights are bounded local vertex multipliers | Pass |
| $N$-uniformity | Only the displayed graph diagonal factors depend on $N$ | Pass |

## Falsification evidence

- The unqualified collision lemma fails exactly as predicted, while its collision-aware repair reaches ratio one in the randomized suite.
- The invalid Hilbert-to-projective inference exhibits the predicted $\sqrt N$ gap through $N=16$; the repaired proof does not use it.
- The global dichotomy was checked on 69,632 exact named-graph placements and 100,000 random component systems.
- All 256 minimal-chain placements across two adaptive nodes pass at $N=2$.
- A spanning strong chain plus a nonspanning weak edge passes all 2,720 admissible $N=2$ two-node placements and six $N=4$ spot placements.
- The all-singleton side passes all relative physical-entry permutations of two spanning chains under repeated adaptive frame draws.
- An exact $N=4$ witness disproves literal multiplicativity after the distinct-label mask; the Walsh expansion repairs it with a diagram constant and no $N$ loss.
- Identity wires exhibit the exact $\sqrt N$ Hilbert--Schmidt gap; the operator-valued frontier regression checks the uniform local update.
- Exact rational support enumeration checks the combined collision/insertion identity through four marks and dose six.
- The adaptive dose recurrence passes 300 random trees through twelve marks and is unchanged after duplicating an outcome branch 100 times.

These computations are falsification evidence only. The pass verdict rests on the two analytic cases and their norm-compatible invariants.

## Residual risk

The greatest residual risk is transcription when this replacement is integrated into a manuscript: reintroducing the old partition into assigned and separately projective components would reintroduce RT-003. The proof must preserve the global dichotomy exactly.
