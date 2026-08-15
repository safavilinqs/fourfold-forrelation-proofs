# Passive four-forrelation consolidation: round 3

## End-of-round status

Round 3 is closed as a research phase.  Start with
ROUND_3_RETROSPECTIVE.md for the consolidated outcome and Round 4 handoff.
The strongest theorem is \(\Omega(N^{1/12})\); passive hard dose greater
than six at \(N=1024\) remains open; and the generic limiting level-twelve
grouped graph norm is now proved exponent-sharp.  The historical plan,
portfolio, and open-problem files below are retained as the research record,
not as the live Round 4 plan.

The live successor is
../open_problem_forr4_passive_floor_consolidation_round_4/README.md.

## Charter

Round 3 is the phase that asks how large, useful, and close to optimal the
proved active--passive separation really is.  It is not defined by the
completion of one signed-permutation ledger.  Read `ROUND_3_CHARTER.md` for
the plain-language mission and `BOUNDARY_MAP.md` for the accepted lower and
upper evidence on both sides of the comparison.

## Mission

Determine how large and practically meaningful the active--passive sensing
separation for four-forrelation really is.

The inherited theorem proved a genuine asymptotic separation:

$$
D_{\mathsf P}^{\rm hard}=\Omega(N^{1/24}),
\qquad
D_{\mathsf A}^{\rm hard}\le 6.
$$

Round 3 now improves the passive side to

$$
D_{\mathsf P}^{\rm hard}=\Omega(N^{1/12}).
$$

Round three should strengthen the science behind that statement rather than
tie success to one inequality or one hard instance.  It has three tracks:

1. **Realistic-size separation.**  Prove passive hard dose greater than six
   near \(N=1024\), beginning with the attenuated signed-permutation plant.
2. **Stronger general theory.**  Improve the asymptotic passive floor or
   extract a reusable contraction principle that explains it.
3. **Test the boundary.**  Search for alternative hard instances, passive
   counterprotocols, and improvements or lower bounds on the active side.

Track A has now cleared every one-batch sector through degree eight and a
substantial set of degree-ten/twelve sectors.  Auditing the newly leading
adjacent witness exposed a much larger improvement in an older accepted
record-three sector: chaining two exact squared-slice tables lowers the
degree-eight conditioned total from about \(0.297666\) to
\(0.248224750\).  With the proved high sectors included, the accepted partial
total is \(0.279758547\).  As a route-selection diagnostic, preserving all
three known physical-orbit charges and assigning \(1/q=1/32\) to every other
open split gives \(0.322669154\), leaving \(0.010664179\).  The current
concrete problem is therefore to prove or classify those coarse \(1/q\)-scale
bounds, not to polish one fragile compound factor.  A complete one-batch
upper bound and the adaptive lift are still required for passive \(>6\).

The first broad classification has now been completed.  Generic
record-incidence and rank--Frobenius bounds prove \(1/q\) for 1,138 of the
6,016 open entries, but all 1,138 are excluded by the exact dose-six
occupation shell; none of the 888 balanced relevant entries is certified.
The live finite-size problem is consequently a chain-aware contraction for
those balanced entries, ordered by actual Perron-ledger impact.  The generic
incidence route is retained as a theorem and stopped as the lead strategy.

The first ten chain-aware balanced-orbit theorems are now complete.  For
\((3,1,1,5):(0,1,0,4)\) and its complement/reversal cuts, a centered
distinct-label Schur factor followed by an exact cubic--Walsh collapse gives
the arbitrary-diagonal coefficient \(0.0934752746\).  This is about
\(2.99121/q\), not the provisional \(1/q\).  Inserting the safe value while
leaving the other open coefficients at their diagnostic targets gives
\(0.326281861<1/3\).

For \((1,1,3,5):(0,0,1,4)\), the exact internally split cubic slice has
energy \((q^2-2q+2)/(q^2(q-1))\).  Combining its square-root Schur factor
with the completed adjacent cubic--quintic link and the sole remaining
\(4|1\) mask gives coefficient \(0.0162724693<1/32\).  With both theorems
inserted, the diagnostic is \(0.325563858\), with \(0.007769475\) slack.
Eight balanced entries are rigorously controlled.

For \((3,1,1,5):(1,0,0,4)\), separating the cubic fixed-pair slice from the
extracted quintic fixed-four slice gives coefficient
\(0.1737428008\).  The normalization is less favorable than an initial
extra-\(1/q\) heuristic because extracting the quintic scalar multiplies its
squared energy by \(N\).  The safe coefficient still beats the generic
two-mask value \(0.2252930474\), but exceeds the provisional \(1/32\).
With all three theorems inserted, the diagnostic was \(0.3326651190\), leaving
only \(0.0006682144\).  Twelve balanced entries were controlled.

The fourth theorem resolves the resulting go/no-go orbit
\((3,1,1,5):(1,1,1,2)\).  Bounding its cubic and quintic endpoints
separately gives \(0.0991470258\) and fails the ledger.  Keeping the common
diagonal law instead permits an exact affine twirl of the cubic endpoint and
a completion-plus-collision bound for the quintic endpoint.  Their shared-law
contraction gives coefficient \(0.0250967461\), safely below the gate
\(0.0450405468\).  With all four theorems inserted, the diagnostic is
\(0.3323683002\), leaving \(0.0009650332\); sixteen of 888 balanced entries
are controlled.  The reranked next priority is
\((3,1,1,5):(0,0,1,4)\), with coefficient gate
\(0.0542506298=1.73602015/q\).

The fifth theorem resolves that orbit without paying its generic
distinct-label factor.  For \((3,1,1,5):(0,0,1,4)\), the exact quintic
fixed-four row energy is \(1-4/N\).  Removing that Schur factor leaves an
\(N\)-row cubic--Hadamard matrix with entries at most \(1/N\), so
rank--Frobenius gives coefficient
\(\sqrt{1-4/N}/\sqrt N=0.0311889051<1/32\).  With all five theorems inserted,
the diagnostic is \(0.3323657941\), leaving \(0.0009675392\); twenty of 888
balanced entries are controlled.  The next priority is
\((1,1,3,5):(0,1,1,3)\), with gate
\(0.0570749885=1.82639963/q\).

The sixth theorem passes that gate by retaining the whole adjacent row until
after the record sectors are separated.  For
\((1,1,3,5):(0,1,1,3)\), exact L-shape incidence controls record one, while
the record-three output is classified by the row pattern of the fixed
quintic triple.  The potentially leading horizontal pattern has exact Walsh
tail \(3/[(q-1)(q-2)]\).  The resulting arbitrary-law coefficient is
\(0.0422410016\).  With all six theorems inserted, the diagnostic is
\(0.3327757792\), leaving \(0.0005575541\); twenty-four of 888 entries are
controlled.  The next priority is \((3,1,1,5):(0,1,1,3)\), with gate
\(0.0484819899\).

The seventh theorem passes that gate by extracting only the split quintic
endpoint.  For \((3,1,1,5):(0,1,1,3)\), the exact fixed-three quintic energy
is \(22365/15872\).  The remaining whole cubic endpoint has XOR-labelled
Walsh columns, so duplicate compression and Schatten Hölder contribute
exactly \(1/q\).  The resulting arbitrary-law coefficient is
\(0.0370952793\).  With all seven theorems inserted, the diagnostic is
\(0.3329643636\), leaving \(0.0003689697\); twenty-eight of 888 entries are
controlled.  The next priority is \((1,3,5,1):(0,2,2,1)\), with gate
\(0.0426269309\).

The eighth theorem passes that gate using a chained row-energy estimate.  For
\((1,3,5,1):(0,2,2,1)\), the two endpoint sums are the exact cubic and
quintic fixed-pair energies.  Taking the universal middle-link maximum only
after preserving those sums gives coefficient
\(0.0285281523<1/32\).  With all eight theorems inserted, the diagnostic is
\(0.3328775891\), leaving \(0.0004557442\); thirty-two of 888 entries are
controlled.  The next priority is \((1,1,3,5):(0,0,3,2)\), with gate
\(0.0454321892\).

The ninth theorem passes that gate while retaining the preceding Walsh
chain.  For \((1,1,3,5):(0,0,3,2)\), normalize the complete fixed-pair
\(M_{35}\) row, then collapse the residual \(H_NM_{13}\) matrix by cubic XOR
and duplicate singleton labels.  Exact record-pattern counts give
coefficient \(0.0250919472<1/32\).  With all nine theorems inserted, the
diagnostic is \(0.3326862124\), leaving \(0.0006471209\); thirty-six of 888
entries are controlled.  The next priority is
\((1,3,5,1):(0,2,3,0)\), with gate \(0.0529177167\).

The tenth theorem passes that gate by keeping both endpoint Walsh factors.
For \((1,3,5,1):(0,2,3,0)\), the scalar cubic--quintic completion row splits
into record-one incidence and record-three endpoint-slice energies; the
residual is a repeated, column-twisted \(H_N\otimes H_N\) matrix with
arbitrary-law coefficient one.  The resulting coefficient is
\(0.0462425962<0.0529177167\), although it is larger than the provisional
\(1/32\).  With all ten theorems inserted, the diagnostic is
\(0.3331326055\), leaving only \(0.0002007278\); forty of 888 entries are
controlled.  The next priority is \((1,1,5,3):(0,1,3,1)\), with gate
\(0.0379251204\).  This narrow margin reinforces the charter: another orbit
is justified only as a bounded gate test, alongside the general-theory and
boundary projects.

The scheduled general-theory audit is also complete.  The accepted proof
already charges a level-\(v\) diagram by \(D^v\) and retains an exact
graph-dependent factor \(N^{-\sigma/2}\).  A checked legal interface family
has \(\sigma=1\) at every level, so existing component bookkeeping alone
cannot improve \(N^{1/24}\).  The corrected \(N^{1/16}\) target needs new
suppression only at levels nine through twelve; a uniform \(N^{-v/8}\)
theorem would now imply \(N^{1/8}\) directly.

The follow-up terminal-image test is also complete.  The earlier relaxed
star violates an exact outer-boundary degree cap and is not terminal.
Nevertheless, six all-new legal Stein transfers produce a genuine
level-twelve sensitive diagram consisting of three disjoint four-layer
paths.  Its local scalar weight is positive, and pairing two marks in one
path gives exact assigned-fiber \(\sigma=1\).  Thus neither terminal-image
exclusion nor scalar branching-sign cancellation improves the exponent.
The next physical gate nevertheless passes.  A safe all-projective
contraction bounds every occupancy-\((2,1,1)\) placement of this forest by
\(N^{-1}\): the paired path costs at most one and each singleton path gives
\(N^{-1/2}\).  This is stronger than the level-twelve \(N^{-3/4}\) target
and never converts a Hilbert auxiliary to projective mass.  Track B now
has classified all true level-nine-through-twelve terminal diagrams using
the better of the assigned and grouped-entry projective contractions.  Of
eight sensitive types, five pass the \(N^{1/16}\) row and three initially
remain joint saturators: two reflected level-nine trees and one level-ten
\(6+4\) forest.  All level-eleven and level-twelve types gain \(N^{-1}\),
first proving \(\Omega(N^{1/20})\).  The level-ten type is now repaired as
well.  Exact enumeration of all 200 coefficient histories finds a forced odd
centered derivative at its degree-three middle vertex.  One more Stein
expansion, followed by exact checks of 282 dangerous and 5,295 extended
physical partitions, gives \(N^{-1}\) for every retained branch.  This first
raises the floor to \(\Omega(N^{1/18})\).  The two reflected level-nine
trees are now repaired too.  All 200 legal histories for the upper-branching
representative force an odd centered endpoint derivative; its exact
expansion checks seven dangerous partitions and 2,906 Type-A branch
placements.  A direct reflected audit checks the retained Type-B branches
and proves that its 1,905 all-fresh placements cancel.  The current general
floor is therefore \(\Omega(N^{1/16})\).  The limiting rung has moved to
the low-level image.  The subsequent complete four-initial-state audit finds
that all four sensitive level-eight types already have \(N^{-1}\) decay.  It
isolates one level-seven and one level-six saturator; all 36 and four of their
respective legal histories force a common odd centered derivative.  Exact
re-expansion repairs every retained branch and cancels the fresh outer
level-six branch.  Thus every sensitive level from five through twelve has
decay \(N^{-1}\), while level four retains \(N^{-1/2}\).  The current general
floor is \(\Omega(N^{1/12})\), now limited by the level-twelve \(N^{-1}\)
row.  That grouped graph norm is exponent-sharp: a legal placement has an
explicit full-distinctness contraction
\(N^{-1}(1-N^{-1})^2(1-2N^{-1})^3\), and all 1,080 histories from its 12
initial triples have positive coefficient.  Further asymptotic progress
therefore needs physical frame structure or a different witness.  A useful
\(N=1024,D=6\) result remains open.

The first active-boundary audit is complete.  An exact tensor-product
endpoint ensemble at \(N=16\), lifted isometrically to \(N=1024\), gives
trace distance \(0.2776778892\) between the two-copy positive and negative
average flag states.  Hence even the optimal collective POVM on two complete
path--mode flags has error \(0.3611610554>1/3\).  Removing the third flag,
changing the classical decoder, or jointly measuring the two completed
flags cannot reduce the active bound below six.  A genuinely interleaved
five-traversal protocol remains open.

The shared output of all three tracks is a characterization of the gap: its
best proved scale, its behavior at realistic sizes, the best lower and upper
evidence for both passive and active sensing, the robustness of the chosen
hard instances, and the physical or mathematical mechanism behind it.

## Start here

1. `ROUND_3_CHARTER.md` states the broad goal in plain language.
2. `BOUNDARY_MAP.md` keeps passive and active lower/upper evidence together.
3. `notes/round3_broadened_mission_guide.md` explains how the charter applies
   to the current technical frontier.
4. `PROGRAM_MAP.md` explains how this round fits the original program.
5. `MISSION_LEDGER.md` maps the broad questions to evidence and deliverables.
6. `ROUND_2_HANDOFF.md` records the inherited facts and finite-size numbers.
7. `GOAL.md` defines the scope and success ladder.
8. `PLAN.md` gives the three-track work program and decision gates.
9. `OPEN_PROBLEMS.md` lists theorem targets across all three tracks.
10. `PORTFOLIO.md` gives bounded live projects and their stop rules.
11. `STATUS.md` records the current lead project and the next portfolio review.
12. `notes/triple_cubic_chained_slice_contraction.md` records the current
   triple-cubic result, and `notes/quintic_slices_and_separated_chain.md`
   records the new cubic--quintic frontier.
13. `notes/opposite_endpoint_orbit_factorization.md` records the exact
    fixed-orbit compound contraction, and
    `notes/opposite_endpoint_mixed_orbit_obstruction.md` records the exact
    physical witness that falsified the old scalar ledger.
14. `notes/euclidean_promise_concentration.md` proves the finite-tilt
    concentration repair and states exactly what remains open.
15. `notes/repaired_open_profile_budget.md` reconstructs the repaired open
    ledger, ranks the next contractions, and states the decision-critical
    \(q=32\) computation.
16. `notes/transposed_dominant_class_and_hybrid_repair.md` records the
    certified second witness, the hybrid concentration theorem, and the
    updated frontier.
17. `notes/adjacent_record_three_chained_repair.md` proves the accepted-sector
    repair, states the new coarse completion target, and records the adjacent
    mixed-law stress tests.
18. `notes/adjacent_cubic_quintic_record_gate.md` preserves the historical
    tensor, exact record slices, and direct-contraction failure.
19. `notes/high_degree_record_incidence_frontier.md` records the rigorous
    generic classification, the dose-six shell obstruction, and the ranked
    balanced frontier.
20. `notes/leading_balanced_disjointness_contraction.md` proves the first
    arbitrary-law bound on the dose-six balanced frontier and records its
    updated ledger cost.
21. `notes/adjacent_balanced_cubic_slice_contraction.md` proves the second
    balanced-orbit theorem using an exact cubic fixed-pair slice and records
    the two-theorem ledger.
22. `notes/separated_balanced_endpoint_slice_contraction.md` proves the third
    balanced-orbit theorem, explains the quintic normalization, and records
    the narrow three-theorem ledger margin.
23. `notes/internal_singleton_shared_law_contraction.md` proves the fourth
    balanced contraction and explains why preserving the common physical law
    succeeds where separate endpoint maxima fail.
24. `notes/column_cubic_quintic_row_contraction.md` proves the fifth balanced
    contraction using exact quintic row energy and a rank--Frobenius collapse.
25. `notes/adjacent_balanced_row_slice_contraction.md` proves the sixth
    balanced contraction by combining whole-row Schur factorization, exact
    record incidence, and horizontal Walsh orthogonality.
26. `notes/whole_cubic_quintic_triple_contraction.md` proves the seventh
    balanced contraction using the exact fixed-three quintic energy and
    XOR compression of the opposite whole cubic endpoint.
27. `notes/middle_cubic_quintic_pair_contraction.md` proves the eighth
    balanced contraction by chaining exact endpoint fixed-pair energies
    through the universal middle-link maximum.
28. `notes/whole_cubic_middle_pair_contraction.md` proves the ninth balanced
    contraction by normalizing the complete middle row and preserving the
    residual XOR/Walsh compression.
29. `notes/double_endpoint_cubic_quintic_row_contraction.md` proves the tenth
    balanced contraction by extracting the scalar completion row while
    preserving a tensor product of both endpoint Walsh factors.
30. `notes/general_theory_exponent_audit.md` derives the corrected exponent
    ladder, exact suppression parameter, interface saturator, and precise
    posterior-stable replacement lemma.
31. `notes/terminal_interpolation_sigma_one_witness.md` resolves the first
    terminal-image fork: the old star is excluded, but an exact positive
    level-twelve three-path forest still has assigned suppression one.
32. `notes/terminal_three_path_projective_repair.md` proves that the same
    forest safely contracts as \(N^{-1}\) in one all-projective norm regime,
    so it is not the next asymptotic saturator.
33. `notes/high_level_terminal_best_of_two_audit.md` enumerates the complete
    high-level terminal image, isolates three joint saturators, and proves
    the improved \(\Omega(N^{1/20})\) passive floor.
34. `notes/level_ten_forest_mean_zero_repair.md` uses the forced centered
    terminal weight to remove the level-ten saturator and proves the intermediate
    \(\Omega(N^{1/18})\) passive floor.
35. `notes/level_nine_tree_centered_repair.md` repairs both reflected
    level-nine trees, closes the high-level obstruction list, and proves the
    intermediate \(\Omega(N^{1/16})\) passive floor.
36. `notes/low_level_terminal_centered_repair.md` enumerates every initial
    collision pattern, repairs the level-seven/six types, and proves the
    current \(\Omega(N^{1/12})\) passive floor.
37. `notes/active_two_flag_collective_obstruction.md` proves that two
    complete active flags fail even under arbitrary collective decoding at
    \(N=1024\), and isolates the genuinely interleaved five-query frontier.
38. `notes/general_theory_loss_ledger.md`,
    `notes/hard_instance_scorecard.md`, and
    `notes/boundary_search_agenda.md` maintain Tracks B and C.
39. `PROJECT_BRIEF_TEMPLATE.md` is the common intake form for a new route or
    calculation.

`BASELINE_RUN.md` records the successful initialization check.  Run the
inherited baseline with:

```sh
./run_round3_checks.sh
```

Round two is a frozen evidence base.  Do not silently weaken or rewrite its
claims in this folder.
