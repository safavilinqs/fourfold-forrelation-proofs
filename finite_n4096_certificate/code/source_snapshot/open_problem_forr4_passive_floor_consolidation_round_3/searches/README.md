# Searches

Place deterministic exploratory calculations here.  Each script must name
the inequality, hard-instance decision, counterprotocol hypothesis, or active
resource claim it is testing.  It must distinguish proof from numerical
evidence and print enough normalization data to reproduce its conclusion.

The first accepted artifact is
alternating_double_endpoint_spectrum.py.  It corrects the unordered-pair
orbit reduction, matches the direct \(q=2,4\) calculation, and prints an
exact integer spectrum certificate at \(q=32\).

weighted_alternating_q2_certificate.py performs the first arbitrary-diagonal
falsification check for OP3-2.  It is explicitly limited to one fixed split
at \(q=2\).

weighted_same_orientation_certificate.py and
mixed_endpoint_weighted_bound.py now control the two relevant alternating
orientations at \(q=32\).  same_middle_weighted_bound.py gives exact
coefficients for the two same-middle classes, and
double_endpoint_occupation_optimization.py combines all 64 cuts against one
law on all 210 dose-six occupation states.

degree_six_joint_occupation_optimization.py is the first A3 ledger.  It
combines the exact endpoint- and middle-cubic fixed-mask bounds with the
refined double endpoint.  The resulting upper overshoots by \(0.1893\), and
more than 99 percent of its degree-six cost comes from the two endpoint-cubic
profiles.  single_cubic_weighted_bound.py records the exact fixed-mask
endpoint constants.  The next target is their shared-law compound
contraction.

occupation_compatible_sector_optimization.py supersedes that coarse barrier
by retaining the exact row/column occupation pairing.  The resulting
objective is a 210-state Perron eigenvalue and now includes every record
sector through degree eight, the proved high endpoint profiles, and the four
leading triple-cubic profiles, plus both separated cubic--quintic reversal
pairs.  It also records exact degree-five endpoint slices and the sharp
middle-quintic four-cycle incidence.  Its accepted partial pass uses
\(\beta=781/1000\).

opposite_endpoint_orbit_scan.py records the first compound-contraction
search.  It proves the xor-labelled cubic/quintic endpoint formulas and
evaluates the critical fixed-orbit block through two length-\(N\) Walsh
spectra.  Its \(q=32\) aligned value is a pure-orbit diagnostic, not an
arbitrary-orbit ledger coefficient.

opposite_endpoint_mixed_orbit_q4.py derives the exact twirled mixed-orbit
formula and searches its small-\(q\) laws.
opposite_endpoint_vertical_mixture_witness.py reduces an explicit physical
law to fourteen frequency classes through \(q=32\), records the old
forced-cut failure, and reoptimizes its finite-tilt promise repair.
opposite_endpoint_repair_scorecard.py ranks the promise, accepted-sector,
and joint-ledger repair targets after that obstruction.
promise_tail_monte_carlo.py is the seeded, explicitly nonrigorous diagnostic
supporting the decision to prioritize a sharper promise theorem.
attenuation_promise_concentration.py implements the global Kearns--Saul
bound, finite-tilt Euclidean certificate, and hybrid theorem that switches
blocks to the global parabola only after they leave the exact branch.  It
also implements the proved two-split extension of the exact Euclidean branch.
endpoint_cubic_physical_orbit.py records the small-\(q\) falsification of the
old unpaired occupation relaxation.

repaired_open_profile_budget.py reconstructs the high-degree ledger after
the physical witness diagnostics and promise repairs.  It preserves the
historical thresholds and now prints the superseding coarse all-open target:
coefficient \(1/q\), optimized total \(0.322669154\), and slack
\(0.010664179\).  This is a proof target, not an arbitrary-law certificate.

high_degree_record_incidence_frontier.py gives the first rigorous broad
classification of that target.  It inventories 43 open profiles, 6,016
profile/splits, and 14,624 record-sector/splits; proves exact one-axis
odd-record incidence by dynamic programming; and combines maximum-entry,
rank, and incidence bounds.  It certifies \(1/q\) for 1,138 entries, but the
exact dose-six shell audit shows that none of the 888 ledger-relevant
balanced entries is among them.  It also ranks the 234 relevant symmetry
orbits by their Perron contribution.  The ranking is diagnostic; the
coefficient classification is rigorous.

leading_balanced_disjointness_contraction.py gives the first rigorous
arbitrary-law upper on the balanced frontier.  For the four-cut
\((3,1,1,5):(0,1,0,4)\) orbit it combines an exact cubic--Walsh collapse
with a centered \(4|1\) distinct-label Schur factor.  The safe \(q=32\)
coefficient is \(0.0934752746\); inserting it raises the provisional
completion diagnostic to \(0.3262818609\), still below \(1/3\) by
\(0.0070514724\).

adjacent_balanced_cubic_slice_contraction.py gives the second rigorous
arbitrary-law upper on the balanced frontier.  For the four-cut
\((1,1,3,5):(0,0,1,4)\) orbit, the exact cubic fixed-pair slice replaces a
failing generic \(1|2\) mask factor; completing the adjacent link and
restoring only the \(4|1\) mask gives coefficient
\(0.0162724693<1/32\).  With both chain-aware theorems inserted, the
diagnostic is \(0.3255638580\), below \(1/3\) by \(0.0077694754\).

separated_balanced_endpoint_slice_contraction.py gives the third rigorous
arbitrary-law upper on the balanced frontier.  For the four-cut
\((3,1,1,5):(1,0,0,4)\) orbit, cubic fixed-pair and quintic fixed-four
slices give coefficient \(0.1737428008\).  The script explicitly records why
extracting the quintic scalar removes a tempting extra \(1/q\), compares the
safe result with the generic two-mask value \(0.2252930474\), and inserts all
three theorems.  The resulting diagnostic is \(0.3326651190\), below
\(1/3\) by only \(0.0006682144\).  It also reranks the frontier and computes
the next-orbit coefficient gate \(0.0450405468=1.44129750/q\).

internal_singleton_shared_law_contraction.py gives the fourth rigorous
arbitrary-law upper on the balanced frontier.  For
\((3,1,1,5):(1,1,1,2)\), it shows that separate endpoint maxima fail, then
computes the exact affine-twirled cubic factor and combines it with a
completion-plus-collision quintic Schur bound.  The coefficient
\(0.0250967461\) gives diagnostic \(0.3323683002\), reranks the remaining
frontier, and computes the next gate \(0.0542506298=1.73602015/q\).  It also
evaluates a symmetric physical-law family as a lower-witness stress test.

column_cubic_quintic_row_contraction.py gives the fifth rigorous
arbitrary-law upper on the balanced frontier.  For
\((3,1,1,5):(0,0,1,4)\), it replaces the failing generic distinctness
multiplier by the exact quintic fixed-four row energy and contracts the
remaining cubic--Hadamard matrix by rank--Frobenius.  The coefficient
\(0.0311889051<1/32\) gives diagnostic \(0.3323657941\), reranks the frontier,
and computes the next gate \(0.0570749885=1.82639963/q\).

adjacent_balanced_row_slice_contraction.py gives the sixth rigorous
arbitrary-law upper on the balanced frontier.  For
\((1,1,3,5):(0,1,1,3)\), it factors the complete \(M_{13}M_{35}\) row,
derives exact record-one and record-three support bounds, checks every
\(q=4\) row and selected \(q=8\) rows, and proves coefficient
\(0.0422410016\).  It inserts all six theorems, obtains diagnostic
\(0.3327757792\), and computes the next gate \(0.0484819899\).

whole_cubic_quintic_triple_contraction.py gives the seventh rigorous
arbitrary-law upper on the balanced frontier.  For
\((3,1,1,5):(0,1,1,3)\), it extracts the exact fixed-three quintic row
energy, compresses the opposite whole cubic endpoint by support XOR, proves
coefficient \(0.0370952793\), inserts all seven theorems to obtain diagnostic
\(0.3329643636\), and computes the next gate \(0.0426269309\).

middle_cubic_quintic_pair_contraction.py gives the eighth rigorous
arbitrary-law upper on the balanced frontier.  For
\((1,3,5,1):(0,2,2,1)\), it chains the exact cubic and quintic fixed-pair
squared energies through the universal \(M_{35}\) maximum.  The coefficient
\(0.0285281523<1/32\) gives diagnostic \(0.3328775891\) and computes the
next gate \(0.0454321892\).

whole_cubic_middle_pair_contraction.py gives the ninth rigorous arbitrary-law
upper on the balanced frontier.  For \((1,1,3,5):(0,0,3,2)\), it normalizes
the complete fixed-pair \(M_{35}\) row, preserves the residual XOR-labelled
\(H_NM_{13}\) chain, and counts both record sectors exactly.  The coefficient
\(0.0250919472<1/32\) gives diagnostic \(0.3326862124\) and computes the next
gate \(0.0529177167\).

double_endpoint_cubic_quintic_row_contraction.py gives the tenth rigorous
arbitrary-law upper on the balanced frontier.  For
\((1,3,5,1):(0,2,3,0)\), it extracts the scalar cubic--quintic completion
row and preserves both endpoint Walsh factors as a repeated
\(H_N\otimes H_N\) residual.  The coefficient
\(0.0462425962<0.0529177167\) gives diagnostic \(0.3331326055\) and computes
the next gate \(0.0379251204\).

active_two_flag_collective_obstruction.py gives the first exact Track C
active-boundary obstruction.  It enumerates the \(N=4\), \(F=\pm1/2\)
factor ensembles, tensors them into an \(N=16\) endpoint pair, constructs an
integer two-copy moment matrix, and factors sixteen exact Gram blocks.  The
resulting trace distance \(0.2776778892<1/3\) rules out arbitrary collective
decoding of two complete flags; a fixed unit tensor factor lifts the result
isometrically to \(N=1024\).

general_theory_exponent_audit.py reconstructs the accepted Track B exponent
ledger with exact rational arithmetic.  It records the retained
\(N^{-\sigma/2}\) graph factor, corrects the stale \(D^{2v}\) research rung,
checks that \(N^{-v/8}\) would now imply \(N^{1/8}\), and constructs a legal
four-layer tree/physical-entry placement with \(\sigma=1\) for every level
four through twelve.  The family limits the current graph interface; it is
not asserted to occur in the terminal interpolation image.

terminal_interpolation_sigma_one_witness.py resolves that terminal-image
question.  It replays six legal all-new Stein transfers that produce three
disjoint four-layer paths at level twelve, checks reflection sensitivity,
component ranks, projective and assigned suppression, and contrasts the true
outer-boundary degree cap with the excluded high-degree star.

terminal_three_path_projective_repair.py proves the next physical-norm gate
for that witness.  It enumerates all fifteen physical-entry partitions of a
four-vertex path, computes exact cut-rank exponents, and certifies that every
occupancy-\((2,1,1)\) placement of the three paths has a safe all-projective
\(N^{-1}\) contraction, stronger than the \(N^{-3/4}\) target.

high_level_terminal_best_of_two_audit.py enumerates the exact four-layer
Stein set dynamics up to layered relabeling.  Among 222 reachable states it
finds 34 terminal types and eight sensitive high-level types, exhausts every
dangerous singleton/pair physical partition, proves \(N^{-1}\) for all
level-11/12 types, isolates three level-9/10 joint saturators, and derives
the improved \(\Omega(N^{1/20})\) theorem.

level_ten_forest_mean_zero_repair.py resolves the unique level-ten
saturator.  It exhausts six initial configurations and 200 legal transfer
histories (four of the six initial triples contribute), proves that every
history has one odd centered derivative at the
same middle vertex, and checks the resulting graphs on 282 dangerous and
5,295 retained extended physical partitions.  Every branch gains
\(N^{-1}\), yielding the intermediate \(\Omega(N^{1/18})\) theorem.

level_nine_tree_centered_repair.py resolves both remaining reflected
level-nine trees.  It enumerates all 200 legal histories of the
upper-branching representative, proves two forced derivative sites in every
profile, and expands the uniform odd outer factor.  Exact grouped-entry
scoring covers seven dangerous partitions, 2,906 Type-A branches, every
retained reflected Type-B branch, and the 1,905 all-fresh Type-B placements
that cancel.  The minimum retained decay is \(N^{-1}\), yielding the current
\(\Omega(N^{1/16})\) theorem.

low_level_terminal_centered_repair.py enumerates the union of all four
initial collision patterns: 236 reachable states, 39 terminal types, and 22
sensitive types.  It proves direct decay \(N^{-1}\) for every sensitive
level-eight type, enumerates 36/four histories for the unique level-seven/six
saturators, and expands their forced centered derivatives.  Exact scoring
covers 2,404 level-seven branches, 31 retained level-six branches, and 144
fresh level-six placements that cancel.  Together with the high-level
repairs this yields the current \(\Omega(N^{1/12})\) theorem.

level_twelve_contraction_sharpness.py constructs eleven grouped unit vectors
for the legal paired-path placement and directly retains every same-layer
distinctness mask.  Its exact value is
\(N^{-1}(1-N^{-1})^2(1-2N^{-1})^3\), matching the accepted \(O(N^{-1})\)
upper exponent.  It also protects the 12 initial triples, all 1,080 positive
all-new histories, and the 0.9922113065 ratio to \(N^{-1}\) at \(N=1024\).

adjacent_cubic_quintic_orbit_witness.py reconstructs the newly leading
physical tensor, checks its direct and orthogonality-reduced \(q=4\) forms,
and derives exact record-one and record-three horizontal slices.  At
\(q=32\) their safe combined local coefficient is \(0.0168662459\), which
quantifies the required shared-law compound gain.

adjacent_combined_repair_scorecard.py identifies the accepted
adjacent-double-cubic family as the cheapest repair.  The audit leads to the
new chained record-three theorem in
occupation_compatible_sector_optimization.py.
adjacent_cubic_quintic_mixed_orbit_q4.py gives an exact twirled Fourier
reduction and direct mixed-law check at \(q=4\).
adjacent_vertical_triple_symmetric_witness.py validates a closed
25-frequency-class formula and shows that its strongest symmetric physical
family falls to \(0.000156598\) at \(q=32\).
adjacent_compound_frame_type_bound.py records a deliberately nonrigorous
route-selection diagnostic showing that the candidate type-block
rank--Frobenius architecture is too loose.

opposite_endpoint_transposed_vertical_witness.py evaluates one physical law
for the former leading unforced split.  Its exact symmetry reduction gives final
blocks of size \(q(q-1)^2\).  Its PSD compression, exact integer numerator
Gram, and explicit roundoff majorants certify the \(q=32\) dominant-class
coefficient lower \(0.0142810242\).

Future searches should be grouped by Track A, B, or C and bounded by an
explicit acceptance or falsification gate from `OPEN_PROBLEMS.md`.
