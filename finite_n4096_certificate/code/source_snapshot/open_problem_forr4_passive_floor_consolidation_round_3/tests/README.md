# Tests

Add regressions here only after an identity is accepted or a counterexample
must remain permanently detectable.  The round-three runner first executes
the complete frozen round-two suite, then every Python file in this folder.

alternating_double_endpoint_spectrum.py protects exact agreement with the
direct \(q=2,4\) values, full block rank, and the exact \(q=32\) spectrum
certificate.

weighted_alternating_q2_certificate.py protects the exact small-orbit Gram
decomposition and the certified numerical upper bound for the \(q=2\)
arbitrary-diagonal stress test.

weighted_double_endpoint_contraction.py validates the mixed Fourier reduction
against direct \(q=2\) matrices and protects both refined \(q=32\)
alternating bounds.  same_middle_weighted_bound.py checks the exact endpoint
orbit spectra and direct \(q=2\) reductions.  The occupation regression
protects the complete 64-cut supporting upper below the isolated
double-endpoint margin.

degree_six_joint_occupation_barrier.py protects the negative A3 result: one
shared occupation law does not yet rescue the four degree-six profiles even
after inserting the refined endpoint- and middle-cubic fixed-mask bounds.
single_cubic_weighted_bound.py independently checks those endpoint formulas
against direct \(q=2\) matrices.

occupation_compatible_sector_optimization.py verifies the exact occupation
pairing, all new incidence and chained-slice tables against direct \(q=4\)
moment enumeration, the degree-eight checkpoint, the four triple-cubic
contractions, the endpoint-quintic slices, both separated cubic--quintic
pairs, the sharp middle-quintic four-cycle saturator, and the updated
remaining degree-ten/twelve target.  It also checks the new xor-labelled
endpoint formulas, the fixed-orbit Gram closed form through \(q=32\), and
the exhaustive \(q=4\) pure-orbit maximum for the opposite endpoint chain.
opposite_endpoint_mixed_orbit_obstruction.py protects the explicit mixed law
and its \(q=32\) coefficient \(0.0395939553\).  It now also verifies that the
superseding accepted-sector repair gives more than \(0.0495\) slack on that
formerly failing restricted diagnostic.  The filename is retained to
preserve the history of the falsification.
attenuation_promise_concentration.py protects the biased-sign moment
generating-function constant, monotone finite-tilt branch, Euclidean packing,
four reverse-martingale coefficient norms, and hybrid exact/global promise
loss.  It also protects the two-split gate and extended optimized promise
loss.

repaired_open_profile_budget.py protects the 43-profile/23-orbit inventory,
the updated repaired/hybrid/extended baselines, preservation of all three
physical-orbit diagnostics, and the coarse all-open \(1/q\) target
\(0.322669154\).

high_degree_record_incidence_frontier.py checks the generalized odd-record
incidence dynamic program against exhaustive \(q=4\) support enumeration,
protects the 6,016/14,624 inventory, records 1,138 generic \(1/q\)
certificates, verifies that zero of 888 dose-six-relevant balanced entries is
certified, and protects the leading Perron route priority.

leading_balanced_disjointness_contraction.py verifies the exact centered
factorization of the \(k\)-set/singleton disjointness mask, stress-tests the
complete chain theorem on constructed small-dimensional tensors with
arbitrary diagonal weights, and protects the \(N=1024\) coefficient and
updated ledger slack.

adjacent_balanced_cubic_slice_contraction.py verifies the exact fixed-pair
cubic slice by complete \(q=4\) signed-permutation enumeration, stress-tests
the abstract adjacent chain under arbitrary diagonal weights, confirms the
generic two-mask failure, and protects the second coefficient and
two-theorem ledger slack.

separated_balanced_endpoint_slice_contraction.py independently enumerates the
\(q=4\) cubic fixed-pair moment slice and the extracted quintic fixed-four
scalar slice, including the essential factor of \(N\).  It stress-tests the
abstract separated-chain factorization under arbitrary diagonal laws and
protects the third coefficient \(0.1737428008\) and the narrow three-theorem
ledger slack.  It also protects the reranked next orbit and its reoptimized
coefficient gate.

internal_singleton_shared_law_contraction.py checks the complete \(q=4\)
cubic endpoint against the exact affine twirl, stress-tests nonsymmetric
laws, enumerates the quintic completion and overlap energy, protects the exact
\(q=32\) spectrum and fourth coefficient \(0.0250967461\), and verifies the
four-theorem ledger, next-orbit gate, and symmetric physical-law diagnostic.

column_cubic_quintic_row_contraction.py enumerates the complete \(q=4\)
quintic fixed-four row energy and compressed base chain, checks rank,
maximum entry, and saturation of the \(1/\sqrt N\) factor, stress-tests exact
target submatrices under correlated diagonal laws, and protects the fifth
coefficient, ledger, reranking, and next gate.

adjacent_balanced_row_slice_contraction.py enumerates every \(q=4\) row,
checks the record-one quintic extension counts and record-three no-even/even
counts through \(q=8\), verifies the exact horizontal Walsh tail for every
record-three cubic, and stress-tests sparse exact occurrence tensors under
correlated laws.  It protects the sixth coefficient \(0.0422410016\), the
six-theorem ledger, reranking, and next gate.

whole_cubic_quintic_triple_contraction.py exhausts the \(q=4\) fixed-three
quintic rows, checks the maximizing \(q=8\) row, verifies exact cubic
XOR-column compression under random laws and a saturating law, and
stress-tests sparse full occurrence tensors.  It protects the seventh
coefficient \(0.0370952793\), seven-theorem ledger, reranking, and next gate.

middle_cubic_quintic_pair_contraction.py constructs every \(q=4\) fixed-pair
row, checks the leading \(q=8\) symmetry representative, separately verifies
both endpoint sums and middle record maxima, and stress-tests sparse exact
occurrence tensors under correlated laws.  It protects the eighth coefficient
\(0.0285281523\), ledger, reranking, and next gate.

whole_cubic_middle_pair_contraction.py constructs every compatible \(q=4\)
middle row, enumerates both fixed-pair geometries through \(q=8\), checks the
record-one and record-three formulas, and stresses the residual Walsh
compression and sparse exact target tensors under correlated laws.  It
protects the ninth coefficient \(0.0250919472\), ledger, reranking, and next
gate.

double_endpoint_cubic_quintic_row_contraction.py constructs every \(q=4\)
scalar completion row, checks the leading \(q=8\) representative, verifies
the repeated \(H_N\otimes H_N\) residual under correlated laws, and
stress-tests sparse exact target tensors.  It protects the tenth coefficient
\(0.0462425962\), ten-theorem ledger, reranking, and next gate.

active_two_flag_collective_obstruction.py protects the complete \(N=4\)
endpoint enumeration, integer Gram spectrum, exact radical trace distance,
\(N=16\) endpoint tensor product, folded-state multiplicativity, and
isometric \(N=1024\) lift for the two-complete-flag obstruction.

general_theory_exponent_audit.py protects the exact \(1/24\) accepted
algebra, the corrected \(N^{-v/8}\Rightarrow N^{1/8}\) implication, the
sufficient high-level second-suppression \(N^{1/16}\) target, and the
\(\sigma=1\) layered interface family at every level four through twelve.

terminal_interpolation_sigma_one_witness.py protects the exact level-twelve
three-path terminal forest, its six all-new transfers, positive local weight,
reflection sensitivity, assigned suppression one, exclusion of the old star,
and all 180 valid transfer orders that produce the same terminal graph.

terminal_three_path_projective_repair.py protects all fifteen path entry
partitions, every exact path flattening at \(N=2,4\), the strong-path norm
one bound, both weak \(N^{-1/2}\) factors, and the combined \(N^{-1}\)
contraction with \(N^{-1/4}\) slack beyond the level-twelve target.

high_level_terminal_best_of_two_audit.py protects the 222-state exact
enumeration, all terminal and sensitivity counts, the unique level-twelve
three-path type, every assigned-sigma-one best-of-two score, the exact list
of three level-9/10 joint saturators, and the improved global exponent
\(1/20\).

level_ten_forest_mean_zero_repair.py protects the 200 exact coefficient
histories, the common centered derivative site, both possible odd derivative
weights, all 282 dangerous original partitions, all 5,295 retained extended
partitions, decay one in every repair branch, and the improved global
exponent \(1/18\).

level_nine_tree_centered_repair.py protects the 200 representative
histories, both forced derivative sites, seven dangerous partitions, exact
Type-A and reflected Type-B branch counts, all retained decay minima, the
all-fresh reflection cancellation class, and the improved global exponent
\(1/16\).

low_level_terminal_centered_repair.py protects all four initial potentials,
the 236-state/39-terminal union, every levelwise sensitive count and safe
minimum, the exact level-seven/six histories, all repair branch counts and
decay minima, the fresh-outer cancellation, and the current global exponent
\(1/12\).

level_twelve_contraction_sharpness.py protects the eleven grouped unit
vectors, direct masked contractions at three powers of two, the exact
\((N-1)^2(N-2)^3/N^6\) formula, the displayed initial triple's 180 orders,
all 1,080 positive histories over the 12 initial triples, matching
upper/lower exponent \(-1\), and the greater-than-99-percent retained ratio
at \(N=1024\).

adjacent_cubic_quintic_orbit_witness.py checks every relevant \(q=4\)
record-one and record-three link formula against complete enumeration,
protects the \(q=4,8\) physical orbit diagnostics, and verifies the exact
\(q=32\) record-slice coefficients and their combined local bound.

adjacent_mixed_orbit_diagnostics.py checks the exact \(q=4\) mixed Fourier
reduction against an independently assembled physical matrix and verifies
that the 25-class closed vertical-triple formula reproduces the complete
\(q=4\) construction.  The occupation-compatible regression separately
protects the chained adjacent record-three coefficient and a direct composed
slice-matrix inequality.

opposite_endpoint_transposed_vertical_witness.py checks the symmetry-reduced
physical witness against a direct \(q=4\) matrix, protects its \(q=8\)
value, and reruns the exact-numerator \(q=32\) dominant-class certificate.

Later tests may protect results from any track: finite-size certificates,
general contraction identities, hard-instance scorecard calculations, or
valid passive/active witnesses.  Each test should cite the theorem target or
decision it protects.
