# Round 3 baseline run

Date: 2026-07-14

Command:

```sh
./run_round3_checks.sh
```

Interpreter: Python 3.13.13 from the recorded Conda audit environment.

Result: PASS.  The command ran the complete round-two reverse-tree and
realistic-size regression suite, including the repaired contraction,
counterexample regressions, joint occurrence packing, degree-eight slice
bounds, the double-endpoint Schur benchmark, and the exact minimal-chain
certificate.  It then reported:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

The newest accepted regression is summarized first; the development history
continues below.

## Latest validated result — sharp level-twelve graph-norm obstruction

The focused regression passed on 2026-07-15.  It constructs eleven unit
grouped-entry vectors for the legal paired-path placement, contracts every
same-layer distinctness mask directly at small powers of two, and reports:

```text
level-twelve contraction sharpness passed:
initials=12,histories=1080,displayed_orders=180,
positive=True,assigned_sigma=1,
upper_exponent=-1,lower_exponent=-1,
N1024_lower=0.000968956353997,
N1024_ratio=0.992211306493,barrier=1/12
```

The exact masked value is
\(N^{-1}(1-N^{-1})^2(1-2N^{-1})^3\).  It matches the accepted
all-projective \(O(N^{-1})\) upper exponent, so the arbitrary grouped graph
norm cannot improve the present \(N^{1/12}\) theorem.  The complete suite
then passed with this regression included:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## Latest validated result — complete low-level image and \(N^{1/12}\) theorem

The focused regression passed on 2026-07-15 after extending the terminal
image to all four initial covariance collision patterns and repairing the
unique low-level saturators.  It reported:

```text
low-level terminal centered repair passed:
reachable=236,terminals=39,sensitive=22,
level8_decay=1,
level7_histories=36,
level7_counts=(20, 106, 586, 1692),
level6_histories=4,
level6_cancelled=144,
global_exponent=1/12
```

All four sensitive level-eight types already have safe decay \(N^{-1}\).
The 36 level-seven histories force one common outer \(\psi''\); all 2,404
repair branches gain at least \(N^{-1}\).  The four level-six histories force
one common internal \(\psi''\); its 31 marked-outer branches gain \(N^{-1}\),
while 144 fresh-outer placements cancel.  Together with the accepted
high-level repairs, this proves the passive floor
\(\Omega(N^{1/12})\).  The bare \(N=1024\) scale is \(1.7817974363\), still
below six.  The complete inherited-plus-Round-3 suite then passed on
2026-07-15 with the new regression included:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## Latest validated result — reflected level-nine repair and \(N^{1/16}\) theorem

The complete suite passed again on 2026-07-15 after removing both reflected
level-nine joint saturators.  The new regression reported:

```text
level-nine reflected-tree centered repair passed:
histories=200,profiles=4,
derivative_sites=(('A', (0, 0)), ('A', (1, 0))),
type_a_counts=(14, 132, 855, 1905),
type_b_counts=(14, 132, 855, 1905),
global_exponent=1/16
```

Every legal history producing the upper-branching representative has a
forced odd centered derivative at the outer endpoint.  The exact expansion
scores all seven dangerous original partitions and 2,906 Type-A branch
placements.  The reflected implementation separately scores every retained
Type-B branch and verifies that its 1,905 all-fresh placements add a fourth
first-layer mark and cancel.  All retained branches gain at least
\(N^{-1}\).

This closes the complete high-level terminal obstruction list and proves the
rigorous passive floor \(\Omega(N^{1/16})\).  It still does not prove passive
dose greater than six at \(N=1024\), where the bare scale is only
\(1024^{1/16}=1.5422108254\).  The run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## Latest validated result — centered level-ten repair and \(N^{1/18}\) theorem

The complete suite passed again on 2026-07-15 after removing the unique
level-ten joint saturator.  The new regression reported:

```text
level-ten forest mean-zero repair passed:
potential_initials=6,contributing_initials=4,
histories=200,profiles=4,
centered_weights=(LocalWeight(kind='gamma', derivative=1),
LocalWeight(kind='stein', derivative=1)),
partitions=282,extended=5295,
duplicate_decay=1,bridge_decay=1,outer_decay=1,
global_exponent=1/18
```

Every legal history producing the level-ten \(6+4\) forest has one forced
odd centered derivative at the same middle vertex.  Applying its exact
Stein identity one more time leaves only reflection-sensitive graphs with
safe decay \(N^{-1}\).  The checker exhausts all 282 dangerous original
physical partitions and all 5,295 retained extensions.  The two reflected
level-nine trees are therefore the only remaining high-level saturators,
and the rigorous passive floor improves to \(\Omega(N^{1/18})\).  This is
still not a passive-dose-greater-than-six theorem at \(N=1024\), where the
bare scale is only \(1024^{1/18}=1.4697344923\).  The run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## Complete high-level audit and \(N^{1/20}\) theorem

The complete suite passed again on 2026-07-15 after repairing the true
terminal three-path witness and exhaustively classifying the reachable
high-level terminal image. The two new regressions reported:

```text
terminal three-path projective repair passed:
partitions=15,occupancy_classes={1: 1, 2: 9, 3: 4, 4: 1},
strong_worst=0,weak=-1/2,combined=-1,accepted=-1/2,target=-3/4,
N1024_bound=0.0009765625,N1024_target=0.00552427172802

high-level terminal best-of-two audit passed:
reachable=222,terminals=34,high=17,sensitive=8,sigma_one=8,
passing=5,joint_saturators=3,failure_levels=(9, 9, 10),
level12_projective=1,global_exponent=1/20
```

The first check proves an all-projective \(N^{-1}\) contraction for the true
level-twelve three-path forest, avoiding any mixed-frame reinterpretation.
The second performs the exact transfer-state enumeration and retains the
better of two separately complete contraction arguments. All level-eleven
and level-twelve terminal types gain \(N^{-1}\). Exactly three
proof-interface saturators remain: two reflected level-nine trees and one
level-ten forest with component sizes \(6+4\). Combining this classification
with the accepted \(N^{-1/2}\) bound at levels four through ten proves the
improved passive floor \(\Omega(N^{1/20})\). This is an asymptotic advance,
not yet a passive-dose-greater-than-six result at \(N=1024\); the bare scale
there is only \(N^{1/20}=\sqrt2\). The run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## True terminal sigma-one witness

The complete suite passed again on 2026-07-15 after resolving the first
terminal-image fork in Track B. The new regression reported:

```text
terminal interpolation sigma-one witness passed:
v=12,e=9,components=3,valid_transfer_orders=180,
first_layer=3,projective_sigma=3,assigned_sigma=1,
old_star_degree=9,true_boundary_cap=2
```

The check replays six legal all-new Stein transfers, verifies the resulting
three disjoint four-layer paths, exhausts all 720 transfer orders and finds
the 180 valid orders producing the same graph, and protects reflection
sensitivity and the physical assignment with suppression one. It also
confirms that the older high-degree star is excluded by the true
outer-boundary degree cap. Thus terminal support alone cannot provide the
extra high-level decay; the next Track B target is a posterior mixed-frame
contraction gaining at least \(N^{-1/4}\) on this explicit forest, or a
physical witness showing that gain is impossible. The run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## Tenth double-endpoint row contraction run

The complete suite passed again on 2026-07-15 after adding the tenth
chain-aware balanced-orbit theorem. The new regression reported:

```text
double-endpoint cubic-quintic row contraction passed:
q4_row=1.58693415638,
q8_row=0.00537738314818,
q32_row=0.00213837770744,
coefficient=0.0462425962446,
ledger_total=0.333132605485,
threshold_slack=0.000200727847845,
residual_worst=0.274541395665,
sparse_ratio=0.000750145028063
```

Complete \(q=4\) scalar-row construction, a direct \(q=8\) representative,
correlated-law Walsh-tensor stress, and sparse exact target tensors protect
the record split, double-endpoint normalization, and
\(H_N\otimes H_N\) residual compression. The new coefficient replaces the
\(1/32\) placeholder on four balanced entries and raises the certified count
to forty of 888 entries. Because \(0.0462425962>1/32\), the diagnostic margin
shrinks to \(0.0002007278\). The remaining 848 provisional coefficients,
physical diagnostic charges, and adaptive lift keep the result from being a
complete passive-dose-six certificate. The run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## Ninth whole-cubic middle-pair contraction run

The complete suite passed again on 2026-07-15 after adding the ninth
chain-aware balanced-orbit theorem. The new regression reported:

```text
whole-cubic middle-pair contraction passed:
q4_record_one=2.75,
q4_record_three=2.5,
q8_record_one=1.17015306122,
q8_record_three_same=0.466666666667,
q8_record_three_distinct=0.25,
residual_ratio=0.295255172788,
coefficient=0.0250919471547,
ledger_total=0.332686212434,
threshold_slack=0.000647120899295,
next_admissible=0.052917716668,
sparse_ratio=0.00237291328662
```

Complete \(q=4\) middle-row construction, pair-geometry enumeration through
\(q=8\), separate record-pattern entry tests, and correlated-law Walsh
stress protect the pair counts and the decisive residual \(1/q\)
compression. Sparse exact target tensors protect the orientation,
disjointness, and normalization. The new coefficient replaces the \(1/32\)
placeholder on four balanced entries and raises the certified count to
thirty-six of 888 entries. The remaining 852 provisional coefficients,
physical diagnostic charges, and adaptive lift keep the result from being a
complete passive-dose-six certificate. The run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## Eighth chained pair-slice contraction run

The complete suite passed again on 2026-07-15 after adding the eighth
chain-aware balanced-orbit theorem. The new regression reported:

```text
middle cubic-quintic pair contraction passed:
q4_row=0.115397805213,
q8_row=0.000289480284345,
q32_row=0.000813855473212,
coefficient=0.0285281522923,
ledger_total=0.332877589131,
threshold_slack=0.000455744202313,
next_admissible=0.0454321892133,
sparse_ratio=0.00238828774266
```

Exact endpoint pair-slice identities and a pointwise middle-link maximum give
the rigorous \(q=32\) row-energy bound. Complete \(q=4\) row enumeration,
selected \(q=8\) rows, endpoint-compatible record checks, and sparse exact
target tensors protect the factorization, incidence cases, and normalization.
The new coefficient replaces the \(1/32\) placeholder on four balanced
entries and raises the certified count to thirty-two of 888 entries. The
remaining 856 provisional coefficients, physical diagnostic charges, and
adaptive lift keep the result from being a complete passive-dose-six
certificate. The run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## Seventh balanced endpoint contraction run

The complete suite passed again on 2026-07-15 after adding the seventh
chain-aware balanced-orbit theorem. The new regression reported:

```text
whole-cubic quintic-triple contraction passed:
q4_slice=0.875,
q8_slice=1.16517857143,
cubic_random_worst=0.272611410757,
coefficient=0.0370952793157,
ledger_total=0.332964363589,
threshold_slack=0.000368969744373,
next_admissible=0.0426269309291,
sparse_ratio=0.0136639102342
```

Complete \(q=4\) fixed-three quintic enumeration and a direct \(q=8\)
maximizing row protect the endpoint energy. Exact random-law and saturating
checks protect the cubic XOR-column compression, while sparse full occurrence
tensors protect the indexing and normalization. The \(q=32\) regression
protects the safe coefficient, seven-theorem ledger, reranking, and next gate.

This gives a rigorous coefficient for four more balanced entries. Together
with the first six theorems, twenty-eight of 888 entries are controlled. The
remaining 860 provisional coefficients, physical diagnostic charges, and
adaptive lift keep the result from being a complete passive-dose-six
certificate. The run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## Sixth balanced whole-row contraction run

The complete suite passed again on 2026-07-15 after adding the sixth
chain-aware balanced-orbit theorem. The new regression reported:

```text
adjacent balanced row-slice contraction passed:
q4_energy=0.376302083333,
q8_horizontal=0.04279649918,
q32_record_one=0.000177699395695,
q32_record_three=0.00160660282258,
coefficient=0.0422410016249,
ledger_total=0.332775779206,
threshold_slack=0.000557554127147,
next_admissible=0.0484819899411,
q4_sparse_ratio=0.00698741756184
```

Complete \(q=4\) enumeration checks all translated whole rows, and selected
\(q=8\) rows protect the exact horizontal Walsh tail. Support-count
regressions reproduce the record-one and record-three incidence formulas;
sparse exact target tensors under correlated diagonal laws protect the Schur
normalization. The \(q=32\) regression protects the safe coefficient,
six-theorem ledger, reranking, and next gate.

This gives a rigorous coefficient for four more balanced entries. Together
with the first five theorems, twenty-four of 888 entries are controlled. The
remaining 864 provisional coefficients, physical diagnostic charges, and
adaptive lift keep the result from being a complete passive-dose-six
certificate. The run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

Bytecode generation is disabled in both runners so the check leaves no
`__pycache__` cleanup artifacts.

## Extended Round 3 run

After resolving OP3-1 and adding the first OP3-2 falsification check, the
same command again passed the complete inherited suite and the two new
Round 3 regressions:

```text
alternating double-endpoint spectrum passed:
q32=0.00035759351713982535053208591764918627335980824852312,
classes=20,rank=1048576

weighted alternating q=2 regression passed:
uniform=0.471591815891,
optimum=0.4718447392358088635,
upper=0.471844739237
```

The first line is an exact integer spectrum certificate before the displayed
decimal is evaluated.  The second is explicitly a certified numerical
fixed-split bound with a conservative \(10^{-12}\) allowance.

The suite was run again after the exact same-middle contraction, refined
mixed coefficient, whole-block cut correction, and complete occupation
optimization were added.  It additionally reported:

```text
weighted double-endpoint contraction passed:
same_q32=0.0109048201874,
mixed_q32=0.020343,
updated_attenuated_ledger=0.0632672427541,
remaining_slack=0.0970908892039

same-middle weighted contraction passed:
equal_q32=0.000947029359879...,
mixed_q32=0.001809666065927...,
hybrid_q32=0.00544009811458...,
attenuated_ledger=0.0825145929133,
slack=0.0778435390447

double-endpoint occupation optimization passed:
objective=0.498692859697,
supporting_upper=0.498694038261,
attenuated_upper=0.115980294719,
margin_slack=0.0443778372385
```

The last regression covers all 64 cuts and all 210 dose-six occupation
states.  It protects the isolated double-endpoint gate, not the complete
degree-twelve or adaptive theorem; the remaining profiles cannot each spend
the same \(0.160358\) margin.

The first OP3-3 degree-six extension is also protected as a negative result:

```text
degree-six joint occupation barrier confirmed:
supporting_upper=0.349674717182,
degree_six=0.245633147213,
double_endpoint=0.104040021372,
overshoot=0.189316585224
```

This run uses the exact endpoint-cubic fixed-mask coefficients
\(0.000324856644\ldots\) and \(0.010395412614\ldots\), the middle-cubic
coefficient \(1/961\), and corrected whole-block cut tables.  The two
endpoint-cubic profiles contribute more than 99 percent of the remaining
degree-six cost.  The next required gain is therefore a compound contraction
that preserves their common physical diagonal law.

The suite was run again after the exact occupation-pairing refinement,
degree-eight incidence ledger, and subgaussian promise theorem were added.
It additionally reported:

```text
attenuation promise concentration passed:
beta=0.7825,
mean=0.374918943789,
sign_proxy=0.743972172822,
chain_proxy=0.00161057653836,
conditioning_loss=0.0157437572845

occupation-compatible sector ledger passed:
degree8_upper=0.281922500044,
promise_loss=0.0157437572845,
partial_total=0.297666257329,
partial_slack=0.0356670760046,
target_1_over_112=0.329941349153,
target_1_over_96=0.334606949048
```

The first regression uses the Kearns--Saul moment-generating-function bound,
not a Gaussian approximation.  The second includes every representation
sector through total degree eight and reduces the occupation maximum to a
210-state Perron eigenvalue.  It is a partial one-batch certificate: degrees
ten and twelve and the adaptive lift remain open.

The suite was run again after the chained-slice contraction for the four
leading triple-cubic profiles.  The occupation regression additionally
reported:

```text
known_high_total=0.324831904411,
known_high_slack=0.00850142892273,
target_1_over_240=0.333045380316,
target_1_over_224=0.333637103977
```

Here `known_high_total` includes all sectors through degree eight, the
already-proved high endpoint profiles, all four triple-cubic profiles, and
the promise loss at \(\beta=25/32\).  It still omits the other high-degree
profiles and is not an adaptive certificate.

The suite was run again after proving the degree-five endpoint slices and
both separated cubic--quintic reversal pairs.  The occupation regression at
\(\beta=781/1000\) now reports:

```text
known_high_total=0.33081366446,
known_high_slack=0.00251966887305,
target_1_over_588=0.333321204449,
target_1_over_584=0.333338429976
```

The same regression reconstructs the full degree-five signed-permutation
moments at \(q=4\), checks every middle-quintic fixed-pair/triple isomorphism
through \(q=8\), and protects the sharp four-cycle saturator.  The accepted
partial ledger now includes both separated pairs.  The opposite-endpoint and
adjacent cubic--quintic contractions and the adaptive lift remain open.

## Repaired open-profile routing run

The complete suite passed again after inserting the physical forced-cut
baseline into the open-profile inventory and adding the leading unforced
physical witness.  The new regressions reported:

```text
transposed opposite-endpoint witness passed:
q4=0.0554262492783,
q8=0.0825349366921,
q4_direct=0.0554262492783,
q8_block_size=392

repaired open-profile budget passed:
open_profiles=43,
reversal_orbits=23,
baseline_slack=0.000997380096829,
common_threshold=0.000735032568437,
top_profile=((3,1,1,5),(5,1,1,3)),
top_split_threshold=0.0139022838042
```

The first regression protects one exact finite physical-law calculation, not
a \(q=32\) contraction.  The second protects a routing diagnostic, not a
complete one-batch certificate.  The run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## Certified second-witness and hybrid-repair run

The complete suite passed on 2026-07-15 after adding the certified
\(q=32\) dominant-class witness, continuing the finite-tilt promise theorem
with its hybrid global branch, and reranking the remaining profiles.  The
decision-critical regressions reported:

```text
attenuation promise concentration passed:
repaired_loss=0.0181887547825,
hybrid_loss=0.0175468358417

transposed opposite-endpoint witness passed:
q4=0.0554262492783,
q8=0.0825349366921,
q32_trace_lower=194,
q32_coefficient_lower=0.0142810242047

repaired open-profile budget passed:
open_profiles=43,
reversal_orbits=23,
certified_obstruction=0.333360945553,
hybrid_total=0.332839308989,
next_split_threshold=0.00902715239134
```

The \(q=32\) check uses an exact integer column Gram plus explicit binary64
roundoff majorants.  It certifies a lower witness from one symmetry class;
it does not extrapolate the nonmonotone \(q=4,8,16\) data.  The hybrid total
is a two-witness scalar diagnostic, not the completed high-degree one-batch
ledger or the adaptive theorem.

The run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## Adjacent-record gate and broadened-mission run

The complete suite passed again on 2026-07-15 after reconstructing the
leading adjacent cubic--quintic tensor, proving its record-one and
record-three slice formulas, extending Euclidean promise packing through the
two-split gate, and aligning the Round 3 program documents around the
resulting compound-frame target.  The decision-critical regressions reported:

```text
adjacent cubic-quintic witness passed:
q4_orbit=0.0295138888889,
q8_orbit=0.00171977970876,
q32_record_one=0.0118822088518,
q32_record_three=0.0119701029103,
q32_combined=0.0168662459036

attenuation promise concentration passed:
extended_loss=0.0200666604394,
same_beta_hybrid_loss=0.0202308591789

repaired open-profile budget passed:
extended_baseline=0.332675996818,
extended_common=0.00051010716305,
extended_threshold=0.011988588293,
combined_adjacent=0.0168662459036
```

The small-orbit values are exact finite diagnostics, not a scaling
extrapolation.  The two record-sector coefficients at \(q=32\) are exact
squared-slice bounds.  Each fits the extended scalar threshold separately,
but their safe direct combination does not.  Thus this run protects a
negative result for the local row-energy contraction and the quantitative
next target: a shared-law compound factor below \(0.7108035992\).  It does
not certify the complete one-batch ledger or the adaptive passive theorem.

The run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## Chained accepted-sector repair and broad-target run

The complete suite passed again on 2026-07-15 after replacing the loose
adjacent double-cubic record-three incidence term by chained exact squared
slices, validating the adjacent mixed-orbit Fourier reduction, and updating
the program documents around the broader \(1/q\)-classification target.  The
new decision-critical regressions reported:

```text
occupation-compatible sector ledger passed:
degree8_upper=0.232480992267,
partial_total=0.248224749551,
known_high_total=0.279758546919,
target_1_over_32=0.326476362740,
target_1_over_24=0.343139196520

adjacent mixed-orbit diagnostics passed:
pure=0.0295138888889,
mixed=0.058361589908,
mixed_direct=0.058361589908,
symmetric_q4=0.185024058902

repaired open-profile budget passed:
forced_entries=12,
coarse_entries=6004,
coarse_beta=0.779698447178,
coarse_total=0.322669154028,
coarse_slack=0.0106641793051
```

The chained record-three coefficient is a proved fixed-split upper bound.
The all-open \(1/q\) ledger is deliberately only a theorem target: its open
coefficients and physical diagnostic charges still require arbitrary-law
upper certificates.  The run does not prove the complete one-batch or
adaptive passive lower bound.

## High-degree shell-audit run

The complete suite passed again on 2026-07-15 after the exact high-degree
record/incidence inventory and dose-six shell audit were added.  The new
regression reported:

```text
high-degree record/incidence frontier passed:
profiles=43,
record_sectors=14624,
certified=1138,
dose_six_relevant=888,
dose_six_certified=0,
leading_contribution=0.00394518332418
```

The incidence dynamic program is checked against exhaustive \(q=4\)
enumeration.  The 1,138 coefficients are rigorous generic upper bounds, but
all lie outside the balanced occurrence shell and therefore contribute zero
to the dose-six ledger.  The leading contribution is a Perron
route-selection diagnostic, not an arbitrary-law coefficient theorem.  The
run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## Corrected Track B exponent-audit run

The complete suite passed again on 2026-07-15 after extracting the exact
graph-dependent suppression from the repaired reverse-tree proof, correcting
the historical theorem ladder, and adding the legal layered interface
saturator.  The new regression reported:

```text
general-theory exponent audit passed:
levels=4-12,
accepted=1/24,
stale_uniform=1/16,
corrected_uniform=1/8,
second_suppression=1/16,
interface_sigma=1
```

The \(1/24\) row and exponent algebra are exact.  The checked
\(\sigma=1\) family is a limitation of the current legal graph/placement
interface; it is not asserted to be a physical posterior witness or a
nonzero terminal interpolation diagram.  The \(1/16\) and \(1/8\) rows
remain theorem targets.  The run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## First balanced chain-aware contraction run

The complete suite passed again on 2026-07-15 after proving the first
arbitrary-law upper bound on the dose-six balanced frontier.  The new
regression reported:

```text
leading balanced disjointness contraction passed:
coefficient=0.0934752745775,
ledger_total=0.326281860918,
threshold_slack=0.00705147241501,
worst_random_ratio=0.382810229268
```

The exact centered factorization of the \(k\)-set/singleton disjointness
mask is checked at several small dimensions.  The full abstract chain tensor
is also constructed under random unit-feature Gram laws and arbitrary
correlated diagonal weights.  The \(N=1024\) coefficient and ledger
optimization are deterministic.

This is a rigorous fixed-split upper for four balanced entries.  The ledger
total still retains provisional \(1/q\) charges and physical lower witnesses
elsewhere, so it is not the complete one-batch certificate or the adaptive
passive lower bound.  The run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## Second balanced cubic-slice contraction run

The complete suite passed again on 2026-07-15 after proving the second
arbitrary-law upper bound on the dose-six balanced frontier.  The new
regression reported:

```text
adjacent balanced cubic-slice contraction passed:
cubic_slice=0.0303049395161,
coefficient=0.0162724692796,
generic_two_mask=0.225293047398,
ledger_total=0.32556385797,
threshold_slack=0.00776947536351,
worst_random_ratio=0.692239897219
```

Complete \(q=4\) signed-permutation enumeration reproduces the exact cubic
fixed-pair slice.  The abstract test constructs the full adjacent chain
under random unit-feature Gram laws and arbitrary correlated diagonal
weights.  It also records that the generic two-mask relaxation is too large
for the current ledger.

This is a rigorous fixed-split upper for four more balanced entries.  Together
with the first chain-aware theorem, eight of 888 balanced entries are
controlled.  The optimized total still retains provisional \(1/q\) charges
and physical lower witnesses elsewhere, so it is not the complete one-batch
certificate or the adaptive passive lower bound.  The run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## Two-complete-flag active-boundary run

The complete suite passed again on 2026-07-15 after adding the exact Track C
obstruction to collective decoding of two complete active flags.  The new
regression reported:

```text
active two-flag collective obstruction passed:
N0=16,
N=1024,
factor_inputs=22528,
rank=100,
trace_distance=0.277677889243,
helstrom_error=0.361161055379,
error_margin=0.0278277220453
```

The regression exhausts the 22,528 \(N=4\), \(F=+1/2\) factor inputs and the
same number with negative sign.  It constructs the exact \(N=16\) endpoint
ensemble, reduces the two-copy hypothesis difference to an integer matrix,
and factors the characteristic polynomials of sixteen 32-dimensional Gram
components.  Tensor multiplicativity and the fixed-unit isometric lift to
\(N=1024\) are checked independently.

This proves an obstruction for arbitrary collective POVMs on two complete
path--mode flag states.  It does not cover a coherent extra one-traversal
query or a genuinely interleaved five-traversal circuit, so it is not a
general active hard-dose-five lower bound.  The run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## Third balanced endpoint-slice run

The complete suite passed again on 2026-07-15 after adding the third
chain-aware balanced-orbit theorem.  The new regression reported:

```text
separated balanced endpoint-slice contraction passed:
cubic_slice=0.0303049395161,
quintic_slice=0.99609375,
coefficient=0.173742800847,
ledger_total=0.332665118954,
threshold_slack=0.000668214378968,
next_contribution=0.00150951796807,
next_admissible=0.0450405467778,
worst_random_ratio=0.579698093621
```

Complete \(q=4\) enumeration independently reproduces the cubic fixed-pair
slice and the extracted quintic fixed-four scalar energy.  In particular it
checks the normalization \(N(1-4/N)\), which rules out a tempting but
incorrect extra factor of \(1/q\).  The abstract regression constructs the
full separated chain under arbitrary diagonal laws.

This gives a rigorous coefficient for four more balanced entries.  Together
with the first two theorems, twelve of 888 entries are controlled.  The
reoptimized diagnostic remains below \(1/3\), but only by
\(0.000668214379\); the remaining provisional coefficients and physical
diagnostic charges keep it from being a complete certificate.  The run ended
with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## Fourth balanced shared-law contraction run

The complete suite passed again on 2026-07-15 after adding the fourth
chain-aware balanced-orbit theorem.  The new regression reported:

```text
internal singleton shared-law contraction passed:
q4_cubic=0.693158531651,
q4_collision_energy=1.91666666667,
q4_random_ratio=0.809441504056,
q32_cubic=0.332653203639,
coefficient=0.0250967461185,
simple_overshoot=0.00385958576689,
ledger_total=0.33236830015,
threshold_slack=0.00096503318286,
next_admissible=0.054250629776,
vertical_mixture=0.00665619197928
```

Complete \(q=4\) enumeration checks the affine-twirled cubic factor, random
nonsymmetric laws, and the quintic completion/collision decomposition.  The
\(q=32\) regression protects the exact Walsh spectrum, safe arbitrary-law
coefficient, reoptimized ledger, and next gate.  The symmetric physical-law
mixture is a lower-witness stress test, not an upper-bound theorem.

This gives a rigorous coefficient for four more balanced entries.  Together
with the first three theorems, sixteen of 888 entries are controlled.  The
remaining provisional coefficients, physical diagnostic charges, and
adaptive lift keep the result from being a complete passive-dose-six
certificate.  The run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```

## Fifth balanced quintic-row contraction run

The complete suite passed again on 2026-07-15 after adding the fifth
chain-aware balanced-orbit theorem.  The new regression reported:

```text
column-cubic quintic-row contraction passed:
q4_row_energy=0.75,
q4_base_rank=16,
q4_base_attained=0.25,
q4_sparse_ratio=0.0159582921714,
coefficient=0.0311889051224,
ledger_total=0.332365794098,
threshold_slack=0.0009675392353,
next_admissible=0.0570749885142
```

Complete \(q=4\) enumeration checks the quintic fixed-four row energy and
the compressed cubic--Hadamard base matrix.  It verifies rank \(N\), maximum
entry \(1/N\), and a law attaining the base \(1/\sqrt N\) factor, then
stress-tests exact target submatrices under correlated diagonal laws.  The
\(q=32\) regression protects the safe coefficient, five-theorem ledger,
reranking, and next gate.

This gives a rigorous coefficient for four more balanced entries.  Together
with the first four theorems, twenty of 888 entries are controlled.  The
remaining provisional coefficients, physical diagnostic charges, and
adaptive lift keep the result from being a complete passive-dose-six
certificate.  The run ended with:

```text
PASS round-two reverse-tree stress checks
PASS round-three initialization and inherited baseline
```
