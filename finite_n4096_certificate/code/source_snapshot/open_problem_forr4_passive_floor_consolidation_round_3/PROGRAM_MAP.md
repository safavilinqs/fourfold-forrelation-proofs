# Program map

Date: 2026-07-15

## The enduring question

How much more experimental exposure does passive sensing need than active
sensing to solve four-forrelation, and what feature of the models causes the
difference?

This question has three layers:

1. **Existence:** is there a rigorous active--passive separation at all?
2. **Strength:** how quickly must passive dose grow with \(N\)?
3. **Relevance and mechanism:** does the gap appear at realistic \(N\), is six
   the right active benchmark, and which structural resource creates it?

The original folder's minimum success criterion was existence: prove any
polynomial passive floor for the exact constant-margin problem, while keeping
the active dose at most six, or else find a subpolynomial passive
counterprotocol.  That criterion has been met.  The remaining program is the
natural second half of the original question: determine how strong, useful,
and close to optimal that separation actually is.

## What each folder contributed

### Original exploration

`../open_problem_forr4_passive_floor`

Developed the model, active protocols, candidate hard distributions,
interpolation and reverse-tree ideas, exact low-dose results, and several
possible routes to a passive floor.  Its role was discovery.

### First consolidation

`../open_problem_forr4_passive_floor_consolidation`

Audited, repaired, and compressed the proof into the established separation

$$
D_{\mathsf P}^{\rm hard}=\Omega(N^{1/24}),
\qquad
D_{\mathsf A}^{\rm hard}\le6.
$$

Its role was to establish that the separation is real and referee-checkable.

### Round two

`../open_problem_forr4_passive_floor_consolidation_round_2`

Stress-tested the reverse-tree contraction, found and repaired four necessary
issues, and showed why the accepted asymptotic theorem does not yet give a
useful (D=6,N=1024) statement.  It then developed the signed-permutation
finite-size route and isolated its weighted compound-frame bottleneck.

Its role was adversarial validation and quantitative diagnosis.

### Round three

This folder.

Round three should determine the strongest meaningful next separation, not
merely finish one inherited calculation.  It combines:

- a focused attempt at passive \(>6\) near \(N=1024\);
- a general program that has improved the inherited \(N^{1/24}\) floor to
  \(N^{1/12}\), closed the complete terminal obstruction list, and now
  targets physical frame structure, unusable constants, and replacement
  witnesses after proving the generic level-twelve graph norm sharp; and
- active tests of alternative instances, passive counterprotocols, and the
  active-dose boundary.

Its role is strengthening, finite-size relevance, and mechanism discovery.
The three live projects and their stop rules are maintained in
`PORTFOLIO.md`; no single hard instance defines the round.

The broadened charter also restores a part of the original program that can
be obscured by proof work: a separation is a comparison of two true
complexities, not just one lower bound and one inherited construction.
Round 3 therefore maintains passive protocol evidence and the active
lower/upper frontier alongside passive lower bounds.  `BOUNDARY_MAP.md`
records this four-sided comparison.

The first complete high-level terminal audit is now a genuine theorem-level
advance.  Exact Stein-state enumeration leaves eight sensitive types at
levels 9--12.  The better of the assigned and all-projective contractions
gives \(N^{-1}\) for every level-11/12 type, while exactly three level-9/10
types retain only \(N^{-1/2}\).  Consequently the general passive floor is
first \(\Omega(N^{1/20})\).  A centered-weight re-expansion now gives
\(N^{-1}\) for the unique level-ten type as well, improving the floor to
\(\Omega(N^{1/18})\).  A second centered endpoint audit now repairs both
reflected level-nine trees and proves \(\Omega(N^{1/16})\).  The high-level
terminal list is closed.  Extending the image to both potential-eight
initial collisions and the potential-four path then shows that every
sensitive level-eight type already has \(N^{-1}\) decay.  Forced centered
derivatives repair the unique level-seven and level-six saturators, proving
\(\Omega(N^{1/12})\).  The complete terminal list is now closed.  An explicit
full-distinctness lower witness matches the limiting \(N^{-1}\) grouped graph
bound, so the next asymptotic work must use physical posterior-frame
structure, preserve more signed coefficient structure before terminal
norms, or test a different hard instance with better scaling or constants.

The first broad finite-size classification illustrates this operating model.
Generic record/incidence estimates rigorously settle 1,138 open entries, but
an exact shell audit shows that all are irrelevant to the dose-six ledger.
Rather than extending that locally successful but mission-irrelevant route,
Round 3 now redirects Track A to 888 balanced entries, schedules one bounded
Track B or C deliverable, and retains an explicit hard-instance pivot.

The first balanced entry class has now yielded to a genuinely chain-aware
argument.  Four complement/reversal cuts receive a rigorous arbitrary-law
coefficient \(0.0934752746\); after paying this conservative value the
finite-size diagnostic still has \(0.0070514724\) slack.  This is the intended
Round 3 workflow: turn a local theorem into a program-level ledger decision,
then move to the next unresolved class rather than treating the local lemma
as the mission.

The second balanced class strengthens that lesson.  Its generic two-mask
bound would fail the current finite-size budget, but retaining the exact
fixed-pair cubic energy proves coefficient \(0.0162724693<1/32\) for four
more cuts.  With both theorems inserted the diagnostic has \(0.0077694754\)
slack, and eight of 888 balanced entries are controlled.  This is evidence
for chain-aware contraction as a mechanism, while the remaining 880 entries,
physical exceptions, and adaptivity keep the broader mission open.

The third balanced class supplies the first sharp warning against extrapolating
from those successes.  Exact endpoint slices prove coefficient
\(0.1737428008\) for \((3,1,1,5):(1,0,0,4)\), better than the generic
two-mask bound but much worse than the provisional \(1/32\).  The resulting
diagnostic has only \(0.0006682144\) slack, with twelve of 888 entries
controlled.  Round 3 therefore changes posture: the next reranked orbit is a
finite-size route decision, while Tracks B and C remain independent ways to
advance the original program if this hard instance cannot close.

That route decision has now passed for structural reasons.  On
\((3,1,1,5):(1,1,1,2)\), separate endpoint maxima would fail, but preserving
the common diagonal law permits an exact cubic twirl and a completed quintic
Schur contraction.  The coefficient \(0.0250967461\) controls four more cuts;
the updated diagnostic is \(0.3323683002\), with \(0.0009650332\) slack.
Sixteen of 888 entries are controlled.  This strengthens the evidence for
shared-law chain structure without narrowing Round 3's mission: the next
orbit is another bounded Track A gate, while asymptotic improvement and
boundary tests remain coequal program outcomes.

The fifth balanced class passes by a complementary mechanism.  For
\((3,1,1,5):(0,0,1,4)\), exact quintic row energy avoids a generic factor
near three, and rank--Frobenius contracts the remaining chain at
\(1/\sqrt N\).  The coefficient \(0.0311889051\) controls four more cuts;
the updated diagnostic is \(0.3323657941\), with \(0.0009675392\) slack.
Twenty of 888 entries are controlled.  This local success keeps Track A
alive but leaves the finite-size theorem, adaptive lift, asymptotic
improvement, and boundary program open.

The sixth balanced class confirms that the useful object can be an entire
chain row rather than a separate endpoint slice.  For
\((1,1,3,5):(0,1,1,3)\), L-shape incidence controls record one and exact
Walsh orthogonality controls the leading horizontal record-three tail.  The
coefficient \(0.0422410016\) passes its gate.  The updated diagnostic is
\(0.3327757792\), with \(0.0005575541\) slack; twenty-four of 888 entries are
controlled.  This preserves Track A but narrows its margin and strengthens
the case for whole-chain mechanisms that may also inform Track B.

The seventh balanced class shows that a full adjacent-row treatment is not
always necessary.  For \((3,1,1,5):(0,1,1,3)\), the exact fixed-three
quintic energy is extracted as a Schur feature, while the opposite whole
cubic endpoint compresses by support XOR to a weighted Walsh matrix.  The
coefficient \(0.0370952793\) passes its gate.  The updated diagnostic is
\(0.3329643636\), with \(0.0003689697\) slack; twenty-eight of 888 entries are
controlled.  Track A remains alive but increasingly margin-limited.

The eighth balanced class shows that even a pointwise middle-link maximum
can pass when it is delayed until after exact endpoint sums are retained.
For \((1,3,5,1):(0,2,2,1)\), chaining the cubic and quintic fixed-pair
energies gives coefficient \(0.0285281523<1/32\).  The updated diagnostic is
\(0.3328775891\), with \(0.0004557442\) slack; thirty-two of 888 entries are
controlled.  This recovers some margin but does not change the need for the
parallel Track B and C programs.

The ninth balanced class preserves more of the remaining chain.  For
\((1,1,3,5):(0,0,3,2)\), a normalized fixed-pair \(M_{35}\) row leaves the
XOR-labelled \(H_NM_{13}\) matrix intact, and duplicate compression supplies
an additional \(1/q\).  Exact record-pattern counts give coefficient
\(0.0250919472<1/32\).  The updated diagnostic is \(0.3326862124\), with
\(0.0006471209\) slack; thirty-six of 888 entries are controlled.  This
reusable normalize-then-compress mechanism keeps Track A alive without
changing the coequal status of Tracks B and C.

The tenth balanced class uses both endpoint Walsh factors at once.  For
\((1,3,5,1):(0,2,3,0)\), extracting the intervening scalar completion row
leaves a repeated \(H_N\otimes H_N\) matrix with coefficient one.  The
record-sector row estimate gives \(0.0462425962<0.0529177167\), but this is
larger than the provisional \(1/32\).  The updated diagnostic is
\(0.3331326055\), with only \(0.0002007278\) slack; forty of 888 entries are
controlled.  This is a useful proof replacement and a warning against
letting local ledger work displace the broader mission.  The next orbit is a
bounded \(0.0379251204\) gate test, while Tracks B and C retain equal claim
on the next substantial block of effort.

The scheduled Track B audit has also been completed.  It corrects a stale
research rung: the accepted proof already pays \(D^v\), so a uniform
\(N^{-v/8}\) contraction would prove \(N^{1/8}\), not \(N^{1/16}\).
The exact existing graph factor is \(N^{-\sigma/2}\), but the theorem
interface permits \(\sigma=1\) through level twelve.  The next general
theorem must therefore use high-level interpolation or posterior physical
structure that the graph invariants discard.

The first high-level image test rules out one tempting shortcut.  The
original relaxed star is absent, but a genuine sensitive level-twelve
diagram of three disjoint four-layer paths is present with positive local
coefficient.  One paired physical entry gives
\(\sigma_{\rm ass}=1\).  Therefore interpolation support alone does not
improve the asymptotic theorem.  But the physical norm test succeeds: a
single all-projective contraction gives \(N^{-1}\) on every
occupancy-\((2,1,1)\) placement, stronger than the required \(N^{-3/4}\),
without invoking the invalid RT-003 conversion.  The remaining Track B
question is now global: whether every true high-level terminal type passes
the better of the assigned and grouped-entry projective bounds.

The first Track C audit also changes the boundary map without overstating a
general lower bound.  At \(N=1024\), an exact endpoint ensemble keeps two
complete folded-chain flags at collective trace distance
\(0.2776778892<1/3\), so their optimal error is \(0.3611610554\).
Deleting the third flag or changing the decoder cannot improve the accepted
six-dose construction.  The only live five-dose continuation must use the
fifth traversal coherently or adopt a genuinely different interleaved word.

For a nontechnical statement of this broadened scope and the integrated end
product, start with `ROUND_3_CHARTER.md`; the technical application to the
current frontier is in `notes/round3_broadened_mission_guide.md`.

## The questions Round 3 inherits

| Program question | Present answer | Round 3 responsibility |
|---|---|---|
| Does an unrestricted active--passive separation exist? | Yes: passive \(\Omega(N^{1/12})\), active at most six | protect the improved proof and continue the boundary audit |
| How large is the passive dose really? | between \(\Omega(N^{1/12})\) and known general upper possibilities | test physical frame restrictions, constants, alternative witnesses, and passive protocols |
| Is the gap visible at realistic size? | not from the asymptotic constants | prove passive \(>6\) near \(N=1024\), or explain quantitatively why a route cannot |
| Is active dose six the right benchmark? | six is an upper bound, not an optimum; two complete flags fail at \(N=1024\) even with arbitrary collective decoding | test a coherent extra query or genuinely interleaved five-traversal circuit |
| What causes the advantage? | coherent active traversal is useful; the sharp passive obstruction is unsettled | isolate a posterior-stable contraction or a counterexample to one |
| Is the current plant the right witness? | promising but not privileged | compare it against alternatives using common gates |

`MISSION_LEDGER.md` turns these questions into maintained evidence records and
end-of-round deliverables.

## How to judge a Round 3 project

A project belongs in this round when it does at least one of the following:

- proves the separation at a realistic size;
- improves the asymptotic passive floor;
- extracts a general principle explaining passive contraction;
- rules out a serious route under the physical model;
- identifies a better hard instance;
- finds a passive counterprotocol that changes the conjectured boundary; or
- sharpens the active side of the comparison.

A calculation that does none of these can still be useful, but it should not
control the program merely because it is the current calculation.
