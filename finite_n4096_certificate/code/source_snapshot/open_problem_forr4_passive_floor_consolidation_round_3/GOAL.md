# Round 3 goal

Date: 2026-07-15

## Plain-language goal

Find out how much more probing passive sensing really needs than active
sensing for four-forrelation.  Prove the difference rigorously, make it
visible at realistic problem sizes, and do not depend on a single candidate
distribution or proof trick.

In one sentence: **determine the size, practical relevance, and cause of the
active--passive gap.**

This requires looking in both directions.  Round 3 should improve passive
lower bounds, but it should also search for passive protocols, test whether
the active cost can fall below six, and compare hard instances.  Otherwise a
successful lower-bound calculation could still leave the true gap largely
unknown.

## Mathematical goal

Determine the strongest trustworthy and quantitatively meaningful
active--passive hard-dose separation available for four-forrelation.

The inherited benchmark is

$$
D_{\mathsf P}^{\rm hard}=\Omega(N^{1/24}),
\qquad
D_{\mathsf A}^{\rm hard}\le6.
$$

Round 3 has now improved the passive exponent:

$$
D_{\mathsf P}^{\rm hard}=\Omega(N^{1/12}).
$$

The complete terminal obstruction list is now closed.  After the high-level
repairs, a four-initial-state audit shows that level eight already has
\(N^{-1}\) decay and repairs the sole level-seven and level-six saturators.
Every sensitive level from five through twelve now has a one-power bound;
level twelve is limiting.  The arbitrary grouped-entry norm at that level is
now known to be exponent-sharp: an explicit same-layer-distinct lower witness
is \(N^{-1}(1-N^{-1})^2(1-2N^{-1})^3\).  A further exponent improvement must
use additional physical frame restrictions, exploit signed coefficients
before terminal norms, or change the hard instance.  The bare scale
\(1024^{1/12}\approx1.782\) is still far below six, so this sharper theorem
does not settle the central realistic-size benchmark.

Round three should advance that benchmark in at least one of three ways:

1. certify passive hard dose greater than six at a realistic size, preferably
   \(N=1024\);
2. improve the general passive exponent or prove a stronger reusable
   contraction theorem; or
3. obtain a decisive structural result about the boundary, such as a better
   hard instance, a passive counterprotocol, or a sharper active-dose result.

The round is therefore a characterization round.  Closing the current
signed-permutation ledger would be an important result, but it would answer
only the finite-size part of the mission.  Conversely, a stronger asymptotic
theorem or a boundary-changing protocol can be a primary Round 3 result even
if that ledger is ultimately abandoned.

The central realistic-size benchmark is \(N=1024,D=6\), not the whole
definition of relevance.  A successful certificate should also be evaluated
at nearby tractable powers of two when possible, and should be compared with
the best valid passive and active protocol evidence rather than reported in
isolation.

## Why this is the next round

The original folder explored whether an active--passive separation existed.
The first consolidation audited and strengthened the asymptotic theorem.
Round two stress-tested the reverse-tree proof and repaired the defects it
found.  The asymptotic separation is therefore the foundation, not the sole
remaining problem.

Round three turns that foundation into a broader research program: make the
separation operational at finite \(N\), strengthen its asymptotic explanation,
and test how close the known bounds are to the true boundary.

## Research tracks

### Track A — realistic-size separation

Try to prove that every classically adaptive passive protocol of hard dose at
most six has transcript total variation below \(1/3\) for an explicit hard
pair near \(N=1024\).  The attenuated signed-permutation plant is the first
route.  A new chained-slice theorem reduces the accepted degree-eight cost by
about \(0.0494\), so a coarse diagnostic that assigns \(1/q\) to every
otherwise open high-degree split now fits below threshold by \(0.010664179\),
even after preserving the three known physical-orbit charges.  The live task
is to control the balanced high-degree chains that actually enter the
dose-six occupation ledger, certify the complete 210-state one-batch bound,
and only then prove the separate adaptive lift.  A first generic
record/incidence classification proves \(1/q\) for 1,138 open splits, but
all are off the dose-six shell; none of the 888 relevant balanced splits is
certified by that generic theorem.  Ten chain-aware results now control
forty of those entries.  The leading \((3,1,1,5):(0,1,0,4)\) orbit has
arbitrary-law coefficient at most \(0.093475275\).  The next
\((1,1,3,5):(0,0,1,4)\) orbit admits a sharper fixed-pair cubic-slice
contraction with coefficient \(0.0162724693<1/32\).  The separated
\((3,1,1,5):(1,0,0,4)\) orbit has a safe endpoint-slice coefficient
\(0.173742801\), substantially above the provisional \(1/32\).  For
\((3,1,1,5):(1,1,1,2)\), separate worst-case endpoint slices fail the live
gate, but a shared-law cubic twirl and a completed quintic Schur factor prove
the much smaller coefficient \(0.0250967461\).  Inserting all four rigorous
values leaves the coarse ledger below \(1/3\) by \(0.000965033\).  The next
\((3,1,1,5):(0,0,1,4)\) orbit admits an exact quintic-row-energy contraction
with coefficient \(0.0311889051<1/32\).  Inserting all five values leaves the
coarse ledger below \(1/3\) by \(0.000967539\).  The sixth
\((1,1,3,5):(0,1,1,3)\) orbit has a whole-row record-slice coefficient
\(0.0422410016\), below its \(0.0570749885\) gate.  Inserting all six values
leaves the ledger below \(1/3\) by \(0.000557554\).  The seventh
\((3,1,1,5):(0,1,1,3)\) orbit factors its fixed-three quintic row and
compresses the whole cubic endpoint by XOR-labelled Walsh columns, giving
coefficient \(0.0370952793<0.0484819899\).  Inserting all seven values leaves
the ledger below \(1/3\) by \(0.000368970\).  The eighth
\((1,3,5,1):(0,2,2,1)\) orbit chains the exact endpoint fixed-pair energies
through the universal middle-link maximum, giving
\(0.0285281523<1/32\).  Inserting all eight values leaves the ledger below
\(1/3\) by \(0.000455744\).  The ninth \((1,1,3,5):(0,0,3,2)\) orbit
extracts the complete middle-quintic row while preserving a separate
Walsh compression, giving \(0.0250919472<1/32\).  Inserting all nine values
leaves the ledger below \(1/3\) by \(0.000647121\).  The tenth
\((1,3,5,1):(0,2,3,0)\) orbit preserves both endpoint Walsh factors and
bounds the intervening scalar completion row, giving
\(0.0462425962<0.0529177167\).  This is larger than the provisional
\(1/32\), so the ten-theorem diagnostic rises to \(0.3331326055\), only
\(0.000200728\) below threshold.  The results do not close the remaining
848 balanced entries, certify the physical exceptions, or address
adaptivity.  Further generic incidence refinement remains stopped.  The
next reranked orbit is \((1,1,5,3):(0,1,3,1)\), with reoptimized acceptance
gate \(0.0379251204\).  That narrow gate makes broader Track B and C work a
coequal priority rather than an optional follow-up.

### Track B — stronger general theory

Revisit the losses in the repaired reverse-tree contraction.  Priority
targets include retaining component-dependent suppression and finding a
posterior-stable Bessel/Carleson or tester inequality that applies beyond the
current plant.  The dose-bookkeeping loss has already been removed: a
uniform \(N^{-v/8}\) theorem would now imply \(\Omega(N^{1/8})\), not
\(\Omega(N^{1/16})\).  The first corrected target is extra suppression only
at interpolation levels nine through twelve.  The terminal-image exclusion
test is now resolved negatively: the old relaxed star is absent, but a true
reflection-sensitive level-twelve forest of three four-layer paths has
nonzero positive branching weight and exact assigned-fiber
\(\sigma=1\).  The subsequent physical gate passes: keeping the whole
frame skeleton projective gives \(N^{-1}\) on every occupancy-pattern
\((2,1,1)\) placement of that forest, stronger than the required
\(N^{-3/4}\), without invoking RT-003.  The next general-theory target is
to enumerate the true terminal image at levels nine through twelve under
the better of the global assigned and grouped-entry projective bounds.  The
complete enumeration first proved \(N^{1/20}\) and isolated two reflected
level-nine trees and one level-ten \(6+4\) forest.  The level-ten forest is
now removed: all 200 legal histories carry a forced odd centered derivative,
and its exact re-expansion has decay \(N^{-1}\) on all 282 dangerous and
5,295 extended physical partitions.  This first proved \(N^{1/18}\).  The
reflected level-nine pair is now removed as well: all 200 representative
histories force an outer centered derivative, and the direct two-orientation
branch audit proves \(N^{-1}\) on every retained branch.  This proves the
intermediate \(N^{1/16}\) theorem.  The complete low-level image is now
audited too: all four sensitive level-eight types already have \(N^{-1}\),
and forced centered derivatives repair the unique level-seven and level-six
saturators.  This proves the current \(N^{1/12}\) theorem.  The next
general-theory gate is a physical-frame restriction that excludes the sharp
grouped-vector optimizer, or coefficient exclusion before terminal norms,
evaluated alongside explicit constants and a replacement witness.

### Track C — boundary and alternatives

Actively look for evidence that changes the problem: passive protocols that
challenge the expected floor, alternative hard instances that close the
finite-size gap, and active protocols or lower bounds that clarify whether
six doses are necessary.  The first active-boundary audit now rules out the
simplest improvement: at \(N=1024\), two complete folded-chain flags have
optimal collective-measurement error at least \(0.361161>1/3\) on an exact
endpoint ensemble.  Thus deleting the third flag or decoding two flags more
cleverly is insufficient.  A genuinely interleaved five-traversal circuit,
including a coherent extra one-query side experiment, remains open.

## Success ladder

- **Required for a healthy round:** a decisive, reproducible verdict on the
  signed-permutation finite-size route and at least one serious result in
  Track B or C.
- **Strong result — asymptotic side achieved:** the passive floor is now
  \(\Omega(N^{1/12})\).  A rigorous passive \(>6\) certificate near
  \(N=1024\) remains open.
- **Best result:** both a realistic-size separation and a stronger general
  theorem, together with a clear account of the mechanism causing the gap.

Round 3 should end with an integrated verdict stating:

1. the strongest rigorous separation now known;
2. whether it is meaningful near \(N=1024\);
3. what is known and unknown about the true passive and active doses;
4. which mechanism is responsible for the proved gap; and
5. which route should be continued, replaced, or stopped.

It should also leave the four-sided bound map in `BOUNDARY_MAP.md` honest:
passive lower evidence, passive upper/protocol evidence, active lower
evidence, and active upper/protocol evidence.  A blank is preferable to
silently treating the absence of a protocol or lower bound as favorable
evidence.

The live bounded deliverables in all three tracks are maintained in
`PORTFOLIO.md`.  This prevents the current finite-size calculation from
silently becoming the definition of success for the round.

A rigorous negative result counts when it rules out a route under the actual
physical constraints and supplies a justified pivot.  An unproved numerical
pattern does not.

## Non-goals

- Polishing exposition without changing the theorem or research decision.
- Re-auditing settled round-two claims without a concrete new defect.
- Treating the signed-permutation plant as mandatory if it fails its gates.
- Reporting uniform-weight numerics as a theorem for physical weights.
- Calling a one-batch or nonadaptive calculation a full passive lower bound.
- Improving an exponent symbolically while ignoring unusable constants and
  the mechanism behind them.
