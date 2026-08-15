# Hard-instance scorecard

Date: 2026-07-15

Status: Track C initialization scaffold.  Use this table before promoting a
new hard instance to the lead finite-size route.

## Common gates

Every candidate must be evaluated on:

1. constant promise separation and explicit conditioning cost;
2. complete one-batch dose-six tester norm at the target \(N\);
3. stability under adaptive outcome selection;
4. exact, analytic, or certifiable finite-size structure;
5. plausible asymptotic passive scaling; and
6. a quantitative advantage over the route it would replace.

## Initial comparison

| Candidate | Promise gate | \(D=6,N=1024\) one-batch gate | Adaptive gate | Asymptotic evidence | Current decision |
|---|---|---|---|---|---|
| Repaired interpolation/reverse-tree plant | established | constants unusable at this scale | established asymptotically | \(\Omega(N^{1/24})\) proved | asymptotic baseline |
| Attenuated signed-permutation plant | exact attenuation theorem plus finite-tilt, hybrid, and extended promise packing | accepted degree-eight conditioned total \(0.248224750\); ten balanced arbitrary-law orbit bounds through \(0.0462425962\); updated coarse completion diagnostic \(0.333132605\), with \(0.000200728\) slack and not yet a theorem | open | potentially stronger, unproved | screen the next ranked orbit against its \(0.0379251204\) gate; in parallel compare a global contraction or replacement hard instance before adaptivity |
| Quadratic-bent exact plant | exact \(F_{4,H}=\pm1\), so no conditioning loss | incomplete; weighted \((5,1)\) link is at most \(2/(N-2)=1/511\), but the complete degree-12 ledger is open | open | larger orbit may cancel decorations; unproved | named pivot candidate if the complete signed-permutation ledger fails |

For the signed-permutation row, the exact mixed-orbit witness closed the old
apparent margin: charging only its four forced cuts gave \(0.334183455\).
The finite-tilt promise theorem repairs that diagnostic to
\(0.332335953\).  A second certified physical witness then exceeds the next
endpoint-tilt gate; hybrid exact/global packing repairs both to
\(0.332839309\).  The two-split Euclidean extension improved the historical
baseline to \(0.332675997\).  A subsequent accepted-sector audit is much
stronger: chained record-three slices lower the degree-eight conditioned
total to \(0.248224750\).  Preserving the three known physical-orbit
diagnostics and assigning \(1/q\) everywhere else gives \(0.322669154\).
The finite-size lead is now a rigorous classification of those coarse
contractions, not the historical adjacent compound bound.

The first balanced classification step is now rigorous.  The four cuts in
the \((3,1,1,5):(0,1,0,4)\) orbit have arbitrary-law coefficient at most
\(0.0934752746\).  This exceeds the provisional \(1/q\), but its exact
insertion still gives \(0.3262818609<1/3\).  The next orbit,
\((1,1,3,5):(0,0,1,4)\), is now bounded at
\(0.0162724693<1/32\) using its exact fixed-pair cubic slice.  Both
insertions give \(0.3255638580<1/3\).  The next unresolved comparison is
\((3,1,1,5):(1,0,0,4)\).  That orbit is now bounded at
\(0.1737428008\): better than its generic two-mask relaxation, but much worse
than the provisional \(1/32\).  All three insertions give
\(0.3326651190<1/3\), leaving only \(0.0006682144\).  The reranked
\((3,1,1,5):(1,1,1,2)\) comparison then passes its coefficient gate:
preserving the common endpoint law proves \(0.0250967461\), whereas separate
endpoint maxima would fail.  All four insertions give
\(0.3323683002<1/3\), leaving \(0.0009650332\).  The next bounded comparison
is \((3,1,1,5):(0,0,1,4)\), with gate
\(0.0542506298=1.73602015/q\), before screening the quadratic-bent plant.

That comparison now passes at \(0.0311889051<1/32\) using exact quintic row
energy rather than a generic distinctness factor.  All five insertions give
\(0.3323657941<1/3\), leaving \(0.0009675392\).  The next bounded comparison
is \((1,1,3,5):(0,1,1,3)\), with gate
\(0.0570749885=1.82639963/q\).

That comparison now passes at \(0.0422410016\).  Keeping the whole adjacent
row exposes an exact horizontal record-three Walsh tail and a separate
L-shape record-one bound.  All six insertions give
\(0.3327757792<1/3\), leaving \(0.0005575541\).  The next bounded comparison
is \((3,1,1,5):(0,1,1,3)\), with gate \(0.0484819899\).

That comparison now passes at \(0.0370952793\).  The exact fixed-three
quintic row energy is extracted before the opposite whole cubic endpoint is
compressed by its XOR-labelled Walsh columns.  All seven insertions give
\(0.3329643636<1/3\), leaving \(0.0003689697\).  The next bounded comparison
is \((1,3,5,1):(0,2,2,1)\), with gate \(0.0426269309\).

That comparison now passes at \(0.0285281523<1/32\).  Exact endpoint
fixed-pair energies make the universal middle-link maximum affordable.  All
eight insertions give \(0.3328775891<1/3\), leaving \(0.0004557442\).  The
next bounded comparison is \((1,1,3,5):(0,0,3,2)\), with gate
\(0.0454321892\).

That comparison now passes at \(0.0250919472<1/32\).  Normalizing the
complete middle row without spending the residual XOR/Walsh chain provides
the decisive extra \(1/q\).  All nine insertions give
\(0.3326862124<1/3\), leaving \(0.0006471209\).  The next bounded comparison
is \((1,3,5,1):(0,2,3,0)\), with gate \(0.0529177167\).

That comparison now passes at \(0.0462425962<0.0529177167\) by extracting a
scalar completion row while retaining both endpoint Walsh factors.  Because
the proved coefficient is larger than the provisional \(1/32\), all ten
insertions give \(0.3331326055<1/3\) with only \(0.0002007278\) margin.  The
next bounded comparison is \((1,1,5,3):(0,1,3,1)\), with gate
\(0.0379251204\).  The scorecard now favors a parallel global-contraction or
quadratic-bent comparison before committing to a long sequence of further
local orbits.

## Pivot rule

A replacement earns lead status only when it has:

- the promise gate in hand;
- a better certified or structurally justified one-batch tester prospect than
  the failed current route; and
- a credible mechanism for surviving adaptive posterior selection.

Clean pointwise moments or attractive uniform-weight spectra are not enough.

## Candidate-generation lanes

Search in bounded lanes tied to a failed gate:

1. modify attenuation or conditioning while keeping the exact
   signed-permutation support;
2. preserve the chain promise but change the hidden matching structure;
3. combine the interpolation witness's adaptive stability with a
   finite-size-friendly exact plant; or
4. use a counterprotocol or saturator to design a witness that removes the
   exploited passive feature.

Each candidate gets one page using `PROJECT_BRIEF_TEMPLATE.md` and one row in
the scorecard.  It must identify which current loss it is designed to improve
before any large moment calculation begins.

## First comparison to run

When the current one-batch gate is decided, compare exactly two rows:

- the signed-permutation plant with its completed or obstructed ledger; and
- one named replacement selected from the lanes above.

For both, record the same promise loss, \(N=1024,D=6\) norm target,
adaptivity obstacle, certification method, and asymptotic prospect.  If no
replacement clears the first two gates on paper, keep Track B—not an
unbounded instance search—as the secondary route.
