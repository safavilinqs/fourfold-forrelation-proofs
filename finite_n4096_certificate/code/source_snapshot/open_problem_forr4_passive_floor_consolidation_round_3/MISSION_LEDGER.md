# Round 3 mission ledger

Date: 2026-07-15

## Purpose

This file keeps local technical work tied to the original program.  It is not
a proof notebook.  It records what question a result answers, how strong the
evidence is, and what decision follows.

The one-sentence mission is:

> Determine the size, practical relevance, and cause of the active--passive
> sensing gap for constant-margin four-forrelation.

## Current program verdict

| Question | Best accepted answer | Important unknown | Next decisive evidence |
|---|---|---|---|
| Does a separation exist? | Yes. Passive hard dose is \(\Omega(N^{1/12})\); active hard dose is at most six. | How close either bound is to optimal. | Use physical frame structure, obtain realistic constants, compare another hard instance, or find a protocol approaching the lower bound. |
| Is it visible near \(N=1024\)? | Not yet.  The accepted degree-eight conditioned total is \(0.248224750\).  Ten chain-aware theorems control forty balanced entries at coefficients \(0.093475275\), \(0.0162724693\), \(0.173742801\), \(0.0250967461\), \(0.0311889051\), \(0.0422410016\), \(0.0370952793\), \(0.0285281523\), \(0.0250919472\), and \(0.0462425962\); the resulting coarse completion diagnostic is \(0.333132605<1/3\), with \(0.000200728\) slack. | Whether chain-aware contraction controls the remaining 848 balanced entries without exhausting the margin, whether the physical witness charges admit upper certificates, and whether the adaptive lift preserves the gain. | Test the next ranked orbit against its \(0.0379251204\) gate only if a credible chain contraction is visible; in parallel, advance a general contraction, protocol search, or hard-instance comparison before further ledger expansion. |
| What is the passive scale? | Polynomial growth \(N^{1/12}\) is now proved.  The complete image audit gives \(N^{-1}\) at level eight and forced centered derivatives remove the unique level-seven/six saturators in addition to the high-level types. | The level-twelve row is sharp for arbitrary grouped graph vectors; whether physical frames or another witness improve it remains open. | Extract realistic constants and compare the present witness against a named alternative and a passive protocol; reopen the contraction only with a concrete frame constraint. |
| What is the active scale? | A hard-dose-six protocol is proved.  At \(N=1024\), two complete folded-chain flags cannot reach error \(1/3\), even with an arbitrary collective POVM: an exact endpoint ensemble gives optimal average error \(0.3611610554\). | Whether a genuinely interleaved five-traversal circuit, or two flags plus one coherent extra query, succeeds; no general active hard-dose-five lower bound is proved. | Optimize the one-extra-query/interleaved-five family against the exact endpoint ensemble. |
| How good can passive protocols be? | No competitive passive upper-bound frontier has yet been consolidated in Round 3; the initialized small-\(N\), hard-dose-six structured search is explicitly exploratory. | Whether a valid passive strategy exploits repeated modes, collective receivers, or posterior adaptation more efficiently than the lower-bound ledgers suggest. | Exhaust the bounded structured small-\(N\) protocol screen and report an actual protocol value separately from relaxation values. |
| What causes the gap? | Active coherent traversal directly measures a folded chain; passive transcripts admit dimension-decaying contractions. | The sharp reusable passive obstruction under adaptive posterior selection. | A Bessel/Carleson, tester, or component-sensitive contraction with a tested sharpness example. |
| Is the witness robust? | The repaired interpolation witness proves the asymptotic theorem; the signed-permutation witness has strong finite-size structure. | Whether either is near-optimal for realistic size and adaptivity. | Common-scorecard comparison with a named alternative hard instance. |

## Required end-of-round deliverables

### 1. Strongest theorem

State the best rigorous active--passive separation, including the physical
model, hard-dose meter, error/margin, asymptotic exponent, and usable
constants.  Keep lower bounds, upper bounds, and conjectures separate.

Current entry: passive \(\Omega(N^{1/12})\), active at most six, for the exact
constant-margin task and unrestricted classically adaptive passive trees.

### 2. Realistic-size verdict

State whether passive hard dose six has been excluded at \(N=1024\).  If not,
name the exact missing theorem or the certified reason the chosen route
cannot close.

Current entry: not established.  A rigorous chained-slice repair lowers the
accepted degree-eight conditioned total to \(0.248224750\), and the proved
high-sector partial total is \(0.279758547\).  A routing calculation that
retains all three known physical-orbit diagnostics and charges every other
open split at \(1/q\) reoptimizes to \(0.322669154\), leaving
\(0.010664179\).  The \(1/q\) charges and physical-witness substitutions are
not yet arbitrary-law upper bounds.  The finite high-sector classification,
full interval certificate, and all adaptivity remain.

New routing evidence: the open inventory contains 888 balanced entries that
can join two total-occupation-six states.  Exact record-incidence bounds prove
\(1/q\) for 1,138 other entries but none of these 888.  This is a rigorous
limitation of the generic architecture, not evidence against the desired
coefficient itself.  Decision: stop generic incidence refinement and test a
chain-aware contraction on the leading balanced orbit.

First chain-aware result: for the four-cut
\((3,1,1,5):(0,1,0,4)\) orbit, exact cubic--Walsh collapse gives \(1/q\)
after completing the quintic link, and the distinct-label \(4|1\) mask costs
the explicit factor \(2.991208786\).  The arbitrary-law coefficient is
\(0.0934752746\).  Substituting this rigorous upper for the provisional
\(1/32\) raises the diagnostic to \(0.3262818609\), leaving
\(0.0070514724\).

Second chain-aware result: for
\((1,1,3,5):(0,0,1,4)\) and its complement/reversal cuts, the exact cubic
fixed-pair energy is
\((q^2-2q+2)/(q^2(q-1))\).  Its square-root Schur factor, the completed
adjacent link, and the quintic \(4|1\) mask give arbitrary-law coefficient
\(0.0162724693<1/32\).  With both theorems inserted, the diagnostic is
\(0.3255638580\), leaving \(0.0077694754\).  Decision: continue the balanced
frontier at \((3,1,1,5):(1,0,0,4)\); do not spend the remaining margin as if
the other provisional coefficients were proved.

Third chain-aware result: for
\((3,1,1,5):(1,0,0,4)\), the cubic fixed-pair slice and extracted quintic
fixed-four slice prove coefficient \(0.1737428008\).  It beats the generic
two-mask coefficient \(0.2252930474\) but is substantially larger than the
provisional \(1/32\).  With all three theorems inserted, the diagnostic is
\(0.3326651190\), leaving only \(0.0006682144\).  Twelve of 888 balanced
entries were controlled.  At that stage, the decision was to use the reranked
\((3,1,1,5):(1,1,1,2)\) orbit as a go/no-go test; if its safe insertion
exceeds the reoptimized gate \(0.0450405468=1.44129750/q\), pivot the
finite-size hard instance rather than assuming the remaining provisional
entries can be repaired.

Fourth chain-aware result: for
\((3,1,1,5):(1,1,1,2)\), separate endpoint maxima would give
\(0.0991470258\) and fail the previous gate.  Joint concavity and affine
twirling reduce the common-law cubic factor to \(0.3326532036\), while a
complete quintic Gram multiplier plus an overlap correction costs at most
\(1+\sqrt2\).  Including the internal singleton link gives arbitrary-law
coefficient \(0.0250967461\).  With all four theorems inserted, the
diagnostic is \(0.3323683002\), leaving \(0.0009650332\), and sixteen of 888
balanced entries are controlled.  Decision: continue with
\((3,1,1,5):(0,0,1,4)\), whose reoptimized gate is
\(0.0542506298=1.73602015/q\), while retaining the hard-instance pivot if a
safe coefficient exceeds its gate.

Fifth chain-aware result: for \((3,1,1,5):(0,0,1,4)\), the generic
distinctness factor would give \(0.0934752746\) and fail the previous gate.
The exact quintic fixed-four row energy is \(1-4/N\), while the remaining
cubic--Hadamard matrix has rank at most \(N\) and entries at most \(1/N\).
Their arbitrary-law coefficient is \(0.0311889051<1/32\).  With all five
theorems inserted, the diagnostic is \(0.3323657941\), leaving
\(0.0009675392\), and twenty of 888 entries are controlled.  Decision:
continue with \((1,1,3,5):(0,1,1,3)\), whose reoptimized gate is
\(0.0570749885=1.82639963/q\).

Sixth chain-aware result: for \((1,1,3,5):(0,1,1,3)\), factoring the full
\(M_{13}M_{35}\) row leaves a repeated Hadamard matrix of weighted trace norm
at most one.  Exact L-shape incidence bounds the record-one row energy, and
record-three extensions split by the fixed triple's row pattern.  The
largest target-order tail is the exact horizontal value
\(3/[(q-1)(q-2)]\).  The resulting arbitrary-law coefficient is
\(0.0422410016<0.0570749885\).  With all six theorems inserted, the
diagnostic is \(0.3327757792\), leaving \(0.0005575541\), and twenty-four of
888 entries are controlled.  Decision: continue with
\((3,1,1,5):(0,1,1,3)\), whose reoptimized gate is \(0.0484819899\).

Seventh chain-aware result: for \((3,1,1,5):(0,1,1,3)\), the exact
fixed-three quintic row energy is \(22365/15872\).  Removing that Schur
feature leaves a repeated cubic--Hadamard matrix.  Duplicate compression and
the XOR-labelled cubic endpoint reduce it to a weighted Walsh matrix with
factor \(1/q\).  The resulting arbitrary-law coefficient is
\(0.0370952793<0.0484819899\).  With all seven theorems inserted, the
diagnostic is \(0.3329643636\), leaving \(0.0003689697\), and twenty-eight of
888 entries are controlled.  Decision: continue with
\((1,3,5,1):(0,2,2,1)\), whose reoptimized gate is \(0.0426269309\).

Eighth chain-aware result: for \((1,3,5,1):(0,2,2,1)\), retain the exact
cubic and quintic fixed-pair squared slices before taking the universal
\(M_{35}\) maximum.  The complete row energy is at most
\(N E_2F_2m_{35}^2\), giving arbitrary-law coefficient
\(0.0285281523<1/32\).  With all eight theorems inserted, the diagnostic is
\(0.3328775891\), leaving \(0.0004557442\), and thirty-two of 888 entries are
controlled.  Decision: continue with \((1,1,3,5):(0,0,3,2)\), whose
reoptimized gate is \(0.0454321892\).

Ninth chain-aware result: for \((1,1,3,5):(0,0,3,2)\), normalize the
complete fixed-pair \(M_{35}\) row but preserve the XOR-labelled
\(H_NM_{13}\) residual.  Exact pair-extension counts and the resulting
Walsh compression give arbitrary-law coefficient
\(0.0250919472<1/32\).  With all nine theorems inserted, the diagnostic is
\(0.3326862124\), leaving \(0.0006471209\), and thirty-six of 888 entries
are controlled.  Decision: continue with \((1,3,5,1):(0,2,3,0)\), whose
reoptimized gate is \(0.0529177167\).

Tenth chain-aware result: for \((1,3,5,1):(0,2,3,0)\), extract the scalar
cubic--quintic completion row while preserving both endpoint Walsh factors.
The residual repeated \(H_N\otimes H_N\) matrix has arbitrary-law
coefficient one.  Record-one incidence and record-three endpoint slices give
coefficient \(0.0462425962<0.0529177167\).  This exceeds the provisional
\(1/32\), so with all ten theorems inserted the diagnostic is
\(0.3331326055\), leaving only \(0.0002007278\), and forty of 888 entries
are controlled.  Decision: the next orbit
\((1,1,5,3):(0,1,3,1)\) is a bounded \(0.0379251204\) gate test, not the
whole program; complete one bounded Track B or C deliverable before opening
a second additional orbit.

### 3. Boundary map

Record the best passive lower bound, passive protocol evidence, active upper
bound, and active lower-bound evidence in one comparison.  Do not describe
the active upper bound six as an optimum.

Current entry: `BOUNDARY_MAP.md` records the accepted four-sided comparison.
The important blank is the passive algorithmic upper frontier.  Filling that
blank with a valid protocol or a certified bound for a precisely stated
protocol class is a program result, not merely a falsification aid.

### 4. Mechanism statement

Separate three levels:

- observed structure in exact or numerical examples;
- a proved instance-specific contraction; and
- a reusable theorem stable under adaptive posterior selection.

Only the third level supports a general explanation of the gap.

### 5. Route decision

For each serious hard instance or proof framework, choose one status:

- continue as lead;
- continue as a bounded secondary project;
- pivot after a named failed gate; or
- stop because a certified obstruction rules out the target.

Current entry: the signed-permutation plant remains the bounded Track A lead
through the complete one-batch gate.  The reverse-tree proof remains the
asymptotic baseline.  Its completed loss audit localizes the next Track B
target first to three exact terminal types at levels nine and ten.  Complete
enumeration proves \(N^{-1}\) for all level-eleven/twelve types.  A second
centered-weight expansion removes the unique level-ten forest and upgrades
the general floor to \(N^{1/18}\).  A second centered endpoint audit removes
both reflected level-nine trees and proves \(N^{1/16}\).  The high-level
terminal list is closed.  The complete low-level image then gives \(N^{-1}\)
at level eight and centered repairs at levels seven and six, proving
\(N^{1/12}\).  The next asymptotic decision is whether the limiting
level-twelve physical contraction or a replacement witness offers a
meaningful gain.
The active complete-flag reuse audit is resolved negatively; the
one-extra-query or interleaved-five family is the next bounded Track C
project.  See `PORTFOLIO.md`.

## Evidence-entry template

Every material result added here should include:

- **Claim:** the mathematical or protocol statement.
- **Label:** theorem, exact finite computation, certified numerical bound,
  exploratory evidence, or conjecture.
- **Scope:** hard instance, \(N\), dose, batch/adaptive model, and physical
  weight class.
- **Consequence:** which mission question changes.
- **Decision:** continue, strengthen, pivot, stop, or no change yet.
- **Artifact:** proof note, checker, certificate, or counterexample.

## Review trigger

Update this ledger after a theorem gate, a falsifying witness, a completed
hard-instance comparison, or a protocol result.  Routine algebra and
unproved numerical improvements stay in their working notes until they
change a program answer.
