# Round 3 charter

Date: 2026-07-15

## Goal in plain language

Work out how real, how large, and how useful the active--passive advantage is
for four-forrelation.

Earlier rounds proved that an advantage exists: passive hard dose must grow
at least as \(\Omega(N^{1/24})\), while an active protocol uses at most six.
Round 3 should now determine what that statement means in practice and how
close it is to the truth.

Current progress: Round 3 has improved the rigorous passive floor to
\(\Omega(N^{1/12})\).  It classified the complete terminal image from every
initial collision pattern and used forced centered local weights to remove
the level-ten, level-nine, level-seven, and level-six obstructions.  The
finite terminal obstruction list is closed.  This remains far from a useful
finite-size separation at \(N=1024\): the bare scale is only about \(1.782\),
before constants.

This is broader than finishing the current \(N=1024\) calculation.  The
current calculation is one important experiment inside the program.

## The five questions

Round 3 should leave the clearest rigorous answer it can to five questions:

1. **How much passive dose is necessary?**  Improve the asymptotic lower
   bound, improve its constants, or identify the precise obstruction to doing
   so with the present proof framework.
2. **Does the advantage appear at useful sizes?**  Decide whether passive
   dose greater than six can be proved around \(N=1024\), and record how the
   answer changes over a small range of realistic powers of two rather than
   relying on one isolated size.
3. **How good can passive sensing actually be?**  Search for valid passive
   protocols and counterexamples to proposed lower-bound mechanisms.  A
   lower bound is informative only when compared with the best upper
   evidence we know.
4. **Is active dose six the true benchmark?**  Seek a protocol below six or
   a genuine lower bound against such protocols.  Six is currently an upper
   bound, not an established optimum.
5. **What causes the gap?**  Separate facts peculiar to one hard distribution
   from a reusable principle that explains why adaptively chosen passive
   measurements lose information that coherent active traversal retains.

The first and fourth questions map the two complexity scales.  The second
makes the result operational.  The third prevents the lower-bound program
from becoming self-confirming.  The fifth determines whether the work teaches
us something beyond one certificate.

## Three coordinated tracks

### Track A — finite-size truth

Determine what can be proved about passive dose six at realistic \(N\).
Continue the signed-permutation plant while its explicit gates pass, but
compare or replace it when a certified obstruction appears.  A one-batch
certificate is an intermediate result; the adaptive lift is required for a
full passive lower bound.

### Track B — asymptotic scale and mechanism

Continue beyond the new \(N^{1/12}\) floor or prove a sharp limitation of
the repaired reverse-tree framework.  Since the complete terminal list is
closed and the arbitrary grouped level-twelve norm now has a matching
\(\Theta(N^{-1})\) lower witness, the next work should test physical frame
restrictions, usable theorem constants, or a replacement witness.  The
preferred output is a reusable, posterior-stable principle or a quantitative
boundary change, not merely a larger catalog or another generic graph norm.

### Track C — algorithms, counterexamples, and alternative witnesses

Challenge both sides of the claimed gap.  Look for passive counterprotocols,
active protocols below six, active lower bounds, and hard instances with
better finite-size or adaptive behavior.  Use the same physical model and
scorecard for every comparison.

These tracks are complementary.  Track A asks whether the gap is visible,
Track B asks how it scales and why, and Track C asks whether the current
algorithms and witnesses are the right ones.

## What counts as progress

A result materially advances Round 3 when it does at least one of the
following:

- proves passive hard dose greater than six at a realistic size;
- improves the asymptotic passive floor or its usable constants;
- proves a reusable passive contraction that survives adaptive posterior
  selection;
- finds a better passive or active protocol;
- gives a meaningful lower bound on active dose;
- identifies a quantitatively better hard instance; or
- rigorously rules out a serious route and supplies a justified pivot.

Local coefficient improvements count only through the program decision they
enable.  Numerical patterns are evidence and falsification tools, not
theorems.

## Success levels

- **Minimum successful round:** produce a decisive verdict on the current
  finite-size route and one independent result that changes the asymptotic,
  active-boundary, passive-protocol, or hard-instance picture.
- **Strong round:** prove either a realistic-size separation or a stronger
  asymptotic theorem, while seriously testing the opposite side of the
  comparison.
- **Best outcome:** give a useful finite-size separation, a stronger general
  theorem, and a credible account of the mechanism and near-optimal boundary.

A negative verdict can satisfy a success condition when it is rigorous,
rules out a route under the actual model, and identifies what should replace
that route.

## End-of-round answer

Round 3 should end with a short synthesis that says:

1. the strongest active--passive separation now proved;
2. what is and is not established near \(N=1024\);
3. the best known lower and upper evidence for passive dose;
4. the best known lower and upper evidence for active dose;
5. which mechanism is proved, which is conjectural, and how adaptivity enters;
6. whether the current hard instance should be continued or replaced; and
7. the single most valuable next theorem or protocol experiment.

`BOUNDARY_MAP.md` maintains the comparison.  `MISSION_LEDGER.md` maintains
the accepted answers.  `PORTFOLIO.md` limits each track to a bounded live
project, and `PLAN.md` defines the gates and pivots.

## Scope boundary

The mission is broadened within constant-margin four-forrelation and its
active/passive hard-dose model.  Round 3 is not an open-ended survey of
unrelated sensing problems, and it is not an exposition-polishing round.
