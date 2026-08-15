# Passive four-forrelation consolidation: round 2

Status: frozen consolidation record.  The repaired asymptotic
\(\Omega(N^{1/24})\) lower bound passed the round-two audit.  A passive
hard-dose-six lower bound near \(N=1024\) was not completed.

Start with [ROUND_2_SUMMARY.md](ROUND_2_SUMMARY.md).  It records the final
verdict, the repairs that are required for the asymptotic theorem, the
realistic-size progress, the rejected approaches, and the round-three
handoff.

## Goal

Stress-test the reverse-tree contraction underlying the $\Omega(N^{1/24})$ passive lower bound until we are confident it is genuinely correct. Prioritize falsification, explicit checks, and proof correctness; do not spend time polishing exposition.

After that claim survives the audit, seek either a stronger contraction or a different hard instance that proves passive dose greater than 6 at realistic sizes, especially $N\sim 10^3$.

## Starting point

- Frozen input: ../open_problem_forr4_passive_floor_consolidation
- Current repaired claim: passive dose lower bound $\Omega(N^{1/24})$ for four-forrelation.
- The reverse-tree contraction is accepted only in the repaired form in
  REPAIRED_REVERSE_TREE_CONTRACTION.md and with the qualifications in
  CONFIDENCE_REPORT.md.
- Numerical evidence may falsify a claim or guide a proof, but it does not replace a rigorous uniform bound.

## Working rule

Every task in this directory should do one of three things: try to break the contraction, certify a delicate step, or improve the quantitative lower bound. Expository cleanup is out of scope until the mathematics is settled.

Round two is now closed.  New mathematical work belongs in the round-three
folder; edits here should be limited to correcting the frozen record or
maintaining its regression suite.

## Success criteria

1. The reverse-tree invariant is stated in a finite, mechanically checkable form with every space, index, norm, and multiplicity explicit.
2. Independent symbolic, exhaustive-small-instance, adversarial numerical, and optimization-based checks find no counterexample—or expose a precise gap.
3. The audit explicitly covers adaptive histories, masks/base fibers, vector-valued bilateral contractions, cross-entry fibers, and all factors depending on $N$, dose, outcomes, and depth.
4. A confidence report records what was proved, what was merely tested, and the residual assumptions.
5. Only after items 1–4 pass, pursue a certified bound that excludes dose 6 near $N=1024$, or replace the hard instance with one that does.
