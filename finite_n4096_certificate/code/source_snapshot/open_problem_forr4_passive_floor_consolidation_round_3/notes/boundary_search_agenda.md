# Boundary-search agenda

Date: 2026-07-14

Status: Track C initialization scaffold.  These are bounded falsification
projects, not an authorization for open-ended parameter sweeps.

## Passive side

### C-P1 — Small valid passive protocols

At the smallest tractable \(N\) and doses, search structured passive protocols
that exploit repeated modes, joint measurements, or adaptive posterior
selection.  Record the actual protocol value separately from PPT, separable,
or other relaxation values.

Decision value: identify a resource missing from the lower-bound norm, or
produce evidence that the norm tracks valid protocols sharply.

### C-P2 — Saturators for the reverse-tree losses

For each proposed Track B improvement, construct a small physical family that
tries to saturate it.  Start with the known two-copy square-function witness,
which falsifies a proof shortcut but is not itself a full counterprotocol to
the theorem.

Decision value: rule out overstrong general inequalities before investing in
their proofs.

## Active side

### C-A1 — Can six doses be reduced?

Write the six-dose active protocol as an exact resource ledger.  Test whether
coherent flags, path reuse, or a different collective measurement can remove
one charged traversal without weakening the constant-margin decision.

Decision value: an explicit protocol below six, or a precise obstruction to
the most credible reuse.

Result: removing the third complete flag and replacing the remaining
measurements by an arbitrary collective POVM is rigorously obstructed at
\(N=1024\).  The exact endpoint ensemble has two-copy trace distance
\(0.2776778892\), hence Helstrom error \(0.3611610554>1/3\).  The result
includes both path and mode registers, but not an extra charged traversal.

### C-A2 — Lower bounds below six

Identify which active resources defeat the known passive proof and formulate
the smallest active lower-bound question that can still be attacked exactly,
beginning with exclusion of very low hard dose.

Decision value: narrow the true active benchmark rather than merely preserve
the current upper bound.

## Selection rule

Run at most one passive and one active boundary search at a time.  Each search
must name its finite search space, certification method, and the program
decision that would change under a positive or negative result.

## Initialized first projects

### Passive project: structured small-\(N\) protocol screen

- **Question:** can a valid passive protocol exploit a feature discarded by
  the current one-batch or reverse-tree ledger?
- **Initial scope:** the smallest exact Sylvester sizes, hard dose at most six,
  and a symmetry-reduced family with repeated modes and at most two adaptive
  batches.
- **Required output:** an actual protocol value and its full hard-dose
  ledger, clearly separated from any relaxation value.
- **Stop rule:** stop after the finite family is exhausted or a protocol
  exposes a named missing resource; do not expand the family without a new
  hypothesis.
- **Program consequence:** challenge the finite-size lower-bound route and
  inform the passive-boundary row of `MISSION_LEDGER.md`.

### Active project: one-extra-query audit

- **Question:** can one extra charged traversal, coherently integrated with
  two folds, overcome the exact two-complete-flag obstruction?
- **Initial scope:** two completed folded-chain flags plus one coherent
  one-query side experiment, or one explicitly specified interleaved
  five-traversal word, with no postselection discount.
- **Required output:** a verified protocol, or an obstruction for the stated
  reuse family—not a general active lower bound unless separately proved.
- **Stop rule:** stop when the one-extra-query/interleaved-five family is
  resolved; broaden only if it produces a quantitatively better margin.
- **Program consequence:** update the active benchmark rather than assuming
  that six is optimal.
