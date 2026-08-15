# Round 3 retrospective and Round 4 handoff

Date: 2026-07-15

Status: end-of-round synthesis.  This file is the shortest authoritative
account of what Round 3 proved, what it only diagnosed, how the work answers
the original program, and why a fourth round is warranted.

## Executive verdict

Round 3 achieved a strong asymptotic result and a sharp framework
limitation:

$$
D_{\mathsf P}^{\rm hard}=\Omega(N^{1/12}),
\qquad
D_{\mathsf A}^{\rm hard}\le6.
\tag{0.1}
$$

It did not prove passive hard dose greater than six at \(N=1024\), did not
find a passive protocol matching six, and did not determine whether active
dose five is possible.

The practical answer is therefore still unknown.  The new asymptotic
theorem has bare scale

$$
1024^{1/12}=1.7817974363,
\tag{0.2}
$$

far below six before theorem constants.  Moreover, the immediate attempt to
improve the exponent again by sharpening the limiting level-twelve grouped
graph norm is closed: a full-distinctness lower witness matches its
\(N^{-1}\) upper power and retains \(0.9922113065\) of that value at
\(N=1024\).

Round 4 should therefore lead with the realistic-size boundary, passive
protocol evidence, and hard-instance comparison.  It should reopen the
asymptotic contraction only when a concrete physical-frame restriction or a
different witness escapes the sharp optimizer.

## 1. What the original program asked

The first folder posed an exact constant-margin problem:

- prove that unrestricted classically adaptive passive hard dose grows as
  \(N^{\Omega(1)}\) for positive-versus-negative Hadamard
  four-forrelation;
- retain the active hard-dose-six protocol; and
- if polynomial passive growth is false, construct and verify an explicit
  \(N^{o(1)}\)-dose passive counterprotocol.

The model was deliberately broad: arbitrary fresh quantum batches,
entanglement, idlers, repeated modes, collective POVMs, vacuum coherence,
classical feed-forward, and branchwise hard photon-pass accounting.
Changing the promise, access class, or dose meter did not count.

That entry goal was achieved before Round 3.  The initial theorem gave
\(\Omega(N^{1/48})\); consolidation removed a duplicate dose charge and
gave \(\Omega(N^{1/24})\).  Round 2 then repaired four genuine contraction
defects and made the \(1/24\) theorem the accepted baseline.

Once a polynomial separation was secure, the program's real questions
became broader:

1. How large is the passive complexity, asymptotically?
2. Is the separation visible at realistic sizes such as \(N=1024\)?
3. Is six the correct active benchmark?
4. Is the current hard instance the right one?
5. Can a passive protocol or counterexample show that the lower-bound
   program is aiming too high?

Round 3 was designed to address all five rather than turn one proof ledger
into the entire research program.

## 2. Round 3 goals versus outcomes

| Original or inherited goal | Round 3 outcome | Verdict |
|---|---|---|
| Preserve the repaired unrestricted adaptive theorem | Complete inherited-plus-Round-3 regression suite passes | achieved |
| Improve the general passive exponent | Improved \(1/24\to1/20\to1/18\to1/16\to1/12\) | achieved |
| Prove passive hard dose \(>6\) at \(N=1024\) | No complete one-batch or adaptive certificate | open |
| Determine whether the current contraction can improve further | Limiting level-twelve grouped graph norm has matching \(\Theta(N^{-1})\) lower witness | achieved as a negative result |
| Test the signed-permutation finite-size witness | Ten chain-aware theorems and a near-threshold completion diagnostic, but 848 entries remain provisional | substantial progress, not a theorem |
| Establish a passive protocol upper frontier | No competitive valid passive protocol consolidated | open |
| Test whether active dose below six is possible | Two complete flags fail even under optimal collective decoding; genuinely interleaved five-traversal protocols remain open | partial |
| Compare alternative hard instances under common gates | Scorecard and pivot rules created; no decisive quantitative comparison completed | open |

## 3. Main theorem achieved in Round 3

### 3.1 Complete terminal image

The reverse-tree image was extended from the single potential-twelve start
to all four initial covariance collision patterns with potentials

$$
(12,8,8,4).
\tag{3.1}
$$

The union contains:

- 236 reachable canonical states;
- 39 terminal types;
- 22 reflection-sensitive terminal types; and
- exact terminal levels four through twelve.

This corrected an important relaxation: every sensitive level-eight type
already has safe decay \(N^{-1}\).  The generic \(N^{-1/2}\) level-eight row
was not attained on the true image.

### 3.2 Exact centered repairs

Four remaining terminal mechanisms were repaired:

- the level-ten \(6+4\) forest;
- both reflected level-nine trees;
- the unique level-seven tree; and
- the unique level-six tree.

Each retained family has a forced odd centered derivative.  Exact
re-expansion and physical-partition enumeration gives \(N^{-1}\) on every
retained branch; the specified all-fresh branches cancel by task
antisymmetry.

After these repairs, the sensitive transcript obeys the schematic theorem

$$
\Delta
\le
C_4(1+D)^4N^{-1/2}
+\sum_{v=5}^{12}C_v(1+D)^vN^{-1}.
\tag{3.2}
$$

Setting \(D=cN^{1/12}\) proves (0.1).

### 3.3 Why the next obvious exponent route stops

The unique sensitive level-twelve type is a forest of three four-layer
Hadamard paths.  A legal physical placement groups two adjacent vertices of
one path and leaves the other ten singleton.  The accepted all-projective
upper bound is \(O(N^{-1})\).

Round 3 constructs eleven grouped unit vectors for which the full
same-layer-distinct tensor has value

$$
N^{-1}(1-N^{-1})^2(1-2N^{-1})^3.
\tag{3.3}
$$

The upper and lower dimension powers therefore match.  Complete coefficient
enumeration gives 12 initial triples and 1,080 all-fresh histories, all with
positive local weight.  The family does not vanish through internal history
cancellation.

This is a limitation of the current graph-norm architecture, not a passive
algorithm and not proof that \(N^{1/12}\) is optimal.  A stronger theorem
must use restrictions on physically realizable posterior frames, cancel
terms before the grouped norm, or change the hard instance.

## 4. Realistic-size work

Round 3 resumed the attenuated signed-permutation route at
\(N=q^2=1024\), \(q=32\).  Its accepted finite-size assets now include:

- the exact attenuated plant variance and promise conditioning;
- the exact minimal one-batch dose-six certificate;
- exact or chain-aware bounds for ten leading balanced orbit families;
- an occupation-compatible degree-eight total
  \(0.2482247496\);
- an accepted partial total through the proved high sectors
  \(0.2797585469\); and
- a coarse completion diagnostic
  \(0.3331326055<1/3\).

That last number is not a theorem.  It has margin only

$$
{1\over3}-0.3331326055
=0.0002007278,
\tag{4.1}
$$

and 848 balanced entries remain provisional.  The diagnostic shows the
route is quantitatively serious, but also fragile.  Continuing one orbit at
a time is justified only when a credible contraction sits below the
reoptimized gate; otherwise the correct response is a global theorem,
hard-instance pivot, or protocol comparison.

No Round 3 result proves:

- one-batch passive dose six impossible at \(N=1024\);
- adaptive passive dose six impossible;
- a finite-size lower bound after promise conditioning; or
- that the signed-permutation plant is better than alternatives.

## 5. Active and passive algorithm boundaries

The active upper bound remains six.  Round 3 proved that simply reducing the
known three complete flags to two cannot work at \(N=1024\), even with an
arbitrary collective measurement.  The exact two-copy endpoint ensemble has

$$
\text{Helstrom error}=0.3611610554>{1\over3}.
\tag{5.1}
$$

This rules out complete-flag reuse as the reason active dose might fall
below six.  It does not rule out a genuinely interleaved five-traversal
circuit or two flags plus one coherent extra query.

The passive upper side remains the least developed part of the boundary.
Round 3 initialized a structured search agenda but did not consolidate a
competitive valid passive protocol.  Absence of such a protocol is not
evidence that the passive lower bound is close to optimal.

## 6. What Round 3 changed conceptually

Round 3 established three research decisions.

1. **The asymptotic separation is stronger and structurally better
   understood.**  The complete terminal list is finite, exact, and repaired.
2. **The practical problem is not solved by reading the exponent literally.**
   The bare \(N=1024\) scale remains far below six, and theorem constants are
   unusable.
3. **Generic grouped graph norms have reached a real barrier.**  The
   level-twelve optimizer survives every distinctness mask nearly unchanged
   at \(N=1024\).  Further norm work needs a physical constraint, not another
   combinatorial catalog.

This justifies a new round.  The next work should be judged by whether it
changes the actual \(N=1024,D=6\) boundary, supplies counterprotocol
evidence, or identifies a hard instance/contraction that escapes the sharp
obstruction.

## 7. Authoritative Round 3 reading order

For the end-of-round state, read:

1. this retrospective;
2. STATUS.md for the chronological result ledger;
3. BOUNDARY_MAP.md for the four-sided active/passive map;
4. notes/low_level_terminal_centered_repair.md for the \(N^{1/12}\) theorem;
5. notes/level_twelve_contraction_sharpness.md for the framework limitation;
6. MISSION_LEDGER.md for proved versus missing outputs; and
7. BASELINE_RUN.md for executable verification.

PLAN.md, PORTFOLIO.md, and OPEN_PROBLEMS.md record the Round 3 research
process and stop rules.  They are retained as archive evidence rather than
the live plan after Round 4 begins.

## 8. Handoff requirements for Round 4

Round 4 should begin from the following frozen facts:

- passive hard dose is rigorously \(\Omega(N^{1/12})\);
- active hard dose is at most six;
- passive \(>6\) at \(N=1024\) is not proved;
- the signed-permutation completion number is diagnostic, not certified;
- the generic level-twelve grouped graph norm is exponent-sharp;
- no competitive passive upper frontier is known; and
- active dose five remains open beyond the failed two-complete-flag family.

It should not spend time polishing the Round 3 exposition or re-enumerating
terminal shapes.  Its lead question is whether the separation is real and
visible at \(N=1024\), with explicit permission to change the hard instance
or find a passive counterprotocol when the present route fails.

The live successor is
../open_problem_forr4_passive_floor_consolidation_round_4/README.md.
