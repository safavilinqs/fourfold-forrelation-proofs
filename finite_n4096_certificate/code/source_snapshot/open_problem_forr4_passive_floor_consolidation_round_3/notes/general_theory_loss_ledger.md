# General-theory loss ledger

Date: 2026-07-15

Status: Track B exponent, terminal-image, physical-norm, complete high-level,
and all centered-weight audits completed.  The current theorem is
\(\Omega(N^{1/12})\).  The complete terminal obstruction list is closed;
the accepted level-twelve \(N^{-1}\) row is now limiting and is sharp for
the arbitrary grouped-entry graph norm.

## Corrected theorem ladder

| Level | Transcript scale | Implied passive floor | Status |
|---|---:|---:|---|
| Accepted theorem | \((1+D)^{12}N^{-1/2}\) | \(\Omega(N^{1/24})\) | proved |
| Exact fixed-diagram refinement | \(D^vN^{-\sigma(G,\phi)/2}\), \(\sigma\ge1\) | still \(N^{1/24}\) in the worst interface case | proved |
| High-level best-of-two theorem | \(D^vN^{-1/2}\) for \(v\le10\), \(D^vN^{-1}\) for \(v=11,12\) | \(\Omega(N^{1/20})\) | proved |
| Centered level-ten repair | old Type C becomes \(((1+D)^{10}+(1+D)^{11})N^{-1}\) | \(\Omega(N^{1/18})\) | proved |
| Centered reflected-tree repair | repaired Type A/B families become finite sums of \((1+D)^vN^{-1}\) | \(\Omega(N^{1/16})\) | proved |
| Minimal level-sensitive improvement | \(D^vN^{-\max\{1/2,v/16\}}\) | \(\Omega(N^{1/16})\) | proved by classification plus centered repairs |
| Natural second-suppression sufficient condition | \(D^vN^{-1/2}\) for \(v\le8\), \(D^vN^{-1}\) for the repaired high-level families | \(\Omega(N^{1/16})\) | proved in the needed terminal image |
| Complete low-level image and repair | \((1+D)^4N^{-1/2}+\sum_{v=5}^{12}(1+D)^vN^{-1}\) | \(\Omega(N^{1/12})\) | proved |
| Level-twelve graph-norm sharpness | upper \(O(N^{-1})\), masked lower \(N^{-1}(1-N^{-1})^2(1-2N^{-1})^3\) | the present grouped-norm route cannot beat \(N^{1/12}\) | proved limitation |
| Uniform component-scale improvement | \(D^vN^{-v/8}\) | \(\Omega(N^{1/8})\) | open |
| Generic passive Fourier growth | posterior-stable bound to be stated | target near \(N^{1/4}/\log N\) | speculative |

The former row \(D^{2v}N^{-v/8}\Rightarrow N^{1/16}\) is obsolete.  The
accepted joint insertion/Bessel audit already removed the duplicate
\(D^v\) charge.  With the current dose ledger, \(N^{-v/8}\) would imply
\(N^{1/8}\) directly.

## Exact loss in the accepted contraction

For a graph component \(C\), write

$$
\ell_C=e_C-v_C+1
$$

for edge surplus and

$$
k_C=\max_E|C\cap E|.
$$

The repaired proof retains the following square-suppression integer before
relaxing it to one:

$$
\sigma(G,\phi)=
\begin{cases}
\sum_C(\ell_C+r_C-1),
  &k_C=1\text{ for every }C,\\[1mm]
\sum_C(\ell_C+k_C-1),
  &\max_C k_C\ge2,
\end{cases}
\tag{2.1}
$$

where \(r_C\) is the natural-cut binary rank in the all-singleton branch.
Every component spans four layers, so \(\sigma\ge1\).  The fixed-diagram
bound is

$$
|\mathfrak C_{G,\phi}|
\le C_G(1+D)^vN^{-\sigma(G,\phi)/2}.
\tag{2.2}
$$

The accepted transcript proof uses only \(\sigma\ge1\).  At level \(v\), the
corresponding dose exponent is \(\sigma/(2v)\).

## Completed interface stress test

For every \(4\le v\le12\), a connected four-layer tree with \(e=v-1\) and
natural-cut rank two has two legal placements—one all-singleton and one with
maximum physical-entry occupancy two—that satisfy

$$
\sigma_{\rm proj}=\sigma_{\rm ass}=1.
$$

It obeys every layered-graph and placement property used by the current
reverse-tree theorem.  Therefore edge surplus, component count, exact
flattening rank, and maximum entry occupancy do not by themselves imply
either \(N^{-v/16}\) or \(N^{-v/8}\).

The particular star-shaped family above is a limitation of the theorem
interface and is not terminal: its outer-boundary degree exceeds the exact
Stein cap.  A subsequent exact replay nevertheless finds a true terminal
replacement at \(v=12\): three disjoint four-layer paths produced by six
all-new transfers.  It is reflection-sensitive, has positive local weight,
and a placement with component occupancies \((2,1,1)\) gives
\(\sigma_{\rm ass}=1\).  Thus a stronger theorem must use graph-tensor
grouping or physical posterior frame structure, not support exclusion or
scalar branching signs.

## Where the \(N^{1/16}\) target began

Levels four through eight already satisfy the \(N^{1/16}\) scaling target
under the accepted \(N^{-1/2}\) contraction.  Only levels nine through
twelve need more dimension decay.  The extra \(N\)-powers required are

| level | extra exponent beyond \(1/2\) |
|---:|---:|
| 9 | \(1/16\) |
| 10 | \(1/8\) |
| 11 | \(3/16\) |
| 12 | \(1/4\) |

At \(N=1024\), these are factors \(0.648420\), \(0.420448\),
\(0.272627\), and \(0.176777\).  Thus the corrected general theorem target
also points directly at the degree-ten/twelve finite-size frontier.

## Precise OP3-7 target

Prove or sharply falsify:

> For every adaptive passive tree, after signed terminal interpolation terms
> are grouped at each level and before terminal absolute values, the
> level-\(v\) contraction is at most
> \(C^v(1+D)^vN^{-\max\{1/2,v/16\}}\), uniformly in outcome width,
> tree depth, posterior-selected descendants, and physical insertion fibers.

The target is now proved.  The first bounded subproblem and its full
high-level follow-up are resolved.
The three-path witness passes with \(N^{-1}\).  Complete enumeration then
proves \(N^{-1}\) for all level-eleven/twelve types and returns three joint
saturators.  The unique level-ten \(6+4\) forest is now repaired by its
forced odd centered local derivative.  Both reflected level-nine trees are
then repaired by their forced outer centered derivative.  Every retained
branch has decay at least \(N^{-1}\), and the sole all-fresh reflected family
cancels.  The subsequent complete low-level audit adds the potential-eight
and potential-four initial states.  All sensitive level-eight types already
have decay one, and centered derivatives repair the unique level-seven/six
saturators.  This proves \(N^{1/12}\).  The next theorem target must improve
the level-twelve one-power row using physical frame restrictions or replace
the witness; more terminal shape or arbitrary grouped-norm work is not
indicated.  The explicit masked sharpness witness matches \(N^{-1}\) up to
the factor \((1-N^{-1})^2(1-2N^{-1})^3\).  The displayed initial triple has
180 legal orders, and complete enumeration gives 1,080 positive histories
over all 12 initial triples.

## Updated loss inventory

| Source | Accepted treatment | Audit verdict | Next evidence |
|---|---|---|---|
| Interpolation degree | levels 4--12 | all 22 sensitive terminal types pass after centered repairs | use a physical-frame restriction or change the witness |
| Marked-node placements | one \(D^v\) charge | already improved; no second charge remains | preserve, do not reopen |
| Graph suppression | \(\sigma\mapsto1\) | best-of-two plus centered expansion closes every terminal type; the limiting graph norm is exponent-sharp | stop generic graph-norm sharpening |
| Complete-frame packing | one global dichotomy | the safe all-projective fallback plus centered repairs prove \(N^{1/12}\) | explicit constants or an optimizer-excluding physical frame identity |
| Adaptive outcomes | operator-valued frontier | no width/depth loss, but signed level grouping is unproved | theorem plus witness tests |
| Promise conditioning | additive \(O(N^{-1})\) | not exponent-limiting | retain accepted bound |

## Stop and promotion rules

- Do not report \(N^{1/16}\) from merely retaining existing component
  factors; the checked interface family forbids that inference.
- Do not pursue terminal-support exclusion further; the true signed image
  already contains a positive level-twelve \(\sigma=1\) diagram.
- Do not call an assigned \(\sigma=1\) diagram a saturator until its whole
  grouped-entry projective norm has also been checked.
- Do not seek another power from the arbitrary grouped-entry norm of the
  level-twelve forest; the same-layer-distinct lower witness is
  \(\Theta(N^{-1})\).
- Stop the level-sensitive target if a physical posterior witness violates
  it.
- Keep a gain instance-specific if it does not survive outcome-selected
  descendants.
- Promote a lemma to the mechanism statement only when its level, dimension,
  dose, and adaptive scope are explicit.

The exponent derivation is in `general_theory_exponent_audit.md`; the exact
terminal construction and its signed survival audit are in
`terminal_interpolation_sigma_one_witness.md`; and its safe projective repair
is in `terminal_three_path_projective_repair.md`.  The complete image audit
and improved intermediate theorem are in
`high_level_terminal_best_of_two_audit.md`.  The current theorem is in
`level_nine_tree_centered_repair.md` and
`low_level_terminal_centered_repair.md`.  The matching lower witness for the
limiting graph norm is in `level_twelve_contraction_sharpness.md`.
