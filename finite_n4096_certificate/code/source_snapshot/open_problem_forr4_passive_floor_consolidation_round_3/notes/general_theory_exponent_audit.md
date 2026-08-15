# Corrected general-theory exponent audit

Date: 2026-07-15

Status: exact audit of the accepted proof and a proved limitation of its
current graph interface.  This note corrects the research ladder and states a
posterior-stable replacement target.  It does not prove a stronger passive
floor.

## 1. The dose ledger was already improved

For a terminal interpolation diagram with \(v\) marked vertices, the accepted
joint insertion/Bessel audit gives

$$
|\mathfrak C_G(\mathcal T)|
\le C_G(1+D)^v N^{-1/2}.
\tag{1.1}
$$

The older \(D^{2v}\) charge was removed before the
\(\Omega(N^{1/24})\) theorem was accepted.  Therefore the inherited research
ladder

$$
D^{2v}N^{-v/8}\Longrightarrow N^{1/16}
$$

is historical.  Under the accepted \(D^v\) ledger, the same dimension
estimate would give

$$
\boxed{D^vN^{-v/8}\Longrightarrow D=\Omega(N^{1/8}).}
\tag{1.2}
$$

There is no remaining dose-bookkeeping step between \(N^{1/16}\) and
\(N^{1/8}\).  An \(N^{1/16}\) theorem now corresponds to a weaker
level-sensitive dimension improvement.

## 2. Exact suppression retained by the repaired proof

The repaired reverse-tree proof contains more information than (1.1).
Write

$$
\ell_C=e_C-v_C+1
$$

for the cyclomatic edge surplus of a connected graph component.

### All-singleton branch

If every component has at most one vertex in every physical amplitude entry,
let \(r_C\) be the binary rank of its natural layer-\((1,3)\) versus
layer-\((2,4)\) cut.  The exact component singular value gives the square
suppression

$$
\sigma_{\rm proj}
=\sum_C(e_C+r_C-v_C)
=\sum_C(\ell_C+r_C-1).
\tag{2.1}
$$

The interpolation audit proves that every component spans all four layers,
so \(r_C\ge2\) and \(\sigma_{\rm proj}\ge1\).

### Assigned-fiber branch

Otherwise, let

$$
k_C=\max_E |C\cap E|
$$

and assign each component to a maximizing physical entry.  The exact
fiber-summed range diagonal contributes

$$
\sigma_{\rm ass}
=\sum_C(e_C-v_C+k_C)
=\sum_C(\ell_C+k_C-1).
\tag{2.2}
$$

Every \(k_C\ge1\), and at least one is at least two, so
\(\sigma_{\rm ass}\ge1\).

Before the final relaxation, the proof therefore establishes

$$
\boxed{
|\mathfrak C_{G,\phi}(\mathcal T)|
\le C_G(1+D)^vN^{-\sigma(G,\phi)/2},
\qquad \sigma(G,\phi)\ge1,
}
\tag{2.3}
$$

with \(\sigma\) given by (2.1) or (2.2).  The accepted theorem replaces every
\(\sigma\) by one.

For one fixed level, (2.3) would force dose exponent

$$
\alpha(G,\phi)={\sigma(G,\phi)\over2v}.
\tag{2.4}
$$

The worst accepted pair \((v,\sigma)=(12,1)\) yields \(1/24\).

## 3. What is actually needed for \(N^{1/16}\)

With only the accepted \(N^{-1/2}\), the per-level dose exponents are

| level \(v\) | accepted exponent |
|---:|---:|
| 4 | \(1/8\) |
| 5 | \(1/10\) |
| 6 | \(1/12\) |
| 7 | \(1/14\) |
| 8 | \(1/16\) |
| 9 | \(1/18\) |
| 10 | \(1/20\) |
| 11 | \(1/22\) |
| 12 | \(1/24\) |

Thus levels four through eight already meet an \(N^{1/16}\) theorem.  Only
levels nine through twelve need improvement.  The mathematically minimal
levelwise target is

$$
|\mathfrak C_v(\mathcal T)|
\le C^v(1+D)^v
N^{-\max\{1/2,v/16\}}.
\tag{3.1}
$$

The required extra \(N\)-powers at levels \(9,10,11,12\) are respectively

$$
{1\over16},\quad {1\over8},\quad {3\over16},\quad {1\over4}.
\tag{3.2}
$$

A more discrete but stronger sufficient target is one additional
\(N^{-1/2}\) for every level at least nine:

$$
|\mathfrak C_v(\mathcal T)|
\le
\begin{cases}
C^v(1+D)^vN^{-1/2},&4\le v\le8,\\
C^v(1+D)^vN^{-1},&9\le v\le12.
\end{cases}
\tag{3.3}
$$

Either (3.1) or (3.3) implies

$$
D=\Omega(N^{1/16}).
$$

The stronger \(N^{1/8}\) floor requires \(N^{-v/8}\) at every level.

## 4. Exact limitation of the current graph interface

For every \(4\le v\le12\), construct one connected layered tree as follows:

1. vertices \(0,1,2,3\) form a path through layers \(1,2,3,4\);
2. the other \(v-4\) vertices are layer-one leaves adjacent to vertex one;
3. for the projective branch, put every vertex in a distinct physical
   amplitude entry; and
4. for the assigned branch, put vertices zero and one in one entry and every
   remaining vertex in a distinct entry.

This family satisfies every graph property used by the accepted contraction:

- it is loopless and uses only adjacent-layer edges;
- it spans all four layers;
- \(e=v-1\ge3v/4\); and
- its natural-cut binary rank is exactly two.

For the all-singleton placement,

$$
e+r-v=1
$$

in the projective branch.  For the paired placement, the maximum entry
occupancy is \(k=2\), so

$$
e-v+k=1.
$$

Both exact retained parameters are therefore

$$
\boxed{\sigma_{\rm proj}=\sigma_{\rm ass}=1}
\tag{4.1}
$$

for every level through twelve.

This is a sharp limitation of any argument using only
\((v,e,r,k,\ell)\) and the present global dichotomy.  The family is a legal
layered graph and physical-entry placement at the theorem interface.  The
particular high-degree star above has since been excluded from the terminal
interpolation image.  However, the follow-up audit constructs a different
true terminal three-path forest with the same assigned parameter one.  That
forest is not a counterexample to a stronger physical or posterior-stable
theorem: the subsequent all-projective audit contracts it as \(N^{-1}\).

## 5. The next replacement lemma

The smallest useful Track B theorem is now:

> **Level-sensitive posterior contraction.** For every adaptive passive tree,
> group the signed terminal interpolation terms at each level before terminal
> absolute values.  Uniformly over posterior-selected descendants and
> physical insertion fibers, the resulting level-\(v\) contraction obeys
> (3.1).

The theorem cannot follow from the current graph invariants alone because of
(4.1).  Terminal support and scalar branching-sign cancellation fail on the
explicit positive three-path witness, but one safe all-projective norm
regime recovers both weak path factors and gives \(N^{-1}\).  A full proof
must show that every remaining high-level terminal type similarly passes
the better of the assigned and grouped-entry projective contractions, or
obtain fractional dimension decay from another physical norm inequality.

This is narrower than asking immediately for generic \(N^{-v/8}\): it targets
only levels nine through twelve and only the extra powers in (3.2).

## 6. Connection to realistic size

At \(N=1024\), the minimal gains in (3.2) are approximately

| level | extra factor |
|---:|---:|
| 9 | \(0.648420\) |
| 10 | \(0.420448\) |
| 11 | \(0.272627\) |
| 12 | \(0.176777\) |

The discrete second-suppression target (3.3) would instead give \(1/32\) at
all four levels.  Degree ten and twelve are exactly the live high-degree
frontier in the signed-permutation \(N=1024\) ledger, so the mechanism is
potentially transferable.  The hard distributions and coefficient
normalizations differ, however; no finite-size bound follows automatically.

## 7. Decision

OP3-6 is complete.  The old \(N^{1/16}\) rung has been corrected, the exact
loss has been localized to high interpolation levels with \(\sigma=1\), and
the current graph interface has a checked saturating family.

The subsequent OP3-7 image audit is now decisive.  The star used in Section
4 is not terminal because it violates the exact outer-boundary degree cap.
However, a true reflection-sensitive level-twelve forest of three
four-layer paths is generated by six all-new Stein transfers.  Its scalar
weight is positive, and assigning two vertices of one path to one physical
entry while leaving all other vertices singleton gives
\(\sigma_{\rm ass}=1\).  Thus terminal support and scalar interpolation
signs do not supply the missing decay.

The subsequent physical gate also passes.  The whole frame skeleton may be
kept projective; the paired path has grouped-entry injective norm at most
one, while the other two paths each retain \(N^{-1/2}\).  Thus the exact
forest contracts as \(N^{-1}\), stronger than the \(N^{-3/4}\) target,
without any RT-003 Hilbert-to-projective conversion.

That enumeration is now complete.  It finds 222 reachable canonical states,
34 terminal types, and eight sensitive types at levels nine through twelve.
Five pass the \(N^{1/16}\) row.  Every level-eleven/twelve type gains
\(N^{-1}\), while two reflected level-nine trees and one level-ten \(6+4\)
forest saturate both safe routes at \(N^{-1/2}\).  Consequently the current
best-of-two theorem improves to \(\Omega(N^{1/20})\).

The next OP3-7 gate is now also complete.  Every one of the 200 legal
histories producing the level-ten forest has a forced odd centered local
derivative.  Its exact re-expansion gives decay \(N^{-1}\) on all 282
dangerous original partitions and 5,295 retained extensions.  This removes
the level-ten obstruction and improves the current full theorem to
\(\Omega(N^{1/18})\).

The final high-level gate is now complete too.  All 200 histories of the
upper-branching level-nine representative force an odd centered outer
derivative.  Exact expansion and direct reflection give decay at least
\(N^{-1}\) on every retained branch; the all-fresh reflected branch cancels.
This removes both level-nine trees and proves \(\Omega(N^{1/16})\).  The
generic level-eight row appeared limiting at this stage.

The complete low-level audit has now corrected that last relaxation.  It
includes both potential-eight initial collision patterns and the
potential-four path, giving 236 reachable states and 39 terminal types.  All
four sensitive level-eight types already have safe decay \(N^{-1}\).  Forced
centered derivatives repair the unique level-seven and level-six saturators.
Together with the high-level repairs this proves \(\Omega(N^{1/12})\).  The
level-twelve one-power row is now limiting.  The follow-up sharpness audit
constructs a full-distinctness lower witness
\(N^{-1}(1-N^{-1})^2(1-2N^{-1})^3\), so arbitrary grouped graph norms cannot
strengthen that row.  The next asymptotic work must use physical frame
restrictions or change the hard instance rather than enumerate more terminal
shapes.

Reproduction:

- `../open_problem_forr4_passive_floor_consolidation_round_2/REPAIRED_REVERSE_TREE_CONTRACTION.md`
  is the authoritative repaired graph/fiber contraction.
- `../open_problem_forr4_passive_floor_consolidation/audits/06_DOSE_LEDGER_STRENGTHENING.md`
  is the accepted single-charge \(D^v\) audit.
- `../open_problem_forr4_passive_floor_consolidation/audits/02_INTERPOLATION_INDEPENDENT_AUDIT.md`
  proves the level and four-layer-component properties used here.
- `searches/general_theory_exponent_audit.py` checks the exact exponent
  algebra and constructs the saturating family at every level.
- `tests/general_theory_exponent_audit.py` protects the corrected ladder,
  graph ranks, suppression parameters, and \(N^{1/16}\) target.
- `terminal_interpolation_sigma_one_witness.md`,
  `searches/terminal_interpolation_sigma_one_witness.py`, and
  `tests/terminal_interpolation_sigma_one_witness.py` resolve the first
  OP3-7 fork by constructing and exhaustively checking the true terminal
  level-twelve witness.
- `terminal_three_path_projective_repair.md`,
  `searches/terminal_three_path_projective_repair.py`, and
  `tests/terminal_three_path_projective_repair.py` prove and check the
  witness-specific \(N^{-1}\) all-projective repair.
- `high_level_terminal_best_of_two_audit.md`,
  `searches/high_level_terminal_best_of_two_audit.py`, and
  `tests/high_level_terminal_best_of_two_audit.py` prove the complete
  high-level classification and the \(N^{1/20}\) theorem.
- `level_ten_forest_mean_zero_repair.md`,
  `searches/level_ten_forest_mean_zero_repair.py`, and
  `tests/level_ten_forest_mean_zero_repair.py` prove the centered repair and
  the intermediate \(N^{1/18}\) theorem.
- `level_nine_tree_centered_repair.md`,
  `searches/level_nine_tree_centered_repair.py`, and
  `tests/level_nine_tree_centered_repair.py` prove both reflected repairs
  and the intermediate \(N^{1/16}\) theorem.
- `low_level_terminal_centered_repair.md`,
  `searches/low_level_terminal_centered_repair.py`, and
  `tests/low_level_terminal_centered_repair.py` prove the complete image,
  low-level repairs, and current \(N^{1/12}\) theorem.
