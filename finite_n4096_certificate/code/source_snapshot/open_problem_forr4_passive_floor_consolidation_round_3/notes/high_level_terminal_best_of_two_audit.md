# Exact high-level terminal audit and the \(N^{1/20}\) theorem

Date: 2026-07-15

Status: rigorous improvement of the general passive lower bound.  Exhaustive
enumeration of the exact four-layer Stein set dynamics, followed by the
better of the global assigned and safe all-projective contractions, proves

$$
D_{\mathsf P}^{\rm hard}=\Omega(N^{1/20}).
\tag{0.1}
$$

This improves the accepted \(\Omega(N^{1/24})\) theorem.  The improvement is
asymptotic; it does not imply passive dose \(>6\) at \(N=1024\).

The same audit identifies exactly three remaining high-level joint
saturators: two reflected connected trees at level nine and one
level-ten forest with component sizes \(6+4\).  These are now the precise
output of this coefficient-independent audit.  Follow-up
`level_ten_forest_mean_zero_repair.md` removes the level-ten type and proves
the intermediate \(\Omega(N^{1/18})\) theorem.  Follow-up
`level_nine_tree_centered_repair.md` removes the reflected level-nine pair
and proves the intermediate \(\Omega(N^{1/16})\) theorem.  The later
complete low-level image audit proves the current \(\Omega(N^{1/12})\)
theorem.

## 1. Why the high-level image is finite and exact

Let \(\beta_0\) be the number of internal blocks whose two initial covariance
endpoints coincide.  The initial branching potential is

$$
\mathcal E_0=4(3-\beta_0).
\tag{1.1}
$$

At a terminal leaf the potential equals the number \(v\) of marked
vertices.  Therefore \(v\ge9\) forces \(\beta_0=0\).  Up to independent
coordinate relabeling in the four layers, there is a unique high-level
initial state: one edge on each boundary and two distinct inherited marks
in each internal layer.

A state is completely determined by its three adjacent-layer edge sets.
For each active internal source, the next Stein transfer chooses:

1. any coordinate already marked in the adjacent layer; or
2. one fresh coordinate, unique up to relabeling.

A fresh transfer preserves \(\mathcal E\).  An existing-coordinate transfer
strictly decreases it.  Canonicalizing after every transfer under independent
within-layer permutations gives a finite breadth-first enumeration with no
heuristic truncation.

The exact counts are:

| quantity | count |
|---|---:|
| reachable canonical states | 222 |
| terminal states at all levels | 34 |
| terminal states at levels 9--12 | 17 |
| reflection-sensitive terminal states at levels 9--12 | 8 |

The levelwise counts are:

| level \(v\) | all terminal types | sensitive types |
|---:|---:|---:|
| 9 | 8 | 3 |
| 10 | 6 | 3 |
| 11 | 2 | 1 |
| 12 | 1 | 1 |

Every enumerated component meets all four layers.  The unique level-twelve
type is the three-path forest from
terminal_interpolation_sigma_one_witness.md.

The enumeration is coefficient-independent in the safe direction: every
actual terminal term belongs to one of these types.  A type that passes the
norm audit is controlled whether or not its signed coefficient later
cancels.

## 2. The two safe contractions

For a connected component \(C\), let

$$
\ell_C=e_C-v_C+1
\tag{2.1}
$$

be its cycle surplus.  For a physical placement, let \(k_C\) be the maximum
number of vertices of \(C\) in one physical amplitude entry.  The global
assigned route gives decay

$$
\delta_{\rm ass}
=\frac12\sum_C(\ell_C+k_C-1).
\tag{2.2}
$$

Only placements with assigned integer

$$
\sigma_{\rm ass}=2\delta_{\rm ass}=1
\tag{2.3}
$$

can miss the level-sensitive target: if \(\sigma_{\rm ass}\ge2\), the
assigned route already gives \(N^{-1}\).

The second route keeps the whole reverse frame skeleton projective.  For a
partition \(\pi_C\) of the vertices of \(C\) into physical-entry blocks,
let

$$
r_{\max}(C,\pi_C)
=\max_{\substack{S\text{ is a union}\\\text{of blocks of }\pi_C}}
\operatorname{rank}_{\mathbb F_2} A_C[S,S^c].
\tag{2.4}
$$

The exact graph flattening formula gives the safe component decay

$$
\delta_{\rm proj}(C,\pi_C)
=\frac{\ell_C+r_{\max}(C,\pi_C)-1}{2}.
\tag{2.5}
$$

Vertical injective-norm multiplicativity adds (2.5) over components.  The
same-layer distinctness mask is expanded into diagonal character unitaries,
so it changes only a diagram constant.  The physical contraction may choose
the better complete proof:

$$
\delta_{\rm best}
=\max\{\delta_{\rm ass},\delta_{\rm proj}\}.
\tag{2.6}
$$

Equation (2.6) does not mix norms inside one induction.  It takes the smaller
of two independently proved upper bounds after both complete contractions
have been established, so RT-003 is not involved.

For \(\sigma_{\rm ass}=1\), there are only two possibilities:

- all components are trees, exactly one has \(k_C=2\), and the rest have
  \(k_C=1\); or
- one component is unicyclic, every component has \(k_C=1\), and all other
  components are trees.

For the \(k_C=2\) component, the executable audit enumerates every partition
into singleton and pair blocks and minimizes (2.5).  For \(k_C=1\), the
partition is all-singleton and every vertex cut is available.

## 3. Complete best-of-two verdict

All eight sensitive high-level types admit assigned suppression one.  Five
pass the required decay and three fail:

| level | sensitive types | pass | joint saturators | worst safe decay |
|---:|---:|---:|---:|---:|
| 9 | 3 | 1 | 2 | \(1/2\) |
| 10 | 3 | 2 | 1 | \(1/2\) |
| 11 | 1 | 1 | 0 | \(1\) |
| 12 | 1 | 1 | 0 | \(1\) |

Thus every sensitive terminal graph at levels eleven and twelve has a safe
fixed-skeleton factor \(N^{-1}\).  Levels nine and ten retain the accepted
\(N^{-1/2}\) bound because of the three types below.

### Type A — upper-branching level-nine tree

Its three boundary edge sets are

$$
\begin{aligned}
E_{12}&=\{(0,0),(0,1)\},\\
E_{23}&=\{(0,0),(0,1),(1,2)\},\\
E_{34}&=\{(0,0),(1,1),(2,2)\}.
\end{aligned}
\tag{3.1}
$$

It is a connected nine-vertex tree with one first-layer vertex.  A
maximum-occupancy-two placement exists for which every entry-respecting cut
has binary rank at most two.  Hence

$$
\delta_{\rm ass}=\delta_{\rm proj}=\frac12
<\frac9{16}.
\tag{3.2}
$$

### Type B — reflected level-nine tree

The layer-reversed type has

$$
\begin{aligned}
E_{12}&=\{(0,0),(1,1),(2,2)\},\\
E_{23}&=\{(0,0),(1,0),(2,1)\},\\
E_{34}&=\{(0,0),(1,0)\}.
\end{aligned}
\tag{3.3}
$$

It has the same two decay exponents as Type A.  Types A and B form one
reflection orbit.

### Type C — level-ten \(6+4\) forest

Its boundaries are

$$
\begin{aligned}
E_{12}&=\{(0,0),(1,1),(2,2)\},\\
E_{23}&=\{(0,0),(1,0),(2,1)\},\\
E_{34}&=\{(0,0),(1,1)\}.
\end{aligned}
\tag{3.4}
$$

One component is a six-vertex tree and the other is a four-layer path.  Put
the paired entry in the six-vertex component.  Its worst projective decay is
zero; the singleton path contributes \(1/2\).  Therefore

$$
\delta_{\rm ass}=\delta_{\rm proj}=\frac12
<\frac{10}{16}=\frac58.
\tag{3.5}
$$

These types are exact reachable, reflection-sensitive terminal set states.
Their signed scalar coefficient grouping has not yet been audited.  They are
therefore precise proof-interface saturators, not yet lower-bound
counterexamples.

## 4. Improved asymptotic theorem

The accepted contraction supplies \(N^{-1/2}\) at levels four through ten.
The exhaustive audit supplies \(N^{-1}\) at levels eleven and twelve.  After
the accepted dose ledger and path-count constants, the task difference is
bounded by

$$
\sum_{v=4}^{10} C_v(1+D)^vN^{-1/2}
+
\sum_{v=11}^{12} C_v(1+D)^vN^{-1}.
\tag{4.1}
$$

Set \(D=cN^{1/20}\).  In the first sum, the worst dimension exponent occurs
at \(v=10\):

$$
N^{10/20}N^{-1/2}=1.
\tag{4.2}
$$

All smaller levels decay.  In the second sum,

$$
N^{v/20}N^{-1}\to0
\qquad(v=11,12).
\tag{4.3}
$$

Choosing the absolute constant \(c>0\) sufficiently small makes the total
distinguishability smaller than the task threshold.  The accepted promise
conditioning loss is \(O(N^{-1})\) and does not change the conclusion.
This proves (0.1).

At \(N=1024\), the bare asymptotic scale is only

$$
1024^{1/20}=\sqrt2.
\tag{4.4}
$$

The theorem therefore improves the genuine asymptotic result but does not
approach the finite-size goal of excluding passive dose six.

## 5. Next decisions

The next Track B work identified by this audit was finite and ranked.

1. Audit the signed terminal coefficient of Type C — completed.  Its forced
   centered local derivative repairs every dangerous partition and improves
   the general exponent to \(1/18\).
2. Audit one representative of the reflected Type A/B pair — completed.
   The forced centered endpoint derivative removes both orientations and the
   mechanism reaches \(1/16\).
3. The high-level list is closed.  Further improvement must address the
   level-eight-or-lower rows, theorem constants, or a different hard
   instance.  Do not reopen generic mixed Hilbert/projective induction.

Reproduction:

- searches/high_level_terminal_best_of_two_audit.py performs the complete
  canonical state enumeration, physical pair-partition search, and exponent
  calculation.
- tests/high_level_terminal_best_of_two_audit.py protects all state counts,
  sensitive types, the exact three-saturator list, the unique level-twelve
  type, and the \(1/20\) theorem exponent.
- terminal_three_path_projective_repair.md supplies the all-projective frame
  lemma used in the scoring step.
