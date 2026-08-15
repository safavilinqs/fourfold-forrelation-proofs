# Complete terminal image, low-level centered repair, and the \(N^{1/12}\) theorem

Date: 2026-07-15

Status: rigorous extension of the terminal audit to every initial covariance
collision pattern, followed by exact centered repairs of the only remaining
safe-norm saturators at levels seven and six.  The four sensitive level-eight
types already have decay \(N^{-1}\); the old generic \(N^{-1/2}\) level-eight
row was not sharp on the true terminal image.  Combining the low-level audit
with the accepted level-nine and level-ten repairs improves the passive floor
from

$$
\Omega(N^{1/16})
\quad\text{to}\quad
\boxed{\Omega(N^{1/12})}.
\tag{0.1}
$$

This remains asymptotic.  At \(N=1024\), its bare scale is

$$
1024^{1/12}=1.7817974363,
\tag{0.2}
$$

so it still does not prove passive hard dose greater than six.

## 1. Why the earlier image was incomplete below level nine

Let \(\beta_0\) be the number of the two internal initial covariance blocks
whose endpoints coincide.  The initial branching potential is

$$
\mathcal E_0=4(3-\beta_0).
\tag{1.1}
$$

For levels nine through twelve, only \(\beta_0=0\) can contribute, which is
why the earlier high-level enumeration had one initial state.  The complete
low-level image must also include both orientations with \(\beta_0=1\) and
the path with \(\beta_0=2\).  Up to independent relabeling inside each
layer, the four initial potentials are

$$
(12,8,8,4).
\tag{1.2}
$$

Fresh transfers preserve the potential and transfers to existing
coordinates strictly lower it.  Canonical breadth-first enumeration from all
four initial states gives 236 distinct reachable states and 39 terminal
types.  Twenty-two terminal types are reflection-sensitive.

The exact level counts are:

| level | all terminal types | sensitive types | pre-repair minimum safe decay |
|---:|---:|---:|---:|
| 4 | 1 | 1 | \(1/2\) |
| 5 | 2 | 2 | \(1\) |
| 6 | 4 | 3 | \(1/2\) |
| 7 | 6 | 4 | \(1/2\) |
| 8 | 9 | 4 | \(1\) |
| 9 | 8 | 3 | \(1/2\) |
| 10 | 6 | 3 | \(1/2\) |
| 11 | 2 | 1 | \(1\) |
| 12 | 1 | 1 | \(1\) |

For each type the safe score is the better of the complete assigned and
grouped-entry projective contractions, capped at one power because larger
values are not needed for this theorem.  The level-nine and level-ten
minimum-\(1/2\) rows are exactly the types already removed by the two earlier
centered repairs.  At low levels there is exactly one minimum-\(1/2\) type at
level seven and exactly one at level six.

The main immediate correction is that level eight was never a true
obstruction.  All four sensitive types have safe decay at least \(N^{-1}\)
without a coefficient argument.

## 2. The unique level-seven saturator

The graph consists of two four-layer paths sharing their first vertex:

$$
a_0-b_0-c_0-d_0,
\qquad
a_0-b_1-c_1-d_1.
\tag{2.1}
$$

It is a seven-vertex tree.  Exactly ten maximum-occupancy-two physical-entry
partitions have entry-respecting cut rank two; both previously accepted safe
routes give only \(N^{-1/2}\) on those partitions.

### 2.1 Exact coefficient histories

Every initial edge triple contained in (2.1) that can reach level seven is
included.  The exact counts are:

| initial potential | initial triples | contributing triples | histories |
|---:|---:|---:|---:|
| 8 | 4 | 4 | 12 |
| 12 | 2 | 2 | 24 |
| total | 6 | 6 | 36 |

The 36 histories give eight terminal local-weight profiles.  Potential-eight
histories have two fresh and one existing-coordinate transfer; potential-
twelve histories have one fresh and two existing-coordinate transfers.
Every history differentiates an existing local factor exactly once, at the
same outer site:

$$
h(a_0)=\gamma'(a_0)=\psi''(a_0).
\tag{2.2}
$$

Thus \(h\) is odd and Gaussian-centered.  The interpolation-time exponents
are \((1,1,1)\) in every profile.  The repair is performed separately on each
history, so no cancellation among unrelated scalar coefficients is assumed.

### 2.2 Centered expansion

Apply the bounded Stein kernel of \(h\) across the first boundary.  The
derivative either hits \(b_0\) or \(b_1\), or it creates a new \(b_2\) and
the resulting raw odd factor transfers right until it hits an occupied
coordinate or the final layer.  Extending every dangerous physical
partition by every new mark gives:

| branch | exact cases | minimum decay exponent |
|---|---:|---:|
| marked \(b_0,b_1\) | 20 | \(1\) |
| new \(b\), marked \(c_0,c_1\) | 106 | \(1\) |
| new \(b,c\), marked \(d_0,d_1\) | 586 | \(3/2\) |
| new \(b,c,d\) | 1,692 | \(1\) |

Every branch remains reflection-sensitive because the number of distinct
first-layer vertices stays one.  Consequently the original level-seven
family satisfies

$$
|\mathfrak C_7|
\le C\sum_{v=7}^{10}(1+D)^vN^{-1}.
\tag{2.3}
$$

## 3. The unique level-six saturator

The graph is

$$
a_0-b_0-c_0-d_0,
\qquad
b_0-c_1-d_1.
\tag{3.1}
$$

It is a six-vertex tree.  Exactly 31 occupancy-two, cut-rank-two partitions
give the old \(N^{-1/2}\) score.

There are two potential-eight initial triples, both contributing.  Their four
legal histories give two profiles.  Every history has one fresh transfer,
one existing-coordinate transfer, and the same centered derivative on the
incoming side of \(b_0\):

$$
h(b_0)=\gamma'(b_0)=\psi''(b_0).
\tag{3.2}
$$

The time exponents are \((0,1,1)\).  Expanding (3.2) to the left has only two
cases:

- the derivative hits the marked \(a_0\); the 31 resulting parallel-edge
  cases all have minimum decay one;
- the derivative creates a fresh \(a_1\); its 144 physical placements have
  two distinct first-layer vertices and cancel from the task difference
  before absolute values.

Therefore

$$
|\mathfrak C_6|
\le C(1+D)^6N^{-1}.
\tag{3.3}
$$

This repair is not strictly required to reach exponent \(1/12\), because the
old \(N^{-1/2}\) level-six row already meets that scaling.  It is included to
close the complete low-level obstruction list and to isolate the true next
frontier.

## 4. The improved global theorem

After combining the complete image audit with all centered repairs, the
sensitive transcript has the schematic bound

$$
\Delta
\le
C_4(1+D)^4N^{-1/2}
+\sum_{v=5}^{12}C_v(1+D)^vN^{-1}.
\tag{4.1}
$$

The constants absorb the finite level shifts produced by each centered
expansion.  Equation (4.1) uses:

- the direct best-of-two \(N^{-1}\) score at levels five and eight;
- the level-six and level-seven repairs above;
- the accepted reflected level-nine and level-ten repairs; and
- the direct \(N^{-1}\) score at levels eleven and twelve.

Set \(D=cN^{1/12}\).  The level-twelve term is made small by choosing the
absolute constant \(c\) sufficiently small.  Every other term decays; for
example,

$$
N^{4/12}N^{-1/2}=N^{-1/6},
\qquad
N^{11/12}N^{-1}=N^{-1/12}.
\tag{4.2}
$$

The accepted promise-conditioning loss is \(O(N^{-1})\).  Hence

$$
D_{\mathsf P}^{\rm hard}=\Omega(N^{1/12}).
\tag{4.3}
$$

## 5. Program decision

The complete terminal obstruction list is now closed through all levels.
The proof is no longer limited by the generic level-eight row.  Its next
asymptotic bottleneck is the \(N^{-1}\) level-twelve contraction, which gives
the ratio \(1/12\).  The follow-up
`level_twelve_contraction_sharpness.md` gives a matching
\(\Theta(N^{-1})\) lower witness for the arbitrary grouped graph norm.  A
gain beyond \(N^{1/12}\) must therefore use additional physical frame
restrictions, remove the positive coefficient family before terminal norms,
or change the hard instance; after that, level eleven becomes the next
likely rung.

The practical conclusion remains more sobering.  Improving the bare
\(N=1024\) scale from \(1.542\) to \(1.782\) still does not approach six, and
the theorem constants are not explicit enough to help.  Round 3 should keep
finite-size certification, passive protocol searches, active-boundary work,
and alternative hard instances coequal with any next asymptotic audit.

Reproduction:

- `searches/low_level_terminal_centered_repair.py` enumerates the 236-state
  union, all 39 terminal types, both exact coefficient-history families, and
  every retained/cancelled repair branch.
- `tests/low_level_terminal_centered_repair.py` protects all image counts,
  levelwise safe minima, coefficient histories, branch counts, decay minima,
  cancellations, and the \(1/12\) exponent.
