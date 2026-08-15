# A true terminal level-twelve sigma-one witness

Date: 2026-07-15

Status: rigorous negative result for the first bounded OP3-7 route. The
specific star used in the earlier graph-interface audit is absent from the
true interpolation image, but a simpler reflection-sensitive level-twelve
forest is present with nonzero coefficient and exact assigned-fiber
suppression \(\sigma=1\). Therefore terminal-image exclusion alone cannot
improve the accepted \(\Omega(N^{1/24})\) exponent.

This does not falsify the desired posterior-stable contraction. It shows
that any improvement must group graph tensors or exploit physical
complete-frame packing; it cannot come from declaring every high-level
\(\sigma=1\) diagram nonexistent.

Follow-up: `terminal_three_path_projective_repair.md` now proves that this
exact forest contracts as \(N^{-1}\) when the entire frame skeleton remains
projective.  Thus the witness is a sharp obstruction to the assigned
parameter alone, but not to the level-twelve \(N^{-3/4}\) target.

## 1. The old relaxed star is not terminal

The prior interface witness attached all \(v-4\) extra layer-one vertices
to one layer-two path vertex. That graph obeys the coarse properties

$$
4\le v\le12,\qquad e\ge {3v\over4},
$$

but the exact Stein dynamics impose an additional outer-boundary degree
cap.

A layer-two vertex can receive one initial edge from layer one. If it is
unmatched toward layer one, it can be the source of exactly one leftward
transfer. That transfer resolves the mismatch permanently. No transfer is
sourced in layer one. Hence every layer-two vertex has at most two
layer-one neighbors. The layer-three/layer-four boundary is symmetric:

$$
\deg_{1\leftrightarrow2}(x)\le2,
\qquad
\deg_{3\leftrightarrow4}(y)\le2.
\tag{1.1}
$$

At \(v=12\), the relaxed star would require boundary degree nine. It is
therefore not in the true terminal image.

This exclusion does not yield extra dimension decay, because the true image
contains the following different witness.

## 2. Exact all-new transfer path

Use layer numbers \(1,2,3,4\). Start from the three covariance-derivative
edges

$$
\begin{aligned}
(1,0)&\mathbin{-}(2,0),\\
(2,3)&\mathbin{-}(3,2),\\
(3,1)&\mathbin{-}(4,1).
\end{aligned}
\tag{2.1}
$$

Both internal layers begin with distinct left- and right-inherited marks.
Apply these six legal Stein transfers:

| step | active source | direction | fresh neighbor |
|---:|---|---|---|
| 1 | \((2,3)\) | left | \((1,1)\) |
| 2 | \((2,0)\) | right | \((3,0)\) |
| 3 | \((3,1)\) | left | \((2,2)\) |
| 4 | \((3,0)\) | right | \((4,0)\) |
| 5 | \((3,2)\) | right | \((4,2)\) |
| 6 | \((2,2)\) | left | \((1,3)\) |

Every neighbor is new when its transfer occurs. The branching potential
therefore stays exactly twelve at every step. At the end, both internal
incidence sets are matched and the terminal graph is

$$
\begin{aligned}
P_0&:(1,0)-(2,0)-(3,0)-(4,0),\\
P_1&:(1,1)-(2,3)-(3,2)-(4,2),\\
P_2&:(1,3)-(2,2)-(3,1)-(4,1).
\end{aligned}
\tag{2.2}
$$

Thus

$$
v=12,\qquad e=9,\qquad c=3,
\tag{2.3}
$$

and every component is a four-layer tree. Labels are distinct within every
layer.

The executable replay checks all active-source conditions and finds 180
valid orders of the six transfers that give the same terminal edge set.

## 3. It survives the signed interpolation

The task difference is interpolated after antisymmetrizing in the first
block. A terminal Fourier support survives exactly when its first-layer
cardinality is odd. Diagram (2.2) has three first-layer vertices, so it is
reflection-sensitive.

Moreover, every transfer in Section 2 creates a new marked vertex. No
transfer differentiates an existing local weight. The terminal scalar
weight therefore contains only undifferentiated Stein-kernel factors and
\(\psi'\) factors. For the capped odd increasing function used by the
accepted plant,

$$
\psi'(x)>0,\qquad
\mathcal S_\psi(x)>0
\tag{3.1}
$$

for every finite \(x\). The first inequality is immediate from the Gaussian
derivative of the error function. For the second, oddness gives zero total
Gaussian integral, while the tail integral defining
\(\mathcal S_\psi\) is positive on either side of zero.

Hence this path has pointwise strictly positive local weight. Its 180 legal
orders add rather than cancel. The first-block antisymmetrization retains
it. The diagram is therefore genuinely present with nonzero coefficient,
not merely reachable in a coefficient-blind set recursion.

## 4. Exact suppression parameter

For every path component \(P_j\),

$$
\ell_j=e_j-v_j+1=0.
\tag{4.1}
$$

Its natural-cut binary rank is two, so the all-singleton projective branch
would retain

$$
\sigma_{\mathrm{proj}}
=\sum_{j=0}^2(\ell_j+r_j-1)
=3.
\tag{4.2}
$$

Now place vertices \((1,0)\) and \((2,0)\) in one physical amplitude entry
and put each of the other ten vertices in a distinct entry. This is a legal
assigned-fiber placement. The maximum component occupancies are

$$
(k_0,k_1,k_2)=(2,1,1).
\tag{4.3}
$$

The exact retained assigned parameter is therefore

$$
\begin{aligned}
\sigma_{\mathrm{ass}}
&=\sum_{j=0}^2(\ell_j+k_j-1)\\
&=(0+2-1)+(0+1-1)+(0+1-1)\\
&=\boxed{1}.
\end{aligned}
\tag{4.4}
$$

The global dichotomy in the repaired proof must use the assigned branch
once any component has a multi-vertex entry. It cannot simultaneously
recover the projective factors of the two weakly assigned paths without a
new mixed Hilbert/projective theorem; that is precisely the inference ruled
out by the existing RT-003 norm counterexample.

At level twelve, (4.4) leaves the accepted factor \(N^{-1/2}\), hence dose
exponent

$$
{1\over2\cdot12}={1\over24}.
\tag{4.5}
$$

The minimal \(N^{1/16}\) target would require \(N^{-3/4}\), equivalently an
additional \(N^{-1/4}\). At \(N=1024\), that is the factor

$$
1024^{-1/4}=0.1767766953.
\tag{4.6}
$$

## 5. Decision for Track B

The first OP3-7 fork is closed:

- the particular relaxed star is absent, but that fact is immaterial;
- a true sensitive level-twelve \(\sigma=1\) forest is present;
- its scalar branching coefficient is positive, so scalar interpolation
  signs do not remove it; and
- the existing graph invariants still yield only \(N^{-1/2}\).

The subsequent bounded target is now complete.  Rather than combine an
assigned Hilbert component with separately projective components, keep the
entire frame skeleton projective.  The paired path has grouped injective
norm at most one and the two singleton paths each give \(N^{-1/2}\), for a
safe total \(N^{-1}\).  This avoids the invalid RT-003 conversion and beats
the required \(N^{-3/4}\).

The next asymptotic target is therefore the complete true terminal image at
levels nine through twelve.  Every diagram/placement must be scored by the
better of its global assigned and all-projective bounds; merely finding
assigned suppression one no longer identifies a genuine obstruction.

Reproduction:

- searches/terminal_interpolation_sigma_one_witness.py replays the exact
  branching path and computes both suppression parameters.
- tests/terminal_interpolation_sigma_one_witness.py checks the graph,
  sensitivity, forest structure, physical placement, 180 valid transfer
  orders, and exponent consequence.
- terminal_three_path_projective_repair.md and its executable certificate
  prove that the physical all-projective bound removes this witness as an
  \(N^{1/16}\) obstruction.
