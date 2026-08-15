# Joint occurrence-profile packing

Date: 2026-07-14

Status: proved the exact unordered-support packing inequality for every
odd four-block profile of total degree at most twelve.  This removes the
separate sum over time, ket/bra, and internal split assignments at the
square-mass level.  It does not by itself control the terminal transcript
norm or mix different Fourier profiles.

## 1. One block

Fix the parity supports in all query-node/frame-copy containers.  For one
physical block, let

$$
w(x)=\#\{\text{containers whose parity support contains }x\},
\qquad M=\sum_xw(x).
$$

Take an unordered Fourier support \(A\) of size \(a\).  The number of
ways to assign each \(x\in A\) to a container containing it, summed over
all such \(A\), is

$$
[z^a]\prod_x(1+w(x)z)=e_a(w).
\tag{1.1}
$$

Replace coordinate \(x\) by \(w(x)\) distinguishable slots.  The left
side of (1.1) counts only slot subsets that never choose two copies of the
same coordinate, while \(\binom Ma\) counts every \(a\)-slot subset.
Therefore

$$
\boxed{e_a(w)\le\binom Ma.}
\tag{1.2}
$$

This is the key normalization point: the Fourier support is an unordered
set.  The right side is a binomial coefficient, not a falling factorial,
so no labeled-mark factorial is introduced.

## 2. Four blocks and two frame copies

For a profile \(a=(a_1,a_2,a_3,a_4)\), apply (1.2) independently in the
four disjoint physical blocks.  Conditional on all realized supports, the
joint assignment count is at most

$$
\prod_{b=1}^4\binom{M_b}{a_b}.
\tag{2.1}
$$

The two complete-frame square sums use two independent support copies.
If the total hard query dose is \(D\), then pointwise

$$
M_1+M_2+M_3+M_4\le2D.
\tag{2.2}
$$

Taking expectations preserves (2.1).  Hence the complete joint
occurrence square mass obeys

$$
\boxed{
\mathcal S_a
\le C_a(D):=
\max_{\substack{M_b\in\mathbb Z_{\ge0}\\\sum_bM_b\le2D}}
\prod_{b=1}^4\binom{M_b}{a_b}.
}
\tag{2.3}
$$

For \(a=(1,1,1,1)\), (2.3) gives \(C_a(6)=3^4=81\), exactly recovering
the earlier minimal-chain square-mass bound.

## 3. Exact dose-six constants

There are seventy odd profiles with \(\sum_ba_b\le12\).  Exhaustive
integer enumeration of (2.3) gives

| total degree | profiles | largest \(C_a(6)\) | largest \(\sqrt{C_a(6)}\) |
|---:|---:|---:|---:|
| 4 | 1 | 81 | 9 |
| 6 | 4 | 160 | \(4\sqrt{10}\) |
| 8 | 10 | 126 | \(3\sqrt{14}\) |
| 10 | 20 | 36 | 6 |
| 12 | 35 | 1 | 1 |

The degree-eight maximum occurs at a single quintic decoration, for
example profile \((5,1,1,1)\) with capacity \((9,1,1,1)\).  A
double-cubic profile has the smaller value \(100\), attained at capacities
\((5,5,1,1)\) up to permutation.  At total degree twelve, every nonzero
term exhausts the capacity and its square mass is one.

Attenuation by \(\beta=5/6\) makes the largest square-root masses by
degree approximately

$$
4.3403,\quad4.2362,\quad2.6106,\quad0.9690,\quad0.1122.
\tag{3.1}
$$

These are per-profile maxima, not a valid sum over representation
sectors.

## 4. What remains

Equation (2.3) solves the combinatorial part of joint occurrence packing:
there is no exponential split count and no hidden factorial.  It also
shows that degrees ten and twelve are automatically occupation-sparse at
dose six.

The missing analytic step is to combine (2.3) with the cut-dependent
signed-permutation contractions before taking a terminal absolute norm.
A universal square-function promotion is impossible by the known
two-node counterexample.  The next contraction must retain the exact
chain/link sector weights; using only the worst fixed-split coefficient
is quantitatively too coarse, especially for \((3,1,1,3)\).

This last point has an explicit dose-six witness.  Put all probe mass on
occupation state \((2,1,1,2)\) and sum the current row/column slice bound
over the sixty-four splits of \((3,1,1,3)\).  The safe coefficient is

$$
2.4563273559,
$$

or \(0.5712632372\) after multiplication by \((5/6)^8\).  By comparison,
the attenuated minimal sector plus promise conditioning already uses
\(0.1729752014\), leaving only \(0.1603581320\) below \(1/3\).  Even if
all other higher sectors vanished, the double-endpoint estimate would
need a factor below \(0.280708\) of its current splitwise sum.

For the adjacent profile, the analogous deterministic state
\((2,2,1,1)\) gives only \(0.0251196\) after attenuation.  The finite-size
priority is therefore a joint two-endpoint Schur contraction across the
central Hadamard, not further polishing of the endpoint-to-L incidence
constant.

Reproduction: `tests/joint_occurrence_profile_packing.py` and
`searches/fixed_split_occupation_barrier.py`.
