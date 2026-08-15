# Two-sided cubic middle-block Schur lift

Date: 2026-07-14

Status: exact uniform-orbit spectrum proved for the balanced cubic split;
all 32 occurrence placements diagonalized directly at \(q=4\).  A coarse
general-diagonal inequality is proved for every fixed occurrence split and
every injective ancillary repetition.  Global packing of overlapping
occurrence assignments remains open.

## 1. Compatible middle supports

Consider the profile \((1,3,1,1)\).  The singleton first and third blocks
force record size one on both sides of the cubic second block.  Its support
must therefore be a three-edge \(2\times2\) L-shape.  Let \(z\) be the
missing corner, equivalently the pair consisting of its odd row and odd
column.  The two adjacent link moments multiply to

$$
M_{1,3}(a,B)M_{3,1}(B,c)
={1\over(q-1)^2}H_N(a,z)H_N(z,c).
\tag{1.1}
$$

Thus two-sided path compatibility removes the isolated cubic norm
\(\sqrt{(q^2+2)/6}\) completely.

## 2. Balanced occurrence split

Put one of the three middle-block coordinates on the ket and two on the
bra.  Put singleton blocks one and three on the ket and block four on the
bra.  The physical kernel has rows \((e,a,c)\), columns
\((\{f,g\},d)\), and is nonzero only when \(\{e,f,g\}\) is an L-shape.

After the two outer Walsh orthogonality sums, the singular spectrum is

$$
{\sqrt2\over q(q-1)^2}
\quad\text{with base multiplicity}\quad
{q(q-1)^2(q+2)\over2},
\tag{2.1}
$$

and

$$
{\sqrt2\over q(q-1)}
\quad\text{with base multiplicity}\quad q(q-1).
\tag{2.2}
$$

Both multiplicities are repeated \(N=q^2\) times by the last singleton
coordinate.  Dividing the nuclear norm by the square root of the two
orbit dimensions gives

$$
\boxed{
{q+4\over q^3\sqrt{q^2-1}}.
}
\tag{2.3}
$$

At \(q=4\), this equals

$$
0.0322748612184,
$$

the largest of all 32 directly diagonalized uniform occurrence
placements.  At \(q=32\), the same balanced orbit has value about

$$
3.434\times10^{-5}.
$$

The asymptotic scale is \(q^{-3}=N^{-3/2}\), far below both the minimal
chain and the cubic endpoint orbit.

## 3. General diagonal weights for one occurrence split

Every nonzero entry of the full four-block cubic-middle tensor has
magnitude

$$
{1\over q^3(q-1)^2}.
\tag{3.1}
$$

There are six distinct marked coordinates in total.  If \(s\) are on the
ket and \(6-s\) on the bra, the two marked basis dimensions are at most
\(N^s\) and \(N^{6-s}\).  Hence every occurrence-split flattening has
rank at most

$$
N^{\min(s,6-s)}\le N^3=q^6.
$$

For arbitrary diagonal masses \(P,Q\), weighted rank--Frobenius now gives

$$
\boxed{
\|D_p^{1/2}KD_q^{1/2}\|_1
\le {1\over(q-1)^2}\sqrt{PQ}
\le {1\over q}\sqrt{PQ}
\quad(q\ge3).
}
\tag{3.2}
$$

The same proof survives arbitrary injective ancillary repetitions of the
marked rows and columns: repetitions do not increase rank after their
diagonal masses are aggregated.  At \(q=32\), the coefficient in (3.2)
is \(1/961\).

## 4. What remains

Equation (3.2) closes arbitrary diagonal weights inside one fixed marked
occurrence assignment.  It does not authorize summing all ket/bra mark
placements with separate unit masses.  A physical probe law is shared by
those placements, and common ket/bra base supports can appear in several
overlapping decompositions.  They must be combined by the distinct-label
mask/Bessel packing used in the repaired reverse contraction, not counted
as independent rank.

The next theorem target is therefore a joint occurrence-packing lemma for
the minimal and cubic sectors.  The local cubic coefficient already has
ample slack: it is \(1/961\) before attenuation at \(q=32\), versus the
allowed \(1/32\) scale.

Reproduction: searches/signed_permutation_middle_cubic_schur_lift.py.
