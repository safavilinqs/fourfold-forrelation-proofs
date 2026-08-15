# Bessel-refined alternating path contraction

Date: 2026-07-14

Status: proved for arbitrary correlated diagonal weights on every
block-coherent alternating flattening.  It removes the growing unweighted
outer-link norms from the two remaining double-cubic path tensors.  The
fixed occurrence-split lifts were subsequently bounded in
double_endpoint_occurrence_contraction.md and
adjacent_double_cubic_occurrence_contraction.md; their joint weighted
packing remains open.

## 1. Refined coefficient

For

$$
K(a,b,c,d)=M_1(a,b)M_2(b,c)M_3(c,d),
$$

define the outer Bessel energies

$$
R_1=\max_a\sum_b|M_1(a,b)|^2,
\qquad
C_3=\max_d\sum_c|M_3(c,d)|^2,
$$

and the middle coherence \(\kappa_2=\max_{b,c}|M_2(b,c)|\).  The
alternating flattening obeys

$$
\boxed{
\|D_p^{1/2}K^{13\mid24}D_q^{1/2}\|_1
\le\sqrt{R_1C_3}\,\kappa_2\sqrt{PQ}.
}
\tag{1.1}
$$

Use the same latent-pair factorization as in the weighted path theorem:

$$
D_p^{1/2}L\,D_{M_2}\,R D_q^{1/2}.
$$

Its two Hilbert--Schmidt norms satisfy

$$
\|D_p^{1/2}L\|_2^2\le R_1P,
\qquad
\|D_{M_2}RD_q^{1/2}\|_2^2
\le\kappa_2^2C_3Q.
$$

Schatten Holder gives (1.1).  Unlike an operator-norm bound, this estimate
is insensitive to coherent multiplicity among different outer rows when
each individual row retains bounded Bessel energy.

## 2. Two decorated endpoints

For a signed-permutation cubic endpoint linked to a singleton, the
singleton character family is orthonormal.  Bessel therefore gives

$$
\sum_b|M_{3,1}(A,b)|^2\le1,
\qquad
\sum_c|M_{1,3}(c,D)|^2\le1.
$$

With \(M_2=H_N\), equation (1.1) proves coefficient \(1/q\) for the
block-coherent alternating cut of \((3,1,1,3)\).  The other whole-block
cuts have at least the same decay by the complete path table.

## 3. Endpoint adjacent to an L-shape

For the record-one sector of \((3,3,1,1)\), restrict the second cubic
support to the compatible L-shape family.  A direct conditional-permutation
count gives

$$
\boxed{
\max_A\sum_{B\in L}|M_{3,3}(A,B)|^2
=1+{2q\over(q-1)^2}.
}
\tag{3.1}
$$

For a type-I endpoint support the energy is one.  For a type-II support,
the exceptional conditional permutation terms occur exactly when two
independent Walsh signs are both equal.  There are
\(q(q/2-1)^2\) exceptional row/character choices; subtracting the baseline
\(1/(q-1)^2\) contribution gives the correction in (3.1).

The next link from the L-shape to a singleton has coherence
\(1/[q(q-1)]\), while the final singleton link has column energy one.
Thus the block-coherent alternating coefficient is

$$
{1\over q(q-1)}
\sqrt{1+{2q\over(q-1)^2}}.
\tag{3.2}
$$

At \(q=32\), (3.2) is about \(0.001041\), far below \(1/q\).

## 4. Remaining interface

The block-coherent result identifies the correct contraction mechanism,
but a cubic support split between ket and bra is an entry mask of the
outer Bessel frame, not merely a repeated row.  The next theorem must show
that distinct-label mask/Bessel packing preserves (1.1) across those
internal splits and their common base fibers.

Reproduction: tests/three_link_weighted_path_contraction.py and
searches/signed_permutation_double_cubic_entries.py.
