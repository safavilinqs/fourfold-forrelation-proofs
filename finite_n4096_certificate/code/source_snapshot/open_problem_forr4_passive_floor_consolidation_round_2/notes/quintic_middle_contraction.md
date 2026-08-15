# Quintic middle-decoration contraction

Date: 2026-07-14

Status: proved for every fixed occurrence split and arbitrary diagonal
weights, with coefficient \(q/(q-1)^2\).  This closes the two
single-decoration degree-eight profiles locally.  Joint occurrence/Bessel
packing remains open.

## 1. Conditional permutation coefficient

Consider the profile \((1,5,1,1)\).  The singleton neighbors force the
five-edge second-block support \(B\) to have exactly one odd row and one
odd column.

For the \((1,5)\) link, condition the hidden permutation on the required
odd-label match.  Any remaining active row label has even degree.  With
five edges there are at most two such active even rows.

For one active even row, the conditional Walsh average over the \(q-1\)
unused preimages has magnitude \(1/(q-1)\).  For two active even rows with
nonzero characters \(\alpha,\gamma\), the normalized injective average is

$$
{1\over(q-1)(q-2)}
\sum_{s\ne t;\ s,t\ne y}
(-1)^{s\cdot\alpha+t\cdot\gamma}.
$$

If \(\alpha=\gamma\), its magnitude is \(1/(q-1)\); otherwise it is
\(2/[(q-1)(q-2)]\), no larger than \(1/(q-1)\) for \(q\ge4\).

There cannot be zero active even rows: five edges in the one odd row would
give five odd columns, contradicting the singleton record on the next
link.  Therefore

$$
\max|M_{1,5}(a,B)|\le {1\over q(q-1)}.
\tag{1.1}
$$

The transposed argument gives the same bound for \(M_{5,1}(B,c)\).

## 2. Four-block entry and rank

Multiplying the two adjacent quintic-link estimates and the last
singleton Hadamard entry gives

$$
|K(a,B,c,d)|
\le {1\over q^3(q-1)^2}.
\tag{2.1}
$$

The profile has eight marked coordinates.  In a fixed ket/bra occurrence
split, the smaller marked side contains at most four coordinates, so its
basis dimension and the matrix rank are at most

$$
N^4=q^8.
$$

For arbitrary diagonal masses \(P,Q\), weighted rank--Frobenius yields

$$
\boxed{
\|D_p^{1/2}KD_q^{1/2}\|_1
\le {q\over(q-1)^2}\sqrt{PQ}.
}
\tag{2.2}
$$

At \(q=32\), the coefficient is

$$
{32\over961}=0.033298647\ldots
=1.065557\ldots\,{1\over q}.
\tag{2.3}
$$

This is slightly above the minimal \(1/q\) scale but far below a constant.
Path reversal closes \((1,1,5,1)\) with the same estimate.

## 3. Scope

The proof is local to one occurrence split.  Injective ancillary
repetitions can be aggregated without increasing rank, as in the cubic
case.  The result does not permit independent unit mass for every split;
the global finite-size sum must preserve their shared passive occupation
law.

At this stage the remaining open degree-eight profiles all had cubic
decorations in two blocks.  Later round-two notes bound all six such
profiles at fixed split.  Joint cut-dependent packing, not another local
single-middle estimate, is the remaining interface.

Reproduction: searches/signed_permutation_quintic_middle_bound.py.
