# Weighted three-link path contraction

Date: 2026-07-14

Status: proved for every block-coherent flattening of a fixed support
sector and every correlated diagonal weight on that flattening.  The
theorem recovers the exact powers of \(N\) in the minimal-chain certificate
and does not duplicate either middle-block mass.  It is not yet the full
higher-sector one-batch Schur contraction, because a multi-coordinate
support inside one physical block can be split between ket and bra.

## 1. Setup

For three independent planted pairs, write their character-moment matrices
as

$$
M_i(u,v)=\mathbb E[\overline{f_{i,u}}g_{i,v}],
\qquad i=1,2,3,
$$

where every feature has modulus one.  The four-block sector tensor is

$$
K(a,b,c,d)=M_1(a,b)M_2(b,c)M_3(c,d).
\tag{1.1}
$$

For a *whole-block* cut \(S\subseteq\{1,2,3,4\}\), flatten \(K\) into
\(K^{S\mid S^c}\).  Put arbitrary nonnegative weights \(p\) and \(q\)
on the complete row and column multi-indices of this flattening, with
total masses \(P,Q\).  In particular, neither weight is assumed to
factor over blocks.

Define

$$
\kappa_i=\max_{u,v}|M_i(u,v)|,
\qquad c_i=\|M_i\|_{\rm op},
\tag{1.2}
$$

and let \(d_2,d_3\) be the cardinalities of the second and third sector
index sets.  The two wedge constants

$$
\tau_{12}=\min\{1,\sqrt{d_2}\,\kappa_1\kappa_2\},
\qquad
\tau_{23}=\min\{1,\sqrt{d_3}\,\kappa_2\kappa_3\}
\tag{1.3}
$$
can be replaced by any sharper weighted nuclear bounds for the two-link
wedge matrices

$$
W_{12}(b;(a,c))=M_1(a,b)M_2(b,c),
\quad
W_{23}(c;(b,d))=M_2(b,c)M_3(c,d).
\tag{1.4}
$$

## 2. Complete cut table

For every cut,

$$
\|D_p^{1/2}K^{S\mid S^c}D_q^{1/2}\|_1
\le B_S\sqrt{PQ},
\tag{2.1}
$$

where \(B_S=B_{S^c}\) and it is enough to list

| cut \(S\) | \(B_S\) |
|---|---:|
| \(\varnothing\) | \(\kappa_1\kappa_2\kappa_3\) |
| \(\{1\}\) | \(\kappa_2\kappa_3\) |
| \(\{2\}\) | \(\kappa_3\tau_{12}\) |
| \(\{3\}\) | \(\kappa_1\tau_{23}\) |
| \(\{4\}\) | \(\kappa_1\kappa_2\) |
| \(\{1,2\}\) | \(\kappa_1\kappa_3\) |
| \(\{2,3\}\) | \(\kappa_2\) |
| \(\{1,3\}\) | \(c_1\kappa_2c_3\) |

No marginal or conditional version of \(p\) or \(q\) is used on two
different links.

## 3. Proofs of the nontrivial cuts

The universal weighted Gram lemma applies after arbitrary repetitions of
a feature.  It immediately proves the endpoint, adjacent, and separated
endpoint rows of the table after the links internal to one side of the cut
are absorbed into that side's diagonal weight.  For example, the
\(\{1,2\}\mid\{3,4\}\) cut leaves \(M_2(b,c)\) as a repeated Gram matrix,
while the two new total masses are at most \(\kappa_1^2P\) and
\(\kappa_3^2Q\).

For a middle singleton, the remaining two-link wedge is itself a Gram
matrix.  This gives the constant one in (1.3).  Independently, it has rank
at most \(d_2\) and entries bounded by \(\kappa_1\kappa_2\), so weighted
rank--Frobenius gives the other term in \(\tau_{12}\).  The third-block
case is identical.

For the alternating cut \(\{1,3\}\mid\{2,4\}\), introduce the latent
pair index \((b,c)\).  The weighted flattening factors as

$$
D_p^{1/2}L\,D_{M_2}\,R D_q^{1/2},
\tag{3.1}
$$

with

$$
L_{(a,c),(b',c')}=\delta_{c,c'}M_1(a,b'),
\qquad
R_{(b',c'),(b,d)}=\delta_{b',b}M_3(c',d).
$$

Therefore

$$
\|D_p^{1/2}L\|_2^2\le c_1^2P,
\qquad
\|D_{M_2}RD_q^{1/2}\|_2^2
\le\kappa_2^2c_3^2Q.
$$

Schatten Holder proves the last row of the table.  This is the place where
three separate link estimates would incorrectly charge the shared labels
more than once.

## 4. Minimal-chain recovery

For the singleton Hadamard sector,

$$
\kappa_1=\kappa_2=\kappa_3=N^{-1/2},
\qquad c_1=c_3=1,
\qquad d_2=d_3=N.
$$

The table gives exactly

$$
B_S=N^{(r_S-3)/2},
$$

where \(r_S\) is the binary cut rank of the three-edge path: \(N^{-3/2}\)
for the empty/full cuts, \(N^{-1}\) for the ten rank-one cuts, and
\(N^{-1/2}\) for the four rank-two cuts.  Thus the result independently
recovers every graph weight used in the exact one-batch dose-six KKT
certificate.

## 5. Remaining finite-size work

For a singleton in each block, every physical ket/bra placement is a
whole-block cut, so the table is exhaustive and independently validates
the minimal-chain graph factors.

For higher sectors it is not exhaustive.  If

$$
R_b=S_b\mathbin\triangle T_b
$$

contains several coordinates, some can lie in the ket support \(S_b\)
and the rest in the bra support \(T_b\).  The physical matrix is then the
Schur lift

$$
K(S_1\triangle T_1,\ldots,S_4\triangle T_4),
$$

which has eight support variables rather than a flattening that assigns
each of the four combined variables wholly to one side.  The cut table
controls the block-coherent subfamilies \(S_b=\varnothing\) or
\(T_b=\varnothing\), but not arbitrary internal splits.

Thus the composition problem has two remaining parts.  First, for the
signed-permutation plant one needs rigorous \(\kappa_i,c_i\), and wedge
constants for every compatible odd record/even-decoration profile of
total physical degree at most twelve.  Second, the proof must lift the
same once-weighted factorization to the eight-sided Schur kernel without
expanding all occurrence placements in \(\ell_1\).
The pure odd-record operator bounds already give

$$
c_i\le {q^{r_i}\over(q)_{r_i}}.
$$

Adaptive outcome-selected frames remain a separate interface after that
full one-batch sector sum is below threshold.

Reproduction: tests/three_link_weighted_path_contraction.py.
