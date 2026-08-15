# Signed-permutation cubic decorations and path compatibility

Date: 2026-07-14

Status: exact degree-one/degree-three link formulas proved.  Cubic even
decorations cause a growing isolated link norm, but middle-block path
compatibility removes that growth.  A cubic endpoint orbit still attains
the minimal \(N^{-1/2}\) scale, so higher sectors cannot be discarded by
claiming an extra power of \(N\).

## 1. Compressed support data

Write \(q=2^m\), \(N=q^2\), and index a physical coordinate by
\((x,y)\in\mathbb F_2^m\times\mathbb F_2^m\).  For a signed permutation
\(P\), the exact pair is

$$
X=K_qP,\qquad Y=PK_q,
$$

where \(K_q\) is the unnormalized Walsh sign matrix.  An \(X\)-support
compresses each column label \(y\) to its parity and the XOR of its
incident \(x\)'s.  A \(Y\)-support has the transposed row-label data.
The sign average matches the odd column labels to the odd row labels;
even labels remain in the permutation phase.

This is why even decorations are invisible to the scalar matching event
but visible to the spectral link norm.

## 2. Exact \(M_{1,3}\) decomposition

Fix one \(X\)-coordinate \((x,y)\).  A three-coordinate \(Y\)-support
with one odd row has two types.

In type I, all three coordinates occupy one row \(u\).  For distinct
\(v_1,v_2,v_3\),

$$
M_{1,3}((x,y);B)
={1\over q}(-1)^{x\cdot u+y\cdot(v_1+v_2+v_3)}.
\tag{2.1}
$$

In type II, one coordinate occupies the odd row \(u\) and two occupy a
different even row.  Then

$$
|M_{1,3}((x,y);B)|={1\over q(q-1)}.
\tag{2.2}
$$

The character columns in each type are repeated normalized Walsh vectors.
Their exact frame operators are

$$
M_{m I}M_{m I}^*
={(q-1)(q-2)\over6}I_N,
\qquad
M_{\rm II}M_{\rm II}^*={q\over2}I_N.
\tag{2.3}
$$

Consequently

$$
\boxed{
M_{1,3}M_{1,3}^*={q^2+2\over6}I_N,
\qquad
\|M_{1,3}\|_{\rm op}=\sqrt{{q^2+2\over6}}.
}
\tag{2.4}
$$

At \(q=32\), the squared norm is exactly \(171\).  Thus the isolated
decorated link really does grow as \(\Theta(\sqrt N)\).

## 3. A cubic middle block

If the same degree-three physical support is used on the right side of
one link and the left side of the next, record size one on both sides
forces one odd row and one odd column.  The support is exactly a
three-edge \(2\times2\) L-shape.  There are

$$
q^2(q-1)^2
\tag{3.1}
$$

such supports.  For each odd-row/odd-column pair, the remaining even row
and column give \((q-1)^2\) copies of one Walsh character.  After the
amplitude \(1/[q(q-1)]\) in (2.2) is included, the restricted frame
operator is the identity:

$$
\boxed{
\|M_{1,3}^{\rm L}\|_{\rm op}
=\|M_{3,1}^{\rm L}\|_{\rm op}=1,
\qquad
\max|M_{1,3}^{\rm L}|={1\over q(q-1)}.
}
\tag{3.2}
$$

At \(q=32\), the coherence is \(1/992\).  In every block-coherent cut of
the \((1,3,1,1)\) path profile, the weighted path theorem therefore gains
at least \(1/[q(q-1)]\), rather than losing the isolated norm in (2.4).

## 4. A cubic endpoint

An endpoint is constrained by only one neighboring link.  Its type-I
orbit consists of three coordinates in one hidden column.  Conditional
on this orbit, the decorated endpoint matrix is a repeated copy of the
minimal Walsh link.  In the alternating \((3,1,1,1)\) path flattening,
uniform physical weights on the two transitive sides give

$$
\boxed{
{\|K_{3,1,1,1}^{13\mid24}\|_1
 \over\sqrt{(\text{row dimension})(\text{column dimension})}}
={1\over q}={1\over\sqrt N}.
}
\tag{4.1}
$$

Thus the growing unweighted link norm does not create a constant-size
passive distinguisher.  It also does not buy an extra \(1/q\): cubic
endpoint decorations can be as dangerous in dimension as the minimal
chain.  For the attenuated plant their only automatic improvement is the
additional factor \(\beta^2\) from total Fourier degree six.

## 5. Exact occurrence-split cubic lift

The most balanced internal split puts one of the three endpoint marks on
the ket, two on the bra, two singleton chain blocks on the ket, and the
last singleton on the bra.  Rows are indexed by \((a,b,c)\); columns by
\((\{e,f\},d)\).  The link moment has the unified form

$$
M_{3,1}(\{a,e,f\},b)
=w(a,e,f)H_N(a\oplus e\oplus f,b),
\tag{5.1}
$$

where \(w=1\) when the three hidden column labels agree,
\(w=-1/(q-1)\) when exactly two agree, and \(w=0\) when all three are
distinct.

In the column Gram matrix, orthogonality of the last two Walsh links
leaves unordered endpoint pairs only.  Split those pairs by their nonzero
XOR \(r\).  If the hidden-label component of \(r\) is zero, the singular
values before the final normalization are

$$
{\sqrt2\over q},qquad
{\sqrt2(q^2-2q+2)\over2q(q-1)}.
$$

If it is nonzero, they are

$$
{\sqrt2\over q},qquad
{\sqrt2\over q(q-1)}.
$$

Their multiplicities sum exactly to the unordered-pair dimension.  The
normalized nuclear norm is

$$
\boxed{
{4(q-1)(q^2-q+1)\over q^4\sqrt{q^2-1}}.
}
\tag{5.2}
$$

It is \(0.15733994844\) at \(q=4\), matching direct diagonalization of all
32 occurrence placements, and

$$
0.003671413\ldots
$$

at \(q=32\).  Its asymptotic scale is \(4/q^2=4/N\), strictly smaller
than the whole-block endpoint value \(1/q\).

## 6. Scope

Equations (3.2), (4.1), and (5.2) now cover the block-coherent cubic path
and every direct \(q=4\) endpoint occurrence split, with the worst balanced
split proved for general \(q\).  A complete degree-six upper bound still
needs general diagonal weights rather than uniform orbit weights, the
cubic middle-block Schur lifts, and simultaneous summation with the
minimal sector.  Any proposed higher-sector theorem must reproduce (4.1),
not claim uniform \(N^{-1}\) decay for every cubic placement.

Reproduction: searches/signed_permutation_full_sector_spectra.py and
searches/signed_permutation_cubic_schur_lift.py.
