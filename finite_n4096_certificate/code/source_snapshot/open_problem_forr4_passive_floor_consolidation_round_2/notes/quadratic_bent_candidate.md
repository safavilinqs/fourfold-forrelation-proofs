# Quadratic-bent exact-plant candidate

Date: 2026-07-14

Status: viable but unresolved.  The unweighted \((5,1)\) endpoint link
norm grows as \(\Theta(\sqrt N)\), which rejects a plain operator-norm
proof.  The physically weighted nuclear norm still decays as \(O(1/N)\),
so this does not reject the hard instance.

## 1. Exact plant

Let \(n\) be even, \(N=2^n\), and sample a Boolean quadratic bent function
\(X\in\{\pm1\}^N\).  Equivalently, its polar alternating form is
nondegenerate.  Then

$$
Y=H_NX\in\{\pm1\}^N.
$$

Three independent pairs give the same exact positive plant as before,

$$
(X_1,\ Y_1X_2,\ Y_2X_3,\ Y_3),
$$

with \(F_{4,H}=1\) pointwise.  Negating the first block gives
\(F_{4,H}=-1\).  Thus this candidate has no promise conditioning loss.

The difference from the signed-permutation plant is the orbit.  It averages
over nondegenerate quadratic forms and their linear shifts rather than over
the Maiorana--McFarland permutation subclass.

## 2. Link operators

For distinct-coordinate supports \(A,B\subseteq[N]\), define

$$
M_{a,b}(A,B)
=\mathbb E\left[\prod_{i\in A}X_i\prod_{j\in B}Y_j\right],
\qquad |A|=a,\quad |B|=b.
$$

These are the complete link operators, including even-pair decorations;
they are not only the pure odd-label compounds.

At \(N=16\), exhaustive enumeration finds all 896 bent functions.  The
operator norms for \(1\le a,b\le6\) are

| \(a\backslash b\) | 1 | 2 | 3 | 4 | 5 | 6 |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 0 | 0.845154 | 0 | 2.360387 | 0 |
| 2 | 0 | 1.142857 | 0 | 1.616244 | 0 | 4.426267 |
| 3 | 0.845154 | 0 | 1.714286 | 0 | 3.833259 | 0 |
| 4 | 0 | 1.616244 | 0 | 6.857143 | 0 | 9.142857 |
| 5 | 2.360387 | 0 | 3.833259 | 0 | 13.714286 | 0 |
| 6 | 0 | 4.426267 | 0 | 9.142857 | 0 | 36.571429 |

Parity gives the checkerboard zeros.

For comparison, the equal-degree signed-permutation values at the same
\(N=16\) are

$$
1,\quad {8\over3},\quad4,\quad16,\quad32,\quad{256\over3}.
$$

The quadratic-bent diagonal is

$$
1,\quad {8\over7},\quad{12\over7},\quad{48\over7},
\quad{96\over7},\quad{256\over7}.
$$

It is strictly smaller at every degree two through six.  The off-diagonal
nonzero sectors improve as well.

## 3. Why this is relevant

The signed-permutation pure odd-record operators are close to unitary at
\(q=32\), but their full support operators acquire substantial
amplification from even-pair decorations.  The quadratic-bent orbit
averages those decorations over a much larger symmetry family.  The
\(N=16\) calculation shows that this changes the complete sector norm, not
merely an entrywise moment.

## 4. Exact obstruction to naive extrapolation

Fix the affine plane

$$
A=\{0,e_1,e_2,e_1+e_2\}\subseteq\mathbb F_2^n.
$$

Every affine 3-flat \(F\) containing \(A\) has eight points.  Put
\(B_F=F\setminus A\).  For every Boolean quadratic polynomial \(q\),
the parity of \(q\) over an affine 3-flat is zero.  Therefore

$$
\prod_{x\in A}(-1)^{q(x)}
=\prod_{x\in B_F}(-1)^{q(x)}.
\tag{4.1}
$$

There are \(N/4-1\) such 3-flat extensions, and their complements are
distinct.  Together with \(A\), this gives \(N/4\) identical columns in
the degree-four feature matrix.  At \(N=1024\), the one-side Gram operator
therefore has norm at least 256.

This does not alone lower-bound the cross operator \(M_{4,4}\), because
Fourier duality can cancel between the two identical-column classes.  It
does prove that the favorable \(N=16\) table cannot be extrapolated using
a design or near-orthogonality argument.  Any revival of this orbit must
diagonalize the cross association scheme and exhibit that cancellation
explicitly.

There is encouraging exact cancellation for the canonical plane class.
Write a nondegenerate alternating matrix in the fixed coordinate basis as
\(A\).  The plane feature of the quadratic phase is
\((-1)^{A_{12}}\), while the corresponding feature after Walsh transform is
\((-1)^{(A^{-1})_{12}}\).  Direct symplectic-form counting gives

$$
\mathbb E_A(-1)^{A_{12}+(A^{-1})_{12}}
={2\over N-2}.
\tag{4.2}
$$

One quick proof uses the Pfaffian cofactor identity
\((A^{-1})_{12}=\operatorname{Pf}(A_{\widehat1\widehat2})\).
Condition on the lower \((n-2)\)-dimensional form being nondegenerate.
The two possible values of the remaining symplectic pairing differ in
count by exactly \(2^{n-2}\); substitution of the standard count of
nondegenerate alternating forms simplifies to (4.2).

Multiplying by the collision-class size \(N/4\) gives

$$
{N\over2(N-2)}\longrightarrow{1\over2},
$$

not a growing cross norm.  Exact enumeration verifies (4.2) at
\(N=16,64\).  Thus the \(N/4\) collision is a barrier to a one-side Gram
proof, but not a counterexample to the cross-sector candidate.

## 5. \((5,1)\) growth and the required physical weighting

The next endpoint sector can also be diagonalized exactly.  For
\(|A|=5\) and one \(Y\)-coordinate, expand

$$
Y_y=N^{-1/2}\sum_z(-1)^{y\cdot z}X_z.
$$

The linear shift average forces \(z=\bigoplus A\).  Rows with the same XOR
are therefore scalar multiples of one orthonormal Walsh row.  The squared
operator norm is the sum of the squared symplectic Fourier coefficients in
one XOR class.

There are exactly two configuration types.  In a fixed \(d\)-dimensional
ambient vector space:

- four nonzero, zero-XOR points spanning dimension three occur seven times
  per 3-subspace and have alternating rank two;
- five nonzero, zero-XOR points spanning dimension four occur 168 times
  per 4-subspace and have alternating rank four.

If

$$
\phi_r={(-1)^r\over
\prod_{j=1}^r(2^{n-2j+1}-1)}
$$

is the Fourier coefficient of the uniform nondegenerate alternating-form
orbit at rank \(2r\), then

$$
\boxed{
\|M_{5,1}\|_{\rm op}^2
=7{n\brack3}_2\phi_1^2
+168{n\brack4}_2\phi_2^2.
}
\tag{5.1}
$$

At \(N=16\), (5.1) is \(39/7\), reproducing the exhaustive value
\(2.360387\).  At \(N=1024\),

$$
\|M_{5,1}\|_{\rm op}^2
={11182413\over64897},
\qquad
\|M_{5,1}\|_{\rm op}\approx13.126697.
\tag{5.2}
$$

The first term in (5.1) is \(\Theta(N)\), so the unweighted norm is
\(\Theta(\sqrt N)\).  This rules out a proof that simply multiplies
unweighted link operator norms.

It does **not** reject the passive hard instance.  Every row has the more
specific form

$$
M_{5,1}(A,y)=\mu_AH_{\oplus A,y},
\qquad |\mu_A|\le {2\over N-2}.
\tag{5.3}
$$

For arbitrary diagonal probe weights \(p_A,q_y\), with total masses
\(P=\sum_Ap_A\) and \(Q=\sum_yq_y\),

$$
\begin{aligned}
\|D_p^{1/2}M_{5,1}D_q^{1/2}\|_F^2
&\le {1\over N}\left({2\over N-2}\right)^2PQ,\\
\operatorname{rank}M_{5,1}&\le N.
\end{aligned}
$$

Therefore

$$
\boxed{
\|D_p^{1/2}M_{5,1}D_q^{1/2}\|_1
\le {2\over N-2}\sqrt{PQ}.
}
\tag{5.4}
$$

The identical-row multiplicity that creates (5.2) cancels against the
physical diagonal normalization.  At \(N=1024\), the coefficient in (5.4)
is \(1/511\).

## 6. Disposition

The finite \(N=16\) improvement was real, while (5.2) shows that
unweighted operator spectra are the wrong extrapolation.  Equation (5.4)
restores the relevant passive scale for the first growing sector.

The next task is to establish analogous weighted nuclear decompositions
for every \(M_{a,b}\) with \(a+b\le12\), then compose three links through
the product middle blocks and adaptive complete frames.  The
quadratic-bent orbit remains a candidate until that weighted program
succeeds or produces a weighted counterexample.

Reproduction: searches/quadratic_bent_sector_spectra.py and
tests/quadratic_bent_collision_barrier.py.
