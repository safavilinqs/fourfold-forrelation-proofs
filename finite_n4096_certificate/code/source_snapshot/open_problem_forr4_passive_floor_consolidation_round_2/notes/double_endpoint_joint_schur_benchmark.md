# Joint double-endpoint Schur benchmark

Date: 2026-07-14

Status: one balanced same-side occurrence slice is diagonalized for all
\(q\), and the alternating slice is diagonalized at \(q=2,4\).  These are
uniform orbit-weight benchmarks, not arbitrary-diagonal contraction
theorems.  They show that the large splitwise row-energy bound destroys a
real tight-frame cancellation.

## 1. Endpoint pair frame

For an endpoint singleton \(i\), an unordered disjoint pair \(E\), and
the neighboring singleton \(b\), put

$$
A(i;E,b)=M_{3,1}(\{i\}\cup E,b).
$$

The exact endpoint slice counts imply the tight-frame identity

$$
\boxed{
AA^*=\lambda_q I_N,
\qquad \lambda_q={q^2+2\over2}.
}
\tag{1.1}
$$

The diagonal is \(q^2E_1=(q^2+2)/2\), using the exact one-endpoint slice
energy.  For two different singled-out coordinates, Walsh orthogonality
cancels the common pair completions, giving the zero off-diagonal.
Direct signed-permutation enumeration verifies (1.1) at \(q=2,4\).

## 2. Both middle singletons on one side

Split both cubic endpoints as one coordinate versus an unordered pair,
and put the two middle singleton coordinates on the pair side.  The full
matrix is

$$
K((i,d),(E,b,c,F))
=A(i;E,b)H_N(b,c)A(d;F,c).
\tag{2.1}
$$

Because every squared Hadamard entry equals \(1/N\), equation (1.1)
gives

$$
KK^*={\lambda_q^2\over N}I_{N^2}.
\tag{2.2}
$$

All \(N^2\) nonzero singular values are \(\lambda_q/q\).  With uniform
row and column masses, the normalized nuclear coefficient is therefore

$$
\boxed{
\gamma_{\mathrm{same}}(q)
={q^2+2\over q^3(q^2-1)}.
}
\tag{2.3}
$$

It equals \(1/4\) at \(q=2\), \(0.01875\) at \(q=4\), and
\(3.06071\times10^{-5}\) at \(q=32\).  The last number is more than three
orders of magnitude below the worst arbitrary-diagonal fixed-split
coefficient.

## 3. Alternating middle singletons

Putting the first middle singleton on the row side and the second on the
column side gives

$$
K((i,b,d),(E,c,F))
=A(i;E,b)H_N(b,c)A(d;F,c).
\tag{3.1}
$$

Its row Gram can be assembled without forming the wide matrix.  In the
\((b,b')\) block it is the Schur product of the endpoint pair Gram and the
Hadamard-twisted endpoint pair Gram.  Exact diagonalization gives

| \(q\) | normalized nuclear coefficient | rank |
|---:|---:|---:|
| 2 | 0.471591815891 | 64 |
| 4 | 0.064200871625 | 4096 |

The alternating slice is the larger uniform benchmark, but it falls by a
factor \(7.35\) from \(q=2\) to \(q=4\), already much faster than the
splitwise \(q^{-1/2}\) estimate.  Two data points do not justify an
asymptotic formula; no \(q=32\) value is claimed.

## 4. Next theorem target

The uniform calculation falsifies the idea that the
\((3,1,1,3)\) sector itself has a large finite-size constant.  The loss
comes from allowing every split to choose unrelated adversarial diagonal
weights and then summing those bounds.

The needed result is a weighted compound-frame inequality for (3.1) that
keeps the two endpoint pair frames and their shared passive occupation
law together.  It must improve the current splitwise occupation sum by a
factor below \(0.280708\).  A proof based only on maximum row and column
energies cannot see (1.1)--(2.2) and is therefore unlikely to suffice.

Reproduction: `searches/double_endpoint_joint_schur_benchmark.py`.
