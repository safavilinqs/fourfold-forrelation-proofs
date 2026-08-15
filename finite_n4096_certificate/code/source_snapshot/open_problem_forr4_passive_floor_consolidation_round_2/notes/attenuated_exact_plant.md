# Attenuated exact plant

Date: 2026-07-14

Status: exact variance lemma proved; promising finite-size candidate, not yet a
passive transcript theorem.

## 1. Construction

Let \(H=H_N\) be the normalized real Hadamard matrix, and take three exact
Boolean pairs

$$
y_r=Hx_r,\qquad x_r,y_r\in\{\pm1\}^N.
$$

The signed-permutation construction supplies such pairs when \(N=q^2\).
Its positive four-chain input is

$$
X=(x_1,\ y_1x_2,\ y_2x_3,\ y_3),
$$

and has \(F_{4,H}(X)=1\).  Negating its first block gives the negative
plant with value \(-1\).

Independently multiply every one of the \(4N\) coordinates by a sign
\(\eta\) with

$$
\mathbb E\eta=\beta.
$$

Call the resulting input \(\widetilde X\).  Multilinearity immediately
gives

$$
\mathbb E[F_{4,H}(\widetilde X)\mid X]
=\pm\beta^4.
\tag{1.1}
$$

Low-degree hypothesis-dependent moments are simultaneously attenuated by
\(\beta^v\) in degree \(v\).

## 2. Exact derivative energy

For a nonempty set \(R\subseteq[4]\) of open blocks, let \(D_R F(X)\) be
the tensor of multilinear derivatives indexed by one coordinate in every
block of \(R\).  Then every exact plant above obeys

$$
\boxed{\|D_R F(X)\|_2^2={1\over N}}
\qquad(R\ne\varnothing).
\tag{2.1}
$$

Square the derivative tensor and sum its open coordinate labels.  At the
first open vertex, viewed from the left, this applies coordinate
dephasing.  If earlier vertices were closed, the exact relations
\(Hx_r=y_r\) propagate a Boolean vector to that vertex; its dephasing is
the identity matrix.  If the first vertex itself is the left endpoint,
summing its open label also starts with the identity.  Every later
Hadamard conjugation, planted diagonal sign, and coordinate dephasing
preserves the identity.  The right boundary therefore evaluates to \(N\).
The two copies of the \(1/N\) normalization in \(F_{4,H}\) give
\(N/N^2=1/N\).

This proof uses only that \(H\) is unitary and the exact trajectory is
Boolean.  It does not rely on averaging over signed permutations.

## 3. Exact conditional variance

Write

$$
\eta=\beta+\sqrt{1-\beta^2}\,\zeta,
$$

where the independent variables \(\zeta\) are centered and have unit
variance.  The orthogonal multilinear expansion, grouped by its nonempty
set \(R\) of selected blocks, and (2.1) give

$$
\begin{aligned}
\operatorname{Var}(F_{4,H}(\widetilde X)\mid X)
&=\sum_{\varnothing\ne R\subseteq[4]}
\beta^{2(4-|R|)}(1-\beta^2)^{|R|}
\|D_RF(X)\|_2^2\\
&={1-\beta^8\over N}.
\end{aligned}
\tag{3.1}
$$

Both the conditional mean and conditional variance are independent of the
latent exact plant.  Hence (1.1) and (3.1) remain exact after mixing the
signed permutations.

## 4. A rational \(N=1024\) promise calculation

Choose \(\beta=5/6\).  Then

$$
\mu=\beta^4={625\over1296}\approx0.482253,
\qquad
\sigma^2={1-\beta^8\over1024}
\approx7.49445\times10^{-4}.
$$

The distance from the promise boundary \(1/4\) is \(a=301/1296\).
Cantelli's inequality bounds the bad-promise probability of either sign by

$$
\epsilon\le {\sigma^2\over\sigma^2+a^2}
={1288991\over94064415}
\approx0.0137033.
\tag{4.1}
$$

Conditioning both hypotheses onto the promise therefore changes their
transcript distance by at most \(2\epsilon\approx0.0274066\).

## 5. Rejected partial transcript budget

Combining the occupation value \(9\), the minimal attenuation \(\beta^4\),
and a hypothetical constant-one two-copy square function would have given

$$
{9\beta^4\over32}+2\epsilon
\approx0.163040<{1\over3}.
\tag{5.1}
$$

That inference is invalid.  The exact witness in
notes/two_copy_square_function_counterexample.md disproves the required
near-unit, dimension-uniform frame lemma.  Equation (5.1) is retained only
as a rejected budget, not as evidence for a lower bound.

The attenuation construction still supplies two valid assets: exact
degree-\(v\) damping by \(\beta^v\), and the exact promise-loss calculation
(4.1).  A viable proof must combine them with a contraction that retains
exact chain entries, cut ranks, or signed-permutation representation
sectors, and must also bound all higher odd sectors.

Reproduction: tests/attenuated_exact_plant_variance.py and
tests/two_copy_square_function_counterexample.py.
