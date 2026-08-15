# Counterexample to the near-unit two-copy square function

Date: 2026-07-14

Status: exact rejection of the dimension-uniform constant required by the
first joint-probe proposal.

## 1. Two valid passive nodes

At \(N=1\), write the four block signs as \(x_1,\ldots,x_4\).  At each of
two nonadaptive dose-one nodes, use the probe law

$$
q(\{b\})={1\over4},\qquad b=1,\ldots,4,
$$

with no vacuum mass, and measure in the real \(4\times4\) Hadamard basis
\(K_4/2\).  The complete-frame coefficients are \(K_4/4\), so outcome
\(y\) has probability

$$
p_y(x)={1\over16}
\left(\sum_{b=1}^4K_{4,yb}x_b\right)^2.
\tag{1.1}
$$

Completeness is immediate from orthogonality of \(K_4/2\).

## 2. Exact minimal-sector output

Pair the two-node transcript with the minimal character
\(x_1x_2x_3x_4\).  Direct Walsh expansion gives the \(4\times4\) signed
output matrix

$$
d_{yz}
=\mathbb E_x[x_1x_2x_3x_4p_y(x)p_z(x)]
=
\begin{cases}
3/32,&y=z,\\
-1/32,&y\ne z.
\end{cases}
\tag{2.1}
$$

An adversarial terminal sign therefore extracts

$$
\sum_{y,z}|d_{yz}|={3\over4}.
\tag{2.2}
$$

This is a nonadaptive example; outcome-selected descendants are not
responsible for the failure.

## 3. Comparison with joint-probe square mass

At either node, a two-block union has local square mass

$$
g_h[\{b,c\}]
=2q(\{b\})q(\{c\})={1\over8}.
$$

Across two nodes, the six ordered complementary two-block masks give

$$
\mathcal S=6\left({1\over8}\right)^2={3\over32}.
\tag{3.1}
$$

Thus the claimed constant-one estimate would compare \(3/4\) with
\(\sqrt{3/32}\), but the exact ratio is

$$
{3/4\over\sqrt{3/32}}=\sqrt6.
\tag{3.2}
$$

The same physical construction embeds into \(N=2\) by using coordinate
zero in every block.  The selected Hadamard-chain entry contributes
\(2^{-3/2}\), while the proposed common graph factor is \(2^{-1/2}\).
The required squared constant becomes

$$
{6\over2^2}={3\over2},
$$

so the required constant is \(\sqrt{3/2}\approx1.224745\), strictly larger
than

$$
{32\over27}\approx1.185185.
$$

## 4. Consequence

The occupation identity

$$
\mathcal S\le(D/2)^4
$$

remains correct.  What fails is the inference from that square mass to the
absolute terminal transcript sum with a near-unit, dimension-uniform
constant.  Complete outcomes can coherently mix the six complementary
block placements.

This small-\(N\) witness does not by itself rule out an
\(N=1024\)-specific inequality that retains exact chain entries, cut ranks,
or representation sectors.  It does rule out using the unweighted
joint-probe square mass as the missing universal lemma, and invalidates the
hypothetical minimal-sector budget based on that lemma.

Reproduction: tests/two_copy_square_function_counterexample.py.
