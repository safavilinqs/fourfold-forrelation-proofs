# The $q=64$ masked four-cubic incidence repair

Date: 2026-07-16

Status: arbitrary-correlated-diagonal one-batch theorem for all 38 residual
$(3,3,3,3)$ entries. Together with the previous repairs, 290 of the 354
affected entries are proved and 64 remain quarantined.

## Record-sector decomposition

On a link between two cubic supports, the common odd record has size one or
three. Split the physical four-block occurrence matrix into the eight disjoint
record triples

$$
\rho=(\rho_1,\rho_2,\rho_3)\in\{1,3\}^3.
$$

The exact signed-permutation permanent formula gives the one-link bounds

$$
a_q=\frac{q+2}{q(q-1)(q-2)}
\quad\text{for record one},
$$

and

$$
b_q=\frac{1}{\binom{q}{3}}
\quad\text{for record three}.
$$

For record three, this is the bound of six on the $3\times3$ odd-block
permanent divided by $(q)_3$. For record one, the exact special-row/special-
column expansion of the even block gives the displayed numerator $q+2$.
These are analytic formulas at $q=64$, not extrapolations from a small-order
screen.

## The masks stay inside the completion degrees

For a cubic family $F$, let $D_F(k)$ be the maximum number of supports in
$F$ that contain a fixed $k$-cell partial support. The five families needed
by a four-cubic chain are:

| family | records seen by the block | $(D_F(0),D_F(1),D_F(2),D_F(3))$ |
|---|---|---|
| $E_1$ | endpoint record one | $(q\binom{q}{3}+q^2(q-1)\binom{q}{2},\binom{q-1}{2}+(q-1)\binom{q}{2}+q(q-1)^2,q^2-2,1)$ |
| $E_3$ | endpoint record three | $(\binom{q}{3}q^3,\binom{q-1}{2}q^2,q(q-2),1)$ |
| $L_{11}$ | middle records $(1,1)$ | $(q^2(q-1)^2,3(q-1)^2,2(q-1),1)$ |
| $S_{13}$ | middle records $(1,3)$ or $(3,1)$ | $(q\binom{q}{3}(3q-2),(q-1)(q-2)(3q-2)/2,q(q-2),1)$ |
| $M_{33}$ | middle records $(3,3)$ | $(6\binom{q}{3}^2,2\binom{q-1}{2}^2,(q-2)^2,1)$ |

Fix a row occurrence partial support of size $s_i$ in block $i$. Every
physical column completion is uniquely the complement of that partial support
inside a three-cell support from the relevant family. It is therefore
disjoint by construction. The number of possible completions in record sector
$\rho$ is at most

$$
R_\rho(s)=\prod_{i=1}^4D_{F_i(\rho)}(s_i).
$$

The analogous column degree is

$$
C_\rho(s)=\prod_{i=1}^4D_{F_i(\rho)}(3-s_i).
$$

No completed unmasked kernel is introduced. The cross-cut distinctness mask
is already present in these completion degrees.

## Arbitrary-law coefficient

Let

$$
m_\rho=\prod_{j=1}^3
\begin{cases}
a_q,&\rho_j=1,\\
b_q,&\rho_j=3.
\end{cases}
$$

Factoring a sector through its physical rows gives coefficient at most
$m_\rho\sqrt{R_\rho(s)}$; transposing gives
$m_\rho\sqrt{C_\rho(s)}$. Thus

$$
c_s\le
\sum_{\rho\in\{1,3\}^3}
m_\rho\sqrt{\min\{R_\rho(s),C_\rho(s)\}}.
$$

This is a row/column feature factorization, so it is uniform over arbitrary
correlated diagonal row and column laws. At $q=64$, all 38 residual
four-cubic entries satisfy

$$
0.00224319523802
\le c_s\le
0.00894228260682<1.
$$

The worst split is $(1,3,0,2)$ and its complement/reversal orbit.

## Exact regression

The regression enumerates all 560 cubic supports at $q=4$ and all 41,664 at
$q=8$, checking the five completion-incidence formulas for every fixed subset
size. At $q=4$ it also evaluates all 157,952 compatible cubic--cubic link
pairs with exact rational moments. The observed record-one and record-three
maxima equal $a_4=b_4=1/4$. Two embedded $q=8$ extremizers attain
$a_8=5/168$ and $b_8=1/56$. The test also checks the 38-entry cover, orbit
closure, outward rounding of every one of the 304 sector coefficients, and
byte-for-byte artifact regeneration.

Reproduce with:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_masked_four_cubic_incidence_repair.py --output artifacts/q64_masked_four_cubic_incidence_repair.json
/opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_masked_four_cubic_incidence_repair.py
```
