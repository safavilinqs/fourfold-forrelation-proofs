# Adjacent double-cubic/quintic-mask insertion at q=64

Date: 2026-07-16

Status: rigorous arbitrary-diagonal one-batch coefficients for eight
four-entry degree-twelve orbits. This raises the q64 theorem count from 348
to 380. It does not prove the remaining 508 entries, intervalize the ledger,
or establish the adaptive lift.

## Result

The canonical cuts and coefficients are

| cut | quintic treatment | coefficient |
|---|---|---:|
| $(3,3,1,5):(0,1,1,4)$ | fixed-four endpoint | $4.70180143564\times10^{-5}$ |
| $(3,3,1,5):(1,0,1,4)$ | fixed-four endpoint | $0.000268011192337$ |
| $(3,3,1,5):(1,3,0,2)$ | fixed-three endpoint | $0.00312706344026$ |
| $(1,3,3,5):(0,2,0,4)$ | generic extreme mask | $0.00137086718895$ |
| $(1,3,3,5):(0,2,3,1)$ | generic extreme mask | $0.000141019599739$ |
| $(1,3,3,5):(0,3,1,2)$ | generic balanced mask | $0.0134380552033$ |
| $(1,3,3,5):(0,3,2,1)$ | generic extreme mask | $0.000803837243791$ |
| $(3,3,1,5):(0,2,0,4)$ | generic extreme mask | $0.00137086718895$ |

Each cut contains a consecutive singleton and two cubics, one whole and one
split, in either path orientation. Complement and reversal generate all 32
entries. Twenty-four have extreme $1|4$ quintic cuts and eight have balanced
$2|3$ cuts.

Honest insertion improves the routing ledger to

$$
P_{\rm total}=0.326050806446,
\qquad
{1\over3}-P_{\rm total}=0.00728252688685.
$$

The margin gain is $0.00129219857109$. There are now 48 quintic entries
left: 32 extreme and 16 balanced. Charging those entries at the local-slice
targets gives total $0.326203188868$ and leaves $0.00613014446556$ beyond
the declared $10^{-3}$ allowance.

## Inherited adjacent double-cubic contraction

The three-block subchain is $M_{33}M_{31}$ or its transpose. Either cubic
may be the endpoint or middle block; the singleton fixes the compatible
record-one support on the middle block. The inherited arbitrary-law
incidence theorem was proved for this subchain before adding any fourth
block. It therefore remains valid after composing a final physical link and
quintic mask.

For an endpoint split with $k$ cells and middle split with $\ell$ cells, put

$$
a_1={ (q+2)^2\over q^2(q-1)^2(q-2)^2}.
$$

The endpoint and middle incidence degrees are

$$
\begin{aligned}
D^E={}&\left(
q\binom q3+q^2(q-1)\binom q2,
\binom{q-1}2+(q-1)\binom q2+q(q-1)^2,
q^2-2,
1\right),\\
D^L={}&\left(q^2(q-1)^2,3(q-1)^2,2(q-1),1\right).
\end{aligned}
$$

The two row/column factorizations give

$$
\Gamma_{k\ell}(q)
\le
\min\left\{
{\sqrt{a_1D^E_kD^L_\ell}\over q-1},
{\sqrt{a_1D^E_{3-k}D^L_{3-\ell}}\over q(q-1)}
\right\}.
$$

Across the eight target orbits, every split pair except $(1,1)$ and $(2,2)$
occurs, each four times. The largest required double-cubic factor is
$\Gamma_{13}=0.00259344785064$ at $q=64$.

## Favorable endpoint orientations

Twelve entries have the singleton on the quintic majority side. Their exact
endpoint slices already include the quintic distinctness mask:

$$
F_4=1-{4\over N}=0.9990234375,
\qquad
F_3=1.45384579613.
$$

The corresponding Schur factors are $\sqrt{F_4}$ and $\sqrt{F_3}$. They
produce the first three coefficients in the table.

## Generic quintic-mask orientations

For the other 20 entries, complete the last signed-permutation link as a
unit cross Gram and then restore the quintic cross-cut distinctness mask.
For an extreme singleton--four-set cut, the centered factor is

$$
\gamma_{N,4}=2.99780260018.
$$

For a balanced pair--triple cut, the direct constant/incidence/containment
factorization gives

$$
\gamma_{2,3}\le1+\sqrt6+\sqrt3=5.18154055035.
$$

Multiplying the appropriate $\Gamma_{k\ell}$ by one of these factors gives
the last five coefficients in the table.

The cubic distinctness mask is already present in the inherited subchain
incidence theorem. The endpoint slices or generic factors insert the
quintic mask exactly once. Schur composition uses each shared physical law
only through the complete occurrence matrix, so no product-law or
invariant-law assumption is made.

## Regression scope

The regression protects all eight four-entry topologies, the eight required
split pairs, the inherited all-$q$ incidence formula including its published
$q=32$ values, the favorable/generic quintic partition, all eight q64
coefficients, the ledger insertion, and deterministic artifact output.

Reproduce with:

    /opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_adjacent_double_cubic_quintic_endpoint_insertion.py --write-artifact
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_adjacent_double_cubic_quintic_endpoint_insertion.py
