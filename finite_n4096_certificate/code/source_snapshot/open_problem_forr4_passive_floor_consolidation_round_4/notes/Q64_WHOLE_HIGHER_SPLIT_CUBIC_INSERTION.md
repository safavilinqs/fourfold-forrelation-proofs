# q64 whole-higher split-cubic insertion

Date: 2026-07-16

Status: the 48-entry mask-aware theorem is preserved, but its upstream
664-entry baseline is invalid. The cumulative 712-entry count, routing total,
and statement that only 176 entries remain are withdrawn. See
`Q64_MASKED_UNIVERSAL_AUDIT.md`.

## Result

The 48 entries split into three finite templates:

| template | entries | coefficient |
|---|---:|---:|
| favorable singleton--split-cubic endpoint | 24 | $0.124035215254$ |
| internal singleton--whole-cubic endpoint | 16 | $0.015625$ |
| complete singleton--quintic--whole-cubic wedge | 8 | $0.00846466875312$ |

Inserting these coefficients after the noncubic/recovered-universal theorem gives

$$
U_{\mathrm{route}}=0.328477421166173,
$$

with raw diagnostic margin

$$
{1\over3}-U_{\mathrm{route}}
=0.004855912167160.
$$

Retaining the declared $10^{-3}$ numerical allowance leaves conditional adaptive additive cap

$$
0.003855912167160,
$$

or multiplicative cap

$$
1.01173874342252.
$$

The remaining 176 entries retain their frozen target $0.124035215254$. Their regenerated common reserve gate is

$$
0.142581909211.
$$

## Favorable cubic endpoint

In 24 entries, the singleton lies on the same side as the pair portion of the split cubic. Complete all other physical links as unit cross-Gram Schur multipliers. The exact cubic endpoint fixed-pair squared slice is

$$
E_2={q^2-2q+2\over q^2(q-1)}.
$$

The arbitrary-law Schur feature therefore gives

$$
\gamma_{\mathrm{fav}}\le\sqrt{E_2}
=0.124035215254
$$

at $q=64$.

## Internal singleton--whole-cubic endpoint

In 16 entries, a singleton and an adjacent whole cubic lie entirely on the same physical side of the cut. Their moment is a row-only or column-only scalar. The shared record is one, and the cubic's other record is one or three, so

$$
|M_{1,3}|\le {1\over q}.
$$

The remaining links form a unit cross-Gram dressing. Thus

$$
\gamma_{\mathrm{int}}\le {1\over q}=0.015625.
$$

## Complete singleton--quintic--cubic wedge

The last eight entries contain a whole singleton, whole quintic, and whole cubic on the alternating sides of a two-link wedge. The remaining split-cubic link is a unit cross-Gram dressing and can be restored after the wedge contraction.

Write

$$
T_{Q,(S,C)}
=M_{1,5}(S,Q)M_{5,3}(Q,C).
$$

For arbitrary nonnegative row weights $p_{S,C}$ and column weights $r_Q$, decompose the weighted matrix into its rank-one columns:

$$
\left\|D_p^{1/2}TD_r^{1/2}\right\|_1
\le
\sum_Q\sqrt{r_Q}\,
\left\|D_p^{1/2}T_{\cdot,Q}\right\|_2.
$$

Cauchy--Schwarz gives

$$
\left\|D_p^{1/2}TD_r^{1/2}\right\|_1
\le
\sqrt{\sum_Qr_Q}\,
\sqrt{
\sum_{S,C}p_{S,C}
\sum_Q|M_{1,5}(S,Q)M_{5,3}(Q,C)|^2
}.
$$

The complete endpoint quintic squared energies, split by its record on the cubic link, are obtained from the fixed-one energies by cell transitivity:

$$
(F_{0,1},F_{0,3},F_{0,5})
=
{q^2\over5}(F_{1,1},F_{1,3},F_{1,5}).
$$

At $q=64$,

$$
(F_{0,1},F_{0,3},F_{0,5})
=(1027,\ 45690,\ 793168).
$$

The record-five sector cannot couple to a cubic. With

$$
m_1={q+2\over q(q-1)(q-2)},
\qquad
m_3={3\over(q-3)\binom q3},
$$

the complete row energy is at most

$$
F_{0,1}m_1^2+F_{0,3}m_3^2.
$$

Therefore

$$
\gamma_{\mathrm{wedge}}
\le
\sqrt{F_{0,1}m_1^2+F_{0,3}m_3^2}
=0.00846466875312.
$$

This derivation uses the same row and column law throughout. It is not a product-law or separate-link optimization.

## Reproduction

Run:

    /opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_whole_higher_split_cubic_insertion.py --write-artifact
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_whole_higher_split_cubic_insertion.py

The regression checks the exact 24/16/8 template partition, endpoint slice identities, direct $q=4$ signed-permutation wedge energy, every q64 coefficient, the remaining class gate, and byte-for-byte artifact regeneration.
