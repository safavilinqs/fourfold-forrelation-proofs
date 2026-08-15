# Decorated adjacent complete-row insertion at q=64

Date: 2026-07-16

Status: rigorous arbitrary-diagonal one-batch coefficients for 16 entries in four complement/reversal orbits. This raises the q64 theorem count from 276 to 292. It does not prove the remaining 596 entries, intervalize the ledger, or establish the adaptive lift.

## Result

Sixteen remaining split-cubic/split-quintic entries contain three consecutive blocks of degrees $(1,3,5)$ or $(5,3,1)$. On the side containing the singleton, their fixed support counts are either one cubic cell and four quintic cells or two cubic cells and three quintic cells. The fourth block is an unsplit cubic at the outer endpoint.

The four canonical cuts are

$$
\begin{aligned}
(1,3,5,3)&:(0,1,2,3),\\
(1,3,5,3)&:(0,2,1,3),\\
(3,1,3,5)&:(0,1,1,4),\\
(3,1,3,5)&:(0,1,2,3).
\end{aligned}
$$

Complement and reversal generate all 16 entries. Eight have extreme quintic splits and eight have balanced quintic splits.

Every entry has the common arbitrary-law coefficient

$$
\boxed{c_{\rm row}=0.0200795672469}.
$$

After insertion, the routing values are

$$
\begin{aligned}
\beta&=0.746098905399,\\
P_{\rm Perron}&=0.311849240525,\\
P_{\rm promise}&=0.0172579672069,\\
P_{\rm total}&=0.329107207732.
\end{aligned}
$$

The routing margin is $0.00422612560146$, an improvement of $0.000776398242664$ over the 276-entry ledger. The remaining quintic inventory is 136 entries: 96 extreme and 40 balanced. Charging those remaining entries at the existing endpoint local-slice coefficients gives total $0.330076650865$ and leaves $0.00225668246805$ beyond the declared $10^{-3}$ allowance.

## Why the inherited complete-row theorem applies

The Round 3 adjacent-row theorem treats a fixed singleton, one fixed cubic cell, and three fixed quintic cells as the row variables. If $E$ is the complementary cubic pair and $G$ the complementary quintic pair, its joint link feature is

$$
L=M_{13}M_{35},
$$

and the complete row energy is bounded at $q=64$ by

$$
R_{64}\le 0.000403189020824.
$$

The arbitrary-law trace-class coefficient is therefore at most $\sqrt{R_{64}}=0.0200795672469$.

The new cuts fix additional cells on the row side. In the $(1,4)$ case, choose any three of the four fixed quintic cells as the old fixed triple; the fourth occupies one of the old complementary coordinates. In the $(2,3)$ case, choose either of the two fixed cubic cells as the old fixed cell; the other occupies one of the old complementary coordinates. In both cases the new complete row is an injectively embedded subrow of an old complete row. Deleting its other coordinates cannot increase squared row energy, so the same $R_{64}$ applies.

The fourth cubic block is unsplit. Its only new link is a signed-permutation moment that separates into row and column characters, hence is a cross Gram of unit vectors. Its trace-class Schur-multiplier norm is at most one. The cubic and quintic within-block distinctness masks are already retained inside the complete $M_{13}M_{35}$ row, while the unsplit outer cubic introduces no cross-cut within-block mask. Schur multiplication by the outer link therefore cannot increase the inherited coefficient.

This argument is uniform over arbitrary correlated diagonal laws. It does not assume product marginals or use the physical lower-witness substitutions.

## Updated acceptance interface

The lead shared-contraction family now has 136 entries. A common coefficient at most $0.151320540139$ is sufficient to retain a $10^{-3}$ one-batch reserve. If the 96 extreme entries are charged at their existing local-slice coefficient $0.123974636390$, the remaining 40 balanced entries may have a common coefficient as large as $0.205605732426$ while retaining that reserve.

These are routing acceptance thresholds, not new contraction theorems. The adaptive additive-overhead allowance under the local two-tier proxy rises to $0.00225668246805$ after reserving $10^{-3}$ for certification and implementation uncertainty.

Reproduce with:

    /opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_decorated_adjacent_row_insertion.py --write-artifact
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_decorated_adjacent_row_insertion.py
