# The $q=64$ block-coherent contraction

Date: 2026-07-16

Status: rigorous arbitrary-diagonal one-batch coefficients for all 70 block-coherent balanced entries in the open degree-ten/twelve profiles. This closes the first part of the $N=4096$ high-sector worklist. It does not close the other 818 entries or the adaptive lift.

## Statement

Let $a=(a_0,a_1,a_2,a_3)$ be one of the open odd profiles of total degree ten or twelve, and let $s$ be a balanced occurrence split with

$$
2|s|=|a|
$$

such that every block lies wholly on one physical side:

$$
s_j\in\{0,a_j\}.
$$

For every arbitrary correlated diagonal probe law, the weighted trace norm of this profile-split entry is at most

$$
C_{a,s}(q)\sqrt{PQ},
$$

where $C_{a,s}(q)$ is the exact rational record-sector sum defined below. At $q=64$, all 70 coefficients obey

$$
\frac1{4096}
\le C_{a,s}(64)
\le
\frac{2609304163}{39728800944}
=0.0656778986780.
$$

The maximum occurs for profile $(3,3,3,3)$ and alternating whole-block cut $(3,0,3,0)$, together with its complement.

## Record-sector calculation

For a compatible odd-record triple $r=(r_1,r_2,r_3)$, define

$$
\kappa_i=\binom q{r_i}^{-1},
\qquad
c_i=\frac{q^{r_i}}{(q)_{r_i}}.
$$

The first quantity is the exact maximum signed-permutation moment entry in that record sector. The second is the accepted pure-record operator bound.

Because the cut is block coherent, no within-block disjointness mask is required. The accepted weighted three-link path theorem applies directly. Its cut table gives:

| whole-block cut, up to complement | sector coefficient |
|---|---:|
| $\varnothing$ | $\kappa_1\kappa_2\kappa_3$ |
| $\{0\}$ | $\kappa_2\kappa_3$ |
| $\{1\}$ | at most $\kappa_3$ |
| $\{2\}$ | at most $\kappa_1$ |
| $\{3\}$ | $\kappa_1\kappa_2$ |
| $\{0,1\}$ | $\kappa_1\kappa_3$ |
| $\{1,2\}$ | $\kappa_2$ |
| $\{0,2\}$ | $c_1\kappa_2c_3$ |

The middle-singleton wedge factors in the fourth and fifth rows are bounded by one. Summing this bound over every compatible record triple gives $C_{a,s}(q)$. There are 196 record-sector bounds across the 70 physical entries. Every number is computed as an exact rational before insertion into the floating Perron ledger.

## Finite-size effect

The earlier sufficient two-tier routing target assigned $0.124035215254$ to every cubic-containing entry and $1/2$ to every other entry. Its total was

$$
0.319181162161,
$$

with margin $0.0141521711721$.

Replacing the 70 block-coherent targets by their exact theorem coefficients gives

$$
\begin{aligned}
\text{Perron upper}&=0.293680914905,\\
\text{promise loss}&=0.0157240921033,\\
\text{total}&=0.309405007008,
\end{aligned}
$$

at

$$
\beta=0.746486947288.
$$

The new margin is

$$
\frac13-0.309405007008
=0.0239283263253.
$$

Thus this theorem recovers $0.00977615515317$ of additional routing reserve.

## Remaining contraction boundary

The 818 balanced entries with an internal block split remain open. The next largest reusable class is the 486 entries with at least one internally split cubic block. A naive product of independent disjointness-mask norms does not pass the ledger; the needed lemma must retain the shared chain or collision completion.

Subsequently, the ten inherited chain-aware theorem families close 40 of
these entries and reduce the live remainder to 778; see
`Q64_CHAIN_AWARE_INSERTION.md`. The 818-entry count above is the boundary of
the block-coherent insertion itself.

The current result is one-batch only. Outward-rounded Perron certification and the adaptive posterior-selection lift remain mandatory before passive hard dose six is excluded.

## Reproduction

Run:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_block_coherent_contraction.py --write-artifact
```

The committed artifact is `artifacts/q64_block_coherent_contraction.json`. The regression verifies every exact coefficient, the worst entry, the ledger insertion, and byte-for-byte artifact regeneration.
