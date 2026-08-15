# Shared-frontier row-orbit screen

Date: 2026-07-15

## Question

Can a native-$q=32$ translation-orbit law supply a high-rank compatible lower witness across the leading 51 unresolved families?

The simplest such law is now tested and closed. This is a negative result about one physical-law family, not about arbitrary passive laws or the shared upper contraction.

## Physical law

Use the 30 occupation states in the current Perron support. For each state $n$, choose one four-block base configuration $A_n$. Apply the same row-XOR translation $u\in\mathbb F_2^5$ to all four blocks and average uniformly over the 32 translated configurations.

Together with the current occupation weights $\rho_n$, this gives the legal diagonal probe law

$$
p(n,u)=\frac{\rho_n}{32}.
$$

Its support contains $30\times32=960$ configurations directly at $q=32$ and $N=1024$. Unlike the earlier embedded-$q=4$ sparse law, it uses the native coordinate geometry.

## Exact screen

Nine deterministic base-configuration candidates were scored using the exact signed-permutation moments. Selection maximized, in order:

1. frontier orbits with a nonzero exact moment;
2. frontier entries with a nonzero exact moment;
3. the sum of absolute exact moments; and
4. the number of nonzero configuration pairs.

The selected trial has 50,112 configuration pairs whose symmetric-difference pattern belongs to the frontier. The signed-permutation parity rule annihilates nearly all of them: only 1,984 pairs, six profile-split entries, and six of the 51 orbits remain nonzero.

The attenuated trace norms are:

| quantity | value |
|---|---:|
| separate sum over the 51 orbit matrices | $4.79416996109\times10^{-13}$ |
| joint frontier matrix | $4.51866434463\times10^{-13}$ |
| joint / separate | $0.942533197884$ |
| current 51-orbit Perron contribution | $0.0282264705537$ |
| separate / current frontier | $1.69846596725\times10^{-11}$ |

All plant moments are exact rationals. The final eigendecompositions use ordinary floating point.

## Decision

Close the common-row-translation family as a lower-witness route. Its rank is higher than the first sparse law and it lives directly at $q=32$, but the exact parity constraint leaves it eleven orders of magnitude below the current frontier scale. Further seed tuning cannot repair that structural mismatch.

This does not show that the signed-permutation witness survives the joint lower test. A useful high-rank law must deliberately populate the odd-domain/odd-codomain parity matches of each signed-permutation link. The remaining candidates are:

- a compressed coordinate-uniform law;
- a mixture of multiple translation orbits chosen by parity type; or
- direct movement to the shared arbitrary-law upper contraction.

Do not infer an upper bound from the small physical witness.

## Reproduction

Run:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/shared_frontier_row_orbit.py --write-artifact
```

The committed output is `artifacts/shared_frontier_row_orbit.json`; the regression regenerates it byte for byte.
