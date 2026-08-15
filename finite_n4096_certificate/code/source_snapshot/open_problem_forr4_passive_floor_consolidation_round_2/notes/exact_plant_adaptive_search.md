# Exact-plant adaptive falsification search

Date: 2026-07-14

Scope: $q=2$, $N=4$ exact signed-permutation plant. The root is a pure passive probe with a binary rank-one effect. Given each root outcome, the child one-batch protocol is optimized by the exact Schur-multiplier SDP. This is a strict protocol slice, not an upper bound over all adaptive protocols.

## Exact hard pair

There are eight signed permutations at $q=2$. Enumerating all $8^3=512$ latent triples constructs the positive law with $F_{4,H}=+1$ and the negative law with $F_{4,H}=-1$, including duplicate physical inputs with their exact multiplicities.

## Results

The one-batch dose-two SDP reproduces

$$
0.50000000224\le\Delta_2\le0.50000001156,
$$

consistent with the exact crossing-pair value $1/q=1/2$.

For a $1+1$ two-stage split, five starts and eight alternating root updates give

$$
0.24999999968\le\operatorname{TV}
\le0.25000002046
$$

for the final fixed root, where the upper number is the repaired child dual bound conditional on that root. The optimized child diagonal weights are uniform over eight one-photon modes in the first two physical blocks.

Thus the search finds no adaptive amplification in this slice; it converges to half the best one-batch dose-two value.

An exploratory $1+2$ run uses a $137\times137$ child SDP and did not complete within the useful diagnostic window. It was terminated without a result and supplies no evidence either way.

## Interpretation

The $1+1$ result is consistent with the idea that temporal splitting loses the crossing-pair coherence that attains $1/\sqrt N$ in one batch. It does not prove this uniformly, because:

- the root receiver is binary and rank one;
- only local alternating updates are used for the root;
- $q=2$ has strong finite-size degeneracies; and
- the fixed-root dual does not optimize over all roots.

Reproduction: searches/exact_plant_adaptive_search.py with `--starts 5 --rounds 8`.
