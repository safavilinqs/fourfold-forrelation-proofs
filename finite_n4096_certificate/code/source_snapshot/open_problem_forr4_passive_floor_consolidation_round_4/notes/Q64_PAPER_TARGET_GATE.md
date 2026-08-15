# The $q=64$ paper-target gate

Date: 2026-07-16

## Result

$N=4096$ is now the most promising fallback size for closing the paper result. Evaluating every inherited lower-sector coefficient at its actual order $q=64$ gives a much more forgiving one-batch target than the $N=1024$ ledger:

| high-sector model | optimized one-batch total | margin to $1/3$ | status |
|---|---:|---:|---|
| common $1/q=1/64$ | $0.113490308350$ | $0.219843024983$ | routing diagnostic |
| common $1/\sqrt q=1/8$ | $0.241894419850$ | $0.091438913483$ | routing diagnostic |
| largest common coefficient at threshold | coefficient $0.199910665542$ | zero | gate |
| largest common coefficient with $10^{-3}$ reserve | coefficient $0.199089176072$ | $0.001$ | project gate |

The high-sector coefficients in this table are proposed envelopes, not proved bounds. The lower degree and known-high-sector formulas are inherited arbitrary-law theorems evaluated at $q=64$.

## A sufficient two-tier target

The 888 balanced entries in the open high-sector profiles divide into:

- 724 entries whose profile contains a cubic block; and
- 164 entries with no cubic block.

For the contraction proof, the finer structural worklist is:

| class | entries | intended mechanism |
|---|---:|---|
| whole-block coherent | 70 | inherited weighted three-link cut table |
| internally split cubic | 486 | one shared cubic-slice/dressing lemma |
| cubic profile, only a higher block split internally | 192 | cubic Gram dressing around a high-support slice |
| no cubic, higher block split internally | 140 | coarse coefficient at most $1/2$ |

The four rows are disjoint and exhaust all 888 balanced entries. The 486-entry internally split cubic class is the first theorem target; proving it one entry at a time would discard the point of the reduction.

The first row is now complete. Exact record-sector application of the weighted three-link cut table closes all 70 whole-block-coherent entries, improves the routing total to $0.309405007008$, and leaves margin $0.0239283263253$. See `Q64_BLOCK_COHERENT_CONTRACTION.md`.

Ten inherited chain-aware theorem families close 40 additional entries,
disjoint from the first row. The combined insertion certifies 110 entries,
improves the routing total to $0.296090867182$, and leaves margin
$0.0372424661512$. See `Q64_CHAIN_AWARE_INSERTION.md`.

Assign the cubic-containing profiles the exact $q=64$ cubic fixed-pair slice scale

$$
c_{\mathrm{cub}}
=
\sqrt{\frac{q^2-2q+2}{q^2(q-1)}}
=0.124035215254,
$$

and assign the remaining profiles the deliberately loose target

$$
c_{\mathrm{noncub}}=\frac12.
$$

The reoptimized occupation ledger is then

$$
0.302847337369
+0.0163338247925
=0.319181162161,
$$

where the second term is the promise-conditioning loss. The optimizing attenuation is

$$
\beta=0.746328499244,
$$

and the margin is

$$
\frac13-0.319181162161
=0.0141521711721.
$$

This is enough room for outward rounding and a modest adaptive-lift loss. It is not permission to spend the margin before those terms are proved.

## Why this changes the plan

At $N=1024$, the unresolved families had to average around $1/q$ and the diagnostic margin was only $2.0\times10^{-4}$. At $N=4096$, even the common $q^{-1/2}$ envelope passes with margin $0.0914$. The two-tier envelope allows the no-cubic minority to be four times larger than $q^{-1/2}$ and still passes.

The immediate mathematical task is therefore no longer “prove 224 near-$1/q$ orbit bounds.” It is:

1. prove one shared contraction at or below $0.124035215254$ for every open profile containing a cubic block;
2. prove the coarse coefficient $1/2$ for the 164 remaining balanced entries;
3. insert the actual coefficient vector and intervalize the Perron calculation; and
4. prove the adaptive posterior-selection lift.

If either coefficient target fails, use the common gate $0.19909$ to trade strength between the two families and reoptimize globally rather than reverting to orbit-by-orbit closure.

## Scope and experimental status

A successful one-batch proof would use

$$
N=4096,
\qquad
M=4N=16{,}384\text{ sign modes}.
$$

No current result establishes that 16,384 sign modes are experimentally credible. This size is inside the predeclared Round 4 feasibility window but must pass the experimental scorecard separately. The active protocol still uses hard dose six; its concrete operations and ancillary modes remain to be itemized.

The displayed calculation also remains one-batch only. It does not exclude unrestricted adaptive passive dose six until the adaptive lift is proved.

## Reproduction

Run:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_paper_target_gate.py --write-artifact
```

The committed output is `artifacts/q64_paper_target_gate.json`. The regression regenerates the calculation, verifies that the inherited module orders are restored afterward, and compares the artifact byte for byte.
