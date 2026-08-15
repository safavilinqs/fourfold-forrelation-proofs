# The $q=64$ chain-aware theorem insertion

Date: 2026-07-16

Status: rigorous arbitrary-diagonal one-batch coefficients for 40 additional balanced entries, obtained by evaluating the ten accepted Round 3 chain-aware theorems at $q=64$. Together with the 70 block-coherent entries, 110 of the 888 open entries now have theorem coefficients. This does not prove the remaining 778 entries or the adaptive lift.

## Result

The Round 3 theorem formulas are dimension-parameterized. Evaluating them at
$q=64$ gives:

| theorem family | entries | coefficient |
|---|---:|---:|
| leading disjointness | 4 | $0.0468406656279$ |
| adjacent cubic slice | 4 | $0.00580989204377$ |
| separated endpoint slice | 4 | $0.123974636390$ |
| internal singleton shared law | 4 | $0.00914383903817$ |
| column cubic--quintic row | 4 | $0.0156173687419$ |
| adjacent balanced row | 4 | $0.0200795672469$ |
| whole cubic--quintic triple | 4 | $0.0188399262557$ |
| middle cubic--quintic pair | 4 | $0.0134819222279$ |
| whole cubic middle pair | 4 | $0.00856128930300$ |
| double-endpoint cubic--quintic row | 4 | $0.0158437199581$ |

The 40 entries are distinct and have no overlap with the 70 block-coherent entries. Every coefficient is an arbitrary-diagonal theorem bound, not a physical-law diagnostic. The largest is the separated-endpoint coefficient

$$
0.123974636390314,
$$

which is just below the initial cubic target $0.124035215254$.

## Finite-size effect

Starting from the block-coherent insertion,

$$
0.309405007008,
$$

inserting the ten theorem families and reoptimizing the attenuation gives

$$
\begin{aligned}
\text{Perron upper}&=0.281033192775,\\
\text{promise loss}&=0.0150576744067,\\
\text{total}&=0.296090867182,
\end{aligned}
$$

at

$$
\beta=0.746667053143.
$$

The new margin is

$$
\frac13-0.296090867182
=0.0372424661512.
$$

This recovers $0.0133141398259$ of reserve beyond the block-coherent result and $0.0230902949790$ beyond the original two-tier target.

## Revised contraction boundary

The remaining 778 entries divide as follows:

| internal structure | entries |
|---|---:|
| exactly one internally split cubic and one split higher block | 280 |
| cubic profile with only a higher block split | 176 |
| no cubic block | 140 |
| two split cubics and one split higher block | 96 |
| exactly one split cubic, no split higher block | 48 |
| two split cubics, no split higher block | 24 |
| three split cubics | 8 |
| four split cubics | 6 |

The first row is now the dominant unresolved class. The enlarged reserve means a shared theorem for this 280-entry class may be substantially looser than the original $0.1240352$ cubic target. The next calculation should determine that exact common-coefficient gate, then test one chain-aware cubic--higher-block contraction against it.

The current result remains one-batch only. Outward-rounded Perron certification and the unrestricted adaptive posterior-selection lift remain mandatory before passive hard dose six is excluded.

## Reproduction

Run:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_chain_aware_insertion.py --write-artifact
```

The committed artifact is `artifacts/q64_chain_aware_insertion.json`. The regression verifies the theorem-family provenance, the disjoint 40-entry partition, the combined ledger insertion, and byte-for-byte artifact regeneration.
