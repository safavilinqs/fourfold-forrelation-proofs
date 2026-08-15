# The $q=64$ remaining-class gate

Date: 2026-07-16

Status: exact generated inventory and floating routing gate after the 110 proved high-sector insertions. The coefficients in this note are proof targets, not arbitrary-law theorems.

## Revised inventory

The 778 remaining balanced entries partition as follows:

| class | entries | current routing target |
|---|---:|---:|
| one split cubic and one split higher block | 280 | $0.124035215254$ |
| only a higher block split, in a cubic profile | 176 | $0.124035215254$ |
| no cubic block | 140 | $0.5$ |
| two split cubics and one split higher block | 96 | $0.124035215254$ |
| one split cubic and no split higher block | 48 | $0.124035215254$ |
| two split cubics and no split higher block | 24 | $0.124035215254$ |
| three split cubics | 8 | $0.124035215254$ |
| four split cubics | 6 | $0.124035215254$ |

The rows are disjoint and exhaust the remaining set. The 140-entry noncubic row contains entries with either one or two internally split higher blocks.

## Lead shared-contraction gate

Leave every other unresolved entry at its existing two-tier routing target and assign one common coefficient $c$ to the 280-entry lead class. Reoptimizing the attenuation gives

$$
c_{\mathrm{threshold}}=0.225536743566
$$

at total $1/3$, and

$$
c_{\mathrm{reserve}}=0.222921146951
$$

at total $1/3-10^{-3}$.

The reserve gate is

$$
\frac{0.222921146951}{0.124035215254}
=1.79724077953
$$

times the original cubic routing target. Thus the next theorem does not need to preserve the near-$1/8$ coefficient that motivated the first $q=64$ screen. It may lose almost a factor of $1.8$ and still preserve the project's declared one-batch reserve.

## Decision

The next mathematical target is:

> Prove or falsify one arbitrary-diagonal chain-aware contraction with common coefficient at most $0.222921146951$ for all 280 entries having exactly one internally split cubic and one internally split higher block.

This is preferable to resolving another isolated orbit. The class contains 184 quintic entries and 96 septimic entries. The quintic fixed-slice factors already visible in the accepted local theorems are below this gate, but extending them across every placement and record sector is not yet proved. The septimic slices have not yet been calibrated.

If the shared contraction passes, insert its actual coefficient vector and recompute the next class gate. If a valid arbitrary-law lower obstruction exceeds the gate, stop this contraction and use the remaining margin on another class or pivot the witness. Interval certification and the adaptive lift remain downstream of a complete one-batch coefficient vector.

## Reproduction

Run:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_remaining_class_gates.py --write-artifact
```

The committed artifact is `artifacts/q64_remaining_class_gates.json`. The regression reconstructs the proved/open partition, reoptimizes both lead gates, and compares the artifact byte for byte.

## Subsequent insertion

The universal coefficient-one Gram bound subsequently closes the 96
septimic entries in the lead row. The live specialized target is therefore
the remaining 184 quintic entries; see
`Q64_UNIVERSAL_SEPTIMIC_INSERTION.md`. The 778-entry partition and
$0.222921146951$ joint gate above remain the correct pre-insertion audit.
