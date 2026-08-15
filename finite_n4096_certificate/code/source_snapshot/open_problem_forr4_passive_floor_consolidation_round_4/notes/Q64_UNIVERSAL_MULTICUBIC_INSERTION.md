# The $q=64$ universal multicubic insertion

Date: 2026-07-16

Status: quarantined historical argument. The inherited cross-Gram proof omits
the cross-cut distinctness masks and does not establish coefficient one for
these 14 physical entries. See `Q64_MASKED_UNIVERSAL_AUDIT.md`.

## Result

The universal cross-Gram lemma proved in `Q64_UNIVERSAL_SEPTIMIC_INSERTION.md` gives coefficient at most one for every fixed occurrence split. Applying it to

- the 8 entries with exactly three split cubic blocks; and
- the 6 entries with four split cubic blocks

closes both classes without another signed-permutation-specific contraction.

The three-split entries do not change the current Perron optimum when raised from their provisional target to one. Raising the four-split entries to one gives

$$
\begin{aligned}
\text{Perron upper}&=0.314565277938,\\
\text{promise loss}&=0.0173705514960,\\
\text{total}&=0.331935829434,
\end{aligned}
$$

at $\beta=0.746071749555$. The margin is

$$
\frac13-0.331935829434
=0.00139750389895.
$$

This remains above the declared $10^{-3}$ reserve, with $0.000397503898949$ of additional routing room.

## Revised boundary

There are now 220 theorem entries and 668 open entries. The universal coefficient-one screen fails for every other remaining structural class, even when each class is tested alone. Thus subsequent progress must improve a coefficient below one; entry counting alone can no longer close a class.

The leading specialized target remains the 184 entries with one split cubic and one split quintic. After both universal insertions, their live common reserve gate is $0.125261095651$. The two inherited local fixed-slice scales are $0.123974636390$ for a $1|4$ quintic split and $0.149556115743$ for a $2|3$ split, so the latter cannot simply be inserted across all placements. Even after assigning the first 104 entries their local scale, the 80 balanced-split entries have reserve gate $0.125681339751$.

## Reproduction

Run:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_universal_multicubic_insertion.py --write-artifact
```

The committed artifact is `artifacts/q64_universal_multicubic_insertion.json`.

The live quintic gates and the failed naive local-slice completion are
recomputed in `Q64_POST_UNIVERSAL_QUINTIC_GATE.md`.
