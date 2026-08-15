# The $q=64$ masked cubic-endpoint repair

Date: 2026-07-16

Status: arbitrary-correlated-diagonal one-batch theorem for 12 actual masked
occurrence entries. Together with the earlier repairs, 246 of 354 affected
entries are proved and 108 remain quarantined.

## Physical endpoint factor

For a split cubic adjacent to a singleton, the better complete-row or
complete-column squared energy is exactly

$$
E_2(q)={q^2-2q+2\over q^2(q-1)}.
$$

The row-vector factorization of the complete physical endpoint kernel gives
Schur coefficient at most $\sqrt{E_2(q)}$ for arbitrary diagonal laws. This
factor already includes the cubic occurrence mask.

The regression enumerates every fixed cubic pair and singleton at $q=4$.
Across 1920 physical rows the only squared energies are $1/24$ and $5/24$,
and the maximum $5/24$ agrees exactly with $E_2(4)$.

## Remaining masks

There are two structural cases.

- Four cubic--septimic entries have a residual $2|5$ septimic mask. Its
  inclusion--exclusion factor is $1+2\sqrt{10}<15/2$.
- Eight recovered cubic/cubic/quintic entries have one residual cubic
  $1|2$ mask and one quintic $1|4$ mask. Their product is below
  $(5/2)\cdot3=15/2$.

After extracting the physical cubic endpoint, the remaining completed moment
kernel is kept as one unit-feature cross Gram. Restoring the listed masks gives

$$
c^2\le {225\over4}E_2(64)<1,
\qquad c=0.930264\ldots.
$$

This proof retains arbitrary correlations between all occurrence blocks and
does not reuse the invalid masked coefficient-one inference.

Reproduce with:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_masked_cubic_endpoint_repair.py --output artifacts/q64_masked_cubic_endpoint_repair.json
/opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_masked_cubic_endpoint_repair.py
```
