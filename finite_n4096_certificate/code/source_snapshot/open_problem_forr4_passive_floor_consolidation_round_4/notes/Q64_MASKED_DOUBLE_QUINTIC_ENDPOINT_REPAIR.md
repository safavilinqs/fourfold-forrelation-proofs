# The $q=64$ masked double-quintic endpoint repair

Date: 2026-07-16

Status: arbitrary-correlated-diagonal one-batch theorem for six actual masked
$(1,5,5,1)$ entries. Together with the previous repairs, 252 of 354 affected
entries are proved and 102 remain quarantined.

## Two physical endpoint factors

Let $E_k$ be the exact physical degree-five endpoint squared slice with $k$
quintic cells on the row occurrence side, and let $s\in\{0,1\}$ say whether
the adjacent singleton is on that side. Factoring through complete physical
rows or columns gives one endpoint squared coefficient

$$
Q_{k,s}=\min\{q^{2(1-s)}E_k,q^{2s}E_{5-k}\}.
$$

The occurrence mask inside that quintic block is already present in the
complete row or column. Composing the two endpoint factorizations while keeping
the middle completed quintic--quintic moment kernel as one unit-feature cross
Gram gives

$$
c_s^2\le Q_{s_2,s_1}Q_{s_3,s_4}.
$$

This composition is a Schur-multiplier factorization and therefore preserves
arbitrary correlations between every row and column occurrence variable.

## Result

The bound is at most one for exactly six of the 18 residual
$(1,5,5,1)$ entries:

- four with one whole quintic endpoint; and
- two with compatibly oriented $1|4$ quintic endpoints.

The maximum exact squared coefficient is

$$
\frac{1023}{1024},
$$

so the maximum outward coefficient is $0.999511599483$. The other 12 entries
are left quarantined; this theorem makes no claim about them.

Reproduce with:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_masked_double_quintic_endpoint_repair.py --output artifacts/q64_masked_double_quintic_endpoint_repair.json
/opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_masked_double_quintic_endpoint_repair.py
```
