# The $q=64$ masked local-Walsh repair

Date: 2026-07-16

Status: arbitrary-correlated-diagonal one-batch theorem for 180 actual
distinctness-masked occurrence entries. It extends the 54-entry quintic
repair and leaves 120 of the original 354 coefficient-one-dependent entries
quarantined.

## Local $q^{-1}$ mechanisms

There are two cases.

First, suppose an internal singleton $z$ lies between odd-degree supports $S$
and $T$; either neighbor may itself be a singleton. The exact
signed-permutation endpoint formula has the
form

$$
M_{d1}(S,z)=v_d(S)H_N(\operatorname{xor}S,z),
\qquad |v_d(S)|\le1.
$$

Multiplying the two links through $z$ gives

$$
M_{d1}(S,z)M_{1e}(z,T)
={v_d(S)v_e(T)\over q}
H_N(\operatorname{xor}S,\operatorname{xor}T).
$$

The residual normalized Walsh matrix has coefficient one after duplicate
compression, so the chain retains an explicit factor $q^{-1}$.

Second, if two adjacent singleton blocks lie wholly on the same occurrence
side, their normalized Walsh link is a row-only or column-only scalar of
magnitude $q^{-1}$. Removing it again leaves completed unmasked link kernels
of coefficient one.

These are statements about the completed kernels only. The physical masks
are restored next.

After extracting the displayed scalar, all endpoint amplitudes and all other
links are kept together as one completed moment kernel. The valid part of the
inherited lemma represents that joint unmasked kernel as a cross Gram of unit
hidden-assignment features. Thus it has coefficient one for arbitrary
correlated diagonal laws. The proof does not multiply independently optimized
link bounds or assume a product occurrence law.

## Every mask is paid explicitly

For an $r$-set and $s$-set, inclusion--exclusion writes the disjointness mask
as a sum of shared-subset Gram kernels and gives

$$
\gamma_{r,s}le
\sum_{t=0}^{\min(r,s)}
\sqrt{{r\choose t}{s\choose t}}.
$$

The theorem uses the exact integer upper bound obtained by replacing each
square root by its ceiling. If $G_s$ is the product of these integer factors
over every split occurrence block, then

$$
c_s\le {G_s\over64}.
$$

Across the 180 entries, $G_s$ ranges from 3 to 54. Hence every coefficient is
at most

$$
{54\over64}=0.84375<1.
$$

No physical mask is dropped, no product law is assumed, and the invalid claim
that masking preserves a coefficient-one Gram contraction is not used.

## Inventory and verification

The original decomposition contains 80 higher--singleton--higher entries and
52 same-side singleton-pair entries. The general internal-singleton identity
also closes 48 opposite-side singleton-pair chains. The theorem is closed under
occurrence complement and path reversal. Together with the previous 54-entry
theorem, 234 of the 354 affected entries are now independently proved and 120
remain quarantined.

At $q=4$, the regression verifies the endpoint-character identity exactly for
every support of degrees 3, 5, 7, and 9, every singleton, and both link
orientations: 889856 exact comparisons. This checks the finite degree list
used by the theorem; the analytic argument supplies arbitrary $q$.

Reproduce with:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_masked_local_walsh_repair.py --output artifacts/q64_masked_local_walsh_repair.json
/opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_masked_local_walsh_repair.py
```
