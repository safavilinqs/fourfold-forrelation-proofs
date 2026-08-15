# q64 same-side whole-link insertion

Date: 2026-07-16

Status: the 96-entry mask-aware theorem is preserved, but its upstream
712-entry baseline is invalid. The cumulative 808-entry count, routing total,
and adaptive allowance are withdrawn. See `Q64_MASKED_UNIVERSAL_AUDIT.md`.

## Result

Every entry in this theorem has exactly one split higher block and an adjacent pair of whole blocks on one physical side. The whole pair is always singleton--singleton, singleton--cubic, or cubic--singleton. Its scalar link has magnitude at most $1/q$. Completing the other physical link kernels and restoring the only split-block mask gives coefficient

$$
c_{d,k}(q)
={1\over q}\sum_{j=0}^{\min(k,d-k)}
\sqrt{\binom{k}{j}\binom{d-k}{j}},
$$

where the split higher block has degree $d$ and smaller side $k$.

At $q=64$ the five templates are

| split | entries | mask factor | coefficient |
|---|---:|---:|---:|
| $1|4$ | 8 | $3$ | $0.046875$ |
| $2|3$ | 28 | $1+\sqrt6+\sqrt3$ | $0.0809615710993$ |
| $1|6$ | 24 | $1+\sqrt6$ | $0.0538982772307$ |
| $2|5$ | 24 | $1+2\sqrt{10}$ | $0.114446176880$ |
| $3|4$ | 12 | $3+\sqrt{12}+\sqrt{18}$ | $0.167292848473$ |

The last coefficient is above the old common target, but the full coefficient vector improves the optimized routing total from

$$
0.328477421166173
$$

to

$$
U_{\mathrm{route}}=0.323362582871308.
$$

The raw diagnostic margin is

$$
{1\over3}-U_{\mathrm{route}}
=0.009970750462025.
$$

Retaining the declared $10^{-3}$ numerical allowance conditionally leaves $0.008970750462025$ for the adaptive lift. Eighty entries remain at their frozen target; their regenerated common reserve gate is

$$
0.190775718804.
$$

## Same-side whole-link scalar

After fixing one matrix side, a whole singleton--singleton link is one normalized Walsh entry and has magnitude $1/q$. A whole singleton--cubic or cubic--singleton link has the endpoint form

$$
M_{13}(a,T)=H_N(a,\xi(T))v_3(T)
$$

or its transpose. Since $|v_3(T)|\le1$ and every normalized Walsh entry has modulus $1/q$, its magnitude is also at most $1/q$.

This factor is row-only or column-only because both blocks lie wholly on the same side of the occurrence cut. It is therefore a diagonal scalar multiplier. Complete every other physical moment kernel as a cross Gram of unit character features; their trace-class Schur-multiplier norms are at most one.

## One split-block mask

Let $A$ and $B$ be the two portions of the unique split higher support, with $|A|=k$ and $|B|=d-k$. Inclusion--exclusion gives

$$
\mathbf 1_{A\cap B=\varnothing}
=\sum_{j=0}^{k}(-1)^j
  \sum_{J:\,|J|=j}
  \mathbf 1_{J\subseteq A}\mathbf 1_{J\subseteq B}.
$$

For level $j$, the row feature has squared multiplicity $\binom{k}{j}$ and the column feature has squared multiplicity $\binom{d-k}{j}$. Optimally scaling the orthogonal direct-sum levels gives Schur factor

$$
\Delta_{k,d-k}
=\sum_j\sqrt{\binom{k}{j}\binom{d-k}{j}}.
$$

There are no other internally split blocks, so there are no other cross-cut distinctness masks. Composing this mask with the $1/q$ same-side scalar and the completed unit link kernels proves the displayed coefficient for an arbitrary correlated diagonal law.

## Regression scope

The regression checks the complete 96/80 inventory split, all 24 symmetry orbits, all five template multiplicities, the inclusion direct-sum identity for every relevant subset size, sampled exact $q=4$ and $q=8$ singleton--cubic endpoint moments, the optimized q64 insertion, the regenerated 80-entry gate, and byte-for-byte artifact reproduction.

Reproduce with:

    /opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_same_side_whole_link_insertion.py --write-artifact
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_same_side_whole_link_insertion.py
