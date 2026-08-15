# The $q=64$ fixed-singleton pair contraction

Date: 2026-07-16

Status: rigorous arbitrary-diagonal one-batch coefficient for another high-impact quintic orbit. It closes four entries, bringing the high-sector theorem count to 236 of 888. It does not prove the other 652 entries or the adaptive lift.

## Theorem

For profile $(1,1,5,3)$ and split $(0,0,3,2)$, pass to the complementary split $(1,1,2,1)$. A complete row then fixes both singleton blocks, a quintic pair, and one cubic cell. The remaining column variables are a quintic triple and a cubic pair.

The normalized Hadamard link contributes squared factor $1/N$. The exact fixed-pair quintic endpoint energy is

$$
F_2(64)=41.3799758185.
$$

There are

$$
\binom{N-1}{2}=8{,}382{,}465
$$

cubic pair completions around the fixed cell. Taking the universal middle-link maximum

$$
m(64)=0.000264016897081
$$

gives complete-row energy

$$
\frac1N F_2(64)\binom{N-1}{2}m(64)^2
=0.00590290071446.
$$

The complete-row Schur-feature argument therefore proves

$$
\left\|D_p^{1/2}KD_r^{1/2}\right\|_1
\le0.0768303372012\sqrt{\left(\sum p\right)\left(\sum r\right)}.
$$

Complement and reversal close the four-entry orbit.

## Finite-size effect

Inserting the coefficient gives total $0.323034695004$ and margin $0.0102986383292$, improving the preceding routing margin by $0.000742085916459$. There are now 236 theorem entries and 652 open entries.

Keeping the four newly proved quintic orbits at their actual coefficients and assigning the other 168 quintic entries their two local slice scales gives diagnostic total $0.325310000395$, with margin $0.00802333293854$. The scalar budget is increasingly comfortable; coverage remains the live issue.

## Reproduction

Run:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_fixed_singleton_pair_contraction.py --write-artifact
```

The committed artifact is `artifacts/q64_fixed_singleton_pair_contraction.json`.
