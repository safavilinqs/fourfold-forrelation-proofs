# The live $q=64$ quintic gate

Date: 2026-07-16

Status: corrected floating routing gate after all 220 current theorem insertions. The local slice values below are proof diagnostics, not all-placement arbitrary-law theorems.

## Live class

The remaining split-cubic/split-quintic class has 184 entries:

| quintic split, up to complement | entries | inherited local slice scale |
|---|---:|---:|
| $1|4$ | 104 | $0.123974636390$ |
| $2|3$ | 80 | $0.149556115743$ |

With all 184 entries assigned one common coefficient, the live coefficient retaining the declared $10^{-3}$ reserve is

$$
c_{\mathrm{common,res}}=0.125261095651.
$$

Thus the $1|4$ local scale fits numerically, while the $2|3$ scale does not.

## Joint local-slice proxy fails

Assign the 104 extreme splits coefficient $0.123974636390$ and reoptimize the remaining 80 balanced splits. Their reserve gate is only

$$
c_{2|3,\mathrm{res}}=0.125681339751.
$$

The inherited fixed-slice product $0.149556115743$ is 19.0 percent above that gate. Assigning both local scales gives

$$
\begin{aligned}
\text{Perron upper}&=0.320532733008,\\
\text{promise loss}&=0.0177153486568,\\
\text{total}&=0.338248081665,
\end{aligned}
$$

which exceeds $1/3$ by

$$
0.00491474833148.
$$

This falsifies the idea that the two already visible local slice factors can simply be declared uniformly across the class after the coarse universal closures.

## Next proof decision

The next useful theorem must do at least one of the following:

1. prove a shared chain-aware coefficient at most $0.125681339751$ for the 80 balanced quintic splits after closing the 104 extreme splits;
2. recover at least $0.00491475$ in the optimized total by replacing some coefficient-one septimic or multicubic insertions with sharper bounds; or
3. prove a nonuniform joint coefficient vector whose reoptimized total retains the reserve.

The second option is attractive because the coefficient-one insertions are deliberately coarse. The first option requires about a 16 percent coefficient improvement over the raw balanced local-slice product. An isolated-orbit result is useful only if its recovered Perron margin is reported against this deficit.

## Reproduction

Run:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_post_universal_quintic_gate.py --write-artifact
```

The committed artifact is `artifacts/q64_post_universal_quintic_gate.json`.

## Subsequent repair

The shifted middle-pair theorem subsequently closes the highest-impact
four-entry orbit and recovers $0.00689076608741$ of margin, more than the
proxy deficit above. See `Q64_SHIFTED_MIDDLE_PAIR_CONTRACTION.md`.
