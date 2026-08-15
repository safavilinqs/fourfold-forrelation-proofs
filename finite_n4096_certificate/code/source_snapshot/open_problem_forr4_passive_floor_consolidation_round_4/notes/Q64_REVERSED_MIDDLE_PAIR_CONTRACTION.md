# The $q=64$ reversed middle-pair contraction

Date: 2026-07-16

Status: rigorous arbitrary-diagonal one-batch coefficient for two high-impact unresolved quintic orbits. It closes eight entries, bringing the high-sector theorem count to 232 of 888. It does not prove the other 656 entries or the adaptive lift.

## Target

For profile $(1,1,5,3)$ and splits $(0,1,3,1)$ and $(0,1,2,2)$, write the kernel as

$$
K=H(a,b)M_{1,5}(b,S)M_{5,3}(S,C),
$$

with $S=F\mathbin{\dot\cup}G$, $|F|=3$, $|G|=2$, and $C=E\mathbin{\dot\cup}\{x\}$, $|E|=2$. The declared rows are $(b,F,x)$ and columns are $(a,G,E)$.

For the first split, use the complete column as a row of the transposed occurrence matrix. For the second split, use the declared complete row. In both cases fix a quintic pair and a cubic pair, then sum the complementary quintic triple and cubic singleton. The normalized Hadamard factor has squared modulus $1/N$. For each fixed adjacent singleton, the exact quintic fixed-pair slice bounds the quintic completion by

$$
F_2(64)=41.3799758185.
$$

Summing over $b$ contributes $N$, exactly canceled by $|H(a,b)|^2=1/N$. There are exactly $N-2=4094$ choices of $x$ disjoint from the fixed cubic pair $E$. The compatible middle-link maximum is

$$
m(64)=0.000264016897081.
$$

Therefore the complete transposed-row energy is at most

$$
(N-2)F_2(64)m(64)^2
=0.0118086844085.
$$

The same complete-row Schur-feature argument as the preceding theorem gives

$$
\left\|D_p^{1/2}KD_r^{1/2}\right\|_1
\le0.108667770790\sqrt{\left(\sum p\right)\left(\sum r\right)}.
$$

Complement and path reversal close both four-entry orbits.

## Finite-size effect

The coefficient is below the live target $0.124035215254$. Inserting it gives

$$
\begin{aligned}
\text{Perron upper}&=0.306868882423,\\
\text{promise loss}&=0.0169078984975,\\
\text{total}&=0.323776780921,
\end{aligned}
$$

with margin $0.00955655241272$. The gain over the previous routing ledger is $0.00126828242636$.

Keeping all three new middle-pair orbits at their theorem coefficients and assigning the other 172 quintic entries their local slice scales gives diagnostic total $0.326448847879$, with margin $0.00688448545451$. The remaining issue is theorem coverage, not the scalar budget.

## Reproduction

Run:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_reversed_middle_pair_contraction.py --write-artifact
```

The committed artifact is `artifacts/q64_reversed_middle_pair_contraction.json`.
