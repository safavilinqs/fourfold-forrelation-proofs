# The $q=64$ shifted middle-pair contraction

Date: 2026-07-16

Status: rigorous arbitrary-diagonal one-batch coefficient for the highest-impact unresolved quintic orbit. It closes four more entries, bringing the high-sector theorem count to 224 of 888. It does not prove the other 664 entries or the adaptive lift.

## Target and theorem

For profile $(1,3,5,1)$ and split $(0,1,3,1)$, write

$$
C=\{x\}\mathbin{\dot\cup}E,
\qquad
S=F\mathbin{\dot\cup}G,
$$

with $|E|=2$, $|F|=3$, and $|G|=2$. Rows are $(x,F,d)$, columns are $(a,E,G)$, and the exact occurrence kernel is

$$
K_{(x,F,d),(a,E,G)}
=M_{1,3}(a,C)M_{3,5}(C,S)M_{5,1}(S,d).
$$

This is the shifted $1|2$ cubic, $3|2$ quintic version of the accepted Round 3 middle-pair row theorem.

Let $E_1(q)$ be the exact cubic endpoint squared slice through one fixed cell and $F_3(q)$ the exact quintic endpoint squared slice through a fixed triple. At $q=64$,

$$
E_1=0.500244140625,
\qquad
F_3=1.45384579613.
$$

The compatible middle link has record-one and record-three maxima

$$
m_1=0.000264016897081,
\qquad
m_3=1.18040341467\times10^{-6},
$$

so $m=m_1$.

Fixing a complete row and taking the middle-link maximum separates the endpoint sums:

$$
\sum_{a,E,G}|K|^2
\le N E_1(q)F_3(q)m(q)^2
=0.000207646085656.
$$

Use the normalized complete kernel row as a Schur feature and the standard basis on the column side, exactly as in the accepted middle-pair theorem. The weighted all-ones base has nuclear norm equal to the geometric mean of the two diagonal masses. Hence every arbitrary correlated diagonal law obeys

$$
\left\|D_p^{1/2}KD_r^{1/2}\right\|_1
\le0.0144099301059\sqrt{\left(\sum p\right)\left(\sum r\right)}.
$$

Complement and path reversal close the four-entry orbit.

## Finite-size effect

Replacing the live target $0.124035215254$ on these four entries gives

$$
\begin{aligned}
\text{Perron upper}&=0.308066639351,\\
\text{promise loss}&=0.0169784239964,\\
\text{total}&=0.325045063347,
\end{aligned}
$$

at $\beta=0.746167078477$. The margin is

$$
0.00828826998635,
$$

an improvement of $0.00689076608741$ over the preceding routing ledger.

This single theorem recovers more than the $0.00491474833148$ deficit of the naive remaining-quintic local-slice proxy. Keeping this proved orbit at its actual coefficient and assigning the other 180 quintic entries their two local slice scales gives diagnostic total $0.329863829155$, with margin $0.00346950417846$. The remaining task is to prove those all-placement coefficients; the numerical obstruction has been removed.

## Reproduction

Run:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_shifted_middle_pair_contraction.py --write-artifact
```

The committed artifact is `artifacts/q64_shifted_middle_pair_contraction.json`.
