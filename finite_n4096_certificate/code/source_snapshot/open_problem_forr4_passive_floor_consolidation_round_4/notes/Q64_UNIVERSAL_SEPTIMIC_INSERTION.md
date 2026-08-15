# The $q=64$ universal septimic insertion

Date: 2026-07-16

Status: quarantined historical argument. The cross-Gram proof below omits the
cross-cut distinctness mask and does not establish coefficient one for these
96 physical entries. See `Q64_MASKED_UNIVERSAL_AUDIT.md`.

## Universal coefficient-one lemma

Let $\nu$ be any probability law on sign vectors $x$, and let the row and column indices of one fixed occurrence split be parity supports $S$ and $T$. Define the unit-modulus character features

$$
f_S(x)=\chi_S(x),
\qquad
g_T(x)=\chi_T(x).
$$

The moment kernel is the cross Gram matrix

$$
K(S,T)
=\mathbb E_{x\sim\nu}\overline{f_S(x)}g_T(x).
$$

For arbitrary nonnegative row and column weights $p_S,q_T$, define operators into $L_2(\nu)$ by

$$
Ue_S=\sqrt{p_S}f_S,
\qquad
Ve_T=\sqrt{q_T}g_T.
$$

Then the weighted kernel is $U^*V$, and Schatten Hölder gives

$$
\left\|D_p^{1/2}KD_q^{1/2}\right\|_1
\le \|U\|_{\mathrm{HS}}\|V\|_{\mathrm{HS}}
=\sqrt{\left(\sum_Sp_S\right)\left(\sum_Tq_T\right)}.
$$

Thus every fixed occurrence split has the universal arbitrary-law coefficient

$$
\gamma\le1.
$$

The ledger's high odd sectors use exactly this unconditioned planted moment normalization; the null odd moment is zero, and promise conditioning remains in the separate concentration term.

## Why the coarse bound is useful here

Among the 778 entries left after the first 110 theorem insertions, 96 have exactly one internally split cubic and one internally split septimic block. Their isolated reserve gate is above one, so no signed-permutation-specific septimic slice calculation is needed.

Replacing their provisional $0.124035215254$ targets by the proved coefficient one gives

$$
\begin{aligned}
\text{Perron upper}&=0.312161099595,\\
\text{promise loss}&=0.0172221220268,\\
\text{total}&=0.329383221622,
\end{aligned}
$$

at

$$
\beta=0.746107587291.
$$

The margin is

$$
\frac13-0.329383221622
=0.00395011171117.
$$

After reserving the declared $10^{-3}$ for certification and the adaptive interface, the routing calculation retains $0.00295011171117$.

## Revised boundary

This insertion closes 96 entries without a new septimic contraction and reduces the live count from 778 to 682. The former 280-entry lead class is now the 184 quintic entries with exactly one split cubic and one split quintic block. Those entries carry substantially more Perron impact and still require a chain-aware theorem; their local fixed-slice factors do not by themselves prove all chain placements.

All other high-sector targets remain provisional. The displayed total is therefore a valid routing insertion, not yet a complete passive theorem.

## Reproduction

Run:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_universal_septimic_insertion.py --write-artifact
```

The committed artifact is `artifacts/q64_universal_septimic_insertion.json`.

## Subsequent insertion

The same coefficient-one lemma subsequently closes the 14 entries with
three or four split cubics; see `Q64_UNIVERSAL_MULTICUBIC_INSERTION.md`.
