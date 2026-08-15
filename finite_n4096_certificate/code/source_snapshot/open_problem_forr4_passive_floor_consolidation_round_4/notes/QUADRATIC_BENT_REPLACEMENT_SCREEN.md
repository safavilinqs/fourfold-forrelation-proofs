# Quadratic-bent replacement screen

Date: 2026-07-15

## Decision

Do not promote the quadratic-bent exact plant to the lead Round 4 route yet.
It removes the signed-permutation promise-conditioning loss, but the current
evidence does not provide a complete one-batch coefficient architecture at
(N=1024). Its finite-size requirement is now explicit.

## What was screened

The quadratic-bent plant has pointwise

$$
F_{4,H}=+1
$$

on the positive hypothesis and (-1) after first-block reflection. Its
natural scalar-ledger screen therefore uses

$$
\beta=1,
\qquad
\text{promise loss}=0.
$$

Retain the current dose-six occupation compatibility exactly. Replace every
nonminimal degree-six, eight, ten, and twelve profile-split coefficient by a
hypothetical common value (c). This does not claim that the physical plant
has common coefficients. It asks how strong a complete theorem would have
to be before the exact promise becomes useful.

The inventory contains:

| degree | dose-compatible profile-split entries |
|---:|---:|
| 6 | 128 |
| 8 | 476 |
| 10 | 920 |
| 12 | 760 |
| total | 2,284 |

There are 69 nonminimal profiles and 7,904 raw profile-split entries before
occupation compatibility removes irrelevant cuts.

## Finite-size gate

Setting every higher-sector coefficient optimistically to zero gives

$$
T_{\mathrm{floor}}=0.281512032891,
$$

leaving

$$
\frac13-T_{\mathrm{floor}}=0.0518213004423.
$$

After reoptimizing the occupation Perron problem, the common coefficient
reaches the threshold at

$$
\boxed{
c_*=0.000529808595115
=\frac{0.542524001398}{N}.
}
$$

Two useful calibration points are:

| hypothetical common coefficient | scalar total | margin / overshoot |
|---:|---:|---:|
| (1/(2N)) | (0.328670235127) | (+0.004663098206) margin |
| (1/N) | (0.387618890853) | (0.0542855575202) overshoot |

Thus “all sectors are (O(1/N))” is not enough. The constants and joint
profile distribution matter. A simple sufficient target is a complete
shared contraction at (1/(2N)), or a nonuniform vector whose reoptimized
total is no larger.

For reference, if only one degree were nonzero, its common-coefficient gate
would be:

| degree | isolated gate |
|---:|---:|
| 6 | (0.00162040162516) |
| 8 | (0.00132809868817) |
| 10 | (0.00229500099696) |
| 12 | (0.00967656775424) |

These isolated numbers cannot be spent simultaneously.

## Relation to the inherited endpoint theorem

The proved weighted quadratic-bent endpoint coefficient for (M_{5,1}) is

$$
\frac{2}{N-2}=\frac1{511}=0.00195694716243.
$$

This is not a failure of the candidate: a four-block profile contains three
links, and endpoint gains may multiply or combine with shared middle-link
contraction. But the endpoint theorem alone cannot populate the scalar
ledger. The unresolved mathematical object remains the arbitrary-split
three-link Schur kernel, especially weighted two-sided internal sectors
(M_{a,b}) with (a,b\ge2).

## Retain-or-pivot consequence

The signed-permutation route remains the finite-size lead because it already
has a complete generated inventory, 40 proved arbitrary-law balanced
entries, and a near-threshold global diagnostic. The quadratic-bent route
currently has a cleaner promise but a larger unpopulated inventory and no
adaptive lift.

Promote the replacement only after it supplies one of:

1. a complete profile-specific coefficient vector whose reoptimized total
   is below (1/3) with at least (10^{-3}) decision reserve; or
2. a shared arbitrary-split contraction at the sufficient (1/(2N)) scale
   or better.

If the signed-permutation shared-contraction attempt fails, the first
quadratic-bent project is now narrowly defined: prove or falsify the
weighted two-sided internal-link contraction needed to reach this gate. Do
not begin by recomputing unweighted spectra.

## Reproduction

Run:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 \
  searches/quadratic_bent_replacement_screen.py --write-artifact
```

The committed result is
`artifacts/quadratic_bent_replacement_screen.json`, and its regression
recomputes every gate and compares the artifact byte for byte.
