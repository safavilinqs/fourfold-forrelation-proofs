# Complete dependency-exact outward $q=64$ ledger

Date: 2026-07-17

Status: certified arbitrary-correlated-diagonal one-batch ledger for probes block diagonal in total signal photon number. The total passes the required $1/3-10^{-3}$ reserve gate. The later direct-sum tree-frontier theorem lifts this exact total through unrestricted classical feed-forward within that number-sector-incoherent class with multiplier one; see `Q64_ADAPTIVE_TREE_FRONTIER_THEOREM.md`.

## Result

The routing inventory contains $6{,}016$ splits on the previously open degree-ten and degree-twelve profiles. Exactly 888 are balanced, all 888 have accepted one-batch coefficients, and none remains open. The other 5,128 routing entries are excluded by the number-sector-incoherent model and are not certified coefficient bounds. This audit found that 272 of them would connect different-total-number states in the at-most-six occupation space, forming 136 undirected edges. They therefore cannot be called irrelevant for a model permitting coherence between signal-number sectors. The final balanced dependencies are the independently certified twelve-entry dual-endpoint theorem and the already certified eighty-entry final residual theorem.

The certificate deliberately does not optimize the attenuation. It fixes the rational value

$$
\beta=\frac{19}{25}=0.76
$$

and obtains

$$
\begin{aligned}
P_{\rm Perron}&\le 0.2587440963855772226792307915292556241,\\
P_{\rm promise}&\le 0.0022251284066307022549501658645609374,\\
P_{\rm one\mbox{-}batch}&\le 0.2609692247922079249341809573938165614.
\end{aligned}
$$

The acceptance threshold is

$$
\frac13-10^{-3}=0.3323333333333333333333333333333333333,
$$

so the certified reserve is at least

$$
0.0713641085411254083991523759395167720.
$$

Therefore

$$
\boxed{P_{\rm one\mbox{-}batch}<\frac13-10^{-3}}.
$$

## Dependency contract

The coefficient map is assembled in theorem dependency order. First, every coefficient-one-dependent entry is reset to its pre-universal value. The nine actual-mask repair families and the joint recovered cubic--quintic theorem are then inserted. Next, the eighty final residual coefficients are inserted from their exact rational squared coefficients. Finally, the twelve independently audited dual-endpoint entries are inserted from

$$
c_{\rm dual}^2=\frac{124116095}{5549064192}.
$$

The resulting balanced status counts are:

| theorem status | entries |
|---|---:|
| inherited nonuniversal proofs | 442 |
| masked quintic slice | 54 |
| masked local Walsh | 180 |
| masked cubic endpoint | 12 |
| masked double-quintic endpoint | 6 |
| masked double-quintic record | 12 |
| masked four-cubic incidence | 38 |
| masked cubic--septimic chain | 12 |
| masked recovered cubic--quintic endpoint row | 28 |
| masked joint recovered cubic--quintic | 12 |
| final residual chain | 80 |
| dual-endpoint Schur | 12 |
| total | 888 |

The artifact records the profile, split, theorem status, source coefficient, and rational outward upper for all 888 number-sector-balanced high-sector entries. It also records the total 6,016-entry routing inventory, the 5,128 excluded unbalanced entries, and the 272-incidence/136-edge scope audit. The regression requires the theorem coefficient map to equal the balanced set exactly, so an unproved routing placeholder cannot silently re-enter.

## Coefficient rounding

Every accepted nonnegative binary64 theorem coefficient is replaced by a strict rational upper on a $10^{-9}$ grid, followed by one additional full grid unit of guard. Thus an accepted source value $c$ is replaced by

$$
\widehat c=\frac{\lceil 10^9c\rceil+1}{10^9}>c.
$$

The largest resulting inflation is below $2\times10^{-9}$. The exact dual-endpoint and final-residual squared coefficients are independently checked to lie below the squares of their ledger values.

## Perron certificate

The occurrence calculation is reconstructed on all 210 at-most-six occupation states using balanced splits only, so the matrix is block diagonal in total signal photon number. Every matrix entry is accumulated with 80-digit `Decimal` arithmetic under `ROUND_CEILING`; every integer square root is advanced to a value whose square is demonstrably no smaller than its argument.

A floating eigensolve supplies only a strictly positive candidate vector $v$. It is not used as a numerical eigenvalue certificate. The committed upper is

$$
\rho(A)\le\max_i\frac{(\widehat A v)_i}{v_i}
$$

by Collatz--Wielandt, where $\widehat A$ is the entrywise outward matrix. This remains valid for any positive candidate vector and therefore does not inherit the eigensolver's rounding error.

## Promise certificate

For a centered biased sign with mean $\beta$, the Kearns--Saul proxy is at most one. The inherited four-link chain proxy is therefore at most

$$
\frac{(1+\beta^2)(1+\beta^4)}{N}.
$$

With $N=4096$ and gap $g=\beta^4-1/4$, the one-sided Chernoff exponent is bounded below by

$$
x=\frac{Ng^2}{2(1+\beta^2)(1+\beta^4)}.
$$

No transcendental floating evaluation is used in the accepted value. For $m=4096$,

$$
e^{-x}\le\left(1+\frac{x}{m}\right)^{-m},
$$

and twice this exact rational upper gives the displayed two-hypothesis promise loss.

## Reproduction and decision

Run:

    /opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_complete_outward_ledger.py --output artifacts/q64_complete_outward_ledger.json
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_complete_outward_ledger.py

The gate passes with much more than the required $10^{-3}$ reserve. This note certifies only the one-batch numerical ledger in the number-sector-incoherent class; `Q64_ADAPTIVE_TREE_FRONTIER_THEOREM.md` supplies the separate outcome-uniform lift and preserves this total exactly. Extending the proof to coherent superpositions of total signal photon number requires certified bounds for the unbalanced physical kernels.
