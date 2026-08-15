# Whole-cubic decorated completion rows at q=64

Date: 2026-07-16

Status: rigorous arbitrary-diagonal one-batch coefficients for 16 entries in four complement/reversal orbits. This raises the q64 theorem count from 304 to 320. It does not prove the remaining 568 entries, intervalize the ledger, or establish the adaptive lift.

## Result

The four canonical cuts and coefficients are

| cut | coefficient |
|---|---:|
| $(3,1,3,5):(0,0,2,4)$ | $0.000196155632204$ |
| $(1,3,5,3):(0,2,4,0)$ | $0.0125539604611$ |
| $(3,1,5,3):(0,0,4,2)$ | $0.00333622329795$ |
| $(1,5,3,3):(0,4,2,0)$ | $0.213518291069$ |

All 16 entries have extreme quintic splits. They currently have zero Perron sensitivity, so inserting their proved values leaves the optimized routing total unchanged:

$$
P_{\rm total}=0.325283979608,
\qquad
{1\over3}-P_{\rm total}=0.00804935372583.
$$

This zero numerical change does not make the theorem optional: the entries must still be controlled before the one-batch ledger is complete. The remaining quintic inventory is 108 entries, with 72 extreme and 36 balanced cuts. Its local-slice proxy remains $0.325988623739$, leaving $0.00634470959400$ beyond the declared $10^{-3}$ allowance.

## Two reusable scalar rows

The first pair of cuts uses the scalar feature

$$
L=v_3(C)M_{35}(C,S),
$$

with a fixed cubic pair and four fixed quintic cells. Its record-sector energies obey

$$
\begin{aligned}
R_1^{(3)}
&={2(q-1)\over(q-1)^2}(q^2-4)m_1^2
=9.05500128879\times10^{-6},\\
R_3^{(3)}
&=NE_2(N-4)m_3^2
=0.000148546921969,
\end{aligned}
$$

where $E_2=0.0153847346230$, $m_1=0.000264016897081$, and $m_3=1/\binom q3$.

The second pair uses

$$
L=v_5(S)M_{53}(S,C),
$$

with four fixed quintic cells and a fixed cubic pair. Its bounds are

$$
\begin{aligned}
R_1^{(5)}
&=(q^2-4)2(q-1)m_1^2
=0.0359393001152,\\
R_3^{(5)}
&=NF_4(N-2)m_3^2
=0.00965076050582,
\end{aligned}
$$

where $F_4=0.9990234375$.

## Why the whole cubic does not spoil the row bound

When the whole cubic precedes a singleton and both lie on the column side, its endpoint moment is a column-only scalar

$$
M_{31}(T,b)=v_3(T)H_N(\xi(T),b),
$$

with magnitude at most $1/q$. The remaining endpoint Walsh link is a repeated unitary matrix. Thus the residual coefficient is $1/q$, giving

$$
{1\over q}\sqrt{R_1^{(3)}+R_3^{(3)}}
=0.000196155632204
$$

and

$$
{1\over q}\sqrt{R_1^{(5)}+R_3^{(5)}}
=0.00333622329795.
$$

When the whole cubic is the trailing block, its link to the adjacent support is a cross Gram of unit character features. Its trace-class Schur-multiplier norm is at most one. Multiplying the scalar completion row by this link therefore gives

$$
\sqrt{R_1^{(3)}+R_3^{(3)}}
=0.0125539604611
$$

and

$$
\sqrt{R_1^{(5)}+R_3^{(5)}}
=0.213518291069.
$$

Both arguments retain the cubic and quintic within-block distinctness masks inside the complete scalar row. No product-law or invariant-law assumption is used.

## Regression scope

Complete $q=4$ enumeration checks both scalar-row orientations, including record-one and record-three maxima separately. The q64 theorem follows from the explicit completion counts, exact endpoint slices, middle-moment bounds, and Schur-feature reductions above.

Reproduce with:

    /opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_whole_cubic_decorated_row_insertion.py --write-artifact
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_whole_cubic_decorated_row_insertion.py
