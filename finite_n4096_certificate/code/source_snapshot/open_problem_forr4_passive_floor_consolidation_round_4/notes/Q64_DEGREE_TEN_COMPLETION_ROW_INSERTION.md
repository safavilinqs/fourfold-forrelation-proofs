# Degree-ten completion-row insertion at q=64

Date: 2026-07-16

Status: rigorous arbitrary-diagonal one-batch coefficients for 12 entries in three complement/reversal orbits. This raises the q64 theorem count from 292 to 304. It does not prove the remaining 584 entries, intervalize the ledger, or establish the adaptive lift.

## Result

The three canonical cuts are

$$
\begin{aligned}
(1,1,3,5)&:(0,0,2,3),\\
(1,1,5,3)&:(0,0,4,1),\\
(1,3,5,1)&:(0,1,4,0).
\end{aligned}
$$

Their arbitrary-law coefficients at $q=64$ are, respectively,

$$
\boxed{
0.00861554231015,qquad
0.0751888832423,qquad
0.0754041939294.
}
$$

Complement and path reversal give all 12 entries. Eight have extreme quintic splits and four have balanced quintic splits.

After insertion, the routing values are

$$
\begin{aligned}
\beta&=0.746146955731,\\
P_{\rm Perron}&=0.308223500018,\\
P_{\rm promise}&=0.0170604795898,\\
P_{\rm total}&=0.325283979608.
\end{aligned}
$$

The routing margin is $0.00804935372583$, an improvement of $0.00382322812438$ over the 292-entry ledger. The remaining quintic inventory is 124 entries: 88 extreme and 36 balanced. Charging those entries at the existing endpoint local-slice values gives total $0.325988623739$ and leaves $0.00634470959400$ beyond the declared $10^{-3}$ allowance.

## Common scalar-row mechanism

For an endpoint support $T$ of odd degree, the signed-permutation endpoint moment has the exact form

$$
M_{1,|T|}(a,T)=v_{|T|}(T)H_N(a,\xi(T)),
$$

where $\xi(T)$ is the support xor. Extract the scalar completion row containing the endpoint amplitude, the middle moment $M_{35}$ or $M_{53}$, and any second endpoint amplitude. The normalized scalar row is a Schur feature. The residual Walsh matrix has arbitrary-diagonal weighted trace norm at most one by duplicate compression and Schatten Hölder.

For the first two cuts, the first two singleton blocks lie entirely on the column side. Their residual factor is

$$
H_N(a,b)H_N(\xi(T),b).
$$

After removing column-only characters, this is $q^{-1}$ times a repeated normalized Walsh matrix. Its arbitrary-diagonal coefficient is therefore $1/q$. The third cut retains the two endpoint Walsh factors and gives the accepted repeated, column-twisted $H_N\otimes H_N$ residual with coefficient one.

The middle moment is split into record-one and record-three sectors. At $q=64$ the bounds are

$$
m_1={q+2\over q(q-1)(q-2)}=0.000264016897081,
$$

and

$$
m_3={1\over\binom q3}=0.0000240015360983.
$$

The second formula is valid throughout the record-three sector for $q\ge8$: the one-even-group correction is at most $3/(q-3)\le1$.

## Two leading singletons, followed by cubic and quintic

For $(1,1,3,5):(0,0,2,3)$, fix a cubic pair and a quintic triple. In record one, there are at most $2(q-1)$ compatible cubic completions, each with endpoint amplitude $1/(q-1)$, and $2(q-2)(2q-1)$ compatible quintic completions. Thus

$$
R_1={2(q-1)\over(q-1)^2}\,2(q-2)(2q-1)m_1^2
=0.0000348480352629.
$$

In record three, the exact cubic fixed-pair endpoint slice is

$$
E_2=0.0153847346230.
$$

Using all $\binom{N-3}{2}$ quintic pair completions is a safe relaxation, so

$$
R_3=NE_2\binom{N-3}{2}m_3^2
=0.304001275810.
$$

Multiplying the row coefficient by the residual factor $1/q$ gives

$$
{1\over q}\sqrt{R_1+R_3}=0.00861554231015.
$$

## Two leading singletons, followed by quintic and cubic

For $(1,1,5,3):(0,0,4,1)$, fix four quintic cells and one cubic cell. In record one, the exact incidence bounds are $q^2-4$ quintic completions and $3(q-1)^2$ cubic completions. In record three, the total completion bounds are $N-4$ and $\binom{N-1}{2}$. Consequently

$$
\begin{aligned}
R_1&=(q^2-4)3(q-1)^2m_1^2
=3.39626386089,\\
R_3&=(N-4)\binom{N-1}{2}m_3^2
=19.7599321357.
\end{aligned}
$$

The same $1/q$ Walsh residual gives

$$
{1\over q}\sqrt{R_1+R_3}=0.0751888832423.
$$

## Double endpoint with fixed counts one and four

For $(1,3,5,1):(0,1,4,0)$, retain both endpoint Walsh factors. The record-one cubic incidence and amplitude combine to $3$, while the fixed-four quintic incidence is $q^2-4$. Hence

$$
R_1=3(q^2-4)m_1^2=0.000855697621791.
$$

The exact endpoint slices are

$$
E_1=0.500244140625,
\qquad
F_4=0.9990234375.
$$

The record-three scalar-row energy is therefore

$$
R_3=N^2E_1F_4m_3^2=0.00483009484036.
$$

The $H_N\otimes H_N$ residual has coefficient one, giving

$$
\sqrt{R_1+R_3}=0.0754041939294.
$$

## Regression scope

Complete $q=4$ enumeration constructs every scalar completion row for all three orientations and checks record-sector maxima separately. A randomized arbitrary-law stress test checks the claimed $1/q$ double-singleton Walsh residual. The q64 theorem itself follows from the explicit incidence counts, endpoint-slice identities, record-sector moment bounds, and Schur-feature factorization above.

Reproduce with:

    /opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_degree_ten_completion_row_insertion.py --write-artifact
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_degree_ten_completion_row_insertion.py
