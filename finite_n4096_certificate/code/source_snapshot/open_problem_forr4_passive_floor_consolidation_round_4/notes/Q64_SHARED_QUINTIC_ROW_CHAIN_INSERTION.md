# q64 shared quintic row/chain insertion

Date: 2026-07-16

Status: rigorous arbitrary-diagonal one-batch coefficients for all 48 previously open degree-twelve split-cubic/split-quintic entries. The theorem raises the q64 high-sector count from 380 to 428 and leaves 460 entries in the four nonquintic classes. It does not prove those 460 entries, intervalize the Perron calculation, or establish the adaptive lift.

## Result

At

$$
q=64,qquad N=q^2=4096,qquad M=4N=16384,
$$

the twelve complement/reversal orbits split into five finite row/chain templates:

| template | orbits | entries | coefficient range |
|---|---:|---:|---:|
| favorable singleton--quintic endpoint, adjacent split cubic | 4 | 16 | $0.0189275747390$--$0.0203737451368$ |
| favorable singleton--quintic endpoint, adjacent whole cubic | 3 | 12 | $0.000263887951092$--$0.000318339767762$ |
| endpoint cubic, quintic, adjacent split cubic | 2 | 8 | $0.0183664518516$--$0.0202103376613$ |
| fixed-one singleton--quintic endpoint, adjacent whole cubic | 2 | 8 | $0.00845369672045$ |
| singleton--whole-cubic endpoint, fixed-four quintic | 1 | 4 | $0.00000535051085909$ |

Thus every coefficient is at most

$$
0.0203737451368,
$$

far below the predeclared common reserve gate $0.410314553367$. With the other 460 open entries retained at their frozen targets, the optimized routing diagnostic becomes

$$
U_{mathrm{route}}=0.323811563171336,qquad {1over3}-U_{mathrm{route}}=0.009521770161998.
$$

This routing value is still not a passive theorem because the other 460 coefficients are targets rather than arbitrary-law bounds.

## Shared arbitrary-law mechanism

Each physical link moment is a cross-Gram kernel of normalized monomial features. After one complete row or chain is normalized by its squared energy, every unused physical link is a unit Schur multiplier. Therefore the argument applies to arbitrary nonnegative row and column diagonal laws; it does not optimize separate laws for different links.

The templates use only three operations:

1. normalize a complete endpoint/chain row by a proved squared-slice bound;
2. decompose incompatible parity records into a finite orthogonal support partition and apply the triangle inequality across sectors; and
3. bound a residual sector by the minimum of its cut-rank, row-incidence, and column-incidence estimates.

Complementation transposes the occurrence matrix and path reversal relabels the features, so each displayed generator proves its full four-entry orbit.

## The decisive fixed-one record split

Let $F_{1,r}$ be the squared singleton--quintic endpoint slice with one fixed quintic cell, after restricting the quintic record on the following cubic link to size $r\in\{1,3,5\}$. Refining the exact endpoint count by column patterns $5$, $4+1$, $3+2$, and $2+2+1$ gives

$$
(F_{1,1},F_{1,3},F_{1,5})
=
(1.253662109375, 55.77392578125, 968.22265625).
$$

Their sum is the inherited exact slice

$$
F_1=1025.250244140625.
$$

The apparent obstruction in the coarse row bound came from charging all of $F_1$ with the record-one cubic--quintic maximum. In fact, the record-five portion cannot couple to a cubic at all. The two compatible middle-link maxima are

$$
m_1={q+2\over q(q-1)(q-2)}=0.000264016897081413,
$$

and

$$
m_3={3\over(q-3)\binom q3}=0.000001180403414671.
$$

For the difficult extreme template this yields

$$
\gamma_{mathrm{ext}}
\le
\sqrt{N\left(F_{1,1}m_1^2+F_{1,3}m_3^2\right)}
=0.018927574738988.
$$

The former universal charge was $0.541036590109$. The improvement is structural: almost all fixed-one endpoint energy lies in a parity sector that either vanishes or couples through the much smaller record-three entry.

For the balanced orientation, the inherited fixed-three endpoint slice $F_3=1.4538457961309523$ already gives

$$
\gamma_{mathrm{bal}}
\le
\sqrt{NF_3m_1^2}
=0.020373745136756.
$$

## The other four templates

When the singleton--quintic endpoint is followed by a whole cubic, no free cubic completion is charged. Fixing the majority endpoint cells gives

$$
\sqrt{F_4}m_1=0.000263887951092
$$

for an extreme cut and

$$
\sqrt{F_3}m_1=0.000318339767762
$$

for a balanced cut. If the singleton lies on the minority side, the safe fixed-one version is

$$
\sqrt{F_1}m_1=0.008453696720453.
$$

For the two endpoint-cubic/quintic/split-cubic orientations, write $r,t\in\{1,3\}$ for the two cubic--quintic link records. The endpoint-cubic amplitude is $1/[q(q-1)]$ when $r=1$ and $1/q$ when $r=3$. For each of the four $(r,t)$ sectors, the code uses the exact one-axis record incidence of each support, relaxes a two-axis incidence to the smaller of its two one-axis incidences, and takes

$$
\min\left\{
q^{\min(|s|,12-|s|)}a_{r,t},
a_{r,t}\sqrt{D_s^{r,t}},
a_{r,t}\sqrt{D_{d-s}^{r,t}}
\right\}.
$$

Summing the four sectors gives $0.0183664518516$ for the extreme orbit and $0.0202103376613$ for the balanced orbit.

For the last singleton--whole-cubic/fixed-four-quintic orbit, the exact record incidences are

$$
D_{4,1}^{(5)}=4092,qquad D_{4,3}^{(5)}=3968.
$$

The two-sector coefficient is

$$
{m_1\sqrt{4092}\over q(q-1)}
+{m_3\sqrt{3968}\over q}
=0.00000535051085909.
$$

## Adaptive requirement and acceptance status

If the remaining 460 entries are eventually proved at their currently frozen targets, retaining the declared $10^{-3}$ numerical allowance leaves the adaptive lift at most

$$
\Delta_{mathrm{adapt}}
\le
{1\over3}-10^{-3}-U_{mathrm{route}}
=0.008521770161998.
$$

Equivalently, a purely multiplicative adaptive comparison would need factor at most

$$
{1/3-10^{-3}\over U_{mathrm{route}}}
=1.026317065637.
$$

These are quantitative design requirements, not an adaptive theorem. The class acceptance gate is now passed only in its one-batch arbitrary-law component. Round 4 completion still requires all of:

- arbitrary-law closure of the remaining 460 entries;
- an outcome-width-, depth-, posterior-, and hard-dose-partition-uniform adaptive recurrence whose evaluated overhead fits the final certified one-batch margin;
- outward-rounded or exact numerical certification; and
- the explicit active six-dose row and experimental resource statement.

## Reproduction

Run:

```text
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_shared_quintic_row_chain_insertion.py --write-artifact
/opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_shared_quintic_row_chain_insertion.py
```

The regression directly enumerates the $q=4$ fixed-one endpoint slice and checks its record-one/three/five decomposition, verifies the exact 48-entry orbit cover, checks the q64 incidences and coefficients, reoptimizes the ledger, and compares the committed artifact byte for byte.
