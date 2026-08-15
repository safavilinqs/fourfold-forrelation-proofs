# Triple-cubic chained-slice contraction

Date: 2026-07-14

Status: proved fixed-split coefficients for the four leading degree-ten
triple-cubic profiles.  Together with the previously proved high endpoint
profiles, the compatible Perron ledger is below \(1/3\) after promise
conditioning at \(N=1024\).  This is still a partial one-batch result:
the other degree-ten/twelve profiles and adaptivity remain open.

## 1. Result

The four profiles are

$$
(3,3,3,1),\quad(1,3,3,3),\quad
(3,1,3,3),\quad(3,3,1,3).
$$

The old maximum-entry/incidence estimate is unusable: on these four profiles
alone it can raise the degree-eight ledger plus promise loss to about
\(0.826\).  The defect is that it charges every record-three output as if it
simultaneously had the largest moment entry and the largest final singleton
coherence.

The replacement keeps squared slice energies through the whole three-link
chain.  At

$$
q=32,\qquad N=q^2=1024,\qquad \beta={25\over32},
$$

the ledger containing

- every profile through degree eight;
- the proved degree-ten/twelve single-endpoint profiles; and
- all four triple-cubic profiles

has Perron upper

$$
F_{\rm known}<0.306046304645.
$$

The subgaussian two-hypothesis promise loss is

$$
2\epsilon_\beta<0.018785599766,
$$

so

$$
\boxed{
F_{\rm known}+2\epsilon_\beta
<0.324831904411
<{1\over3}.}
\tag{1.1}
$$

The remaining margin is \(0.008501428922\).

## 2. Record classification

Write \(r_i\in\{1,3\}\) for the odd record size on hidden link \(i\).
Singleton neighbors force the adjacent record to equal one.

For \((3,3,3,1)\), the possible record triples and support families are

| record triple | block families |
|---|---|
| \(111\) | endpoint record one, L-shape, L-shape, singleton |
| \(311\) | endpoint record three, \((3,1)\)-star, L-shape, singleton |
| \(131\) | endpoint record one, \((1,3)\)-star, \((3,1)\)-star, singleton |
| \(331\) | endpoint record three, matching cubic, \((3,1)\)-star, singleton |

For \((3,1,3,3)\), only \(111\) and \(113\) occur.  Reversal gives the
other two profiles.

The matching-cubic incidence degrees are

$$
D_k^{M}=\left(
6\binom q3^2,
2\binom{q-1}2^2,
(q-2)^2,
1
\right)_k.
\tag{2.1}
$$

Direct \(q=4\) enumeration verifies the record classification and all
incidence degrees.

## 3. Three local squared-slice tables

### 3.1 Endpoint next to a singleton

If \(k\) endpoint cells are fixed and the other endpoint cells are summed,
the exact squared slice is

$$
E_k=\left(
{q^2+2\over6},
{q^2+2\over2q^2},
{q^2-2q+2\over q^2(q-1)},
{1\over q^2}
\right)_k.
\tag{3.1}
$$

This is the earlier endpoint-slice theorem.

### 3.2 A record-\((1,3)\) star next to a singleton

A star is type A when all three record-one labels coincide and type B when
one label is odd and a second label is repeated.  Its link entry has squared
magnitude \(q^{-2}\) in type A and \(q^{-2}(q-1)^{-2}\) in type B.

Counting extensions of a fixed \(k\)-cell subset gives

$$
T_k=\left(
{q^2-4\over6},
{q^2-4\over2q^2},
{q-2\over q(q-1)},
{1\over q^2}
\right)_k.
\tag{3.2}
$$

For example, a fixed cell has

$$
\binom{q-1}{2}
$$

type-A extensions and

$$
{3\over2}(q-1)^2(q-2)
$$

type-B extensions.  Weighting these two counts gives \(T_1\).

### 3.3 A pure record-three output

For a fixed cubic input with record three, let \(C_k\) be the largest
squared moment sum over all record-three cubic outputs containing a fixed
\(k\)-cell subset.  The pure signed-permutation moment is

$$
M(A,B)={1\over(q)_3}
\sum_{\pi\in S_3}\prod_{i=1}^3
\chi_{x_i}(u_{\pi(i)})\chi_{v_{\pi(i)}}(y_i).
\tag{3.3}
$$

Expanding the square and summing the free Walsh labels gives unnormalized
squared-permanent sums

$$
q^4(q-1)(q-2),\quad
3q^2(q-1)(q-2),\quad
12q(q-2),\quad
36
\tag{3.4}
$$

when zero, one, two, or three output cells are fixed.  Division by
\((q)_3^2\) yields

$$
C_k=\left(
{q^2\over(q-1)(q-2)},
{3\over(q-1)(q-2)},
{12\over q(q-1)^2(q-2)},
{36\over q^2(q-1)^2(q-2)^2}
\right)_k.
\tag{3.5}
$$

The regression enumerates the exact \(q=4\) moment matrices.  Independent
direct evaluations at \(q=8,16\) reproduce (3.5).

## 4. Retaining the final singleton suppression

The output of a record-three link can be:

- a type-A star, whose summed singleton-link squared energy is one;
- a type-B star, whose energy is \((q-1)^{-2}\); or
- a matching cubic, whose singleton link is zero.

The type-A star incidence degrees are

$$
A_k=\left(q\binom q3,\binom{q-1}2,q-2,1\right)_k,
\tag{4.1}
$$

and every record-three entry is at most

$$
\kappa_3={1\over\binom q3}.
$$

Consequently a record-three slice followed by its actual final singleton
link is bounded by

$$
\boxed{
W_k={C_k\over(q-1)^2}
+\left(1-{1\over(q-1)^2}\right)\kappa_3^2A_k.}
\tag{4.2}
$$

This is the decisive gain.  At \(q=32\),

$$
(C_0,C_1,C_2,C_3)
=(1.10107527,\,0.00322581,\,1.30073\!\times10^{-5},
4.06478\!\times10^{-8}),
$$

whereas

$$
(W_0,W_1,W_2,W_3)
=(0.00759066,\,2.22383\!\times10^{-5},
1.23170\!\times10^{-6},\,4.06478\!\times10^{-8}).
$$

The earlier bound replaced every \(W_k\) by \(C_k\).

## 5. The endpoint-to-star first link

Let \(P_\ell\) bound the squared output slice over stars containing a fixed
\(\ell\)-cell subset, for one fixed endpoint support.  Type-A star entries
are at most \(1/q\).  Type-B entries are at most

$$
b_1^{1/2}={q+2\over q(q-1)(q-2)}.
$$

Separating their incidences gives

$$
P_\ell=A_\ell/q^2+B_\ell b_1,
\tag{5.1}
$$

for \(\ell<3\), with \(P_3=1/q^2\), where

$$
B_\ell=\left(
q^2(q-1)\binom{q-1}2,
{3\over2}(q-1)^2(q-2),
(q-1)(q-2),0
\right)_\ell.
\tag{5.2}
$$

In the reverse direction, let \(Q_k\) bound endpoint inputs containing a
fixed \(k\)-cell subset for one fixed star.  A type-A fixed star gives the
exact table \(E_k\).  For a type-B fixed star, the full energy is

$$
Q_0^{B}=D_0^E b_0+X(b_1-b_0),
\quad
b_0={1\over q^2(q-1)^2},
\quad
X={q^3(q-2)^2\over8}.
\tag{5.3}
$$

Here \(X\) counts endpoint supports satisfying the two Walsh orthogonality
conditions that produce the exceptional entry.  Thus

$$
Q_k=\max\left\{E_k,
\min\{Q_0^B,b_1D_k^E\}\right\}.
\tag{5.4}
$$

For fixed endpoint and star parts, the complete first-link squared sum is
therefore at most

$$
J_{k\ell}=\min\{D_k^E P_\ell,D_\ell^{\rm star}Q_k\}.
\tag{5.5}
$$

At \(q=32\), (5.4) simplifies numerically to \(Q_k=E_k\), but the code
uses the general maximum in (5.4).

## 6. The two chained contractions

### 6.1 Separated profile \((3,1,3,3)\)

For the \(113\) record sector and split \(s\), condition on the singleton in
block two.  The other three squared slices separate.  If that singleton is
not selected, it is summed over all \(q^2\) cells.  Hence the row energy is

$$
R_{113}(s)=q^{2(1-s_2)}E_{s_1}T_{s_3}C_{s_4}.
\tag{6.1}
$$

The arbitrary-diagonal coefficient is

$$
\gamma_{113,s}
\le\sqrt{\min\{R_{113}(s),R_{113}(a-s)\}}.
\tag{6.2}
$$

The \(111\) sector is added separately using its endpoint--singleton--L--
endpoint maximum entry and cut-dependent incidence degrees.

### 6.2 End-chain profile \((3,3,3,1)\)

For the \(131\) record sector, first sum the endpoint-to-star link by (5.5),
then sum the record-three link together with the final singleton by (4.2).
The row energy is

$$
R_{131}(s)=J_{s_1s_2}W_{s_3}q^{-2s_4},
\tag{6.3}
$$

and

$$
\gamma_{131,s}
\le\sqrt{\min\{R_{131}(s),R_{131}(a-s)\}}.
\tag{6.4}
$$

The \(111\), \(311\), and \(331\) sectors already have adequate
entry/incidence bounds.  Their maximum entries are respectively

$$
{q+2\over q^3(q-1)^2(q-2)},\qquad
{1\over\binom q3q^2(q-1)},\qquad
{1\over q\binom q3^2}.
\tag{6.5}
$$

Adding (6.4) and these three disjoint record sectors gives the fixed-split
coefficient.  Path reversal proves the other orientation.

## 7. What remains

The next step after (1.1) is recorded in
`notes/quintic_slices_and_separated_chain.md`.  It proves both separated
cubic--quintic reversal pairs, derives the exact endpoint-quintic slice
table, and keeps the enlarged partial ledger below threshold.

After charging both pairs at their proved coefficients, a common provisional
coefficient \(1/588\) for every still-open profile passes and \(1/584\)
fails.  The transition is near \(1/585\).  The new note also records the
sharp local obstruction for the opposite-endpoint orientation.

This is a historical pre-obstruction checkpoint.  The first repaired common
threshold after one physical forced cut was about \(1/1360.48\).  After the
second certified witness and hybrid promise repair, the current common
target is about \(1/2604.1\).  See `repaired_open_profile_budget.md` and
`transposed_dominant_class_and_hybrid_repair.md`.

The current leading targets are maintained in `OPEN_PROBLEMS.md`; older
local contraction failures remain documented in the quintic note.

Even a complete one-batch ledger would not yet prove the passive lower
bound.  The adaptive lift remains a separate theorem gate.

Reproduction:

- `searches/occupation_compatible_sector_optimization.py`;
- `tests/occupation_compatible_sector_optimization.py`; and
- `./run_round3_checks.sh`.
