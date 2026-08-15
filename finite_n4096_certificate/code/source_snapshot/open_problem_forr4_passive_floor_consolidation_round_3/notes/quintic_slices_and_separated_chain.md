# Quintic slices and the separated cubic chain

Date: 2026-07-14

Status: proved endpoint-quintic squared slices, sharp middle-quintic
incidences from two through five fixed cells, and fixed-split coefficients
for the two degree-ten reversal pairs

$$
(5,1,3,1)/(1,3,1,5),\qquad
(3,1,5,1)/(1,5,1,3).
$$

The accepted compatible Perron ledger remains below \(1/3\) after adding
both pairs.  The same calculation also gives a sharp obstruction to the
direct row-energy treatment of the next opposite-endpoint orientation.

## 1. Degree-five endpoint moments

Let \(Q\) be a degree-five endpoint support whose hidden column record has
size one, and let the neighboring singleton match that odd column.  Partition
the five endpoint cells by their column multiplicities.  After the odd label
is matched, injective Walsh averaging gives three possible squared moment
weights:

$$
w_0={1\over q^2},\qquad
w_1={1\over q^2(q-1)^2},\qquad
w_2={4\over q^2(q-1)^2(q-2)^2}.
\tag{1.1}
$$

The column patterns and their weights are:

| multiplicities | condition | squared weight |
|---|---|---:|
| \(5\) | none | \(w_0\) |
| \(4+1\) | four-row XOR zero | \(w_0\) |
| \(4+1\) | otherwise | \(w_1\) |
| \(3+2\) | always | \(w_1\) |
| \(2+2+1\) | the two pair XORs agree | \(w_1\) |
| \(2+2+1\) | distinct pair XORs | \(w_2\) |

This follows from the normalized injective character average over zero, one,
or two even labels.  With one even label, a zero character averages to one
and a nonzero character to \(-1/(q-1)\).  With two distinct nonzero
characters the magnitude is \(2/[(q-1)(q-2)]\).

## 2. Exact endpoint squared slices

Let \(F_k\) be the largest squared moment sum over compatible endpoint
supports containing a fixed \(k\)-cell subset, with the singleton coordinate
also optimized.  Define

$$
Z_4={q(q-1)(q-2)\over24}.
$$

The full-slice counts are

$$
\begin{aligned}
n_5&=q\binom q5,\\
n_{41}&=q^2(q-1)\binom q4,
&h_{41}&=q^2(q-1)Z_4,\\
n_{32}&=q(q-1)\binom q3\binom q2,\\
n_{221}&=q^2\binom{q-1}2\binom q2^2,
&h_{221}&={q^4(q-1)\over4}\binom{q-1}2.
\end{aligned}
\tag{2.1}
$$

Here \(h_{41}\) and \(h_{221}\) count the exceptional supports receiving
weight \(w_0\) and \(w_1\), respectively.  Therefore

$$
\begin{aligned}
F_0={}&(n_5+h_{41})w_0\\
&+(n_{41}-h_{41}+n_{32}+h_{221})w_1\\
&+(n_{221}-h_{221})w_2,
\\
F_1={}&{5F_0\over q^2}.
\end{aligned}
\tag{2.2}
$$

The maximizing fixed pair and triple lie in one column.  Direct extension
counting gives

$$
\begin{aligned}
F_2={}&
\binom{q-2}3w_0\\
&+q(q-1)\left[
\left({q\over2}-1\right)w_0+
\left(\binom{q-2}2-{q\over2}+1\right)w_1
\right]\\
&+\left[(q-2)(q-1)\binom q2+(q-1)\binom q3\right]w_1\\
&+q(q-1)(q-2)\left[
{q\over2}w_1+
\left(\binom q2-{q\over2}\right)w_2
\right],
\\[2mm]
F_3={}&
\binom{q-3}2w_0+
q(q-1)\left[w_0+(q-4)w_1\right]+
(q-1)\binom q2w_1,
\\
F_4={}&1-{4\over q^2},
\qquad
F_5={1\over q^2}.
\end{aligned}
\tag{2.3}
$$

For pairs there are three equality types: one column, one row, or disjoint
row and column.  For triples there are six types obtained from column
partitions \(3\), \(2+1\), and \(1+1+1\) and the possible row equalities.
Substitution in each orbit proves the stated maxima.  At \(q=32\), the
nonmaximal pair types both give \(3249/7936\), while the one-column type
gives \(159457/7936\).  The nonmaximal triple types give at most
\(8931/492032\), while the one-column type gives \(22365/15872\).

The target-order table is

$$
\boxed{
(F_0,\ldots,F_5)=
\left(
52685,\,
{263425\over1024},\,
{159457\over7936},\,
{22365\over15872},\,
{255\over256},\,
{1\over1024}
\right).}
\tag{2.4}
$$

Exact signed-permutation enumeration at \(q=4\) reproduces

$$
\left(17,{85\over16},{25\over12},{7\over8},{3\over4},{1\over16}\right).
\tag{2.5}
$$

## 3. Middle-quintic incidences

A degree-five middle support between singleton neighbors must have record one
on both axes.  If \(D_k^{(5)}\) is the maximum number of such supports
containing a fixed \(k\)-cell subset, then

$$
\boxed{
\begin{aligned}
D_2^{(5)}&=(q-2)(q-1)(5q-4),\\
D_3^{(5)}&=2(q-2)(2q-1),\\
D_4^{(5)}&=q^2-4,\\
D_5^{(5)}&=1.
\end{aligned}}
\tag{3.1}
$$

Fixed pairs maximize in one row or column.  Fixed triples maximize on an
L-shape.  The complete pair/triple isomorphism count is

| fixed cells | bipartite type | compatible extensions |
|---:|---|---:|
| 2 | one row or one column | \((q-2)(q-1)(5q-4)\) |
| 2 | matching | \((q-2)(15q-22)\) |
| 3 | three-star | \(3(q-1)\) |
| 3 | L-shape | \(2(q-2)(2q-1)\) |
| 3 | two-star plus disjoint edge | \(7q-11\) |
| 3 | matching | \(9\) |

These formulas follow by classifying the equality pattern of the remaining
three or two edges.  The test suite enumerates every extension of each
isomorphism type exactly for \(q=4,\ldots,8\), in addition to enumerating
every \(q=4\) support.  The displayed comparison proves the maxima for
\(q\ge4\).  Fixed four cells maximize on a four-cycle.  At \(q=4\),
exhaustive support enumeration gives

$$
(D_0^{(5)},\ldots,D_5^{(5)})
=(1008,315,96,28,12,1),
\tag{3.2}
$$

and reproduces (3.1).  The zero- and one-cell source bounds deliberately use
the larger one-record incidence family; they are safe but not asserted
sharp.

The four-cycle value is a genuine local obstruction.  Fixing a rectangle
leaves \(q^2-4\) possible fifth edges, every one of which has record one on
both axes.  For each completion, both neighboring singleton moment entries
attain \(1/[q(q-1)]\).  Thus neither a smaller local entry bound nor a smaller
\(D_4^{(5)}\) can improve that cut.

## 4. The resolved reversal pair

For

$$
a=(5,1,3,1),
$$

the first block is an endpoint record-one quintic and the third block is a
forced L-shape.  The three moment entries are bounded by

$$
{1\over q},\qquad {1\over q(q-1)},\qquad
{1\over q(q-1)}.
\tag{4.1}
$$

Let \(D_k^E\) be the exact endpoint-quintic one-record incidence, \(D_k^L\)
the L-shape incidence, and \(D_k^S=(q^2,1)_k\) the singleton incidence.  For
a fixed split \(s\), the proved coefficient is

$$
\gamma_{a,s}\le
\min\left\{
{q^{\min(|s|,10-|s|)}\over q^3(q-1)^2},
{ \sqrt{\prod_j D_{s_j}^{(j)}}\over q^3(q-1)^2},
{ \sqrt{\prod_j D_{a_j-s_j}^{(j)}}\over q^3(q-1)^2}
\right\}.
\tag{4.2}
$$

Path reversal gives \((1,3,1,5)\).  At \(q=32\), the largest coefficient is

$$
\gamma_{\max}=0.0246729924519
\tag{4.3}
$$

at split \((4,1,0,0)\) and its complements/reversals.

## 5. Updated finite-size ledger

At the target size

$$
q=32,\qquad N=1024,
$$

the slightly better rational attenuation is \(\beta=781/1000\).  Adding the
four resolved profiles to all previously accepted sectors gives

$$
F_{\rm known}<0.311361620285.
\tag{5.1}
$$

The \((5,1,3,1)/(1,3,1,5)\) contributions are
\(0.001722015303\) each and the
\((3,1,5,1)/(1,5,1,3)\) contributions are \(0.001345386692\) each.
The subgaussian promise loss is \(0.019452044176\), so

$$
\boxed{
F_{\rm known}+2\epsilon_\beta
<0.330813664461
<{1\over3}.}
\tag{5.2}
$$

The remaining margin is \(0.002519668873\).  A local attenuation scan places
the optimum near \(\beta=0.781\), which is why the current certificate uses
the exact rational \(781/1000\).

For every still-open high-degree profile, a common provisional coefficient
\(1/588\) passes with total \(0.333321204449\), while \(1/584\) fails with
total \(0.333338429976\).  The transition is near \(1/585\).  These are
targets, not proved coefficients.

Update: this paragraph is the pre-obstruction checkpoint.  After inserting
the physical mixed-orbit cuts and finite-tilt promise repair, the comparable
common threshold is \(0.000735032568\), about \(1/1360.48\).  See
`repaired_open_profile_budget.md` for the current target.

## 6. A route obstruction and a second resolved pair

### 6.1 Opposite cubic and quintic endpoints

For \((3,1,1,5)\), directly pairing the cubic endpoint slices \(E_k\) with
the quintic slices \(F_\ell\) gives

$$
\max_{k,\ell}
\sqrt{\min\{E_kF_\ell,E_{3-k}F_{5-\ell}\}}
=0.780328871525.
\tag{6.1}
$$

The maximum occurs at \((k,\ell)=(2,2)\).  This rules out a direct repeat of
the double-endpoint row-energy contraction.  The profile needs a compound
weighted spectral or frame argument.

That compound calculation has now started.  The exact endpoint moments are
xor-labelled Walsh rows, and every fixed pair/triple translation-orbit block
has a factored row Gram.  The aligned full block is \(0.00344731635\) at
\(q=32\).  The follow-up exact mixture calculation finds a valid physical
law with coefficient \(0.0395939553\), which makes the independent scalar
ledger fail by \(0.0008501212\).  See
`opposite_endpoint_orbit_factorization.md` and
`opposite_endpoint_mixed_orbit_obstruction.md`.  The route now needs a
promise/accepted-sector repair, a joint contraction, or a hard-instance
pivot.

### 6.2 Cubic endpoint and middle quintic

For \(b=(3,1,5,1)\), use the endpoint-cubic record-one incidences, a
singleton, the middle-quintic incidences from Section 3, and the final
singleton.  The same link-entry product \(1/[q^3(q-1)^2]\) and the same
rank/incidence minimum as (4.2), with these four families and profile \(b\),
give maximum coefficient

$$
0.0332010766872
\tag{6.2}
$$

at split \((2,0,4,0)\).  The four-cycle saturator shows why further local
incidence counting cannot improve the decisive cut.  Nevertheless, the
sharp coefficient fits the actual compatible Perron ledger, so this pair
is resolved.  A chain-aware contraction would only be needed later to
recover margin for other sectors.

## 7. Consequence

Two leading reversal pairs are closed and the partial theorem still fits.
The next useful work is not another independent row-energy estimate.  It is
one of:

1. a compound contraction for \((3,1,1,5)/(5,1,1,3)\);
2. a record-aware contraction for the adjacent
   \((1,1,3,5)/(5,3,1,1)\) chain; or
3. an alternative hard-instance calculation if these physical saturators
   force the complete ledger above threshold.

Reproduction:

- `searches/occupation_compatible_sector_optimization.py`;
- `tests/occupation_compatible_sector_optimization.py`; and
- `./run_round3_checks.sh`.
