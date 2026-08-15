# Chained repair of the adjacent record-three sector

Date: 2026-07-15

Status: the squared-slice contraction in Sections 2--4 is proved and tested.
The all-open \(1/q\) ledger in Section 6 is a theorem target, not a claim that
the remaining coefficients already obey that bound.

## 1. Why this became the lead result

The adjacent cubic--quintic investigation initially appeared to require a
delicate shared-law improvement by a factor below \(0.710804\).  A repair
scorecard showed that a reduction of less than one half percent in the
already accepted adjacent double-cubic family would instead recover the old
scalar ledger.  Auditing that family exposed a much larger avoidable loss:
its record-three sector still used an entry-incidence estimate even though
the two consecutive squared-slice tables were already exact.

The chained estimate below replaces that loose term.  It lowers the accepted
degree-eight ledger enough that the former adjacent compound gate is no
longer decision-critical.

## 2. Canonical record-three chain

Consider the record-three sector of the degree-eight profile

$$
(1,1,3,3),
$$

in canonical order.  The first two blocks are singletons, the third support
is a record-\((1,3)\) cubic star, and the fourth is a pure record-three cubic
endpoint.  Its kernel is

$$
K(a,b,S,E)=M_{11}(a,b)M_{13}(b,S)M_{33}(S,E).
$$

For a cut \(s=(s_0,s_1,s_2,s_3)\), let \(s_j\) be the number of selected
cells in block \(j\), and put \(u=s_0+s_1\).

The exact singleton--star squared-slice energies are

$$
T_0={q^2-4\over6},\quad
T_1={q^2-4\over2q^2},\quad
T_2={q-2\over q(q-1)},\quad
T_3={1\over q^2}.
\tag{2.1}
$$

For a fixed input cubic, the exact pure record-three output slices are

$$
C_0={q^2\over(q-1)(q-2)},\quad
C_1={3\over(q-1)(q-2)},
$$

$$
C_2={12\over q(q-1)^2(q-2)},\quad
C_3={36\over q^2(q-1)^2(q-2)^2}.
\tag{2.2}
$$

These tables were already proved independently and checked by complete
\(q=4\) moment enumeration.

## 3. Chained squared-slice theorem

Fix the selected partial supports.  First sum the squared \(M_{33}\) entries
over endpoint completions.  Equation (2.2) bounds this uniformly in the
completed star \(S\).  The remaining weighted sum over star completions is
then exactly controlled by (2.1).  Finally,

$$
|M_{11}(a,b)|^2={1\over q^2}
$$

and summing or fixing the two singleton cells contributes
\(q^{2(1-u)}\).  Therefore the row squared energy is at most

$$
R_s=q^{2(1-u)}T_{s_2}C_{s_3}.
\tag{3.1}
$$

Applying the same argument to the complementary cut gives

$$
C_s=q^{2(u-1)}T_{3-s_2}C_{3-s_3}.
\tag{3.2}
$$

The two Schatten factorizations yield the fixed-split coefficient

$$
\boxed{
\gamma_s^{\rm chain}
\le \sqrt{\min\{R_s,C_s\}}.
}
\tag{3.3}
$$

The implementation takes the minimum of (3.3) and the inherited
entry-incidence bound.  Reversing the chain proves the same formula for
\((3,3,1,1)\).

## 4. Direct validation and the repaired cut

At \(q=4\), the test constructs the actual squared singleton--star matrix,
the actual squared record-three matrix, and every relevant endpoint
incidence slice.  Their composed slice is bounded by the product
\(T_{s_2}C_{s_3}\), independently validating the only nontrivial chaining
step.

At \(q=32\), the formerly dominant split orbit contains

$$
(1,1,3,3):(0,1,2,1)
$$

and its cut/complement/reversal images.  On that split,

$$
\gamma^{(1)}=0.001957449200,
$$

while the old record-three incidence term was

$$
\gamma^{(3)}_{\rm old}=0.040812609606.
$$

The chained theorem gives instead

$$
\boxed{
\gamma^{(3)}_{\rm chain}=0.002545238140,
}
$$

so the complete record-one plus record-three coefficient becomes

$$
0.004502687340.
\tag{4.1}
$$

This is not a numerical fit; it is the evaluation of (2.1)--(3.3).

## 5. Consequence for the accepted partial ledger

With every sector through degree eight and the standard promise theorem, at
\(\beta=313/400\) the repaired compatible-occupation ledger is

$$
F_{\le8}=0.232480992267,
$$

$$
F_{\le8}+2\epsilon_\beta
=0.248224749551.
\tag{5.1}
$$

The earlier accepted value was about \(0.297666\).  With the already proved
high endpoint, triple-cubic, and separated cubic--quintic sectors at
\(\beta=781/1000\), the partial total is now

$$
0.279758546919,
\tag{5.2}
$$

leaving \(0.053574786414\) before the unresolved high-degree profiles.

## 6. New finite-size completion target

For route selection, preserve the three known physical-orbit diagnostic
values, including the conservative adjacent combined coefficient
\(0.016866245904\), and assign coefficient \(1/q=1/32\) to every other
unresolved degree-ten/twelve split.  Under the extended promise theorem and
optimized attenuation, this gives

$$
\beta=0.779698447178,
$$

$$
\boxed{
F_{\rm coarse}+2\epsilon_\beta
=0.322669154028<{1\over3},
}
\tag{6.1}
$$

with slack

$$
0.010664179305.
\tag{6.2}
$$

Equation (6.1) is a diagnostic theorem target.  It does **not** prove that
the 6,004 unforced open split entries are at most \(1/q\), and the physical
witness values inserted into this routing ledger are not arbitrary-law upper
bounds.  Its consequence is strategic: a coarse classification theorem is
now more valuable than a delicate improvement of one adjacent witness.

The immediate one-batch problem is therefore:

1. identify which remaining record chains already satisfy a fixed-split
   \(1/q\) theorem by Bessel, rank--Frobenius, or chained slices;
2. isolate and treat any exceptions, with the opposite-endpoint witness as a
   known warning that a universal \(1/q\) statement is false;
3. insert rigorous arbitrary-physical-law bounds for every open split; and
4. intervalize the resulting 210-state Perron certificate.

The adaptive lift remains a separate theorem after this one-batch gate.

## 7. Adjacent mixed-law stress test

The investigation that led to the repair remains useful negative evidence.
An exact \(q=4\) Fourier reduction was independently checked against direct
mixed occurrence matrices.  Its optimized symmetric vertical-triple family
has coefficient \(0.185024058902\) at \(q=4\), but a closed 25-frequency-class
formula gives only

$$
0.000156597909
$$

at \(q=32\).  It is therefore not a physical obstruction at realistic size.
A candidate type-block rank--Frobenius architecture was also too loose:
its sampled \(q=32\) diagnostic was about \(0.02555\).  Because that value is
sampled and uses unproved within-type averaging, it is evidence to stop that
architecture, not an upper bound.

Reproduction:

- `searches/occupation_compatible_sector_optimization.py`;
- `searches/repaired_open_profile_budget.py`;
- `searches/adjacent_combined_repair_scorecard.py`;
- `searches/adjacent_cubic_quintic_mixed_orbit_q4.py`;
- `searches/adjacent_vertical_triple_symmetric_witness.py`;
- `searches/adjacent_compound_frame_type_bound.py`;
- `tests/occupation_compatible_sector_optimization.py`;
- `tests/adjacent_mixed_orbit_diagnostics.py`; and
- `tests/repaired_open_profile_budget.py`.
