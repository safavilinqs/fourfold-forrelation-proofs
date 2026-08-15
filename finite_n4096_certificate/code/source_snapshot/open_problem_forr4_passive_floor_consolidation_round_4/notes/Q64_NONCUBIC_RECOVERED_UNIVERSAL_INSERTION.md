# q64 noncubic and recovered-universal insertion

Date: 2026-07-16

Status: mixed and partially quarantined. The specialized 16-entry middle
degree-seven endpoint theorem is preserved. The other 124 noncubic claims and
all 96 recovered-universal claims rely on the invalid unmasked coefficient-one
step. The 664-entry count and downstream routing values are withdrawn. See
`Q64_MASKED_UNIVERSAL_AUDIT.md`.

## Result

At

$$
q=64,\qquad N=4096,\qquad M=16384,
$$

the noncubic class has:

- 16 entries with one middle degree-seven block and three singletons;
- 16 entries with one middle degree-nine block and three singletons; and
- 108 entries with two degree-five blocks and two singletons.

The degree-seven entries receive the new coefficients

$$
0.0170899096903
\quad\text{or}\quad
0.0382305883153.
$$

The other 124 noncubic entries receive the previously proved universal arbitrary-diagonal cross-Gram coefficient one. With these 140 coefficients inserted, the routing diagnostic improves from

$$
0.323811563171336
\quad\text{to}\quad
0.286902076794188.
$$

That gain makes the universal coefficient-one theorem affordable on the 96-entry two-split-cubic/one-split-higher class. After inserting it, the live diagnostic is

$$
U_{\mathrm{route}}=0.328938230122941,
$$

with raw margin

$$
{1\over3}-U_{\mathrm{route}}
=0.004395103210392.
$$

Retaining the declared $10^{-3}$ numerical allowance leaves conditional adaptive additive cap

$$
0.003395103210392,
$$

or multiplicative cap

$$
1.01032140049250.
$$

These caps remain requirements conditional on the 224 open entries retaining their frozen targets.

## Endpoint injective-character formula

Let $S$ be a degree-seven support with one odd row and one odd column. The two singleton endpoint links force those record-one sectors. After matching the odd label, the signed-permutation endpoint moment has the form

$$
M_{1,7}={1\over q}\,\omega\,J_x(\alpha_1,\ldots,\alpha_t),
$$

where $|\omega|=1$, the $\alpha_i$ are the nonzero XOR labels of the even support groups, and

$$
J_x(\alpha_1,\ldots,\alpha_t)
=
{1\over(q-1)_t}
\sum_{\substack{y_1,\ldots,y_t\ne x\\\text{distinct}}}
\prod_{i=1}^t\chi_{\alpha_i}(y_i).
$$

The reversed endpoint has the same formula for the active even groups on the other axis.

For $t\le3$, injective inclusion-exclusion gives

$$
|J_x(\alpha)|={1\over q-1},
$$

$$
|J_x(\alpha,\beta)|\le {1\over q-1},
$$

and

$$
|J_x(\alpha,\beta,\gamma)|
\le {3\over(q-1)(q-3)}.
$$

For the third line, expand the injective sum over set partitions of three indices. A nonzero Walsh character sums to a phase of magnitude one on $\mathbb F_q\setminus\{x\}$. If all three labels agree, the numerator has magnitude $3(q-2)$ and equality holds. If exactly two agree, if three distinct labels XOR to zero, or if their total XOR is nonzero, the numerator is no larger. The regression enumerates every label tuple and excluded point at $q=4,8$.

## Why the two endpoint residuals cannot both be large

Call an even row active when the XOR of its neighboring column labels is nonzero, and define active columns symmetrically. Let their counts be $t$ and $u$.

The only dangerous case would have no active group on one axis and fewer than three on the other. This cannot occur for a simple seven-edge bipartite support with one odd vertex on each side.

Suppose $t=0$, and let the odd row have degree $\delta$.

- If $\delta=7$, all seven columns are odd.
- If $\delta=5$, the remaining degree-two row is necessarily active.
- If $\delta=3$, the remaining four edges form one degree-four even row with zero neighbor XOR. Having exactly one odd column forces the odd row's three neighbors to lie inside that four-set. Those three shared columns have degree two and are all active, so $u=3$.
- If $\delta=1$, the remaining six edges either include an active degree-two row or form one degree-six row; the latter gives five or seven odd columns.

Thus $t=0$ forces $u=3$. By symmetry, $u=0$ forces $t=3$. If both are nonzero, the first two injective bounds give product at most $1/(q-1)^2$. Consequently every compatible degree-seven support satisfies

$$
|J_{\mathrm{row}}J_{\mathrm{column}}|
\le
A_q,
\qquad
A_q:={3\over(q-1)(q-3)}.
$$

At $q=64$,

$$
A_{64}=0.000780640124902420.
$$

The test enumerates all $2{,}928$ bidegree-$(1,1)$ seven-edge supports at $q=4$. Their active-group pairs are exactly

$$
(0,3),(1,3),(2,2),(2,3),
(3,0),(3,1),(3,2),(3,3).
$$

Direct enumeration of every signed permutation at $q=4$ also attains

$$
\max_S\max_{b,d}
|M_{1,7}(b,S)M_{7,1}(S,d)|
={A_4\over4^2}
={1\over16}.
$$

## Arbitrary-law occurrence coefficient

The full four-block chain has one additional singleton Hadamard link of magnitude $1/q$. Therefore the degree-seven endpoint-product lemma improves the universal record-one entry cap from $q^{-3}$ to

$$
{A_q\over q^3}.
$$

For each occurrence split, use the inherited arbitrary-law rank/incidence theorem. A singleton block has incidence $q^2$ or $1$, and the middle degree-seven block is constrained to record one on both axes. Its incidence is safely bounded by the smaller of the two exact one-axis record-one extension counts. Taking the minimum of cut rank, row incidence, and column incidence gives the inherited geometric coefficients

$$
21.8921743133
\quad\text{or}\quad
48.9733836319
$$

before the joint endpoint improvement. Multiplication by $A_{64}$ gives the two displayed theorem coefficients.

The proof is for arbitrary correlated row and column diagonal laws. It does not optimize the two endpoint links under separate laws.

## Universal recovery

For any fixed occurrence split, the complete physical moment kernel is a cross Gram matrix of unit-modulus features. Arbitrary diagonal weighting therefore has normalized nuclear coefficient at most one. This is the universal theorem already used for the q64 septimic and multicubic insertions.

It applies immediately to:

- the 16 middle degree-nine and 108 double-quintic noncubic entries; and
- the 96 two-split-cubic/one-split-higher entries.

Coefficient one was previously too expensive on the latter class. The degree-seven improvement creates enough routing margin to insert all 96 while preserving the declared reserve.

## Remaining gates

Two classes remain:

| class | entries | frozen target | common reserve gate |
|---|---:|---:|---:|
| higher-split-only cubic profile | 176 | $0.124035215254$ | $0.140343030565$ |
| one split cubic, no split higher block | 48 | $0.124035215254$ | $0.344887217413$ |

The 176-entry class is the live tight gate. The 48-entry class is also a possible margin-recovery project because its gate is much looser.

## Reproduction

Run:

    /opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_noncubic_recovered_universal_insertion.py --write-artifact
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_noncubic_recovered_universal_insertion.py

The regression checks the exact class inventory, character-average maxima, exhaustive degree-seven support geometry, direct signed-permutation endpoint moments, coefficient map, both Perron insertions, remaining class gates, and byte-for-byte artifact regeneration.
