# Adjacent balanced whole-row contraction

Date: 2026-07-15

Status: rigorous arbitrary-diagonal bound for the sixth balanced orbit. At
\(N=1024\), the coefficient is

$$
0.0422410016249<0.0570749885142,
$$

so it passes the live route gate. With all six chain-aware orbit theorems
inserted, the coarse completion diagnostic is \(0.332775779206<1/3\).
This is still not a complete one-batch or adaptive passive lower bound.

## 1. Exact target

For profile

$$
(1,1,3,5)
$$

and split \((0,1,1,3)\), write the first singleton labels as \(a,b\), the
cubic support as \(C=\{x\}\mathbin{\dot\cup}E\), with \(|E|=2\), and the
quintic support as \(S=F\mathbin{\dot\cup}G\), with \(|F|=3\) and
\(|G|=2\). Rows are \((b,x,F)\), columns are \((a,E,G)\), and the exact occurrence
kernel is

$$
K_{(b,x,F),(a,E,G)}
=H_N(b,a)M_{1,3}(b,C)M_{3,5}(C,S),
\tag{1.1}
$$

with the within-block disjointness conditions \(x\notin E\) and
\(F\cap G=\varnothing\).

The complement/reversal orbit is

$$
\begin{aligned}
(1,1,3,5)&:(0,1,1,3),\\
(1,1,3,5)&:(1,0,2,2),\\
(5,3,1,1)&:(3,1,1,0),\\
(5,3,1,1)&:(2,2,0,1).
\end{aligned}
\tag{1.2}
$$

The older \((0,1,2,2)\) adjacent slice cannot be reused directly: in the
present cut the first Hadamard link crosses the cut, and \(b\) shares a row
with \((x,F)\).

## 2. Reduce the theorem to one chain row energy

Define

$$
L_{(b,x,F),(E,G)}
=M_{1,3}(b,\{x\}\cup E)
 M_{3,5}(\{x\}\cup E,F\cup G).
\tag{2.1}
$$

For each row \(r=(b,x,F)\), use the vector

$$
u_r=(L_{r,(E,G)})_{E,G},
$$

and for a column \(c=(a,E,G)\), use the standard basis vector
\(v_c=e_{(E,G)}\). Thus \(L_{r,(E,G)}=\langle u_r,v_c\rangle\), so Schur
multiplication by \(L\) has trace-class norm at most

$$
\sqrt{R_q},\qquad
R_q=\max_{b,x,F}\sum_{E,G}|L_{(b,x,F),(E,G)}|^2.
\tag{2.2}
$$

After removing \(L\), the remaining matrix is the repeated Hadamard matrix
\(A_{(b,x,F),(a,E,G)}=H_N(b,a)\). For arbitrary nonnegative row and column
laws \(p,r\), duplicate compression gives masses

$$
P_b=\sum_{x,F}p_{b,x,F},\qquad
Q_a=\sum_{E,G}r_{a,E,G}.
$$

Because \(H_N\) is unitary, Schatten Hölder gives

$$
\left\|D_P^{1/2}H_ND_Q^{1/2}\right\|_1
\le
\left\|D_P^{1/2}\right\|_2
\left\|H_ND_Q^{1/2}\right\|_2
=\sqrt{\left(\sum P\right)\left(\sum Q\right)}.
\tag{2.3}
$$

Consequently the arbitrary-law coefficient of (1.1) is at most
\(\sqrt{R_q}\). It remains to bound the full \(M_{13}M_{35}\) row, including
both record sectors.

## 3. Record-one part

If the cubic has record one on both axes, it is an L-shape.  There are

$$
3(q-1)^2
\tag{3.1}
$$

such cubics through a fixed cell \(x\). Their preceding singleton moment
has constant squared modulus

$$
|M_{1,3}|^2={1\over q^2(q-1)^2}.
\tag{3.2}
$$

For a fixed triple \(F\), the number of degree-five row-record-one supports
containing \(F\) is maximized when the row multiplicities of \(F\) are
\(2+1\). Directly placing the two new cells gives

$$
D_{5,3}^{(1)}
={q-2\over2}(3q^2+q-6).
\tag{3.3}
$$

For comparison, the \(3\) and \(1+1+1\) row patterns give

$$
\binom{q-3}{2}+(q-1)\binom q2+q(q-1)(q-3)
$$

and \(3(q-1)^2\), respectively. Formula (3.3) is the largest for \(q\ge4\).

The universal record-one \(M_{35}\) entry bound is

$$
|M_{3,5}|
\le {q+2\over q(q-1)(q-2)}.
\tag{3.4}
$$

Multiplying (3.1)--(3.4) gives

$$
\boxed{
R_q^{(1)}
\le
{3D_{5,3}^{(1)}(q+2)^2
 \over q^4(q-1)^2(q-2)^2}.}
\tag{3.5}
$$

At \(q=32\),

$$
D_{5,3}^{(1)}=46470,
\qquad
R_{32}^{(1)}\le
{447661\over2519203840}
=0.000177699395695.
\tag{3.6}
$$

## 4. Record-three part

### 4.1 Cubic weight

A record-three cubic through \(x\) has row record one and three distinct
column labels.  There are

$$
D_*={(q-1)(q-2)(3q-2)\over2}
\tag{4.1}
$$

such cubics.  Of these,

$$
A=\binom{q-1}{2}
\tag{4.2}
$$

lie in one row and have \(|M_{13}|^2=1/q^2\). The other \(D_*-A\) have row
pattern \(2+1\) and squared moment \(1/[q^2(q-1)^2]\). Hence

$$
W_q^{(3)}
=\sum_C|M_{1,3}(b,C)|^2
={A\over q^2}+{D_*-A\over q^2(q-1)^2}.
\tag{4.3}
$$

This value is independent of \(b\): changing \(b\) changes only Walsh
phases. At \(q=32\),

$$
W_{32}^{(3)}={255\over512}=0.498046875.
\tag{4.4}
$$

### 4.2 Three fixed cells in one row

Suppose \(F\) lies in one row. A compatible record-three quintic adds one
cell in each of two distinct other rows. For a fixed record-three cubic,
expand \(|M_{35}|^2\) as a sum over two permutations of its three odd columns.
Summing the two new column labels over all of \(\mathbb F_q\) kills every
cross term between distinct permutations.  The remaining diagonal terms,
summed over the unordered pair of new rows, give the exact tail

$$
T_3={3\over(q-1)(q-2)}.
\tag{4.5}
$$

This identity holds for every record-three cubic, irrespective of its row
xors.  It is the decisive orthogonality gain over maximum-entry incidence.

### 4.3 The other row patterns

If the row pattern of \(F\) is \(2+1\), the compatible quintic extensions
split into

$$
n_0=q(q-2)^2,
\qquad
n_1=q^2\binom{q-2}{2},
\tag{4.6}
$$

with zero or one even row group. If the pattern is \(1+1+1\), the counts
are

$$
\begin{aligned}
n_0&=3\binom{q-1}{2},\\
n_1&=(q-3)\binom q2+3q(q-1)(q-3).
\end{aligned}
\tag{4.7}
$$

With no even row, \(|M_{35}|\le1/\binom q3\). With one even row, the
conditional character average has magnitude at most \(3/(q-3)\). Therefore

$$
T\le {n_0\over\binom q3^2}
 +{9n_1\over(q-3)^2\binom q3^2}.
\tag{4.8}
$$

At \(q=32\), the three tail bounds for row patterns \(3,2+1,1+1+1\) are

$$
0.00322580645161,\qquad
0.00136441924719,\qquad
0.000100502294216.
\tag{4.9}
$$

Thus (4.5) is the largest at the target order, and

$$
\boxed{
R_{32}^{(3)}
\le W_{32}^{(3)}T_3
={51\over31744}
=0.00160660282258.}
\tag{4.10}
$$

## 5. Arbitrary-law theorem at \(N=1024\)

The two record sectors occupy disjoint entries, so their squared row
energies add.  Equations (3.6) and (4.10) give

$$
R_{32}
\le R_{32}^{(1)}+R_{32}^{(3)}
={4495021\over2519203840}
=0.00178430221828.
\tag{5.1}
$$

Combining (2.2), (2.3), and (5.1), for arbitrary nonnegative diagonal laws,

$$
\boxed{
\|D_p^{1/2}KD_r^{1/2}\|_1
\le0.0422410016249
\sqrt{\left(\sum p\right)\left(\sum r\right)}.}
\tag{5.2}
$$

In particular,

$$
\boxed{
0.0422410016249
<0.0570749885142
=1.82639963245/q.}
\tag{5.3}
$$

Transpose and path reversal give (5.2) on all four cuts in (1.2).

## 6. Exact stress tests

Complete \(q=4\) enumeration evaluates all \(8960\) translated row types.
The largest exact energy is

$$
0.376302083333334,
\tag{6.1}
$$

attained at \(b=x=0\) and a vertical fixed triple. At \(q=8\), the exact
horizontal row has

$$
R_8^{(1)}=0.00931435632289,
\qquad
R_8^{(3)}=0.0334821428571.
\tag{6.2}
$$

The record-three value in (6.2) exactly matches (4.3)--(4.5). Exhaustive
extension counts at \(q=4,8\) reproduce (3.3), (4.6), and (4.7), and every
record-three cubic at those orders reproduces the horizontal identity
(4.5).  Random correlated diagonal laws on sparse submatrices of the exact
\(q=4\) occurrence tensor remain below the exact row-energy coefficient.

These computations protect the indexing, support restrictions, and
normalization. The \(q=32\) theorem follows from the explicit counts,
entry bounds, Walsh orthogonality, and Schur factorization above.

## 7. Ledger consequence

Insert (5.2) with the preceding five chain-aware orbit theorems.  Reoptimizing
the exact 210-state diagnostic and extended promise loss gives

$$
\beta=0.779325491007,
\qquad
\operatorname{TV}_{\rm diagnostic}=0.332775779206.
\tag{7.1}
$$

Thus

$$
\boxed{
{1\over3}-\operatorname{TV}_{\rm diagnostic}
=0.000557554127147.}
\tag{7.2}
$$

Twenty-four of the 888 balanced entries now have chain-aware arbitrary-law
coefficients.  The reranked next provisional orbit is generated by

$$
(3,1,1,5):(0,1,1,3),
\tag{7.3}
$$

with Perron contribution \(0.00100671864913\). Varying only that four-cut
orbit gives the next safe-coefficient gate

$$
\boxed{c_{\rm next}<0.0484819899411.}
\tag{7.4}
$$

## 8. Reproduction and scope

- `searches/adjacent_balanced_row_slice_contraction.py` computes the analytic
  coefficient, exact \(q=4\) table, selected \(q=8\) rows, six-theorem
  ledger, reranking, and next gate.
- `tests/adjacent_balanced_row_slice_contraction.py` checks all support counts,
  the horizontal Walsh identity, exact finite rows, target constants, and
  sparse arbitrary-law occurrence tensors.
- `./run_round3_checks.sh` runs the complete inherited and Round 3 suite.

The theorem does not certify the remaining 864 balanced entries, replace
the physical lower witnesses by arbitrary-law upper bounds, close the
one-batch ledger, or prove the adaptive lift.  It establishes a reusable
whole-row mechanism: retain the full adjacent chain until record
orthogonality and fixed-support incidence can be combined before applying
the arbitrary-law Schur bound.
