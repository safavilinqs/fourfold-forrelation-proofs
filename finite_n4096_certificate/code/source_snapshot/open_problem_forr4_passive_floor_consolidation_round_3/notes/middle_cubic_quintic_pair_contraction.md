# Middle cubic--quintic pair contraction

Date: 2026-07-15

Status: rigorous arbitrary-diagonal bound for the eighth balanced orbit. At
\(N=1024\), the coefficient is

$$
0.0285281522923<0.0426269309291,
$$

and is also below the provisional \(1/32\). With all eight chain-aware orbit
theorems inserted, the coarse completion diagnostic is
\(0.332877589131<1/3\). This is still not a complete one-batch or adaptive
passive lower bound.

## 1. Exact target

For profile

$$
(1,3,5,1)
$$

and split \((0,2,2,1)\), write

$$
C=E\mathbin{\dot\cup}\{x\},
\qquad
S=F\mathbin{\dot\cup}G,
\qquad
|E|=|F|=2,\quad |G|=3.
$$

Rows are \((E,F,d)\), columns are \((a,x,G)\), and the exact occurrence
kernel is

$$
K_{(E,F,d),(a,x,G)}
=M_{1,3}(a,C)M_{3,5}(C,S)M_{5,1}(S,d).
\tag{1.1}
$$

The complement/reversal orbit is

$$
\begin{aligned}
(1,3,5,1)&:(0,2,2,1),\\
(1,3,5,1)&:(1,1,3,0),\\
(1,5,3,1)&:(0,3,1,1),\\
(1,5,3,1)&:(1,2,2,0).
\end{aligned}
\tag{1.2}
$$

Both decorated middle blocks are split, but the selected side fixes two
cells in each. This makes the complete row amenable to a chained endpoint
slice bound.

## 2. Endpoint fixed-pair energies

Let \(E_2(q)\) be the exact cubic endpoint squared slice through a fixed
pair, with the adjacent singleton fixed. The accepted endpoint
classification gives

$$
E_2(q)
={q^2-2q+2\over q^2(q-1)}.
\tag{2.1}
$$

Because the first singleton \(a\) is summed in a row of (1.1),

$$
\sum_{\substack{a,x\\x\notin E}}
|M_{1,3}(a,E\cup\{x\})|^2
\le N E_2(q).
\tag{2.2}
$$

Let \(F_2(q)\) be the exact quintic endpoint squared slice through a fixed
pair, with the final singleton fixed. Put

$$
w_0={1\over q^2},
\quad
w_1={1\over q^2(q-1)^2},
\quad
w_2={4\over q^2(q-1)^2(q-2)^2}.
\tag{2.3}
$$

The maximizing pair lies in one hidden column, and direct extension
classification gives

$$
\begin{aligned}
F_2(q)={}&
\binom{q-2}{3}w_0\\
&+q(q-1)\left[
\left({q\over2}-1\right)w_0+
\left(\binom{q-2}{2}-{q\over2}+1\right)w_1
\right]\\
&+\left[(q-2)(q-1)\binom q2+(q-1)\binom q3\right]w_1\\
&+q(q-1)(q-2)\left[
{q\over2}w_1+
\left(\binom q2-{q\over2}\right)w_2
\right].
\end{aligned}
\tag{2.4}
$$

Therefore, for fixed \(F,d\),

$$
\sum_{\substack{G:\ |G|=3\\G\cap F=\varnothing}}
|M_{5,1}(F\cup G,d)|^2
\le F_2(q).
\tag{2.5}
$$

At \(q=32\),

$$
E_2(32)={481\over15872},
\qquad
F_2(32)={159457\over7936}.
\tag{2.6}
$$

## 3. Universal middle-link maximum

The \(M_{3,5}\) link has only record-one and record-three sectors compatible
with the endpoint singletons. The accepted exact moment formulas give

$$
m_1(q)
={q+2\over q(q-1)(q-2)}
\tag{3.1}
$$

in record one. In record three, summing the six odd-label permutations and
allowing the single even row group gives

$$
m_3(q)
={3\over(q-3)\binom q3}
={18\over q(q-1)(q-2)(q-3)}.
\tag{3.2}
$$

For \(q\ge8\),

$$
{m_1(q)\over m_3(q)}
={(q+2)(q-3)\over18}>1.
\tag{3.3}
$$

Thus at the target order every compatible middle entry satisfies

$$
|M_{3,5}(C,S)|
\le m_1(32)
={17\over14880}
=0.00114247311828.
\tag{3.4}
$$

## 4. Chained row-energy theorem

Fix a row \((E,F,d)\). Taking the pointwise middle maximum in (1.1) separates
the two exact endpoint sums:

$$
\begin{aligned}
\sum_{a,x,G}|K_{(E,F,d),(a,x,G)}|^2
&\le
m_1(q)^2
\left(\sum_{a,x}|M_{1,3}(a,C)|^2\right)
\left(\sum_G|M_{5,1}(S,d)|^2\right)\\
&\le
N E_2(q)F_2(q)m_1(q)^2.
\end{aligned}
\tag{4.1}
$$

Use the complete kernel row as a Schur feature and the standard basis on
the column side. The weighted all-ones base has trace norm equal to the
geometric mean of the row and column masses. Hence, for arbitrary
nonnegative diagonal laws,

$$
\boxed{
\left\|D_p^{1/2}KD_r^{1/2}\right\|_1
\le
\sqrt{N E_2(q)F_2(q)}\,m_1(q)
\sqrt{\left(\sum p\right)\left(\sum r\right)}.}
\tag{4.2}
$$

At \(q=32\), the squared coefficient is exactly

$$
R_{32}
={22165958113\over27235742515200}
=0.000813855473212,
\tag{4.3}
$$

so

$$
\boxed{
\sqrt{R_{32}}
=0.0285281522923
< {1\over32}
<0.0426269309291.}
\tag{4.4}
$$

Transpose and path reversal give (4.2) on every cut in (1.2).

## 5. Exact stress tests

Complete \(q=4\) enumeration constructs the full weighted
cubic--quintic moment table and every fixed-pair incidence row. The maximum
exact squared row energy is

$$
0.115397805213,
\tag{5.1}
$$

with a horizontal fixed cubic pair and a vertical fixed quintic pair. The
record-one and record-three maxima are respectively

$$
0.00617283950617,
\qquad
0.110939643347.
\tag{5.2}
$$

At \(q=8\), the corresponding symmetry representative has exact squared row
energy

$$
0.000289480284345.
\tag{5.3}
$$

The tests separately enumerate the two endpoint sums in (2.2) and (2.5),
check the record-sector middle maxima, and stress sparse submatrices of the
full \(q=4\) tensor under correlated diagonal laws. These checks protect the
cut indexing, orientation of both endpoint records, within-support
disjointness, and normalization.

## 6. Ledger consequence

Insert (4.4) with the preceding seven chain-aware orbit theorems.
Reoptimizing the exact 210-state diagnostic and extended promise loss gives

$$
\beta=0.779322677591,
\qquad
\operatorname{TV}_{\rm diagnostic}=0.332877589131.
\tag{6.1}
$$

Thus

$$
\boxed{
{1\over3}-\operatorname{TV}_{\rm diagnostic}
=0.000455744202313.}
\tag{6.2}
$$

Thirty-two of the 888 balanced entries now have chain-aware arbitrary-law
coefficients. The reranked next provisional orbit is

$$
\begin{aligned}
(1,1,3,5)&:(0,0,3,2),\\
(1,1,3,5)&:(1,1,0,3),\\
(5,3,1,1)&:(2,3,0,0),\\
(5,3,1,1)&:(3,0,1,1),
\end{aligned}
\tag{6.3}
$$

with Perron contribution \(0.000981156614081\). Varying only this four-cut
orbit gives the next safe-coefficient gate

$$
\boxed{c_{\rm next}<0.0454321892133.}
\tag{6.4}
$$

## 7. Reproduction and scope

- searches/middle_cubic_quintic_pair_contraction.py computes the endpoint
  slices, middle maximum, eighth coefficient, ledger, reranking, and next
  gate.
- tests/middle_cubic_quintic_pair_contraction.py checks complete \(q=4\)
  rows, the leading \(q=8\) row, endpoint factors, middle maxima, sparse
  occurrence tensors, ledger, and next gate.
- ./run_round3_checks.sh runs the complete inherited and Round 3 suite.

The theorem does not certify the remaining 856 balanced entries, replace the
physical lower witnesses by arbitrary-law upper bounds, close the one-batch
ledger, or prove the adaptive lift. It shows that a pointwise middle-link
bound can be effective when it is applied only after exact endpoint slice
energies have been retained on both sides.
