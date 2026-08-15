# Whole-cubic middle-pair contraction

Date: 2026-07-15

Status: rigorous arbitrary-diagonal bound for the ninth balanced orbit. At
\(N=1024\), the coefficient is

$$
0.0250919471547<0.0454321892133,
$$

and is also below the provisional \(1/32\). With all nine chain-aware orbit
theorems inserted, the coarse completion diagnostic is
\(0.332686212434<1/3\). This is still not a complete one-batch or adaptive
passive lower bound.

## 1. Exact target

For profile

$$
(1,1,3,5)
$$

and split \((0,0,3,2)\), let \(a,b\) be the singleton labels, let \(C\) be
the complete cubic support, and write

$$
S=F\mathbin{\dot\cup}G,
\qquad |C|=|G|=3,\quad |F|=2.
$$

Rows are \((C,F)\), columns are \((a,b,G)\), and the exact occurrence kernel
is

$$
K_{(C,F),(a,b,G)}
=H_N(a,b)M_{1,3}(b,C)M_{3,5}(C,F\cup G).
\tag{1.1}
$$

The complement/reversal orbit is

$$
\begin{aligned}
(1,1,3,5)&:(0,0,3,2),\\
(1,1,3,5)&:(1,1,0,3),\\
(5,3,1,1)&:(2,3,0,0),\\
(5,3,1,1)&:(3,0,1,1).
\end{aligned}
\tag{1.2}
$$

The useful feature is that the cubic is wholly fixed in a row. This permits
the complete \(M_{35}\) extension row to be extracted without discarding the
Walsh compression in the preceding \(H_NM_{13}\) chain.

## 2. Normalize the complete middle row

The endpoint link \(M_{13}(b,C)\) forces the parity record of \(C\) on that
side to have size one. Its record on the \(M_{35}\) side therefore has size
\(\rho=1\) or \(3\). The exact endpoint formula is

$$
M_{1,3}(b,C)
=v_\rho(C)H_N(\xi(C),b),
\qquad
\xi(C)=\bigoplus_{x\in C}x,
\tag{2.1}
$$

with

$$
|v_1(C)|={1\over q-1},
\qquad
|v_3(C)|=1.
\tag{2.2}
$$

Let

$$
S_\rho(q)=
\max_{\substack{C,F\\\operatorname{rec}(C)=\rho}}
\sum_{\substack{G:\ |G|=3\\G\cap F=\varnothing}}
|M_{3,5}(C,F\cup G)|^2.
\tag{2.3}
$$

For a row \(r=(C,F)\) in sector \(\rho\), use the Schur feature

$$
u_r=
{1\over\sqrt{S_\rho(q)}}
\bigl(M_{3,5}(C,F\cup G)\bigr)_G,
\tag{2.4}
$$

and for a column \((a,b,G)\), use \(e_G\). Every row feature has norm at most
one. After this contraction, the residual matrix is

$$
B_{(C,F),(a,b,G)}
=\sqrt{S_\rho(q)}H_N(a,b)M_{1,3}(b,C).
\tag{2.5}
$$

## 3. Preserve the extra Walsh factor

Let \(p_{C,F}\) and \(r_{a,b,G}\) be arbitrary nonnegative row and column
laws. Compress the duplicate labels:

$$
P_C=\sum_Fp_{C,F},
\qquad
R_{a,b}=\sum_Gr_{a,b,G}.
\tag{3.1}
$$

Using (2.1), collapse all cubic rows with the same XOR label into

$$
Z_x=
\sum_{\substack{C:\ \xi(C)=x}}
P_C S_{\rho(C)}(q)|v_{\rho(C)}(C)|^2.
\tag{3.2}
$$

For each \(b\), all columns indexed by \(a\) are proportional. Since every
normalized Hadamard entry has squared modulus \(1/N=1/q^2\), the second
duplicate compression gives

$$
W_b=\sum_aR_{a,b}|H_N(a,b)|^2
={1\over N}\sum_aR_{a,b}.
\tag{3.3}
$$

The remaining matrix is \(D_Z^{1/2}H_ND_W^{1/2}\). Schatten Hölder and
unitarity of \(H_N\) give

$$
\|D_Z^{1/2}H_ND_W^{1/2}\|_1
\le\sqrt{\left(\sum_xZ_x\right)\left(\sum_bW_b\right)}.
\tag{3.4}
$$

Consequently, the arbitrary-law coefficient of (1.1) is at most

$$
\boxed{
c(q)\le
\max\left\{
{\sqrt{S_1(q)}\over q(q-1)},
{\sqrt{S_3(q)}\over q}
\right\}.}
\tag{3.5}
$$

The factor \(1/q\) in (3.5) is decision-critical. Treating \(M_{13}M_{35}\)
as one row feature would spend the cubic Walsh structure and lose this gain.

## 4. Record-one slice

The universal record-one middle entry satisfies

$$
|M_{3,5}|
\le m_1(q)
={q+2\over q(q-1)(q-2)}.
\tag{4.1}
$$

For a fixed pair in one physical row, the number of record-one quintic
extensions is

$$
\begin{aligned}
D_1(q)={}&
\binom{q-2}{3}+(q-1)\binom q3
 +(q-2)(q-1)\binom q2\\
&+(q-1)q\left[
\binom{q-2}{2}+(q-2)\binom q2
\right].
\end{aligned}
\tag{4.2}
$$

This classifies the three added cells by row multiplicity \(3\) or \(2+1\)
and requires exactly one final odd row. If the fixed pair lies in distinct
rows, the corresponding count is

$$
2\binom{q-1}{3}
+2(q-1)\left[\binom{q-1}{2}+(q-2)\binom q2\right]
+(q-1)^2(q-2)q.
\tag{4.3}
$$

The difference between (4.2) and (4.3) is

$$
{q(q-2)(q+1)(3q^2-14q+14)\over6}>0
\tag{4.4}
$$

for \(q\ge4\). Thus

$$
S_1(q)\le D_1(q)m_1(q)^2.
\tag{4.5}
$$

At \(q=32\),

$$
D_1=15{,}811{,}580,
\qquad
S_1\le {228477331\over11070720}
=20.6379829857,
\tag{4.6}
$$

and the record-one contribution to (3.5) is only

$$
c_1(32)\le0.00457954101408.
\tag{4.7}
$$

## 5. Record-three slice

A compatible record-three cubic has three odd column labels. A quintic with
three odd row labels has one of two row-multiplicity patterns:

- \(3+1+1\), with no even row; or
- \(2+1+1+1\), with one even row.

The exact record-three moment formula gives the pointwise bounds

$$
m_{3,0}(q)={1\over\binom q3},
\qquad
m_{3,1}(q)={3\over(q-3)\binom q3},
\tag{5.1}
$$

respectively. The second factor is the injective average of the one even-row
Walsh label over the \(q-3\) unused columns.

Classifying the three added cells gives the exact extension counts

| fixed-pair rows | \(3+1+1\) | \(2+1+1+1\) |
|---|---:|---:|
| same row | \(\binom{q-1}{2}(q-2)q^2\) | \(\binom{q-1}{3}q^3\) |
| distinct rows | \((q-2)\binom q3+2(q-2)q\binom{q-1}{2}\) | \((q-2)q(q-3)\binom q2+2(q-1)\binom{q-2}{2}q^2\) |

Hence

$$
\begin{aligned}
S_{3,\mathrm{same}}
&\le
n_{\mathrm{same},0}m_{3,0}^2
+n_{\mathrm{same},1}m_{3,1}^2,\\
S_{3,\mathrm{distinct}}
&\le
n_{\mathrm{distinct},0}m_{3,0}^2
+n_{\mathrm{distinct},1}m_{3,1}^2.
\end{aligned}
\tag{5.2}
$$

At \(q=32\), the four counts are

$$
14{,}284{,}800,\quad147{,}292{,}160,\quad
1{,}041{,}600,\quad41{,}425{,}920,
\tag{5.3}
$$

and therefore

$$
S_{3,\mathrm{same}}
\le {2898\over4495}
=0.644716351502,
\qquad
S_{3,\mathrm{distinct}}
\le {4341\over71920}
=0.0603587319244.
\tag{5.4}
$$

The same-row class dominates. Substitution in (3.5) gives

$$
\boxed{
c(32)=c_3(32)
\le\sqrt{{1449\over2301440}}
=0.0250919471547
<{1\over32}
<0.0454321892133.}
\tag{5.5}
$$

Transpose and path reversal prove (5.5) on every cut in (1.2).

## 6. Exact stress tests

Complete \(q=4\) construction of every compatible cubic and every fixed-pair
\(M_{35}\) row gives maximum exact squared slices

$$
2.75\quad\text{in record one},
\qquad
2.50\quad\text{in record three}.
\tag{6.1}
$$

After multiplying by the exact cubic endpoint energy, these become
\(0.305555555556\) and \(2.50\). At \(q=8\), the leading checked
representatives give

$$
S_1^{\rm exact}=1.17015306122,
\quad
S_{3,\rm same}^{\rm exact}=0.466666666667,
\quad
S_{3,\rm distinct}^{\rm exact}=0.25.
\tag{6.2}
$$

The regression enumerates both pair geometries for every \(q=4,\ldots,8\),
verifies the four counts in Section 5, checks the residual Walsh compression
under correlated diagonal laws, and stresses sparse exact \(q=4\) target
tensors. These checks protect record orientation, pair disjointness, sector
normalization, and the extra \(1/q\).

## 7. Ledger consequence

Insert (5.5) with the preceding eight chain-aware orbit theorems. Reoptimizing
the exact 210-state diagnostic and extended promise loss gives

$$
\beta=0.779331115015,
\qquad
\operatorname{TV}_{\rm diagnostic}=0.332686212434.
\tag{7.1}
$$

Thus

$$
\boxed{
{1\over3}-\operatorname{TV}_{\rm diagnostic}
=0.000647120899295.}
\tag{7.2}
$$

Thirty-six of the 888 balanced entries now have chain-aware arbitrary-law
coefficients. The reranked next provisional orbit is

$$
\begin{aligned}
(1,3,5,1)&:(0,2,3,0),\\
(1,3,5,1)&:(1,1,2,1),\\
(1,5,3,1)&:(0,3,2,0),\\
(1,5,3,1)&:(1,2,1,1),
\end{aligned}
\tag{7.3}
$$

with Perron contribution \(0.000924037577192\). Varying only this four-cut
orbit gives the next safe-coefficient gate

$$
\boxed{c_{\rm next}<0.0529177166680.}
\tag{7.4}
$$

## 8. Reproduction and scope

- searches/whole_cubic_middle_pair_contraction.py computes the pair counts,
  sector bounds, ninth coefficient, ledger, reranking, and next gate.
- tests/whole_cubic_middle_pair_contraction.py checks exact \(q=4\) slices,
  \(q=8\) representatives, all pair-count formulas, Walsh compression,
  sparse occurrence tensors, ledger, and next gate.
- ./run_round3_checks.sh runs the complete inherited and Round 3 suite.

The theorem does not certify the remaining 852 balanced entries, replace the
physical lower witnesses by arbitrary-law upper bounds, close the one-batch
ledger, or prove the adaptive lift. It does identify a reusable mechanism:
extract a normalized local row without spending the XOR-labelled endpoint
matrix that remains behind it.
