# Whole-cubic quintic-triple contraction

Date: 2026-07-15

Status: rigorous arbitrary-diagonal bound for the seventh balanced orbit. At
\(N=1024\), the coefficient is

$$
0.0370952793157<0.0484819899411,
$$

so it passes the live route gate. With all seven chain-aware orbit theorems
inserted, the coarse completion diagnostic is \(0.332964363589<1/3\).
This remains a target diagnostic, not a one-batch or adaptive passive lower
bound.

## 1. Exact target

For profile

$$
(3,1,1,5)
$$

and split \((0,1,1,3)\), let \(Q\) be the cubic support, \(b,c\) the two
singleton coordinates, and

$$
S=F\mathbin{\dot\cup}G,
\qquad |F|=3,\quad |G|=2.
$$

Rows are \((b,c,F)\), columns are \((Q,G)\), and the exact occurrence
kernel is

$$
K_{(b,c,F),(Q,G)}
=M_{3,1}(Q,b)H_N(b,c)M_{1,5}(c,F\cup G).
\tag{1.1}
$$

The complement/reversal orbit is

$$
\begin{aligned}
(3,1,1,5)&:(0,1,1,3),\\
(3,1,1,5)&:(3,0,0,2),\\
(5,1,1,3)&:(2,0,0,3),\\
(5,1,1,3)&:(3,1,1,0).
\end{aligned}
\tag{1.2}
$$

The useful asymmetry is that the quintic endpoint is split \(3|2\), whereas
the cubic endpoint is wholly on the column side. This allows the exact
quintic row energy to be extracted before compressing the complete cubic
endpoint.

## 2. Exact fixed-three quintic energy

For each row \(r=(b,c,F)\), define the Schur feature

$$
u_r=\bigl(M_{1,5}(c,F\cup G)\bigr)_G,
$$

and for each column \(d=(Q,G)\), use the standard basis vector \(v_d=e_G\).
The Schur multiplier norm is therefore at most

$$
\sqrt{F_3(q)},
\qquad
F_3(q)=
\max_{c,F}\sum_{\substack{G:\ |G|=2\\G\cap F=\varnothing}}
|M_{1,5}(c,F\cup G)|^2.
\tag{2.1}
$$

The accepted endpoint-quintic classification shows that the maximizing
fixed triple lies in one hidden column. Put

$$
w_0={1\over q^2},
\qquad
w_1={1\over q^2(q-1)^2}.
\tag{2.2}
$$

Directly classifying the two-cell extensions by column multiplicity and XOR
type gives

$$
\boxed{
F_3(q)=
\binom{q-3}{2}w_0
+q(q-1)\bigl[w_0+(q-4)w_1\bigr]
+(q-1)\binom q2w_1.}
\tag{2.3}
$$

At \(q=32\),

$$
F_3(32)
={22365\over15872}
=1.409085181451613.
\tag{2.4}
$$

## 3. Compress the remaining cubic--Hadamard matrix

After removing the quintic Schur feature, the remaining matrix is

$$
B_{(b,c,F),(Q,G)}
=H_N(b,c)M_{3,1}(Q,b).
\tag{3.1}
$$

Let \(p_{b,c,F}\) and \(r_{Q,G}\) be arbitrary nonnegative row and column
laws, and define their compressed masses

$$
P_b=\sum_{c,F}p_{b,c,F},
\qquad
R_Q=\sum_G r_{Q,G}.
\tag{3.2}
$$

For fixed \(b\), all rows indexed by \(c,F\) in (3.1) are proportional, and
\(|H_N(b,c)|=1/q\). For fixed \(Q\), all columns indexed by \(G\) are
duplicates. Row and column isometries therefore give

$$
\left\|D_p^{1/2}BD_r^{1/2}\right\|_1
={1\over q}
\left\|D_P^{1/2}M_{3,1}D_R^{1/2}\right\|_1.
\tag{3.3}
$$

The cubic endpoint moment has the exact XOR-labelled form

$$
M_{3,1}(Q,b)
=v_3(Q)H_N(\xi(Q),b),
\qquad
\xi(Q)=\bigoplus_{x\in Q}x,
\qquad
|v_3(Q)|\le1.
\tag{3.4}
$$

Collapse all proportional columns with the same XOR label by setting

$$
Z_x=\sum_{\substack{Q:\ \xi(Q)=x}}R_Q|v_3(Q)|^2.
\tag{3.5}
$$

This preserves the nonzero singular values and yields

$$
\left\|D_P^{1/2}M_{3,1}D_R^{1/2}\right\|_1
=
\left\|D_P^{1/2}H_N^TD_Z^{1/2}\right\|_1.
\tag{3.6}
$$

Since \(H_N\) is unitary, Schatten Hölder gives

$$
\begin{aligned}
\left\|D_P^{1/2}H_N^TD_Z^{1/2}\right\|_1
&\le
\left\|D_P^{1/2}\right\|_2
\left\|H_N^TD_Z^{1/2}\right\|_2\\
&=\sqrt{\left(\sum_bP_b\right)
         \left(\sum_xZ_x\right)}\\
&\le
\sqrt{\left(\sum p\right)\left(\sum r\right)}.
\end{aligned}
\tag{3.7}
$$

Thus the complete base contributes exactly the safe factor \(1/q\).

## 4. Arbitrary-law theorem at \(N=1024\)

Combining (2.1), (3.3), and (3.7) gives, for arbitrary nonnegative diagonal
laws,

$$
\boxed{
\left\|D_p^{1/2}KD_r^{1/2}\right\|_1
\le
{\sqrt{F_3(q)}\over q}
\sqrt{\left(\sum p\right)\left(\sum r\right)}.}
\tag{4.1}
$$

At \(q=32\),

$$
\boxed{
{\sqrt{F_3(32)}\over32}
=0.0370952793157
<0.0484819899411.}
\tag{4.2}
$$

Transpose and path reversal give the same bound on all four cuts in (1.2).
The coefficient is larger than the provisional \(1/32\), but the live gate
was computed after reoptimizing the complete ledger and is the relevant
acceptance test.

## 5. Exact stress tests

Complete \(q=4\) enumeration of every fixed triple and singleton gives

$$
F_3(4)={7\over8},
\tag{5.1}
$$

with a maximizing triple in one hidden column. The same direct endpoint
calculation at \(q=8\) gives

$$
F_3(8)={261\over224}.
\tag{5.2}
$$

The cubic compression identity (3.6) is checked directly at \(q=4\) under
random row and column laws. It also has a saturating law: choose one
unit-weight cubic for every XOR label and make both compressed laws uniform.
The resulting weighted cubic trace norm is exactly one.

Random correlated laws on sparse submatrices of the full \(q=4\) tensor
(1.1) stay below (4.1). These computations protect the occurrence indexing,
within-support disjointness, XOR compression, and all normalizations.

## 6. Ledger consequence

Insert (4.2) with the preceding six chain-aware orbit theorems. Reoptimizing
the exact 210-state diagnostic and extended promise loss gives

$$
\beta=0.779319212599,
\qquad
\operatorname{TV}_{\rm diagnostic}=0.332964363589.
\tag{6.1}
$$

Hence

$$
\boxed{
{1\over3}-\operatorname{TV}_{\rm diagnostic}
=0.000368969744373.}
\tag{6.2}
$$

Twenty-eight of the 888 balanced entries now have chain-aware arbitrary-law
coefficients. The reranked next provisional orbit is

$$
\begin{aligned}
(1,3,5,1)&:(0,2,2,1),\\
(1,3,5,1)&:(1,1,3,0),\\
(1,5,3,1)&:(0,3,1,1),\\
(1,5,3,1)&:(1,2,2,0),
\end{aligned}
\tag{6.3}
$$

with Perron contribution \(0.000999555962451\). Varying only this four-cut
orbit gives the next safe-coefficient gate

$$
\boxed{c_{\rm next}<0.0426269309291.}
\tag{6.4}
$$

## 7. Reproduction and scope

- searches/whole_cubic_quintic_triple_contraction.py computes the exact
  slice formula, target coefficient, seven-theorem ledger, reranking, and
  next gate.
- tests/whole_cubic_quintic_triple_contraction.py checks the endpoint slices,
  exact XOR compression and saturation, sparse occurrence tensors, ledger,
  and next gate.
- ./run_round3_checks.sh runs the complete inherited and Round 3 suite.

The theorem does not certify the remaining 860 balanced entries, replace the
physical lower witnesses by arbitrary-law upper bounds, close the one-batch
ledger, or prove the adaptive lift. It supplies another reusable mechanism:
extract a split high-degree endpoint as a Schur row, then compress the
opposite whole endpoint by its XOR-labelled Walsh columns before applying
Schatten Hölder.
