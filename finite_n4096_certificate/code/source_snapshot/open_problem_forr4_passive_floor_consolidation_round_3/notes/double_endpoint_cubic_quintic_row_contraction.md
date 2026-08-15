# Double-endpoint cubic--quintic row contraction

Date: 2026-07-15

Status: rigorous arbitrary-diagonal bound for the tenth balanced orbit. At
\(N=1024\), the coefficient is

$$
0.0462425962446<0.0529177166680.
$$

With all ten chain-aware orbit theorems inserted, the coarse completion
diagnostic is

$$
0.333132605485<1/3.
$$

This is still a diagnostic inside the repaired moment-ledger route, not a
complete one-batch or adaptive passive lower bound. The new coefficient is
also larger than the provisional \(1/32\); its value is that it replaces an
assumption by a proof while remaining below the live acceptance gate.

## 1. Exact target

For profile

$$
(1,3,5,1)
$$

and split \((0,2,3,0)\), write

$$
C=E\mathbin{\dot\cup}\{x\},
\qquad
S=F\mathbin{\dot\cup}G,
\qquad
|E|=|G|=2,\quad |F|=3.
$$

Rows are \((E,F)\), columns are \((a,x,G,d)\), and the exact occurrence
kernel is

$$
K_{(E,F),(a,x,G,d)}
=M_{1,3}(a,C)M_{3,5}(C,S)M_{5,1}(S,d).
\tag{1.1}
$$

The complement/reversal orbit is

$$
\begin{aligned}
(1,3,5,1)&:(0,2,3,0),\\
(1,3,5,1)&:(1,1,2,1),\\
(1,5,3,1)&:(0,3,2,0),\\
(1,5,3,1)&:(1,2,1,1).
\end{aligned}
\tag{1.2}
$$

The useful feature is that both endpoint links expose Walsh characters.
Keeping both characters, instead of absorbing either into a row norm,
leaves a tensor product of two unitary Hadamards.

## 2. Extract one scalar completion row

The endpoint formulas are

$$
\begin{aligned}
M_{1,3}(a,C)&=v_3(C)H_N(\xi(C),a),\\
M_{5,1}(S,d)&=v_5(S)H_N(\xi(S),d),
\end{aligned}
\qquad
\xi(T)=\bigoplus_{z\in T}z.
\tag{2.1}
$$

For a fixed row \(r=(E,F)\), define its scalar completion vector by

$$
L_r(x,G)
=v_3(E\cup\{x\})
 M_{3,5}(E\cup\{x\},F\cup G)
 v_5(F\cup G).
\tag{2.2}
$$

Let

$$
R(q)=\max_{E,F}\sum_{x,G}|L_{E,F}(x,G)|^2.
\tag{2.3}
$$

The Schur features \(L_r/\sqrt{R(q)}\) and the coordinate vectors
\(e_{x,G}\) have norm at most one. After extracting their inner product,
the residual matrix is

$$
B_{(E,F),(a,x,G,d)}
=\sqrt{R(q)}
H_N(\xi(E)\oplus x,a)
H_N(\xi(F)\oplus\xi(G),d).
\tag{2.4}
$$

## 3. The residual coefficient is one

Walsh multiplication gives

$$
\begin{aligned}
H_N(\xi(E)\oplus x,a)
 &=\chi(x,a)H_N(\xi(E),a),\\
H_N(\xi(F)\oplus\xi(G),d)
 &=\chi(\xi(G),d)H_N(\xi(F),d).
\end{aligned}
\tag{3.1}
$$

The two characters in (3.1) depend only on the column and therefore do not
change a weighted nuclear norm. For arbitrary row and column probability
laws, compress rows sharing \((\xi(E),\xi(F))\), and compress columns sharing
\((a,d)\). The remaining matrix is

$$
D_P^{1/2}(H_N\otimes H_N)D_Q^{1/2},
\qquad
\operatorname{tr}D_P=\operatorname{tr}D_Q=1.
\tag{3.2}
$$

Because \(H_N\otimes H_N\) is unitary, Schatten Hölder gives

$$
\left\|D_P^{1/2}(H_N\otimes H_N)D_Q^{1/2}\right\|_1
\le
\|D_P^{1/2}\|_2\|D_Q^{1/2}\|_2
=1.
\tag{3.3}
$$

Thus the arbitrary-diagonal coefficient of (1.1) is at most

$$
\boxed{c(q)\le\sqrt{R(q)}}.
\tag{3.4}
$$

## 4. Record-one sector

The middle moment is nonzero only when its two parity records have the same
odd size. In the record-one sector, a cubic through a fixed pair \(E\) has
at most

$$
D_C(q)=2(q-1)
\tag{4.1}
$$

compatible L-shape completions. Each has

$$
|v_3(C)|={1\over q-1}.
\tag{4.2}
$$

For a quintic through a fixed triple \(F\), the exact maximum number with
both records one is

$$
D_Q(q)=2(q-2)(2q-1),
\tag{4.3}
$$

and \(|v_5(S)|\le1\). The universal record-one middle bound is

$$
m_1(q)={q+2\over q(q-1)(q-2)}.
\tag{4.4}
$$

Consequently,

$$
R_1(q)
\le
{2\over q-1}\,
2(q-2)(2q-1)\,
\left({q+2\over q(q-1)(q-2)}\right)^2.
\tag{4.5}
$$

At \(q=32\),

$$
D_Q=3780,
\qquad
R_1={6069\over19066240}
=0.000318311318855.
\tag{4.6}
$$

## 5. Record-three sector

In the no-even-group record-three class, the exact permanent formula gives

$$
|M_{3,5}|\le{1\over\binom q3}.
\tag{5.1}
$$

The one-even-group correction is at most \(3/(q-3)\) times the same bound.
For \(q\ge8\), this factor is at most one, so (5.1) holds throughout the
record-three sector.

Let \(E_2(q)\) be the exact cubic endpoint squared slice through a fixed pair
and let \(F_3(q)\) be the exact quintic endpoint squared slice through a
fixed triple. Removing the fixed endpoint singleton from the normalized
Hadamard moment multiplies these energies by \(N=q^2\):

$$
\sum_x|v_3(E\cup\{x\})|^2\le N E_2(q),
\qquad
\sum_G|v_5(F\cup G)|^2\le N F_3(q).
\tag{5.2}
$$

Therefore

$$
R_3(q)
\le
{N^2E_2(q)F_3(q)\over\binom q3^2}.
\tag{5.3}
$$

At \(q=32\), the exact endpoint slices are

$$
E_2={481\over15872},
\qquad
F_3={22365\over15872},
\tag{5.4}
$$

and hence

$$
R_3={2151513\over1182106880}
=0.00182006638858.
\tag{5.5}
$$

The two record sectors occur in the same scalar row, so their squared
energies add rather than taking a maximum:

$$
\begin{aligned}
R(32)
&\le R_1+R_3\\
&=0.00213837770744,\\
c(32)
&\le\sqrt{R_1+R_3}
=0.0462425962446.
\end{aligned}
\tag{5.6}
$$

## 6. Independent checks

The regression tests perform four checks that are independent of the final
\(q=32\) arithmetic:

1. Exact signed-permutation enumeration at \(q=4\) constructs all scalar
   rows and finds maximum squared energy \(1.58693415638\).
2. Direct \(q=8\) enumeration of the horizontal-pair/vertical-triple
   representative gives record-one energy
   \(4.19680575\times10^{-6}\) and record-three energy
   \(0.00537318634\).
3. Random correlated laws stress both the residual
   \(H_N\otimes H_N\) compression and sparse exact target tensors.
4. Exact endpoint-slice formulas and the ten-theorem ledger reproduce every
   displayed \(q=32\) number.

The executable certificate is
searches/double_endpoint_cubic_quintic_row_contraction.py and its regression
is tests/double_endpoint_cubic_quintic_row_contraction.py.

## 7. Ledger consequence

The theorem controls four more profile-split entries, bringing the
chain-aware total to \(40\) of \(888\) relevant entries. At the reoptimized

$$
\beta=0.779315347679,
\tag{7.1}
$$

the current diagnostic is

$$
\mathcal D_{10}=0.333132605485,
\qquad
{1\over3}-\mathcal D_{10}
=0.000200727847845.
\tag{7.2}
$$

The next sensitivity-ranked unresolved orbit is led by

$$
(1,1,5,3):(0,1,3,1).
\tag{7.3}
$$

Its present contribution is \(0.000936859407035\), and its live admissible
coefficient is

$$
c_{\mathrm{next}}<0.0379251204234.
\tag{7.4}
$$

This narrow gate is the important planning signal: future work should not
continue mechanically through the ledger unless the next tensor exposes a
clear chain-aware contraction below (7.4). Parallel effort on a stronger
global contraction, a different hard instance, or a passive counterexample
remains part of the Round 3 mission.
