# Exact obstruction to two complete active flags

Date: 2026-07-15

Status: rigorous obstruction to the simplest lower-than-six modification of
the accepted active protocol.  At \(N=1024\), two complete folded-chain flags
cannot reach worst-case error \(1/3\), even with an arbitrary collective
POVM on both path and mode registers.  This is not a general active
hard-dose-five lower bound: a genuinely interleaved five-traversal circuit,
or two flags coherently combined with an extra one-traversal query, remains
open.

## 1. The reuse family being tested

For an input \(x\), the accepted active construction prepares

$$
|\Psi_x\rangle
={|0\rangle|L_x\rangle+|1\rangle|R_x\rangle\over\sqrt2},
\qquad
\langle L_x|R_x\rangle=F_{4,H}(x).
\tag{1.1}
$$

One complete flag has deterministic hard dose two.  The published
six-dose protocol prepares three independent copies, measures path \(X\),
and takes a majority vote.

The first reuse question is stronger than deleting the third Bernoulli
sample: prepare only two complete copies of (1.1), retain every path and mode
register, apply arbitrary sample-free joint quantum processing, and perform
an arbitrary collective two-outcome POVM.  Any classical reprocessing,
entangling decoder, delayed path measurement, or joint measurement of the
two completed flags is a special case.

A hard-dose-five protocol that still obtains information only by preparing
complete copies has at most two copies.  The unused fifth traversal cannot
create a third complete flag.  Section 6 states precisely which uses of that
traversal are outside the present family.

## 2. Exact endpoint ensemble at \(N=16\)

Start at \(N=4\).  Complete sign-cube enumeration gives 22,528 inputs with

$$
F_{4,H_4}=+{1\over2}
$$

and the same number with value \(-1/2\).  Let \(\mathcal A_+\) and
\(\mathcal A_-\) be the uniform distributions on those two sets.

For two independent factor inputs \(x,y\), define their tensor-product input
blockwise:

$$
(x\otimes y)^{(a)}=x^{(a)}\otimes y^{(a)}.
\tag{2.1}
$$

Because \(H_{16}=H_4\otimes H_4\),

$$
F_{4,H_{16}}(x\otimes y)
=F_{4,H_4}(x)F_{4,H_4}(y),
\tag{2.2}
$$

and the folded states tensor:

$$
L_{x\otimes y}=L_x\otimes L_y,
\qquad
R_{x\otimes y}=R_x\otimes R_y.
\tag{2.3}
$$

Define \(\mathcal E_+\) by choosing a uniform sign \(s\), then drawing
\(x\sim\mathcal A_s\) and \(y\sim\mathcal A_s\).  Define
\(\mathcal E_-\) by drawing \(x\sim\mathcal A_s\) and
\(y\sim\mathcal A_{-s}\).  Every input is exactly on the promise endpoint:

$$
F_{4,H_{16}}=+{1\over4}
\quad\hbox{under }\mathcal E_+,
\qquad
F_{4,H_{16}}=-{1\over4}
\quad\hbox{under }\mathcal E_-.
\tag{2.4}
$$

For \(k\) complete flags, let

$$
\rho_\pm^{(k)}
=\mathbb E_{z\sim\mathcal E_\pm}
\left(|\Psi_z\rangle\langle\Psi_z|\right)^{\otimes k}.
\tag{2.5}
$$

The optimal equal-prior average error of every collective POVM is

$$
P_{\rm Hel}^{(k)}
={1-\operatorname{TD}(\rho_+^{(k)},\rho_-^{(k)})\over2}.
\tag{2.6}
$$

Thus one exact average-state calculation gives a worst-case obstruction by
Yao's principle.

## 3. Two-copy moment reduction

For one \(N=4\) factor, write \(V_0=L\), \(V_1=R\).  For path word
\(p=(p_1,p_2)\in\{0,1\}^2\), define the exact conditional moment block

$$
M_s[p,q]
=\mathbb E_{x\sim\mathcal A_s}
\left(
V_{p_1}\otimes V_{p_2}
\right)
\left(
V_{q_1}\otimes V_{q_2}
\right)^{\mathsf T}.
\tag{3.1}
$$

The four path words are ordered \(00,01,10,11\).  Tensor factorization and
the sign-product mixture give every block of
\(\Delta=\rho_+^{(2)}-\rho_-^{(2)}\) exactly:

$$
\Delta[p,q]
={1\over8}\sum_{s=\pm1}
\left(
M_s[p,q]\otimes M_s[p,q]
-M_s[p,q]\otimes M_{-s}[p,q]
\right).
\tag{3.2}
$$

All same-parity path blocks vanish.  In the even/odd decomposition,

$$
\Delta=
\begin{pmatrix}
0&B\\
B^{\mathsf T}&0
\end{pmatrix},
\qquad
\operatorname{TD}(\rho_+^{(2)},\rho_-^{(2)})
=\|B\|_1.
\tag{3.3}
$$

The folded \(N=4\) state entries lie in
\(\{0,\pm1/2,\pm1\}\).  Consequently all of (3.1)--(3.3) can be formed by
integer counting.  The checker proves that

$$
\widehat B=247808\,B
\tag{3.4}
$$

is integral.

## 4. Exact Gram certificate

The nonzero graph of \(\widehat B\widehat B^{\mathsf T}\) has sixteen
connected components of size 32.  Exact characteristic-polynomial
factorization on those integer blocks gives the squared singular-value
spectrum

| squared singular value of \(\widehat B\) | multiplicity |
|---:|---:|
| \(14992384\) | 1 |
| \(2478080\) | 12 |
| \(1115136\) | 6 |
| \(671744\) | 18 |
| \(184320\) | 36 |
| \(147456\) | 18 |
| \(82944\) | 9 |
| \(0\) | 412 |

The rank is 100.  Summing the exact singular values and dividing by
247808 gives

$$
\boxed{
\operatorname{TD}(\rho_+^{(2)},\rho_-^{(2)})
={7\over88}
+{15\sqrt5\over242}
+{9\sqrt{41}\over968}
=0.277677889242648
<{1\over3}.}
\tag{4.1}
$$

Therefore

$$
\boxed{
P_{\rm Hel}^{(2)}
=0.361161055378676
>{1\over3}}
\tag{4.2}
$$

by \(0.0278277220453\).  If a two-complete-flag measurement had worst-case
error at most \(1/3\), its average error on \(\mathcal E_\pm\) would also be
at most \(1/3\), contradicting (4.2).

This obstruction includes the full mode registers.  The familiar calculation
with two classical path-\(X\) outcomes alone gives error \(3/8\), but that
weaker calculation is not used.

## 5. Isometric lift to \(N=1024\)

The explicit \(N=4\) input

$$
\begin{aligned}
z^{(1)}&=(-1,-1,-1,-1),\\
z^{(2)}&=(-1,-1,-1,-1),\\
z^{(3)}&=(-1,-1,-1,+1),\\
z^{(4)}&=(-1,-1,-1,+1)
\end{aligned}
\tag{5.1}
$$

has \(F_{4,H_4}(z)=1\), hence \(L_z=R_z\).  Tensoring three copies gives a
fixed \(N=64\) unit instance.  Tensor every \(N=16\) endpoint input from
Section 2 with this fixed factor.  Since \(16\cdot64=1024\), the result is an
exact \(N=1024\) endpoint ensemble.

On the complete flag state, this lift is the isometry

$$
|0\rangle|v\rangle\mapsto
|0\rangle|v\rangle|L_z\rangle^{\otimes3},
\qquad
|1\rangle|v\rangle\mapsto
|1\rangle|v\rangle|R_z\rangle^{\otimes3}.
\tag{5.2}
$$

Because \(L_z=R_z\), it is visibly a common tensor factor; even without that
equality, path orthogonality makes the controlled map isometric.  Trace
distance is preserved.  Equations (4.1)--(4.2) therefore hold unchanged at
\(N=1024\).

## 6. Decision and remaining active frontier

The following modifications of the accepted protocol are now ruled out at
the realistic target size:

1. delete the third complete flag and keep the original measurements;
2. replace majority vote by any classical rule on two complete flags;
3. delay all measurements and perform an arbitrary collective POVM on both
   complete path--mode states; or
4. add arbitrary sample-free coherent processing before that POVM.

Thus hard dose four, and hard dose five spent only on complete-flag
preparations plus sample-free decoding, cannot reach error \(1/3\).

The result does **not** cover:

- an extra one-traversal query whose output is kept coherently correlated
  with the two flags;
- a circuit that shares a charged traversal by ceasing to factor into
  complete flags;
- a new interleaved five-query algorithm; or
- a general active hard-dose-five lower bound.

The next bounded Track C question is therefore precise: can one additional
sample traversal, coherently integrated with two folds, overcome the exact
\(0.02783\) error deficit?  Any claimed active improvement must answer that
question with a complete five-traversal circuit, rather than remeasuring the
existing two flags.

## 7. Reproduction

- `searches/active_two_flag_collective_obstruction.py` enumerates the exact
  factor ensembles, constructs the integer moment matrix, factors all Gram
  blocks, and prints the Helstrom obstruction.
- `tests/active_two_flag_collective_obstruction.py` protects the exact
  spectrum, radical expression, endpoint tensor product, folded-state
  multiplicativity, and \(N=1024\) lift.
