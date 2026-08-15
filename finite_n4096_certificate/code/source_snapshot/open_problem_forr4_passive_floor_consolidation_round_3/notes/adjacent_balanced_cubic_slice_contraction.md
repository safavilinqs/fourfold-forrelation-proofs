# Adjacent balanced cubic-slice contraction

Date: 2026-07-15

Status: rigorous arbitrary-diagonal fixed-split upper bound for the second
previously unforced dose-six-relevant orbit.  The proved coefficient
\(0.01627246928\) is below the provisional \(1/32\).  With both
balanced-orbit theorems inserted, the coarse completion ledger is
\(0.32556385797\),
still not a complete one-batch or adaptive certificate.

## 1. Exact tensor and cut

Take profile

$$
(1,1,3,5)
$$

and split \((0,0,1,4)\).  Write the first two singleton labels as \(a,b\),
the cubic support as \(\{x\}\mathbin{\dot\cup}E\), with \(|E|=2\), and the
quintic support as \(F\mathbin{\dot\cup}\{e\}\), with \(|F|=4\).  Rows are
\((x,F)\), columns are \((a,b,E,e)\), and the occurrence tensor is

$$
K[(x,F),(a,b,E,e)]
=H_N(a,b)
 M_{1,3}(b,\{x\}\cup E)
 M_{3,5}(\{x\}\cup E,F\cup\{e\}),
\tag{1.1}
$$

on \(x\notin E\) and \(e\notin F\), and is zero otherwise.  This is a new
reshuffling of the adjacent cubic--quintic tensor.  The earlier
\((0,1,2,2)\) squared-slice result cannot be inserted without reconstructing
the present row and column fibers.

All bounds below allow arbitrary nonnegative diagonal weights on these row
and column labels.

## 2. Why the generic two-mask theorem is insufficient

Completing both link moments gives a unit-feature Gram Schur multiplier.
Restoring the cubic \(1|2\) distinct-label mask and quintic \(4|1\) mask
separately would cost

$$
\gamma_{N,2}\gamma_{N,4},
$$

where

$$
\gamma_{N,k}
=1-{k\over N}
+\sqrt{k(1-k/N)(1-1/N)}.
$$

The first Walsh link in (1.1) lies entirely on the column side and has
modulus \(1/q\).  Thus the generic coefficient would be

$$
{ \gamma_{N,2}\gamma_{N,4}\over q}
=0.225293047398
\tag{2.1}
$$

at \(q=32,N=1024\).  Inserting (2.1), after retaining the first balanced
orbit theorem, gives optimized total

$$
0.337715432160,
\tag{2.2}
$$

above \(1/3\) by \(0.004382098827\).  The largest coefficient this orbit can
tolerate in that diagnostic is \(0.158677132585\), so the generic
architecture would need an unexplained factor below \(0.704315\).

This is a failure of a relaxation, not of the signed-permutation plant.

## 3. Exact cubic fixed-pair slice

Define the internally split cubic symbol

$$
A[x;(b,E)]
=\mathbf 1_{x\notin E}\,
M_{1,3}(b,\{x\}\cup E).
\tag{3.1}
$$

For every fixed singleton \(b\) and pair \(E\), the exact signed-permutation
cubic slice table gives

$$
\sum_{x\notin E}|A[x;(b,E)]|^2
\le T_2,
\qquad
T_2={q^2-2q+2\over q^2(q-1)}.
\tag{3.2}
$$

Equation (3.2) includes the distinct-label restriction.  It follows directly
from

$$
M_{1,3}(b,C)=v_3(C)H_N(b,\operatorname{xor}C)
$$

and the exact enumeration of same-label, two-plus-one, and three-distinct
hidden-label cubic supports.

Regard the row vector in (3.1) as a vector in \(\ell_2([N])\).  Factoring
\(A[x;(b,E)]\) as the inner product of the basis vector \(u_x\) and this
column vector proves the Schur multiplier bound

$$
\boxed{\|S_A\|_{S_1\to S_1}\le\sqrt{T_2}.}
\tag{3.3}
$$

The extra row label \(F\) and column labels \(a,e\) merely repeat the same
factor vectors.  At \(q=32\),

$$
T_2=0.0303049395161,
\qquad
\sqrt{T_2}=0.174083139666.
\tag{3.4}
$$

This is the decisive gain over \(\gamma_{1024,2}=2.41018866664\).

A complete \(q=4\) signed-permutation enumeration independently gives

$$
\max_{b,E}
\sum_{x\notin E}|M_{1,3}(b,\{x\}\cup E)|^2
={10\over48},
\tag{3.5}
$$

exactly matching (3.2).

## 4. Complete the adjacent link and restore the quintic mask

Before imposing the within-block distinctness restrictions, the adjacent
moment has the Gram representation

$$
\widetilde M_{3,5}((x,F);(E,e))
=\mathbb E\!\left[
  (X_xY_F)(X_EY_e)
\right].
\tag{4.1}
$$

Both features have unit modulus, so (4.1) is a Schur multiplier of norm at
most one.  Cubic collisions need no second mask: (3.1) is already zero when
\(x\in E\).  Restoring only \(e\notin F\) costs \(\gamma_{N,4}\).

Finally, \(H_N(a,b)\) is a column-only factor of modulus \(1/q\).  Combining
(3.3), (4.1), and the quintic mask gives

$$
\boxed{
\|D_p^{1/2}KD_r^{1/2}\|_1
\le {\gamma_{N,4}\sqrt{T_2}\over q}
\sqrt{\left(\sum p\right)\left(\sum r\right)}.}
\tag{4.2}
$$

The complement and path-reversal cuts have the same bound.  At the target
size,

$$
\boxed{
{\gamma_{1024,4}\sqrt{T_2}\over32}
=0.0162724692796354
<{1\over32}.}
\tag{4.3}
$$

No record-one/record-three triangle sum is needed.  Both adjacent record
sectors remain inside the completed unit-feature Gram multiplier; the exact
cubic slice and the sole remaining distinctness mask carry the quantitative
cost.

## 5. Finite-size consequence

Insert the first balanced-orbit coefficient \(0.093475274578\), replace the
four provisional charges in the present orbit by (4.3), retain the three
physical-orbit diagnostics, and leave every other open split at its
provisional \(1/32\).  Reoptimizing gives

$$
\beta=0.779586351033,
\qquad
\operatorname{TV}_{\rm diagnostic}=0.325563857970,
\tag{5.1}
$$

with

$$
\boxed{0.007769475364}
\tag{5.2}
$$

below \(1/3\).  Because (4.3) is sharper than the old provisional charge,
this second theorem recovers about \(0.000718\) of the margin spent by the
first conservative theorem.

Eight of the 888 balanced entries now have rigorous arbitrary-law upper
bounds from the two chain-aware theorems.  The remaining provisional charges
and the physical lower-witness substitutions still prevent (5.1) from being
a complete one-batch theorem, and no adaptive lift is claimed.

After removing the already proved and explicitly forced orbits, the largest
remaining provisional orbit is

$$
(3,1,1,5):(1,0,0,4),
$$

followed closely by \((3,1,1,5):(1,1,1,2)\).  The next Track A theorem should
continue at that frontier.

## 6. Reproduction

- `searches/adjacent_balanced_cubic_slice_contraction.py` evaluates the
  generic failure, exact coefficient, and updated ledger.
- `tests/adjacent_balanced_cubic_slice_contraction.py` checks (3.2) by full
  \(q=4\) moment enumeration, constructs complete abstract chain tensors
  under arbitrary diagonal laws, and protects the \(N=1024\) constants.
