# Leading balanced distinct-label contraction

Date: 2026-07-15

Status: rigorous arbitrary-diagonal fixed-split upper bound for the first
previously unforced dose-six-relevant orbit.  The coefficient is
\(0.09347527458\), not the provisional \(1/32\), but the current coarse
completion ledger remains below \(1/3\) after the replacement.

## 1. The exact occurrence tensor

Take profile

$$
(3,1,1,5)
$$

and split \((0,1,0,4)\).  Write the cubic endpoint support as \(Q\), the two
middle singleton labels as \(b,c\), and the quintic endpoint support as
\(F\mathbin{\dot\cup}\{e\}\), where \(|F|=4\).  Rows are \((b,F)\), columns
are \((Q,c,e)\), and the exact distinct-label tensor is

$$
K[(b,F),(Q,c,e)]
=M_{3,1}(Q,b)H_N(b,c)M_{5,1}(F\cup\{e\},c),
\qquad e\notin F.
\tag{1.1}
$$

This corrects the possible temptation to reuse the already studied
\((2,0,1,2)\) opposite-endpoint reshuffling.  The two cuts have different
row and column fibers and different physical weighted norms.

The theorem below allows arbitrary nonnegative diagonal row and column
weights \(p_{b,F}\) and \(r_{Q,c,e}\).  Thus it is an upper bound for every
physical passive probe law, not a uniform-orbit calculation.

## 2. The distinct-label Schur factor

Let the ground-set size be \(N\), let \(|F|=k\), and define

$$
D(F,e)=\mathbf 1_{e\notin F}.
$$

Put

$$
\alpha=1-{k\over N},\qquad
x_F=\mathbf 1_F-{k\over N}\mathbf 1,\qquad
y_e=\mathbf 1_e-{1\over N}\mathbf 1.
$$

Then

$$
D(F,e)=\alpha-\langle x_F,y_e\rangle,
\tag{2.1}
$$

with

$$
\|x_F\|^2=k\alpha,
\qquad
\|y_e\|^2=1-{1\over N}.
$$

For any \(s>0\), factor (2.1) as the inner product of

$$
u_F=(\sqrt\alpha,sx_F),
\qquad
v_e=(\sqrt\alpha,-s^{-1}y_e).
$$

Choosing

$$
s^4={1-1/N\over k\alpha}
$$

minimizes the product of the two uniform norm bounds and gives

$$
\boxed{
\gamma_{N,k}
=1-{k\over N}
+\sqrt{k\left(1-{k\over N}\right)
              \left(1-{1\over N}\right)}.}
\tag{2.2}
$$

Consequently Schur multiplication by \(D\) has trace-class norm at most
\(\gamma_{N,k}\).  This is the standard Hilbert-space factorization bound
for a Schur multiplier: a symbol
\(D_{ij}=\langle u_i,v_j\rangle\) acts on \(S_1\) with norm at most
\(\max_i\|u_i\|\max_j\|v_j\|\).

At \(N=1024,k=4\),

$$
\gamma_{1024,4}=2.991208786480118.
\tag{2.3}
$$

## 3. Complete the quintic link, then collapse the chain

Before imposing \(e\notin F\), the quintic moment has the Gram form

$$
\widetilde M_{5,1}(F;(e,c))
=\mathbb E\!\left[X_F(X_eY_c)\right].
\tag{3.1}
$$

Both features in (3.1) have unit modulus.  Its Schur multiplier norm is
therefore at most one.  Equations (2.1)--(3.1) show that the actual
distinct-label quintic factor in (1.1) has Schur multiplier norm at most
\(\gamma_{N,4}\).

Remove that factor temporarily.  The remaining weighted matrix is constant
over \(F\) and \(e\), so the duplicate rows and columns compress isometrically
to weights

$$
P_b=\sum_Fp_{b,F},
\qquad
R_{Q,c}=\sum_er_{Q,c,e}.
$$

For the signed-permutation endpoint,

$$
M_{3,1}(Q,b)=v_3(Q)H_N(\operatorname{xor}Q,b),
\qquad |v_3(Q)|\le1.
\tag{3.2}
$$

Since \(H_N(x,b)H_N(b,c)=q^{-1}H_N(b,x\oplus c)\), with
\(q=\sqrt N\), the compressed matrix is a column-duplicated Hadamard
matrix.  Its effective column weights are

$$
W_z={1\over q^2}
\sum_{Q,c:\,\operatorname{xor}Q\oplus c=z}
|v_3(Q)|^2R_{Q,c},
\qquad
\sum_zW_z\le{1\over q^2}\sum_{Q,c}R_{Q,c}.
\tag{3.3}
$$

Schatten Holder and unitarity of \(H_N\) give

$$
\|D_P^{1/2}H_ND_W^{1/2}\|_1
\le\sqrt{\left(\sum_bP_b\right)
          \left(\sum_zW_z\right)}.
\tag{3.4}
$$

Restoring the distinct quintic labels using (2.2) proves

$$
\boxed{
\|D_p^{1/2}KD_r^{1/2}\|_1
\le {\gamma_{N,4}\over\sqrt N}
\sqrt{\left(\sum p\right)\left(\sum r\right)}.}
\tag{3.5}
$$

The same bound holds for the complement and path-reversal cuts.  At the
target size its coefficient is

$$
\boxed{
{\gamma_{1024,4}\over32}
=0.09347527457750368.}
\tag{3.6}
$$

The proof is chain-aware: the \(1/q\) comes from the exact cubic--Walsh
collapse, while the factor near three is solely the price of restoring the
multilinear \(4|1\) distinct-label mask.

## 4. Finite-size consequence

The previous route-selection ledger charged all four cuts in this orbit at
the unproved value \(1/32\).  Replace those four charges by the rigorous
upper (3.6), retain the three earlier physical-orbit diagnostics, and leave
every other open split at the provisional \(1/32\) target.  Reoptimizing the
attenuation and extended promise bound gives

$$
\beta=0.779557316457838,
\qquad
\operatorname{TV}_{\rm diagnostic}=0.326281860918323,
\tag{4.1}
$$

so the ledger retains

$$
\boxed{0.007051472415011}
\tag{4.2}
$$

below \(1/3\).  Equation (4.1) is still not a complete one-batch theorem,
because the other open \(1/32\) charges and the physical-witness charges are
not all arbitrary-law upper bounds.  It proves something narrower but
decision-relevant: the first unforced balanced orbit can be bounded
rigorously without killing the realistic-size route.

Under the updated diagnostic, this proved orbit itself has the largest
Perron contribution because of its conservative factor near three.  The
largest still-unresolved provisional orbit is now

$$
(1,1,3,5):(0,0,1,4),
$$

followed closely by the remaining internally split
\((3,1,1,5)\) orbits.  Improving (3.6) is optional unless later rigorous
charges consume the remaining margin; the next new theorem should attack
the unresolved orbit rather than polish this factor immediately.

## 5. Scope and decision

- This is an arbitrary-diagonal fixed-split upper bound, not a physical
  lower witness or a uniform-weight spectrum.
- It controls four of the 888 balanced dose-six-relevant entries.
- It does not prove the provisional \(1/q\) coefficient; the exact safe
  coefficient is about \(2.99121/q\).
- It does not address adaptive posterior selection.
- It preserves the feasibility of the one-batch \(N=1024,D=6\) program with
  more than \(0.00705\) diagnostic slack.

Reproduction:

- `searches/leading_balanced_disjointness_contraction.py` evaluates (2.2),
  (3.6), and the reoptimized ledger (4.1).
- `tests/leading_balanced_disjointness_contraction.py` checks the exact
  centered factorization, constructs random complete chain tensors at small
  dimension, verifies (3.5), and protects the \(N=1024\) constants.
