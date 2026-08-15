# Column-cubic/quintic-row contraction

Date: 2026-07-15

Status: rigorous arbitrary-diagonal bound for the fifth balanced orbit.  The
coefficient (0.03118890512) is below the provisional (1/32) and the live
route gate (0.05425062978).  The proof replaces a failing generic
distinct-label multiplier by the exact row energy of the quintic endpoint.

## 1. Exact target

For profile ((3,1,1,5)) and split ((0,0,1,4)), let (Q) be the cubic
support, let (b,c) be the singleton labels, and write the quintic support as
(F\mathbin{\dot\cup}\{e\}), where (|F|=4).  Rows are ((c,F)), columns are
((Q,b,e)), and the exact occurrence tensor is

$$
K_{(c,F),(Q,b,e)}
=M_{3,1}(Q,b)H_N(b,c)M_{1,5}(c,F\cup\{e\}),
\qquad e\notin F.
\tag{1.1}
$$

The desired theorem permits arbitrary correlated nonnegative diagonal row
and column laws.  Uniform-orbit estimates are therefore insufficient.

## 2. Use the quintic row as a Schur feature

Define

$$
B_{(c,F),e}=\mathbf 1_{\{e\notin F\}}
M_{1,5}(c,F\cup\{e\}).
\tag{2.1}
$$

The exact fixed-four endpoint slice gives

$$
\max_{c,F}\sum_{e\notin F}|B_{(c,F),e}|^2
=T_4=1-{4\over N}.
\tag{2.2}
$$

Take (u_{c,F}=(B_{(c,F),e})_e) and take the column feature for
((Q,b,e)) to be the standard basis vector (v_{Q,b,e}=\mathbf e_e).
Then (B_{(c,F),e}=\langle u_{c,F},v_{Q,b,e}\rangle), so Schur
multiplication by the repeated quintic factor has trace-class norm at most

$$
\boxed{\sqrt{T_4}=\sqrt{1-4/N}.}
\tag{2.3}
$$

This factorization already includes the condition (e\notin F).  The
generic centered-disjointness factor would instead cost

$$
\gamma_{N,4}=1-{4\over N}
+\sqrt{4(1-4/N)(1-1/N)},
\tag{2.4}
$$

which is about (2.99121) at (N=1024) and is too expensive for this gate.

## 3. Rank--Frobenius collapse of the remaining chain

Remove (2.1) temporarily.  The remaining matrix is independent of (F,e):

$$
A_{c,(Q,b)}=M_{3,1}(Q,b)H_N(b,c).
\tag{3.1}
$$

For arbitrary row and column laws (p_{c,F}) and (r_{Q,b,e}), duplicate
compression produces

$$
P_c=\sum_Fp_{c,F},
\qquad
R_{Q,b}=\sum_er_{Q,b,e}.
\tag{3.2}
$$

The signed-permutation endpoint formula is

$$
M_{3,1}(Q,b)=v_3(Q)H_N(\operatorname{xor}Q,b),
\qquad |v_3(Q)|\le1.
\tag{3.3}
$$

Both Hadamard factors have modulus (1/\sqrt N).  Hence every entry of
(A) has modulus at most (1/N), while (A) has only (N) rows.  Therefore

$$
\begin{aligned}
\|D_P^{1/2}AD_R^{1/2}\|_1
&\le\sqrt N\,\|D_P^{1/2}AD_R^{1/2}\|_2\\
&\le {1\over\sqrt N}
\sqrt{\left(\sum P\right)\left(\sum R\right)}.
\end{aligned}
\tag{3.4}
$$

At (q=4), the complete compressed matrix has rank (N=16).  A fixed
vertical cubic support with uniform (b,c) laws attains the factor
(1/\sqrt N=1/4), confirming that (3.4) itself cannot be improved without
using the quintic law jointly.

## 4. Arbitrary-law whole-chain theorem

Restore the quintic Schur factor from (2.3).  Equations (2.3) and (3.4) give

$$
\boxed{
\|D_p^{1/2}KD_r^{1/2}\|_1
\le {\sqrt{1-4/N}\over\sqrt N}
\sqrt{\left(\sum p\right)\left(\sum r\right)}.
}
\tag{4.1}
$$

At (N=1024),

$$
\boxed{
c_{(3,1,1,5):(0,0,1,4)}
=0.0311889051224
<{1\over32}
<0.0542506297760.
}
\tag{4.2}
$$

Complement and path reversal prove the same value on the other three cuts.
For comparison, using (2.4) would give (0.09347527458), above the live
gate.  The exact quintic row energy is therefore decision-critical.

## 5. Stress tests

Complete (q=4) enumeration gives

$$
T_4={3\over4},
\qquad
\operatorname{rank}A=16,
\qquad
\max|A_{ij}|={1\over16}.
\tag{5.1}
$$

The fixed-vertical law attains the base coefficient (1/4).  Direct sparse
submatrices of the exact tensor (1.1), with correlated random diagonal laws,
remain below (4.1).  These checks protect the indexing and normalization;
the arbitrary-law theorem follows from the Schur and rank--Frobenius
arguments, not from sampling.

## 6. Ledger consequence

Insert (4.2) together with the four earlier chain-aware theorems.  The exact
210-state diagnostic and extended promise loss reoptimize to

$$
\beta=0.779338939024,
\qquad
\operatorname{TV}_{\rm diagnostic}=0.332365794098.
\tag{6.1}
$$

Thus

$$
\boxed{
{1\over3}-\operatorname{TV}_{\rm diagnostic}
=0.000967539235.
}
\tag{6.2}
$$

Twenty of the 888 balanced entries now have chain-aware arbitrary-law
coefficients.  The largest remaining provisional orbit is

$$
(1,1,3,5):(0,1,1,3),
\tag{6.3}
$$

with Perron contribution (0.001161909690).  Reoptimizing while varying only
that four-cut orbit gives the next acceptance gate

$$
\boxed{
c_{(1,1,3,5):(0,1,1,3)}
<0.0570749885142
=1.82639963245/q.
}
\tag{6.4}
$$

## 7. Reproduction and scope

- `searches/column_cubic_quintic_row_contraction.py` computes the coefficient,
  five-theorem ledger, reranked orbit, and next gate.
- `tests/column_cubic_quintic_row_contraction.py` enumerates the (q=4)
  quintic row energy and base chain, checks saturation of (3.4), and
  stress-tests the exact tensor under correlated diagonal laws.
- `./run_round3_checks.sh` runs the complete inherited and Round 3 suite.

This theorem does not certify the remaining 868 balanced entries, replace
the physical lower witnesses by arbitrary-law upper bounds, or prove the
adaptive lift.  It passes one more finite-size gate and identifies exact
one-sided row energy as a reusable alternative to generic distinctness
multipliers.
