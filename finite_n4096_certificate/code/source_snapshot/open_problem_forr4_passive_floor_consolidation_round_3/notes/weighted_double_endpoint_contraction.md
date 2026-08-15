# Weighted double-endpoint fixed-split contraction

Date: 2026-07-14

Status: the two quantitatively relevant alternating occurrence orientations
have \(q=32\) arbitrary-diagonal bounds.  The mixed result in this note has
since been sharpened to \(0.020343\), and the same-middle masks have been
solved exactly.  The complete 64-cut occupation ledger now fits below its
allocated margin; see `same_middle_and_full_double_endpoint_ledger.md`.

## 1. Translation twirl

For any fixed occurrence flattening \(K\), put

$$
\Phi_K(p,r)=\|D_p^{1/2}KD_r^{1/2}\|_1,
\qquad \sum p=\sum r=1.
\tag{1.1}
$$

This is the root fidelity between \(D_p\) and \(KD_rK^*\), so it is jointly
concave in \((p,r)\).  Independent translations of the four physical
coordinate blocks act on \(K\) by row/column permutations and diagonal
Walsh signs.  Averaging a diagonal law over these translations therefore
cannot decrease (1.1).

For an endpoint pair \(E\), its XOR \(x=\sigma(E)\ne0\) is unchanged by
translation.  The twirled law is uniform over the \(N/2\) unordered pairs
with a fixed difference.  Thus:

- for rows \((i,b,d)\) and columns \((E,F,c)\), the remaining variable is a
  joint law \(Q_{x,y}\) on two pair differences;
- for rows \((i,b,F)\) and columns \((E,c,d)\), the remaining variables are
  separate row and column difference laws \(P_y,Q_x\).

These reductions retain arbitrary diagonal correlations allowed by the
fixed occurrence split.

## 2. Same endpoint orientation

Let \(U_i(E)=(q-1)w_i(E)\) and

$$
C(h,x)=
\sum_{\substack{E\ {\rm unordered}\\\sigma(E)=x}}
U_0(E)U_h(E),
\qquad
D(\mu,x)=\sum_hC(h,x)\chi_\mu(h).
\tag{2.1}
$$

For the twirled joint difference law \(Q\), the complete block spectrum is

$$
\lambda_{\alpha,\beta}(Q)
={4\over N^8(q-1)^4}
\sum_{x,y\ne0}
Q_{x,y}D(\beta,x)D(\alpha\oplus x,y),
\tag{2.2}
$$

and

$$
\Phi(Q)=N\sum_{\alpha,\beta}
\sqrt{\lambda_{\alpha,\beta}(Q)}.
\tag{2.3}
$$

Equations (2.2)--(2.3) use exact integer Walsh data before square roots are
taken.  At \(q=32\), the optimizer is uniform on the 465 pairs of nonzero
same-hidden-column differences satisfying
\(\langle x,y\rangle=0\).  The direct supporting gradient is constant on
the support and below it elsewhere by \(2.13\times10^{-5}\).  With a
conservative \(10^{-12}\) allowance,

$$
\boxed{
\sup_Q\Phi(Q)<0.010905.
}
\tag{2.4}
$$

For comparison, the values at \(q=4,8,16,32\) are

| \(q\) | arbitrary-diagonal coefficient |
|---:|---:|
| 4 | 0.171034270767 |
| 8 | 0.074160451046 |
| 16 | 0.029266144596 |
| 32 | 0.010904820188 |

This is the orientation with one endpoint coordinate on the row side at
both decorated endpoints.  It is not the mixed orientation that attained
the old worst coefficient.

## 3. Mixed endpoint orientation

Consider rows \((i,b,F)\) and columns \((E,c,d)\).  Replace each unordered
pair of difference \(x\) by its two ordered representatives.  Splitting its
diagonal mass equally preserves the weighted singular values.

After Walsh transforms in the pair bases and endpoint singleton labels, the
matrix separates into blocks.  If

$$
\widehat f_x(\mu)
=\sum_z w_0(\{z,z\oplus x\})\chi_\mu(z),
\tag{3.1}
$$

then

$$
|\widehat f_x(\mu)|^2
={2D(\mu,x)\over(q-1)^2}.
\tag{3.2}
$$

The weighted trace norm is exactly

$$
\Phi(P,Q)
={1\over N^3}
\sum_{\mu,\nu}
\left\|
D_{a^\nu}W_N^{(\ne0)}D_{b^\mu}
\right\|_1,
\tag{3.3}
$$

where

$$
a^\nu_y=\sqrt{P_y}\,|\widehat f_y(\nu)|,
\qquad
b^\mu_x=\sqrt{Q_x}\,|\widehat f_x(\mu)|,
\tag{3.4}
$$

and \(W_N^{(\ne0)}\) is the unnormalized Walsh matrix with the zero
difference column removed.  Direct \(q=2\) matrices match (3.3) for random
correlated diagonal laws.

## 4. Three difference types

Write a nonzero difference as \(x=(r,c)\).  There are three types:

$$
V=(r\ne0,c=0),\qquad
H=(r=0,c\ne0),\qquad
D=(r\ne0,c\ne0),
\tag{4.1}
$$

with cardinalities \(n,n,n^2\), where \(n=q-1\).

For a probability law supported on one type \(s\), define

$$
A_s=
\sup_P\sum_\mu
\sqrt{\sum_{x\in s}P_x|\widehat f_x(\mu)|^2}.
\tag{4.2}
$$

The incidence action is transitive within each type, and the expression is
concave in \(P\), so the supremum occurs at the uniform law.  Exact endpoint
Fourier counts give, with

$$
a={q^2-2q+2\over n},\qquad
s=\sqrt{{2(q-2)\over n}},\qquad
d={\sqrt{2(q^2-2q+2)}\over n^2},
\tag{4.3}
$$

$$
A_V=2+na+n s+n^2s,
\tag{4.4}
$$

$$
A_H=2+n s+2+n s,
\tag{4.5}
$$

$$
A_D=2+n s+s+n^2d.
\tag{4.6}
$$

For a block from row type \(s\) to column type \(t\), rank--Frobenius in
(3.3) gives

$$
B_{s,t}
={\sqrt{\min(|s|,|t|)}\over N^3}A_sA_t.
\tag{4.7}
$$

If \(u_s=\sqrt{P(s)}\) and \(v_t=\sqrt{Q(t)}\), triangle inequality across
the nine type blocks gives

$$
\Phi(P,Q)\le u^TBv
\le\|B\|_{\rm op}
\le\max_s\sum_tB_{s,t}.
\tag{4.8}
$$

At \(q=32\),

$$
B=
\begin{pmatrix}
0.0284923593142&0.00109705412969&0.00109859688729\\
0.00109705412969&0.000042240368732&0.000042299770222\\
0.00109859688729&0.000042299770222&0.000235846351798
\end{pmatrix}.
\tag{4.9}
$$

Therefore

$$
\boxed{
\Phi(P,Q)
\le0.0306880103312.
}
\tag{4.10}
$$

This was the first analytic bound.  Exact treatment of its dominant
vertical/vertical type block, followed by a rational Collatz certificate,
now improves (4.10) to

$$
\Phi(P,Q)<0.020343.
\tag{4.11}
$$

The old mixed fixed-split coefficient was

$$
0.1232155202300.
$$

The new bound is a factor \(0.249060\) of the old one, below the
\(0.280708\) local improvement target.

## 5. Superseded fixed-split ledger

For the barrier occupation \((2,1,1,2)\), only endpoint splits one and two
have nonzero occurrence mass.  Use:

- \(1/q\) when the two middle singleton marks lie on the same side;
- (2.4) for alternating middle marks and equal endpoint orientations; and
- (4.10) for alternating middle marks and mixed endpoint orientations.

The resulting safe local triangle ledger is

$$
0.8327440826492
\tag{5.1}
$$

before attenuation, and

$$
\left({5\over6}\right)^8
0.8327440826492
=0.1936696585915.
\tag{5.2}
$$

The available margin is \(0.160358131958\), so the sharpened ledger still
overshoots by

$$
0.033311526634.
\tag{5.3}
$$

This obstruction is superseded.  The exact same-middle coefficients are
\(0.00094703\) and \(0.00180967\), not \(1/q\), and the refined alternating
mixed coefficient is (4.11).  On the same diagnostic state the attenuated
ledger becomes \(0.0632673\).  More importantly, optimizing all 64 cuts over
all 210 dose-six occupation states gives the conservative numerical upper
\(0.115981\), below the allocated \(0.160358\) margin.  The additional gain
comes from the one-balanced/one-unsplit same-middle bound in the companion
note.

## 6. Next theorem target

The remaining finite-size task is to combine this profile with the other
degree-six through degree-twelve profiles without allocating the full TV
margin to each profile separately.  Intervalization of the occupation
supporting hyperplane remains worthwhile before a final theorem claim.

Reproduction:

- searches/weighted_same_orientation_certificate.py;
- searches/mixed_endpoint_weighted_bound.py;
- searches/same_middle_weighted_bound.py;
- searches/double_endpoint_occupation_optimization.py; and
- the corresponding regressions in `tests/`.
