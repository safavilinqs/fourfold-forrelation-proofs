# Occupation compatibility and the degree-eight one-batch ledger

Date: 2026-07-15

Status: proved the occupation-pairing refinement and a sharper promise
concentration theorem.  With all degree-four, degree-six, and degree-eight
signed-permutation sectors included, the resulting (N=1024) one-batch
numerical certificate is below (1/3) after promise conditioning.  This is
not yet a passive dose-six lower bound: degree-ten and degree-twelve sectors
and the adaptive lift remain open.

The global quadratic promise theorem in Section 4 remains valid, but it is
no longer the sharpest bound in the attenuation window used by the current
forced-cut ledger.  The finite-tilt Euclidean refinement is proved in
`euclidean_promise_concentration.md`.

## 1. The compatibility discarded by the old moment relaxation

Fix an odd Fourier profile (a=(a_1,ldots,a_4)) and an occurrence split
(sle a).  For a physical matrix entry, let (S_b) and (T_b) be the row
and column parity supports in block (b).  The selected marks are exactly
(S_b\setminus T_b), while the complementary marks are
(T_b\setminus S_b).  Therefore

$$
|S_b|-s_b=|S_b\cap T_b|
=|T_b|-(a_b-s_b).
\tag{1.1}
$$

If (n_b=|S_b|) and (m_b=|T_b|), every nonzero block consequently obeys

$$
\boxed{m=n+a-2s.}
\tag{1.2}
$$

The earlier relaxation bounded the same split by

$$
\gamma_{a,s}
\sqrt{
\left(\sum_n\rho_n\prod_b\binom{n_b}{s_b}\right)
\left(\sum_m\rho_m\prod_b\binom{m_b}{a_b-s_b}\right)
}.
\tag{1.3}
$$

It pairs every row occupation with every column occupation, including pairs
that cannot be the two sides of one physical Schur entry.  This is why its
endpoint-cubic optimizer can put all its mass on parity-incompatible states;
the direct (q=2) physical matrix for that law is identically zero.

For fixed (s), row occupations (n) and their partners
(m=n+a-2s) form disjoint row/column occupation blocks.  Applying the
arbitrary-diagonal coefficient inside each block and summing their nuclear
norms gives the sharper valid contribution

$$
\boxed{
\gamma_{a,s}
\sum_{n:\,m=n+a-2s}
\sqrt{
\rho_n\rho_m
\prod_b\binom{n_b}{s_b}\binom{m_b}{a_b-s_b}
}.}
\tag{1.4}
$$

Equation (1.4) is no larger than (1.3) by Cauchy--Schwarz, but it retains the
exact zero pattern of the physical matrix.

## 2. The optimization is an exact Perron problem

Let (x_n=\sqrt{\rho_n}).  After summing (1.4) over profiles and splits,
the complete ledger has the form

$$
F(\rho)=\sum_{n<m}c_{nm}\sqrt{\rho_n\rho_m}
=x^TBx,
\qquad \|x\|_2=1,quad x\ge0,
\tag{2.1}
$$

where (B_{nm}=c_{nm}/2) is a nonnegative symmetric matrix on the 210
dose-six occupation states.  Perron--Frobenius therefore gives

$$
\boxed{\sup_\rho F(\rho)=\lambda_{\max}(B).}
\tag{2.2}
$$

This replaces the nonlinear occupation solve by one (210\times210)
eigenvalue calculation.  A positive Perron vector also gives a
Collatz--Wielandt upper.  The current displayed decimal uses floating-point
coefficients with a conservative allowance; intervalization remains final
certificate cleanup.

## 3. Cut-dependent incidence improvements

The degree-eight pass uses the already proved endpoint Fourier and path
coefficients, plus row/column incidence bounds that retain the actual support
families.  The relevant maximum incidence degrees at order (q) are:

$$
D_k^{L}=\bigl(q^2(q-1)^2,\ 3(q-1)^2,\ 2(q-1),\ 1\bigr)_k
\tag{3.1}
$$

for cubic L-shapes,

$$
D_k^{E_3}=\left(
q\binom q3+q^2(q-1)\binom q2,
\binom{q-1}2+(q-1)\binom q2+q(q-1)^2,
q^2-2,1\right)_k
\tag{3.2}
$$

for cubic endpoint record-one supports, and

$$
D_k^{\mathrm{star}}
=\left(q\binom q3(3q-2),
{(q-1)(q-2)(3q-2)\over2},q(q-2),1\right)_k
\tag{3.3}
$$

for a cubic support with record sizes one and three on its two sides.  A
record-three endpoint has

$$
D_k^{E_{3,r=3}}
=\left(\binom q3q^3,\binom{q-1}2q^2,q(q-2),1\right)_k.
\tag{3.4}
$$

For a degree-five endpoint record-one family, the exact incidence degree is
computed by a five-cell column-parity generating function.  At (q=32),

$$
(D_0,ldots,D_5)
=(120731415552,589508865,15811580,46470,1020,1).
\tag{3.5}
$$

For a support family with degrees (D^{(b)}_{s_b}), maximum entry
(\kappa), and complementary split (a-s), row/column energy gives

$$
\gamma_{a,s}
\le\kappa\min\left\{
\sqrt{\prod_bD^{(b)}_{s_b}},
\sqrt{\prod_bD^{(b)}_{a_b-s_b}}
\right\}.
\tag{3.6}
$$

Direct \(q=4\) enumeration verifies every incidence formula used here.
Equation (3.6) sharply improves the blanket balanced-rank charge for middle
cubics, middle quintics, separated double cubics, and both record sectors of
the adjacent double-cubic profiles.

## 4. A sharper promise theorem

Let a biased sign \(\eta\in\{\pm1\}\) have mean \(\beta\).  The Kearns--Saul
inequality gives

$$
\mathbb E e^{t(\eta-\beta)}
\le e^{\kappa_\beta t^2/2},
\qquad
\kappa_\beta={2\beta\over
\log((1+\beta)/(1-\beta))}.
\tag{4.1}
$$

Condition on the two middle noise blocks.  The four-chain statistic is a
bilinear form

$$
{1\over N}\eta_1^TA\eta_4
\tag{4.2}
$$

with \(A\) orthogonal.  Conditioning once more on one endpoint and applying
(4.1) twice shows that its endpoint fluctuation has proxy

$$
{\kappa_\beta(1+\beta^2)\over N}.
\tag{4.3}
$$

The endpoint conditional mean is \(\beta^2Q\), where exact planted
identities reduce \(Q\) to a two-block bilinear Hadamard form.  The same
argument gives (4.3) for \(Q-\beta^2\).  Combining the conditional moment
generating functions proves

$$
\boxed{
F_{4,H}(\widetilde X)-\beta^4
\text{ is subgaussian with proxy }
V_\beta={\kappa_\beta(1+\beta^2)(1+\beta^4)\over N}.}
\tag{4.4}
$$

Thus the bad-promise probability for either sign is at most

$$
\epsilon_\beta
\le\exp\left[-{(\beta^4-1/4)^2\over2V_\beta}\right].
\tag{4.5}
$$

This improves substantially on the earlier variance-only Cantelli bound.

## 5. The degree-eight partial budget

Choose the rational attenuation

$$
\beta={313\over400}=0.7825.
\tag{5.1}
$$

At (N=1024), (4.5) gives

$$
2\epsilon_\beta<0.015743758.
\tag{5.2}
$$

The Perron upper for every profile of total degree four, six, or eight,
including the central and adjacent record-three sectors, is

$$
F_{\le8}<0.281922501.
\tag{5.3}
$$

Consequently

$$
\boxed{
F_{\le8}+2\epsilon_\beta
<0.297666259
<{1\over3},}
\tag{5.4}
$$

leaving (0.03566707) for degrees ten and twelve.  This simultaneously
resolves the former endpoint-cubic and double-endpoint numerical blockers;
neither is the current lead obstruction.

## 6. High-degree continuation

At this checkpoint, the 20 degree-ten and 35 degree-twelve profiles were not
included in (5.3).  The first target diagnostic passed at a common open
coefficient \(1/112\) and failed at \(1/96\), identifying the four
triple-cubic degree-ten profiles as the leading objects.

That target has now been acted on.  Chained squared-slice bounds include all
four triple-cubic profiles together with the already-proved high endpoint
profiles.  At \(\beta=25/32\), the enlarged known-sector Perron upper plus
promise loss is

$$
0.324831904411<{1\over3}.
$$

The next cubic--quintic step proves both separated reversal pairs
\((5,1,3,1)/(1,3,1,5)\) and
\((3,1,5,1)/(1,5,1,3)\).  At \(\beta=781/1000\), the current
known-sector total is \(0.330813664460\), leaving \(0.002519668873\).
The common-coefficient diagnostic for the still-open profiles now passes at
\(1/588\) and fails at \(1/584\), with transition near \(1/585\).  See
`triple_cubic_chained_slice_contraction.md` and
`quintic_slices_and_separated_chain.md` for the proofs and revised frontier.
Those diagnostics predate the physical mixed-orbit cuts and the subsequent
accepted-sector repair.  Chaining the adjacent double-cubic record-three
slices now lowers the degree-eight Perron upper to \(0.232480992267\) and the
conditioned total to \(0.248224749551\).  With the proved high sectors, the
partial total is \(0.279758546919\).  See
`adjacent_record_three_chained_repair.md` for the proof and current
\(1/q\)-classification target.  Do not spend time on adaptivity until the
remaining one-batch profiles are inserted or this route receives a
quantitative obstruction.

Reproduction:

- `searches/occupation_compatible_sector_optimization.py`;
- `searches/attenuation_promise_concentration.py`;
- `searches/endpoint_cubic_physical_orbit.py`;
- `tests/occupation_compatible_sector_optimization.py`; and
- `tests/attenuation_promise_concentration.py`.
