# Euclidean finite-tilt promise concentration

Date: 2026-07-15

Status: proved analytic repair for the first opposite-endpoint forced-cut
diagnostic.  Section 6 records the subsequent hybrid extension needed by
the second physical witness, and Section 7 records the later two-split
extension.  None completes the one-batch ledger or the adaptive passive
lower bound.

## 1. Decision result

The mixed-orbit witness in
opposite_endpoint_mixed_orbit_obstruction.md made the old independent scalar
ledger exceed \(1/3\) by

$$
0.000850121165945.
$$

The smallest proposed repair was a \(0.9448\%\) improvement in the promise
proxy.  A finite-tilt use of the exact biased-sign moment generating
function supplies more than this.  With the four forced cuts fixed at the
physical \(q=32\) witness value and attenuation reoptimized, the diagnostic
becomes

$$
\beta=0.779512334639,
\qquad
\text{Perron upper}+\text{promise loss}
=0.332335953237,
$$

which is below \(1/3\) by

$$
0.000997380096833.
\tag{1.1}
$$

Decision: the explicit mixed-orbit witness no longer forces a
signed-permutation pivot.  The scalar route is viable again, but its margin
is narrow and the remaining high-degree profiles are still absent.

## 2. The scalar monotone branch

Let \(\eta\in\{-1,+1\}\) have mean \(\beta\in(0,1)\), let
\(\xi=\eta-\beta\), and define the lower-tail log-MGF

$$
\psi_\beta(s)
=\log\mathbb E e^{-s\xi},
\qquad s\ge0.
\tag{2.1}
$$

Set

$$
L_\beta=\log{1+\beta\over1-\beta}.
\tag{2.2}
$$

Then

$$
{\psi_\beta(s)\over s^2}
\quad\text{is increasing on }(0,L_\beta],
\tag{2.3}
$$

and its endpoint is the Kearns--Saul equality point:

$$
\psi_\beta(L_\beta)=\beta L_\beta.
\tag{2.4}
$$

For completeness, write \(p=(1+\beta)/2\) and
\(B=(\eta+1)/2\sim\operatorname{Ber}(p)\).  Then

$$
{\psi_\beta(s)\over s^2}
=4g_p(-2s),
$$

where \(g_p(t)=t^{-2}\log\mathbb E e^{t(B-p)}\).  Schlemm proves that
\(g_p\) is strictly unimodal with unique maximum
\(2\log((1-p)/p)=-2L_\beta\).  On \(s\in[0,L_\beta]\), this gives (2.3)
directly.  See Eckhard Schlemm, “The Kearns--Saul inequality for Bernoulli
and Poisson-binomial distributions,” Theorem 1.1,
[arXiv:1405.4496](https://arxiv.org/abs/1405.4496),
[published DOI](https://doi.org/10.1007/s10959-014-0564-x).

## 3. Exact Euclidean packing at finite tilt

Let \(a=(a_i)\) be real coefficients with \(A=\lVert a\rVert_2\).  For
\(0\le tA\le L_\beta\), independence and (2.3) give

$$
\log\mathbb E\exp\left[-t\sum_i a_i\xi_i\right]
\le \psi_\beta(tA).
\tag{3.1}
$$

Indeed, the lower-tail MGF dominates the upper-tail MGF for a positive-mean
sign, so the \(i\)-th term is at most
\(\psi_\beta(t|a_i|)\).  Monotonicity of the ratio then gives

$$
\psi_\beta(t|a_i|)
\le {a_i^2\over A^2}\psi_\beta(tA),
$$

and summing proves (3.1).  This is the gain over the previous proof: all
independent linear terms in one block are aggregated at their Euclidean norm
before paying the exact finite-tilt MGF.

## 4. The four-block martingale

Let \(F\) be the four-chain statistic of the independently attenuated exact
plant.  Average the noise blocks from one end and write the reverse Doob
decomposition

$$
F-\beta^4=D_0+D_1+D_2+D_3.
\tag{4.1}
$$

After the first \(j\) blocks have been replaced by their mean, \(D_j\) is a
linear form in the next centered noise block.  The exact planted identities
collapse the averaged prefix to a sign vector.  Conditional on all later
blocks, the remaining signed diagonal and normalized Hadamard factors are
isometries, so its coefficient norm is exactly

$$
A_j={\beta^j\over\sqrt N}.
\tag{4.2}
$$

Iteratively integrate the four independent blocks using (3.1).  For every
\(0\le t\le L_\beta\sqrt N\),

$$
\log\mathbb E e^{-t(F-\beta^4)}
\le
\sum_{j=0}^3
\psi_\beta\left({t\beta^j\over\sqrt N}\right).
\tag{4.3}
$$

Taking the explicit endpoint tilt \(t=L_\beta\sqrt N\) gives the rigorous
one-sided promise failure bound

$$
\epsilon_\beta^{\rm E}
\le
\exp\left[
-L_\beta\sqrt N(\beta^4-1/4)
+\sum_{j=0}^3\psi_\beta(L_\beta\beta^j)
\right].
\tag{4.4}
$$

The negative hypothesis has the same bound by the sign symmetry of the exact
plant, so conditioning costs at most \(2\epsilon_\beta^{\rm E}\).

This is a finite-tilt statement.  It need not improve the global quadratic
Kearns--Saul bound for every \(\beta\); it improves it in the narrow
attenuation window relevant to the \(N=1024,D=6\) ledger.

## 5. Numerical comparison and scope

At the old forced-ledger optimizer \(\beta=0.780899845855\), (4.4) gives

$$
2\epsilon_\beta^{\rm E}=0.018188754783,
$$

compared with the old quadratic loss \(0.01972476027\).  Reoptimizing the
Perron-plus-promise objective gives (1.1).

What is now proved:

- the exact scalar MGF packing (3.1);
- the four-block promise bound (4.4); and
- a passing forced-cut diagnostic for the explicit mixed-orbit witness.

What is not proved:

- a complete degree-twelve one-batch bound;
- a joint treatment of every still-open split and representation record;
- an interval or exact certificate for the final 210-state optimization; or
- any adaptive finite-size passive lower bound.

Reproduction:

- searches/attenuation_promise_concentration.py implements (2.1) and (4.4);
- tests/attenuation_promise_concentration.py checks the scalar branch, exact
  finite signed-weight products, all four conditional coefficient norms at
  \(q=4\), and the finite-\(N\) value;
- searches/opposite_endpoint_vertical_mixture_witness.py reoptimizes the
  repaired forced ledger; and
- tests/opposite_endpoint_mixed_orbit_obstruction.py protects the old
  obstruction and the new passing diagnostic side by side.

## 6. Subsequent hybrid extension

The transposed physical witness later forced the finite-tilt optimum just
beyond the endpoint allowed in (4.3).  This does not invalidate the scalar
route.  For each martingale block separately, retain the exact branch while
\(s\beta^j\le L_\beta\), and use the global Kearns--Saul parabola afterward:

$$
\Phi_j(s)=
\begin{cases}
\psi_\beta(s\beta^j),&s\beta^j\le L_\beta,\\
\frac12\kappa_\beta(s\beta^j)^2,&s\beta^j>L_\beta,
\end{cases}
\qquad
\kappa_\beta={2\beta\over L_\beta}.
\tag{6.1}
$$

Iterated conditioning now gives the valid bound

$$
\log\mathbb P(F-\beta^4\le-(\beta^4-1/4))
\le
-s\sqrt N(\beta^4-1/4)+\sum_{j=0}^3\Phi_j(s)
\tag{6.2}
$$

for every \(s\ge0\).  With both certified witness orbits inserted, its
reoptimized diagnostic is

$$
\beta=0.779767976462,\qquad
\text{total}=0.332839308989,
$$

below \(1/3\) by \(0.000494024345\).  Exactly one martingale block uses the
global branch at the optimizer.  The full certificate and the resulting
new frontier are in
`notes/transposed_dominant_class_and_hybrid_repair.md`.

## 7. Two-split Euclidean extension

The hybrid theorem is globally valid but pays its quadratic branch as soon
as a block passes \(L_\beta\).  Exact Euclidean packing actually remains
valid farther.  For

$$
a_\beta(s)={\psi_\beta'(s)\over s},
$$

direct differentiation shows that its derivative changes sign only once.
A pairwise rotation at fixed Euclidean norm therefore reduces the extremal
coefficient vector to \(k\) equal nonzero magnitudes.  If \(G_\beta\) is the
first root above \(L_\beta\) of

$$
\psi_\beta(G_\beta)
=2\psi_\beta(G_\beta/\sqrt2),
$$

then \(k=2\) is the strongest competitor before that root, and

$$
\sum_i\psi_\beta(|y_i|)
\le\psi_\beta(\lVert y\rVert_2),
\qquad \lVert y\rVert_2\le G_\beta.
\tag{7.1}
$$

Applying (7.1) blockwise and keeping the two physical witnesses fixed gives

$$
\beta=0.779868743331,qquad
\text{total}=0.332675996818,
$$

below \(1/3\) by \(0.000657336515\).  This was the promise baseline before
the subsequent accepted-sector chained-slice repair; the promise theorem
itself remains in force.
The proof, derivative calculation, and adjacent-profile consequence are in
`notes/adjacent_cubic_quintic_record_gate.md`; the implementation and
regression remain in `searches/attenuation_promise_concentration.py` and
`tests/attenuation_promise_concentration.py`.
