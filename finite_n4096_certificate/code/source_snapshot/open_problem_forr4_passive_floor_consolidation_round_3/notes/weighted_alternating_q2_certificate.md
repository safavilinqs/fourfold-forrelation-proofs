# Weighted alternating \(q=2\) certificate

Date: 2026-07-14

Status: certified numerical bound for arbitrary correlated diagonal row and
column laws in one fixed alternating occurrence split at \(q=2\).  This is
a falsification check for OP3-2, not the physical joint-split theorem and not
an extrapolation to \(q=32\).

## 1. Question

For the complete \(64\times144\) alternating matrix \(K\), define

$$
\Phi(p,r)=
\left\|D_p^{1/2}KD_r^{1/2}\right\|_1,
\qquad
\sum p=\sum r=1.
\tag{1.1}
$$

The uniform value is

$$
\Phi_{\rm uniform}=0.4715918158911433.
\tag{1.2}
$$

The test asks whether arbitrary diagonal concentration already destroys this
uniform cancellation at the smallest dimension.

## 2. Concavity and orbit reduction

Equation (1.1) is the root fidelity between \(D_p\) and
\(KD_rK^*\), hence is jointly concave in \((p,r)\).  The exact \(q=2\)
symmetries make the row law uniform at a stationary point.  They split the
144 columns \((E,F,c)\) into:

- 48 high-orbit columns with
  \(\langle\sigma(E),\sigma(F)\rangle=0\); and
- 96 low-orbit columns with nonzero binary inner product.

Let \(t\) be total mass on the high orbit.  Integer Gram matrices for the
two orbits commute.  Their exact joint eigenvalue pairs, after the
normalizations in the search artifact, have multiplicities

$$
(0,1/24576)^{[24]},\quad
(1/12288,0)^{[12]},\quad
(1/12288,1/12288)^{[12]},
$$

$$
(1/4096,0)^{[4]},\quad
(1/12288,1/6144)^{[12]}.
\tag{2.1}
$$

Therefore

$$
\Phi(t)=
{3\over8\sqrt6}\sqrt{1-t}
+{\sqrt3+1\over16}\sqrt t
+{\sqrt3\over16}
+{\sqrt3\over16}\sqrt{2-t}.
\tag{2.2}
$$

The unique stationary point is

$$
t_*=0.375486496922672\ldots
\tag{2.3}
$$

and gives

$$
\Phi(t_*)=
0.4718447392358088\ldots.
\tag{2.4}
$$

## 3. Global diagonal check

At the law in (2.3), the direct polar-factor derivative is constant across
all 64 row coordinates and all 144 column coordinates to
\(7\times10^{-16}\).  Joint concavity then supplies a supporting-hyperplane
upper bound.  With a conservative \(10^{-12}\) numerical allowance, the
machine-checked result is

$$
\boxed{
\sup_{p,r}\Phi(p,r)<0.471845.
}
\tag{3.1}
$$

The uniform-to-optimal increase is only \(0.0537\%\).  The safe
row/column-energy coefficient at \(q=2\) is \(\sqrt{3/8}=0.612372\ldots\),
so the actual arbitrary-diagonal optimum is only \(0.77052\) of that bound.

## 4. Interpretation and limitation

This check finds no small-\(q\) diagonal concentration catastrophe.  It
supports continuing the weighted route and shows that the local row/column
bound is already substantially loose at \(q=2\).

It does not establish OP3-2:

- it treats one occurrence split, not all splits sharing one probe law;
- it uses \(q=2\), where the type-I cubic orbit is absent;
- the final supporting bound is computer-assisted rather than interval
  arithmetic; and
- no scaling claim is made.

The next useful calculation is the corresponding symmetry-reduced \(q=4\)
problem.  Forming the full \(4096\times230400\) matrix is unnecessary; the
orbit Gram should be assembled from endpoint incidences.

Reproduction:

- searches/weighted_alternating_q2_certificate.py; and
- tests/weighted_alternating_q2_certificate.py.
