# Degree-six joint-occupation barrier

Date: 2026-07-14

Status: historical intermediate obstruction.  Exact cut-dependent bounds
remove most of the first OP3-3 obstruction, but the moment relaxation below
still pairs physically incompatible row and column occupations.  The later
occupation-pairing theorem in
`occupation_compatibility_and_promise_concentration.md` supersedes this
barrier and clears every sector through degree eight.  The calculations here
remain useful fixed-mask regressions.

## 1. Profiles included

The diagnostic uses one common distribution on the 210 occupation states

$$
\{n\in\mathbb Z_{\ge0}^4:\ |n|\le6\}.
\tag{1.1}
$$

It includes the refined double-endpoint profile \((3,1,1,3)\) and all four
degree-six profiles

$$
(3,1,1,1),\quad(1,3,1,1),\quad
(1,1,3,1),\quad(1,1,1,3).
\tag{1.2}
$$

For a profile \(a\) and cut \(s\le a\), its occupation term is

$$
\beta^{|a|}\gamma_{a,s}
\sqrt{m_s(\rho)m_{a-s}(\rho)},
\qquad
m_s(\rho)=\mathbb E_\rho\prod_b\binom{n_b}{s_b}.
\tag{1.3}
$$

The double-endpoint constants are those in
`same_middle_and_full_double_endpoint_ledger.md`.

## 2. Refined single-cubic constants

Write \(q=32\), \(n=q-1\), and \(A=q^2-2q+2\).  Translation twirling and
the exact endpoint weight-matrix spectrum give the following
arbitrary-diagonal coefficients for an internal \(1|2\) endpoint-cubic
split:

$$
\gamma_{\mathrm{extreme}}
=\frac{2\sqrt{n(A^2+q^3)}}{q^5}
=0.000324856644\ldots,
\tag{2.1}
$$

when zero or all three chain singletons accompany the selected endpoint
mark, and

$$
\gamma_{\mathrm{balanced}}
=\frac{2\sqrt{n(A^2+q^3)}}{q^4}
=0.010395412614\ldots,
\tag{2.2}
$$

when one or two do.  The pair-difference column Grams are orthogonal in the
balanced case.  Direct \(q=2\) matrices attain the corresponding mask-class
maxima and provide an independent normalization check.

For a cubic decoration on either middle block, the two singleton records
force the L-shape kernel.  Rank--Frobenius over its six marks gives

$$
\gamma_{\mathrm{middle}}=\frac{1}{(q-1)^2}
=\frac1{961}.
\tag{2.3}
$$

Whole-block cuts use the actual path coherences.  Ordinary three-Hadamard
profiles have coefficient \(q^{-3}\) for the empty cut, \(q^{-2}\) for any
singleton or adjacent outer pair, and \(q^{-1}\) for the remaining pair
classes.  A middle-cubic profile instead has two L-shape links of coherence
\(1/[q(q-1)]\); the optimizer uses the resulting complete mask table rather
than a blanket \(1/q\) charge.

These are arbitrary-diagonal bounds for each fixed mask.  They do not yet
exploit that all masks come from the same physical probe law, which is the
remaining opportunity for a compound contraction.

## 3. Concave optimization and result

The sum of (1.3) is concave in \(\rho\).  An interior numerical candidate
supplies a global tangent upper over all 210 states.  At \(q=32\),
\(\beta=5/6\), the result is

$$
F_{6+\mathrm{de}}<0.349675.
\tag{3.1}
$$

At the candidate law, the contributions are

$$
F_{\mathrm{degree}\ 6}=0.2456331472,
\qquad
F_{(3,1,1,3)}=0.1040400214.
\tag{3.2}
$$

The available higher-sector margin is only

$$
0.160358131958,
\tag{3.3}
$$

so the certified numerical upper overshoots it by \(0.189317\).  The two
endpoint-cubic profiles contribute \(0.12176943\) each.  The two
middle-cubic profiles contribute only \(0.00104714\) each.  The endpoint
profiles therefore account for more than 99 percent of the degree-six cost.

The leading individual endpoint terms include whole-block cuts and balanced
internal masks.  This matters: replacing only one advertised coefficient
will not close the ledger unless the proof controls the coupled family of
masks.

## 4. Quantitative next target

If the double-endpoint cost remained at its present candidate value, degree
six could spend at most

$$
0.160358132-0.104040021=0.056318111.
\tag{4.1}
$$

Thus the current degree-six contribution would need an effective factor
below about \(0.229\).  After reserving the already-small middle-cubic cost,
the two endpoint profiles could spend about \(0.0542\), versus their current
\(0.2435\).  The exact target can move because the shared occupation
optimizer will reallocate mass, so these ratios are guidance rather than a
separate proof obligation.

The uniform endpoint-cubic Schur benchmark is

$$
\frac{4(q-1)(q^2-q+1)}{q^4\sqrt{q^2-1}}
=0.003671413\ldots
\quad(q=32).
\tag{4.2}
$$

This is about \(0.353\) of the arbitrary-diagonal balanced coefficient in
(2.2).  It confirms that substantial physical-law structure is available,
but a successful theorem must preserve that structure across the whole
endpoint-cubic mask family.  The concrete next task is therefore:

> Prove or falsify a shared-physical-law contraction for the endpoint-cubic
> profile strong enough that the joint upper in (3.1) drops below
> \(0.160358132\).

Only after that gate passes should the program extend the ledger to the
remaining degree-eight, degree-ten, or degree-twelve profiles, and then to
adaptivity.

Reproduction:

- `searches/single_cubic_weighted_bound.py`;
- `tests/single_cubic_weighted_bound.py`;
- `searches/degree_six_joint_occupation_optimization.py`; and
- `tests/degree_six_joint_occupation_barrier.py`.
