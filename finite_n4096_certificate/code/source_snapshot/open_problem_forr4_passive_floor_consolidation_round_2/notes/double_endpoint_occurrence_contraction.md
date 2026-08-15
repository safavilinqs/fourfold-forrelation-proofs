# Double-endpoint occurrence-split contraction

Date: 2026-07-14

Status: proved for every fixed unordered occurrence split of
\((3,1,1,3)\) and arbitrary correlated diagonal weights, before the joint
sum over shared occurrence assignments.  The worst local coefficient is
\(0.1232\) at \(q=32\), or about \(0.0287\) after degree-eight attenuation.

## 1. Exact one-endpoint slice energies

For a cubic endpoint linked to a singleton, let \(E_k\) be the largest
unweighted squared row norm when \(k\) of the three endpoint coordinates,
together with the neighboring singleton when it is on that side, are held
fixed and the other \(3-k\) endpoint coordinates are summed.  Exact
support counting and the two cubic link amplitudes give

$$
\boxed{
E_0={q^2+2\over6},\quad
E_1={q^2+2\over2q^2},\quad
E_2={q^2-2q+2\over q^2(q-1)},\quad
E_3={1\over q^2}.
}
\tag{1.1}
$$

For example, with one endpoint coordinate fixed, the opposite pair can
either share its hidden label or form a type-II decoration.  Summing their
squared amplitudes gives

$$
{\binom{q-1}{2}\over q^2}
+{(q-1)\binom q2+q(q-1)^2\over q^2(q-1)^2}
={q^2+2\over2q^2}.
$$

The other formulas follow from the same two orbit types.  Complete
\(q=4\) enumeration verifies all four values.

## 2. Two endpoints around opposite singletons

Put the two middle singleton coordinates on opposite physical sides.  If
the first and fourth cubic endpoints contribute \(k\) and \(\ell\) marks
to the ket side, the squared row norm of the full four-block kernel is at
most

$$
E_kE_\ell.
\tag{2.1}
$$

This follows by summing the endpoint supports first and then using
\(\sum_c|H_N(b,c)|^2=1\) on the central singleton link.  Reversing ket and
bra gives squared column norm at most

$$
E_{3-k}E_{3-\ell}.
\tag{2.2}
$$

For any matrix \(W\) and diagonal masses \(P,Q\), Schatten Holder gives

$$
\|D_p^{1/2}WD_q^{1/2}\|_1
\le\min\{\sqrt{R},\sqrt{C}\}\sqrt{PQ},
\tag{2.3}
$$

when every unweighted row and column has squared norm at most \(R,C\).
Applying (2.1)--(2.2), the worst internal split is \((k,\ell)=(1,2)\) or
its reversal:

$$
\boxed{
\gamma_{3,1,1,3}(q)
=\sqrt{E_1E_2}
=\sqrt{
{(q^2+2)(q^2-2q+2)\over2q^4(q-1)}
}.
}
\tag{2.4}
$$

If the middle singleton coordinates lie on the same side, their Hadamard
entry is an internal factor \(1/q\) and the remaining independent endpoint
links form a Gram dressing.  If both cubic supports are unsplit, the
Bessel-refined whole-block path theorem gives \(1/q\).  Thus (2.4), which
is larger for \(q\ge4\), covers every fixed occurrence placement.

At \(q=32\),

$$
\gamma_{3,1,1,3}(32)=0.1232\ldots,
\qquad
\beta^8\gamma_{3,1,1,3}(32)=0.0287\ldots
$$

for \(\beta=5/6\).

## 3. Scope

The proof uses the actual unordered endpoint supports and therefore pays
no labeled-mark factorial.  It permits arbitrary correlated diagonal
weights inside one split.  Common ancillary repetitions can be aggregated
when their marked incidence map is injective.

It does not sum all split masks against one passive occupation law.  That
joint Bessel step is still required before (2.4) can enter a transcript
budget.  The coefficient is nevertheless small enough that the
degree-eight attenuation leaves useful finite-size slack.

Reproduction: searches/double_endpoint_slice_energies.py.
