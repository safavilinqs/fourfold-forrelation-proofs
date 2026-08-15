# Universal weighted link Gram contraction

Date: 2026-07-14

Status: proved for every support sector and every planted-pair law.  This
closes the single-link even-decoration problem for both exact-plant
candidates.

## 1. Statement

Let \((\Omega,\nu)\) be any probability space.  Let
\(\{f_A\}_{A\in\mathcal A}\) and
\(\{g_B\}_{B\in\mathcal B}\) be unit-modulus features and define

$$
M(A,B)=\mathbb E_{\omega\sim\nu}
\left[\overline{f_A(\omega)}g_B(\omega)\right].
$$

For arbitrary nonnegative diagonal weights \(p_A,q_B\), put

$$
P=\sum_Ap_A,\qquad Q=\sum_Bq_B.
$$

Then

$$
\boxed{
\|D_p^{1/2}MD_q^{1/2}\|_1\le\sqrt{PQ}.
}
\tag{1.1}
$$

No orthogonality, representation decomposition, or bound on the
unweighted operator norm is required.

## 2. Proof

Define operators into \(L_2(\Omega,\nu)\) by

$$
Ue_A=\sqrt{p_A}\,f_A,
\qquad
Ve_B=\sqrt{q_B}\,g_B.
$$

The weighted moment matrix is \(U^*V\).  Since all features have modulus
one,

$$
\|U\|_{\rm HS}^2=P,\qquad
\|V\|_{\rm HS}^2=Q.
$$

Schatten Hölder immediately gives

$$
\|U^*V\|_1
\le\|U^*\|_2\|V\|_2
=\sqrt{PQ}.
$$

This also explains why a large unweighted feature Gram norm can be
irrelevant to a passive probe: the same feature multiplicity is charged
inside the diagonal probability mass.

## 3. Exact-plant consequences

For one signed-permutation or quadratic-bent pair \((X,Y)\), take

$$
f_A=\prod_{i\in A}X_i,\qquad
g_B=\prod_{j\in B}Y_j.
$$

Equation (1.1) applies simultaneously to every \(M_{a,b}\), including all
odd records and even-pair decorations through total degree twelve.  It
strictly subsumes the endpoint calculation in
notes/quadratic_endpoint_weighted_contraction.md.

The pure-sector bosonic-compound norms and the quadratic affine-flat
collisions remain useful diagnostics of unweighted coherence, but neither
can introduce a factor after the correct physical diagonal normalization
is retained on one link.

## 4. Remaining nontrivial interface

The four-block exact plant has moment tensor

$$
m(A,B,C,D)=M_1(A,B)M_2(B,C)M_3(C,D).
$$

The middle support labels \(B,C\) occur in two neighboring link kernels
and in one physical product block.  A passive probe law can correlate its
four block supports, and an adaptive child law can depend on the complete
parent outcome.  Therefore three separate applications of (1.1) do not
automatically multiply: doing so would duplicate the middle diagonal
weights.

The exact remaining theorem is a three-link weighted composition
inequality that assigns each middle physical mass once and survives
outcome-selected preparations.  This is narrower than the earlier
representation-sector program: all single-link spectra are already
controlled.

Reproduction: tests/weighted_link_gram_contraction.py.
