# Operator-valued reverse-frontier invariant

Date: 2026-07-14

Status: replacement for the insufficient literal unit-vector formulation of
the all-assigned reverse induction.

## 1. Why a vector invariant is too narrow

At an intermediate frontier, a component input may belong to an ancestor
frame that has not yet been removed. Its coordinate is an open identity
wire. On an \(N\)-dimensional wire,

$$
\|I_N\|_{\rm op}=1,
\qquad
\|I_N\|_{\rm HS}=\sqrt N.
$$

Thus one cannot vectorize the entire frontier and claim Hilbert norm one.
Doing so would introduce a spurious dimension factor. The Hilbert-vector
argument is valid only after a unit vector on the open input boundary has
been fixed.

## 2. Correct invariant

At a reverse frontier, curry all still-open ancestor coordinates into an
input Hilbert space \(J_h\), and retain the unremoved auxiliary coordinates
in \(K_h\). After dividing by the accumulated Bessel mass, the residual has
a factorization

$$
R_h=L_h Z_h,
\qquad
Z_h:J_h\longrightarrow H_{\rm pre,h}\otimes K_h,
\qquad
\|Z_h\|_{\rm op}\le1.
\tag{2.1}
$$

Here \(L_h\) is the tensor product of the sliced component maps currently
exposed at the frontier, with identity maps on free selected coordinates.
Identity wires are part of \(Z_h\), where they cost operator norm one.

The scalar terminal boundary and the final root boundary are the special
case \(J_h=K_h=\mathbb C\).

## 3. One marked reverse step

Fix a unit \(\xi\in J_h\). Then

$$
z_\xi=Z_h\xi,
\qquad
\|z_\xi\|_2\le1.
$$

Reshape \(z_\xi\) across the current selected preimage coordinates and the
retained auxiliary coordinates. The resulting operator \(S_\xi\) obeys

$$
\|S_\xi\|_{\rm op}
\le\|S_\xi\|_{\rm HS}
=\|z_\xi\|_2
\le1.
$$

Hence the current selected operator satisfies the required covariance
majorant

$$
W_\xi W_\xi^*
\preceq
(L_{\rm ket}L_{\rm ket}^*)\otimes
(L_{\rm bra}L_{\rm bra}^*),
$$

with identities on free current coordinates. Apply the Hilbert-valued
collision-aware frame lemma. It returns a vector in \(K_h\) with norm at
most the local bound, uniformly for every unit \(\xi\).

The reverse update is linear (or conjugate-linear, depending on the fixed
ket/bra convention) in \(\xi\). Taking the supremum over unit \(\xi\)
therefore proves that the returned frontier map has operator norm at most
the same local bound. This is (2.1) at the parent.

The same argument covers nodes with two, one, or zero assigned output
entries. In the zero-output case both graph majorants are identities; open
identity wires remain in the operator boundary and are never measured in
Hilbert--Schmidt norm.

## 4. Adaptive outcomes and fibers

The child frontier map may depend arbitrarily on the current outcome. The
Hilbert-valued local lemma already permits outcome-dependent selected
tensors. Its bound is uniform after \(\xi\) is fixed, so the operator
supremum introduces no outcome-width factor.

Cross-entry graph coordinates are still fixed before the covariance is
formed and summed with the sliced diagonal afterward. The operator
boundary changes neither that fiber sum nor the insertion/Bessel ledger.

## 5. Audit consequence

The sentence "reshape the unit preimage tensor" is sound only when read as
"fix a unit open-boundary input, then reshape its image." The invariant is
operator-valued before that input is fixed. The regression
tests/open_frontier_operator_stress.py checks both the identity-wire norm
gap and the uniform operator-valued local update.
