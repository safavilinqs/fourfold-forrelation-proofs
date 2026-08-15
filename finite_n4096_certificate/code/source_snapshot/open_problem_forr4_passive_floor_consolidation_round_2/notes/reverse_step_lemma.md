# Abstract reverse step and the repaired ancillary interface

Date: 2026-07-14

Status: the Hilbert preimage step is proved. RT-003 shows why a mixed
Hilbert/projective induction fails; the global dichotomy avoids that mixed
regime, and notes/operator_frontier_invariant.md supplies the required
operator-valued formulation for unresolved identity wires.

## 1. Hilbert-valued reverse step

Let $E=X\otimes Z$ and let $K$ be an auxiliary Hilbert space. For every $(y,r,s)$ let

$$
W_{y,r,s}:K\longrightarrow E
$$

satisfy

$$
W_{y,r,s}W_{y,r,s}^*\preceq P_s\otimes Q_r.
$$

For coefficients $|\alpha_{y,r,s}|\le1$, define the reverse update

$$
z=\sum_{y,r,s}\alpha_{y,r,s}W_{y,r,s}^*(U_{y,r}\otimes V_{y,s})\in K.
$$

Then the injective fiberwise lemma gives $\|z\|_2\le\sqrt{\Lambda\mathrm M}$, or the collision-aware version gives

$$
\|z\|_2\le
\sqrt{\kappa_\tau\kappa_\upsilon A_{\rm diag}B_{\rm diag}}.
$$

Indeed, for every unit $\xi\in K$,

$$
(W_{y,r,s}\xi)(W_{y,r,s}\xi)^*
\preceq W_{y,r,s}W_{y,r,s}^*
\preceq P_s\otimes Q_r.
$$

Testing $z$ against $\xi$, applying the scalar lemma, and taking the supremum over unit $\xi$ proves the claim. This is the precise Hilbert-duality argument that was implicit in round one.

If $A:K\to H$ is the tensor product of the graph maps remaining above the frontier, the physical output $Az$ stays in $\operatorname{ran}A$ with preimage norm bounded by the same quantity. Thus a marked outcome is used once and arbitrary dependence of descendants on $y$ is harmless.

## 2. Why the PSD majorant follows from a Hilbert preimage

Suppose the current residual before frame elimination is

$$
(L\otimes A)\xi,
$$

where $L$ is the tensor product of the current sliced graph maps, $A$ contains the remaining ancestor maps, and $\|\xi\|_2\le1$. Reshape $\xi$ into an operator $Z:K\to E_{\rm pre}$. Then

$$
\|Z\|_{\rm op}\le\|Z\|_{\rm HS}=\|\xi\|_2\le1.
$$

The selected operator is $W=LZ$, so

$$
WW^*=LZZ^*L^*\preceq LL^*.
$$

When the current ket and bra component families are disjoint, $L=L_{\rm ket}\otimes L_{\rm bra}$ and $LL^*=P_s\otimes Q_r$. This proves the majorant required by the local lemma, provided the global residual really has the displayed Hilbert preimage factorization.

Unmarked stochastic refinements only take convex combinations, so they do not change this conclusion or introduce outcome-width and depth factors.

## 3. The remaining gap: singleton ancillary factors use a different norm

Round one also carries components for which every physical entry contains at most one component vertex. Their graph tensors are bounded in Hilbert injective norm, and the corresponding frame skeleton must therefore be bounded in the dual Hilbert projective norm.

The Hilbert-valued reverse step above controls only the ordinary Hilbert norm of its auxiliary output. It does not control its projective norm across the singleton component parties. These norms can differ by a dimension factor.

For the normalized Hadamard edge $H$, viewed as a two-party tensor,

$$
\|H\|_\varepsilon=\|H\|_{\rm op}=1.
$$

But the auxiliary vector

$$
z={\operatorname{vec}H\over\sqrt N}
$$

has $\|z\|_2=1$ and projective norm $\sqrt N$. Its contraction is

$$
\langle\operatorname{vec}H,z\rangle=\sqrt N.
$$

Therefore the sentence that a singleton component can equivalently be absorbed as an arbitrary Hilbert auxiliary output is false without an additional structural invariant. The vertical injective-norm multiplicativity theorem does not by itself convert Hilbert control into projective control.

## 4. What closes the full induction

One of the following is needed:

1. A mixed Hilbert/projective reverse invariant proving that singleton-component coordinates remain in a projective decomposition of controlled mass while assigned-component coordinates use sliced Hilbert preimages.
2. A Banach-valued version of fiberwise packing whose output norm is exactly the required singleton projective norm.
3. A simultaneous global cut or component reassignment that treats every carried singleton component as an operator contraction throughout the temporal order, with no boundary at which only its Hilbert norm is known.
4. A counterexample showing that none of these is possible for the current hard instance.

The round-two global dichotomy implements the third route in a stronger
form: in the all-singleton case the projective invariant is kept from the
start, while in the other case every component is assigned and the
induction is operator/Hilbert-valued throughout. The two norm regimes are
never mixed. This avoids the counterexample rather than trying to bound its
projective norm.
