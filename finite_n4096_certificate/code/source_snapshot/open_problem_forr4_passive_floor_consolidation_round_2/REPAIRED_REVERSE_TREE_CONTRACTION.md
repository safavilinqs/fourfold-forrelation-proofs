# Repaired reverse-tree contraction

Date: 2026-07-14

Status: validated round-two theorem interface after the line-item audit in
AUDIT_REPAIRED_CONTRACTION.md.  This proof replaces the mixed-component
step from round one; it does not reuse the invalid Hilbert-to-projective
ancillary inference.  It must be integrated together with the collision,
distinct-label-mask, operator-frontier, and interpolation repairs.

## Theorem

Fix a marked passive skeleton, a ket/bra split, its insertion fibers, and a layered Hadamard multigraph $G$. Let $\mathcal E$ be the finite set of physical amplitude entries, one ket and one bra entry at every marked node. Each graph vertex is placed in exactly one $E\in\mathcal E$.

After normalizing by the one joint insertion/Bessel mass $\mathcal B_\phi$, the graph contraction satisfies

$$
|\mathfrak C_{G,\phi}|\le C_GN^{-1/2}\mathcal B_\phi,
$$

where $C_G$ depends only on the finite diagram and not on $N$, depth, outcome width, or dose.

## Facts used

1. Every graph component $C$ is nontrivial, so $e_C\ge v_C-1$.
2. At least one component spans all four layers.
3. For $A_C\subseteq C$, $B_C=C\setminus A_C$, and a fixed opposite-entry fiber $s_C$ on $F_C\subseteq B_C$,

   $$
   \operatorname{diag}(M_{C,s_C}M_{C,s_C}^*)
   =N^{|B_C|-|F_C|-e_C}\mathbf1.
   $$

   Summing all $N^{|F_C|}$ fibers gives $N^{|B_C|-e_C}\mathbf1$.
4. The collision-aware bilateral frame lemma and its unilateral special case hold in Hilbert-valued form.
5. A complete marked frame has joint insertion masses at most $t^a$ and $t^b$ for $a$ ket and $b$ bra marks. Within-support ordering collisions cost at most $a!b!$ and are diagram constants.
6. Unmarked histories are stochastic contractions.

Before either case below, impose the distinct-coordinate condition on the
marked Fourier labels. For Sylvester labels,

$$
{\bf1}\{x\ne y\}
=\left(1-{1\over N}\right)
-{1\over N}\sum_{\chi\ne1}\chi(x)\overline{\chi(y)}.
$$

Expanding all same-layer distinctness factors has coefficient mass at most
$2^{P_G}$, where $P_G=\sum_r\binom{|V_r|}{2}$. Each summand only applies
local diagonal character unitaries to the component tensors. It therefore
preserves all injective norms and sliced covariance diagonals below. We
prove the estimate termwise and absorb $2^{P_G}$ into $C_G$. This step is
necessary: the masked tensor is not literally a vertical product (see
`notes/distinct_label_mask_repair.md`).

## Case I: all components are singleton in every physical entry

Assume $|C\cap E|\le1$ for every $C$ and $E$.

View each component tensor as a tensor over the common party list $\mathcal E$, inserting one-dimensional dummy parties. Because a component contributes at most one vertex to a party, its injective norm over $\mathcal E$ is bounded by its natural layer-cut operator norm:

$$
\|T_C\|_{\varepsilon,\mathcal E}
\le
\begin{cases}
N^{-1/2},&C\text{ spans four layers},\\
1,&\text{otherwise}.
\end{cases}
$$

Termwise vertical tensor-product multiplicativity gives

$$
\left\|\mathop{\boxtimes}_CT_C\right\|_{\varepsilon,\mathcal E}
\le N^{-s_0/2}\le N^{-1/2}.
$$

After summing the mask expansion, the left side acquires only the diagram
constant $2^{P_G}$.

At a marked node, hold every base and insertion fiber fixed. Its ket and bra coefficients are Hilbert vectors in the two grouped physical-entry spaces. Complete-frame Cauchy--Schwarz bounds the sum of the products of their norms by the local joint Bessel mass. This remains true when the vectors are joint, entangled coordinates from several graph components.

Reverse outcome summation and projective-norm subadditivity therefore express the full frame skeleton as a projective tensor over $\mathcal E$ of mass at most $C_G\mathcal B_\phi$. Pairing it with the graph tensor and using injective/projective duality proves the theorem in Case I.

No Hilbert auxiliary is later reinterpreted as a projective tensor: the projective norm is maintained from the first reverse step to the final contraction.

## Case II: some component has a multi-vertex physical entry

For every component $C$, choose an entry $E(C)$ maximizing

$$
k_C=|C\cap E(C)|\ge1,
$$

and assign $C$ to that entry. Thus every component is assigned, including components with $k_C=1$. Put $A_C=C\cap E(C)$ and $B_C=C\setminus A_C$.

The fiber-summed range diagonal obeys

$$
N^{|B_C|-e_C}
=N^{v_C-k_C-e_C}
\le N^{1-k_C}.
$$

Hence a weak assignment $k_C=1$ costs at most one, while a strong assignment $k_C\ge2$ supplies at least $N^{-1/2}$ after the frame square root. Case II contains at least one strong assignment.

It remains to justify simultaneous reverse summation. In a fixed
character-expansion term, orient every component tensor as its sliced map
from $B_C$ to $A_C$ and tensor all these maps. The local character unitaries
do not change the sliced diagonal identity. At a reverse frontier, after
fixing cross-entry fibers, maintain the following operator-valued invariant:

> Curry every still-open ancestor coordinate into an input space $J_h$.
> After division by accumulated Bessel mass, the residual factors through
> the tensor product of the currently exposed sliced maps as
> $R_h=L_hZ_h$, where
> $Z_h:J_h\to H_{\rm pre,h}\otimes K_h$ has operator norm at most one.
> Coordinates emitted by maps already removed and unresolved identity
> wires are retained in $K_h$ or $J_h$; they are not vectorized and charged
> in Hilbert--Schmidt norm.

The invariant starts at terminal scalars. Suppose a marked node is removed. Combine all assigned outputs in its ket entry into $L_{\rm ket}$ and all assigned outputs in its bra entry into $L_{\rm bra}$. All other coordinates present in those entries are free Hilbert coordinates or inputs of components assigned elsewhere.

Fix a unit vector $\xi\in J_h$. Its preimage $Z_h\xi$ has Hilbert norm at
most one. Reshape this vector across the current selected coordinates and
the retained auxiliary coordinates. The resulting operator has operator
norm at most its Hilbert--Schmidt norm, hence at most one. Therefore the
current selected operator $W_\xi$ satisfies

$$
W_\xi W_\xi^*\preceq
(L_{\rm ket}L_{\rm ket}^*)\otimes
(L_{\rm bra}L_{\rm bra}^*),
$$

with an identity factor on any free selected coordinates. Inputs of a current component that occur in the opposite physical entry are fixed before this majorant is formed and are summed as graph fibers afterward.

The Hilbert-valued collision-aware bilateral lemma removes the complete
outcome and returns a vector in the retained auxiliary Hilbert space with
the required bound, uniformly for every unit $\xi$. Taking the supremum
over $\xi$ proves the same operator-norm bound for the returned frontier
map. If only one entry has assigned outputs, use the unilateral version. If
neither does, the same lemma with identity majorants is ordinary
complete-frame Bessel contraction. Thus every node type preserves the
operator invariant and every physical frame is used once. The identity-wire
point is detailed in notes/operator_frontier_invariant.md.

Tensoring the exact fiber diagonals over components gives the graph factor

$$
\prod_C N^{(v_C-k_C-e_C)/2}
\le N^{\frac12\sum_C(1-k_C)}
\le N^{-1/2},
$$

because all $k_C\ge1$ and at least one $k_C\ge2$.

This proves the theorem in Case II. The proof never creates a separately projective ancillary component, so RT-003 does not apply.

## Restoring the dose ledger

At a node with $a+b$ marks, collision-aware packing and joint completeness cost at most

$$
C_{a,b}t^{(a+b)/2},
$$

where $C_{a,b}$ is a factorial depending only on the diagram. Opposite-entry graph fibers are already included in the sliced diagonal sum. Insertion-slot fibers remain in the same falling-factorial Bessel sum, so there is no second placement selection.

Explicitly, for one entry the base/ordered-mark map satisfies

$$
\sum_{R,x}q(R\cup\{x_1,\ldots,x_a\})
=\sum_Sq(S)(|S|)_a\le t^a,
$$

while every fixed $(R,S)$ fiber has size at most $a!$. Thus the factorial
collision loss and the insertion mass occur in the same sum, as detailed in
notes/combined_collision_bessel_identity.md.

Summing ket/bra splits and marked times on a branch gives

$$
|\mathfrak C_G(\mathcal T)|
\le C_G(1+D)^{v(G)}N^{-1/2}.
$$

For $v(G)\le12$, this is the reverse-tree input required for the transcript bound $C(1+D)^{12}/\sqrt N$.
