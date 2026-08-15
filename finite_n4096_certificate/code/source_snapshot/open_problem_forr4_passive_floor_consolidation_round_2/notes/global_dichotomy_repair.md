# Global dichotomy repair for mixed component norms

Date: 2026-07-14

Status: proved and integrated into the repaired contraction.  The
independent line-by-line audit passed; see
AUDIT_REPAIRED_CONTRACTION.md and CONFIDENCE_REPORT.md.

## 1. The problematic split is unnecessary

Round one simultaneously treated components with a multi-vertex amplitude entry by sliced Hilbert range maps and components with only singleton entries by projective/injective norms. That creates the Hilbert-versus-projective interface exposed by RT-003.

For the theorem, only one factor $N^{-1/2}$ is needed. This permits a global dichotomy that never mixes the two norm regimes.

## 2. Case A: every component is singleton at every amplitude entry

Assume

$$
|C\cap E|\le1
$$

for every graph component $C$ and every physical ket or bra amplitude entry $E$.

Group the vertex Hilbert spaces by physical amplitude entry. Each component tensor has one-dimensional dummy parties where it is absent. Its injective norm under this grouping is at most its natural-cut operator norm, hence at most $N^{-1/2}$ if it spans all four layers and at most one otherwise.

Vertical tensor-product multiplicativity therefore gives

$$
\left\|\mathop{\boxtimes}_C T_C\right\|_\varepsilon
\le N^{-s_0/2},
$$

where $s_0$ is the number of four-layer components. At least one component spans all four layers, so $s_0\ge1$.

At each marked physical node, its ket/bra frame atom is one Hilbert vector in each of the two grouped entry spaces. Complete-frame Cauchy--Schwarz bounds the sum of the products of their Hilbert norms by the single joint insertion/Bessel mass. Reverse stochastic summation therefore produces a projective decomposition over physical amplitude entries, not separately over graph components. Duality with the displayed injective norm gives $N^{-1/2}$.

Arbitrary cross-component entanglement inside one physical entry is already allowed because that entire grouped entry is one Hilbert party.

## 3. Case B: at least one component has a multi-vertex entry

For every component $C$, choose one amplitude entry $E(C)$ maximizing

$$
k_C=|C\cap E(C)|.
$$

Unlike round one, assign components with $k_C=1$ as well as those with $k_C\ge2$. Put $A_C=C\cap E(C)$ and $B_C=C\setminus A_C$. For cross-entry fibers $F_C\subseteq B_C$, use the sliced graph map $M_{C,s_C}$.

The exact diagonal sum is unchanged:

$$
\sum_{s_C}\operatorname{diag}(M_{C,s_C}M_{C,s_C}^*)
=N^{|B_C|-e_C}\mathbf1.
$$

Because $e_C\ge v_C-1$,

$$
|B_C|-e_C=v_C-k_C-e_C\le1-k_C.
$$

Thus a weak assignment with $k_C=1$ costs at most one, while every strong assignment with $k_C\ge2$ contributes at least $N^{-1/2}$. By the hypothesis of Case B, at least one strong assignment exists.

Tensor the sliced maps of all assigned components. There are now no separately projective graph components. Every unselected graph coordinate is either:

1. an input coordinate of one of these sliced Hilbert maps;
2. a fixed cross-entry fiber; or
3. a free Hilbert coordinate emitted after its assigned map has been eliminated.

The reverse invariant is purely Hilbert-valued. At a node with assigned component outputs in both entries, apply collision-aware bilateral packing. With outputs in only one entry, use its unilateral specialization. With no assigned outputs, complete-frame Bessel packing acts directly on the current free Hilbert coordinates and costs at most one after normalization.

The operator-majorant hypothesis follows from the reshape argument in notes/reverse_step_lemma.md. Because every component was assigned, no later step asks an ordinary Hilbert auxiliary output to satisfy a projective norm bound. The one strong component supplies $N^{-1/2}$; all weak components preserve the norm.

## 4. Dose and collision ledger

The collision-aware lemma keeps base and insertion-slot labels inside its diagonal square sums. Within-support ordering collisions cost at most a diagram factorial. Opposite-entry graph fibers cancel against the sliced diagonal sum. The remaining insertion sums are

$$
\sum_Sq(S)(|S|)_a\le t^a,
\qquad
\sum_Sq(S)(|S|)_b\le t^b.
$$

Therefore the local cost remains an absolute diagram constant times $t^{(a+b)/2}$. The global marked-time sum remains $C_G(1+D)^{v(G)}$ rather than acquiring a second slot count.

## 5. Consequence, after the physical index audit

The two cases are exhaustive:

- If every component is singleton in every entry, Case A gives $N^{-1/2}$ from a four-layer component.
- Otherwise Case B gives $N^{-1/2}$ from a strong assignment, regardless of which layers it spans.

This restores the finite-skeleton contraction without the invalid
ancillary-norm inference in RT-003.

## 6. Independent checks used for acceptance

The acceptance audit verified:

1. Case A frame atoms really form a projective decomposition over grouped physical entries under adaptive reverse summation.
2. In Case B, the pure Hilbert preimage invariant covers nodes with zero, one, or two assigned entries and arbitrary free coordinates from previously eliminated maps.
3. Cross-entry fibers of weak $k_C=1$ assignments obey the same diagonal cancellation.
4. The all-assigned construction does not reuse a complete physical frame or a marked insertion fiber.
5. Small physical optimizations fail to exploit the change of component orientation across temporal frontiers.
