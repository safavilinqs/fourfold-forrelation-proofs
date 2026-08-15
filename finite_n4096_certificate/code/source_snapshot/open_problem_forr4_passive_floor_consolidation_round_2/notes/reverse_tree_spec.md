# Reverse-tree contraction: typed audit specification

Date: 2026-07-14

Status: frozen audit interface.  It is not itself a proof, but every
acceptance row was discharged by the repaired contraction and the
line-item audit.

## 1. Claim under audit

For each fixed marked skeleton, fixed ket/bra split, and layered Hadamard multigraph $G$, the normalized contraction must satisfy

$$
|\mathfrak C_{G,\phi}|\le N^{-1/2}\mathcal B_\phi,
$$

uniformly over the adaptive tree, its outcome alphabets, all passive complete frames, all base supports, and all terminal signs. Here $\mathcal B_\phi$ is the single joint insertion/Bessel mass that is later summed in the dose ledger.

The statement must not hide a factor depending on $N$, tree depth, outcome width, the number of fibers, or a second selection of insertion slots.

## 2. Finite data that must be explicit

A mechanically checkable instance consists of:

1. A finite rooted history tree. Each marked node $h$ has outcome set $Y_h$ and may choose descendants as a function of its outcome.
2. A complete rank-one frame $c_{h,y}\in\mathbb C^{\Omega_h}$ satisfying

   $$\sum_{y\in Y_h}c_{h,y}c_{h,y}^*=D_{q_h},\qquad q_h\ge0,\qquad\operatorname{Tr}D_{q_h}=1.$$

3. A fixed labeled-slot lift of every marked derivative. The ket and bra insertion maps, including their base labels, must be written as actual maps into $\Omega_h$.
4. A layered Hadamard multigraph $G$ with each vertex assigned to exactly one marked ket or bra entry.
5. For every connected component assigned to an entry, a selected set $A_C$ with $|A_C|\ge2$, the complement $B_C$, the subset $F_C\subseteq B_C$ occurring in the opposite entry of the same frame, and every sliced graph map $M_{C,s_C}$.
6. The remaining singleton-projective components, their four-layer status, and the precise physical party grouping used for their injective norm.
7. The terminal functional, bounded in modulus by one, and every unmarked stochastic kernel between marked nodes.

Without these types, the phrase selected joint vector does not determine which norm or range condition is being asserted.

## 3. Correct local fiberwise lemma

Let

$$
U_{y,r}(x)=c_y(\tau(r,x)),\qquad
V_{y,s}(z)=\overline{c_y(\upsilon(s,z))}.
$$

The diagonal formulas used in the round-one proof are valid if $x\mapsto\tau(r,x)$ is injective for every fixed $r$, and $z\mapsto\upsilon(s,z)$ is injective for every fixed $s$. Under those hypotheses, define

$$
\Lambda=\max_\omega\sum_{r,x:\tau(r,x)=\omega}\sum_s(P_s)_{xx},
$$

and define $\mathrm M$ symmetrically. If

$$
W_{y,r,s}W_{y,r,s}^*\preceq P_s\otimes Q_r,
$$

then

$$
\sum_{y,r,s}|\langle U_{y,r}\otimes V_{y,s},W_{y,r,s}\rangle|
\le\sqrt{\Lambda\mathrm M}.
$$

The same bound is Hilbert-valued if $W_{y,r,s}$ is viewed as a map from an auxiliary output space and satisfies $W_{y,r,s}W_{y,r,s}^*\preceq P_s\otimes Q_r$.

The injectivity hypotheses are essential. If collisions are retained, the exact first square sum contains

$$
\sum_{r,s}\sum_{x,x':\tau(r,x)=\tau(r,x')}
(P_s)_{x x'}q(\tau(r,x)),
$$

not only the diagonal of $P_s$. A safe collision-aware version may either use this exact positive block sum or pay an explicit collision multiplicity. That multiplicity must then be reconciled with the joint insertion mass rather than charged a second time.

## 4. Exact reverse invariant

At every reverse frontier $h$ and for every fixed value $f$ of the already
exposed graph fibers, curry the open ancestor coordinates into an input
space $J_{h,f}$. The residual tensor must admit

$$
R_{h,f}=A_{h,f}Z_{h,f},
\qquad
Z_{h,f}:J_{h,f}\to H_{\rm pre,h,f}\otimes K_{h,f},
\qquad
\|Z_{h,f}\|_{\rm op}\le b_{h,f}.
$$

Here $A_{h,f}$ is the tensor product of precisely the sliced graph maps that
remain assigned below that frontier, with identities on free current
coordinates. The scalar $b_{h,f}$ is the unused portion of the one joint
insertion/Bessel mass. Open identity wires are measured in operator norm,
not Hilbert--Schmidt norm.

The global dichotomy is applied before this invariant: the all-singleton
case stays projective from start to finish, while the all-assigned case uses
this operator-valued invariant and has no separately projective ancillary
component. Distinct-label masks are Walsh-expanded before either case.

Equivalently, whenever one current ket slice and one current bra slice are exposed, their selected residual must obey

$$
W_{y,r,s}W_{y,r,s}^*\preceq P_s\otimes Q_r,
$$

with $P_s$ and $Q_r$ the corresponding sliced range covariances. These
majorants must be independent of the current outcome $y$. To obtain them,
fix a unit $\xi\in J_{h,f}$, reshape the Hilbert vector $Z_{h,f}\xi$, and
then apply the Hilbert-valued local lemma.

The induction step must prove more than a scalar estimate. Its Hilbert bound
must hold uniformly for every unit $\xi$; taking the supremum then returns
the operator-valued parent invariant with the updated
$b_{\operatorname{parent}(h),f'}$.

## 5. Obligations at one marked elimination

For each marked node, a complete proof must verify all of the following in one formula:

1. Labeled-slot completeness gives the joint ket/bra square masses $t^a$ and $t^b$ before any graph-coordinate estimate.
2. Cross-entry graph coordinates are fixed as fibers before forming range covariances.
3. The maps from the remaining selected coordinates into the labeled frame basis are injective. If the parity quotient is used instead, its collision operator or multiplicity is explicit.
4. After a unit open-boundary input is fixed, the descendant residual is the image of a unit Hilbert preimage under the tensor product of the current sliced maps. Uniformity in that input returns the operator-valued parent invariant.
5. Components assigned to the ket and bra entries are disjoint tensor factors even if their coordinate values coincide numerically.
6. The global dichotomy has selected one norm regime: projective in the all-singleton case or operator/Hilbert in the all-assigned case. They are never mixed.
7. Hilbert-valued bilateral packing removes the physical outcome once and returns the parent preimage invariant.
8. The fibers summed here are the same insertion fibers counted by the dose ledger; no second slot or collision factor is introduced.

## 6. Current evidence classification

| Item | Current evidence | Audit status |
|---|---|---|
| Sliced graph diagonal identity | Exact algebra plus deterministic and randomized finite checks | Proved and independently reproduced |
| Unsliced vectorized majorant | Explicit rank-$N$ counterexample | Disproved; must never be used |
| Fiberwise scalar lemma with injective maps | Direct Cauchy--Schwarz proof | Proved |
| Fiberwise lemma as stated for arbitrary maps | Two-coordinate collision counterexample | Disproved as stated |
| Hilbert-valued local lemma under an operator majorant | Duality proof; adversarial regression | Proved and regression-backed |
| PSD majorant for every residual created by the adaptive reverse induction | Operator-valued frontier invariant; fix a unit input and reshape | Proved after RT-007 repair |
| Compatibility of collision handling with the nonduplicated $t^{k/2}$ insertion mass | Exact base/ordered-mark identity and rational enumeration | Proved and reproduced |
| Singleton ancillary preservation under mixed norm regimes | Global dichotomy avoids every mixed regime | Obsolete obligation; RT-003 avoided |
| Distinct-label mask before component tensoring | Bounded-mass Walsh expansion; exact omission witness | Proved after RT-006 repair |

## 7. Acceptance gate

The repaired interface now has no row classified as unresolved. Acceptance
still depends on preserving the Walsh expansion, global dichotomy,
operator-valued frontier, and combined collision/Bessel identity together;
the regressions are falsification evidence, not substitutes for those
analytic steps.
