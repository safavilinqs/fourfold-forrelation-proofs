# Reverse-tree counterexample and gap log

Date: 2026-07-14

## RT-001 — Collision defect in the fiberwise lemma statement

Status: confirmed local counterexample to the round-one statement;
repaired by collision-aware PSD packing and the combined
collision/insertion identity.  The repaired theorem does not use the
false diagonal-only lemma.

Take one frame coordinate and one frame outcome, with $c=1$ and $q=1$. Let $X=\{0,1\}$ and let both selected coordinates collide under $\tau$:

$$
\tau(0,0)=\tau(0,1)=\omega.
$$

Let the other side be one-dimensional, take

$$
P=\begin{pmatrix}1&1\\1&1\end{pmatrix},\qquad Q=(1),\qquad W=(1,1)^\mathsf T.
$$

Then $WW^*=P\otimes Q$, so the claimed PSD-majorant hypothesis holds. The left side of the round-one Lemma 4.1 is $2$, whereas its diagonal-only right side is $\sqrt2$. The exact collision-aware block sum is $4$ and gives the correct bound $2$.

The earlier translated-frame Lemma 1.1 in Note 18 included fixed-base
injectivity, but the repaired bilateral statement in audit 04 and the
consolidated proof omitted it.  Round two supplied the missing combined
formula: within-base collisions cost only the appropriate factorial, and
that fiber is part of the same falling-factorial Bessel sum rather than a
second insertion charge.

Reproduction: tests/fiberwise_bilateral_stress.py.

## RT-002 — Reverse invariant is not yet a checkable proposition

Status: resolved by the typed operator-valued frontier invariant and the
one-node Hilbert-valued update.  The literal unit-vector formulation was
false for open identity wires, as recorded in RT-007.

The round-one text says the residual selected vector has preimage norm at most one under the remaining sliced graph maps. It does not define the residual spaces or write the marked-node update that returns the same range/preimage condition at the parent. Consequently the PSD hypothesis of the local bilateral lemma cannot currently be checked from the written induction.

Resolution: notes/operator_frontier_invariant.md curries open inputs into
a boundary space, applies the one-node Hilbert update uniformly on unit
boundary vectors, and then takes the operator supremum.  The global
dichotomy removes the mixed Hilbert/projective regime.

## RT-003 — Hilbert-valued control does not preserve singleton projective mass

Status: confirmed norm-level counterexample to a round-one proof
inference.  The global dichotomy repair avoids the inference entirely and
passed the physical-index audit; no physical counterprotocol was found.

The local bilateral lemma returns an auxiliary output with controlled Hilbert norm. Round one then says singleton-projective components can equivalently be absorbed into that auxiliary output and later bounded by their injective graph norm. This inference needs projective control, not Hilbert control.

For a normalized $N\times N$ Hadamard edge $H$, $\|H\|_{\rm op}=1$. Nevertheless

$$
z=\operatorname{vec}(H)/\sqrt N
$$

has Hilbert norm one and projective norm $\sqrt N$, and $\langle\operatorname{vec}H,z\rangle=\sqrt N$. Thus arbitrary Hilbert ancillary absorption can lose a growing dimension factor even for the simplest nonspanning singleton component.

Resolution: the exhaustive global dichotomy uses a projective invariant
only when every component is singleton in every physical entry, and a
pure Hilbert/operator invariant after assigning every component
otherwise.  It never converts the displayed Hilbert auxiliary into a
projective one.

Reproduction: tests/singleton_ancilla_norm_gap.py and notes/reverse_step_lemma.md.

## RT-004 — Entrywise signed-permutation recording budget

Status: falsified as a realistic-size theorem target.

The first exact-plant program proposed bounding transcript distance by an absolute constant times the product of three hidden-label match probabilities. For one occurrence per block this predicts $2/q^3$.

A single passive dose-two batch using the adjacent-pair split of the minimal chain attains TV $1/N=1/q^2$, while either crossing-pair split attains TV $1/\sqrt N=1/q$. The crossing flattening has nuclear norm $N^{3/2}$, and uniform half-mass on its two pair sectors achieves the stated value. The ratio to the entrywise record budget is $q^2/2$, so no dimension-independent recording constant can make that target true.

This does not defeat the exact plant: at $N=1024$ the sharp explicit value is $1/32$. It shows both that the replacement bound must retain spectral norms and that a uniform $C(D)/N$ target is impossible.

Reproduction: searches/minimal_chain_recording_counterexample.py.

## RT-005 — Termwise adaptive block-occupation ledger

Status: rejected as a finite-size proof method.

Exact block-occupation masses and cut ranks prove the one-batch minimal-chain coefficient is below the $N=1024$ threshold. They do not help if every assignment of four marked blocks to adaptive ket/bra entries is optimized and summed separately.

For the dose partition $(1,1,1,1,1,1)$, exact enumeration gives coefficient $8730$ relative to $N^{-1/2}$, versus the required $32/3$. The overshoot is about $818$. This is an upper-ledger failure, not a passive counterprotocol.

Required resolution: combine temporal placements before absolute values using a square-function, martingale, or exact tester norm. Any proof retaining the termwise $\ell_1$ marked-time sum cannot reach realistic sizes.

Reproduction: searches/rejected_termwise_adaptive_ledger.py.

## RT-006 — Distinct-label masks break literal vertical multiplicativity

Status: confirmed proof-level counterexample; repaired with a bounded-mass
character expansion.

Fourier supports contain distinct coordinates within each block. The
corresponding all-distinct mask couples disconnected graph components, so
the masked tensor is not their vertical product. For two four-layer chains
at $N=4$, with relative physical-entry placement $(0,1,3,2)$, an exact
product-vector witness has value

$$
{15\over32\sqrt3}>{1\over4},
$$

where $1/4=1/N$ is the product of the two unmasked component norms.

For Sylvester labels, expand every inequality indicator in Walsh
characters. Each pair factor has coefficient mass at most two, and every
expanded term is a true vertical product modified only by local diagonal
unitaries. The total loss is therefore at most the diagram constant

$$
2^{\sum_r\binom{|V_r|}{2}},
$$

with no power of $N$. The $N^{-1/2}$ contraction survives, while any claim
of exact multiplicativity for the masked tensor must be rejected.

Reproduction: `tests/all_singleton_masked_graph_norm.py` and
`notes/distinct_label_mask_repair.md`.

## RT-007 — Open identity wires are not unit Hilbert vectors

Status: confirmed invariant-level defect; repaired by an operator-valued
frontier.

At an intermediate reverse frontier, unresolved graph coordinates carry
identity wires. Although $\|I_N\|_{\rm op}=1$, vectorizing such a wire gives
Hilbert--Schmidt norm $\sqrt N$. Thus the phrase "the residual is the image
of a unit Hilbert vector" is false before all open inputs have been fixed.

Curry those inputs into a boundary space $J_h$ and maintain a contraction
$Z_h:J_h\to H_{\rm pre,h}\otimes K_h$. For each unit $\xi\in J_h$,
$Z_h\xi$ is a unit Hilbert vector, so the existing reshape and
Hilbert-valued local lemma apply. The estimate is uniform in $\xi$; taking
the supremum preserves the operator norm of the returned frontier map.

Reproduction: `tests/open_frontier_operator_stress.py` and
`notes/operator_frontier_invariant.md`.
