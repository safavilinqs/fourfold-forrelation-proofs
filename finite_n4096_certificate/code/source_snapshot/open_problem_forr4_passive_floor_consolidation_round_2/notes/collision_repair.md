# Collision-aware repair of fiberwise frame packing

Date: 2026-07-14

Status: proved local repair, integrated with the combined insertion/Bessel
identity, and accepted as part of the audited global contraction.  The
separate RT-002 frontier issue was later resolved by the operator-valued
invariant.

## 1. Collision-aware lemma

Retain the complete frame and maps

$$
U_{y,r}(x)=c_y(\tau(r,x)),\qquad
V_{y,s}(z)=\overline{c_y(\upsilon(s,z))}.
$$

Do not assume injectivity. Define the maximum within-base fiber sizes

$$
\kappa_\tau=\max_{r,\omega}|\{x:\tau(r,x)=\omega\}|,
\qquad
\kappa_\upsilon=\max_{s,\omega}|\{z:\upsilon(s,z)=\omega\}|.
$$

For PSD families $P_s$ and $Q_r$, put

$$
A_{\rm diag}=
\sum_{r,s,x}(P_s)_{xx}q(\tau(r,x)),
\qquad
B_{\rm diag}=
\sum_{r,s,z}(Q_r)_{zz}q(\upsilon(s,z)).
$$

If

$$
W_{y,r,s}W_{y,r,s}^*\preceq P_s\otimes Q_r,
$$

then

$$
\boxed{
\sum_{y,r,s}|\langle U_{y,r}\otimes V_{y,s},W_{y,r,s}\rangle|
\le\sqrt{\kappa_\tau\kappa_\upsilon A_{\rm diag}B_{\rm diag}}.}
$$

The same statement holds for Hilbert-valued outputs.

## 2. Proof

The operator majorant and Cauchy--Schwarz reduce the left side to the product of

$$
\left(\sum_{y,r,s}\langle U_{y,r},P_sU_{y,r}\rangle\right)^{1/2}
$$

and its bra analogue. Frame completeness turns the first squared quantity into

$$
\sum_{r,s,\omega}q(\omega)
\sum_{x,x'\in F_{r,\omega}}(P_s)_{xx'},
\qquad
F_{r,\omega}=\{x:\tau(r,x)=\omega\}.
$$

For any PSD matrix $P$ and any set $F$,

$$
\sum_{x,x'\in F}P_{xx'}
=\langle\mathbf1_F,P\mathbf1_F\rangle
\le |F|\operatorname{Tr}P[F,F]
\le\kappa_\tau\sum_{x\in F}P_{xx}.
$$

Therefore the ket square sum is at most $\kappa_\tau A_{\rm diag}$. The bra side is at most $\kappa_\upsilon B_{\rm diag}$. This proves the scalar bound, and testing a Hilbert-valued sum against a unit auxiliary vector proves the vector form.

RT-001 attains this collision factor: $P$ is the two-dimensional all-ones matrix and $|F|=2$.

## 3. Consequence for fixed marked diagrams

For a node receiving $a$ labeled ket marks and $b$ labeled bra marks, hold the base configuration and insertion-slot choice fixed. A physical parity support determines the marked coordinate tuple up to permutations of equal-block marked vertices. Hence

$$
\kappa_\tau\le a!,\qquad\kappa_\upsilon\le b!,
$$

and in fact the product of within-block factorials is enough. Since every diagram has at most twelve vertices, this is an absolute diagram constant.

Crucially, the insertion-slot label is part of $r$ or $s$ and is not included in $\kappa$. Its sum remains inside $A_{\rm diag}$ or $B_{\rm diag}$. For a fixed sliced graph covariance, the diagonal is constant in the selected graph coordinate. Summing its opposite-entry graph fibers gives the round-one identity

$$
\sum_{s_C}\operatorname{diag}(M_{C,s_C}M_{C,s_C}^*)
=N^{|B_C|-e_C}\mathbf1.
$$

The remaining $q$-weighted insertion sum is precisely the falling-factorial Bessel mass

$$
\sum_Sq(S)(|S|)_a\le t^a
$$

on the ket and at most $t^b$ on the bra. Thus the repaired local estimate costs

$$
\sqrt{a!b!}\,t^{(a+b)/2}
$$

times the graph suppression. It creates no power of $N$, no second insertion-slot choice, and no change to the dose exponent.

## 4. Integration checks subsequently discharged

To instantiate this repair, the physical skeleton must explicitly separate:

1. base and insertion-slot indices, which are summed in the falling-factorial mass;
2. opposite-entry graph-coordinate fibers, which cancel against sliced diagonals; and
3. within-support ordering collisions, which pay only the factorial above.

This resolves the omitted injectivity condition locally.  The separate
adaptive-descendant PSD-majorant question is resolved in
notes/operator_frontier_invariant.md and audited in
AUDIT_REPAIRED_CONTRACTION.md.
