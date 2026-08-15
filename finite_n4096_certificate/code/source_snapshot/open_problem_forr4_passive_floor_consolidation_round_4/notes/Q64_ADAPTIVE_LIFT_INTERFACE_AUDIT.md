# $q=64$ adaptive-lift interface audit

Date: 2026-07-18

Status: **SUPERSEDED — the missing global normalization is proved with multiplier one in `Q64_ADAPTIVE_TREE_FRONTIER_THEOREM.md`.**

This note records the interface audit that preceded the adaptive theorem. Its arithmetic cap, rejection of nodewise posterior reuse, and individual Schur-amplification lemma remain valid. The later direct-sum induction supplies the common normalized feature family that was missing here and closes the finite-size adaptive lower bound without changing the 888/888 registry.

## Decision

The adaptive lift cannot be obtained by applying the one-batch theorem independently after each outcome. Conditioning changes the hard-instance law, and a child posterior need not satisfy the original unconditioned moment bounds. Nor can all measurements simply be deferred: coherent control of later probes enlarges the passive protocol to an inter-batch coherent strategy and can admit the active behavior that the theorem is intended to separate.

The viable route is an operator-valued frontier argument. It must keep the complete outcome tree together until the physical occurrence kernels have been contracted. The missing issue is global normalization across all nodes and leaves, not amplification of an individual physical kernel.

## 1. Exact numerical gate

Write the certified complete total as

$$
U_1=P+E,
$$

where

$$
\begin{aligned}
P&\le 0.2666330066532943971677563163721884844,\\
E&\le 0.0022251284066307022549501658645609374,\\
U_1&\le 0.2688581350599250994227064822367494217.
\end{aligned}
$$

The outward acceptance threshold is

$$
T=\frac13-10^{-3}=0.3323333333333333333333333333333333333.
$$

The additive adaptive overhead must therefore be strictly below

$$
T-U_1=0.0634751982734082339106268510965839116\ldots.
$$

If the promise loss is paid once and an adaptive frontier multiplies only the occurrence/Perron term by $C_{\rm front}$, the sufficient gate is

$$
C_{\rm front}<\frac{T-E}{P}=1.2380620429185861374326799922207765185\ldots.
$$

If a theorem instead multiplies the entire one-batch total, it must satisfy

$$
C_{\rm total}<\frac{T}{U_1}=1.2360917896691517626640466898161393920\ldots.
$$

Thus a $6/5$ occurrence multiplier fits, while $5/4$, $\sqrt{2}$, two, and any depth-linear bound fail. The project needs a genuinely near-isometric frontier theorem, but it does not need literal zero overhead.

## 2. What the arbitrary-diagonal coefficient already proves

For a finite physical occurrence matrix $K=(K_{ij})$, define

$$
\Gamma(K)=\sup_{p,w}
\left\|
\operatorname{diag}(\sqrt p)K\operatorname{diag}(\sqrt w)
\right\|_1,
$$

where $p,w$ range over probability laws on the physical row and column supports. The finite-dimensional Hilbert factorization theorem gives

$$
\Gamma(K)=
\inf_{K_{ij}=\langle u_i,v_j\rangle}
\left(\max_i\|u_i\|\right)
\left(\max_j\|v_j\|\right).
$$

The easy direction follows by writing the weighted matrix as the product of the weighted row-feature and column-feature matrices and applying Schatten Hölder. The reverse direction is the dual finite-dimensional factorization theorem for the weighted trace norm. Consequently every accepted Round-4 statement of the form

$$
\left\|
\operatorname{diag}(\sqrt p)K\operatorname{diag}(\sqrt w)
\right\|_1
\le c\sqrt{\left(\sum p\right)\left(\sum w\right)}
$$

is exactly a Hilbert factorization bound $\Gamma(K)\le c$ on that physical support family.

This automatically survives one contractive outcome amplification. If

$$
G_{ij}=\langle a_i,b_j\rangle,
\qquad
\|a_i\|,\|b_j\|\le1,
$$

then

$$
(K\circ G)_{ij}
=\langle u_i\otimes a_i,v_j\otimes b_j\rangle,
$$

and hence

$$
\boxed{\Gamma(K\circ G)\le\Gamma(K).}
$$

The same proof covers $G_{ij}=\langle a_i,Zb_j\rangle$ for any contraction $Z$. Therefore an auxiliary outcome register, quantum side-information register, or leaf-sign contraction does not inflate an individual one-batch coefficient when it enters as one normalized cross-Gram factor. This conclusion applies uniformly to all 888 accepted physical entries; their individual proofs do not need to be repeated with matrix-valued diagonal weights.

## 3. Why this is not yet the adaptive theorem

An adaptive tree generally produces one continuation for every outcome history. Applying the preceding inequality independently to histories gives a sum of frontier normalizations. Nothing in the one-batch registry proves that this sum is at most one, or even at most $1.238062$. A complete proof must show that instrument completeness and branchwise dose combine all histories into one global cross-Gram factor or into an operator square function with the same bounded normalization.

The danger can already be seen in a two-bit classical toy model. Under the positive hypothesis let $(A,B)$ be uniform on $\{(-1,-1),(1,1)\}$, and under the negative hypothesis let it be uniform on $\{(-1,1),(1,-1)\}$. Observing $A$ has total variation zero, but after conditioning on either value of $A$, the posterior laws of $B$ are disjoint. Thus an unconditional child bound cannot be re-used after the root outcome. The joint two-coordinate experiment still has total variation one, so the example does not obstruct a globally grouped proof; it rejects only the nodewise-posterior shortcut.

The Round-3 warning that a universal near-unit two-copy temporal square function is false points in the same direction. Any proposed frontier lemma must be tested against that witness before promotion. This audit does not reconstruct that Round-2 witness because the present work is restricted to Round 4 and its imported Round-3 interface.

## 4. The precise missing lemma

The narrow sufficient statement is:

> **Global passive frontier lemma.** Fix any unrestricted classically adaptive passive tree with arbitrary fresh batches, idlers, entanglement, vacuum coherence, collective batch POVMs, classical feed-forward, and branchwise hard photon-pass dose at most six. After choosing the optimal signs of the terminal transcript differences and expanding the unconditioned attenuated hard-instance moments, all degree-four through degree-twelve physical occurrence terms admit one common operator-valued frontier such that: (i) each occurrence matrix is Schur-multiplied by a contractive cross-Gram symbol; (ii) the squared frontier masses define one law on the 210 occupation states; (iii) the exact pairing $m=n+a-2s$ is preserved; (iv) the total frontier normalization is at most $C_{\rm front}$ independently of outcome width, depth, and dose partition; and (v) $C_{\rm front}<1.238062042918586$.

If this lemma holds, the already certified factorization coefficients and the same 210-state Collatz matrix give

$$
\operatorname{TV}_{\rm adaptive}
\le C_{\rm front}P+E<T.
$$

The particularly clean case $C_{\rm front}=1$ would reuse the complete outward total without any adaptive overhead. It must not be asserted from deferred measurement or from the scalar registry alone.

The proof must also explain why promise conditioning is paid once after the unconditioned adaptive comparison. Reconditioning independently at child nodes would duplicate the promise loss and is not allowed.

## 5. Next bounded test

Before attempting an arbitrary-depth proof, construct the exact two-batch operator frontier with dose partitions $(1,5)$, $(2,4)$, and $(3,3)$. Keep the root instrument effects and all child continuations as operators. The test must determine whether their combined occurrence symbol is one cross-Gram contraction or, failing that, compute the sharp square-function normalization. It must include the inherited two-node witness and posterior choices with rare root outcomes.

Pass if the worst two-batch normalization is proved below the dynamic cap $1.238062042918586$ and the proof exposes an induction invariant. Fail if a physical two-batch passive tree exceeds that cap. A failure of a proposed frontier norm is a proof-architecture obstruction, not a passive protocol for the original promise problem unless the full transcript advantage is also evaluated.

## Reproduction

Run:

    /opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_adaptive_interface_gate.py --output artifacts/q64_adaptive_interface_gate.json
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_adaptive_interface_gate.py

The regression checks the exact caps against the committed outward ledger, rejects nodewise posterior reuse by direct enumeration, and stress-tests the contractive Schur-amplification identity on independently generated factorizations. Those checks protect this interface audit; they do not supply the missing global frontier lemma.
