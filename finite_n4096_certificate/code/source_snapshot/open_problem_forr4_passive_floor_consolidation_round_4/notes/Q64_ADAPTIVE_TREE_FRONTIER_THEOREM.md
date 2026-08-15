# $q=64$ classically adaptive tree-frontier theorem

VERDICT: CERTIFIED

Date: 2026-07-18

The complete $q=64$ balanced one-batch bound extends with no multiplicative loss to every finite classically adaptive passive tree whose fresh probes are block diagonal in total signal photon number. The occurrence-frontier multiplier is exactly one. Therefore the outward transcript bound is

$$
0.26096922479220792493418095739381656137869710501965562750914672535608982983539778
<
\frac13-10^{-3},
$$

with certified reserve

$$
0.071364108541125408399152375939516771954636228313677705824186607977243503497935553.
$$

This closes the adaptive-lift obligation for classical feed-forward between number-sector-incoherent fresh quantum batches. It does not cover coherence between different total signal-number sectors, coherent quantum memory, or coherent control carried between batches.

## 1. Model and signed transcript functional

Fix a finite passive tree. At each node the protocol prepares a fresh signal--idler state that commutes with total signal photon number; sends the signal through the diagonal sign oracle; performs an arbitrary collective POVM on that batch and its idlers; records the classical outcome; and chooses the next fresh batch from the classical history. The state may have arbitrary idlers, repeated physical modes, and within-sector entanglement and coherence. Every root-to-leaf branch has total photon-pass dose at most six.

Let $P_+$ and $P_-$ be the two unconditioned hard-instance transcript laws. Their total variation is the supremum of the signed transcript functional over terminal weights $h_\tau\in[-1,1]$, with the same conventional factor $1/2$ used in the one-batch ledger. It is therefore enough to bound the expectation difference for an arbitrary fixed collection of terminal signs.

A number-sector-incoherent mixed probe has a spectral decomposition into pure states supported in fixed total signal-number sectors. By convexity, reveal the decomposition label to the controller and include it in an orthogonal direct sum; this can only enlarge the strategy class. Randomized controllers are treated the same way. It is enough to prove the factorization for pure fixed-number probes and finite POVMs. The direct-integral version gives the same argument for general outcome spaces, or equivalently finite coarse grainings converge to the transcript total variation.

## 2. Normalized strategy-factorization lemma

For a subtree, index a complete temporal ket history by $I$ and a complete temporal bra history by $J$, including every Fock occupation label selected along that branch. After removing the oracle characters, let $D(I,J)$ be the coefficient of the signed transcript functional. The induction invariant is

$$
D(I,J)=\langle v_J,u_I\rangle,
\qquad
\sum_I\lVert u_I\rVert^2\le1,
\qquad
\sum_J\lVert v_J\rVert^2\le1.
\tag{2.1}
$$

This invariant includes the outcome-selected probe amplitudes; it is not the rejected unweighted temporal square function.

Conditional on the refined classical label, write the first fixed-number fresh probe as

$$
|\psi\rangle=\sum_i |i\rangle|\alpha_i\rangle,
\qquad
\sum_i\lVert\alpha_i\rVert^2=1,
\tag{2.2}
$$

where $i$ is the first-batch Fock occupation, all occupied $i$ have the same total signal photon number, and $|\alpha_i\rangle$ is an arbitrary idler vector. Let $\{E_y\}_y$ be the root POVM and define

$$
a_{y,i}=E_y^{1/2}(|i\rangle|\alpha_i\rangle).
\tag{2.3}
$$

POVM completeness gives

$$
\sum_y\lVert a_{y,i}\rVert^2=\lVert\alpha_i\rVert^2,
\qquad
\sum_{y,i}\lVert a_{y,i}\rVert^2=1.
\tag{2.4}
$$

For each outcome $y$, apply the induction hypothesis to the arbitrarily chosen child subtree, including its $y$-dependent probe and all later terminal signs:

$$
D_y(I,J)=\langle v_{y,J},u_{y,I}\rangle,
\qquad
\sum_I\lVert u_{y,I}\rVert^2\le1,
\qquad
\sum_J\lVert v_{y,J}\rVert^2\le1.
\tag{2.5}
$$

Because the child is fresh and is selected only by the classical outcome, the complete signed coefficient is

$$
D((i,I),(j,J))
=
\sum_y
\langle a_{y,j},a_{y,i}\rangle
\langle v_{y,J},u_{y,I}\rangle.
\tag{2.6}
$$

Define the direct-sum features

$$
U_{i,I}=\bigoplus_y a_{y,i}\otimes u_{y,I},
\qquad
V_{j,J}=\bigoplus_y a_{y,j}\otimes v_{y,J}.
\tag{2.7}
$$

They reproduce (2.6), and their total squared masses satisfy

$$
\begin{aligned}
\sum_{i,I}\lVert U_{i,I}\rVert^2
&=
\sum_y
\left(\sum_i\lVert a_{y,i}\rVert^2\right)
\left(\sum_I\lVert u_{y,I}\rVert^2\right)\\
&\le
\sum_{y,i}\lVert a_{y,i}\rVert^2
=1.
\end{aligned}
\tag{2.8}
$$

The column calculation is identical. At a leaf, $|h_\tau|\le1$ supplies the scalar base factorization. This proves (2.1) for arbitrary width, depth, outcome probabilities, outcome-selected child probes, and branch-dependent dose partitions. Rare outcomes do not introduce inverse probabilities because the proof never normalizes a posterior law.

## 3. Pullback to the one-batch occurrence kernels

For a temporal history $I=(i_1,\ldots,i_d)$, let

$$
f(I)=i_1+\cdots+i_d
\tag{3.1}
$$

be its aggregate physical Fock occupation. Diagonal sign characters multiply, so

$$
\chi_{i_1}(x)\cdots\chi_{i_d}(x)=\chi_{f(I)}(x).
\tag{3.2}
$$

Repeated physical modes are retained in the integer occupation $f(I)$; their parity is taken automatically by the sign character. Orthogonality of the refined sector labels makes the cross-Gram coefficient zero unless ket and bra have equal total signal photon number in every batch. Thus every surviving aggregate pair has $|f(I)|=|f(J)|$. The branchwise hard-dose condition implies $|f(I)|\le6$. Histories of different depths may be padded by vacuum batches.

Consequently every averaged hard-instance moment in the adaptive expansion is the pullback

$$
K_{f(I),f(J)}
\tag{3.3}
$$

of the same balanced degree-four through degree-twelve physical occurrence kernels certified in the complete one-batch registry. Distinct temporal decompositions of one aggregate occupation merely duplicate a row or column. Such pullback does not increase a Hilbert factorization norm: a factorization $K_{s,t}=\langle b_t,c_s\rangle$ pulls back by replacing $s,t$ with $f(I),f(J)$.

Set $p_I=\lVert U_I\rVert^2$ and $w_J=\lVert V_J\rVert^2$. After normalizing nonzero features, (2.1) writes the adaptive occurrence block as

$$
\operatorname{diag}(\sqrt p)
\bigl(K\circ G\bigr)
\operatorname{diag}(\sqrt w),
\tag{3.4}
$$

where $G$ is one contractive cross-Gram symbol and both $p$ and $w$ have total mass at most one. Tensoring the certified factorization of $K$ with the normalized temporal features proves that every accepted arbitrary-diagonal coefficient is unchanged. This is precisely the individual Schur-amplification lemma from the interface audit, now supplied with one global normalized feature family for the whole tree.

This also resolves the earlier Round-2 nonadaptive/adaptive boundary correctly. Outcome-selected child probes do prevent a product-law XOR pushforward, but the present $p$ and $w$ are not asserted to be product laws. They are single normalized, generally correlated laws on complete temporal histories. The completed $q=64$ coefficient registry was proved for arbitrary correlated diagonal row and column laws, which is exactly the extra scope needed here.

## 4. The 210-state Perron interface accepts two laws

Push $p$ and $w$ forward by the four block-occupation vector of $f(I)$ and $f(J)$. For one profile and one split, exact occupation compatibility remains

$$
m=n+a-2s.
\tag{4.1}
$$

The resulting complete ledger is a bilinear form $x^TBy$, where $x_n=\sqrt{p_n}$, $y_m=\sqrt{w_m}$, and $B$ is the same nonnegative symmetric 210-state matrix used by the outward certificate, block diagonal in total signal photon number. The two feature-mass bounds imply $\lVert x\rVert_2,\lVert y\rVert_2\le1$. Therefore

$$
x^TBy
\le
\lVert B\rVert_{2\to2}
=
\lambda_{\max}(B).
\tag{4.2}
$$

For a symmetric entrywise nonnegative matrix, the spectral norm equals its Perron root. Thus the existing outward Collatz--Wielandt value bounds independent row and column laws as well as the common law displayed in the original one-batch derivation. No temporal multiplier is spent.

## 5. Promise conditioning and numerical consequence

The factorization applies to the unconditioned attenuated hard pair. Promise conditioning is then paid once at the two distribution endpoints by the same total-variation conditioning inequality used in the complete ledger. It is not applied at internal posteriors. Hence

$$
\operatorname{TV}_{\rm adaptive}
\le
0.25874409638557722267923079152925562403311987598701005323368022563713944348400092
+
0.00222512840663070225495016586456093734557722903264557427546649971895038635139684472693787993,
$$

which is bounded outward by the certified total at the start of this note. The adaptive frontier multiplier is $C_{\rm front}=1$, strictly inside the former cap $1.238062042918586$.

## 6. Independent checks and the inherited witnesses

The independent regression reconstructs the old two-copy Hadamard witness exactly. Its signed output matrix has diagonal $3/32$, off-diagonal $-1/32$, and transcript mass $3/4$. It then constructs the full 64-coordinate direct-sum strategy features and verifies exact feature mass $1/16$ on each of the sixteen temporal histories, total row and column masses one, and exact recovery of the transcript functional. The old ratio $\sqrt6$ compares transcript mass with a discarded joint-square proxy of mass $3/32$; it is not an adaptive multiplier over the physical occurrence kernel.

The regression also checks eighty independently generated complex POVM/child-factorization instances, a root outcome of mass below $7\times10^{-15}$ without posterior division, direct repeated-mode character pullbacks for $(1,5)$, $(2,4)$, and $(3,3)$, one hundred duplicated-history pullbacks, and one hundred independent-law Perron inequalities. These computations validate the algebra and catch conjugation, outcome-normalization, rare-outcome, dose-partition, duplicate-history, and common-law assumptions; the proof itself is the exact argument in Sections 2 through 4.

Run:

    /opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_adaptive_tree_frontier.py --output artifacts/q64_adaptive_tree_frontier.json
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_adaptive_tree_frontier.py

## 7. Scope and remaining paper work

The theorem covers unrestricted classical adaptivity between number-sector-incoherent fresh passive batches. It permits arbitrary child probes satisfying that restriction and selected by the complete classical history, and it does not assume product laws, uniform outcomes, exchangeability, fixed dose partitions, or lower bounds on outcome probabilities.

The theorem does not cover coherence between distinct total signal-number sectors, coherent quantum memory passed between batches, or coherent control of later probes. The broader number-sector-coherent extension has 272 unbalanced high-sector split/state incidences on 136 occupation edges whose physical coefficients remain unproved. The theorem also does not by itself establish experimental credibility for $16{,}384$ sign modes.
