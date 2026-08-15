# Audit of the interpolation-to-reverse-tree handoff

Date: 2026-07-14

Status: analytic integration audit for the repaired fixed-skeleton
contraction.

## 1. Form of one terminal interpolation term

Condition on the three interpolation times, the latent Gaussian variables,
the fresh-noise boundary variables, and one terminal Stein-transfer path.
For a terminal graph \(G=(V,E)\), the label sum has the form

$$
\sum_{\substack{\iota:V\to[N]\\
\iota\ {\rm injective\ within\ each\ layer}}}
\left(\prod_{\{u,v\}\in E}H_{\iota(u),\iota(v)}\right)
\left(\prod_{v\in V}w_v(\iota(v))\right)
\partial_{\iota(V)}F_{\mathcal T}(\mu).
\tag{1.1}
$$

Every \(w_v\) is a derivative of the fixed functions \(\psi'\) or
\(\mathcal S_\psi\) evaluated at the Gaussian variable attached to \(v\).
Multiple derivatives at one marked coordinate remain a single local
factor. The branching potential permits at most twelve vertices and at
most six differentiations of existing weights, so

$$
\|w_v\|_\infty\le C_{\tau,\beta}.
$$

The expectation over latent Gaussians is outside (1.1). Therefore the
weights need not be independent random variables: after conditioning, they
are local diagonal multipliers of norm at most \(C_{\tau,\beta}\). They
cannot act as an arbitrary joint sign tensor that cancels the Hadamard
phases.

The distinct-label indicators in (1.1) are expanded by the Walsh repair
before components are tensorized. Local vertex weights and Walsh
characters are then absorbed into diagonal contractions on the graph
vertex spaces. They change only the diagram constant.

## 2. Uniformity of the fixed-skeleton estimate

The repaired reverse-tree theorem is uniform over:

1. every bias/base value with modulus at most one;
2. every terminal transcript sign of modulus at most one;
3. every adaptive outcome alphabet and outcome-selected descendant;
4. every fixed assignment of marked vertices to physical nodes and ket/bra
   entries; and
5. every local diagonal vertex multiplier of bounded norm.

Consequently (1.1) is bounded before the interpolation times and latent
variables are integrated. Tonelli/triangle inequality then introduces only
the bounded time volume and the finite path multiplicity \(C^v\).

Antisymmetrizing the transcript statistic before absolute values ensures
that every retained terminal graph is hypothesis-sensitive. The
interpolation audit proves the stronger task-specific fact that every
connected component spans all four layers.

## 3. Adaptive marked-time sum

Let \(S_h\) be the branch potential

$$
S_h=\sqrt{t_h}+\max_y S_{hy},
\qquad S_{\rm leaf}=0.
$$

Because every positive dose is an integer and
\(t_h+\max_yD_{hy}\le D_h\),

$$
S_h\le t_h+\max_yD_{hy}\le D_h.
$$

For \(k\) labeled marks remaining at node \(h\), choose a subset of \(i\)
marks for the current node. Ket/bra splitting and joint Bessel mass cost at
most \((2\sqrt{t_h})^i\), up to the fixed factorial diagram constant. The
other \(k-i\) marks enter one outcome-selected child. Complete-frame
packing or unmarked stochasticity removes the current outcome without an
alphabet-size factor. Inductively,

$$
\begin{aligned}
F_h(k)
&\le\sum_{i=0}^k {k\choose i}
(2\sqrt{t_h})^i
\max_y(2S_{hy})^{k-i}\\
&\le(2S_h)^k
\le(2D_h)^k.
\end{aligned}
\tag{3.1}
$$

This is the tree-level version of the branchwise multinomial identity. It
shows explicitly why adaptive width does not turn the marked-time sum into
a sum over all nodes of the tree.

## 4. Result of the handoff

For every terminal graph with \(v\le12\), the repaired fixed-skeleton
estimate and (3.1) give

$$
|\mathfrak C_G(\mathcal T)|
\le C_G(1+D)^vN^{-1/2}.
$$

The finite Stein-path sum, bounded local weights, time integration, Walsh
mask expansion, ket/bra splits, collision factors, and conditioning
constants are all absorbed into an \(N\)- and \(D\)-independent constant.
Taking \(v\le12\) yields the claimed transcript bound

$$
\operatorname{TV}\le C(1+D)^{12}N^{-1/2}+O(N^{-1}).
$$

The integration audit found no additional power of \(N\), dose, tree
depth, or outcome width. The numerical regression
tests/adaptive_mark_assignment_ledger.py checks the exact recurrence (3.1)
on adversarial finite trees; the proof is the binomial induction above.
