# Joint-probe square mass for the minimal chain

Date: 2026-07-14

Status: proved occupation identity and finite-size constant. The proposed
near-unit temporal/frame square-function inference is false; see
notes/two_copy_square_function_counterexample.md.

## 1. Local square-mass polynomial

The minimal chain has one marked coordinate in each of four blocks. At a
node \(h\), let \(q_h\) be the probe's parity-support law and

$$
n_b(S)=|S\cap\text{block }b|.
$$

If the block set \(A\) is placed in the ket entry and the disjoint block
set \(B\) in the bra entry, complete-frame packing gives amplitude mass

$$
\sqrt{m_{h,A}m_{h,B}},
\qquad
m_{h,A}=\mathbb E_{S\sim q_h}\prod_{b\in A}n_b(S).
$$

Its square is \(m_{h,A}m_{h,B}\). Introduce square-free variables
\(z_1,\ldots,z_4\). The local generating polynomial for all ket/bra
placements is

$$
g_h(z)=
\sum_{\substack{A,B\subseteq[4]\\A\cap B=\varnothing}}
m_{h,A}m_{h,B}z_{A\cup B}.
\tag{1.1}
$$

Draw independent \(S_h,T_h\sim q_h\). Expanding block by block gives the
exact identity

$$
g_h(z)
=\mathbb E_{S_h,T_h}
\prod_{b=1}^4
\left[1+\bigl(n_b(S_h)+n_b(T_h)\bigr)z_b\right].
\tag{1.2}
$$

Thus one probe distribution is shared by every placement. Optimizing every
placement's moment separately, as in the rejected ledgers, discards this
constraint.

## 2. Sum over all temporal placements

For a fixed branch with nodes \(h\), the sum of squared Bessel masses over
every assignment of the four block marks to nodes and ket/bra entries is
the coefficient

$$
\mathcal S=[z_1z_2z_3z_4]\prod_hg_h(z).
$$

Using (1.2), set

$$
M_b=\sum_h\bigl(n_b(S_h)+n_b(T_h)\bigr).
$$

Coefficient extraction gives

$$
\mathcal S=\mathbb E\prod_{b=1}^4M_b.
\tag{2.1}
$$

Every support at node \(h\) has size at most its hard dose \(t_h\).
Consequently, pointwise,

$$
\sum_{b=1}^4M_b
\le2\sum_ht_h
\le2D.
$$

Arithmetic--geometric mean in (2.1) proves

$$
\boxed{
\mathcal S\le\left({D\over2}\right)^4.
}
\tag{2.2}
$$

At \(D=6\),

$$
\sqrt{\mathcal S}\le9,
\qquad
{9\over\sqrt{1024}}={9\over32}=0.28125<{1\over3}.
\tag{2.3}
$$

The estimate is independent of the number and sizes of the batches. Exact
minimal-chain graph weights are at most the factored \(N^{-1/2}\) scale, so
retaining them can only lower \(\mathcal S\).

## 3. What this does and does not prove

Equation (2.2) is a rigorous joint-probe occupation bound. It removes the
factor \(7.97\) obstruction created by optimizing each placement
independently.

It is not yet a passive transcript bound. The remaining theorem must show
that the reverse adaptive contraction combines marked-time placements with
the square mass \(\mathcal S\), up to a constant below

$$
{32/3\over9}={32\over27}\approx1.18519.
$$

A generic triangle inequality fails.  The exact two-node Hadamard-frame
witness in notes/two_copy_square_function_counterexample.md requires
constant \(\sqrt6\) at \(N=1\) and still exceeds \(32/27\) after its
one-coordinate embedding at \(N=2\).  Therefore (2.2) remains useful
occupation algebra, but it cannot be promoted by the proposed universal
near-unit frame lemma.

Reproduction: tests/joint_probe_square_mass.py.
