# Dimension-sensitive two-batch contraction

Date: 2026-07-14

Status: proved for two fixed nonadaptive dose-one probes and the minimal
four-chain sector.  This is a replacement benchmark, not an adaptive
dose-six theorem.

## 1. Averaged two-probe operator

Give each of two probes one photon.  A mode is a pair \((b,i)\), where
\(b\in[4]\) is the input block and \(i\in[N]\).  For fixed probe laws
\(q_1,q_2\), average the two-probe density matrix against the one-sided
minimal moment tensor

$$
T_{ijkl}=H_{ij}H_{jk}H_{kl}.
$$

Only matrix entries whose four ket/bra modes meet every block exactly once
survive.  Call the resulting Hermitian operator \(\Omega(q_1,q_2)\).
Every separate or joint measurement on the two probes has minimal-sector
transcript mass at most

$$
\|\Omega(q_1,q_2)\|_1.
\tag{1.1}
$$

This combines all temporal placements before the terminal absolute value.

## 2. Three exact flattenings

The operator splits as a direct sum over the three unordered \(2+2\)
partitions of the four blocks.  The two crossing chain flattenings have
rank \(N^2\) and nuclear norm \(N^{3/2}\); the adjacent flattening has rank
\(N\) and nuclear norm \(N\).

Let \(Q_{h,b}\) be probe \(h\)'s total one-photon mass in block \(b\).
Weighted Frobenius-to-nuclear estimates give

$$
\begin{aligned}
\|\Omega(q_1,q_2)\|_1
\le 2\sum_{\{A,A^c\}}
w_A\bigg[
&\sqrt{\prod_{b\in A}Q_{1,b}\prod_{b\in A^c}Q_{2,b}}\\
+&\sqrt{\prod_{b\in A^c}Q_{1,b}\prod_{b\in A}Q_{2,b}}
\bigg],
\end{aligned}
\tag{2.1}
$$

where the two crossing weights are \(w_A=N^{-1/2}\) and the adjacent
weight is \(w_A=N^{-1}\).

## 3. Exact block-mass maximum

Put \(a_b=\sqrt{Q_{1,b}}\) and \(b_b=\sqrt{Q_{2,b}}\).  Group the four
crossing terms relative to the adjacent split \(01\mid23\).  Their squared
pair-feature mass is

$$
(Q_{h,0}+Q_{h,1})(Q_{h,2}+Q_{h,3})\le{1\over4}.
$$

If \(u_h=Q_{h,0}+Q_{h,1}\) and
\(x_h^2=u_h(1-u_h)\), the adjacent pair-feature mass is at most

$$
{1\over4}-{x_h^2\over2}.
$$

Cauchy--Schwarz across the two probes now gives the sharp maximum of
(2.1):

$$
\boxed{
\|\Omega(q_1,q_2)\|_1
\le {1\over2\sqrt N}+{1\over4N}.
}
\tag{3.1}
$$

Uniform block masses attain the block-mass bound.  With uniform coordinate
weights, direct diagonalization attains (3.1) itself.

## 4. Comparison with the rejected square function

For the uniform probes, the joint occupation square mass is
\(\mathcal S=3/32\).  Hence the exact constant relative to
\(N^{-1/2}\sqrt{\mathcal S}\) is

$$
{1/2+1/(4\sqrt N)\over\sqrt{3/32}}
\longrightarrow\sqrt{8/3}\approx1.632993.
\tag{4.1}
$$

This explains both sides of the earlier falsification:

- a universal constant one is false;
- retaining the exact chain flattenings gives a finite constant well below
  the crude \(\sqrt6\) small-\(N\) value.

At \(N=1024\), (3.1) equals

$$
{1\over64}+{1\over4096}
\approx0.0158691.
$$

## 5. Nonadaptive consolidation and the adaptive boundary

More generally, fixed nonadaptive probes with parity-support laws
\(q_1,\ldots,q_m\) consolidate exactly.  Map a tuple of supports to its XOR
\(R=S_1\triangle\cdots\triangle S_m\), and let \(Q_R\) be the pushforward
of the product law.  The weighted moment matrix on support tuples is an
isometric lift of

$$
D_Q^{1/2}C D_Q^{1/2},
\qquad C(R,R')=c(R\triangle R').
$$

Since \(|R|\le\sum_h|S_h|\le D\), every fixed nonadaptive schedule is
bounded by the one-batch dose-\(D\) Schur-symbol optimum.

Outcome-selected child probe laws break this XOR pushforward: there is no
single product law before the root outcome is measured.  Thus (3.1) and
nonadaptive consolidation do not settle unrestricted connected
adaptivity.  The next theorem must retain the dimension-sensitive
flattening norm through outcome-selected preparations rather than revert
to the false unweighted square function.

Reproduction: tests/two_batch_trace_norm.py.
