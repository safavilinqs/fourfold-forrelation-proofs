# Barrier for a placement-only temporal square function

Date: 2026-07-14

Status: exact rejection of a generic \(\ell_2\) repair of the adaptive
minimal-chain ledger.

## 1. Worst temporal partition

Take six sequential dose-one nodes and the four labeled minimal-chain
marks. At a dose-one node, a nonzero term can place at most one mark in its
ket entry and at most one in its bra entry. A node receiving one mark has
local mass one; a node receiving one ket and one bra mark has optimal local
mass \(1/2\).

Every nonzero placement is all-singleton, so its graph factor relative to
the target \(N^{-1/2}\) is one. Exact enumeration of the \(12^4\) entry
placements gives

$$
\sum_p b_p=8730,
\qquad
\sum_p b_p^2={14445\over2}.
$$

Thus even an ideal square-sum over all temporal placements yields

$$
\left(\sum_p b_p^2\right)^{1/2}
=\sqrt{14445/2}
\approx84.9853.
$$

At \(N=1024\), the threshold coefficient is \(32/3\approx10.6667\).
The placement-only square function misses it by a factor

$$
{3\over32}\sqrt{14445/2}\approx7.96737.
$$

## 2. Consequence

A theorem that merely replaces the rejected temporal \(\ell_1\) sum by an
unweighted \(\ell_2\) sum cannot certify dose greater than six. The needed
finite-size argument must also exploit at least one of:

1. orthogonality or cancellation between signed-permutation representation
   sectors;
2. incompatibility of independently optimized local probe masses across
   nodes and placements;
3. a stronger tester norm that keeps the entire adaptive comb coupled; or
4. a different hard instance with a substantially smaller dose-one
   temporal constant.

This barrier does not describe an achievable protocol. Each \(b_p\)
optimizes local masses independently, so the square sum remains an upper
ledger. It quantifies how much more structure a successful contraction
must retain.

Reproduction: tests/temporal_square_function_barrier.py.
