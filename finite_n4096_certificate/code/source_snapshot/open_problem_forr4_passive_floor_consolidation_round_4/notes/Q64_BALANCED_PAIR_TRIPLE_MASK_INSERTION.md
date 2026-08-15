# Balanced pair--triple mask insertion at q=64

Date: 2026-07-16

Status: rigorous arbitrary-diagonal one-batch coefficient for all eight
remaining entries with a favorable split-cubic endpoint but no favorable
quintic endpoint. This raises the q64 theorem count from 340 to 348. It does
not prove the remaining 540 entries, intervalize the ledger, or establish
the adaptive lift.

## Result

The two canonical cuts are

$$
\begin{aligned}
(1,3,3,5)&:(0,1,3,2),\\
(3,1,3,5)&:(1,0,3,2).
\end{aligned}
$$

Complement and reversal generate all eight entries. Every quintic split is
balanced $2|3$. The common arbitrary-law coefficient is

$$
\boxed{0.642693497508}.
$$

After honest insertion, the routing ledger is

$$
P_{\rm total}=0.327343005018,
\qquad
{1\over3}-P_{\rm total}=0.00599032831576.
$$

Thus the theorem spends $0.00257667806424$ of routing margin but retains
more than the declared $10^{-3}$ reserve. There are now 80 quintic entries
left: 56 extreme and 24 balanced.

## Exact cubic endpoint

In every target cut, the split cubic is adjacent to the singleton on its
two-cell side. Assigning the cubic distinctness mask to that endpoint gives
the exact fixed-pair squared slice

$$
E_2=0.0153847346230,
\qquad
\gamma_2(A_{\rm cubic})\le\sqrt{E_2}=0.124035215254.
$$

This is the same arbitrary-law endpoint factor used by the dual-endpoint
theorem.

## Direct balanced quintic mask

Let $F$ be the two quintic cells on one side of the cut and $G$ the three
cells on the other. Their cross-cut distinctness symbol obeys the exact
identity

$$
\mathbf 1_{F\cap G=\varnothing}
=1-|F\cap G|+\mathbf 1_{F\subseteq G}.
$$

The three terms factor respectively through a constant coordinate,
singleton-incidence coordinates, and two-subset coordinates. Their row
squared multiplicities are

$$
(1,2,1),
$$

and their column squared multiplicities are

$$
(1,3,3).
$$

For positive direct-sum scales $x_i$, the product of the maximum feature
norms is

$$
\sqrt{
\left(\sum_i a_ix_i\right)
\left(\sum_i{b_i\over x_i}\right)}.
$$

Cauchy--Schwarz shows that the minimum is attained at
$x_i=\sqrt{b_i/a_i}$ and equals

$$
\boxed{
\gamma_{2,3}
\le\sum_i\sqrt{a_ib_i}
=1+\sqrt6+\sqrt3
=5.18154055035.
}
$$

This direct factor is substantially smaller than the sequential
singleton-versus-triple charge $\gamma_{N,3}^2=7.45548065275$ used as a
generic fallback in the preceding endpoint theorem.

## Complete arbitrary-law contraction

Temporarily complete every signed-permutation link moment across repeated
coordinates. Each completed link is a cross Gram of unit physical
character features and therefore has trace-class Schur-multiplier norm at
most one. The unsplit cubic introduces no cross-cut within-support mask.

Restore the cubic mask through its exact endpoint factor and the balanced
quintic mask through the direct factor above. Schur-multiplier
submultiplicativity gives

$$
\begin{aligned}
\gamma_2(K)
&\le \sqrt{E_2}(1+\sqrt6+\sqrt3)\\
&=0.642693497508.
\end{aligned}
$$

No product-law or invariant-law assumption is used. All physical links and
both support masks are composed once on the complete occurrence matrix.

## Regression scope

The regression protects the two-orbit topology, exact coefficient and
ledger insertion, deterministic artifact, and explicit feature
factorization for every pair/triple at dimensions five and seven. The
cubic endpoint slice and completed-link Gram principle are inherited proved
dependencies.

Reproduce with:

    /opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_balanced_pair_triple_mask_insertion.py --write-artifact
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_balanced_pair_triple_mask_insertion.py
