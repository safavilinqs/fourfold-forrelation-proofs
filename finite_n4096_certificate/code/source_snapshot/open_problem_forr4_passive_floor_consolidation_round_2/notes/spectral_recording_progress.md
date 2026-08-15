# Spectral recording progress for the exact plant

Date: 2026-07-14

## 1. Pure match sectors are bosonic Hadamard compounds

Fix one signed-permutation pair $(X,Y)$. Consider a moment sector with $r$ distinct odd column labels on the $X$ side and $r$ distinct odd row labels on the $Y$ side, with one physical occurrence per odd label.

Index the left modes by $(x_i,y_i)$ with distinct $y_i$, and the right modes by $(u_j,v_j)$ with distinct $u_j$. The exact permutation formula gives

$$
M_r(\mathbf x,\mathbf y;\mathbf u,\mathbf v)
={1\over(q)_r}
\sum_{\sigma\in S_r}
(-1)^{\sum_i x_i\cdot u_{\sigma(i)}
             +\sum_i v_{\sigma(i)}\cdot y_i}.
$$

Let $H_N=H_q\otimes H_q$, $N=q^2$. The permanent compound of the corresponding $r\times r$ submatrix satisfies

$$
\operatorname{perm}(H_N[\mathrm L,\mathrm R])
=q^{-r}
\sum_{\sigma\in S_r}(-1)^{\cdots}.
$$

Therefore

$$
\boxed{
M_r={q^r\over(q)_r}
P_{\rm right}\,\operatorname{Sym}^r(H_N)\,P_{\rm left},}
$$

where the $P$'s restrict to collision-free odd-label sectors and normalized distinct-mode bosonic basis vectors are understood.

Since $\operatorname{Sym}^r(H_N)$ is unitary,

$$
\boxed{\|M_r\|_{\rm op}\le {q^r\over(q)_r}.}
$$

The small-$q$ calculation shows equality for every tested sector. At $q=32$ the factor remains modest for all $r$ accessible at dose six:

| $r$ | $q^r/(q)_r$ |
|---:|---:|
| 1 | 1 |
| 2 | 1.03226 |
| 3 | 1.10108 |
| 4 | 1.21498 |
| 5 | 1.38855 |
| 6 | 1.64569 |

## 2. Interpretation

The hidden-label probability $1/\binom qr$ is not the correct operator scale. Coherent character summation promotes the pure record sector to an almost-unitary bosonic transform. This exactly explains RT-004: for $r=1$, the link is the full Hadamard unitary, not an entrywise $1/q$ contraction.

The useful fact is different: sector amplification is controlled by the explicit falling-factorial ratio above, which is close to one at $q=32$, $r\le6$. A three-link sector therefore carries a small explicit representation-theoretic constant rather than an unknown dimension loss.

## 3. Candidate spectral contraction

Decompose every signed-permutation moment into orthogonal odd-label sectors $(r_1,r_2,r_3)$ and even-pair decorations. For each pure sector:

1. replace each hidden-permutation link by the compressed bosonic Hadamard transform $M_{r_a}$;
2. retain its exact rank and flat diagonal under every physical ket/bra grouping;
3. apply the repaired global reverse-tree dichotomy to this compound network; and
4. square-sum inequivalent permutation-representation sectors before summing adaptive histories.

The crossing minimal-chain sector rules out a uniform $1/N$ prefactor. The sharp possible form is

$$
\operatorname{TV}
\le {1\over\sqrt N}
\sum_{r_1,r_2,r_3}
c_{r_1,r_2,r_3}(D)
\prod_{a=1}^3{q^{r_a}\over(q)_{r_a}},
$$

with the exact integer-dose-six sum below $32/3$.

The adjacent-pair minimal sector attains $1/N$, but either crossing-pair sector attains $1/\sqrt N$. The latter is the correct uniform obstruction and must be retained in every finite-size estimate.

## 4. Remaining obstacles

- Moments may contain even repeated labels in addition to the odd record set. These must be shown to be contractions or incorporated as lower-sector multiplicities.
- Block-two and block-three product masks couple the outgoing sector of one permutation to the incoming sector of the next.
- Adaptive complete outcomes can mix representation sectors. The needed square-sum must follow from frame completeness, not from assuming orthogonality after outcome selection.
- The final constants must be evaluated at the exact occurrence profiles with total size at most twelve; an unspecified polynomial in $D$ is insufficient.

The next falsification target is a two-link or three-link sector in which outcome-selected frame weights destroy the expected orthogonal square-sum.
