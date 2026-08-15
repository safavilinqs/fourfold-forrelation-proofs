VERDICT: CERTIFIED

# Joint shared-quintic contraction for the twelve recovered entries

Date: 2026-07-17

The twelve quarantined $(1,3,5,3)$ entries are proved at $q=64$. The proof never multiplies independently optimized cubic--quintic maxima. It first fixes the shared physical quintic, classifies its joint row/column multiplicity shape, and bounds the left endpoint--cubic--quintic factor and the right quintic--cubic factor within that same shape. Exact physical completion incidences then give arbitrary-correlated-diagonal row, column, and rank bounds.

The three outward-rounded canonical coefficients are

$$
c_{(0,1,3,2)}=0.33828697324447987,
\qquad
c_{(0,2,2,2)}=0.11863696369021283,
\qquad
c_{(0,2,3,1)}=0.3144332243430522.
$$

All are strictly below one. Complement and reversal therefore prove all twelve entries.

## 1. Physical sector matrix

Let $q=2^m\ge8$, let $\Omega=\mathbb F_2^m\times\mathbb F_2^m$, and let $M(S,T)$ be the exact signed-permutation link moment defined in the independent audit. For profile $d=(1,3,5,3)$ and split $t$, write the four simple full supports as

$$
U_i=A_i\sqcup B_i,
\qquad
|A_i|=t_i,
\qquad
|B_i|=d_i-t_i.
$$

The disjoint union is the physical cross-cut mask: every $U_i$ is a simple support, so $A_i\cap B_i=\varnothing$. In record sector $\rho=(1,r,s)$, the occurrence matrix is

$$
\mathcal O_t^\rho(A,B)
=
\mathbf 1\{A_i\cap B_i=\varnothing\text{ for all }i\}
M(U_0,U_1)M(U_1,U_2)M(U_2,U_3),
$$

where the cubic $U_1$ has row record one and column record $r$, the shared quintic $U_2$ has row record $r$ and column record $s$, and the terminal cubic $U_3$ has row record $s$.

For arbitrary probability laws $p$ and $w$ on the complete row and column tuples, including arbitrary correlations among blocks, the target is

$$
\left\|
\operatorname{diag}(\sqrt p)\,
\mathcal O_t^\rho\,
\operatorname{diag}(\sqrt w)
\right\|_1.
$$

No product-law, exchangeability, translation-invariance, or uniform-law assumption is used.

## 2. Link constants used inside one fixed quintic shape

Define the exact rational quantities

$$
e_1=\frac1{q(q-1)},
\qquad
e_3=\frac1q,
\qquad
b=\frac{q+2}{q(q-1)(q-2)},
$$

and

$$
g=\binom q3^{-1},
\qquad
d=\frac3{(q-3)\binom q3}.
$$

For $q\ge8$, $d<g$.

The endpoint factors $e_1,e_3$ are the exact maxima from the previous audit. A cubic with endpoint record one must be an L shape, and its singleton link has magnitude $e_1$. A cubic with endpoint record three is bounded by $e_3$, attained by a horizontal triple.

### Record-three quintic side

If the relevant quintic axis has multiplicity pattern $311$, it has three odd groups and no even group. The record-three odd-block permanent has at most six terms of modulus one, so

$$
|M_{35}|\le g.
$$

If the pattern is $2111$, the even pair has a nonzero xor. For each odd matching, its remaining character average has magnitude at most $3/(q-3)$, hence

$$
|M_{35}|\le d.
$$

The same statements apply to $M_{53}$ after transposing the cubic and quintic.

### Record-one quintic side

An endpoint-compatible record-one cubic is an L shape and has one nonzero even-column xor. Condition on the unique odd row of the quintic mapping to the unique odd column of the cubic. The remaining permutation is uniform on $q-1$ points. For quintic pattern $32$, its even pair has nonzero xor and the resulting two-point character average is bounded by

$$
|M_{35}|\le b.
$$

For pattern $41$, the even four-group may have zero xor. If it does, only the cubic's nonzero even xor remains and the conditional character average has magnitude $1/(q-1)$, giving $1/[q(q-1)]\le b$ after the odd match. If the four-group xor is nonzero, the same two-point character calculation gives at most $b$. Thus $b$ covers both xor subclasses of $41$; it does not assume that the four-group is active.

The remaining $221$ quintic pattern has two nonzero even-pair xors, which may be equal or distinct. Put $n=q-1$. If the cubic has no even group, the conditional record-one average is at most

$$
\frac{n+1}{n(n-1)},
$$

before the final factor $1/q$, which is below $b$. If the cubic has one even group, the two single-even terms contribute at most $4/[n(n-1)]$. The remaining pair term contributes at most

$$
\frac{3(n+3)}{n(n-1)(n-2)}.
$$

Thus the conditional average is at most

$$
\frac{7n+1}{n(n-1)(n-2)}.
$$

For $q\ge8$, equivalently $n\ge7$,

$$
7n+1\le(n+3)(n-2)
$$

because the difference is $(n-7)(n+1)$. After the final factor $1/q$, this again gives $b$.

For the terminal record-one link, transpose the link so that the terminal cubic is on the left. A quintic column pattern $32$ or $221$ is therefore bounded by $b$. A pattern $41$ can have a zero xor in its even four-group, so only the universal value $1/q$ is safe there. The simple-support shape table below is what prevents that dangerous $41$ value from being combined with an independently optimized left link.

## 3. Complete joint quintic-shape table

The notation $221$, $32$, $41$, $2111$, and $311$ records the multiplicities of equal row or column labels among the five simple cells. Enumerating the two restricted-growth set partitions of the five labeled cells and rejecting repeated row--column pairs gives exactly the following 15 feasible joint shapes in the four sectors.

The xor information retained by the estimates is:

- every even pair in patterns $32$, $221$, and $2111$ has nonzero xor because the physical cells are distinct;
- the two even-pair xors in $221$ may coincide or differ, and the $221$ estimate covers both cases;
- the even four-group in $41$ may have zero or nonzero xor, and the proof uses $b$ on the endpoint-compatible left link but the safe value $1/q$ on a terminal record-one link;
- pattern $311$ has no even group.

| quintic row pattern | quintic column pattern | $(r,s)$ | left endpoint--middle bound | right bound |
|---|---|---:|---:|---:|
| $221$ | $221$ | $(1,1)$ | $e_1b$ | $b$ |
| $221$ | $32$ | $(1,1)$ | $e_1b$ | $b$ |
| $32$ | $221$ | $(1,1)$ | $e_1b$ | $b$ |
| $221$ | $2111$ | $(1,3)$ | $e_1b$ | $d$ |
| $221$ | $311$ | $(1,3)$ | $e_1b$ | $g$ |
| $32$ | $2111$ | $(1,3)$ | $e_1b$ | $d$ |
| $41$ | $2111$ | $(1,3)$ | $e_1b$ | $d$ |
| $2111$ | $221$ | $(3,1)$ | $e_3d$ | $b$ |
| $2111$ | $32$ | $(3,1)$ | $e_3d$ | $b$ |
| $2111$ | $41$ | $(3,1)$ | $e_3d$ | $1/q$ |
| $311$ | $221$ | $(3,1)$ | $e_3g$ | $b$ |
| $2111$ | $2111$ | $(3,3)$ | $e_3d$ | $d$ |
| $2111$ | $311$ | $(3,3)$ | $e_3d$ | $g$ |
| $311$ | $2111$ | $(3,3)$ | $e_3g$ | $d$ |
| $311$ | $311$ | $(3,3)$ | $e_3g$ | $g$ |

This table incorporates the within-quintic simple-support constraint. In particular:

- sector $(1,1)$ permits only column patterns $32$ and $221$, so its terminal link is bounded by $b$, not the unsafe universal $1/q$;
- a column pattern $41$ in sector $(3,1)$ forces row pattern $2111$, so the left middle link receives the smaller $d$ bound;
- the two record-three links can both take the larger $g$ bound only on the joint $311/311$ plus shape.

Taking the maximum within each sector gives the valid complete-chain entry bounds

$$
J_{11}=e_1b^2,
\qquad
J_{13}=e_1bg,
$$

$$
J_{31}=e_3d\frac1q,
\qquad
J_{33}=e_3g^2.
$$

For $J_{31}$, the other candidate is $e_3gb$. The displayed choice is larger for $q\ge8$ because

$$
\frac{d/q}{gb}
=
\frac{3(q-1)(q-2)}{(q-3)(q+2)}
>1.
$$

At $q=64$, the exact reduced fractions are

| sector | $J_{r,s}$ |
|---|---:|
| $(1,1)$ | $121/6999104028672$ |
| $(1,3)$ | $11/6999104028672$ |
| $(3,1)$ | $1/3470000128$ |
| $(3,3)$ | $1/111096889344$ |

These are joint bounds through one shared $U_2$. They are not products of separately optimized $M_{35}$ and $M_{53}$ maxima.

## 4. Complete-chain row, column, and rank energies

For a split $t$ and sector $(1,r,s)$, let

$$
I_0(k)=
\begin{cases}
q^2,&k=0,\\
1,&k=1,
\end{cases}
$$

let $I_1^{1,r}(k)$ be the two-axis simple cubic completion incidence, let $I_2^{r,s}(k)$ be the two-axis simple quintic completion incidence, and let $I_3^s(k)$ be the one-axis terminal cubic incidence. The two-axis quantities are the minimum of the two one-axis relaxations. Each one-axis count still chooses distinct cells within every occupied coordinate group, so taking their minimum is an upper bound on the physical intersection and never deletes within-block distinctness.

Define

$$
R_{r,s}(t)
=
\prod_{i=0}^3 I_i(t_i),
\qquad
C_{r,s}(t)
=
\prod_{i=0}^3 I_i(d_i-t_i).
$$

Every nonzero entry has magnitude at most $J_{r,s}$, and every fixed row has at most $R=R_{r,s}(t)$ physical completions. To see the arbitrary-law row bound explicitly, set $X=\operatorname{diag}(\sqrt p)\mathcal O_t^{(1,r,s)}\operatorname{diag}(\sqrt w)$ and pair it with an arbitrary matrix $Z$ satisfying $\|Z\|\le1$. If $S_i$ is the support of row $i$, then Cauchy--Schwarz first within each row and then over rows gives

$$
\begin{aligned}
|\langle Z,X\rangle|
&\le J_{r,s}\sum_i\sqrt{p_i}\sum_{j\in S_i}\sqrt{w_j}|Z_{ij}|\\
&\le J_{r,s}\sqrt R\left(\sum_jw_j\sum_i|Z_{ij}|^2\right)^{1/2}\\
&\le J_{r,s}\sqrt R.
\end{aligned}
$$

The last inequality uses that every column of an operator-norm-one matrix has Euclidean norm at most one and $\sum_jw_j=1$. Trace-norm duality therefore gives

$$
\left\|
\operatorname{diag}(\sqrt p)\,
\mathcal O_t^{(1,r,s)}\,
\operatorname{diag}(\sqrt w)
\right\|_1
\le
J_{r,s}\sqrt{R_{r,s}(t)}.
$$

Applying the same argument to $X^\top$ gives $J_{r,s}\sqrt{C_{r,s}(t)}$. Both estimates use the same arbitrary correlated laws $p,w$; neither replaces them with product laws or independently optimized link laws.

There are at most $q^{2|t|}$ row features and $q^{2(12-|t|)}$ column features. Because $\sum_{i,j}p_iw_j=1$, the weighted Frobenius norm is at most $J_{r,s}$. The weighted matrix has rank at most $q^{2\min(|t|,12-|t|)}$, so the trace-norm/Frobenius inequality gives

$$
q^{\min(|t|,12-|t|)}J_{r,s}.
$$

Consequently the exact rational squared sector coefficient is

$$
c_{t,r,s}^2
=
\min\left\{
q^{2\min(|t|,12-|t|)}J_{r,s}^2,\,
J_{r,s}^2R_{r,s}(t),\,
J_{r,s}^2C_{r,s}(t)
\right\}.
$$

This is the requested complete joint-chain energy formula. Each of the three terms bounds the same sector matrix for every $p,w$, so their minimum is valid. Sign parity makes every physical support tuple belong to exactly one record sector, so the four sector matrices have disjoint entry support and sum to the full physical occurrence matrix. Orthogonality is not required: the trace-norm triangle inequality gives

$$
\sup_{p,w}
\left\|
\operatorname{diag}(\sqrt p)\,
\mathcal O_t\,
\operatorname{diag}(\sqrt w)
\right\|_1
\le
\sum_{r,s\in\{1,3\}}c_{t,r,s}.
$$

## 5. Certified q64 values

All incidences and sector squares remain integers or Fraction values. Each displayed square root and each accumulated sector sum is rounded upward.

| canonical split $t$ | $(1,1)$ | $(1,3)$ | $(3,1)$ | $(3,3)$ | outward sum |
|---|---:|---:|---:|---:|---:|
| $(0,1,3,2)$ | $0.015641999334377113$ | $0.0024257518540529172$ | $0.2567034840429827$ | $0.06351573801306697$ | $0.33828697324447987$ |
| $(0,2,2,2)$ | $0.004095996697862384$ | $0.0017018123143495655$ | $0.06827900524431511$ | $0.044560149433685706$ | $0.11863696369021283$ |
| $(0,2,3,1)$ | $0.015641999334377113$ | $0.0013999466656331846$ | $0.2607473182634366$ | $0.036643960079605124$ | $0.3144332243430522$ |

The maximum coefficient is below one by more than $0.6617$.

Complementing a split exchanges $A_i$ and $B_i$ inside the same simple support $U_i$ and transposes the occurrence matrix after swapping $p$ and $w$; all physical masks and coefficients are unchanged. Reversal uses $M(S,T)=M(T^\top,S^\top)$, preserves simplicity, and reverses the record triple. Because the proof includes every $(r,s)$ sector, reversal is a bijection of the sector sum. These two operations carry the three canonical splits to all twelve registry entries without an omitted record orientation.

## 6. Independent physical validation

The regression tests/q64_joint_recovered_cubic_quintic_contraction.py does not use the production moment evaluator for its physical checks.

At $q=4$, it constructs all $4!2^4=384$ signed permutations and all physical cubic and quintic supports. It independently recovers all 15 feasible quintic shape pairs and computes the exact shared-quintic maxima in all four sectors:

$$
J_{11}^{(q=4)}=J_{13}^{(q=4)}=\frac1{192},
\qquad
J_{31}^{(q=4)}=J_{33}^{(q=4)}=\frac1{64}.
$$

At $q=8$, it directly sums all $8!$ permutations for representatives of every one of the 15 shape-table rows. These checks cover the $32$, $221$, $41$, $2111$, and $311$ geometries, the active and zero-xor cases, and all four record sectors. The rejected old-proof counterexample remains exactly

$$
\frac1{17920}.
$$

It is now safely contained in the $(3,1)$ joint shape bound. The $311/311$ plus shape attains the new $q=8$ $(3,3)$ bound exactly:

$$
J_{33}^{(q=8)}=\frac1{25088}.
$$

The artifact regression evaluates every exact rational sector square at 100-digit decimal precision and verifies that every committed float dominates the exact sum.

## 7. Registry decision and scope

The joint theorem repairs the twelve previously quarantined entries. The dependency-exact registry may therefore return to:

- 354 of 354 coefficient-one-dependent entries repaired;
- zero coefficient-one quarantine;
- 808 provisionally supported and 80 open before the separate dual-endpoint caveat;
- 796 conservatively supported and 92 open after withholding the 12 dual-endpoint entries.

This result repairs only the twelve $(1,3,5,3)$ entries. It does not alter the rejection of the old independent-maxima proof, does not address the later 80 residual entries, does not certify the full routing ledger, and does not prove the adaptive lift.
