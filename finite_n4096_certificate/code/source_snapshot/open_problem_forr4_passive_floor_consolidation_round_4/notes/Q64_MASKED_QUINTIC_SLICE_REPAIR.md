# The $q=64$ masked quintic slice repair

Date: 2026-07-16

Status: proved arbitrary-correlated-diagonal one-batch theorem for 54 entries
of the actual distinctness-masked occurrence matrix. This is a partial repair
of the 354-entry quarantine, not a proof of the universal coefficient-one
lemma.

## The row-energy lemma

Let $K=(K_{ij})$ be the complete physical occurrence matrix after every
cross-cut distinctness mask has been imposed. If

$$
\sup_i\sum_j |K_{ij}|^2\le R,
$$

then for arbitrary probability laws $p$ and $q$ on its rows and columns,

$$
\left\|D_p^{1/2}KD_q^{1/2}\right\|_1\le\sqrt R.
$$

Indeed, take $u_i=(K_{ij})_j$ and $v_j=e_j$ in the column Hilbert space, so
$K_{ij}=\langle u_i,v_j\rangle$. The associated weighted feature maps have
Hilbert--Schmidt norms at most $\sqrt R$ and $1$. Schatten Hölder gives the
claim. Applying the same argument to $K^*$ gives the better of the maximum
row and maximum column energies.

This factorization never completes the masked matrix to an unmasked Gram
kernel. For a fixed occurrence row, every surviving column is literally a
disjoint completion of the fixed partial supports.

## Exact quintic slices

Let $E_k$ be the accepted exact squared-moment slice for a degree-five support
next to a singleton, when a fixed $k$-cell subset of the quintic is retained.
The singleton coordinate is optimized. At $q=64$,

| $k$ | $E_k$ |
|---:|---:|
| 0 | $839885$ |
| 1 | $1025.250244140625$ |
| 2 | $41.37997581845239$ |
| 3 | $1.4538457961309523$ |
| 4 | $0.9990234375$ |
| 5 | $0.000244140625$ |

The theorem module evaluates the defining combinatorial formulas as exact
rational numbers. Decimal values here are only for scale. A singleton--
singleton link has squared magnitude $q^{-2}$. A fixed singleton--quintic
entry has squared magnitude at most $q^{-2}$, while summing over a variable
singleton costs at most one.

### Two-sided quintic tail

An internal quintic in the profiles $(1,5,1,5)$ and $(5,1,5,1)$ has record
one on both adjacent axes. The exact record-one link formula shows that a
degree-five support can attain singleton amplitude $q^{-1}$ only in the
column-occupancy types $(5)$ or $(4,1)$ with zero xor on the four-cell
column. Type $(5)$ has row record five. In the exceptional type $(4,1)$,
the singleton row either belongs to the four-cell column, leaving row record
three, or does not, leaving row record five. Neither can have row record one.

Every bi-record-one quintic therefore has singleton-link entries at most
$1/[q(q-1)]$. Summing their squares over the $q^2$ singleton supports gives
the uniform tail bound

$$
T_5^{\rm bi}\le {1\over(q-1)^2}.
$$

The regression independently enumerates all 1008 bi-record-one quintics at
$q=4$ and evaluates every singleton moment as an exact fraction. Every tail
energy equals $1/9$.

## Chained physical row bounds

Write $s=(s_1,s_2,s_3,s_4)$ for the occurrence split. For profile
$(5,1,1,5)$, the central singleton--singleton link separates the two
quintic slices and gives

$$
R_s=q^{2(1-s_2-s_3)}E_{s_1}E_{s_4}.
$$

For profile $(1,5,1,5)$, sum the endpoint quintic slice first. On the first
two links, two fixed singletons contribute $q^{-2}$. With one variable
singleton, the bi-record tail contributes $(q-1)^{-2}$. With two variable
singletons, summing the other singleton and its quintic slice costs $q^2$,
while the retained bi-record tail contributes $(q-1)^{-2}$. Thus

$$
R_s=E_{s_2}E_{s_4}
\begin{cases}
q^2/(q-1)^2,&s_1+s_3=0,\\
1/(q-1)^2,&s_1+s_3=1,\\
q^{-2},&s_1+s_3=2.
\end{cases}
$$

The profile $(5,1,5,1)$ follows by path reversal. In every case the physical
coefficient is bounded by

$$
\gamma_s\le
\sqrt{\min\{R_s,R_{d-s}\}},
$$

where $d-s$ is the complementary occurrence split. These are direct sums
over disjoint completions; no universal unmasked coefficient is used.

## Separated-endpoint tensor factorization

The remaining profile $(5,1,1,5)$ admits a second masked factorization. For
one quintic--singleton endpoint link, let the quintic occurrence size be $k$
and let $s\in\{0,1\}$ say whether the singleton lies on the row. Its maximum
physical row and column energies are bounded by

$$
q^{2(1-s)}E_k
\quad\text{and}\quad
q^{2s}E_{5-k},
$$

respectively. Factoring through complete rows or complete columns therefore
gives squared feature-norm product

$$
Q_{k,s}=\min\{q^{2(1-s)}E_k,q^{2s}E_{5-k}\}.
$$

Tensoring the two endpoint feature factorizations preserves arbitrary
correlations between the two endpoint occurrence laws. If the central
singletons lie on the same side, their link is a row or column phase of
magnitude $q^{-1}$. If they lie on opposite sides, the link is a unit-feature
Gram factor and its Schur product costs one. Hence

$$
T_s=Q_{s_1,s_2}Q_{s_4,s_3}
\begin{cases}
q^{-2},&s_2=s_3,\\
1,&s_2\ne s_3.
\end{cases}
$$

For this profile the final squared coefficient is bounded by
$\min\{R_s,R_{d-s},T_s\}$. This proves six entries that the scalar row-energy
bound alone misses.

## Result

All 54 quarantined entries in the three profiles above have a displayed
squared coefficient at most one. Their eight exact squared
coefficients are

$$
{1023\over4194304},\quad {81\over4096},\quad {81\over256},\quad
{55638713045\over3670705963008},\quad
{1432003925\over5549064192},\quad
{4199425\over16257024},\quad
{262577446975\over721554505728},\quad
{791806327225\over1894080577536}.
$$

Thus their outward-rounded coefficients lie between
$0.015617368742$ and $0.646562122163$. The repair is closed under occurrence
complement and path reversal. It raises the safely supported registry from
454 to 508 entries and leaves 300 of the original coefficient-one-dependent
entries quarantined. With the separate 12-entry dual-endpoint caveat withheld,
the conservative supported count is 496.

This result also clarifies why the rejected record-one shortcut failed. The
proof does not assign norm one to decorated links. It retains the exact
quintic fixed-slice energies and the two-sided record constraint through the
full three-link path, then uses the physical complement relation to choose the
lower-energy orientation. The final six candidates require the central
Walsh-chain completion below.

## Central Walsh-chain completion

For every record-one quintic support $S$ and singleton $z$, the exact endpoint
formula is

$$
M_{51}(S,z)=v_5(S)H_N(\operatorname{xor}S,z),\qquad |v_5(S)|\le1.
$$

Consequently, pointwise in the shared singleton,

$$
M_{51}(S,a)H_N(a,b)
={v_5(S)\over q}H_N(a,\operatorname{xor}S\mathbin\oplus b).
$$

Duplicate compression makes the residual normalized Walsh matrix a
coefficient-one Schur factor, leaving the explicit factor $q^{-1}$.

Restore the occurrence mask at each split quintic by inclusion--exclusion.
For an $r$-set and an $s$-set this gives

$$
\gamma_{r,s}\le
\sum_{t=0}^{\min(r,s)}\sqrt{{r\choose t}{s\choose t}}.
$$

Thus $gamma_{1,4}\le3$ and
$\gamma_{2,3}\le1+\sqrt6+\sqrt3<6$. Completing both unmasked endpoint
kernels and then restoring both physical masks gives

$$
c_s\le {\Gamma_{s_1}\Gamma_{s_4}\over q},\qquad
\Gamma_1=\Gamma_4=3,\quad \Gamma_2=\Gamma_3=6.
$$

At $q=64$, the two new coefficient bounds are $9/64$ and $9/16$. Both are
below one. This repairs the final six entries while explicitly paying for
the masks omitted by the rejected universal lemma. The regression checks the
endpoint identity for all 4368 quintic supports and all 16 singleton supports
in both orientations at $q=4$, totaling 139776 exact comparisons.

## Residual exact-moment screens

The last six entries form three complement pairs. Their representative
splits are

$$
(1,1,0,4),\qquad (2,0,1,3),\qquad (2,1,0,3).
$$

The full-group $q=4$ canonical pure-law screen gives coefficients
$0.015625$, $0.010416666667$, and $0.0226153926$, respectively. A targeted
$q=8$ screen evaluates ten deterministic 64-element translation-subspace
orbits for each representative using exact signed-permutation moments and
every physical distinctness mask. The largest optimized nuclear lower
witnesses are, in the same order,

$$
6.7250649\mathbin{\times}10^{-5},\qquad
7.53097278296\mathbin{\times}10^{-5},\qquad
6.0607366\mathbin{\times}10^{-5}.
$$

No screened coefficient exceeds one. These numbers are lower witnesses and
orbit-restricted tangent values, not arbitrary-law upper bounds. They did not
repair the six entries by themselves; the Walsh-chain completion above is the
arbitrary-law proof.

## Reproduction

Run:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_masked_quintic_slice_repair.py --write-artifact
/opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_masked_quintic_slice_repair.py
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q8_masked_separated_quintic_residual_screen.py --output artifacts/q8_masked_separated_quintic_residual_screen.json
/opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q8_masked_separated_quintic_residual_screen.py
```

The deterministic artifacts are
`artifacts/q64_masked_quintic_slice_repair.json` and
`artifacts/q8_masked_separated_quintic_residual_screen.json`.
