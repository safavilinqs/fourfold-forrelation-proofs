# Opposite-endpoint mixed-orbit obstruction

Date: 2026-07-15

Status: exact physical witness. The fixed-orbit contraction does not extend
to arbitrary orbit mixtures, and the current independent profile/split
Perron ledger cannot pass at $N=1024$ without an additional improvement.
This obstructs the ledger, not the signed-permutation hard instance itself.

## 1. Why the fixed-orbit result was insufficient

For the critical split $(2,0,1,2)$ of $(3,1,1,5)$, a single aligned
translation orbit has coefficient

$$
0.00344731635
$$

at $q=32$. That is below the pair-specific uniform-coefficient target
$0.00624122$. Concavity, however, permits the diagonal probe laws to mix
selected-pair differences and complement-triple shapes. The pure-orbit
maximum is therefore only a lower diagnostic, not an upper certificate.

The mixed calculation below shows that this distinction is decisive.

## 2. Exact mixed-orbit Fourier reduction

Let $G=\mathbb F_2^{\log_2N}$. Translation twirling reduces the row law to
a distribution $p_{x,y}$ on the nonzero XORs of the selected cubic and
quintic endpoint pairs. It reduces the column law to a distribution
$r_\tau$ on translation orbits of the complement triple.

For a base pair of difference $x$, let $f^3_x(s)$ be its cubic endpoint
amplitude when the complement singleton is translated by $s$. For a
quintic pair of difference $y$ and triple representative $\tau$, define
$f^5_{y,\tau}(s)$ analogously. Put

$$
a_x(\alpha)=\widehat f^3_x(\alpha),
\qquad
g_{y,\tau}(\delta)
=(-1)^{\langle\operatorname{xor}\tau,\delta\rangle}
\widehat f^5_{y,\tau}(\delta).
\tag{2.1}
$$

The phase in $g$ makes it independent of the chosen translated
representative of $\tau$.

Use oriented pair representatives; assigning one quarter of the physical
row weight to each of the four duplicate orientations is an isometry. Walsh
transform the two pair translations and the free middle singleton. The
remaining matrix splits into $N^2$ blocks. The exact arbitrary-orbit
objective is

$$
\boxed{
\Phi(p,r)=N^{-3}\sum_{\alpha,\gamma\in G}
\|L_{\alpha,\gamma}(p,r)\|_1.}
\tag{2.2}
$$

Here $L_{\alpha,\gamma}$ has rows indexed by $(x,y)$, columns indexed by
$\tau$, and entries

$$
L_{\alpha,\gamma}[(x,y),\tau]
=\sqrt{p_{x,y}r_\tau}\,
a_x(\alpha)
g_{y,\tau}(\gamma\mathbin\oplus x).
\tag{2.3}
$$

For point masses $p,r$, (2.2) reduces exactly to the factored fixed-orbit
formula. At $q=4$, a genuinely mixed two-difference law also agrees with
the direct $2048$-by-$4096$ physical matrix to
$2.1\times10^{-15}$. Thus cross-orbit terms, not merely the diagonal
blocks, are retained.

## 3. An explicit physical mixed-orbit witness

Take $p$ uniform over every ordered pair of nonzero vertical differences.
Take $r$ uniform over translation orbits of triples contained in one hidden
column. Both are valid physical diagonal laws.

This family has a fourteen-class frequency reduction. Nonzero row
frequencies $\alpha_R,\gamma_R$ are classified by the bilinear value
$\langle\alpha_R,\gamma_R\rangle$, while their hidden-column frequencies
only distinguish zero from nonzero. Consequently (2.2) through $q=32$
requires only fourteen small nuclear norms.

The exact numerical values are:

| $q$ | $N$ | witness coefficient |
|---:|---:|---:|
| 4 | 16 | 0.067658234671 |
| 8 | 64 | 0.127869435555 |
| 16 | 256 | 0.085309750834 |
| 32 | 1024 | 0.039593955295 |

The nonmonotonic small-$q$ behavior is real. At the target size, the
witness is $6.34$ times the former pair-specific target. Therefore no
arbitrary-diagonal coefficient below $0.0395939552946$ can hold for the
critical split.

## 4. Forced failure of the old scalar ledger

The critical split, its complement, and the two reversed cuts all have the
same tester norm. Insert only those four forced coefficients into the
otherwise accepted partial ledger; set every other split of the two profiles
to zero. Reoptimizing the attenuation gives

$$
\beta_*=0.780899845855,
\qquad
F_{\rm forced}+2\epsilon_{\beta_*}
=0.334183454499.
\tag{4.1}
$$

Thus even this deliberately favorable ledger exceeds $1/3$ by

$$
\boxed{0.000850121166.}
\tag{4.2}
$$

No better upper bound for the other occurrence splits can repair (4.2).
The independent profile/split scalar ledger is quantitatively obstructed.

This is not a passive counterprotocol and does not prove that the attenuated
signed-permutation plant fails. A joint contraction can retain cancellation
between cuts or profiles, and an improvement elsewhere in the accepted
budget can create enough room for the forced norm.

## 5. Ranked repair targets

The overshoot is narrow. Reoptimizing $\beta$ after uniformly scaling one
accepted profile family gives the following sufficient reductions:

| repair family | required reduction |
|---|---:|
| adjacent double-cubic pair $(1,1,3,3)/(3,3,1,1)$ | 1.4273% |
| endpoint-cubic degree-six pair | 1.6366% |
| endpoint-quintic degree-eight pair | 2.1453% |
| all four triple-cubic profiles | 3.7660% |
| central double cubic $(1,3,3,1)$ | 4.4512% |
| double endpoint $(3,1,1,3)$ | 5.1303% |

The promise route has an even smaller analytic target. At $\beta_*$, the
proved one-sided failure bound is $0.00986238014$, coming from subgaussian
proxy

$$
V_{\beta_*}=0.001607508491.
$$

Holding the Perron side fixed, it suffices to replace this by

$$
V_{\rm target}=0.001592321190,
$$

a $0.9448\%$ proxy reduction. Equivalently, the two-hypothesis promise
loss need only fall by $0.000850122$, or $4.31\%$.

A seeded Monte Carlo diagnostic over $10^5$ independently attenuated
samples from 25 random exact plants found three promise failures, an empirical
rate $3\times10^{-5}$. This is not proof, but the gap from the current
$9.86\times10^{-3}$ bound makes sharper promise concentration the leading
repair project. The adjacent double-cubic $1.43\%$ target is the best
backup within the existing tester ledger.

## 6. Program decision

The next finite-size work should not seek a mixed-orbit upper below the
explicit witness. It should proceed in this order:

1. use the now-proved finite-tilt promise repair and rerank every still-open
   profile against its narrow remaining margin;
2. if the complete ledger fails, test whether the adjacent double-cubic
   triangle/incidence bound can be tightened by $1.43\%$ while retaining its
   physical law;
3. if that bounded repair fails, replace the independent scalar ledger by a
   joint cut/profile contraction or move the finite-size lead to the
   quadratic-bent exact plant, which removes promise loss entirely;
4. keep the repaired reverse-tree theorem as the asymptotic baseline and do
   not infer anything about adaptivity from this one-batch obstruction.

This is exactly the Round 3 pivot rule: record a quantitative failure, name
the amount that must be recovered, and compare bounded alternatives rather
than continuing to polish a falsified inequality.

## 7. Subsequent repair

The first target in Section 5 has been met rigorously.  Strict unimodality of
the centered Bernoulli log-MGF implies an exact Euclidean packing inequality
on the finite tilt interval used by the four-block reverse martingale.
Reoptimizing (4.1) with that promise theorem gives

$$
\beta=0.779512334639,\qquad
F_{\rm forced}+2\epsilon_\beta^{\rm E}
=0.332335953237,
$$

below \(1/3\) by \(0.000997380097\).

Thus this witness remains a valid falsification of the old concentration
ledger but is no longer a forced obstruction to every scalar completion.
The complete remaining high-degree ledger is still open.  See
euclidean_promise_concentration.md for the proof and exact scope.

A second physical witness subsequently exceeded this repaired endpoint
ledger.  Extending the same theorem blockwise with the global Kearns--Saul
branch repairs the combined two-witness diagnostic.  See
transposed_dominant_class_and_hybrid_repair.md; the current frontier should
be read there rather than inferred from the historical numbers above.

## 8. Reproduction

- searches/opposite_endpoint_mixed_orbit_q4.py implements (2.2), validates
  pure orbits, and searches $q=4$ mixed laws.
- searches/opposite_endpoint_vertical_mixture_witness.py evaluates the
  fourteen-class witness and the optimized forced ledger.
- searches/opposite_endpoint_repair_scorecard.py reproduces the ranked
  scalar repair percentages and promise-proxy target.
- searches/promise_tail_monte_carlo.py reproduces the explicitly
  nonrigorous seeded promise-tail diagnostic.
- tests/opposite_endpoint_mixed_orbit_obstruction.py protects the four
  witness values, (4.1)--(4.2), and the repaired forced-cut diagnostic.
- searches/occupation_compatible_sector_optimization.py accepts explicit
  profile/split coefficients for obstruction and future joint diagnostics.
