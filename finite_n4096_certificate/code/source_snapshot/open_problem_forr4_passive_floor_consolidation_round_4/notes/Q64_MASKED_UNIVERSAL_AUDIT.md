# q64 masked-universal audit

Date: 2026-07-17

Status: the universal coefficient-one proof remains withdrawn, and the old final recovered cubic--quintic product proof remains rejected. All 354 dependent entries now have independent actual-mask proofs: 54 quintic, 180 local-Walsh, 12 cubic-endpoint, 6 double-quintic endpoint, 12 double-quintic record-resolved, 38 four-cubic incidence, 12 cubic--septimic chain, 28 recovered cubic--quintic endpoint-row, and 12 replacement joint shared-quintic entries.

## The proof failure

The completed moment kernel

$$
K_0(S,T)=\mathbb E\overline{f_S}g_T
$$

is a unit-feature cross Gram and has arbitrary-diagonal coefficient at most one. The physical fixed-profile occurrence matrix is instead

$$
K(S,T)=\mathbf 1_{S\cap T=\varnothing}K_0(S,T)
$$

for every internally split block, with one such mask per split block. Schur multiplication by this mask is not contractive in trace class.

For example, if the completed kernel is constant and both sides are singletons, the physical matrix is $J-I$. Under uniform laws on $N$ labels,

$$
\left\|\frac{J-I}{N}\right\|_1
=2\left(1-\frac{1}{N}\right)>1.
$$

Therefore the written cross-Gram proof does not imply coefficient one. Any repair must use cancellation or geometry specific to the signed-permutation moments, or pay an explicit distinctness-mask factor.

## Exact dependency registry

The affected inventories are disjoint:

| dependency | entries |
|---|---:|
| universal septimic | 96 |
| universal three-/four-cubic | 14 |
| universal double-cubic | 24 |
| universal noncubic | 124 |
| recovered universal | 96 |
| total coefficient-one dependent | 354 |
| repaired by actual-mask quintic factorizations | 54 |
| repaired by actual-mask local-Walsh factorizations | 180 |
| repaired by physical cubic-endpoint factorizations | 12 |
| repaired by double-quintic endpoint factorizations | 6 |
| repaired by double-quintic record-resolved rows | 12 |
| repaired by four-cubic record-sector incidences | 38 |
| repaired by cubic--septimic record-compatible chains | 12 |
| repaired by recovered cubic--quintic endpoint rows | 28 |
| repaired by joint recovered cubic--quintic shapes | 12 |
| total independently repaired | 354 |
| remaining quarantined | 0 |

Replacing every coefficient-one-dependent claim by its independent repair, adding the separate final residual theorem, and applying the independently audited dual-endpoint theorem gives 888 certified one-batch entries and zero open entries. The generated artifact records one explicit status row for every one of the 888 entries.

The registry module's complete floating cross-check, after inserting the final-80 coefficients, is $0.239911667770$. This optimizer output is not the accepted numerical certificate. The separate dependency-exact outward ledger fixes $\beta=19/25$ and certifies total $0.268858135059926<1/3-10^{-3}$ by directed rounding and Collatz--Wielandt. See `Q64_COMPLETE_OUTWARD_LEDGER.md`.

The repair uses exact rational singleton--quintic fixed-slice energies, a two-sided record-one tail, a separated-endpoint tensor factorization, and an exact central Walsh-chain completion with both physical masks restored by inclusion--exclusion. Its 54 coefficients are at most $0.646562122163$ and never invoke the invalid masked Gram inference. See `Q64_MASKED_QUINTIC_SLICE_REPAIR.md`.

The second repair covers 180 entries having an internal singleton. This
includes higher--singleton--higher chains and both orientations of adjacent
singleton pairs. The local chain supplies an exact $q^{-1}$ and
integer-ceiled inclusion--exclusion factors restore every physical mask. Its
maximum coefficient is $54/64=0.84375$. See
`Q64_MASKED_LOCAL_WALSH_REPAIR.md`.

The third repair uses an exact physical split-cubic endpoint factor and
rational upper bounds for every remaining mask. It closes 12 entries with
coefficient $0.930264114402$. See
`Q64_MASKED_CUBIC_ENDPOINT_REPAIR.md`.

The fourth repair composes two complete physical quintic endpoint
factorizations while retaining the central quintic--quintic moment kernel as
one cross Gram. Both endpoint occurrence masks and both singleton-side
orientation factors are included in their exact slices. It closes six
$(1,5,5,1)$ entries, with worst coefficient
$\sqrt{1023/1024}=0.999511599482$. See
`Q64_MASKED_DOUBLE_QUINTIC_ENDPOINT_REPAIR.md`.

The fifth repair splits every $(3,3,3,3)$ physical kernel into eight disjoint
record triples. Exact record-one and record-three link maxima, together with
five exact cubic completion-incidence families, bound complete physical rows
and columns without removing a mask. It closes all 38 residual four-cubic
entries with worst coefficient $0.00894228260681$. See
`Q64_MASKED_FOUR_CUBIC_INCIDENCE_REPAIR.md`.

The sixth repair retains the cubic endpoint record and the two compatible
degree-seven bidegree sectors. Exact physical completion incidences keep the
septimic mask, while the fixed-one cubic endpoint split prevents its large
record-three energy from being charged at the record-one central moment. For
the far-singleton rows, a zero-active/injective-Walsh split bounds the
dangerous record-one endpoint energy: at most $5(q-1)$ completions lack a
Walsh residual, and every other completion gains $1/(q-1)^2$. It closes all
twelve residual cubic--septimic entries with coefficients from
$0.0884219637995$ to $0.538626546351$. See
`Q64_MASKED_CUBIC_SEPTIMIC_CHAIN_REPAIR.md`.

The seventh repair now retains only the 28 recovered endpoint-row entries. Exact singleton--cubic/quintic endpoint energies followed by physical block incidences give seven orbit coefficients from $0.0162888571820$ to $0.703615181088$. The former 12-entry four-sector chain is rejected: its separate cubic--quintic link maxima do not have the claimed physical scope, and an exact sector-$(1,3,1)$ chain violates the product maximum. See `Q64_MASKED_RECOVERED_CUBIC_QUINTIC_INCIDENCE_REPAIR.md` and `Q64_RECOVERED_CUBIC_QUINTIC_INDEPENDENT_AUDIT.md`.

The eighth repair resolves the central odd record in all five remaining $(1,5,5,1)$ orbits. Exact physical endpoint sums retain both quintic masks, and the central record-$r$ permanent is bounded by $1/\binom qr$. Exact affine-shape sums close ten entries; a parity bound closes the final variable-four pair. The twelve coefficients range from $0.000260884928748$ to $0.198263998705$, removing the entire double-quintic family from quarantine. The regression independently checks all $q=4$ rows, selected $q=8$ rows, and the full finite $q=64$ pair/triple affine-shape enumeration. See `Q64_MASKED_DOUBLE_QUINTIC_RECORD_REPAIR.md`.

The ninth repair replaces the rejected recovered-chain argument. It enumerates all 15 feasible simple quintic row/column multiplicity pairs and bounds the endpoint--cubic--quintic factor and quintic--cubic factor within the same physical quintic shape. Exact row, column, and rank energies give the three canonical coefficients $0.338286973244$, $0.118636963690$, and $0.314433224343$. Full $q=4$ physical enumeration and direct $q=8$ permutation checks of every shape row pass. See Q64_JOINT_RECOVERED_CUBIC_QUINTIC_CONTRACTION.md.

The separate final residual theorem closes the later 80 entries. Forty-eight use the accepted local-Walsh mechanism with the sole higher-block mask restored explicitly. Thirty-two use exact terminal-high, double-endpoint, endpoint--cubic--high, and endpoint--high--cubic scalar features. Their maximum coefficient is $0.447299757774$. Independent $q=4$ rows, direct $q=8$ permutation sums, and exhaustive generation of the exceptional zero-active septimic family pass. See Q64_FINAL_RESIDUAL_CHAIN_CONTRACTION.md.

## Exact physical $q=2$ screen

The audit constructs the complete physical signed-permutation occurrence matrix for profile $(3,3,3,3)$ in three representative split geometries. Every cross-cut within-block distinctness condition is imposed before evaluating the exact chain moment. Eight times each matrix has entries in $\{-1,0,1\}$.

| split | split cubics | matrix | nonzero | uniform norm | tangent upper |
|---|---:|---:|---:|---:|---:|
| $(0,1,2,3)$ | 2 | $96\times96$ | 2,304 | $0.416666667$ | $0.416666667$ |
| $(1,1,1,3)$ | 3 | $256\times216$ | 6,912 | $0.519214888$ | $0.584366971$ |
| $(1,1,2,2)$ | 4 | $576\times576$ | 20,736 | $0.610633422$ | $0.680474822$ |

The weighted trace norm is jointly concave in the two diagonal laws. A
polar-factor tangent at the uniform laws gives the displayed floating global
upper estimate over both probability simplexes. A $10^{-10}$ numerical
allowance is included, but this is not an outward-rounded interval
certificate. Thus coefficient one survives this numerical screen on the
complete exact $q=2$ cubic representatives.

This is not evidence that coefficient one holds at $q=64$. It only rules out the easiest small-order cubic counterexample and validates that the exact test includes the masks correctly.

## Generic formal repair fails numerically

Inclusion--exclusion gives a rigorous Schur factor for disjoint $k$- and $l$-subsets,

$$
\Delta_{k,l}
=\sum_j\sqrt{\binom{k}{j}\binom{l}{j}}.
$$

Composing one such factor for every split block repairs the missing logical step without using cancellation. Across the 354 affected entries, the resulting coefficients range from $3$ to $33.9705627485$. Their q64 routing total is $3.40344112205$, far above $1/3$. The generic repair is therefore formally valid but useless for the paper target.

## Decision

The foundational masked-contraction, final-residual, dual-endpoint, and outward-ledger gates are passed. The dependency-valid count is 354 of 354 affected entries repaired, all 80 later residual entries proved, and all 12 dual-endpoint entries independently certified. The registry is 888 certified and zero open. The later direct-sum tree-frontier theorem lifts the outward total through unrestricted classical feed-forward with multiplier one; see `Q64_ADAPTIVE_TREE_FRONTIER_THEOREM.md`.

## Acceptance criteria met by the joint repair

The masked-contraction project passes only if all of the following hold:

1. The occurrence kernel includes every cross-cut distinctness mask before
   any Gram, Schur, symmetry, or completion step.
2. The 12 recovered chain entries are covered by one joint-link theorem whose 15 structural templates have a machine-checked union and simple-support feasibility classification.
3. The coefficient bound is uniform over arbitrary correlated row and column
   diagonal laws. Uniform-law, product-law, or invariant-law evidence is not
   enough.
4. Every analytic mask or cancellation factor is exact or outward rounded.
   Floating small-$q$ optimization is a falsification screen, not the proof.
5. A falsification is decisive only if an exact or interval-certified
   compatible physical law gives a masked coefficient strictly above one.
6. After either outcome, one command regenerates all 888 registry rows and
   rejects any routing artifact that depends on a quarantined entry.
7. Adaptive work resumes only if the rebuilt, fully proved one-batch ledger is
   at most $1/3-10^{-3}$. Otherwise pivot the contraction or hard instance.

These criteria serve the paper goal: they prevent another attractive routing
total from outrunning its physical operator proof, while preserving a fast
exit if the signed-permutation witness cannot yield a credible finite-size
separation.

Reproduce with:

    /opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_masked_universal_audit.py --write-artifact
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_masked_universal_audit.py
