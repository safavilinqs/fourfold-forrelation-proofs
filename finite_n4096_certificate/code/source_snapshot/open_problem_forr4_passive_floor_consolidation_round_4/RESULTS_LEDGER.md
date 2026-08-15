# Round 4 results ledger

Date initialized: 2026-07-15

This file tracks evidence for the paper-ready finite-size separation.  Update
it only with scoped, reproducible evidence.

## Paper headline target

| field | current value | status |
|---|---|---|
| chosen size $N_*$ | $4096$ | mathematically certified |
| sign modes $M=4N_*$ | $16{,}384$ | exact resource count |
| passive hard dose six | excluded in the classically adaptive model | proved plus certified numeric |
| adaptive margin to $1/3-10^{-3}$ | $0.063475198273408$ | certified numeric |
| adaptive and conditioning scope | multiplier-one tree lift; promise paid once | proved |
| active hard dose | six; error $81/256$ | proved protocol |
| experimental feasibility | ideal resources complete; present-day credibility not certified | decided negative on reviewed evidence |

The q64 registry has all 354 actual-mask repairs, all 80 later residual entries, and all 12 independently audited dual-endpoint entries proved. It supports all 888 balanced entries with zero open. The old universal-lemma route remains withdrawn. The dependency-exact outward total is $0.268858135059926$. The direct-sum adaptive frontier has multiplier one, so passive hard dose six is excluded at $N=4096$ in the declared classically adaptive model. The theorem package is complete; the hardware screen finds no reviewed demonstration of the required $4096$-dimensional coherent transform at combined contrast above $0.904294855157$.

## Foundational masked-universal audit

- The five disjoint affected inventories have sizes 96, 14, 24, 124, and 96,
  totaling 354 entries.
- The abstract completed-cross-Gram argument is insufficient: a constant
  singleton--singleton kernel becomes (J-I), whose uniformly weighted trace
  norm is (2(1-1/N)>1).
- Complete exact physical matrices for the three cubic split geometries at
  (q=2) have numerical arbitrary-law tangent uppers (0.416667),
  (0.584367), and (0.680475). Coefficient one survives this screen but is
  not proved at q64.
- The generic inclusion--exclusion mask repair is rigorous, with affected
  coefficients from 3 to (33.9705627485), but its q64 routing total is
  (3.40344112205), so it fails the paper target.
- Evidence: `notes/Q64_MASKED_UNIVERSAL_AUDIT.md` and
  `artifacts/q64_masked_universal_audit.json`, including all 888 registry rows.

### Actual-mask quintic slice repair

- A direct row-vector factorization bounds arbitrary correlated diagonal laws
  by the square root of the maximum complete physical row or column energy.
- Exact rational degree-five singleton slices control that energy without
  removing any cross-cut distinctness mask.
- The two-sided record-one constraint excludes the exceptional quintics with
  singleton-tail energy one, improving the internal tail to $(q-1)^{-2}$.
- A tensor product of the two endpoint feature factorizations and the central
  singleton Gram/phase factor closes six more separated-endpoint entries.
- The final six use the exact quintic endpoint character, central Walsh-chain
  collapse, and explicit inclusion--exclusion factors for both physical masks.
- This proves 54 entries across $(1,5,1,5)$, $(5,1,1,5)$, and $(5,1,5,1)$.
  Their coefficients range from $0.015617368742$ to $0.646562122163$.
- The original 354-entry dependency set is unchanged as provenance, but only
  none of those entries remain quarantined after all eight repairs.
- Evidence: `notes/Q64_MASKED_QUINTIC_SLICE_REPAIR.md` and
  `artifacts/q64_masked_quintic_slice_repair.json`.

### Actual-mask local-Walsh repair

- An internal singleton between two higher odd-degree blocks collapses its
  two endpoint characters to a normalized Walsh factor and an explicit
  $q^{-1}$.
- A same-side adjacent singleton pair supplies the same scalar gain directly.
- Integer-ceiled inclusion--exclusion factors restore every split-block mask.
  The maximum coefficient is $54/64=0.84375$.
- This proves 180 entries: 80 higher--singleton--higher cases, 52 same-side
  singleton-pair cases, and 48 opposite-side singleton-pair chains.
  Together the actual-mask theorems prove 234 of 354 affected entries, leaving
  120 quarantined.
- Evidence: `notes/Q64_MASKED_LOCAL_WALSH_REPAIR.md` and
  `artifacts/q64_masked_local_walsh_repair.json`.

### Actual-mask cubic-endpoint repair

- The exact physical split-cubic endpoint squared factor is
  $(q^2-2q+2)/(q^2(q-1))$.
- Rational inclusion--exclusion bounds cost at most $15/2$ for all remaining
  masks, giving coefficient $0.930264114402$.
- This proves four cubic--septimic and eight recovered
  cubic/cubic/quintic entries, leaving 108 quarantined.
- Evidence: `notes/Q64_MASKED_CUBIC_ENDPOINT_REPAIR.md` and
  `artifacts/q64_masked_cubic_endpoint_repair.json`.

### Actual-mask double-quintic endpoint repair

- The exact degree-five physical endpoint slice gives the squared factor
  $Q_{k,s}=\min\{q^{2(1-s)}E_k,q^{2s}E_{5-k}\}$, where $k$ counts quintic
  cells on the row side and $s$ records the adjacent singleton side.
- Composing two endpoint factors while retaining the middle
  quintic--quintic moment kernel as one cross Gram preserves arbitrary
  correlated diagonal laws.
- This proves six of the 18 residual $(1,5,5,1)$ entries. The worst outward
  coefficient is $\sqrt{1023/1024}=0.999511599482$; the other 12 stay
  quarantined.
- Evidence: `notes/Q64_MASKED_DOUBLE_QUINTIC_ENDPOINT_REPAIR.md` and
  `artifacts/q64_masked_double_quintic_endpoint_repair.json`.

### Actual-mask double-quintic record repair

- Five orbits containing all twelve residual $(1,5,5,1)$ entries are split by the central odd record $r\in\{1,3,5\}$.
- Exact endpoint row sums keep the pair--triple distinctness mask in the first quintic and the triple--pair mask in the second.
- The central record-$r$ permanent is bounded by $1/\binom qr$, giving exact squared complete-row coefficient $5901977909483/1291082725544951808$.
- Exact affine-shape sums close ten entries; the parity bounds $B_1\le2q$ and $B_3,B_5\le q^2$ close the last variable-four pair. The outward coefficients range from $0.000260884928748$ to $0.198263998705$.
- Evidence: `notes/Q64_MASKED_DOUBLE_QUINTIC_RECORD_REPAIR.md` and `artifacts/q64_masked_double_quintic_record_repair.json`.

### Actual-mask four-cubic incidence repair

- The physical $(3,3,3,3)$ kernel is split into its eight disjoint
  record-one/record-three triples.
- Exact one-link moment maxima are combined with completion incidences for
  endpoint record-one, endpoint record-three, middle $(1,1)$, middle
  $(1,3)$, and middle $(3,3)$ cubic families.
- Every column completion is the disjoint complement of a fixed row partial
  support inside one of those families, so the physical masks remain present.
- This proves all 38 residual four-cubic entries with coefficients between
  $0.00224319523802$ and $0.00894228260681$.
- Evidence: `notes/Q64_MASKED_FOUR_CUBIC_INCIDENCE_REPAIR.md` and
  `artifacts/q64_masked_four_cubic_incidence_repair.json`.

### Actual-mask cubic--septimic chain repair

- The two singleton endpoints force record one, and the central
  cubic--septimic link is split into records one and three.
- Exact degree-seven bidegree completion incidences retain the septimic
  occurrence mask; the complete cubic endpoint slice retains the cubic mask.
- Splitting the fixed-one cubic endpoint energy by central record aligns its
  large record-three part with the smaller $1/\binom q3$ link moment.
- A zero-active/injective-Walsh split controls the far-singleton rows: at
  most $5(q-1)$ completions lack a Walsh residual and every other completion
  gains $1/(q-1)^2$ in squared energy.
- This proves all twelve residual cubic--septimic entries in three orbits,
  with coefficients from $0.0884219637995$ to $0.538626546351$.
- Evidence: `notes/Q64_MASKED_CUBIC_SEPTIMIC_CHAIN_REPAIR.md` and
  `artifacts/q64_masked_cubic_septimic_chain_repair.json`.

### Actual-mask recovered cubic--quintic incidence repair

- The retained profiles are $(1,3,3,5)$ and $(1,5,3,3)$, with four common-record sectors $(1,r,s)$ for $r,s\in\{1,3\}$.
- Exact equality-shape completion incidences retain the occurrence masks in
  the quintic and both cubic blocks.
- Complete endpoint rows close 28 recovered entries in seven orbits, with coefficients from $0.0162888571820$ to $0.703615181088$.
- The former twelve-entry four-sector chain is rejected: an exact $q=8$ sector-$(1,3,1)$ physical entry is $1/17920$, above the claimed maximum $1/25088$.
- Evidence: `notes/Q64_RECOVERED_CUBIC_QUINTIC_INDEPENDENT_AUDIT.md`, `notes/Q64_MASKED_RECOVERED_CUBIC_QUINTIC_INCIDENCE_REPAIR.md`, and `artifacts/q64_masked_recovered_cubic_quintic_incidence_repair.json`.

### Final residual chain contraction

- The final residual inventory contains 80 entries in 20 complement/reversal orbits, each with one split quintic or septimic block and no split cubic.
- Forty-eight entries in 12 orbits inherit the accepted local-Walsh theorem and have coefficient at most $12/64=0.1875$.
- Thirty-two entries in eight orbits use four scalar-feature families: terminal split higher block, double endpoint, endpoint--cubic--higher, and endpoint--higher--cubic.
- Exact physical incidences retain the sole higher-block mask. The exceptional record-$(3,1)$ septimic family has $15119032320$ zero-active supports and exact fixed-cell incidence $25838190$ at $q=64$.
- The eight chain-orbit coefficients range from $0.018839926256$ to $0.447299757774$.
- Independent checks use all $4!2^4$ signed permutations at $q=4$, direct $q=8$ permutation sums for the exceptional link geometries, and exhaustive generation of all 18,816 zero-active septimic supports at $q=8$.
- Evidence: `notes/Q64_FINAL_RESIDUAL_CHAIN_CONTRACTION.md` and `artifacts/q64_final_residual_chain_contraction.json`.

### Masked translation reduction

- Exact signed-permutation link covariance was proved and checked on 400 link
  and 400 chain/occurrence identities at each of $q=4$ and $q=8$.
- The translation character separates across disjoint row and column
  occurrence supports, so the physical masked matrix changes only by
  permutations and diagonal signs.
- Joint concavity therefore permits arbitrary diagonal laws to be twirled
  without decreasing the coefficient. This is a global reduction, not an
  invariant-law assumption.
- The original 354 affected entries reduce to 97 complement/reversal templates; all 354 entries now have independent actual-mask proofs.
- Exact cocycle row reduction further classifies them into normalized
  projective commutator ranks $0,4,8$, with $21,26,50$ templates. The same
  classification holds at $q=4,8,64$, and every complementary row/column
  sector has the identical cocycle, not merely matching rank.
- The exact full-group Clifford formula reduces every mixed orbit-shape law
  to matrix-valued twisted Fourier symbols while retaining all cross-shape
  blocks.
- One complete $q=4$ orbit for all 97 templates has maximum
  $0.0836986669474$. A 210-shape focused search raises the pure maximum to
  $0.176776695297$.
- Selected complete-group mixed simplexes for normalized ranks $0,4,8$ have
  tangent uppers $0.0223813195039$, $0.166741363473$, and
  $0.176776695297$, respectively.
- Exact $q=4$ translation-subspace screens cover all 97; the maximum tested
  pure-orbit coefficient is $0.0220970869121$.
- Selected $q=8$ pure-orbit and three-shape mixture screens have maxima
  $6.01796939448\times10^{-5}$ and $0.000279017854455$, respectively.
- Shape mixing is visible but remains far below one. The missing theorem is a
  uniform bound over every physical shape-indexed symbol matrix.
- A proposed 118-entry record-one singleton-anchor shortcut is invalid:
  record-one operator norm one applies only before even decorations are
  restored. Exact $q=4$ full-sector norms are $\sqrt3$ for degrees $1,3$ and
  $4$ for degrees $3,3$. The next Bessel proof must retain two-sided
  decoration compatibility or fixed-slice energies.
- Evidence: `notes/MASKED_TRANSLATION_REDUCTION.md`, five exact screen
  artifacts, and three exact cocycle-inventory artifacts.

## Historical q32 compatible-law baseline

- All 16 leading unresolved families are now represented simultaneously by
  one legal configuration-level law.
- The selected order-four law uses 180 basis configurations across the 30
  Perron-support occupation states and activates 26 of the 64 leading
  profile-split entries.
- Its attenuated joint/separate trace-norm ratio is
  $0.746162307493$. The exact-moment order-32 embedding has ratio
  $0.854097529510$.
- The calculation proves that cancellation matters and supplies an exact
  order-32 moment engine. It is not a theorem coefficient or a decision on
  the hard instance.
- Evidence: notes/JOINT_PHYSICAL_LAW_DIAGNOSTIC.md and
  artifacts/joint_impact_sparse_diagnostic.json.

### Closed native row-translation screen

- One base configuration for each of the 30 Perron-support occupation states,
  averaged over 32 common row translations, gives a legal native-$q=32$
  law on 960 configurations.
- Among 50,112 combinatorially matching pairs, exact parity leaves 1,984
  nonzero pairs in only 6 frontier orbits.
- The separate attenuated mass is $4.79417\times10^{-13}$, only
  $1.70\times10^{-11}$ of the current 51-orbit frontier contribution.
- Decision: close common-row-translation laws as a lower-witness route. This
  is not an arbitrary-law upper bound.
- Evidence: notes/SHARED_FRONTIER_ROW_ORBIT_SCREEN.md and
  artifacts/shared_frontier_row_orbit.json.

### Closed leading-orbit non-invariant screen

- A valid vertical--horizontal physical law was checked against the
  unreduced order-four formula before use at order 32.
- It improves the invariant coefficient from $0.0395996495754$ to
  $0.0396118487001$.
- It remains $0.00185046959644$ below the
  $0.0414623182965$ scalar kill gate.
- Decision: stop local witness tuning and proceed to the global shared
  contraction. This does not exclude all arbitrary physical laws.
- Evidence: notes/OPPOSITE_ENDPOINT_TWO_AXIS_SCREEN.md and
  artifacts/opposite_endpoint_two_axis_screen.json.

## Historical initial Round 4 decision

The q32 independent orbit-by-orbit scalar ledger was rejected as the lead
project. Its $0.0002007278$ diagnostic margin was too small to justify 224
separate orbit bounds. The initial retain-or-pivot test used one compatible
physical law over this impact frontier:

| frontier | unresolved Perron contribution | role |
|---:|---:|---|
| leading orbit | 13.6% | included diagnostic |
| leading 16 orbits | 51.1% | first lower-witness stage |
| leading 51 orbits | 90.5% | shared-contraction target |

These percentages remain routing diagnostics, not theorem bounds. A
compatible lower-witness vector reaching $1/3$ would have closed the scalar
proof architecture, not the hard instance itself. This decision framework
led to the live q64 route described below.

### Exact shared-frontier factorization

- The leading 51 orbits contain 198 centrally balanced entries but only six
  unordered profile patterns.
- Their compatible dose-six action uses 125 occupation states and 241
  undirected edges in eight connected components.
- The 100 degree-ten entries each have five compatible occupation terms; the
  98 degree-twelve entries each have one.
- The `(5,3,1,1)` family carries 62.1 percent of the frontier's current
  Perron contribution; adding `(5,3,3,1)` raises coverage to 81.8 percent.
- This is an exact structural reduction, not an arbitrary-law coefficient
  theorem. It records the superseded q32 shared-operator target.
- Evidence: notes/SHARED_FRONTIER_STRUCTURE.md and
  artifacts/shared_frontier_structure.json.

## Inherited mathematical boundary

| Question | Best accepted result | Status | Missing evidence |
|---|---|---|---|
| Passive asymptotic lower | $\Omega(N^{1/12})$ | proved | better physical contraction or hard instance |
| Passive lower at $N=4096$ | hard dose $>6$ | proved plus certified numeric | coherent inter-batch quantum memory remains outside the declared model |
| Active upper | hard dose six, error $81/256$ | proved | current-platform implementation not certified |
| Current hard instance | attenuated signed-permutation q64 route | certified lower-bound witness | experimental credibility is a separate negative decision |
| Asymptotic mechanism | complete terminal repair plus $N^{-1}$ decay | proved | optimizer-excluding physical-frame principle |
| Framework obstruction | level-twelve grouped graph norm is $\Theta(N^{-1})$ | proved | whether optimizer is physically realizable |

## Historical q32 signed-permutation ledger

| Quantity | Current value | Label |
|---|---:|---|
| conditioned degree-eight total | $0.2482247496$ | proved |
| partial total through proved high sectors | $0.2797585469$ | proved |
| coarse completion total | $0.3331326055$ | diagnostic |
| coarse completion margin | $0.0002007278$ | diagnostic |
| unresolved balanced entries | 848 | proved count |
| unresolved complement/reversal orbits | 224 | generated count |
| largest unresolved physical coefficient | $0.0395939553$ | physical lower witness |
| reoptimized gate for that orbit | $0.0414623183$ | diagnostic scalar gate |
| best invariant correlated lower witness | $0.0395996496$ | physical diagnostic |
| next provisional-$1/q$ gate | $0.0379251204$ | diagnostic scalar gate |

## Candidate hard-instance scorecard

| Candidate | promise | active cost | finite-size evidence | adaptive stability | level-twelve behavior | decision |
|---|---|---|---|---|---|---|
| capped-Gaussian interpolation | constant margin after conditioning | six | constants unusable | proved asymptotically | positive sharp three-path family | asymptotic baseline |
| attenuated signed permutation | exact orbit plus attenuation/conditioning | six-compatible | q64: 888/888 certified; outward adaptive total $0.268858135059926$ | multiplier-one tree lift proved | finite-sector route | mathematical theorem at $N=4096$; present-day experimental credibility not certified |
| quadratic-bent exact plant | exact $F_{4,H}=\pm1$; no conditioning loss | six-compatible | 2,284 compatible higher entries; common gate $0.542524/N$; $1/N$ common proxy fails and $1/(2N)$ passes | open | larger orbit may cancel decorations; unproved | quantified first fallback; not promoted |

The quadratic-bent replacement screen gives an optimistic zero-higher-sector
floor of $0.281512032891$. A hypothetical common $1/(2N)$ contraction
would give $0.328670235127$, leaving $0.004663098206$, while common
$1/N$ gives $0.387618890853$. This is a candidate promotion gate, not a
coefficient theorem. Evidence:
notes/QUADRATIC_BENT_REPLACEMENT_SCREEN.md.

## Finite-size implementation window

| $N$ | $q=\sqrt N$ | current witness status |
|---:|---:|---|
| 256 | 16 | geometry supported; full ledger not calibrated |
| 512 | — | unsupported by current $N=q^2$ construction |
| 1024 | 32 | complete diagnostic routing ledger; not a theorem |
| 2048 | — | unsupported by current $N=q^2$ construction |
| 4096 | 64 | complete mathematical separation certified; hardware screen decided not yet experimentally credible |

### Historical pre-theorem $N=4096$ routing target

- Evaluating the inherited lower-sector formulas at $q=64$ gives a common
  unresolved-coefficient threshold of $0.199910665542$, or $0.199089176072$
  with the project reserve of $10^{-3}$.
- A common $1/\sqrt q=1/8$ envelope would give total $0.241894419850$.
- A more permissive two-tier target assigns $0.124035215254$ to the 724
  cubic-containing entries and $1/2$ to the other 164. It gives the
  diagnostic total $0.319181162161$ and margin $0.0141521711721$.
- These high-sector values were proof targets, not arbitrary-law theorems. The later dependency-exact ledger and adaptive theorem supersede them.
- The associated resource count is $M=16{,}384$ sign modes; experimental
  credibility has not been assessed.
- Evidence: notes/Q64_PAPER_TARGET_GATE.md and
  artifacts/q64_paper_target_gate.json.

### Proved $q=64$ block-coherent sector

- All 70 block-coherent balanced entries in the open profiles now have
  arbitrary-diagonal one-batch coefficients from the accepted weighted
  three-link path theorem.
- The calculation sums 196 compatible odd-record sectors using exact
  rational arithmetic.
- The coefficient range is $1/4096$ through
  $2609304163/39728800944=0.0656778986780$.
- Inserting these coefficients into the remaining two-tier target gives
  total $0.309405007008$ and diagnostic margin $0.0239283263253$.
- The other 818 internally split entries, interval certification, and the
  adaptive lift remain open.
- Evidence: notes/Q64_BLOCK_COHERENT_CONTRACTION.md and
  artifacts/q64_block_coherent_contraction.json.

### Proved $q=64$ chain-aware sectors

- The ten accepted Round 3 arbitrary-diagonal theorems close 40 further
  balanced entries at their actual order $q=64$.
- These entries are disjoint from the 70 block-coherent entries, so 110 of
  the 888 open entries now have theorem coefficients.
- The coefficient range is $0.00580989204377$ through
  $0.123974636390$.
- Their combined ledger insertion gives total $0.296090867182$ and margin
  $0.0372424661512$, an improvement of $0.0133141398259$ over the
  block-coherent insertion.
- The other 778 entries, interval certification, and the adaptive lift remain
  open. The largest remaining reusable class has 280 entries with exactly
  one split cubic and one split higher block.
- Evidence: notes/Q64_CHAIN_AWARE_INSERTION.md and
  artifacts/q64_chain_aware_insertion.json.

### Next $q=64$ shared-contraction gate

- The 778 remaining entries divide into structural classes of sizes 280,
  176, 140, 96, 48, 24, 8, and 6.
- The lead 280-entry class has exactly one internally split cubic and one
  internally split higher block; it contains 184 quintic and 96 septimic
  entries.
- Assigning that class one common coefficient, while leaving all other open
  entries at their routing targets, gives threshold coefficient
  $0.225536743566$.
- The coefficient $0.222921146951$ retains the declared $10^{-3}$ reserve.
  This is $1.79724077953$ times the original cubic routing target.
- These are floating proof-allocation gates, not coefficient theorems.
- Evidence: notes/Q64_REMAINING_CLASS_GATES.md and
  artifacts/q64_remaining_class_gates.json.

### Proved $q=64$ universal septimic sector

- Every fixed occurrence-split moment matrix is a cross Gram matrix of
  unit-modulus global character features, so arbitrary diagonal weights have
  universal coefficient at most one.
- The 96 entries with one split cubic and one split septimic block have an
  isolated reserve gate above one. Assigning them the proved universal
  coefficient therefore closes the entire septimic part of the former lead
  class without a signed-permutation-specific slice theorem.
- Together with the first 110 theorem entries, 206 of the 888 open entries
  now have arbitrary-law coefficients and 682 remain open.
- The resulting routing total is $0.329383221622$, with margin
  $0.00395011171117$ and $0.00295011171117$ left after the declared
  $10^{-3}$ allowance.
- The next specialized target is the 184-entry split-cubic/split-quintic
  class.
- Evidence: notes/Q64_UNIVERSAL_SEPTIMIC_INSERTION.md and
  artifacts/q64_universal_septimic_insertion.json.

### Proved $q=64$ universal multicubic tail

- The universal cross-Gram coefficient one also closes all 8 entries with
  three split cubics and all 6 entries with four split cubics.
- The first class does not move the current Perron optimum; adding the second
  gives total $0.331935829434$ and margin $0.00139750389895$.
- There are now 220 theorem entries and 668 open entries. Every other
  remaining structural class fails the reserve test at coefficient one.
- The live 184-entry split-cubic/split-quintic class has post-insertion common
  reserve gate $0.125261095651$.
- Evidence: notes/Q64_UNIVERSAL_MULTICUBIC_INSERTION.md and
  artifacts/q64_universal_multicubic_insertion.json.

### Corrected live quintic decision

- The 184 split-cubic/split-quintic entries divide into 104 extreme $1|4$
  splits and 80 balanced $2|3$ splits.
- Their live common reserve gate after all 220 theorem insertions is
  $0.125261095651$.
- The inherited local scales are $0.123974636390$ and $0.149556115743$,
  respectively. After assigning the first scale, the balanced-split reserve
  gate is $0.125681339751$.
- Assigning both local scales gives total $0.338248081665$, exceeding $1/3$
  by $0.00491474833148$.
- Decision: the local factors need an additional shared-chain contraction,
  or a sharper bound must recover the same margin from a coefficient-one
  class. These remain routing diagnostics, not new coefficient theorems.
- Evidence: notes/Q64_POST_UNIVERSAL_QUINTIC_GATE.md and
  artifacts/q64_post_universal_quintic_gate.json.

### Proved $q=64$ shifted middle-pair orbit

- The highest-impact unresolved quintic orbit has profile $(1,3,5,1)$ and
  split $(0,1,3,1)$, up to complement and reversal.
- A complete-row Schur-feature argument using exact cubic fixed-singleton
  and quintic fixed-triple slices gives arbitrary-law coefficient
  $0.0144099301059$ on all four entries.
- There are now 224 theorem entries and 664 open entries.
- The insertion gives total $0.325045063347$, margin
  $0.00828826998635$, and margin gain $0.00689076608741$.
- Keeping this theorem coefficient while assigning the other 180 quintic
  entries their local slice scales gives diagnostic total $0.329863829155$
  and margin $0.00346950417846$. The class-wide theorem remains open, but
  its prior numerical deficit is gone.
- Evidence: notes/Q64_SHIFTED_MIDDLE_PAIR_CONTRACTION.md and
  artifacts/q64_shifted_middle_pair_contraction.json.

### Proved $q=64$ reversed middle-pair orbits

- Two unresolved quintic orbits have profile $(1,1,5,3)$ and splits
  $(0,1,3,1)$ and $(0,1,2,2)$, up to complement and reversal.
- Using the complete column as a transposed row, Hadamard flatness cancels
  the summed singleton. The exact quintic fixed-pair slice and $N-2$ cubic
  completions give arbitrary-law coefficient $0.108667770790$.
- There are now 232 theorem entries and 656 open entries.
- The routing total is $0.323776780921$, with margin
  $0.00955655241272$ and gain $0.00126828242636$.
- Keeping all three new middle-pair theorem coefficients while assigning the
  other 172 quintic entries their local scales gives diagnostic total
  $0.326448847879$ and margin $0.00688448545451$.
- Evidence: notes/Q64_REVERSED_MIDDLE_PAIR_CONTRACTION.md and
  artifacts/q64_reversed_middle_pair_contraction.json.

### Proved $q=64$ fixed-singleton pair orbit

- For profile $(1,1,5,3)$ and split $(0,0,3,2)$, the complementary row
  fixes both singletons, a quintic pair, and one cubic cell.
- Hadamard flatness, the exact quintic pair slice, and
  $\binom{N-1}{2}$ cubic completions give arbitrary-law coefficient
  $0.0768303372012$.
- There are now 236 theorem entries and 652 open entries.
- The routing total is $0.323034695004$, with margin
  $0.0102986383292$ and gain $0.000742085916459$.
- Assigning the other 168 quintic entries their local scales gives
  diagnostic total $0.325310000395$ and margin $0.00802333293854$.
- Evidence: notes/Q64_FIXED_SINGLETON_PAIR_CONTRACTION.md and
  artifacts/q64_fixed_singleton_pair_contraction.json.

### Proved $q=64$ universal double-cubic class

- All 24 entries with exactly two split cubic blocks and no split
  higher-degree block have arbitrary-law cross-Gram coefficient one.
- There are now 260 theorem entries and 628 open entries.
- The routing total is $0.330563353867$, with margin
  $0.00276997946657$ and $0.00176997946657$ beyond the declared allowance.
- Assigning the other 168 quintic entries their local scales gives
  diagnostic total $0.333103976654$ and margin $0.000229356680$. A shared
  quintic improvement of at least $0.000770643321$ is needed to restore the
  allowance on that proxy.
- Evidence: notes/Q64_UNIVERSAL_DOUBLE_CUBIC_INSERTION.md and
  artifacts/q64_universal_double_cubic_insertion.json.

### Proved $q=64$ fixed-pair adjacent-row orbit

- For profile $(1,1,3,5)$ and split $(0,1,2,2)$, retain the complete
  $M_{13}M_{35}$ row as one Schur feature.
- A fixed cubic pair has at most $2(q-1)$ record-one L completions, $q-2$
  high-amplitude horizontal record-three completions, and
  $(q-1)(q-2)$ remaining record-three completions.
- Combining these counts with the proved fixed-pair $M_{35}$ tail bounds
  gives arbitrary-law coefficient $0.0699895941730$ on all four entries.
- There are now 264 theorem entries and 624 open entries. The routing total
  is $0.329248254134$, with margin $0.00408507919913$.
- Assigning local-slice scales to the other 164 quintic entries gives
  diagnostic total $0.331066966509$, margin $0.00226636682480$, and
  $0.00126636682480$ beyond the declared allowance.
- Evidence: notes/Q64_FIXED_PAIR_ADJACENT_ROW_CONTRACTION.md and
  artifacts/q64_fixed_pair_adjacent_row_contraction.json.

### Proved $q=64$ dual endpoint-slice Schur class

- Twelve entries in three orbits have favorable singleton placements for
  both internally split blocks.
- Exact cubic fixed-pair and quintic fixed-triple slice factorizations absorb
  the two distinctness masks. The remaining link is a unit cross-Gram Schur
  multiplier.
- Their arbitrary-law coefficient is $0.149556115743$.
- There are now 276 theorem entries and 612 open entries. The routing total
  is $0.329883605975$, with margin $0.00344972735879$.
- The remaining quintic inventory has 152 entries: 104 extreme and 48
  balanced.
- Evidence: notes/Q64_DUAL_ENDPOINT_SCHUR_INSERTION.md and
  artifacts/q64_dual_endpoint_schur_insertion.json.

### Live shared-quintic and adaptive acceptance gate

- The decorated adjacent complete-row theorem closes 16 degree-twelve
  entries in four orbits at coefficient $0.0200795672469$.
- Three degree-ten scalar completion-row theorems close 12 further entries
  at coefficients $0.00861554231015$, $0.0751888832423$, and
  $0.0754041939294$.
- Four whole-cubic decorated rows close another 16 degree-twelve entries,
  with coefficients from $0.000196155632204$ to $0.213518291069$.
- The final degree-ten orbit closes at coefficient $0.0910312181521$.
- The internal whole-cubic endpoint theorem closes 16 degree-twelve entries
  at coefficient $0.113036239514$; its balanced schema coefficient is
  $0.281119075921$.
- The balanced pair--triple mask theorem closes eight entries at coefficient
  $0.642693497508$ using the direct factor $1+\sqrt6+\sqrt3$.
- The adjacent double-cubic/quintic-mask theorem closes 32 entries at
  coefficients from $4.70180143564\times10^{-5}$ to $0.0134380552033$.
- The shared row/chain theorem closes the final 48 quintic entries in twelve orbits using five templates. Its coefficients range from $5.35051085909\times10^{-6}$ to $0.0203737451368$.
- There are now 428 theorem entries and 460 open entries. The routing diagnostic is $0.323811563171336$, with margin $0.009521770161998$.
- Retaining the allowance conditionally permits adaptive additive overhead at most $0.008521770161998$, equivalently multiplicative amplification at most $1.026317065637$.
- The residual classes have sizes 176, 140, 96, and 48, with common reserve gates $0.155710812601$, $0.535855735188$, $0.349193343122$, and $0.557126634930$.
- The 48-entry coefficients are theorems. The combined routing value, residual-class gates, and adaptive cap remain diagnostic requirements because 460 entries retain targets and no adaptive recurrence has been proved.
- Evidence: notes/Q64_DECORATED_ADJACENT_ROW_INSERTION.md,
  artifacts/q64_decorated_adjacent_row_insertion.json,
  notes/Q64_DEGREE_TEN_COMPLETION_ROW_INSERTION.md,
  artifacts/q64_degree_ten_completion_row_insertion.json,
  notes/Q64_WHOLE_CUBIC_DECORATED_ROW_INSERTION.md,
  artifacts/q64_whole_cubic_decorated_row_insertion.json,
  notes/Q64_LAST_DEGREE_TEN_CHAIN_INSERTION.md,
  artifacts/q64_last_degree_ten_chain_insertion.json,
  notes/Q64_INTERNAL_WHOLE_CUBIC_ENDPOINT_INSERTION.md,
  artifacts/q64_internal_whole_cubic_endpoint_insertion.json,
  notes/Q64_BALANCED_PAIR_TRIPLE_MASK_INSERTION.md,
  artifacts/q64_balanced_pair_triple_mask_insertion.json,
  notes/Q64_ADJACENT_DOUBLE_CUBIC_QUINTIC_ENDPOINT_INSERTION.md,
  artifacts/q64_adjacent_double_cubic_quintic_endpoint_insertion.json,
  notes/Q64_SHARED_QUINTIC_AND_ADAPTIVE_ACCEPTANCE.md, and
  artifacts/q64_shared_quintic_acceptance_gate.json,
  notes/Q64_SHARED_QUINTIC_ROW_CHAIN_INSERTION.md, and
  artifacts/q64_shared_quintic_row_chain_insertion.json.

### Historical conditional noncubic and recovered-universal insertion

- The degree-seven endpoint-product lemma proves joint injective-character factor $3/[(q-1)(q-3)]$.
- Sixteen middle degree-seven entries receive coefficients from $0.0170899096903$ to $0.0382305883153$.
- The other 124 noncubic entries receive the universal cross-Gram coefficient one, closing the entire 140-entry class.
- The intermediate routing total is $0.286902076794188$, which makes coefficient one affordable on the 96-entry two-split-cubic/one-split-higher class.
- There are now 664 theorem entries and 224 open entries. The live routing diagnostic is $0.328938230122941$, with margin $0.004395103210392$.
- The remaining class gates are $0.140343030565$ for 176 higher-split-only cubic-profile entries and $0.344887217413$ for 48 one-split-cubic/no-split-higher entries.
- Retaining the declared allowance conditionally permits adaptive additive overhead $0.003395103210392$, or multiplier $1.01032140049250$.
- Evidence: notes/Q64_NONCUBIC_RECOVERED_UNIVERSAL_INSERTION.md and artifacts/q64_noncubic_recovered_universal_insertion.json.

### Historical conditional whole-higher split-cubic insertion

- The 48 one-split-cubic/no-split-higher entries reduce to 24 favorable cubic endpoints, 16 internal singleton--whole-cubic endpoints, and eight complete singleton--quintic--whole-cubic wedges.
- Their arbitrary-law coefficients are $0.124035215254$, $0.015625$, and $0.00846466875312$.
- The historical conditional calculation reported 712 entries, 176 open entries, and routing diagnostic $0.328477421166173$. These cumulative values are withdrawn.
- The final higher-split-only cubic-profile class has common reserve gate $0.142581909211$.
- Retaining the declared allowance conditionally permits adaptive additive overhead $0.003855912167160$, or multiplier $1.01173874342252$.
- Evidence: notes/Q64_WHOLE_HIGHER_SPLIT_CUBIC_INSERTION.md and artifacts/q64_whole_higher_split_cubic_insertion.json.

### Preserved same-side whole-link theorem

- A mask-aware class theorem covers 96 entries using one same-side whole link
  and the exact inclusion--exclusion factor for the sole split higher block.
- Its coefficients range from (0.046875) to (0.167292848473).
- The theorem remains reusable. Its historical cumulative 808-entry routing total and adaptive allowance remain withdrawn; the dependency and final-residual registries have since been rebuilt independently.
- Evidence: `notes/Q64_SAME_SIDE_WHOLE_LINK_INSERTION.md` and
  `artifacts/q64_same_side_whole_link_insertion.json`.

## Paper resource row

| Side | Size | hard dose | error/advantage | resource status | theorem status |
|---|---:|---:|---:|---|---|
| active | $4096$ target | 6 | error $81/256$; margin $13/768$ | three folded flags; ideal resource row complete | proved |
| passive | $4096$ | exclude 6 | transcript reserve $0.063475198273408$ | mathematical scope complete; experimental feasibility is an active-side implementation issue | proved plus certified numeric |

The active row uses three independent single photons. Each photon is coherent
over a binary path flag and 4096 mode labels and crosses two of the four sign
blocks. The controlled words are $D_1HD_2$ and $D_4HD_3H$; a path-$X$
receiver yields one binary flag, and three flags are majority decoded. No
postselection is used. Evidence: notes/ACTIVE_SIX_DOSE_RESOURCE_ROW.md and
artifacts/active_six_resource_row.json.

The exact scalar active robustness gate is flag expectation
$\mu_*=0.226073713789$. At the promise boundary this requires combined
multiplicative contrast above $0.904294855157$, or per-pass power
transmission above $0.950944191400$ if the only imperfection is two equal
lossy sample traversals. This is a resource allocation, not a device-specific
noise theorem. Evidence: notes/ACTIVE_SIX_ROBUSTNESS_GATE.md and
artifacts/active_six_robustness_gate.json.

The dated primary-evidence screen gives verdict `NOT YET EXPERIMENTALLY CREDIBLE`. At nominal detector efficiency 98%, all non-detector imperfections share only $0.3491602$ dB of geometric-mean loss budget. Existing reviewed time-bin and spatial demonstrations do not supply a coherent $H_{4096}$ at this threshold. Evidence: notes/EXPERIMENTAL_FEASIBILITY_DECISION.md and artifacts/q64_experimental_feasibility_gate.json.

Active dose five, passive upper bounds, and the full protocol frontier are
preserved in the long-horizon backlog but are not Round 4 completion gates.

## Result-entry template

For each new result record:

- date and file;
- exact model and hard instance;
- dimension and dose;
- analytic or certified numerical statement;
- theorem/protocol/diagnostic label;
- regression command;
- effect on the boundary; and
- next decision.
