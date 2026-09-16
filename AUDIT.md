# Standalone exposition and current finite scope

The canonical [self-contained PDF](fourfold_forrelation.pdf) includes every mathematical ingredient used for its stated results. It also proves a finite dose-two obstruction at `N=4096` and an all-split, 210-state sufficient spectral criterion for arbitrary hard-dose-six parallel probes, including photon-number coherence. The needed dose-six spectral estimate remains unresolved.

The detailed audit below records the earlier finite coefficient ledger and its balanced-parity restriction. It remains useful for diagnosing those attempted estimates; the new PDF does not depend on that ledger or on this audit. The all-split reduction defines exact kernel norms directly and retains the edges excluded by the balanced restriction. No earlier unproved numerical separation is promoted to a theorem.

---

# Proof audit: current mathematical status

Original repository base: `95c4c006d7566c7a01f68f0ee2ed465714dadf59`.
Independent-review checkpoint: `2b7a3cb72e4f54b4989d66e213d4436333ad4554`.

The asymptotic parallel theorem and the two-pass construction survive review. The finite lower bound is **conditional, not established**: its previous photon-number scope was incorrect, and some coefficient derivations contain false intermediate bounds. The arithmetic certificate is a certificate for the stored matrix, conditional on valid analytic inputs. Both adaptive extensions remain unproved.

| Argument | Current assessment |
|---|---|
| Parallel Fourier bound | Valid, and strengthened by a factor of two using `E − I/2` |
| Gaussian moments, variance and promise transfer | Valid; direct promise-mass accounting gives a much shorter closing argument |
| Two-pass protocol, exact plant and finite concentration | Valid |
| Finite physical-number-to-balanced-cut step | False; parity-support size is the relevant quantity |
| Complete-kernel-to-count-matrix reduction | Proved in the revised note under an explicit coefficient hypothesis |
| Completed local-Walsh contraction | False in the claimed generality, including at q=64 |
| Four-cubic maximum-entry premise | False for q≥8 |
| Several inherited/quintic/residual coefficient steps | Unproved as supplied; see the family guide |
| Directed matrix arithmetic and positive Collatz candidate | Valid conditional calculation; cannot repair false or missing analytic premises |
| Adaptive norm arguments | Still incomplete; the finite abstract norm inference has a counterexample |

## Photon number is not parity-support size

This error was introduced in the first shortened presentation: it called the ledger's support counts physical Fock totals. The earlier source [occupation note](finite_n4096_certificate/code/source_snapshot/open_problem_forr4_passive_floor_consolidation_round_3/notes/occupation_compatibility_and_promise_concentration.md) explicitly defines them as parity-support cardinalities.

For an occupation ν, put `S_b(ν)={i:ν_bi is odd}` and `d_b(ν)=sum_i floor(ν_bi/2)`. Its physical totals are `N_b(ν)=|S_b(ν)|+2d_b(ν)`. Thus for a pair of occupations the physical identity is

`N(μ)=N(ν)+a−2t+2(d(μ)−d(ν))`,

where `a_b=|S_b(ν) △ S_b(μ)|` and `t_b=|S_b(ν) \ S_b(μ)|`. Equal physical totals do not imply the balanced condition `2|t|=|a|`.

For example, superpose three photons in coordinate zero of block one with one photon in coordinate zero of each of blocks two, three and four. Both components have physical photon number three; their parity sizes are one and three. The surviving quartic profile is `(1,1,1,1)` with split `(1,0,0,0)`. Its unconditioned half-difference moment is

`beta^4/q^3 = 130321/102400000000 > 0` at q=64.

The balanced ledger drops this entry. This refutes the deletion step, not the final proposed TV bound. More generally, adding photon pairs realizes every equal-parity pair of supports of size at most six in one physical-number sector. Physical number diagonality therefore cannot justify deleting these unbalanced parity edges.

The corrected conditional theorem requires `[rho,Q]=0`, where `Q` counts odd-occupied modes, together with hard signal dose at most six. Grouping all Fock amplitudes with the same parity support includes repeated modes and idlers exactly. If signal occupations are collision-free, Q equals physical photon number; in general it does not.

## An exact completed-Walsh counterexample

The [local-Walsh note](finite_n4096_certificate/code/source_snapshot/open_problem_forr4_passive_floor_consolidation_round_4/notes/Q64_MASKED_LOCAL_WALSH_REPAIR.md) first uses an incorrect phase identity. The product of the two Walsh phases through a singleton z is proportional to `H_N(z, xor(S) xor xor(T))`, not `H_N(xor(S),xor(T))`; it depends on z.

There is also a substantive norm error. Let `M` denote the completed signed-permutation link moment, let x range over individual cells, and let E range over unordered two-cell sets. Put

`V[x,E] = q M({x} △ E, {0})`.

Repeated hidden signs cancel on overlaps, so the symmetric difference is essential. Exact permutation averaging gives:

| Pair E | Position of x | V[x,E] |
|---|---|---:|
| Same column b | Column b | 1 |
| Same column b | Another column | −1/(q−1) |
| Distinct columns b,c | x belongs to E | 1 |
| Distinct columns b,c | Another cell in b or c | −1/(q−1) |
| Distinct columns b,c | Another column | 0 |

The row space splits into constants, column-constant vectors of total sum zero, and vectors summing to zero in every column. On these spaces `VVᵀ` has eigenvalues

`0`, `q^4/[2(q−1)]`, `q^3/(q−1)`,

with multiplicities `1`, `q−1`, `q(q−1)` respectively. To see this directly, each same-column pair contributes the column vector

`u_b = q/(q−1) (1_b − 1/q · 1)`,

repeated `binom(q,2)` times. Each cross-column pair contributes `a_rb+a_sc`, with

`a_rb = q/(q−1) (e_rb − 1/q · 1_b)`.

Summing over r,s cancels the cross terms. These two contributions act on the two indicated nonconstant subspaces and give the displayed eigenvalues.

Consequently, under uniform row and column laws,

`||V / sqrt(N binom(N,2))||_1 = (1+sqrt(2q))/sqrt(q+1) > 1`.

At q=8 this is exactly `5/3`; at q=64 it is `(1+8sqrt(2))/sqrt(65)`. Thus the completed amplitude is not a unit Schur factor.

For the completed two-link wedge

`W[(x,z),(E,t)] = M({x} △ E,z) H_N(z,t)`,

a row phase and a column permutation identify `qW` with `V ⊗ H_N`. Its uniform weighted trace coefficient is therefore

`(1+sqrt(2q))/(q sqrt(q+1)) > 1/q`.

This disproves the completed contraction used by several accepted proofs. It does **not** by itself disprove their final masked four-block coefficients: a different joint argument could still establish those. The [coefficient guide](finite_n4096_certificate/code/COEFFICIENTS.md) identifies affected entries, including four retained inherited rows, eight quintic rows, 128 local-Walsh rows and 48 residual rows.

## An exact four-cubic counterexample

Let `V={(0,0),(1,0),(2,0)}` and `H={(0,0),(0,1),(0,2)}`. Sign averaging gives

`M(V,H)=1/q`, `M(H,V)=1/binom(q,3)`.

The first average survives precisely when the permutation fixes zero; the second survives precisely when it maps `{0,1,2}` onto itself as a set. Every surviving Walsh phase is +1.

The completely physical four-block tuple `(V,H,V,H)` is in record sector `(1,3,1)`. Its entry is `g/q^2`, where `g=1/binom(q,3)`. The [four-cubic implementation](finite_n4096_certificate/code/source_snapshot/open_problem_forr4_passive_floor_consolidation_round_4/searches/q64_masked_four_cubic_incidence_repair.py) instead uses `b^2 g`, with `b=(q+2)/[q(q−1)(q−2)]`, as a universal maximum. For q≥8 this is too small. At q=8,

`actual = 1/3584 > 25/1580544 = asserted maximum`.

The old exhaustive q=4 check misses the issue because b=1/q there. Its selected q=8 examples show attainment of particular values, not a universal maximum.

This refutes the maximum-entry premise supporting 38 coefficients, not those final norm coefficients themselves. Simply replacing b by the valid 1/q bound does not recover the recorded certificate: the corresponding incidence bounds exceed one. A stronger joint shape analysis is needed.

## Other coefficient obligations

The residual orbit generated by `(1,5,3,3):(0,3,3,0)` uses an unspecified normalized Walsh dressing. The supplied fixed-endpoint slice calculation bounds the relevant row energy by `F_3`; it does not supply the claimed `F_3/q^2`. Four coefficients remain unproved for that reason. Their tests check separate factors, not the missing joint normalization.

Other reviewed physical row-energy and completed-character reductions are sound in their stated scope. The [family-by-family record](finite_n4096_certificate/code/COEFFICIENTS.md) distinguishes these from inherited scalar enumerations that were not exhaustively re-proved. Correct reductions do not justify the unrelated failed ones.

## The finite adaptive error

Let `K=J_n`, the all-ones matrix, and let every scalar feature be `U_i=V_i=1/sqrt(n)`. Both feature families have squared mass one. The kernel has Schur factorization norm one. With `p_i=w_i=1/n`, the weighted matrix is `J_n/n`, whose trace norm is one. But the entry sum is

`sum_ij K_ij <V_j,U_i> = n`.

Thus normalized feature masses and a contractive Schur symbol do **not** imply the needed scalar bound. Pairing a matrix with the all-ones observable can cost its operator norm, which is `n`. The old proof controls the matrix in its equation (3.4), then passes to the entry sum without controlling that pairing. The old regression checks the induction and individual Schur/Perron facts, but never tests this missing implication.

This is a counterexample to the abstract inference, not a counterexample to the desired physical adaptive lower bound. A repair must use additional causal/measurement structure that the invariant has discarded.

## The asymptotic adaptive gap

The original fresh-cut lemma introduces input/cut spaces and asserts a contraction after a zig-zag rearrangement. It does not give compatible explicit maps that establish the rearranged operator norm. In particular, an instrument isometry on its original input-output split is not automatically contractive across a new tensor split. Outcome-dependent preparations must be accounted for in that same norm.

The new proof avoids this issue for a parallel probe. Its input and effect are fixed, and a record kernel is explicitly the product of a pulled-back small coefficient matrix, a residual-parity equality kernel and two one-sided filters. This proves the same Fourier bound for that class. No counterexample to the broader adaptive Fourier inequality is asserted; it is recorded as an open obligation rather than an established lemma.

## What the arithmetic establishes

The unchanged backend reproduces the 888 accepted balanced high-sector rows and the inherited lower-degree inputs, excludes unbalanced rows, reconstructs the 210-state matrix with directed arithmetic, evaluates the positive Collatz candidate and checks the rational promise relaxation. Its `proved`/`certified` field names are historical metadata.

Upward rounding certifies domination of the supplied floating values. Exact squared checks additionally certify domination of specified rational formulas for several families. Neither establishes that those formulas bound the complete physical kernels. The new regression independently reproduces the parity, Walsh and four-cubic counterexamples instead of importing the frozen coefficient code.

## Scope after independent review

- Proved: the parallel asymptotic theorem, with arbitrary idlers, repeated modes and number coherence; the two-pass dose-six construction; the exact finite plant and concentration bound.
- Conditional: the finite TV bound under valid complete-kernel inequalities and diagonality in parity-support size. The revised note supplies the complete residual-support lifting proof.
- Open: the intended unconditional finite separation, the physical-number-only claim with repeated modes, and both adaptive extensions. Mean-dose bounds and multipass optimality are not established.

The frozen 469-file snapshot remains unchanged. The external sensing manuscript was not modified. See [REVIEW.md](REVIEW.md) for adopted reviewer inputs and [VALIDATION.md](VALIDATION.md) for executable checks.
