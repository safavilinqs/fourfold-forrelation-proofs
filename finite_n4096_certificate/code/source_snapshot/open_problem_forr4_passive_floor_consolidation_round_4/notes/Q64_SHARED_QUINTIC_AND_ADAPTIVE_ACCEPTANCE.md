# q64 shared-quintic and adaptive acceptance gate

Date: 2026-07-16

Status: historical acceptance checkpoint for the valid 48-entry
shared-quintic theorem. All cumulative counts, routing values, and adaptive
allowances in this note are withdrawn because later totals depend on the
unmasked universal lemma. The current gate is
`Q64_MASKED_UNIVERSAL_AUDIT.md`.

## Decision

The proposed shared $q=64$ contraction was the correct bounded project. Since the 168-entry checkpoint, the fixed-pair, dual-endpoint, decorated-row, degree-ten completion-row, whole-cubic decoration, final degree-ten chain, internal whole-cubic endpoint, balanced pair--triple mask, adjacent double-cubic/quintic-mask, and final shared row/chain theorems have closed the entire class. The last theorem covers all 48 remaining degree-twelve entries: 32 extreme $1|4$ quintic cuts and 16 balanced $2|3$ cuts.

The finite-orientation interpretation was essential: five row/chain templates cover twelve complement/reversal orbits. The largest coefficient is $0.0203737451368$, and the optimized routing diagnostic is $0.323811563171336$.

## Class-level routing acceptance

The contraction is evaluated against the following criteria:

1. It covers all 48 entries for arbitrary correlated diagonal physical
   laws, including complement and reversal, without an invariant-law or
   product-law assumption.
2. Every coefficient is connected to a proved row, chain, Gram, or shared
   operator inequality and a deterministic regression.
3. Reoptimizing the q64 routing ledger with the other 460 open entries held
   at their frozen target values and rounding outward gives

   $$
   U_{\rm route}\le\frac13-10^{-3}=0.332333333333\ldots.
   $$

4. It regenerates gates and a closure plan for the 460 entries left in the
   four residual classes:

   - 176 higher-split-only cubic-profile entries;
   - 140 noncubic-profile entries;
   - 96 two-split-cubic/one-split-higher entries; and
   - 48 one-split-cubic/no-split-higher entries.

Criteria 1, 2, and 4 pass. Criterion 3 passes in floating arithmetic with visible slack, but outward-rounded interval certification is still open. Thus the arbitrary-law coefficient theorem is accepted while the complete numerical one-batch gate is not yet certified.

The displayed $1/3-10^{-3}$ condition is passed with diagnostic slack $0.008521770161998$ after retaining the allowance. It remains a class-level routing gate while those 460 entries retain targets and becomes a one-batch theorem gate only after every residual coefficient is proved.

Two sufficient coefficient targets are now known:

- one common coefficient at most $0.410314553367$ on all 48 entries; or
- coefficient at most $0.123974636390$ on the 32 extreme entries and at
  most $0.149556115743$ on the 16 balanced entries.

After fixing the extreme value, the balanced reserve gate is actually
$1.09336558289$, so the two-tier target has substantial room.

The accepted coefficient map improves the earlier two-tier proxy and gives

$$
U_{\rm route}^{\rm accepted}=0.323811563171336,
$$

with raw diagnostic margin $0.009521770161998$.

The recomputed common reserve gates for the four residual classes are:

- $0.155710812601$ for the 176 higher-split-only cubic-profile entries;
- $0.535855735188$ for the 140 noncubic-profile entries;
- $0.349193343122$ for the 96 two-split-cubic/one-split-higher entries; and
- $0.557126634930$ for the 48 one-split-cubic/no-split-higher entries.

Against the frozen targets $0.124035215254$, $0.5$, $0.124035215254$, and $0.124035215254$, respectively, the noncubic class has the least relative headroom and is the next falsification target.

## Simultaneous adaptive requirement

Before investing in a long class proof, write the candidate adaptive theorem
in one of the following quantitative forms:

$$
\operatorname{TV}_{\rm adaptive}(D\le6)
\le U_1+\Delta_{\rm ad}(N,6),
$$

or

$$
\operatorname{TV}_{\rm adaptive}(D\le6)
\le C_{\rm ad}(N,6)U_1.
$$

The recurrence must be uniform in outcome width, depth, dose partition,
posterior-selected child laws, fresh batches, and branchwise hard-dose-six
trees. A calculation for a fixed batching pattern does not qualify.

At the accepted 48-entry insertion, retaining the declared $10^{-3}$ allowance requires

$$
\boxed{\Delta_{\rm ad}(4096,6)\le0.008521770161998}
$$

or, in multiplicative form,

$$
\boxed{C_{\rm ad}(4096,6)\le1.026317065637}.
$$

These are design requirements inferred from a routing calculation with 460 frozen targets, not certified adaptive limits. After all 460 entries are proved, recompute them from the outward-rounded complete one-batch value $U_1$. Without retaining the allowance, the present additive limit is $0.009521770161998$. The likely successful theorem is posterior-stable with very little amplification, or the remaining one-batch contractions must create more margin.

## Kill and pivot conditions

Stop this project and reconsider the witness or contraction architecture if:

- a legal physical law exceeds any proposed class envelope;
- the proved coefficient map cannot reach the outward-rounded
  $1/3-10^{-3}$ routing gate;
- closing the quintic class leaves no credible contraction program for the
  four residual classes;
- every plausible outcome-uniform adaptive recurrence exceeds the dynamic
  cap $1/3-10^{-3}-U_1$; or
- the only route is a placement-only temporal square function already ruled
  out by the Round 2 barrier.

An adaptive recurrence that fails the current cap does not automatically
kill the hard instance. It requires either a stronger one-batch contraction
that creates more room, a posterior-stable norm with smaller amplification,
or the predeclared hard-instance pivot.

## Deliverables

The bounded project deliverables now stand as follows:

1. complete: a machine-readable 48-entry coverage map;
2. complete: one shared proof schema with five finite row/chain templates;
3. floating complete, interval certification open: class insertion, residual-class gates, and routing margin;
4. open: a precise adaptive recurrence candidate and evaluated $N=4096,D=6$ overhead;
5. open: a pass/fail decision against the dynamic adaptive cap; and
6. open: an updated retain-or-pivot decision after the 460-entry and adaptive tests.

Reproduce the numerical gates with:

    /opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_shared_quintic_acceptance_gate.py --write-artifact
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_shared_quintic_acceptance_gate.py
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_shared_quintic_row_chain_insertion.py --write-artifact
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_shared_quintic_row_chain_insertion.py
