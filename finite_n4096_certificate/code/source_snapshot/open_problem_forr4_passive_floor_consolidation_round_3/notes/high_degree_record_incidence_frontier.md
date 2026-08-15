# High-degree record/incidence frontier

Date: 2026-07-15

Status: rigorous coarse coefficient bounds, exact finite inventory, and a
route-selection verdict.  This is not the complete one-batch certificate.

## 1. Question

After the chained degree-eight repair, the current finite-size target assigns
coefficient \(1/q=1/32\) to the otherwise open degree-ten/twelve splits.  The
purpose of this audit is to determine how much of that target follows from
generic maximum-entry, rank--Frobenius, and incidence arguments before doing
more specialized chain algebra.

There are 43 open profiles, 6,016 profile/split entries, 92 compatible record
triples across the profiles, and 14,624 record-sector/split entries.

## 2. Universal record-sector entry bound

Let \(r_j\) be the number of odd labels on link \(j\).  Matching one prescribed
odd-label set through a uniform hidden permutation has probability

$$
{1\over {q\choose r_j}}.
$$

The conditional phase moment of every even decoration has magnitude at most
one.  Therefore every entry in record sector
\(r=(r_0,r_1,r_2)\) obeys

$$
|M_{S,T}|
\le
\prod_{j=0}^2 {1\over {q\choose r_j}}.
\tag{2.1}
$$

Every compatible record is a positive odd integer no larger than the degrees
of its two incident blocks.

## 3. Exact one-axis incidence

Fix an \(s\)-cell selected support in one block and write its multiplicities
on one link axis as \(f_1,\ldots,f_q\).  A degree-\(a\) extension chooses
\(e_i\) new cells from the \(q-f_i\) unused cells with that label.  Hence the
number of extensions with final odd-label record \(r\) is exactly

$$
\sum_{
\substack{
0\le e_i\le q-f_i\\
\sum_i e_i=a-s\\
\sum_i(f_i+e_i\bmod2)=r
}}
\prod_{i=1}^q {q-f_i\choose e_i}.
\tag{3.1}
$$

A two-coordinate dynamic program evaluates (3.1).  Maximizing it over the
integer partitions of \(s\) gives the exact uniform one-axis incidence
\(I_q(a,s,r)\).  For a middle block constrained on both axes, the number of
compatible extensions is safely at most the smaller of the two corresponding
one-axis maxima.  The checker exhaustively reproduces these maxima at \(q=4\)
for degrees one, three, and five.

For one record sector, combine (2.1) with the smaller of:

1. the crude square-root rank \(q^k\), where
   \(k=\min(\sum s_b,L-\sum s_b)\);
2. the square root of the product of row-extension incidences; and
3. the corresponding column-extension quantity.

Record sectors may be triangle-summed.  Independently, the entire profile can
be bounded in one step using maximum entry \(q^{-3}\) and the all-support
incidence

$$
{q^2-s\choose a-s}.
$$

The implemented coefficient is the minimum of the whole-profile bound and
the record-sector triangle bound.

## 4. Exact frontier at \(q=32\)

The generic estimate proves coefficient at most \(1/q\) for:

| class | certified | total |
|---|---:|---:|
| degree ten | 256 | 896 |
| degree twelve | 882 | 5,120 |
| all open profile/splits | 1,138 | 6,016 |

This looks like useful progress until occupation compatibility is imposed.

For a profile of total degree \(L\), let \(k=\sum_b s_b\).  The exact pairing
condition in the occupation ledger maps a row state of total occupation six
to a column state of total occupation

$$
6+L-2k.
$$

Both sides can occur in the dose-six ledger only when

$$
k={L\over2}.
\tag{4.1}
$$

There are 144 such degree-ten entries and 744 such degree-twelve entries,
for 888 dose-six-relevant entries in 234 reversal/complement orbits.  The
generic estimate proves \(1/q\) for

$$
\boxed{0\text{ of the }888\text{ dose-six-relevant entries}.}
$$

All 1,138 generic successes are off shell and contribute exactly zero to the
current six-dose Perron ledger.  Replacing their provisional \(1/q\) charges
by the sharper proved values consequently leaves the optimized diagnostic
unchanged at \(0.322669154028\).

## 5. Which balanced orbits matter first

At the optimized all-open diagnostic, the Perron vector gives an exact local
sensitivity to each provisional coefficient.  Grouping by complement and
path reversal gives the following leading route priorities:

| rank | representative | current charge | Perron contribution |
|---:|---|---:|---:|
| 1 | \((3,1,1,5):(1,1,0,3)\) | \(0.0395939553\) | \(0.0039451833\) |
| 2 | \((3,1,1,5):(0,1,0,4)\) | \(1/32\) | \(0.0017081466\) |
| 3 | \((1,1,3,5):(0,0,1,4)\) | \(1/32\) | \(0.0015451842\) |
| 4 | \((3,1,1,5):(1,1,1,2)\) | \(1/32\) | \(0.0015394748\) |
| 5 | \((3,1,1,5):(1,0,0,4)\) | \(1/32\) | \(0.0014627498\) |

Rank one is the already known opposite-endpoint physical orbit.  The
ranking is a derivative/decomposition of the current numerical target, not
an upper-bound theorem.  It prevents the next proof attempt from being
selected by raw coefficient size or by the number of record sectors.

Subsequent result: rank two has now been controlled rigorously, though not at
the provisional \(1/q\).  Its complement/reversal orbit has arbitrary-law
coefficient at most \(0.0934752746\).  Paying this conservative value raises
the diagnostic to \(0.3262818609\) but leaves \(0.0070514724\) below
threshold.  Rank three is now controlled at \(0.0162724693<1/q\) by an
exact cubic fixed-pair slice.  Inserting both theorems gives diagnostic
\(0.3255638580\), with \(0.0077694754\) slack.

The former rank-five orbit \((3,1,1,5):(1,0,0,4)\) is now controlled by
exact cubic fixed-pair and quintic fixed-four slices at coefficient
\(0.1737428008\).  This is below the generic two-mask value but above
\(1/q\).  With all three theorems inserted the diagnostic becomes
\(0.3326651190\), with only \(0.0006682144\) slack.

After removing the first three proved and explicitly forced orbits and
reoptimizing, the largest remaining provisional orbit at that stage was
\((3,1,1,5):(1,1,1,2)\), with Perron contribution \(0.00150951797\).
Its reoptimized acceptance gate is coefficient
\(0.0450405468=1.44129750/q\).
It is followed by \((3,1,1,5):(0,0,1,4)\), with contribution
\(0.00128274051\).  See
`leading_balanced_disjointness_contraction.md`,
`adjacent_balanced_cubic_slice_contraction.md`, and
`separated_balanced_endpoint_slice_contraction.md`.

Subsequent result: the \((3,1,1,5):(1,1,1,2)\) orbit passes that gate at
coefficient \(0.0250967461\).  The proof preserves the common law across its
endpoints; separate worst-case slice bounds would instead give
\(0.0991470258\) and fail.  With four chain-aware theorems inserted, the
diagnostic is \(0.3323683002\), leaving \(0.0009650332\), and sixteen of 888
balanced entries are controlled.  Reoptimization makes
\((3,1,1,5):(0,0,1,4)\) the largest remaining provisional orbit, with
Perron contribution \(0.00128192174\) and acceptance gate
\(0.0542506298=1.73602015/q\).  See
`internal_singleton_shared_law_contraction.md`.

Subsequent fifth result: the \((3,1,1,5):(0,0,1,4)\) orbit has coefficient
\(0.0311889051<1/32\).  Its exact quintic row energy avoids the generic
distinctness factor, and the remaining cubic--Hadamard chain collapses by
rank--Frobenius.  With five chain-aware theorems inserted, the diagnostic is
\(0.3323657941\), leaving \(0.0009675392\), and twenty of 888 entries are
controlled.  Reoptimization makes \((1,1,3,5):(0,1,1,3)\) the largest
remaining provisional orbit, with Perron contribution \(0.00116190969\) and
acceptance gate \(0.0570749885=1.82639963/q\).  See
`column_cubic_quintic_row_contraction.md`.

Subsequent sixth result: the \((1,1,3,5):(0,1,1,3)\) orbit has coefficient
\(0.0422410016\).  Factoring the complete adjacent row avoids separate
endpoint masks; exact L-shape incidence controls record one and a horizontal
Walsh identity controls the leading record-three triple pattern.  With six
chain-aware theorems inserted, the diagnostic is \(0.3327757792\), leaving
\(0.0005575541\), and twenty-four of 888 entries are controlled.
Reoptimization makes \((3,1,1,5):(0,1,1,3)\) the largest provisional orbit,
with Perron contribution \(0.00100671865\) and gate \(0.0484819899\).  See
`adjacent_balanced_row_slice_contraction.md`.

Subsequent seventh result: the \((3,1,1,5):(0,1,1,3)\) orbit has coefficient
\(0.0370952793\).  Extracting the exact fixed-three quintic row energy leaves
a whole cubic endpoint whose proportional columns compress by support XOR to
a weighted Walsh matrix.  With seven chain-aware theorems inserted, the
diagnostic is \(0.3329643636\), leaving \(0.0003689697\), and twenty-eight of
888 entries are controlled.  Reoptimization makes
\((1,3,5,1):(0,2,2,1)\) the largest provisional orbit, with Perron
contribution \(0.000999555962\) and gate \(0.0426269309\).  See
whole_cubic_quintic_triple_contraction.md.

Subsequent eighth result: the \((1,3,5,1):(0,2,2,1)\) orbit has coefficient
\(0.0285281523<1/32\).  Exact cubic and quintic fixed-pair squared sums are
retained before taking the universal middle-link maximum.  With eight
chain-aware theorems inserted, the diagnostic is \(0.3328775891\), leaving
\(0.0004557442\), and thirty-two of 888 entries are controlled.
Reoptimization makes \((1,1,3,5):(0,0,3,2)\) the largest provisional orbit,
with Perron contribution \(0.000981156614\) and gate \(0.0454321892\).  See
middle_cubic_quintic_pair_contraction.md.

Subsequent ninth result: the \((1,1,3,5):(0,0,3,2)\) orbit has coefficient
\(0.0250919472<1/32\).  A normalized complete \(M_{35}\) row leaves the
XOR-labelled \(H_NM_{13}\) chain available for a separate \(1/q\)
compression.  With nine chain-aware theorems inserted, the diagnostic is
\(0.3326862124\), leaving \(0.0006471209\), and thirty-six of 888 entries are
controlled.  Reoptimization makes \((1,3,5,1):(0,2,3,0)\) the largest
provisional orbit, with Perron contribution \(0.000924037577\) and gate
\(0.0529177167\).  See whole_cubic_middle_pair_contraction.md.

Subsequent tenth result: the \((1,3,5,1):(0,2,3,0)\) orbit has coefficient
\(0.0462425962<0.0529177167\).  A scalar completion-row Schur feature leaves
a repeated, column-twisted \(H_N\otimes H_N\) residual.  With ten
chain-aware theorems inserted, the diagnostic is \(0.3331326055\), leaving
\(0.0002007278\), and forty of 888 entries are controlled.  Reoptimization
makes \((1,1,5,3):(0,1,3,1)\) the largest provisional orbit, with Perron
contribution \(0.000936859407\) and gate \(0.0379251204\).  See
double_endpoint_cubic_quintic_row_contraction.md.

## 6. Decision

The generic record/incidence route is accepted as a theorem but stopped as
the lead strategy for the realistic-size result.  Its successes lie entirely
outside the dose-six shell.  The next Track A theorem must exploit balanced
chain structure that this argument discards: squared-slice energy, a
whole-chain Bessel/Carleson inequality, or an exact physical-law reduction
for a leading orbit.

This was also a portfolio review trigger.  The bounded Track B exponent audit
was completed before the next balanced calculation, and the first ten
chain-aware theorems have preserved the finite-size route, still with little
margin.  Continue for one bounded orbit in reoptimized ranked order only if
its chain architecture can plausibly pass the \(0.0379251204\) gate.  A
certified failure of a chain-aware gate
triggers the named alternative-hard-instance comparison rather than another
sweep of generic incidence refinements.

## 7. Reproduction

- `searches/high_degree_record_incidence_frontier.py` implements the exact
  inventory, incidence dynamic program, rigorous coefficient bounds,
  dose-six shell filter, symmetry orbits, and Perron priority diagnostic.
- `tests/high_degree_record_incidence_frontier.py` compares the incidence
  dynamic program with direct \(q=4\) enumeration and protects every frontier
  count and the leading route priority.
