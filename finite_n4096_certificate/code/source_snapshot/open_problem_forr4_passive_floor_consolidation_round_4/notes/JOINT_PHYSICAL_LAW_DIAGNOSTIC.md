# Joint physical-law diagnostic

Date: 2026-07-15

## Question

Can the leading unresolved families be evaluated under one actual passive
probe law, rather than as unrelated scalar orbit constants?

The first exact calculation now exists. It is a baseline for the joint-law
route, not yet the retain-or-pivot certificate.

## The compatible object

Let (n=(n_0,n_1,n_2,n_3)) be one of the dose-six occupation states and let
(ho_n) be its weight. For every (n), choose a conditional distribution
(mu_n) on four-block coordinate supports

$$
A=(A_0,A_1,A_2,A_3),
\qquad |A_j|=n_j.
$$

The physical probe law is

$$
p(n,A)=ho_nmu_n(A).
$$

For two configurations (A,B), put (Q_j=A_j\mathbin\triangle B_j). The
plus-plant Fourier moment factors over the three independent hidden
signed-permutation links:

$$
\mathcal M(Q_0,Q_1,Q_2,Q_3)
=
\prod_{j=0}^{2}
\mathbb E_P\left[
\chi_{Q_j}(KP)\chi_{Q_{j+1}}(PK)
\right].
$$

The joint Hermitian kernel has entries

$$
K_{(n,A),(m,B)}
=
\sqrt{p(n,A)p(m,B)}
\beta^{\sum_j|Q_j|}
\mathcal M(Q_0,Q_1,Q_2,Q_3),
$$

restricted here to the leading 16 unresolved orbit families. This single
matrix retains compatibility and cancellation. Computing a trace norm for
each orbit separately and adding the answers is the triangle relaxation
that the joint program is meant to test.

## Exact link arithmetic

For one link, averaging over the random signs first imposes a parity
matching condition: the hidden permutation must map the odd column degrees
of the left support to the odd row degrees of the right support. The
remaining permutation average is a permanent.

The new evaluator computes this moment exactly. The odd-parity permanent is
small because the Fourier supports have dose at most six. The even-parity
block differs from the all-ones matrix only in touched rows and columns, so a
special-row/special-column expansion avoids enumerating all (q!)
permutations. The implementation was checked against all signed
permutations for 240 randomly selected support pairs at (q=2,4), and
against an exact (1/32) case at (q=32).

## Why the old direct model cannot answer the question

Every one of the leading 16 unresolved orbits has a block degree at least
five. At (q=2), each block contains only four coordinates. Therefore none
of the 16 families exists in the old direct physical-law model.

At (q=4), all leading profiles exist, but the full coordinate-uniform law
on the 30 Perron-support occupation states has dimension

$$
\sum_{n\in\operatorname{supp}\rho}
\prod_{j=0}^3 {16\choose n_j}
=18{,}904{,}064.
$$

A dense direct calculation is consequently the wrong representation.

## First sparse compatible law

As a tractable exact baseline, take six sampled coordinate configurations
for each of the 30 occupation states and search 24 deterministic seeded
candidates. Select by:

1. number of activated leading orbits;
2. number of activated profile-split entries; and
3. attenuated joint trace norm at (q=4).

The selected law has a 180-dimensional basis. It activates 26 of the 64
profile-split entries and all 16 leading orbits. The same configurations are
then embedded into the upper-left (4\times4) coordinates of each
(32\times32) block and reevaluated with exact (q=32) plant moments.

| calculation | separate attenuated norms | joint attenuated norm | joint / separate |
|---|---:|---:|---:|
| (q=4), (N=16) | (1.95832159156\times10^{-4}) | (1.46122575757\times10^{-4}) | (0.746162307493) |
| embedded (q=32), (N=1024) | (1.11522671396\times10^{-9}) | (9.52512381241\times10^{-10}) | (0.854097529510) |

Thus shared-law cancellation removes about 25.4 percent of the separate
sum at (q=4), and 14.6 percent for this sparse (q=32) embedding.

## What this decides

It establishes three facts:

- the 16 leading families can coexist under one legal physical law;
- their joint contribution is materially smaller than the sum of their
  separately evaluated contributions; and
- exact configuration-level (q=32) moments are now computationally
  available without signed-permutation enumeration.

It does not yet decide the hard instance. The sparse law is deliberately
low rank, the occupation weights are already built into its matrix, and its
trace norm is not a vector of scalar theorem coefficients that may be
inserted into the current Perron ledger. The final eigendecomposition is
also ordinary floating point, not interval arithmetic.

In particular, the small absolute (q=32) value is evidence only that this
sparse embedding is a weak witness. It is not evidence that all compatible
physical laws are weak.

## Next mathematical target

Replace the sparse conditional law by a compressed high-rank law. The next
calculation should reproduce the coordinate-uniform (q=4) answer without
forming its 18.9-million-dimensional matrix, then use the same verified
representation at (q=32). It should cover the leading 51 orbits, or a
natural family carrying at least 90 percent of the unresolved Perron
contribution.

Only that high-rank calculation can trigger the declared gate:

- if a certified compatible lower witness closes the (1/3) margin, pivot;
- if it stays below the gate, attempt the shared upper contraction; and
- do not return to independent orbit-by-orbit polishing.

## Reproduction

Run:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 \
  searches/joint_impact_sparse_q4.py --write-artifact
```

The committed result is
`artifacts/joint_impact_sparse_diagnostic.json`. The Round 4 test suite
regenerates it and compares it byte for byte.
