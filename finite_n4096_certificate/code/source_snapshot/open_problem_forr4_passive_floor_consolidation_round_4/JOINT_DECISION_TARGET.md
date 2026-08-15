# Joint hard-instance decision target

Date: 2026-07-15

Status: historical initial decision framework. It retained the
signed-permutation route provisionally and led to the q64 program. The live
gate is notes/Q64_SHARED_QUINTIC_AND_ADAPTIVE_ACCEPTANCE.md: $N=4096$, 264
proved entries, and 164 entries in the next shared class.

## One-sentence target

Decide whether the attenuated signed-permutation hard instance is still a
viable route to a rigorous passive-dose-greater-than-six result at a credible
finite size by testing the dominant unresolved families under one compatible
physical law and one shared contraction; retain it only if that global test
leaves enough margin for certification and adaptivity, otherwise pivot
immediately to a replacement witness.

## Why this replaces orbit-by-orbit closure

The current $N=1024$ total is only a diagnostic:

$$
0.333132605485488 < \frac13
$$

by $0.000200727847845$. Only 40 of 888 balanced entries have
arbitrary-law upper bounds. The other 848 entries form 224 unresolved
complement/reversal orbits. Resolving the four-entry opposite-endpoint orbit
would still leave 844 entries in 223 orbits.

The unresolved Perron contribution is concentrated enough to support a
bounded joint test:

- the leading orbit accounts for 13.6 percent;
- the leading 16 orbits account for 51.1 percent; and
- the leading 51 orbits account for 90.5 percent.

These percentages are routing diagnostics computed from the present
coefficient map and Perron sensitivities. They are not theorem bounds.

## The joint lower-witness test

Construct one legal diagonal probe law whose latent choices are shared across
the leading unresolved families. For each induced cut, compute a rigorous
lower bound on the corresponding physical coefficient. Reoptimize the scalar
ledger with the whole compatible coefficient vector, not with independently
chosen per-orbit witnesses.

Start with the leading 16 unresolved orbits and extend to the leading 51 if
the first stage does not decide the route.

The first sparse compatible-law baseline is complete. It activates all 16
leading orbits and proves that joint cancellation is numerically material,
but its low rank makes it too weak to apply this gate. The operative next
target is the compressed high-rank law in
notes/JOINT_PHYSICAL_LAW_DIAGNOSTIC.md, not further tuning of the sparse
sample.

The last local kill screen found a valid non-invariant two-axis coefficient
of \(0.0396118487001<0.0414623182965\). This closes the planned
opposite-endpoint lower-witness variants without deciding arbitrary laws.
Proceed to the shared contraction rather than opening another orbit-local
search.

The first native-\(q=32\) higher-rank symmetry class is closed as well. A
960-configuration common-row-translation law has nonzero exact moments in
only 6 of the 51 frontier orbits and is eleven orders of magnitude below the
current frontier scale. Any further compatible lower law must enforce the
signed-permutation parity matches by construction; random translation-orbit
tuning is not a live route.

If a certified compatible coefficient vector forces the scalar ledger to at
least $1/3$, the independent scalar-ledger proof architecture is closed.
This would not prove that the hard instance itself is easy for passive
protocols; it is the predeclared signal to stop using it for the finite-size
paper result and move to the replacement scorecard.

## The shared upper-contraction test

If the lower-witness test stays below $1/3$, attempt one joint contraction
covering at least the leading 51 orbits, or another natural family carrying
at least 90 percent of the unresolved Perron contribution. The contraction
must retain the common physical law and cancellation between cuts; summing
independently optimized scalar orbit bounds does not qualify.

The target now has an exact finite occupation representation. The 51-orbit
frontier collapses to six unordered profile patterns on 241 edges among 125
occupation states, split into eight connected components. The first two
patterns, `(5,3,1,1)` and `(5,3,3,1)`, carry 81.8 percent of the frontier's
current Perron contribution. This does not prove a contraction, but it
replaces the apparent 198-entry problem by one bounded componentwise operator
problem. See notes/SHARED_FRONTIER_STRUCTURE.md.

Retain the signed-permutation route only if the resulting one-batch bound is
rigorous or interval-certifiable and leaves at least $10^{-3}$ reserve below
$1/3$, unless a proved adaptive-lift estimate justifies a different reserve.
The $10^{-3}$ value is a project decision threshold, not part of the final
theorem.

## Retain or pivot

Retain the current witness only after the shared contraction passes. Then:

1. replace floating-point Perron optimization by outward-rounded intervals;
2. complete the residual low-impact families globally;
3. prove the adaptive posterior-selection lift; and
4. write the active six-dose protocol and experimental resource row.

Pivot if the compatible lower witness reaches $1/3$, the shared contraction
cannot cover the impact frontier with reserve, or the surviving size is not
experimentally credible. The first named replacement is the quadratic-bent
exact plant because it removes the promise-conditioning loss; it must be
instantiated and scored rather than left as a placeholder.

That first score now exists. The exact promise leaves an optimistic
\(0.0518213\) higher-sector budget, but a common \(1/N\) coefficient vector
fails. The replacement needs a complete vector averaging roughly
\(0.542524/N\), with \(1/(2N)\) as a simple sufficient target. It remains the
first fallback but is not promoted before a weighted arbitrary-split
contraction reaches that scale.

## Non-goals during the decision

Do not spend the decision phase on:

- another long sequence of independent orbit theorems;
- improving the inherited $N^{1/12}$ exponent;
- active dose five;
- paper exposition or Obsidian reformatting; or
- unsupported claims that 4096 sign modes are experimentally feasible.

The broader scientific questions about passive-frame realizability and
passive protocols remain valuable, but they are pursued here only when they
decide the finite-size witness or supply its replacement.
