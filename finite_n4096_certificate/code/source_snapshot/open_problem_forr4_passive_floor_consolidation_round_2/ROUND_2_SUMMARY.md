# Round 2 summary and handoff

Date: 2026-07-14

## Executive outcome

Round two had two goals.

1. Stress-test the reverse-tree contraction behind the passive
   \(\Omega(N^{1/24})\) four-forrelation lower bound.
2. After that, seek a quantitatively stronger route that excludes passive
   hard dose six near \(N=1024\).

The first goal succeeded, but only after replacing four false or incomplete
round-one steps.  The second goal did not reach a lower bound.  It produced
an exact hard-instance program, several sharp local contractions, an exact
one-batch benchmark, and a precise remaining obstruction.

The correct final claims are:

- The asymptotic \(\Omega(N^{1/24})\) lower bound is supported in the
  repaired form audited here.
- The original round-one proof must not be cited unchanged.
- No result in this folder proves passive hard dose greater than six at
  \(N=1024\).

The full regression command is:

```sh
./run_round2_checks.sh
```

## 1. What changed in the asymptotic proof

The audit found four real defects.

| Defect | Falsification | Required repair |
|---|---|---|
| Collision-free diagonal packing was applied to noninjective fibers | A two-coordinate collision loses \(\sqrt2\) | Collision-aware PSD packing, with its factorial fiber inside the original insertion/Bessel sum |
| A Hilbert auxiliary was later treated as though it had projective norm one | A normalized Hadamard edge has a \(\sqrt N\) Hilbert/projective gap | Exhaustive global dichotomy: all-projective singleton case or all-assigned Hilbert/operator case |
| Distinct Fourier labels were omitted before multiplying graph components | An exact \(N=4\) masked product exceeds literal multiplicativity | Walsh-expand the distinct-label mask; absorb only a finite diagram constant |
| Open identity wires were vectorized as unit Hilbert vectors | Their Hilbert--Schmidt norm is \(\sqrt N\) | Carry an operator-valued frontier and fix a unit boundary input before applying the local Hilbert lemma |

The interpolation handoff and adaptive marked-time ledger were then
re-audited.  They give the transcript estimate

$$
\operatorname{TV}\le C(1+D)^{12}N^{-1/2},
$$

and hence the passive floor \(\Omega(N^{1/24})\).  The proof depends on all
four repairs simultaneously.  The authoritative files are:

- `REPAIRED_REVERSE_TREE_CONTRACTION.md`;
- `AUDIT_REPAIRED_CONTRACTION.md`;
- `CONFIDENCE_REPORT.md`;
- `notes/interpolation_handoff_audit.md`; and
- the falsification history in `notes/counterexample_log.md`.

## 2. Why the asymptotic proof cannot solve the realistic-size problem

At \(D=6\), \(N=1024\), the repaired asymptotic ledger contains

$$
{7^{12}\over32}\approx4.33\times10^8
$$

before its absolute constant.  No constant cleanup can make this less than
the \(1/3\) transcript-distance threshold.  A qualitatively different
contraction or hard instance is required.

## 3. Realistic-size assets proved in round two

### Exact and attenuated signed-permutation plant

The Maiorana--McFarland signed-permutation orbit gives exact Boolean
four-forrelation inputs at \(N=q^2\).  Independent sign attenuation by
\(\beta=5/6\) gives mean \(\pm\beta^4\), exact conditional variance

$$
{1-\beta^8\over N},
$$

and conditioning loss at most \(0.0274066\) at \(N=1024\).

### Exact minimal one-batch benchmark

For one passive batch, the minimal four-vertex sector has exact dose-six
coefficient

$$
F_6={2337\over256}+{3\sqrt2\over8}=9.659236\ldots.
$$

Without attenuation this gives \(F_6/32=0.301851\ldots<1/3\).  With the
attenuated plant, the minimal contribution plus promise conditioning uses
\(0.172975\ldots\), leaving \(0.160358\ldots\) for all higher sectors and
adaptivity.

### Weighted contraction toolkit

Round two proved:

- universal diagonal-weighted Gram-link contraction;
- block-coherent weighted three-link path contraction;
- Gram-dressed endpoint contraction for every endpoint degree;
- cubic and quintic middle-decoration bounds;
- fixed-split bounds for all six degree-eight double-cubic profiles; and
- exact slice energies for the formerly hard endpoint and adjacent
  double-cubic cases.

Every degree-eight profile is now locally bounded at fixed occurrence
split.  Nineteen of the seventy odd profiles through total degree twelve
have explicit local classifications; fifty-one degree-ten or degree-twelve
profiles remain locally incomplete.

### Joint occurrence packing

For an odd degree profile \(a=(a_1,a_2,a_3,a_4)\), the unordered joint
occurrence square mass satisfies

$$
\mathcal S_a\le
\max_{\sum_bM_b\le12}\prod_b\binom{M_b}{a_b}.
$$

The largest square masses at degrees \(4,6,8,10,12\) are respectively

$$
81,\quad160,\quad126,\quad36,\quad1.
$$

This removes an exponential split count and any labeled-mark factorial.
It does not promote automatically to terminal trace norm.

## 4. Approaches falsified or rejected

These failures should remain in the record because they constrain any
future proof.

- Entrywise hidden-label matching is false as an operator theorem.  A
  crossing minimal-chain cut has scale \(N^{-1/2}\), not \(N^{-1}\).
- Optimizing and summing marked-time assignments in \(\ell_1\) gives
  coefficient \(8730\) at dose six.
- A universal near-unit two-copy square-function lemma is false.  The exact
  two-node witness needs constant \(\sqrt6\) at \(N=1\) and
  \(\sqrt{3/2}\) after the relevant \(N=2\) embedding.
- Unweighted link operator norms can grow even for promising alternative
  hard orbits.  Physical diagonal weights must remain in the theorem.
- The first fixed-split row/column-energy sum for \((3,1,1,3)\) is far too
  coarse: after attenuation it spends \(0.571263\) on one safe occupation
  ledger, versus only \(0.160358\) total remaining TV margin.

## 5. Most important unfinished calculation

For the balanced double-endpoint slice, define

$$
A(i;E,b)=M_{3,1}(\{i\}\cup E,b).
$$

The endpoint pair frame is exactly tight:

$$
AA^*={q^2+2\over2}I_N.
$$

When the two middle singleton coordinates are on the same side, the full
uniform Schur coefficient is exactly

$$
{q^2+2\over q^3(q^2-1)},
$$

only \(3.06\times10^{-5}\) at \(q=32\).  The alternating placement has
uniform coefficients \(0.471592\) at \(q=2\) and \(0.0642009\) at \(q=4\).
This is strong evidence of a compound tight-frame cancellation, but no
arbitrary-diagonal or \(q=32\) alternating theorem was proved.

The best immediate target is therefore:

> Prove or falsify a weighted compound-frame contraction for the
> alternating \((3,1,1,3)\) Schur lift that preserves the shared passive
> occupation law and improves the current splitwise sum by a factor below
> \(0.280708\).

## 6. Round-three handoff

The next round should proceed in this order:

1. diagonalize the alternating double-endpoint Gram through its
   translation/association scheme and test the candidate weighted norm;
2. if it survives, combine cut-dependent degree-six and degree-eight
   constants with the exact joint occupation law in a one-batch
   \(N=1024\) budget;
3. extend only the quantitatively relevant degree-ten and degree-twelve
   sectors, exploiting their small occurrence masses;
4. attack the adaptive outcome-selection interface after the one-batch
   budget is genuinely below threshold; and
5. switch hard instances only if the weighted double-endpoint theorem is
   falsified or the complete finite-size ledger cannot fit below \(1/3\).

The round-three initialization documents contain explicit success gates,
stop conditions, and the inherited evidence map.
