# Round 2 handoff

Date: 2026-07-14

The authoritative full summary is
`../open_problem_forr4_passive_floor_consolidation_round_2/ROUND_2_SUMMARY.md`.

## Inherited theorem status

- The passive \(\Omega(N^{1/24})\) result is accepted only with the
  collision-aware packing lemma, distinct-label Walsh expansion, global
  dichotomy, operator-valued frontier, and interpolation handoff in the
  round-two folder.
- The original round-one mixed-component proof is false as written.
- The repaired asymptotic constant is unusable at \(D=6,N=1024\).

The frozen mathematical baseline is round-two commit `788c826`; subsequent
round-two edits only consolidate status and handoff documentation.

## Inherited finite-size numbers

For the attenuated signed-permutation plant at \(N=1024\), \(q=32\), and
\(\beta=5/6\):

- promise-conditioning loss: \(0.02740656\);
- attenuated minimal one-batch contribution: \(0.14556864\);
- minimal plus promise: \(0.17297520\);
- remaining margin below \(1/3\): \(0.16035813\).

Every degree-eight profile has a fixed-split local bound.  The unordered
joint occurrence square-mass maxima at degrees \(4,6,8,10,12\) are
\(81,160,126,36,1\).

## Main obstruction

The splitwise row/column estimate for \((3,1,1,3)\) spends \(0.571263\)
after attenuation on one safe occupation ledger.  It must improve by a
factor below \(0.280708\), even if every other higher sector vanished.

This loss appears artificial.  The balanced endpoint pair frame satisfies

$$
AA^*={q^2+2\over2}I_N.
$$

One same-side joint slice has exact uniform coefficient

$$
{q^2+2\over q^3(q^2-1)}=3.0607\times10^{-5}
$$

at \(q=32\).  The alternating uniform slice is \(0.471592\) at \(q=2\)
and \(0.0642009\) at \(q=4\), but its general-\(q\) and weighted behavior
remain open.

## Known invalid shortcuts

- Entrywise matching probabilities do not control operator norm.
- Termwise adaptive \(\ell_1\) summation overshoots by hundreds.
- A universal near-unit two-copy square function is false.
- The first proposed Walsh-convolution reduction of the alternating Gram
  does not reproduce the exact \(q=2,4\) spectrum; see
  `notes/alternating_double_endpoint_target.md`.
