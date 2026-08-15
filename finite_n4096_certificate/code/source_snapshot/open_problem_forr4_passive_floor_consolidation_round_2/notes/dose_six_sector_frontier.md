# Exact dose-six sector frontier

Date: 2026-07-14

Status: bookkeeping frontier proved.  Every profile through total Fourier
degree eight now has a local fixed-split contraction.  The three formerly
hard record-one double-cubic profiles are controlled by slice-energy
estimates.  This does not yet combine the sector contractions with the
joint occurrence mass or the adaptive lift.

## 1. Exact profile list

Hypothesis sensitivity forces an odd signed-permutation record through all
three links.  In every block, the Fourier symmetric-difference degree is
therefore odd.  Write

$$
n_b=1+2k_b,
\qquad k_b\ge0.
$$

A hard-dose-six ket/bra pair has total Fourier degree at most twelve, so

$$
k_1+k_2+k_3+k_4\le4.
$$

There are exactly

$$
{8\choose4}=70
\tag{1.1}
$$

degree profiles.  Their counts by total degree are

| degree | profiles |
|---:|---:|
| 4 | 1 |
| 6 | 4 |
| 8 | 10 |
| 10 | 20 |
| 12 | 35 |

Odd link-record sizes are bounded independently by the smaller adjacent
block degree.  Keeping all such choices produces 130 record triples
across the 70 profiles.

## 2. Closed local profiles

The current contractions close or explicitly bound nineteen degree profiles locally:

- the minimal \((1,1,1,1)\) profile;
- every profile with a single decorated endpoint, of degrees
  \(3,5,7,9\), by the Gram-dressed tail theorem; and
- the two cubic middle profiles \((1,3,1,1)\) and \((1,1,3,1)\), by the
  L-shape rank--Frobenius bound; and
- the two quintic middle profiles \((1,5,1,1)\) and \((1,1,5,1)\), by the
  conditional-permutation rank--Frobenius bound; and
- three double-cubic profiles, by compatible-entry rank--Frobenius bounds;
  and
- the other three double-cubic profiles, by endpoint/L-shape occurrence
  slice energies.

Consequently all five profiles of total degree at most six have a local
coefficient no larger than \(1/q\) for every fixed occurrence split, and
all ten degree-eight profiles are locally bounded.  The record-three
sector of \((1,3,3,1)\) still has the coarse coefficient about \(6/q\),
and the double-endpoint record-one coefficient is too large for a final
finite-size sum, but neither is structurally open at fixed split.

At higher degree, 16 of the 20 degree-ten profiles and 31 of the 35
degree-twelve profiles remain locally open.  Record sizes three and above
carry much smaller matching coherences, but record-one sectors exist in
this frontier and must be treated explicitly.

## 3. Why this is not yet a finite-size budget

The 70 profile count is not a factor that may simply multiply \(1/q\).
Different occurrence assignments share one passive probe law, and their
common base supports are packed by complete-frame Bessel identities.  The
minimal profile alone demonstrates the issue: its 16 fixed cuts combine
to the sharp occupation coefficient

$$
{2337\over256}+{3\sqrt2\over8},
$$

not to either one cut or sixteen independent unit masses.

The joint unordered-support packing inequality now bounds each profile's
complete occurrence square mass by

$$
\max_{\sum_bM_b\le12}\prod_b\binom{M_b}{a_b}.
$$

The next calculation must combine that mass with the cut-dependent local
contractions before taking terminal absolute values.  Re-expanding the
already closed cases entrywise, or replacing every cut by its worst local
coefficient, discards the required gain.

Reproduction: searches/dose_six_sector_frontier.py and
tests/joint_occurrence_profile_packing.py.
