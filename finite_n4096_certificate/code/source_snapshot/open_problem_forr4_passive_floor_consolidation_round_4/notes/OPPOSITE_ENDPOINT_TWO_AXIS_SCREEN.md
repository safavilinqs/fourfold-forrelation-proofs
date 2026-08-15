# Opposite-endpoint two-axis screen

Date: 2026-07-15

## Decision

Close the local physical-witness search on the leading opposite-endpoint
orbit. The best valid two-axis non-invariant law remains below the scalar
kill gate. Do not spend another Round 4 project tuning this orbit; return to
the global shared-contraction target.

## Why this screen was needed

The best row-symmetry-invariant law has coefficient

$$
0.0395996495754,
$$

while the leading-orbit scalar gate is

$$
0.0414623182965.
$$

At (q=4), the unrestricted mixed-orbit optimizer is much larger than the
invariant vertical law. This left open the possibility that a simple valid
non-invariant direction at (q=32) would cross the gate and force an
immediate pivot.

## Physical family

Keep the cubic selected-pair difference vertical and the triple translation
orbit vertical. In the quintic selected-pair difference, mix:

- the correlated vertical law, whose equal-difference mass is optimized;
  and
- a uniform horizontal difference with total mass (t).

This is one legal diagonal physical law. It is not the invalid fourteen-class
extension to arbitrary nonsymmetric laws. Independent row-label symmetries
inside the vertical and horizontal axes reduce the exact frequency sum to
25 blocks.

The reduction is checked at (q=4) against the inherited unreduced formula
for a genuinely non-invariant law. For equal mass (2/5) and horizontal
mass (1/20), both calculations give

$$
0.0837940294383.
$$

## Result at (q=32)

The optimized family gives:

| quantity | value |
|---|---:|
| invariant equal-difference mass | (0.0192515087131) |
| invariant coefficient | (0.0395996495754) |
| horizontal mass | (0.000647887695545) |
| best two-axis coefficient | (0.0396118487001) |
| absolute improvement | (0.0000121991247051) |
| relative improvement | (0.0308061%) |
| remaining gate headroom | (0.00185046959644) |

The horizontal component is real but tiny. It improves the invariant witness
by only about three hundredths of one percent and closes less than one
percent of the remaining gap to the kill gate.

## Scope

This does not prove that no arbitrary physical law crosses the gate. It does
rule out the most immediate non-invariant direction suggested by the exact
(q=4) optimizer, using a reduction that has been validated without a
symmetry assumption error.

Round 4 therefore applies its stop rule:

1. the invariant leading-orbit search did not kill the scalar route;
2. the first valid non-invariant extension did not materially change it;
3. the compatible 16-orbit sparse law showed genuine cross-family
   cancellation; and
4. the next work is one high-rank shared upper contraction, not another
   leading-orbit lower-witness variant.

## Reproduction

Run:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 \
  searches/opposite_endpoint_two_axis_mixture.py --write-artifact
```

The committed artifact is
`artifacts/opposite_endpoint_two_axis_screen.json`. The regression checks the
unreduced (q=4) identity and regenerates the (q=32) optimization.
