# Confidence report: reverse-tree contraction

Date: 2026-07-14

Verdict: **validated only with the round-two repairs in this directory**. The round-one mixed-component inference is not correct, and literal vertical multiplicativity also fails after distinct Fourier labels are imposed. REPAIRED_REVERSE_TREE_CONTRACTION.md uses an exhaustive global dichotomy plus a bounded-mass Walsh expansion of the distinct-label mask to restore the finite-skeleton $N^{-1/2}$ bound without changing the dose exponent.

This verdict is about mathematical correctness of the reverse-tree interface. It does not mean the bound is quantitatively useful at realistic sizes.

## 1. What the stress test found

1. The frozen round-one deterministic suite reproduces under its recorded Python 3.13.13 environment.
2. The repaired bilateral lemma omitted an injectivity or collision hypothesis. A two-coordinate collision violates the stated diagonal-only bound by $\sqrt2$.
3. Collision-aware PSD packing repairs that defect with a within-base fiber factor. For diagrams of at most twelve vertices, it is only an absolute factorial and changes neither the $N$ exponent nor the single insertion/Bessel charge.
4. The Hilbert-valued local reverse update is correct once its operator majorant and preimage spaces are stated explicitly.
5. Round one incorrectly treated ordinary Hilbert auxiliary control as equivalent to the projective control required by separately carried singleton components. A normalized Hadamard edge gives a concrete $\sqrt N$ norm gap.
6. The global dichotomy repair avoids that conversion entirely.
7. Distinct-coordinate masks couple disconnected components and invalidate literal vertical multiplicativity. An exact $N=4$ product witness exceeds the claimed $1/N$ two-component value. Walsh expansion costs only a finite diagram constant and restores termwise multiplicativity.
8. A literal Hilbert-vector frontier would charge unresolved identity wires by $\sqrt N$. The reverse induction must be operator-valued until a unit open-boundary input is fixed; the local Hilbert lemma then applies uniformly and preserves operator norm.
9. Exact parity-support enumeration confirms that the factorial collision fiber and the falling-factorial insertion mass are two aspects of one sum, so the collision repair does not duplicate the dose charge.
10. Conditioning on the interpolation variables leaves only bounded local vertex multipliers, not an arbitrary joint label weight. A tree-level binomial recurrence bounds all adaptive marked-time assignments by $(2D)^v$ without an outcome-width factor.

## 2. The replacement argument

There are exactly two cases.

### All-singleton case

First Walsh-expand the same-layer distinct-coordinate masks. If every graph component contributes at most one vertex to every physical amplitude entry, group graph coordinates by physical entry and keep a projective frame invariant throughout. Vertical injective norms multiply in each expanded term, and a four-layer component supplies $N^{-1/2}$ through the natural layer cut. Cross-component entanglement inside an entry is already included in the grouped party; summing the mask expansion costs only $C_G$.

### All-assigned Hilbert case

Otherwise, assign every component to a maximum-occupancy entry, including components whose maximum occupancy is one. For $k_C$ selected vertices,

$$
|B_C|-e_C=v_C-k_C-e_C\le1-k_C.
$$

Weak $k_C=1$ assignments cost at most one. At least one strong $k_C\ge2$ assignment exists and supplies $N^{-1/2}$. Because every component is now oriented as a sliced Hilbert map, the reverse induction remains in the Hilbert/operator regime and never asks a Hilbert auxiliary to obey a projective norm bound. Open ancestor coordinates are carried as operator-norm-one boundary wires, not vectorized in Hilbert--Schmidt norm.

Collision-aware complete-frame packing removes every physical outcome once. Opposite-entry graph fibers cancel against sliced covariance diagonals, while insertion-slot fibers remain inside the original falling-factorial mass.

## 3. Evidence classification

| Requirement | Evidence | Status |
|---|---|---|
| Exact graph spectra and sliced diagonals | Algebra plus frozen checks | Proved and reproduced |
| Layer-cut $N^{-1/2}$ for a four-layer component | Forest/rank proof plus exact/random checks | Proved and reproduced |
| Collision-aware frame packing | PSD proof plus a saturating regression | Proved |
| Operator-valued reverse frontier | Fix a unit boundary input, apply Hilbert duality/reshape, then take the operator supremum; identity-wire and adaptive stresses | Proved after frontier repair |
| Distinct-label mask | Exact Walsh expansion with coefficient mass at most $2^{P_G}$ | Proved; exact regression detects omission |
| All-singleton projective case | Grouped-entry projective induction and termwise vertical injective multiplicativity | Proved after mask expansion |
| Strong-plus-weak all-assigned case | Sliced diagonal inequality and pure Hilbert induction | Proved |
| Dichotomy exhaustiveness | Immediate logic; 69,632 exact placements and 100,000 random systems | Proved and stress-tested |
| No outcome-width or depth factor | Hilbert/projective reverse recurrences and stochastic unmarked nodes | Proved |
| Single insertion/Bessel charge | Exact base/ordered-mark identity plus exhaustive rational support enumeration | Proved and reproduced |
| Fixed-skeleton contraction $C_GN^{-1/2}$ | Two exhaustive analytic cases | Proved after repair |
| Interpolation handoff | Conditioned local vertex factors, antisymmetrization, and exact adaptive marked-time recurrence | Proved and independently re-audited |
| Transcript bound $C(1+D)^{12}/\sqrt N$ | Interpolation handoff plus repaired contraction | Supported after integration audit |
| Passive floor $\Omega(N^{1/24})$ | Algebraic consequence of the transcript bound | Supported |

The computational suite also covers arbitrary entangled selected tensors, two adaptive levels, all 256 minimal-chain placements, 2,720 mixed chain/edge placements, all relative all-singleton placements of two chains, and selected $N=4$ mixed placements. These are falsification checks, not the basis of the universal proof.

## 4. Important qualification

The original round-one proof should not be cited unchanged. It contains four real local defects:

- the unqualified collision statement; and
- the Hilbert-to-projective singleton ancillary inference;
- the omission of the distinct-label mask before vertical tensor multiplication; and
- the literal Hilbert-vector formulation at frontiers with unresolved identity wires.

The result is supported only with the collision-aware lemma, the distinct-label Walsh expansion, and the global dichotomy replacement in this directory. Integrating the repair must preserve the two cases exactly; reverting to a simultaneous assigned/projective component split reopens the gap.

## 5. Realistic-size verdict

The validated asymptotic bound is quantitatively useless near $N=1024$. At dose six,

$$
{(1+6)^{12}\over\sqrt{1024}}
={7^{12}\over32}
\approx4.33\times10^8
$$

before the unknown absolute constant. No honest constant extraction from the present proof can turn that into transcript distance below $1/3$.

The next phase must therefore seek a qualitatively stronger contraction or a different hard instance. Minor bookkeeping improvements are not enough.
