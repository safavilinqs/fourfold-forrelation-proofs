# Realistic-size program: exact plant after the first falsification

Date: 2026-07-14

Status: the entrywise matching-budget theorem proposed in the first version is false. The exact signed-permutation plant remains promising, but it requires a spectral recording inequality.

Round-two disposition: closed as a chronological research log, not as the
current plan.  Sections intentionally retain falsified intermediate
targets.  ROUND_2_SUMMARY.md gives the authoritative end-of-round state;
continuing work is initialized in the round-three folder.

## 1. Counterexample to the entrywise target

Take the minimal four-layer chain and split its vertices between two amplitude sectors of a single passive dose-two batch.

For the exact plant, the positive-minus-negative minimal moment is twice the normalized Hadamard chain tensor $T$. The adjacent-pair flattening has nuclear norm $N$, but either crossing-pair flattening has rank $N^2$, singular values $N^{-1/2}$, and nuclear norm $N^{3/2}$.

Give each sector total diagonal mass $1/2$, uniformly across its $N^2$ coordinate pairs. The exact one-batch Helstrom values are

$$
\operatorname{TV}_{\rm adjacent}={1\over N}={1\over q^2},
\qquad
\operatorname{TV}_{\rm crossing}={1\over\sqrt N}={1\over q}.
$$

The entrywise three-permutation record budget for one occurrence per block is

$$
{2\over q^3}.
$$

The crossing ratio to the entrywise budget is $q^2/2$, equal to $512$ at $q=32$. Therefore no dimension-independent constant can turn the entrywise event budget into a transcript bound.

The missing factor comes from coherent summation over the character coordinates that are not recorded by the hidden permutation label sets. This is exactly the operator amplification warned about in the initial program.

## 2. What remains useful

The combinatorial search still identifies the dangerous block occurrence allocation $(2,4,4,2)$ at dose six and quantifies the hidden-label part of the problem. It is a diagnostic, not a certificate.

The adjacent placement improves to $1/N$, but the crossing placement saturates the generic $N^{-1/2}$ contraction. At $N=1024$, the explicit passive dose-two advantage is therefore at least $1/32$.

## 3. Revised theorem target

Replace the scalar match union bound by a spectral record decomposition. For each triple of match sizes $(r_1,r_2,r_3)$ and block occurrence profile, form the exact character-weighted matching operator from the permutation moment formula. The desired theorem is

$$
\operatorname{TV}(Q_+^{\mathcal T},Q_-^{\mathcal T})
\le
\sum_{r_1,r_2,r_3\ge1}
\mathcal B_{r_1,r_2,r_3}(D,q),
$$

where each $\mathcal B$ is an operator or completely bounded norm after passive complete-frame weighting, not an entrywise probability.

The finite-size milestone is unchanged: the rigorously summed right side must be below $1/3$ for $D=6$, $q=32$.

## 4. Exact one-batch benchmark

For the minimal four-vertex chain, keeping each ket/bra cut's exact binary rank reduces the dose-six one-batch problem to a concave optimization over four block occupations. It has the rigorous optimum

$$
F_6={2337\over256}+{3\sqrt2\over8}
\approx9.65923633589.
$$

The corresponding transcript contribution at $N=1024$ is

$$
F_6/32\approx0.301851135497<1/3.
$$

Global optimality is certified by an exact KKT supporting hyperplane over all 210 integer occupation states of total dose at most six. This establishes the desired realistic-size inequality for one batch and the minimal sector only.

The margin, about $0.03148$, is too small to add adaptive and higher-sector errors by a coarse triangle inequality.

## 5. Next concrete calculations

1. Diagonalize the one-link size-$r$ character-weighted matching operator under the action of the affine group. The $r=1$ sector must reproduce the exact Hadamard unitary rather than the scalar $1/q$ event probability.
2. Retain all minimal-chain groupings: adjacent gives $1/N$, while crossing gives the sharp $1/\sqrt N$ obstruction.
3. Extend the rank-weighted occupation certificate across temporal partitions without summing marked-time assignments in $\ell_1$.
4. Enumerate all occurrence profiles with total size at most twelve and retain spectral ranks and singular values, not only event counts.
5. Search two-batch passive seesaws at $q=2,4$ for coherent amplification between record sectors.
6. Seek a reverse-tree inequality that square-sums orthogonal representation sectors across adaptive histories. An entrywise triangle inequality cannot meet the target.

## 6. Rejected adaptive ledger

Optimizing every marked-time placement independently and then summing absolute values is hopeless even after exact block-occupation masses are inserted. For six one-photon nodes, this safe ledger has coefficient exactly $8730$, overshooting the required $32/3$ by a factor about $818$.

This does not describe an achievable protocol. It isolates the needed theorem: temporal placements must be combined before absolute values, through an orthogonal square function, a martingale argument, or an exact tester norm.

## 7. Small-$q$ adaptive falsification

For the exact $q=2$, $N=4$ plant, the one-batch dose-two SDP gives TV $=1/2$. A restricted two-stage $1+1$ seesaw with exact SDP children converges to TV $=1/4$, with matching repaired child dual bounds for its final root. No adaptive amplification is seen in that slice.

This is encouraging but noncertifying. The $1+2$ child SDP is already large, and a global tester relaxation will need symmetry reduction before it is useful.

## 8. Quantitative bar

A uniform $C(D)/N$ bound is false. A bound of the sharp form $C(D)/\sqrt N$ can prove the desired result at $N=1024$ only if

$$
C(6)<{32\over3}\approx10.667.
$$

This is now an explicit constant problem, not an exponent problem.

Candidate milestones are therefore:

- a square-function bound whose exact dose-six value is below $10.667$;
- preferably, $C(D)=O(D)$ or $O(\sqrt{\binom D2})$ with a modest constant; or
- a sectorwise bound that is substantially smaller than its worst polynomial envelope at the specific integer dose six.

## 9. Placement-only square functions are insufficient

For the six one-photon-node partition, exact minimal-chain placement weights
satisfy

$$
\sum_pb_p=8730,
\qquad
\left(\sum_pb_p^2\right)^{1/2}
=\sqrt{14445/2}
\approx84.9853.
$$

Thus even replacing the rejected temporal $\ell_1$ ledger by an ideal
unweighted $\ell_2$ square sum misses the required coefficient $32/3$ by a
factor about $7.97$. A successful realistic-size proof must retain
representation-sector cancellation, joint probe incompatibility, or the
full adaptive tester norm; a generic placement martingale alone is not
enough.

Reproduction: tests/temporal_square_function_barrier.py and
notes/temporal_square_function_barrier.md.

## 10. Joint-probe square mass: valid algebra, false frame inference

The preceding placement-only calculation optimizes each marked placement
independently and is therefore not the right square function.  If all
placements retain their shared passive probe law, their exact summed square
mass is

$$
\mathcal S=\mathbb E\prod_{b=1}^4 M_b,
\qquad \sum_bM_b\le2D.
$$

Consequently

$$
\mathcal S\le(D/2)^4.
$$

At dose six this gives \(\sqrt{\mathcal S}\le9\), hence \(9/32<1/3\) at
\(N=1024\).  This is a rigorous occupation bound uniform over every dose
partition.

It does not imply the corresponding transcript bound.  Two nonadaptive
dose-one probes with Hadamard measurements give an exact ratio \(\sqrt6\)
between minimal transcript mass and \(\sqrt{\mathcal S}\) at \(N=1\).
The one-coordinate \(N=2\) embedding still requires
\(\sqrt{3/2}>32/27\).  Thus the proposed dimension-uniform near-unit frame
lemma is false.

Reproduction: notes/joint_probe_square_mass.md and
tests/joint_probe_square_mass.py.

## 11. Attenuated exact plant

Independently bias-flipping every coordinate of the exact plant with mean
\(\beta\) attenuates a degree-\(v\) hypothesis-dependent moment by
\(\beta^v\).  Every nonempty block derivative tensor of the exact
four-chain has squared norm exactly \(1/N\), which gives the exact formula

$$
\operatorname{Var}(F_{4,H}\mid X)={1-\beta^8\over N}.
$$

For \(N=1024\) and the rational choice \(\beta=5/6\), conditioning each
hypothesis onto its promised side costs at most \(0.013704\).

The previously displayed \(0.163041<1/3\) partial budget depended on the
false near-unit frame lemma and is rejected.  What survives is the exact
variance and promise-loss calculation.  The candidate now requires a
dimension-sensitive contraction retaining cut ranks or representation
sectors, plus a bound on all higher odd sectors.

This is a candidate hard instance, not yet a lower bound.  Its open tasks
are now sharply isolated: prove or falsify the temporal/frame square
contraction in the correct sector norm, then sum the attenuated higher odd
sectors.

Reproduction: notes/attenuated_exact_plant.md and
tests/attenuated_exact_plant_variance.py.

## 12. Dimension-sensitive two-batch replacement

For two fixed nonadaptive dose-one probes, all six minimal temporal
placements are entries of one averaged two-probe Hermitian operator.
Splitting that operator into the three \(2+2\) chain flattenings and
optimizing the two probes' block masses gives the sharp bound

$$
\|\Omega\|_1\le {1\over2\sqrt N}+{1\over4N}.
$$

Uniform block and coordinate masses attain the formula.  Relative to the
joint occupation square mass \(3/32\), its constant tends to
\(\sqrt{8/3}\), so exact cut geometry repairs much of the failed
constant-one proposal.

The same argument consolidates every fixed nonadaptive multi-batch
schedule into the one-batch dose-\(D\) Schur symbol by XOR-pushing forward
the product probe law.  It does not consolidate an outcome-selected child
law.  Preserving this dimension-sensitive trace norm through adaptive
preparation is now the precise minimal-sector boundary.

Reproduction: notes/dimension_sensitive_two_batch_contraction.md and
tests/two_batch_trace_norm.py.

## 13. Alternate exact orbit: quadratic bent functions

The uniform quadratic-bent orbit gives another pointwise exact plant:
\(Y=H_NX\) is Boolean, three independent pairs produce \(F_{4,H}=1\), and
first-block reflection produces \(-1\).

At \(N=16\), exhaustive enumeration of all 896 bent functions shows that
every nonzero complete link norm through six occurrences is smaller than
for the signed-permutation orbit.  The equal-degree norms improve from

$$
1,\ {8\over3},\ 4,\ 16,\ 32,\ {256\over3}
$$

to

$$
1,\ {8\over7},\ {12\over7},\ {48\over7},\
{96\over7},\ {256\over7}.
$$

This includes even-pair decorations, which were the unresolved part of the
signed-permutation sector proposal.

There is an exact growing collision.  A four-point affine plane
and the four-point complements inside all of its 3-flat extensions define
the same feature on every quadratic function.  The resulting
degree-four feature class has size \(N/4\), equal to 256 at \(N=1024\).
The canonical \((4,4)\) cross class cancels this multiplicity.  The next
endpoint sector has a growing *unweighted* norm: exact
association-scheme counting gives

$$
\|M_{5,1}\|_{\rm op}^2
={11182413\over64897}
$$

at \(N=1024\), hence norm \(13.126697\).  In general this unweighted norm
grows as \(\Theta(\sqrt N)\), rejecting a plain link-operator proof.

The physical diagonal weighting changes the conclusion.  Rows obey
\(M_{5,1}(A,y)=\mu_AH_{\oplus A,y}\) with
\(|\mu_A|\le2/(N-2)\).  Weighted Frobenius-to-nuclear conversion therefore
gives

$$
\|D_p^{1/2}M_{5,1}D_q^{1/2}\|_1
\le {2\over N-2}\sqrt{PQ}.
$$

Thus the affine multiplicity disappears at the passive Schur-symbol
level.  The quadratic-bent orbit remains viable, but only a fully weighted
sector analysis can certify it.

The weighted endpoint statement in fact holds for every degree:

$$
M_{a,1}(A,y)=\mu_AH_{\oplus A,y},
\qquad
\|D_p^{1/2}M_{a,1}D_q^{1/2}\|_1\le\sqrt{PQ}.
$$

It also holds for \(M_{1,a}\).  Therefore every outer link with a
singleton adjacent block is controlled uniformly, including deterministic
affine-flat relations.  The smallest unresolved quadratic-bent object is
the weighted two-sided internal sector \(M_{a,b}\) with \(a,b\ge2\).

Reproduction: notes/quadratic_bent_candidate.md,
notes/quadratic_endpoint_weighted_contraction.md,
searches/quadratic_bent_sector_spectra.py,
tests/quadratic_bent_collision_barrier.py, and
tests/quadratic_endpoint_weighted_bound.py.

## 14. Universal weighted link contraction

The endpoint proof is a special case of a general Gram factorization.  For
any planted pair and any support families,

$$
M(A,B)=\mathbb E[\overline{f_A}g_B]
$$

obeys

$$
\|D_p^{1/2}MD_q^{1/2}\|_1\le\sqrt{PQ}.
$$

This controls every odd-record and even-decoration sector on one link,
for both the signed-permutation and quadratic-bent plants, with constant
one.  Large unweighted sector norms are absorbed by the physical diagonal
mass.

The unresolved object is now specifically the three-link product

$$
M_1(A,B)M_2(B,C)M_3(C,D),
$$

where the shared middle labels must be charged once, not once per adjacent
Gram bound, and the inequality must remain valid under outcome-selected
child preparations.

Reproduction: notes/weighted_link_gram_contraction.md and
tests/weighted_link_gram_contraction.py.

## 15. Block-coherent weighted three-link path contraction

The three-link composition inequality is now proved for a block-coherent
flattening of an arbitrary fixed support sector and an arbitrary correlated
diagonal weight on that flattening.  If

$$
K(a,b,c,d)=M_1(a,b)M_2(b,c)M_3(c,d),
$$

then every one of its 16 whole-block flattenings has a weighted nuclear bound
that charges each full row or column multi-index once.  For the alternating
cut, the sharp factorization is

$$
(I\otimes M_1)\,\operatorname{diag}(M_2)\,(I\otimes M_3),
$$

which gives coefficient
\(\|M_1\|_{\rm op}\max|M_2|\|M_3\|_{\rm op}\).  The adjacent and
separated cuts follow from the universal link Gram lemma, and the middle
singletons use a two-link wedge rank--Frobenius bound.

For \(M_i=H_N\), the complete table exactly reproduces the minimal-chain
cut factors \(N^{-3/2}\), \(N^{-1}\), and \(N^{-1/2}\).  Thus the
three-link tensor geometry itself is no longer an unknown source of loss
in the singleton sector.

This does not yet prove the full higher-sector Schur contraction.  A
multi-coordinate symmetric-difference support in one physical block can
have some marks on the ket and others on the bra.  Such an occurrence-level
placement is an eight-variable Schur lift, not one of the 16 whole-block
cuts.  The remaining one-batch task is therefore both structural and
plant-specific: lift the once-weighted factorization to those internal
splits, enumerate
the compatible signed-permutation odd-record/even-decoration sectors,
insert their exact entry coherences, compound operator norms, and wedge
ranks, and sum the physical occupation masses at dose six.  Only after
that sum passes the \(N=1024\) threshold should the same contraction be
lifted through outcome-selected complete frames.

Reproduction: notes/three_link_weighted_path_contraction.md and
tests/three_link_weighted_path_contraction.py.

## 16. Cubic-decoration compatibility

The first complete signed-permutation decorated sector has now been
diagonalized.  An isolated \((1,3)\) record-one link obeys

$$
M_{1,3}M_{1,3}^*={q^2+2\over6}I_N.
$$

Its operator norm therefore grows as \(q/\sqrt6\); even decorations are a
real spectral effect, not merely duplicate physical labels.

Path compatibility changes the conclusion.  A degree-three middle block
with one odd row and one odd column must be a \(2\times2\) L-shape.  The
compatible restricted links have operator norm one and coherence
\(1/[q(q-1)]\).  Hence the block-coherent \((1,3,1,1)\) path gains an
extra \(1/(q-1)\) relative to the minimal chain.

At an endpoint there is no two-sided compatibility.  The orbit in which
all three coordinates share one hidden label gives a normalized
\((3,1,1,1)\) alternating flattening of exactly \(1/q\).  Thus cubic
endpoint decorations do not yield a constant-size distinguisher, but they
remain at the sharp \(N^{-1/2}\) scale.  In the attenuated candidate their
guaranteed improvement is \(\beta^2\), not another dimension factor.

The endpoint occurrence-split Schur lift has also been diagonalized.  The
worst balanced \(1|2\) split has normalized nuclear norm

$$
{4(q-1)(q^2-q+1)\over q^4\sqrt{q^2-1}},
$$

equal to \(0.003672\) at \(q=32\) and asymptotic to \(4/N\).  Direct
\(q=4\) diagonalization of all 32 endpoint placements confirms that none
exceeds the whole-block \(1/q\) orbit.

The next required calculation is the general-diagonal and middle-block
Schur lift, followed by a joint minimal-plus-cubic occupation sum;
whole-block compatibility alone is not an upper bound for every passive
ket/bra placement.

Reproduction: notes/signed_permutation_decoration_compatibility.md and
searches/signed_permutation_full_sector_spectra.py.

## 17. Gram-dressed endpoint theorem

The general diagonal-weight problem for endpoint decorations is now
closed.  If a base kernel has weighted nuclear coefficient \(\gamma\),
then repeating its rows/columns and Schur-multiplying it by any
unit-feature cross Gram kernel preserves \(\gamma\).  This follows by
writing the Gram kernel as an average of left/right diagonal unitary
dressings and applying trace-norm convexity.

The last two singleton Hadamard links have coefficient \(1/\sqrt N\) for
every physical placement: either one link is internal and has magnitude
\(1/\sqrt N\), or the alternating two-link wedge has rank \(N\) and entry
magnitude \(1/N\).  The arbitrary first link, including any internally
split decorated endpoint support, is exactly a Gram dressing of this tail.
Therefore

$$
(a,1,1,1),\ (1,1,1,a)
\quad\Longrightarrow\quad
\text{weighted coefficient}\le {1\over\sqrt N}
$$

for every support degree \(a\), every occurrence split, and every
correlated diagonal probe weight.  The cubic whole-block orbit attains the
bound, so the dimension power is sharp.

The unresolved higher-sector core is now a two-sided middle decoration,
starting with \((1,3,1,1)\) and its reversal: that support participates in
two planted links and cannot be absorbed into a one-sided Gram dressing.

Reproduction: notes/gram_dressed_tail_contraction.md and
tests/gram_dressed_tail_contraction.py.

## 18. Two-sided cubic middle benchmark

The \((1,3,1,1)\) middle decoration has now been screened at the actual
occurrence-split Schur level.  Singleton neighbors force the cubic support
to be an L-shape with missing corner \(z\), and its adjacent moments reduce
exactly to

$$
{1\over(q-1)^2}H_N(a,z)H_N(z,c).
$$

For the balanced \(1|2\) middle split, incidence diagonalization gives
normalized nuclear norm

$$
{q+4\over q^3\sqrt{q^2-1}},
$$

which is \(0.032275\) at \(q=4\), about \(3.4\times10^{-5}\) at
\(q=32\), and asymptotic to \(N^{-3/2}\).  Direct \(q=4\)
diagonalization of all 32 occurrence placements finds this to be the
largest uniform-orbit value.

For a fixed occurrence split, a general-diagonal theorem follows already
from rank--Frobenius.  Every nonzero entry has magnitude
\(1/[q^3(q-1)^2]\), and the smaller marked side has dimension at most
\(N^3=q^6\).  Therefore its weighted coefficient is at most

$$
{1\over(q-1)^2}\le {1\over q},
$$

equal to \(1/961\) at \(q=32\).  This is stable under injective ancillary
repetitions.  The remaining task is joint Bessel packing across overlapping
occurrence assignments and common base fibers, not another local
L-shape norm estimate.

Reproduction: notes/middle_cubic_schur_lift.md and
searches/signed_permutation_middle_cubic_schur_lift.py.

## 19. Exact dose-six sector frontier

Every hypothesis-sensitive block degree is odd.  Writing
\(n_b=1+2k_b\), the dose-six constraint is \(\sum_bk_b\le4\), giving
exactly 70 Fourier degree profiles through total degree twelve.  The
current local contractions close or explicitly bound nineteen:

- the minimal profile;
- all eight profiles with one decorated endpoint of degree \(3,5,7,9\);
  and
- the two cubic and two quintic middle profiles; and
- three double-cubic profiles with entry-rank bounds; and
- the other three double-cubic profiles with occurrence-slice bounds.

In particular, every total-degree-four or total-degree-six profile is now
locally controlled at coefficient at most \(1/q\) for a fixed occurrence
split, and every degree-eight profile now has an explicit local bound.
There are 51 open
profiles overall and 130 compatible odd record triples before symmetry
reduction.

This is a routing result, not a transcript budget.  The next pass must
combine the degree-eight constants with the shared occupation/Bessel mass,
improve the coarse record-three constant in the adjacent-middle case, and
then move to degrees ten and twelve; multiplying the 70 profile count by
\(1/q\) would undo the finite-size gain.

Reproduction: notes/dose_six_sector_frontier.md and
searches/dose_six_sector_frontier.py.

## 20. Quintic middle decoration

The two single-decoration degree-eight profiles are now locally closed.
A five-edge middle support compatible with singleton neighbors has one odd
row and one odd column.  Conditioning on either odd match leaves at least
one and at most two active even labels; their exact injective Walsh average
has magnitude at most \(1/(q-1)\).  Hence

$$
|M_{1,5}|,|M_{5,1}|\le {1\over q(q-1)}
$$

on the compatible sector, and the full four-block entry is at most
\(1/[q^3(q-1)^2]\).

There are eight marks, so a fixed occurrence-split flattening has rank at
most \(N^4=q^8\).  Weighted rank--Frobenius gives coefficient

$$
{q\over(q-1)^2}.
$$

At \(q=32\), this is \(32/961=0.033299\), only \(1.066\) times the
minimal \(1/q\) scale.  The only locally open degree-eight profiles now
have cubic decorations in two blocks.

Reproduction: notes/quintic_middle_contraction.md and
searches/signed_permutation_quintic_middle_bound.py.

## 21. Three double-cubic profiles

Three of the six double-cubic degree-eight profiles now have explicit local
bounds.  The separated profiles \((3,1,3,1)\) and \((1,3,1,3)\) contain
a cubic L-shape between singleton neighbors, so their entry magnitude is
at most \(1/[q^3(q-1)^2]\) and their fixed-split coefficient is at most
\(q/(q-1)^2\).

For \((1,3,3,1)\), the central record-one sector has the same bound.  Its
record-three sector obeys the coarser estimate

$$
{q^2\over\binom q3},
$$

equal to \(0.206452\) at \(q=32\).  This sector is locally finite but its
constant still needs the pure record-three compound contraction before a
joint dose-six sum.

The remaining hard degree-eight profiles are exactly
\((3,3,1,1)\), its reversal, and \((3,1,1,3)\).  Their record-one sectors
require a Bessel/compound argument rather than entrywise rank.

Reproduction: notes/double_cubic_entry_contractions.md and
searches/signed_permutation_double_cubic_entries.py.

## 22. Bessel-refined double-cubic path

The alternating block-coherent path estimate can use outer row/column
energies instead of unweighted operator norms:

$$
\|K^{13\mid24}\|_{1,p,q}
\le\sqrt{R_1C_3}\,\max|M_2|\sqrt{PQ}.
$$

For two cubic endpoints, both Bessel energies are one, so the central
Hadamard coherence gives the sharp \(1/q\) coefficient despite the
growing isolated cubic link norms.

For the record-one \((3,3,1,1)\) sector, restricting the middle cubic
support to an L-shape gives the exact row energy

$$
1+{2q\over(q-1)^2}.
$$

Together with the next-link coherence \(1/[q(q-1)]\), the alternating
coefficient at \(q=32\) is only about \(0.001041\).

This closes the whole-block flattenings of the three hard degree-eight
profiles.  The remaining degree-eight issue is specifically the
occurrence-split lift: prove that distinct-label mask/Bessel packing
preserves these outer energies when a cubic support is divided between
ket and bra.

Reproduction: notes/bessel_refined_path_contraction.md,
tests/three_link_weighted_path_contraction.py, and
searches/signed_permutation_double_cubic_entries.py.

## 23. Double-endpoint occurrence split

The hard \((3,1,1,3)\) occurrence lift is now bounded without factorial
collisions.  The exact cubic endpoint slice energies are

$$
E_0={q^2+2\over6},\quad
E_1={q^2+2\over2q^2},\quad
E_2={q^2-2q+2\over q^2(q-1)},\quad
E_3=q^{-2}.
$$

For endpoint splits \((k,\ell)\), the full kernel has row and column
energies at most \(E_kE_\ell\) and
\(E_{3-k}E_{3-\ell}\).  The better of the two Schatten factorizations is
worst at the \(1|2\) pairing and gives

$$
\gamma_{3,1,1,3}(q)=
\sqrt{{(q^2+2)(q^2-2q+2)\over2q^4(q-1)}}.
$$

At \(q=32\), this is about \(0.1232\), and degree-eight attenuation by
\((5/6)^8\) reduces it to about \(0.0287\).  This is valid for arbitrary
diagonal weights in each fixed unordered occurrence split.  Only the joint
packing of all splits against one passive probe law remains.

Reproduction: notes/double_endpoint_occurrence_contraction.md and
searches/double_endpoint_slice_energies.py.

## 24. Adjacent double-cubic occurrence split

The record-one \((3,3,1,1)\) sector and its reversal are now locally
controlled for every fixed unordered occurrence split.  Endpoint
record-one cubic supports and middle L-shapes have maximum incidence
degrees

$$
\begin{aligned}
D^E&=\left(
q\binom q3+q^2(q-1)\binom q2,
\binom{q-1}2+(q-1)\binom q2+q(q-1)^2,
q^2-2,1\right),\\
D^L&=\left(q^2(q-1)^2,3(q-1)^2,2(q-1),1\right).
\end{aligned}
$$

The largest squared endpoint-to-L entry is

$$
a_1={ (q+2)^2\over q^2(q-1)^2(q-2)^2}.
$$

Combining the incidence slice with the two row/column Schatten
factorizations gives

$$
\Gamma_{k\ell}(q)
\le\min\left\{
{\sqrt{a_1D^E_kD^L_\ell}\over q-1},
{\sqrt{a_1D^E_{3-k}D^L_{3-\ell}}\over q(q-1)}
\right\}.
$$

At \(q=32\), the largest of all sixteen bounds is \(0.009277\), attained
by the \((2,2)\) placement; degree-eight attenuation reduces it to
\(0.002158\).  Exact exceptional-pair counting improves this particular
placement to \(0.008715\).  Complete enumeration at \(q=4,8\) confirms
the formulas and the identity of the worst split, while the stated
\(q=32\) result uses the uniform incidence bound.

Thus all formerly hard record-one degree-eight profiles have small local
fixed-split coefficients.  The next bottleneck is joint occupation-mask
packing, followed by the record-three central sector and the adaptive
lift.

Reproduction: notes/adjacent_double_cubic_occurrence_contraction.md,
searches/adjacent_double_cubic_slice_energies.py, and
searches/adjacent_double_cubic_orbit_scaling.py.

## 25. Joint unordered occurrence packing

For a fixed realization of all query-node/frame-copy supports, let
\(w_b(x)\) be the number of containers in block \(b\) containing physical
coordinate \(x\), and set \(M_b=\sum_xw_b(x)\).  The joint assignment
count for an unordered degree-\(a_b\) Fourier support is

$$
e_{a_b}(w_b)=[z^{a_b}]\prod_x(1+w_b(x)z)
\le\binom{M_b}{a_b}.
$$

Across four blocks and the two complete-frame copies, \(\sum_bM_b\le12\)
at hard dose six.  Therefore every odd profile \(a\) has complete joint
occurrence square mass

$$
\mathcal S_a\le
\max_{\sum_bM_b\le12}\prod_b\binom{M_b}{a_b}.
$$

Exact enumeration of all seventy profiles gives maximum square masses by
total degree

$$
81,\quad160,\quad126,\quad36,\quad1
$$

for degrees \(4,6,8,10,12\), respectively.  A double-cubic degree-eight
profile has the sharper value \(100\).  This proves there is no separate
exponential split factor and no labeled-mark factorial; it also shows
that the highest degrees are strongly occupation-sparse.

This is a square-mass packing theorem, not yet a terminal trace-norm
bound.  The known two-node counterexample prevents a universal
square-function promotion.  The next contraction must combine these
profile masses with the exact cut-dependent sector weights.  Replacing
all cuts by their worst fixed-split coefficient is already quantitatively
too coarse for \((3,1,1,3)\).

Indeed, the deterministic occupation state \((2,1,1,2)\) makes the sum of
the current sixty-four double-endpoint split bounds equal to \(2.45633\),
or \(0.57126\) after degree-eight attenuation.  The attenuated minimal
sector and promise conditioning leave only \(0.16036\) below the TV
threshold.  Thus a successful two-endpoint contraction must improve this
splitwise sum by a factor below \(0.28071\), even before reserving room for
other sectors.  The analogous adjacent-double-cubic number is only
\(0.02512\), so it is no longer the main finite-size obstacle.

Reproduction: notes/joint_occurrence_profile_packing.md and
tests/joint_occurrence_profile_packing.py, with the quantitative barrier in
searches/fixed_split_occupation_barrier.py.

## 26. Joint double-endpoint Schur benchmark

The first balanced joint slice of \((3,1,1,3)\) has a hidden tight frame.
For a singled-out endpoint coordinate \(i\), complementary pair \(E\),
and neighboring singleton \(b\), define

$$
A(i;E,b)=M_{3,1}(\{i\}\cup E,b).
$$

Walsh orthogonality and the exact endpoint slice count give

$$
AA^*={q^2+2\over2}I_N.
$$

When both middle singleton coordinates lie on the pair side, the full
two-endpoint Schur matrix is also tight.  Its uniform normalized nuclear
coefficient is exactly

$$
{q^2+2\over q^3(q^2-1)},
$$

equal to \(0.01875\) at \(q=4\) and \(3.06071\times10^{-5}\) at \(q=32\).
For the alternating middle placement, exact row-Gram diagonalization gives
\(0.471592\) at \(q=2\) and \(0.0642009\) at \(q=4\).  No \(q=32\)
extrapolation is used for this second slice.

These are uniform orbit-weight benchmarks, not arbitrary-diagonal
theorems, but they locate the missing gain: the splitwise row-energy
method discards a strong compound tight-frame cancellation.  The next
target is a weighted version for the alternating matrix that preserves
the shared passive occupation law and gains at least the required factor
\(0.280708\).

Reproduction: notes/double_endpoint_joint_schur_benchmark.md and
searches/double_endpoint_joint_schur_benchmark.py.
