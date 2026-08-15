# Sharpness of the limiting level-twelve contraction

Date: 2026-07-15

Status: rigorous limitation of the current grouped-entry norm route.  The
unique sensitive level-twelve terminal type has a legal physical placement
for which the accepted all-projective upper bound is \(O(N^{-1})\), while an
explicit unit-vector contraction of the full same-layer-distinct tensor is

$$
N^{-1}(1-N^{-1})^2(1-2N^{-1})^3.
\tag{0.1}
$$

Thus the \(N^{-1}\) power is sharp for this graph-norm architecture.  This
does not prove that the true passive lower bound stops at \(N^{1/12}\), and
it is not a passive protocol.  It proves that a stronger exponent cannot
come from improving the same grouped injective norm while treating its frame
parties as arbitrary unit vectors.  Further progress requires additional
physical frame restrictions, cancellation before terminal norms, or a
different hard instance.

## 1. The terminal tensor and legal placement

The unique level-twelve type is the three-path forest

$$
G(a,b,c,d)
=\prod_{j=0}^2
 H_{a_jb_j}H_{b_jc_j}H_{c_jd_j},
\tag{1.1}
$$

with \(a_0,a_1,a_2\) distinct and likewise for the \(b,c,d\) layers.  Here
\(H\) is the normalized Sylvester matrix.  The physical placement groups
\((a_0,b_0)\) into one amplitude entry and leaves the other ten vertices in
singleton entries.  This is the accepted legal placement with assigned
suppression one.

The all-projective upper proof is already exact at the exponent level.  The
paired path has grouped injective norm at most one, and each singleton path
has norm at most \(N^{-1/2}\).  Hence

$$
\|G\|_{\varepsilon,\mathcal E}
\le C_GN^{-1},
\tag{1.2}
$$

where the same-layer distinctness mask contributes only a diagram constant
in that upper proof.

## 2. Explicit grouped unit vectors

Choose distinct coordinates \(r_0,r_1,r_2\) and distinct coordinates
\(s_1,s_2\).  For the paired path set

$$
U_0(a,b)=H_{ab}H_{br_0},\qquad
C_0=e_{r_0},\qquad
D_0=He_{r_0}.
\tag{2.1}
$$

For singleton path \(j=1,2\), set

$$
A_j=He_{s_j},\qquad
B_j=e_{s_j},\qquad
C_j=e_{r_j},\qquad
D_j=\operatorname{sgn}(H_{s_jr_j})He_{r_j}.
\tag{2.2}
$$

Every displayed vector has Euclidean norm one.  For the grouped vector,

$$
\sum_{a,b}|U_0(a,b)|^2
=\sum_{a,b}|H_{ab}|^2|H_{br_0}|^2
=1.
\tag{2.3}
$$

Without the distinctness masks, the paired path contracts to one and each
singleton path contracts to \(N^{-1/2}\).  The combined value is therefore
exactly \(N^{-1}\).

## 3. Retaining every same-layer distinctness mask

The chosen \(c\)-coordinates are already distinct.  In every remaining
layer the vector signs make all surviving summands nonnegative, so the mask
cost can be counted exactly.

- In the \(a\) layer, the fraction of distinct triples is
  \((N)_3/N^3=(1-N^{-1})(1-2N^{-1})\).
- In the \(b\) layer, \(b_1=s_1\) and \(b_2=s_2\), while the dense \(b_0\)
  coordinate avoids those two values.  The factor is \(1-2N^{-1}\).
- The \(c\)-layer factor is one.
- The \(d\)-layer factor is again
  \((1-N^{-1})(1-2N^{-1})\).

Multiplying by the unmasked value gives

$$
\begin{aligned}
\langle G_{\rm distinct},
 U_0,C_0,D_0,A_1,B_1,C_1,D_1,A_2,B_2,C_2,D_2\rangle
&=N^{-1}(1-N^{-1})^2(1-2N^{-1})^3\\
&={ (N-1)^2(N-2)^3\over N^6}.
\end{aligned}
\tag{3.1}
$$

This is a lower bound on the grouped injective norm.  Its ratio to \(N^{-1}\)
tends to one, so (1.2) and (3.1) match in their dimension exponent.

At \(N=1024\), the exact value and ratio are

$$
0.000968956353997,\qquad
{0.000968956353997\over1024^{-1}}
=0.992211306493.
\tag{3.2}
$$

The distinct-coordinate mask therefore removes less than one percent from
this witness at the realistic-size benchmark.  It cannot supply either a
better \(N\)-power or a useful constant-scale rescue of this proof route.

## 4. The coefficient family does not internally cancel

The terminal graph is produced by six fresh transfers.  The initial triple
displayed in the earlier witness has 180 legal orders producing the same
labeled support.  Complete enumeration of the canonical type gives 12
initial triples and 1,080 legal histories.  Every transfer in every history
is fresh, so none differentiates an existing local weight; all local factors
are positive Stein kernels or positive \(\psi'\) factors.  The terminal
support has three first-layer marks and survives the task antisymmetrization.
Consequently all 1,080 histories have the same scalar sign.

This eliminates two simple escape routes:

- same-layer distinctness does not improve the graph exponent; and
- the explicit level-twelve coefficient does not vanish by cancellation of
  its legal transfer histories.

The result remains scoped.  A proof that uses restrictions on the physical
frame factors, rather than arbitrary grouped unit vectors, could still gain.
A different hard-instance mixture could also cancel or suppress this family
while preserving the active signal.

## 5. Program decision

The current asymptotic theorem remains

$$
D_{\mathsf P}^{\rm hard}=\Omega(N^{1/12}),
\tag{5.1}
$$

but the immediate proposal to improve it by sharpening the level-twelve
grouped graph norm is closed.  Round 3 should now prioritize one of:

1. an explicit physical-frame restriction that excludes the optimizer in
   Section 2 and is stable under adaptive posterior selection;
2. an explicit-constant finite-size evaluation of the accepted transcript
   bound at \(N=1024,D=6\); or
3. a common-scorecard alternative plant whose active signal survives while
   the all-fresh three-path coefficient is cancelled or attenuated.

The realistic-size track should remain the lead: the present lower witness
shows that even the distinctness mask is almost saturated at \(N=1024\), so
another graph-norm refinement is less promising than a physical or
hard-instance change.

Reproduction:

- `searches/level_twelve_contraction_sharpness.py` constructs the eleven
  grouped unit vectors, directly contracts all distinctness masks at small
  powers of two, counts all 12 initial triples and 1,080 positive histories,
  and evaluates (3.1).
- `tests/level_twelve_contraction_sharpness.py` protects the exact formula,
  unit norms, matching upper/lower exponents, positivity, history count, and
  the \(N=1024\) ratio.
