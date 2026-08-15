# All-projective repair of the terminal three-path witness

Date: 2026-07-15

Status: rigorous positive result for the second bounded OP3-7 test.  The
true terminal forest from `terminal_interpolation_sigma_one_witness.md`
does not saturate the best safe physical contraction.  Although the global
assigned dichotomy retains only \(N^{-1/2}\), a separate all-projective
contraction gives \(N^{-1}\).  This is stronger than the \(N^{-3/4}\)
level-twelve target and avoids the invalid RT-003 Hilbert-to-projective
conversion.

This result removes the explicit three-path forest as an obstruction to
\(N^{1/16}\).  It does not yet prove \(N^{1/16}\) for the full interpolation
image, and its diagram constants do not give passive dose \(>6\) at
\(N=1024\).

Follow-up: the complete enumeration in
`high_level_terminal_best_of_two_audit.md` uses this lemma to prove the new
\(\Omega(N^{1/20})\) theorem and isolates three different level-nine/ten
joint saturators.  `level_ten_forest_mean_zero_repair.md` subsequently
removes the level-ten type and improves the current theorem to
\(\Omega(N^{1/18})\).

## 1. The safe norm choice

Fix a physical placement \(\phi\) and lift every marked physical support to
its ordered coordinate fibers.  Let \(\mathcal E\) be the set of physical
ket and bra amplitude entries.  The reverse complete-frame argument has two
logically separate parts:

1. after division by the joint insertion/Bessel mass
   \(\mathcal B_\phi\), the frame skeleton has projective mass at most a
   diagram constant over the grouped entry spaces
   \(\bigotimes_{E\in\mathcal E}\mathcal H_E\); and
2. the graph tensor must be bounded in the dual grouped-entry injective
   norm.

The first statement does not require every graph component to be singleton
in every entry.  A multi-mark entry is one joint Hilbert coordinate.  At a
marked node, complete-frame Cauchy--Schwarz bounds the product of the ket
and bra Hilbert norms by the same joint Bessel mass; reverse stochastic
summation and projective subadditivity preserve the bound.  Ordered-support
collisions cost only a diagram factorial.

The singleton hypothesis in Case I of the repaired contraction is used to
obtain a universal graph-norm estimate, not to construct the projective
frame decomposition.  Therefore any fixed placement with a directly small
grouped-entry graph norm may use this all-projective route, even if one
component has two marks in an entry.

This is crucial: the proof stays projective from the terminal leaves to the
final graph pairing.  It never first produces a Hilbert auxiliary and then
asks for its projective norm.

## 2. Exact flattening norm of one path

For one four-layer path write

$$
T(a,b,c,d)=H_{ab}H_{bc}H_{cd},
\tag{2.1}
$$

where \(H\) is the normalized Sylvester matrix.  For a vertex cut
\(S\mid S^c\), let \(r(S)\) be the binary rank of the path adjacency across
the cut.  Internal edges contribute only row or column diagonal unitaries.
The bicharacter matrix across the cut has exact operator norm

$$
\left\|T_{S\mid S^c}\right\|_{\rm op}
=N^{(4-3-r(S))/2}.
\tag{2.2}
$$

There are three cases relevant here.

### Every vertex is in a distinct entry

The cut \(S=\{a,c\}\), \(S^c=\{b,d\}\) respects the physical entries and
has rank two.  Hence

$$
\|T\|_{\varepsilon,\mathcal E}
\le N^{-1/2}.
\tag{2.3}
$$

Indeed the corresponding \(N^2\times N^2\) flattening \(M\) satisfies

$$
MM^*={1\over N}I_{N^2}.
\tag{2.4}
$$

### The largest entry contains two vertices

The path still occupies at least two physical entries.  Choose any
nonempty proper union of its entry blocks.  Because the path is connected,
at least one edge crosses, so \(r(S)\ge1\).  Equation (2.2) gives

$$
\|T\|_{\varepsilon,\mathcal E}\le1.
\tag{2.5}
$$

For the explicit adjacent pair \((a,b)\), this bound can also be seen
directly.  For unit \(u_{ab},v_c,w_d\), successive contraction by \(H\)
gives vectors of norm at most one, and therefore

$$
\left|
\sum_{a,b,c,d}H_{ab}H_{bc}H_{cd}u_{ab}v_cw_d
\right|\le1.
\tag{2.6}
$$

### All four vertices are in one entry

There is no nontrivial entry cut.  The grouped tensor is one Hilbert vector
with norm \(N^{(4-3)/2}=N^{1/2}\).  This case is not present in the terminal
witness, but it explains why the all-projective route is not a universal
replacement for the assigned contraction.

## 3. Apply the bound to the true terminal forest

Let \(P_0,P_1,P_2\) be the three disjoint paths.  For every placement with
maximum component occupancies

$$
(k_0,k_1,k_2)=(2,1,1),
\tag{3.1}
$$

the component bounds are

$$
\|T_{P_0}\|_{\varepsilon,\mathcal E}\le1,
\qquad
\|T_{P_1}\|_{\varepsilon,\mathcal E},
\|T_{P_2}\|_{\varepsilon,\mathcal E}\le N^{-1/2}.
\tag{3.2}
$$

Termwise vertical injective-norm multiplicativity gives

$$
\|T_{P_0}\boxtimes T_{P_1}\boxtimes T_{P_2}\|_
{\varepsilon,\mathcal E}
\le N^{-1}.
\tag{3.3}
$$

The same-layer distinct-coordinate mask couples the three components, but
the accepted character expansion resolves it.  There are three vertices in
each of four layers, hence

$$
P_G=4\binom32=12.
\tag{3.4}
$$

The expansion has coefficient mass at most \(2^{12}\).  Each term applies
only local diagonal character unitaries, including a joint diagonal unitary
on the paired entry, so (3.3) is unchanged termwise.  Thus the normalized
fixed-skeleton contraction obeys

$$
|\mathfrak C_{G,\phi}|
\le C_G\,N^{-1}\mathcal B_\phi,
\tag{3.5}
$$

where \(C_G\) is independent of \(N\), depth, outcome width, and dose.

## 4. Comparison with the required gain

The global assigned dichotomy uses

$$
\sigma_{\rm ass}
=(2-1)+(1-1)+(1-1)=1
\tag{4.1}
$$

and therefore returns only \(N^{-1/2}\).  Equation (3.5) supplies an extra
\(N^{-1/2}\), while the level-twelve \(N^{1/16}\) row requires only an
extra \(N^{-1/4}\):

$$
N^{-1}
\le N^{-3/4}
\le N^{-1/2}.
\tag{4.2}
$$

At \(N=1024\), the three graph powers are

$$
N^{-1}=0.0009765625,
\qquad
N^{-3/4}=0.0055242717,
\qquad
N^{-1/2}=0.03125.
\tag{4.3}
$$

These numbers compare only the dimension powers.  The safe mask and path
constants have not been optimized, so (4.3) is not a realistic-size
passive-dose certificate.

## 5. Why RT-003 does not apply

RT-003 starts with an arbitrary Hilbert auxiliary of norm one and then
tries to bound its projective mass.  The maximally entangled Hadamard vector
shows a factor \(\sqrt N\) can be lost.

No such conversion occurs here.  The frame skeleton is expanded in the
grouped-entry projective norm from the start, including the paired entry as
one Hilbert party.  The graph is paired with that decomposition only after
the projective mass is controlled.  The proof is therefore the same safe
norm regime as the accepted all-singleton case, with a placement-specific
graph estimate replacing the universal singleton estimate.

## 6. Track B decision

The explicit terminal three-path family passes the requested gate with more
decay than required.  It is not a counterexample to the \(N^{1/16}\)
mechanism.

The subsequent obstruction search is complete.  It enumerates the true
terminal image at levels nine through twelve and evaluates, for every
physical placement, the better of:

1. the global assigned exponent; and
2. the safe grouped-entry all-projective exponent.

Exactly three types keep both routes at only \(N^{-1/2}\): two reflected
level-nine trees and one level-ten \(6+4\) forest.  This classification
proves \(\Omega(N^{1/20})\).  The subsequent centered-weight audit removes
the level-ten type and proves \(\Omega(N^{1/18})\); the reflected level-nine
pair is now the next signed and physical contraction target.

Reproduction:

- `searches/terminal_three_path_projective_repair.py` enumerates all fifteen
  physical-entry partitions of one path and computes the exact safe cut
  exponents.
- `tests/terminal_three_path_projective_repair.py` checks every path cut
  directly at \(N=2,4\), protects the occupancy classes, and verifies the
  combined \(N^{-1}\) bound and \(N^{-1/4}\) target slack.
