# Separated balanced endpoint-slice contraction

Date: 2026-07-15

Status: rigorous arbitrary-diagonal bound for the third balanced orbit.  The
bound keeps the realistic-size diagnostic below \(1/3\), but only by
\(0.0006682144\).  This is evidence that exact chain structure works; it is
also a warning that the present completion route is fragile.

## 1. Target and exact cut

Consider profile \((3,1,1,5)\) and split \((1,0,0,4)\).  Write the cubic
endpoint support as

$$
\{x\}\mathbin{\dot\cup}E,\qquad |E|=2,
$$

and the quintic endpoint support as

$$
F\mathbin{\dot\cup}\{e\},\qquad |F|=4.
$$

With rows \((x,F)\) and columns \((E,b,c,e)\), the exact occurrence tensor is

$$
K_{(x,F),(E,b,c,e)}
=M_{3,1}(\{x\}\cup E,b)H_N(b,c)M_{1,5}(c,F\cup\{e\}),
\tag{1.1}
$$

on \(x\notin E\) and \(e\notin F\).  Here \(N=q^2\) and \(H_N\) is the
normalized Sylvester matrix, so every entry has modulus \(N^{-1/2}=q^{-1}\).
Complement and path reversal give three additional cuts with the same bound.

## 2. The generic relaxation is too expensive

Treating the \(1|2\) and \(4|1\) distinct-label masks independently gives

$$
{\gamma_{N,2}\gamma_{N,4}\over q}
=0.225293047398
\tag{2.1}
$$

at \(N=1024\).  This is safe but unusable in the current Perron ledger.  The
proof therefore has to retain the exact endpoint slices instead of replacing
both by generic disjointness multipliers.

## 3. Cubic fixed-pair slice

For every fixed pair \(E\) and singleton \(b\), exact signed-permutation
moments give

$$
\sum_{x\notin E}
|M_{3,1}(\{x\}\cup E,b)|^2
\le T_2,
\qquad
T_2={q^2-2q+2\over q^2(q-1)}.
\tag{3.1}
$$

Regard the values indexed by \(x\) as a feature vector for the column
\((E,b)\).  The standard feature factorization of a Schur multiplier then
gives norm at most \(\sqrt{T_2}\).  At \(q=32\),

$$
T_2=0.0303049395161,\qquad
\sqrt{T_2}=0.174083139666.
\tag{3.2}
$$

## 4. Collapse the central link into the quintic slice

The endpoint moment has the exact form

$$
M_{1,5}(c,Q)=v_5(Q)H_N(c,\operatorname{xor}Q).
\tag{4.1}
$$

Consequently,

$$
H_N(b,c)M_{1,5}(c,Q)
={v_5(Q)\over q}
 H_N(c,b\oplus\operatorname{xor}Q).
\tag{4.2}
$$

The residual Hadamard symbol in (4.2) is a unit-feature Schur multiplier:
its dependence on \(F\) is one normalized Hadamard row, while its dependence
on \((b,c,e)\) is a phase times the basis vector indexed by \(c\).

For fixed \(F\), the exact fixed-four quintic moment slice is

$$
\sum_{e\notin F}|M_{1,5}(c,F\cup\{e\})|^2
\le T_4,
\qquad T_4=1-{4\over N}.
\tag{4.3}
$$

Because \(|H_N|^2=1/N\), extracting \(v_5\) from (4.1) multiplies this
squared energy by \(N\):

$$
\sum_{e\notin F}|v_5(F\cup\{e\})|^2\le NT_4.
\tag{4.4}
$$

Thus the scalar \((F,e)\) multiplier has norm at most \(\sqrt{NT_4}\).
The explicit \(1/q=1/\sqrt N\) in (4.2) cancels the \(\sqrt N\), leaving
only \(\sqrt{T_4}\).

This normalization is essential.  Dividing \(\sqrt{T_2T_4}\) by another
factor of \(q\) would be wrong: the apparently missing gain has already been
spent when the normalized Hadamard factor is removed from \(M_{1,5}\).
Complete \(q=4\) enumeration independently gives quintic scalar energy
\(12=N(1-4/N)\), protecting this point.

## 5. The theorem

Multiplying the cubic, scalar-quintic, and unit residual factors proves

$$
\boxed{
\|K\|_{\mathrm{Schur}}
\le \sqrt{T_2T_4}
=\sqrt{
{q^2-2q+2\over q^2(q-1)}
\left(1-{4\over q^2}\right)}
}.
\tag{5.1}
$$

The feature proof is stable under arbitrary nonnegative row and column
diagonal weights, so (5.1) is an arbitrary-physical-law coefficient, not a
uniform-weight diagnostic.  At \(N=1024\),

$$
\boxed{c_{(3,1,1,5):(1,0,0,4)}=0.173742800847.}
\tag{5.2}
$$

This improves the generic value \(0.225293047398\), but it is much larger
than the provisional \(1/32\).

## 6. Ledger consequence

Insert (5.2) on its four complement/reversal cuts together with the two
previous balanced theorems.  Reoptimizing the exact 210-state diagnostic and
the extended promise loss gives

$$
\beta=0.779328831259,\qquad
\operatorname{TV}_{\rm diagnostic}=0.332665118954.
\tag{6.1}
$$

Hence

$$
\boxed{{1\over3}-\operatorname{TV}_{\rm diagnostic}
=0.000668214379.}
\tag{6.2}
$$

Twelve of the 888 balanced entries are now controlled by chain-aware
arbitrary-law theorems; 876 still carry provisional charges.  The route
survives, but the remaining margin is only about \(8.6\%\) of the margin
after the second theorem.  No remaining provisional coefficient may be
treated as harmless.

After reoptimization and removal of proved and explicitly forced orbits, the
largest unresolved provisional orbit is

$$
(3,1,1,5):(1,1,1,2),
$$

whose current Perron contribution is \(0.00150951797\).  This already exceeds
the total remaining slack.  Reoptimizing while varying only this four-cut
orbit gives the numerical acceptance gate

$$
\boxed{c_{(3,1,1,5):(1,1,1,2)}<0.04504054678
=1.441297497/q.}
\tag{6.3}
$$

Its rigorous coefficient is therefore the next go/no-go test.  The following
orbit is \((3,1,1,5):(0,0,1,4)\), with contribution \(0.00128274051\).

## 7. Reproduction and scope

- `searches/separated_balanced_endpoint_slice_contraction.py` evaluates the
  exact slices, safe and generic coefficients, reoptimized ledger, reranked
  frontier, and next-orbit acceptance gate.
- `tests/separated_balanced_endpoint_slice_contraction.py` independently
  enumerates the \(q=4\) cubic and extracted-quintic slices, checks the
  normalization \(NT_4\), and stress-tests the abstract factorization under
  arbitrary diagonal laws.
- `./run_round3_checks.sh` runs the complete inherited and round-three suite.

This theorem does not complete the high-sector ledger, certify the remaining
physical diagnostic charges, or lift the one-batch result through adaptivity.
It establishes one more local contraction and sharply changes the route
decision: continue for one reranked orbit at a time, with an immediate pivot
if its coefficient exceeds \(0.0450405468\).

## 8. Subsequent fourth-orbit result

The next-orbit gate in (6.3) has now passed.  A shared-law contraction for
\((3,1,1,5):(1,1,1,2)\) gives coefficient \(0.0250967461\), even though
separate endpoint slice maxima would give \(0.0991470258\).  With this fourth
theorem inserted, the diagnostic is \(0.3323683002\), leaving
\(0.0009650332\); sixteen entries are controlled and 872 remain provisional.
The next reranked orbit is \((3,1,1,5):(0,0,1,4)\), with gate
\(0.0542506298\).  See `internal_singleton_shared_law_contraction.md` for the
proof and stress tests.

## 9. Subsequent fifth-orbit result

That next gate has also passed.  Exact quintic row energy plus
rank--Frobenius gives coefficient \(0.0311889051<1/32\) for
\((3,1,1,5):(0,0,1,4)\).  The five-theorem diagnostic is
\(0.3323657941\), leaving \(0.0009675392\).  The next reranked orbit is
\((1,1,3,5):(0,1,1,3)\), with gate \(0.0570749885\).  See
`column_cubic_quintic_row_contraction.md`.
