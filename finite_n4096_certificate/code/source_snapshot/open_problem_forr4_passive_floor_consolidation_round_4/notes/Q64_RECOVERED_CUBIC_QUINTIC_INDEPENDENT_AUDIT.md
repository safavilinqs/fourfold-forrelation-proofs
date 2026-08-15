VERDICT: REJECTED

# Independent audit of the recovered cubic--quintic repair

Date: 2026-07-17

This audit rejects the claim that all forty recovered cubic--quintic entries are proved. The twenty-eight endpoint-row entries do admit a complete arbitrary-correlated-diagonal proof after one unsafe cubic--cubic record-one use is replaced by the generic physical bound. The twelve entries assigned to `four_sector_physical_chain` do not: an exact endpoint-compatible $q=8$ physical chain has entry $1/17920$, strictly larger than the claimed maximum-entry bound $1/25088$.

At the time of this audit, the registry consequence was that 28 recovered entries remained proved and the 12 canonical-profile $(1,3,5,3)$ entries returned to the coefficient-one quarantine. The global actual-mask count at that checkpoint was therefore 342 of 354, not 354 of 354. Section 8 records the later replacement theorem and current registry state.

## 1. Physical occurrence matrices and arbitrary-law norm

Let $q=2^m$ and let $\Omega=\mathbb F_2^m\times\mathbb F_2^m$ be the $q^2$ cells of one Fourier block. A physical degree-$d$ support is a simple subset $U\subseteq\Omega$ with $|U|=d$; simplicity means that no cell occurs twice. For a uniformly random signed permutation $P$ and the unnormalised Walsh matrix $H$, define

$$
X=HP,
\qquad
Y=PH,
$$

and define the exact link moment

$$
M(S,T)=\mathbb E_P\left[\prod_{(a,c)\in S}X_{a,c}\prod_{(b,d)\in T}Y_{b,d}\right].
$$

For a four-block profile $d=(d_0,d_1,d_2,d_3)$ and a split $t=(t_0,t_1,t_2,t_3)$, a row index is a tuple $A=(A_0,A_1,A_2,A_3)$ with $|A_i|=t_i$, and a column index is a tuple $B=(B_0,B_1,B_2,B_3)$ with $|B_i|=d_i-t_i$. The physical occurrence matrix is

$$
\mathcal O_t(A,B)
=
\mathbf 1\{A_i\cap B_i=\varnothing\text{ for every }i\}
\prod_{i=0}^2 M(A_i\cup B_i,A_{i+1}\cup B_{i+1}).
$$

Thus every cross-cut distinctness mask is the literal condition $A_i\cap B_i=\varnothing$ inside one simple support $U_i=A_i\cup B_i$. There is no product-law assumption among the four blocks.

For arbitrary probability laws $p$ on all row tuples and $w$ on all column tuples, including laws correlated across blocks, the coefficient to be bounded is

$$
\left\|
\operatorname{diag}(\sqrt p)\,
\mathcal O_t\,
\operatorname{diag}(\sqrt w)
\right\|_1.
$$

No step below assumes uniformity, exchangeability, translation invariance, product structure across blocks, or separately chosen laws on different links.

### Link records

For a support $S$, let $R(S)$ be the set of row coordinates having odd multiplicity and let $C(S)$ be the set of column coordinates having odd multiplicity. Sign averaging forces a link $M(S,T)$ to vanish unless $|C(S)|=|R(T)|$. The common odd cardinality is the link record.

For the profiles in this audit, the record triple is

$$
\rho=(1,r,s),
\qquad
r,s\in\{1,3\}.
$$

The first record belongs to the singleton--first-higher-degree link, $r$ belongs to the next link, and $s$ belongs to the final link. In particular, `endpoint_energy(..., records[1])` correctly receives the outgoing record $r$ of the first higher-degree block; passing `records[0]` would incorrectly pass the forced singleton record.

### Generic link bound

After sign averaging, the hidden permutation must map the $r$ odd domain coordinates bijectively to the $r$ odd codomain coordinates. Both the odd and even permutation sums have entries of modulus one. Therefore

$$
|M(S,T)|
\le
\frac{r!(q-r)!}{q!}
=
\binom qr^{-1}.
$$

This bound is universal over all simple physical supports in the fixed record sector.

## 2. The endpoint-row mechanism

The endpoint mechanism applies to the seven complement/reversal orbits with canonical profile $(1,3,3,5)$ or $(1,5,3,3)$, for 28 entries total.

Fix $\rho=(1,r,s)$. For a row split, let $E_r(k,\epsilon)$ be the maximum, over the fixed singleton/first-higher partial support, of the complete squared row sum over the singleton complement and first-higher complement, restricted to incoming record one and outgoing record $r$:

$$
E_r(k,\epsilon)
=
\max_{A_0,A_1}
\sum_{B_0,B_1}
\mathbf 1\{A_j\cap B_j=\varnothing\}
\mathbf 1\{|R(U_1)|=1,\ |C(U_1)|=r\}
|M(U_0,U_1)|^2.
$$

Here $|A_1|=k$, $|A_0|=\epsilon$, and $U_j=A_j\cup B_j$. This energy already includes the singleton--first-higher link and no other link.

For blocks 2 and 3, let $D_i^\rho(t_i)$ be the maximum number of simple physical completions of a fixed $t_i$-cell partial support in the prescribed incoming and outgoing record sector. Since the two remaining link moments are bounded pointwise by $m_r$ and $m_s$, a complete row satisfies

$$
\sum_B|\mathcal O_t^\rho(A,B)|^2
\le
E_r(t_1,t_0)
\,D_2^\rho(t_2)D_3^\rho(t_3)
\,(m_rm_s)^2.
$$

The adjacent factors are ordinary multiplication; the displayed line is the precise justification for the code expression

$$
E_r(k,\epsilon)D_3^\rho(t_3)D_4^\rho(t_4)(m_rm_s)^2.
$$

The endpoint link is neither omitted nor double-counted: it occurs once inside $E_r$, while the product contains exactly links 1 and 2 of the remaining chain.

### Why a complete-row energy controls arbitrary laws

Let $L_\rho$ be the displayed complete-row bound. Factor the sector matrix as $\mathcal O_t^\rho=XY^*$ with $X_A=(\mathcal O_t^\rho(A,B))_B$ and $Y_B=e_B$. Then

$$
\left\|\operatorname{diag}(\sqrt p)X\right\|_F^2
\le L_\rho,
\qquad
\left\|\operatorname{diag}(\sqrt w)Y\right\|_F^2=1.
$$

The Schatten Hölder inequality gives the arbitrary-law bound

$$
\left\|\operatorname{diag}(\sqrt p)\mathcal O_t^\rho\operatorname{diag}(\sqrt w)\right\|_1
\le\sqrt{L_\rho}.
$$

Applying the same argument to the transpose gives the complete-column bound. Both bounds apply to the same fixed record sector and the same arbitrary laws, so their minimum is valid inside that sector. The four record sectors are disjoint at the level of full physical supports; even without using orthogonality, the trace-norm triangle inequality justifies summing their four independently bounded norms.

### Endpoint-energy scope

The imported quintic endpoint energies have the following exact scope at $q=64$: fixed-one energies are record-resolved; variable-singleton fixed-pair and fixed-singleton fixed-triple energies are maxima over every simple affine partial-support type; fixed-singleton values are divided by $q^2$ when obtained from a variable-singleton sum; and the variable-singleton fixed-four values are stated rational upper bounds. For partial supports of size at most three over $\mathbb F_2^m$, the equality type determines the affine type because two distinct nonzero differences are linearly independent. The enumeration therefore covers the physical shapes it claims to cover.

The cubic endpoint implementation is record-resolved only at the fixed-one slice. At the other slices, `cubic_endpoint_slice_energies` is the total over outgoing records, and the code safely uses that total as an upper bound for each individual record. The independent $q=4$ enumeration confirms this distinction: for a fixed cubic pair the direct record-one and record-three energies are $1/24$ and $1/6$, while the code charges the safe unrefined total $5/24$ to either record. Those slices must not be described as record-resolved.

Whenever the singleton lies on the completion side, there are exactly $q^2$ possible singleton cells. Multiplying a uniform fixed-singleton row bound by $q^2$ is therefore safe. It is a counting inequality; equality of the maximising partial shape under every translation is not required. The independent $q=4$ enumeration obtains the exact factor $16=q^2$ at the fixed-one cubic slice.

### The cubic--cubic record-one scope correction

The smaller cubic--cubic record-one estimate

$$
b_{11}=\frac{q+2}{q(q-1)(q-2)}
$$

is not universal. A vertical cubic followed by a horizontal cubic has exact moment $1/q>b_{11}$ for $q\ge8$. The estimate is valid when the left cubic also has exterior record one. In that case the left cubic must be an L shape and has one nonzero even-column xor. Conditioning on the forced odd match leaves $n=q-1$ coordinates. If the right cubic has no even row, the remaining character sum has magnitude $1/[q(q-1)]$. If it has one even row, separating whether the special domain maps to the special codomain gives a numerator bounded by

$$
(n-1)+4=q+2
$$

over $q(q-1)(q-2)$ possible normalisation, which proves $b_{11}$.

For profile $(1,3,3,5)$, the relevant left cubic is the endpoint cubic and has exterior incoming record one, so $b_{11}$ is in scope. For profile $(1,5,3,3)$ and sector $(1,3,1)$, the left cubic has incoming record three, so the old use of $b_{11}$ was invalid. The audited code now uses the universal $1/q$ bound in that sector. With that correction, all 28 endpoint entries remain below one; their seven orbit coefficients are

$$
0.0162888571820,
\quad 0.0264913195629,
\quad 0.0331908274616,
\quad 0.0506481680416,
\quad 0.110898445736,
\quad 0.413910743027,
\quad 0.703615181088.
$$

### Complement and reversal

Complementing a split swaps $A_i$ and $B_i$ inside the same full support, preserves every record and every distinctness mask, and transposes the occurrence matrix. Reversal sends

$$
(U_0,U_1,U_2,U_3)
\longmapsto
(U_3^\top,U_2^\top,U_1^\top,U_0^\top).
$$

The change of variables $P\mapsto P^\top$ gives $M(S,T)=M(T^\top,S^\top)$, so reversal preserves the physical kernel and reverses the record triple. Since the proof sums all $(r,s)$ sectors, this is a bijection of sectors rather than an omitted orientation assumption. The independent test evaluates nonzero physical $q=4$ chains and their transposed reversals for all ten claimed orbits, not merely the tuple transformation.

## 3. The four-sector chain mechanism is invalid

The twelve rejected entries have canonical profile $(1,3,5,3)$. The endpoint-cubic bounds themselves are correct:

$$
e_1=\frac1{q(q-1)},
\qquad
e_3=\frac1q.
$$

To prove them directly, condition the singleton's odd coordinate to map to the cubic's unique odd row. If the cubic has one even row, its nonzero xor is averaged over the remaining $q-1$ preimages and gives magnitude $1/[q(q-1)]$. If the cubic consists of three cells in its odd row, no even character remains and the magnitude is $1/q$. These cases cover every simple cubic support with incoming record one, and they respectively realise the record-one and record-three maxima.

The claimed cubic--quintic values

$$
m_1=\frac{q+2}{q(q-1)(q-2)},
\qquad
m_3=\frac3{(q-3)\binom q3}
$$

are not universal over the physical record sectors used by the chain.

For record one, take a cubic with three cells in one column and a quintic with cells $(0,0),(0,1),(0,2),(0,3),(1,0)$. The four cells in the even quintic row have xor zero, the link record is one, the other quintic record is three, and the exact moment is $1/q$. Thus $m_1$ already fails at $q=8$ in a geometry allowed at the terminal link.

For record three, take a cubic with cells $(0,0),(0,1),(0,2)$ and a quintic with cells $(0,0),(0,1),(0,2),(1,0),(2,1)$. The cubic has exterior record one, the link record is three, the other quintic record is one, and direct permutation summation gives

$$
|M|=\frac1{3\binom q3}.
$$

At $q=64$ this exceeds $m_3$ because $1/3>3/(q-3)$. Hence $m_3$ also lacks the claimed endpoint-compatible scope.

### Smallest explicit failure of the complete maximum-entry product

At $q=8$, use flattened cell coordinates and the four simple supports

$$
U_0=\{0\},
\qquad
U_1=\{0,1,2\},
\qquad
U_2=\{0,1,8,16,24\},
\qquad
U_3=\{0,1,2\}.
$$

Their block records are

$$
(|R(U_1)|,|C(U_1)|)=(1,3),
\quad
(|R(U_2)|,|C(U_2)|)=(3,1),
\quad
(|R(U_3)|,|C(U_3)|)=(1,3),
$$

so this is an endpoint-compatible sector $(1,3,1)$. Direct averaging over all $8!$ permutations, with the sign parity imposed exactly, gives

$$
M(U_0,U_1)=\frac18,
\qquad
M(U_1,U_2)=-\frac1{280},
\qquad
M(U_2,U_3)=-\frac18.
$$

The physical chain entry is therefore

$$
\left|\prod_{i=0}^2M(U_i,U_{i+1})\right|
=
\frac1{17920}.
$$

The claimed sector maximum is

$$
a_{3,1}=e_3m_3m_1
=
\frac18\cdot\frac3{280}\cdot\frac5{168}
=
\frac1{25088}.
$$

Hence

$$
\frac1{17920}>\frac1{25088},
$$

and the first invalid inequality is exactly

$$
|\mathcal O_t^\rho(A,B)|\le a_{r,s}=e_rm_rm_s.
$$

At $q=64$, the same support family violates the claimed product by the exact factor $217/11$. Because any full support can be partitioned into a row subset of any requested size and its disjoint complement, this same physical chain supplies an occurrence entry for every split in the three $(1,3,5,3)$ complement/reversal orbits. All twelve entries are therefore affected.

### What remains valid in the rank and incidence argument

Conditional on a valid maximum-entry bound $a$, the three norm estimates in the code are correct.

There are at most $(q^2)^{|t|}=q^{2|t|}$ row indices and at most $q^{2(12-|t|)}$ column indices. Thus

$$
\operatorname{rank}(\mathcal O_t^\rho)
\le q^{2\min(|t|,12-|t|)}.
$$

The diagonally weighted Frobenius norm is at most $a$, so the trace norm is at most

$$
q^{\min(|t|,12-|t|)}a.
$$

The exponent is correct even though a cell has $q^2$ values: the number of features is bounded by $q^{2|t|}$, and the trace-norm/Frobenius conversion takes its square root.

If every row has at most $R$ physical completions, factor the matrix by its column-coordinate basis. The first factor has weighted Frobenius norm at most $a\sqrt R$ and the second has weighted Frobenius norm one, giving the arbitrary-law row-incidence bound $a\sqrt R$. The transposed factorisation gives $a\sqrt C$. A maximum-entry estimate alone is not used as a nuclear-norm estimate.

`odd_record_incidence` counts a genuine one-axis relaxation: for each occupied coordinate group it chooses distinct cells from the $q$ cells on that axis and tracks the final odd multiplicity. It therefore retains within-block cell distinctness. A physical completion with both prescribed axis records belongs to each one-axis relaxation, so

$$
|\text{physical intersection}|\le\min\{\text{left relaxation},\text{right relaxation}\}.
$$

Taking the minimum of rank, row, and column bounds is valid within a fixed record sector because all three bound the same weighted sector matrix for all laws. Summing the four sector bounds is valid by the triangle inequality. None of these valid steps repairs the false $a_{r,s}$ input.

## 4. Independent small-order validation

The new regression `tests/q64_recovered_cubic_quintic_independent_audit.py` does not call `coefficient`, `endpoint_coefficient`, `chain_coefficient`, `chain_sector_squared_coefficient`, `nonfavorable_adjacent_split_coefficient`, or a second closed-form coefficient implementation.

At $q=4$, it explicitly constructs all $4!2^4=384$ signed permutations, forms the physical $HP$ and $PH$ coordinates, and evaluates moments by direct averaging. It checks 20 complete endpoint energies, 12 exact two-axis incidence families, four exact one-axis relaxation families, all four $(r,s)$ chain sectors, and all 40 complement/reversal representatives of the ten orbits. Endpoint checks retain exact extremal rows, so the validation is not a uniform-law calculation.

At $q=8$, it directly averages selected links over all $8!$ permutations after exact sign averaging. The selected supports cover the endpoint record-one L shape, endpoint record-three row triple, compatible and incompatible cubic record-one geometries, quintic row patterns $5$, $4+1$ with zero and nonzero xor, $3+2$, $2+2+1$ with equal and distinct even-group xors, and the full sector-$(1,3,1)$ counterexample. The direct result is compared with the proposed analytic bounds, not with a production moment helper.

The largest $q=4$ dense object is a $560\times4368$ `int16` moment table, below 5 MB. The largest feature table is below 3.4 MB. No giant occurrence matrix is constructed, the run stays far below 8 GB, and the focused test completes in seconds rather than approaching the 30-minute limit.

## 5. Numerical certification of the retained endpoint theorem

Every endpoint-sector square in the regenerated artifact is stored as an integer numerator and denominator. The production code retains `Fraction` values through every product and minimum, rounds each square root upward with `nextafter`, and rounds every accumulated sector sum upward.

The independent test parses the committed rational sector squares, evaluates their square roots with 100-digit `Decimal` arithmetic, and verifies that every displayed artifact coefficient dominates the exact high-precision sum. The largest retained coefficient is

$$
0.7036151810879134<1,
$$

leaving margin greater than $0.296$, vastly larger than floating-point rounding error. No rejected chain entry receives an artifact theorem coefficient.

## 6. Exact rejected entries

The three rejected complement/reversal orbits are represented by

$$
((1,3,5,3),(0,1,3,2)),
\quad
((1,3,5,3),(0,2,2,2)),
\quad
((1,3,5,3),(0,2,3,1)).
$$

Their twelve entries are:

- $((1,3,5,3),(0,1,3,2))$
- $((1,3,5,3),(1,2,2,1))$
- $((3,5,3,1),(1,2,2,1))$
- $((3,5,3,1),(2,3,1,0))$
- $((1,3,5,3),(0,2,2,2))$
- $((1,3,5,3),(1,1,3,1))$
- $((3,5,3,1),(1,3,1,1))$
- $((3,5,3,1),(2,2,2,0))$
- $((1,3,5,3),(0,2,3,1))$
- $((1,3,5,3),(1,1,2,2))$
- $((3,5,3,1),(1,3,2,0))$
- $((3,5,3,1),(2,2,1,1))$

## 7. Registry decision and narrow missing lemma at audit time

At this rejected-proof checkpoint, the registry recorded:

- 342 of 354 coefficient-one-dependent entries repaired;
- recovered-chain quarantine count: 12;
- 796 provisionally supported and 92 open before the separate dual-endpoint caveat;
- 784 conservatively supported and 104 open after withholding the 12 dual-endpoint entries.

The narrow missing result at that checkpoint was a joint endpoint-compatible contraction for the two cubic--quintic links sharing the quintic block. It had to control their product with the full support geometry and record orientation retained; separate scalar bounds $m_rm_s$ are false. A safe generic replacement did not close the two extreme $(1,3,5,3)$ orbits, so this was not a rounding or bookkeeping repair.

This audit does not address the later 80 residual entries, the 12-entry dual-endpoint caveat, the outward-rounded full routing ledger, or the adaptive lift. Its immediate next goal was to prove or falsify the missing joint-link contraction for these twelve quarantined entries before returning to any later residual family.

## 8. Follow-up resolution

The verdict at the top of this note remains the correct verdict on the old independent-maxima proof. A separate replacement theorem now supplies the missing joint contraction by classifying all 15 feasible shared-quintic row/column shapes. Its three canonical q64 coefficients are $0.338286973244$, $0.118636963690$, and $0.314433224343$, so the twelve entries are no longer quarantined. The current registry is 354 of 354 repaired, 808 provisionally supported and 80 open before the dual-endpoint caveat, or 796 conservatively supported and 92 open after withholding it. See Q64_JOINT_RECOVERED_CUBIC_QUINTIC_CONTRACTION.md.
