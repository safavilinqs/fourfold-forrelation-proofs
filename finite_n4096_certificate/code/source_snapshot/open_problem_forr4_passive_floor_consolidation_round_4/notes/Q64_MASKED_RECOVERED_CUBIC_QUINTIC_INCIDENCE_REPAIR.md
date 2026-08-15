# The $q=64$ masked recovered cubic--quintic endpoint repair

Date: 2026-07-17

Status: arbitrary-correlated-diagonal one-batch theorem for 28 recovered cubic--quintic endpoint-row entries. The former claim for twelve $(1,3,5,3)$ four-sector chain entries is rejected in this module; a separate joint shared-quintic theorem now proves them. The independent verdict on the old proof is in Q64_RECOVERED_CUBIC_QUINTIC_INDEPENDENT_AUDIT.md, and the replacement is in Q64_JOINT_RECOVERED_CUBIC_QUINTIC_CONTRACTION.md.

## Physical record sectors

Every canonical profile begins with an endpoint singleton. The first link forces record one, and the other two common link records lie in $\{1,3\}$. Thus every entry splits into four physical sectors

$$
\rho=(1,r,s).
$$

All support extensions are counted inside their prescribed row and column records. A completion is the disjoint complement of the fixed occurrence partial inside one simple bipartite support, so every cross-cut distinctness mask is retained.

## Complete endpoint rows

The retained entries have canonical profile $(1,3,3,5)$ or $(1,5,3,3)$. For the first higher block, let $E_r(k,\epsilon)$ be the complete singleton--cubic or singleton--quintic squared row energy in outgoing record $r$, where $k$ higher-block cells and $\epsilon$ singleton cells are fixed.

The cubic fixed-one slice is record-resolved:

$$
E_1(1,1)=\frac3{q^2},
\qquad
E_3(1,1)=\frac{q^2-4}{2q^2}.
$$

The other cubic slices use the safe unrefined total energy for either record; they are not described as record-resolved. When the singleton is variable, summing over its exactly $q^2$ possible cells contributes a safe factor $q^2$. The quintic energies are exact rational record-resolved endpoint maxima for the enumerated slices or stated rational upper bounds for the fixed-four variable-singleton slice.

For the two remaining blocks, write $D_i^\rho(t_i)$ for the exact physical completion incidence at the selected subset size. After removing the endpoint link already contained in $E_r$, the remaining two pointwise moment bounds are $m_r$ and $m_s$. A complete row has squared coefficient at most

$$
E_r(k,\epsilon)D_3^\rho(t_3)D_4^\rho(t_4)(m_rm_s)^2.
$$

Transposition gives the complementary complete-column bound. The theorem takes the smaller exact rational value inside each fixed record sector, applies an outward-rounded square root, and then sums the four sector bounds.

## Corrected cubic--cubic record-one scope

The refined cubic--cubic record-one value

$$
\frac{q+2}{q(q-1)(q-2)}
$$

requires the left cubic's exterior incoming record also to be one. It is valid for the record-one cubic--cubic link in profile $(1,3,3,5)$ and in the compatible sectors of $(1,5,3,3)$. In sector $(1,3,1)$ of profile $(1,5,3,3)$ that exterior record is three, so the audited implementation now uses the universal physical value $1/q$. The resulting seven orbit coefficients are

$$
0.0162888571820,
\quad 0.0264913195629,
\quad 0.0331908274616,
\quad 0.0506481680416,
\quad 0.110898445736,
\quad 0.413910743027,
\quad 0.703615181088.
$$

All remain strictly below one.

## Rejected four-sector chain

The former twelve-entry chain used

$$
a_{r,s}=e_rm_rm_s
$$

as a maximum physical entry. At $q=8$, the exact sector-$(1,3,1)$ supports

$$
\{0\},
\quad
\{0,1,2\},
\quad
\{0,1,8,16,24\},
\quad
\{0,1,2\}
$$

have chain entry $1/17920$, while the claimed maximum is $1/25088$. The maximum-entry premise is therefore false even though the subsequent rank and incidence norm estimates are valid conditional on such a premise. These twelve entries are not included in `repaired_entries()` or the theorem artifact.

## Exact regression

The independent regression directly constructs the physical signed-permutation moments and masks at small order. It verifies all ten complement/reversal orbits at $q=4$, all relevant endpoint energies, exact physical and relaxed incidences, selected $q=8$ link geometries, the rejecting four-block counterexample, and independent high-precision outward rounding for all retained artifact coefficients.

Reproduce with:

    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_recovered_cubic_quintic_independent_audit.py
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_masked_recovered_cubic_quintic_incidence_repair.py

## Replacement theorem

This note continues to certify only the 28 endpoint-row entries and continues to reject its old independent-maxima chain mechanism. The separate joint shared-quintic theorem now proves the twelve formerly quarantined $(1,3,5,3)$ entries with canonical coefficients $0.338286973244$, $0.118636963690$, and $0.314433224343$. See Q64_JOINT_RECOVERED_CUBIC_QUINTIC_CONTRACTION.md.
