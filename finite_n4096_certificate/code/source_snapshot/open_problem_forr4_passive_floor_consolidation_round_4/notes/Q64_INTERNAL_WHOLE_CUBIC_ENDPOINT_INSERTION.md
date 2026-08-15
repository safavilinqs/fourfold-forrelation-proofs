# Internal whole-cubic endpoint insertion at q=64

Date: 2026-07-16

Status: rigorous arbitrary-diagonal one-batch coefficients for 16 degree-twelve entries in four complement/reversal orbits. This raises the q64 theorem count from 324 to 340. It does not prove the remaining 548 entries, intervalize the ledger, or establish the adaptive lift.

## Result

The theorem applies whenever the unique unsplit cubic is adjacent to the singleton and lies entirely on the same side of the cut. The four canonical cuts are

$$
\begin{aligned}
(1,3,3,5)&:(0,0,2,4),\\
(1,3,5,3)&:(0,0,4,2),\\
(3,1,3,5)&:(1,1,3,1),\\
(3,3,1,5)&:(1,3,1,1).
\end{aligned}
$$

All 16 cuts have an extreme $1|4$ quintic split. Their arbitrary-law coefficient is

$$
\boxed{0.113036239514}.
$$

The same proof schema would give coefficient $0.281119075921$ for a balanced $2|3$ quintic cut, although no remaining cut with this topology is balanced.

These degree-twelve entries have zero current Perron sensitivity. Insertion therefore leaves

$$
P_{\rm total}=0.324766326953,
\qquad
{1\over3}-P_{\rm total}=0.00856700638000.
$$

The remaining quintic inventory is 88 entries: 56 extreme and 32 balanced. Its local-slice proxy remains $0.325037612931$, leaving $0.00729572040278$ beyond the declared $10^{-3}$ allowance.

## Complete the physical link kernels

Ignore the within-support distinctness conditions temporarily. Every signed-permutation link moment has the form

$$
\mathbb E[U_rV_c]
$$

for unit-modulus physical character features on the two sides of the cut. It is therefore a cross Gram whose trace-class Schur-multiplier norm is at most one. The product of the three completed link kernels also has norm at most one.

The unsplit cubic introduces no cross-cut within-support mask. Only the split cubic and split quintic masks must be restored.

## Split-cubic mask

For a cubic split as one cell against a pair, the mask is

$$
D(x,E)=\mathbf 1_{x\notin E}.
$$

The centered-vector factorization for a singleton against a fixed two-set gives

$$
\gamma_{N,2}
=1-{2\over N}
+\sqrt{2\left(1-{2\over N}\right)
              \left(1-{1\over N}\right)}.
$$

At $N=4096$,

$$
\gamma_{N,2}=2.41320737011.
$$

## Split-quintic masks

For an extreme split, the quintic mask is a singleton against a four-set. Its factor is

$$
\gamma_{N,4}=2.99780260018.
$$

For a balanced split, write the selected pair in deterministic order $(e_1,e_2)$. Then

$$
\mathbf 1_{\{e_1,e_2\}\cap F=\varnothing}
=\mathbf 1_{e_1\notin F}\mathbf 1_{e_2\notin F}.
$$

Each singleton-versus-triple mask costs $\gamma_{N,3}$, so the balanced factor is at most

$$
\gamma_{N,3}^2=7.45548065275.
$$

This sequential factorization is valid for unordered support pairs because the ordering is fixed deterministically before the Schur features are assigned.

## Same-side whole-cubic endpoint

Let $T$ be the unsplit cubic and $a$ the adjacent singleton. Since both lie on one side of the cut, their endpoint moment is a row-only or column-only scalar:

$$
M_{31}(T,a)=v_3(T)H_N(\xi(T),a).
$$

The endpoint amplitude satisfies $|v_3(T)|\le1$, and the normalized Walsh entry has modulus $1/q$. This supplies the scalar factor $1/q$ independently of the other two completed link kernels.

Combining the endpoint scalar with the distinctness factors proves

$$
\boxed{
c_{1|4}
\le{\gamma_{N,2}\gamma_{N,4}\over q}
=0.113036239514
}
$$

and

$$
\boxed{
c_{2|3}
\le{\gamma_{N,2}\gamma_{N,3}^2\over q}
=0.281119075921.
}
$$

No physical law is duplicated: the completed links and masks are composed as Schur multipliers on the full occurrence matrix. The proof is therefore uniform over arbitrary correlated diagonal laws.

## Regression scope

The regression protects the four-orbit topology partition, both coefficient identities, deterministic artifact output, and the exact $1/q$ magnitude of every sampled whole-cubic endpoint moment at $q=4$ and $q=8$. The centered-mask formula and completed-link Gram principle are inherited proved dependencies.

Reproduce with:

    /opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_internal_whole_cubic_endpoint_insertion.py --write-artifact
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_internal_whole_cubic_endpoint_insertion.py
