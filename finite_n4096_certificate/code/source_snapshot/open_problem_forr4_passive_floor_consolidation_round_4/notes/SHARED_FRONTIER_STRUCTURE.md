# Shared-frontier occupation structure

Date: 2026-07-15

## Result

The leading 51 unresolved symmetry orbits are not 51 independent contraction problems. They carry 90.5467 percent of the unresolved Perron contribution and collapse exactly to one dose-six occupation kernel with:

- 198 profile-split entries;
- six unordered block-degree patterns;
- 125 participating occupation states;
- 241 distinct undirected occupation edges; and
- eight connected components.

This is an exact combinatorial reduction of the target. It is not an arbitrary-law coefficient bound and does not improve the passive lower bound by itself.

## Why the reduction occurs

Every entry in the frontier is centrally balanced. Write its odd block-degree profile as $a=(a_0,a_1,a_2,a_3)$ and its selected occurrence split as $s$. Every compatible occupation pair has the form

$$
n=s+k,
\qquad
m=(a-s)+k,
$$

where $k$ is the shared unmarked occupation.

For the 100 degree-ten entries, $|s|=5$, so the hard-dose-six condition leaves only $|k|=0$ or $1$. Each entry therefore has exactly five compatible terms: one on the dose-five layer and four on the dose-six layer. For the 98 degree-twelve entries, $|s|=6$, so $k=0$ and each entry has exactly one compatible term on the dose-six layer.

Thus the nominal entry inventory produces only $500+98=598$ oriented compatible terms. Complement and reversal symmetries, together with overlap between distinct profiles, merge them into 241 undirected edges. Twenty-five edges are shared by both the degree-ten and degree-twelve families.

## Six profile families

| unordered profile | degree | orbits | entries | occupation edges | current Perron contribution |
|---|---:|---:|---:|---:|---:|
| $(5,3,1,1)$ | 10 | 21 | 84 | 189 | $0.0175265235$ |
| $(5,3,3,1)$ | 12 | 15 | 60 | 30 | $0.0055625540$ |
| $(7,3,1,1)$ | 12 | 4 | 16 | 8 | $0.0015747806$ |
| $(7,1,1,1)$ | 10 | 4 | 16 | 40 | $0.0015165587$ |
| $(5,5,1,1)$ | 12 | 5 | 16 | 8 | $0.0013234466$ |
| $(3,3,3,3)$ | 12 | 2 | 6 | 3 | $0.0007226071$ |

The first pattern alone carries 62.1 percent of this frontier's current Perron contribution. The first two patterns carry 81.8 percent. These percentages are routing diagnostics from the non-theorem coefficient ledger, not certified bounds.

## Component structure

The occupation graph has four dose-five components, each on 11 states, and four dose-six components on 24, 20, 19, and 18 states. Consequently a shared contraction can be assembled and checked component by component. It does not require a dense $210\times210$ argument, much less 198 independent scalar theorems.

This is the correct interface between the configuration-level signed-permutation problem and the scalar Perron certificate:

1. derive one shared physical or operator bound for the six profile families;
2. push its edge weights into these eight occupation components;
3. take the largest certified component Perron value; and
4. insert the resulting joint frontier bound into the full finite-size ledger.

The contraction must preserve cancellation among entries that share an occupation edge. Replacing it by 198 independently optimized coefficients would discard the reduction.

## Next target

The immediate target is the dominant $(5,3,1,1)$ family, followed by the $(5,3,3,1)$ family. A useful compressed representation must reproduce the full coordinate-uniform $q=4$ operator for these families and expose the same eight-component edge weights at $q=32$.

Once those two families are controlled, the same machinery should absorb the remaining four patterns. If it cannot produce a joint bound with visible reserve, the signed-permutation witness fails the declared Round 4 retain gate and the project pivots.

## Reproduction

Run:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/shared_frontier_structure.py --write-artifact
```

The committed machine-readable audit is `artifacts/shared_frontier_structure.json`. The Round 4 regression regenerates it and compares it byte for byte.
