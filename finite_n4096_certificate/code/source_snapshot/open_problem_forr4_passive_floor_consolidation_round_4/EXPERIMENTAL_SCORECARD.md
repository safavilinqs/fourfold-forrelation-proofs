# Experimental feasibility scorecard

Complete this file for every candidate $N$. A mathematical pass is not automatically an experimental pass.

VERDICT: NOT YET EXPERIMENTALLY CREDIBLE

Status: the mathematical active--passive separation and the active ideal-resource row are complete at $N=4096$. The dated primary-evidence screen is also complete and does not find a current demonstration of the required $4096$-dimensional coherent Sylvester transform at combined contrast above $0.904294855157$. See notes/EXPERIMENTAL_FEASIBILITY_DECISION.md.

## Resource table

| field | required entry |
|---|---|
| $N$ | power-of-two problem size |
| sign modes $M$ | $4N$ |
| ancillary modes | count and role |
| active sample traversals | six, itemized |
| passive dose excluded | exact hard cap |
| state preparation | explicit states and coherence requirements |
| transformation network | Hadamard/controlled operations required |
| receiver | measurement or POVM description |
| repetitions | if any, with total charged dose |
| promise preparation rate | including conditioning/postselection policy |
| loss model | ideal or quantified |
| phase/sign error model | ideal or quantified |
| detector assumptions | efficiency, dark counts, resolution |
| certified decision error | both promise sides |
| theorem margin | distance from the $1/3$ threshold |

### Current $N=4096$ active row

| field | current value |
|---|---|
| $N$, sign modes | $4096$; $M=16{,}384$ in four sign banks |
| ancillary modes | one binary path flag per photon; three photons may run sequentially |
| active traversals | two sign-mask traversals per photon; six total |
| state preparation | $|+\rangle_{\rm path}\otimes|u\rangle_{\rm mode}$ |
| transformation | branch words $D_1HD_2$ and $D_4HD_3H$ |
| receiver | path-$X$ port, mode-insensitive detection, three-flag majority |
| postselection | none; no-click must remain in the record |
| ideal error | $81/256$ |
| scalar robustness | require $g/4-b>0.226073713789$ |
| combined zero-bias contrast | $g>0.904294855157$ |
| nominal 98% detector screen | remaining non-detector contrast $>0.922749852201$; geometric-mean loss $<0.3491602$ dB |
| transform screen | one $H_{4096}$ on the left, two on the right; no reviewed platform demonstration at the required contrast |
| status | resource-complete ideal protocol; not yet experimentally credible |

Here $g=\xi\sqrt{T_LT_R}\,\nu\cos\delta$ combines detection/retention, geometric-mean arm power transmission, task-mode overlap, and differential phase; $b$ is worst-signed additive flag bias. This scalar allocation is not a substitute for a measured platform noise model.

## Feasibility verdict

For each candidate size, label:

- mathematically certified;
- resource-complete but experimentally ambitious;
- plausible with named current capabilities;
- requires an additional robustness corollary; or
- not credible for this paper.

The current $N=4096$ row receives the label `not yet experimentally credible`. This is not a physical impossibility claim. The blocking resource is the measured coherent transform and total contrast, not addressable pixels or detector count.

## Comparison template

| $N$ | $M$ | passive $D=6$ bound | margin | conditioning | active resources | feasibility |
|---:|---:|---:|---:|---:|---|---|
| 256 | 1024 | open | open | open | six traversals | evaluate |
| 512 | 2048 | current witness unsupported | — | — | six traversals | needs different witness geometry |
| 1024 | 4096 | diagnostic only | $0.0002007278$ diagnostic | inherited | six traversals | preferred target |
| 2048 | 8192 | current witness unsupported | — | — | six traversals | needs different witness geometry |
| 4096 | 16384 | classically adaptive passive hard dose six excluded; outward transcript bound $0.268858135059926$ | $0.063475198273408$ below $1/3-10^{-3}$ | $0.002225128406631$ paid once | three photons; six charged traversals; error $81/256$; require $g>0.9042949$ at zero bias | mathematically certified; ideal resources explicit; not yet experimentally credible |

The unsupported labels concern the present $N=q^2$ signed-permutation witness, not the definition of the four-forrelation task at those sizes. The $N=1024$ row remains the smaller experimental upgrade target but needs a new dependency-exact one-batch theorem and its own platform certificate.
