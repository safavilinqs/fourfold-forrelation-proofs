# $N=4096$ experimental feasibility decision

VERDICT: NOT YET EXPERIMENTALLY CREDIBLE

Date: 2026-07-18

The active protocol is resource-complete as an ideal experiment, but the reviewed primary evidence does not demonstrate the required $4096$-dimensional coherent Sylvester transform with enough end-to-end contrast. This is an evidence-based paper-framing decision, not a physical impossibility theorem.

## What the experiment actually needs

Each of three single-photon flags needs a coherent superposition of 4096 logical modes and one binary path ancilla. The left path applies one public $H_{4096}$ between sign banks 1 and 2. The right path applies two public $H_{4096}$ transforms between sign banks 4 and 3 and before recombination. Four independently specified 4096-coordinate phase-sign patterns must be addressed, every photon must cross exactly two charged sign banks, the paths must remain mutually coherent, and the final detector must retain no-click events rather than retrying them away.

The detector resolves only the output path port; it does not need 16,384 pixels. Conversely, having 16,384 addressable pixels does not implement a coherent $H_{4096}$ or certify its loss and mode overlap.

## Quantitative budget

For zero additive bias, the three-flag majority decoder requires

$$
g=\xi\sqrt{T_LT_R}\,\nu\cos\delta
>
0.904294855156832,
$$

where $\xi$ is detection and retention efficiency, $T_L,T_R$ are complete arm power transmissions, $\nu$ is task-mode overlap, and $\delta$ is differential phase error.

A fiber-coupled SNSPD has demonstrated nominal $98.0\%\pm0.5\%$ system detection efficiency at 1550 nm [in this NIST/Optica experiment](https://www.nist.gov/publications/superconducting-nanowire-single-photon-detectors-98-system-detection-efficiency-1550-nm). Using the nominal $\xi=0.98$ as an optimistic benchmark and setting bias, phase error, and mode mismatch to zero leaves

$$
\sqrt{T_LT_R}\,\nu\cos\delta
>
0.922749852200849.
$$

Thus every non-detector imperfection combined receives less than

$$
-10\log_{10}(0.922749852200849)
=0.349160156943327\ \mathrm{dB}
$$

of geometric-mean power-loss budget. This is only about $0.35$ dB and already assumes perfect sign masks, routing, overlap, phase, and zero additive bias.

A butterfly realization of $H_{4096}$ has twelve stages. With one transform on the left and two on the right, the geometric-mean transform depth is eighteen stages. If every stage had the same power transmission and all other components were perfect, each stage would need

$$
\eta_{m stage}>0.995543454598465,
$$

or less than $0.0193977865$ dB loss per stage. This is an illustrative allocation, not a claim that every physical Sylvester implementation must use that architecture.

## Primary-evidence screen

| platform evidence | demonstrated scale | relevance to this protocol | decision |
|---|---|---|---|
| [Ultrafast time-bin circuits, PRL 2024](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.133.090601) | 362 programmable unitaries through dimension 8; passive networks through 36 modes; reported fidelities above 97% | Time bins are the cleanest conceptual match for sequential sign phases and stable collinear interference | Strong candidate architecture, but direct coherent dimension is more than two orders of magnitude short; fidelity is not the required end-to-end contrast |
| [Borealis time-domain processor, Nature 2022](https://www.nature.com/articles/s41586-022-04725-x) | 216 temporal modes with three-loop nonuniversal connectivity | Shows that hundreds of coherent time bins, switching, and detection are real | Does not implement a dense or Sylvester $H_{4096}$; the paper explicitly avoids universal depth because loss accumulates |
| [Complex-medium spatial circuits, Nature Physics 2024](https://www.nature.com/articles/s41567-023-02319-6) | mixer dimension about 200; characterized Fourier circuits through dimension 7, with 81.4% fidelity at dimension 7 | Shows programmable high-dimensional spatial mixing | Physical mixer dimension is not a certified circuit dimension; no $4096$-mode transform or compatible contrast is demonstrated |
| [Spatial MPLC characterization](https://doi.org/10.1016/j.yofte.2016.09.005) | 10-mode multiplexers, $-26$ dB crosstalk, insertion loss reported below 4 dB | Shows mature mode conversion | The published loss statement is too coarse to certify the entire roughly 0.35 dB remaining budget, and the mode count is far short |
| [105-mode SLM communication](https://www.nature.com/articles/srep27674) | 35 spatial modes across three wavelengths, detected as 105 carriers | Shows many addressable classical carriers | It is modal encoding and projection, not one coherent $4096$-dimensional unitary on a single photon |
| [400,000-pixel SNSPD camera](https://www.nist.gov/publications/superconducting-nanowire-single-photon-camera-400000-pixels) | 400,000 detector pixels | Shows that raw pixel count need not be the limiting number | The protocol needs a two-port coherent receiver, not a large camera; the transform and contrast remain the bottleneck |

The strongest plausible mapping is a custom time-bin device with fast coordinate-wise phase modulation and fixed low-loss delay/butterfly networks. The reviewed evidence makes that a research direction, not a platform-certified implementation today. Frequency bins and spatial modes remain possible alternatives, but no reviewed result closes both the dimension and contrast requirements.

## What would change the verdict

A platform proposal passes only if it supplies all of the following with measured uncertainty:

1. A heralded or deterministic single photon in a coherent uniform superposition over 4096 modes, with any heralding completed before charged sample access.
2. Four independently programmable or fabricated 4096-coordinate sign banks and an itemized optical path proving exactly two charged crossings per retained flag.
3. A measured transfer matrix or task-specific certificate for one $H_{4096}$ use on the left and two on the right, including loss, leakage, and phase.
4. A path-$X$ interferometer and mode-insensitive two-port receiver whose no-click and ambiguous outcomes are retained and assigned according to the declared rule.
5. A conservative lower confidence bound satisfying $g/4-b>0.226073713789208$ after detector inefficiency, complete arm transmissions, task-mode overlap, differential phase, dark counts, multiphoton contamination, and sign errors are combined.
6. Three independent flags or a justified correlated-noise analysis proving majority error below $1/3$ under total hard dose six.

Raw mode count, simulation, component-by-component nominal efficiency, or fidelity without throughput does not pass this gate.

## Paper decision

The $N=4096$ result should be presented as a rigorous finite-size theory separation with a fully specified ideal active protocol and a quantitative hardware target. It should not be advertised as an experimentally feasible or demonstrated separation on current evidence.

If present-day experimental credibility is mandatory for the headline, Round 4 has not met that part of the goal. The narrow next project is then an $N=1024$ upgrade: construct a dependency-exact outward $q=32$ one-batch ledger below $1/3-10^{-3}$ and combine it with a platform design that passes the same measured contrast inequality. The existing $q=32$ total is only diagnostic and misses the reserve gate by about $7.99\times10^{-4}$ even before theoremization. The multiplier-one adaptive theorem would apply once a valid one-batch ledger exists.

## Reproduction

Run:

    /opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_experimental_feasibility_gate.py
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_experimental_feasibility_gate.py

The artifact `artifacts/q64_experimental_feasibility_gate.json` records the arithmetic and the primary-evidence snapshot. Its hardware entries are a dated engineering screen, not theorem dependencies.
