# Bibliography notes

Date checked: 2026-07-19

Every conventional bibliography citation in the manuscript is resolved below. A source supports only the stated role; none of the external papers is a dependency of the finite passive certificate.

## `AaronsonAmbainis2015`

Scott Aaronson and Andris Ambainis, “Forrelation: A Problem that Optimally Separates Quantum from Classical Computing,” STOC 2015, 307–316, DOI [10.1145/2746539.2746547](https://doi.org/10.1145/2746539.2746547), arXiv:1411.5729.

Claim scope: origin and query-complexity interpretation of forrelation. This source does not establish the passive sensing theorem or the finite $N=4096$ constants.

## `BansalSinha2021`

Nikhil Bansal and Makrand Sinha, “$k$-Forrelation Optimally Separates Quantum and Classical Query Complexity,” STOC 2021, 1303–1316, DOI [10.1145/3406325.3451040](https://doi.org/10.1145/3406325.3451040), arXiv:2008.07003.

Claim scope: the $k$-fold forrelation lineage and abstract coherent-query context. Its access model and bounds are not substituted for the fresh-batch sensing proof.

## `SafaviNaeini2026`

Amir Safavi-Naeini, “A Hard-Dose Floor for Passive Sensing of Four-Fold Forrelation: Six Coherent Passes Beat Every Fresh-Probe Strategy,” companion manuscript, 2026. The internal source package is preserved at `scratch/20260718_forr4_floor_paper`.

Claim scope: the provisional asymptotic passive $D\ge(2/15)N^{1/8}$ hard-dose result for the full fresh-probe model, the active dose-six upper bound, and the access-model definitions. The finite $N=4096$ theorem for number-sector-incoherent fresh batches is established by the present source snapshot.

## `GirishEtAl2024`

Uma Girish, Makrand Sinha, Avishay Tal, and Kewen Wu, “The Power of Adaptivity in Quantum Query Algorithms,” STOC 2024, 1488–1497, DOI [10.1145/3618260.3649621](https://doi.org/10.1145/3618260.3649621), arXiv:2311.16057.

Claim scope: conceptual comparison with query-adaptive algorithms. The paper uses a different query model and does not prove the adaptive lift used here.

## `AharonovCotlerQi2022`

Dorit Aharonov, Jordan Cotler, and Xiao-Liang Qi, “Quantum Algorithmic Measurement,” Nature Communications 13, 887 (2022), DOI [10.1038/s41467-021-27922-0](https://doi.org/10.1038/s41467-021-27922-0), arXiv:2101.04634.

Claim scope: conceptual comparison between coherent and incoherent experimental access. Its QUALM definitions are not identified with the passive and active classes of the present theorem.

## `ChenEtAl2022`

Sitan Chen, Jordan Cotler, Hsin-Yuan Huang, and Jerry Li, “Exponential Separations between Learning with and without Quantum Memory,” FOCS 2021 proceedings published in 2022, 574–585, DOI [10.1109/FOCS52979.2021.00063](https://doi.org/10.1109/FOCS52979.2021.00063), arXiv:2111.05881.

Claim scope: conceptual evidence that coherent quantum memory can change learning power. The learning problems and memory model are distinct from the four-bank sensing task.

## `ChenEtAl1988`

Yan-Song Chen, Shi-Hai Zheng, Bi-Zhen Dong, De-Hua Li, and Guo-Zhen Yang, “Optical Walsh–Hadamard Transform for Orders 32 and 64,” Applied Optics 27, 2608–2611 (1988), DOI [10.1364/AO.27.002608](https://doi.org/10.1364/AO.27.002608).

Claim scope: a coherent optical Walsh–Hadamard implementation at orders 32 and 64 using a holographic mask and Fourier lenses. This is transform-principle evidence and not an $H_{4096}$ benchmark.

## `ReddyEtAl2020`

Dileep V. Reddy, Robert R. Nerem, Sae Woo Nam, Richard P. Mirin, and Varun B. Verma, “Superconducting Nanowire Single-Photon Detectors with 98% System Detection Efficiency at 1550 nm,” Optica 7, 1649–1653 (2020), DOI [10.1364/OPTICA.400751](https://doi.org/10.1364/OPTICA.400751).

Claim scope: the nominal $98.0\pm0.5$ percent system-detection-efficiency benchmark used in the illustrative scalar loss allocation. The detector number is not an end-to-end active-protocol efficiency.

## `BouchardEtAl2024`

Frédéric Bouchard, Kate Fenwick, Kent Bonsma-Fisher, Duncan England, Philip J. Bustard, Khabat Heshami, and Benjamin Sussman, “Programmable Photonic Quantum Circuits with Ultrafast Time-Bin Encoding,” Physical Review Letters 133, 090601 (2024), DOI [10.1103/PhysRevLett.133.090601](https://doi.org/10.1103/PhysRevLett.133.090601).

Claim scope: 362 programmed unitaries through dimension eight, a passive network through 36 modes, reported fidelity above 97 percent, and multi-day passive phase stability. The reported circuits are not a coherent $H_{4096}$ implementation.

## `MadsenEtAl2022`

Lars S. Madsen et al., “Quantum Computational Advantage with a Programmable Photonic Processor,” Nature 606, 75–81 (2022), DOI [10.1038/s41586-022-04725-x](https://doi.org/10.1038/s41586-022-04725-x).

Claim scope: a programmable time-domain Gaussian processor with 216 squeezed modes and three loop layers. The architecture, state class, and transformation are different from the single-photon Sylvester circuit.

## `GoelEtAl2024`

Suraj Goel et al., “Inverse Design of High-Dimensional Quantum Optical Circuits in a Complex Medium,” Nature Physics 20, 232–239 (2024), DOI [10.1038/s41567-023-02319-6](https://doi.org/10.1038/s41567-023-02319-6).

Claim scope: high-dimensional spatial-mode circuits embedded in an ambient multimode mixer, with quantum-gate demonstrations through dimension seven. The source does not certify $H_{4096}$ or the required retained contrast.

## `LiuEtAl2026`

Hua-Liang Liu, Hao Su, Yu-Hao Deng et al., “Gaussian Boson Sampling with 1,024 Squeezed States in 8,176 Modes,” Nature 653, 687–692 (2026), DOI [10.1038/s41586-026-10523-6](https://doi.org/10.1038/s41586-026-10523-6).

Claim scope: a recent raw optical scale benchmark in a Gaussian boson-sampling architecture. The presence of 8,176 modes and 4096 connectivity does not establish a coherent single-photon $H_{4096}$ transform at the contrast required here.
