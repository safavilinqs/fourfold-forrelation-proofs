# A two-pass construction and a conditional finite lower bound

Read [the PDF](output/pdf/forr4_n4096_advantage.pdf), built from [main.tex](main.tex).

Three single photons, each traversing two masks, solve the promise at hard dose six with error `81/256`. This construction is proved.

The proposed lower bound is **conditional**. The stored matrix gives `TV <= 0.260969224792207925` only if the complete occurrence-kernel bounds are valid and the parallel probe is block diagonal in parity-support size. Independent review found false premises in several coefficient derivations. The recorded arithmetic therefore does not yet certify a physical lower bound.

Fixed physical photon number is insufficient for the balanced-cut restriction when modes repeat. For example, three photons in one mode have parity degree one, while three photons in distinct blocks have parity degree three. Their superposition has fixed physical photon number but an omitted unbalanced parity cut.

[The audit](../AUDIT.md) gives exact counterexamples and the unresolved adaptive step. [The coefficient guide](code/COEFFICIENTS.md) distinguishes sound reductions, unproved inputs and false intermediate claims. The short note supplies an explicit residual-support lifting argument for its corrected conditional theorem.

Run `make check` and `make build` at the repository root. The frozen source snapshot remains unchanged for reproduction; its historical `proved` and `CERTIFIED` labels are not current theorem verdicts.
