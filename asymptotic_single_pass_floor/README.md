# Asymptotic parallel-query bound

Read [main.pdf](main.pdf), built from [main.tex](main.tex).

The proof establishes hard dose `D >= (2/15) N^(1/8)` for powers of two `N >= 2^30`, for one parallel quantum interrogation. Idlers, repeated modes, mixtures and coherence between different signal-number sectors are allowed. Predetermined batches can be combined into one interrogation.

The presentation replaces the former cut construction by an explicit matrix factorization. Section 6 states exactly the Fourier inequality still needed for outcome-dependent fresh batches. The former adaptive proof was incomplete; this is not a proof that its conclusion is false. See [the audit](../AUDIT.md).

Run `python3 verify_constants.py` for exact arithmetic, or `make check` at the repository root for all supported checks. These computations do not independently verify Gaussian or operator lemmas.
