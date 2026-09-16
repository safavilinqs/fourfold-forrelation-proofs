# Asymptotic parallel-query bound

**For the current self-contained exposition, read [the complete PDF](../fourfold_forrelation.pdf).** This directory preserves the earlier separate note and its verification history.

Read [main.pdf](main.pdf), built from [main.tex](main.tex).

The proof establishes hard dose `D >= (2/15) N^(1/8)` for powers of two `N >= 2^30`, for one parallel quantum interrogation. Idlers, repeated modes, mixtures and coherence between different signal-number sectors are allowed. Batches whose inputs are jointly prepared before any output can influence an unqueried signal can be combined into one interrogation.

The proof uses an explicit matrix factorization and a centered measurement effect. Only even Fourier levels contribute, giving an average acceptance gap below `1/12`; correctness on the promises requires a gap above `1/5`. Section 6 states a sufficient Fourier inequality still unproved for outcome-dependent fresh batches. See [the audit](../AUDIT.md) for the gap in the former adaptive argument.

Run `python3 verify_constants.py` for exact arithmetic, or `make check` at the repository root for all supported checks. These computations do not independently verify Gaussian or operator lemmas.
