# Active six-dose resource row at $N=4096$

Date: 2026-07-16

Status: proved ideal active protocol. This closes the active mathematical row for the current $N=4096$ paper target. It does not establish laboratory robustness or the passive lower bound.

## Exact protocol

Let $H=H_N$ be the normalized Sylvester Hadamard, let

$$
|u\rangle=N^{-1/2}\sum_i|i\rangle,
$$

and let $D_a=\operatorname{diag}(x^{(a)})$ be the four unknown sign masks. Define

$$
|L_x\rangle=D_2HD_1|u\rangle,
\qquad
|R_x\rangle=HD_3HD_4|u\rangle.
$$

Prepare one photon in

$$
\frac{|0\rangle|L_x\rangle+|1\rangle|R_x\rangle}{\sqrt2}.
$$

The two controlled branch words, written chronologically, are:

| branch | operations | charged sign blocks | public Hadamards | hard dose |
|---|---|---|---:|---:|
| left | $D_1,H,D_2$ | 1, 2 | 1 | 2 |
| right | $D_4,H,D_3,H$ | 4, 3 | 2 | 2 |

Every branch crosses exactly two unknown sample regions. Public Hadamards, coherent routing, delay, and the final path interferometer do not add sample dose.

The folded-state overlap is exactly

$$
\langle L_x|R_x\rangle=F_{4,H}(x).
$$

Measuring Pauli $X$ on the path therefore returns a binary flag $X$ with

$$
\mathbb E[X|x]=F_{4,H}(x).
$$

The mode outcome may be ignored: physically, the receiver recombines the two paths and records the output port while summing over the mode register.

## Error and dose

At either promise boundary $|F_{4,H}|=1/4$, one flag is correct with probability $5/8$. Run three independent flags and take their majority. The worst-case promised error is

$$
\left(\frac38\right)^3
+3\left(\frac58\right)\left(\frac38\right)^2
=\frac{81}{256}
=0.31640625.
$$

The exact error margin is

$$
\frac13-\frac{81}{256}
=\frac{13}{768}
=0.0169270833333.
$$

Each flag has deterministic hard dose two, so three flags have deterministic branchwise hard dose six. Every binary result is retained; there is no postselection or heralding discount.

## $N=4096$ resource ledger

| field | value |
|---|---|
| problem size | $N=4096$ |
| unknown sign blocks | 4 |
| sign modes | $M=4N=16{,}384$ |
| logical mode dimension per flag | 4096 |
| logical path dimension per flag | 2 |
| photons / independent flags | 3 |
| charged traversals per photon | 2 |
| total hard dose | 6 |
| input | $|+\rangle_{\rm path}\otimes|u\rangle_{\rm mode}$ |
| receiver | path-$X$ interferometer, mode-insensitive port detection, majority vote |
| postselection | none |
| worst promised error | $81/256$ |
| error margin | $13/768$ |

The three flags can be run sequentially with the same mode network, so the theorem does not require three copies of the 16,384 sign modes. Parallel execution is optional and does not change the hard-dose accounting.

## Experimental nonclaims

The theorem assumes ideal coherent sign masks, lossless public Hadamards and routing, correct relative phase between arms, ideal single-photon preparation, and ideal detection. It currently supplies no tolerated loss, phase noise, detector inefficiency, dark-count rate, or source impurity. Those quantities must be added before calling $N=4096$ experimentally feasible.

The subsequent scalar allocation in `ACTIVE_SIX_ROBUSTNESS_GATE.md` requires
combined retained contrast above $0.904294855157$ when additive bias is zero.
It quantifies the available budget but is not a device-specific robustness
theorem.

The protocol is pointwise for every promised input and does not depend on the passive hard distribution or its conditioning probability.

## Reproduction

Run:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/active_six_resource_row.py --write-artifact
```

The artifact is `artifacts/active_six_resource_row.json`. The original folded-overlap regression remains `../open_problem_forr4_passive_floor/active_six_check.py`.
