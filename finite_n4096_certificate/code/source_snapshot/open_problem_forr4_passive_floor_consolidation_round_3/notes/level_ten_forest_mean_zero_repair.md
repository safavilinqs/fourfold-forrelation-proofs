# Centered-weight repair of the level-ten forest

Date: 2026-07-15

Status: rigorous removal of the unique level-ten joint saturator and the
intermediate \(\Omega(N^{1/18})\) theorem.  The later
`level_nine_tree_centered_repair.md` removes both remaining trees and proves
the intermediate \(\Omega(N^{1/16})\) theorem; the later complete-image
audit proves the current \(\Omega(N^{1/12})\) theorem.  The
existing-coordinate Stein transfer in every history producing the
six-plus-four forest differentiates an even local weight at the same middle
vertex.  The resulting odd centered factor may be expanded once more by an
exact Gaussian Stein identity.  Every reflection-sensitive graph produced
by that expansion has grouped-entry decay at least (N^{-1}).

Consequently the general passive floor improves from

$$

\Omega(N^{1/20})

\quad\text{to}\quad

\boxed{\Omega(N^{1/18})}.
\tag{0.1}

$$

This remains an asymptotic theorem.  At (N=1024), its bare scale is only

$$
1024^{1/18}=1.4697344923,
\tag{0.2}
$$

so it does not prove passive dose greater than six.

## 1. The exact obstruction

Use layer labels (1,2,3,4).  The level-ten graph from
`high_level_terminal_best_of_two_audit.md` is

$$
\begin{aligned}
T_6:
&\quad a_0-b_0-c_0-d_0,
\qquad a_1-b_1-c_0,\\
P_4:
&\quad a_2-b_2-c_1-d_1.
\end{aligned}
\tag{1.1}
$$

Thus (T_6) is a six-vertex tree and (P_4) is a disjoint four-layer
path.  There are eight Hadamard edges.

Only assigned-suppression-one placements need repair.  The path must be
singleton in every physical entry.  On (T_6), the exact rank-one
entry-respecting partitions are

$$
\begin{aligned}
\pi_0
&=\{\{a_0,b_0\},\{a_1,b_1\},\{c_0\},\{d_0\}\},\\
\pi_1
&=\{\{a_0,b_0\},\{a_1,b_1\},\{c_0,d_0\}\}.
\end{aligned}
\tag{1.2}
$$

Path vertices may share entries with distinct blocks in (1.2), but two path
vertices may not share an entry.  Up to entry relabeling, this gives exactly
282 global dangerous partitions.  On each one,

$$
\delta_{\rm ass}=\delta_{\rm proj}=\frac12.
\tag{1.3}
$$

## 2. The coefficient audit

Write

$$
\gamma=\psi',
\qquad
\mathcal S=\mathcal S_\psi,
\tag{2.1}
$$

for the positive even derivative of the capped odd rounding function and
its even Stein kernel.  The local-function update is the one in the
integration-by-parts branching lemma: a transfer to a new incidence adds a
(gamma) factor, while a transfer to an already occupied incidence
differentiates the factor already there.  This is the exact update in
Lemma 5.2 and Claim 5.3 of
[Bansal--Sinha](https://arxiv.org/abs/2008.07003), specialized to the capped
rounding used here.

The labeled history enumeration for (1.1) is finite and exact:

| quantity | value |
|---|---:|
| potential-twelve initial edge triples | 6 |
| initial triples that reach this terminal graph | 4 |
| complete legal transfer histories | 200 |
| terminal local-weight profiles | 4 |
| fresh transfers in every history | 4 |
| existing-coordinate transfers in every history | 1 |
| derivatives of an existing local weight | 1 |

The differentiated factor is always the incoming-side factor at (c_0).
Depending on whether its first incoming incidence was initial or created by
a transfer, it is exactly one of

$$
h=\gamma'=\psi'',
\qquad
h=\mathcal S'.
\tag{2.2}
$$

Both (gamma) and (mathcal S) are even.  Hence both choices of (h)
are odd and Gaussian-centered:

$$
\mathbb E h(Z)=0,
\qquad Z\sim N(0,1).
\tag{2.3}
$$

All other local factors in the four profiles are undifferentiated
(gamma)'s and (mathcal S)'s.  The transfer coefficient has the same
Hadamard monomial and the same interpolation-time powers

$$
(t_1^2,t_2^2,t_3^1)
\tag{2.4}
$$

for every history.  Thus there is no unsupported scalar-sign cancellation
claim: every history is repaired separately using its centered factor.

## 3. One more exact Stein expansion

For either (h) in (2.2), define its bounded Stein kernel

$$
\mathcal S_h(x)
=e^{x^2/2}\int_x^\infty h(s)e^{-s^2/2}\,ds.
\tag{3.1}
$$

Oddness and smooth Gaussian decay imply that (mathcal S_h) and the fixed
number of derivatives used below are bounded.  If (B) and
((U_b)_{b\in[N]}) are jointly standard Gaussian with

$$
\operatorname{Cov}(B,U_b)=t_2H_{bc_0},
\tag{3.2}
$$

then, whenever (Q) has no explicit dependence on (B),

$$
\mathbb E[h(B)Q]
=\sum_{b=1}^N t_2H_{bc_0}
\mathbb E[\mathcal S_h(B)\,\partial_{U_b}Q].
\tag{3.3}
$$

Here (B) is the latent variable carrying the incoming (c_0) factor.
The terminal transcript factor is (partial_JF(\mu)), with (c_0\in J),
so multilinearity makes it independent of (mu_{c_0}).  The outgoing
local factor at (c_0) uses the other, independent latent variable.
Therefore the remaining integrand (Q) has no explicit (B)-dependence,
and (3.3) applies exactly.

There are two cases.

### 3.1 The neighbor (b) is already marked

The marked second-layer labels are (b_0,b_1,b_2).  Since
(partial_JF) is independent of their biases, (partial_{U_b}) hits
only the existing bounded local factor at (b).  Equation (3.3) adds one
Hadamard edge (b-c_0) and no transcript mark.

- For (b=b_0) or (b=b_1), the edge is parallel to an existing edge.
  The six-vertex component becomes unicyclic.  Its projective decay is at
  least (1/2), including the even-parity cancellation of the parallel
  edge in the binary cut matrix, while (P_4) contributes another
  (1/2).
- For (b=b_2), the new edge joins (T_6) to (P_4).  The resulting
  ten-vertex graph is a tree.  On all 282 dangerous partitions, either the
  assigned occupancy is at least three or an entry-respecting cut has
  binary rank at least three.

Thus every marked-neighbor branch satisfies

$$
\delta_{\rm best}\ge1.
\tag{3.4}
$$

### 3.2 The neighbor (b) is new

Differentiating through the unmarked output bias gives

$$
\partial_{U_b}\partial_JF(\mu)
=\psi'(U_b)\psi(V_b)
\partial_{J\cup\{b\}}F(\mu).
\tag{3.5}
$$

The raw factor (psi(V_b)) is transferred once toward the first layer:

$$
\mathbb E[\psi(V_b)R]
=\sum_{a=1}^N t_1H_{ab}
\mathbb E[\mathcal S_\psi(V_b)\,\partial_{U_a}R].
\tag{3.6}
$$

If (a) is new, the graph has four distinct first-layer vertices and is
reflection-insensitive.  It cancels before absolute values are taken.  If
(a\in\{a_0,a_1,a_2\}), the derivative hits the existing outer
(gamma) factor and the diagram remains reflection-sensitive.  It has one
new second-layer vertex and therefore level eleven.

For completeness, the reflection sign can be checked locally.  If an outer
vertex has boundary degree (d), its local factor is
(psi^{(d)}).  Under first-block reflection, the (d) boundary
Hadamard factors contribute ((-1)^d), while
(psi^{(d)}(-x)=(-1)^{d+1}psi^{(d)}(x)).  The product is (-1) per
distinct outer vertex.  Hence the whole diagram has sign

$$
(-1)^{|V_1|},
\tag{3.7}
$$

which rigorously justifies the first-layer-cardinality sensitivity test,
including repeated outer coordinates.

The new mark (b) can occupy any old physical entry or a new one.  Extending
all 282 dangerous partitions gives 5,295 exact cases for each collection of
the three outer choices.  The two choices (a_0,a_1) create a unicyclic
component and have minimum decay (3/2).  The choice (a_2) joins the two
old components into an eleven-vertex tree and has minimum decay one.  Thus

$$
\delta_{\rm best}\ge1
\tag{3.8}
$$

for every retained new-neighbor branch.

## 4. Exact physical-partition lemma

The finite audit behind (3.4) and (3.8) uses only the following exact
quantities.  For a component (C), physical-entry partition (pi_C),
cycle surplus (ell_C), maximum occupancy (k_C), and maximum
entry-respecting binary cut rank (r_C), it computes

$$
\delta_{\rm ass}
=\frac12\sum_C(\ell_C+k_C-1),
\qquad
\delta_{\rm proj}
=\frac12\sum_C(\ell_C+r_C-1).
\tag{4.1}
$$

The score is the better of these two separately complete contractions.  No
mixed Hilbert/projective induction is used.  Parallel edges are retained in
(ell_C) and reduced modulo two only in the cut adjacency matrix, exactly
as required by the Hadamard graph formula.

The exhaustive counts and minima are:

| branch | exact partitions | minimum (delta_{\rm best}) |
|---|---:|---:|
| duplicate (b_0-c_0) or (b_1-c_0) | (2\cdot282) | (1) |
| bridge (b_2-c_0) | (282) | (1) |
| new (b), existing (a_0,a_1,a_2) | (5{,}295) | (1) |

The enumeration is complete because (1.2) lists every induced rank-one
partition of (T_6), every path vertex is injected into at most one old
block or its own new block, and the new (b) is then inserted into every
old block or a fresh block.  Coordinate and entry relabelings do not alter
any quantity in (4.1).

## 5. The improved transcript theorem

Combining the two cases in Section 3 bounds the original Type-C family by

$$
|\mathfrak C_{\rm C}|
\le C
\left((1+D)^{10}+(1+D)^{11}\right)N^{-1}.
\tag{5.1}
$$

All other level-ten types already have decay at least (N^{-5/8}).  The
complete high-level audit gives (N^{-1}) at levels eleven and twelve.
Levels four through nine retain the accepted (N^{-1/2}), and the two
level-nine trees are now the only limiting high-level types.  Therefore the
task difference obeys

$$
\begin{aligned}
\Delta
\le{}&
\sum_{v=4}^{9}C_v(1+D)^vN^{-1/2}
+C_{10}(1+D)^{10}N^{-5/8}\\
&+C_{\rm C}
\left((1+D)^{10}+(1+D)^{11}\right)N^{-1}
+\sum_{v=11}^{12}C_v(1+D)^vN^{-1}.
\end{aligned}
\tag{5.2}
$$

Set (D=cN^{1/18}).  The level-nine term is constant at the dimension
scale and is made small by choosing (c) sufficiently small.  Every other
term decays; for example,

$$
N^{10/18}N^{-5/8}=N^{-5/72},
\qquad
N^{11/18}N^{-1}=N^{-7/18}.
\tag{5.3}
$$

The accepted promise-conditioning loss remains (O(N^{-1})).  This proves
(0.1).

## 6. Decision for Track B

The level-ten obstruction is removed.  It is neither a scalar-cancellation
counterexample nor a surviving physical-norm obstruction once its forced
centered local derivative is used.

The next asymptotic target identified here was one representative of the
reflected level-nine tree pair.  That audit is now complete and raises the
present mechanism to \(\Omega(N^{1/16})\).  The realistic-size objective
remains separate: the
constants in (5.2) are not controlled sharply enough to exclude passive
dose six at (N=1024).

Reproduction:

- `searches/level_ten_forest_mean_zero_repair.py` enumerates all 200 legal
  coefficient histories, all 282 dangerous physical partitions, and all
  5,295 retained extended partitions.
- `tests/level_ten_forest_mean_zero_repair.py` protects the common centered
  factor, exact branch counts, decay-one minima, and the (1/18) exponent.
