# Centered endpoint repair of the reflected level-nine trees

Date: 2026-07-15

Status: rigorous removal of the final two joint saturators in the
level-nine-through-twelve audit and the intermediate \(N^{1/16}\) theorem.
The later `low_level_terminal_centered_repair.md` proves the current
\(N^{1/12}\) theorem.  Every legal coefficient history
for the upper-branching tree contains a forced odd centered derivative at its
outer endpoint.  Expanding that factor by one exact Gaussian Stein identity
gives grouped-entry decay at least \(N^{-1}\) in every retained branch.  The
layer-reflected tree has the same retained minima, while its one all-fresh
branch cancels from the planted-minus-null task difference.

Consequently the general passive floor improves from

$$

\Omega(N^{1/18})

\quad\text{to}\quad

\boxed{\Omega(N^{1/16})}.
\tag{0.1}

$$

This is still an asymptotic theorem.  Its bare scale at \(N=1024\) is

$$
1024^{1/16}=1.5422108254,
\tag{0.2}
$$

so it does not by itself prove passive hard dose greater than six.

## 1. The two exact obstructions

Use layer labels \(A,B,C,D\).  The upper-branching Type-A tree is

$$
\begin{aligned}
E_{AB}&=\{a_0b_0,a_0b_1\},\\
E_{BC}&=\{b_0c_0,b_0c_1,b_1c_2\},\\
E_{CD}&=\{c_0d_0,c_1d_1,c_2d_2\}.
\end{aligned}
\tag{1.1}
$$

It is a connected nine-vertex tree with eight Hadamard edges.  Type B is
its layer reflection:

$$
\begin{aligned}
E_{AB}&=\{a_0b_0,a_1b_1,a_2b_2\},\\
E_{BC}&=\{b_0c_0,b_1c_0,b_2c_1\},\\
E_{CD}&=\{c_0d_0,c_1d_0\}.
\end{aligned}
\tag{1.2}
$$

Only assigned-suppression-one placements can miss the target
\(N^{-9/16}\).  Maximum occupancy is then two.  For Type A, the exact
entry-respecting cut-rank-two partitions are

$$
\begin{aligned}
&\{a_0\},\{b_0\},\{b_1\},\{c_0,d_0\},\{c_1,d_1\},\{c_2,d_2\};\\
&\{a_0\},\{b_0,b_1\},\{c_0,d_0\},\{c_1,d_1\},\{c_2,d_2\};\\
&\{a_0,b_0\},\{b_1\},\{c_0,d_0\},\{c_1,d_1\},\{c_2,d_2\};\\
&\{a_0,b_1\},\{b_0\},\{c_0,d_0\},\{c_1,d_1\},\{c_2\},\{d_2\};\\
&\{a_0,b_1\},\{b_0\},\{c_0,d_0\},\{c_1,d_1\},\{c_2,d_2\};\\
&\{a_0,b_1\},\{b_0,c_2\},\{c_0,d_0\},\{c_1,d_1\},\{d_2\};\\
&\{a_0,b_1\},\{b_0,d_2\},\{c_0,d_0\},\{c_1,d_1\},\{c_2\}.
\end{aligned}
\tag{1.3}
$$

Both the global assigned bound and the all-projective grouped-entry bound
give only \(N^{-1/2}\) on each partition.  Reflecting (1.3) gives all seven
dangerous Type-B partitions.

## 2. Exact coefficient histories

Write

$$
\gamma=\psi',
\qquad
\mathcal S=\mathcal S_\psi,
\tag{2.1}
$$

where \(\psi\) is the capped odd rounding function.  A transfer to a new
incidence adds a \(\gamma\) factor, while a transfer to an incidence already
occupied at that coordinate differentiates its local factor.  This is the
same coefficient update used in the accepted level-ten repair.

The labeled Type-A enumeration is finite and exact:

| quantity | value |
|---|---:|
| potential-twelve initial edge triples | 6 |
| initial triples that reach Type A | 4 |
| complete legal transfer histories | 200 |
| terminal local-weight profiles | 4 |
| fresh transfers in every history | 3 |
| existing-coordinate transfers in every history | 2 |
| differentiated local weights in every history | 2 |

Every profile has exactly the same derivative sites:

$$
(A,a_0)
\quad\text{and}\quad
(A,b_0).
\tag{2.2}
$$

At \(a_0\), the factor is always

$$
h=\gamma'=\psi''.
\tag{2.3}
$$

At \(b_0\), it is either \(\gamma'\) or \(\mathcal S'\).  Both sites are
useful, but the proof needs only the uniform outer factor (2.3).  It is odd
and Gaussian-centered:

$$
\mathbb E h(Z)=0,
\qquad Z\sim N(0,1).
\tag{2.4}
$$

Every history has interpolation-time exponents \((1,2,2)\).  The repair is
performed history by history, so the argument does not assume cancellation
between unrelated scalar coefficients.

## 3. Exact centered expansion of Type A

Define the bounded Stein kernel

$$
\mathcal S_h(x)
=e^{x^2/2}\int_x^\infty h(s)e^{-s^2/2}\,ds.
\tag{3.1}
$$

If \(X_{a_0}\) is the Gaussian carrying (2.3), then the marked-coordinate
multilinearity of the transcript derivative makes the rest of the
integrand independent of its bias.  Thus the exact Stein identity applies:

$$
\mathbb E[h(X_{a_0})Q]
=\sum_b t_1H_{a_0b}
\mathbb E[\mathcal S_h(X_{a_0})\,\partial_{Y_b}Q].
\tag{3.2}
$$

The derivative either hits a marked neighbor or creates a fresh mark and a
raw odd factor.  In the latter case that factor is transferred right until
it hits an occupied coordinate or reaches layer \(D\).  This gives exactly
four branch families:

1. \(b=b_0\) or \(b_1\): add a parallel edge and stay at level nine;
2. \(b\) new and \(c\in\{c_0,c_1,c_2\}\): level ten;
3. \(b,c\) new and \(d\in\{d_0,d_1,d_2\}\): level eleven;
4. \(b,c,d\) all new: level twelve.

New marks may occupy any existing physical entry or a fresh entry.  Extending
all seven partitions in (1.3) gives the following exhaustive audit:

| Type-A branch | exact scored cases | minimum decay exponent |
|---|---:|---:|
| marked \(b_0,b_1\) | 14 | \(1\) |
| new \(b\), marked \(c_0,c_1,c_2\) | 132 | \(1\) |
| new \(b,c\), marked \(d_0,d_1,d_2\) | 855 | \(3/2\) |
| new \(b,c,d\) | 1,905 | \(1\) |

The exponent is the better of two separately valid complete contractions:

$$
\delta_{\rm ass}
=\frac12\sum_K(\ell_K+k_K-1),
\qquad
\delta_{\rm proj}
=\frac12\sum_K(\ell_K+r_K-1).
\tag{3.3}
$$

Here \(\ell_K\) is cycle surplus, \(k_K\) is maximum physical-entry
occupancy, and \(r_K\) is maximum entry-respecting binary cut rank in a
connected component \(K\).  No Hilbert/projective regimes are mixed.
Consequently

$$
|\mathfrak C_{\rm A}|
\le C\sum_{v=9}^{12}(1+D)^vN^{-1}.
\tag{3.4}
$$

## 4. Direct reflected audit of Type B

Layer reversal is an automorphism of the four-layer Gaussian plant and maps
each legal Type-A transfer history to a legal Type-B history with the same
bounded local factors and reversed interpolation-time tuple.  In particular,
the forced outer \(\psi''\) factor maps from \(a_0\) to \(d_0\).

The implementation does not merely cite symmetry.  It reflects every edge
and every physical-entry block, checks that the result is exactly (1.2), and
then scores the retained left-moving expansion directly.

The first three branch families have the same exact case counts and minima:

| Type-B branch | exact scored cases | minimum decay exponent |
|---|---:|---:|
| marked \(c_0,c_1\) | 14 | \(1\) |
| new \(c\), marked \(b_0,b_1,b_2\) | 132 | \(1\) |
| new \(c,b\), marked \(a_0,a_1,a_2\) | 855 | \(3/2\) |

If the expansion also creates a new \(a_3\), the number of distinct
first-layer vertices changes from three to four.  The task coefficient has
reflection sign

$$
(-1)^{|V_A|}.
\tag{4.1}
$$

Indeed a degree-\(d\) outer vertex contributes \((-1)^d\) from its boundary
Hadamard factors and \((-1)^{d+1}\) from
\(\psi^{(d)}(-x)\), hence one minus sign per distinct outer vertex.  The
1,905 all-fresh physical placements are therefore reflection-even and
cancel before absolute values.  Thus

$$
|\mathfrak C_{\rm B}|
\le C\sum_{v=9}^{11}(1+D)^vN^{-1}.
\tag{4.2}
$$

## 5. The \(N^{1/16}\) transcript theorem

The complete earlier audit already gives the target \(N^{-v/16}\) decay for
all other level-nine and level-ten types and \(N^{-1}\) at levels eleven and
twelve.  The unique level-ten forest was separately repaired to \(N^{-1}\).
Levels four through eight retain the accepted \(N^{-1/2}\).  Combining those
facts with (3.4) and (4.2) gives

$$
\begin{aligned}
\Delta\le{}&
\sum_{v=4}^{8}C_v(1+D)^vN^{-1/2}
+\sum_{v=9}^{10}\widetilde C_v(1+D)^vN^{-v/16}\\
&+C_{\rm repair}\sum_{v=9}^{12}(1+D)^vN^{-1}
+\sum_{v=11}^{12}\widehat C_v(1+D)^vN^{-1}.
\end{aligned}
\tag{5.1}
$$

Set \(D=cN^{1/16}\).  The level-eight term and the already-passing
high-level rows are made small by choosing the absolute constant \(c\)
sufficiently small.  All repaired terms decay; their worst dimension power
is

$$
N^{12/16}N^{-1}=N^{-1/4}.
\tag{5.2}
$$

The accepted promise-conditioning loss is \(O(N^{-1})\).  Therefore

$$
D_{\mathsf P}^{\rm hard}=\Omega(N^{1/16}).
\tag{5.3}
$$

At this stage the generic low-level \(v=8\) row appeared limiting.  The later
complete-image audit proves that every sensitive level-eight type actually
has \(N^{-1}\) decay and repairs the level-seven/six types.

## 6. Program decision

The finite high-level obstruction list is closed.  The subsequent low-level
audit closes the complete terminal list and proves \(N^{1/12}\).  The later
sharpness audit matches the generic level-twelve \(N^{-1}\) graph bound.  The
next asymptotic gain therefore needs physical frame restrictions, better use
of signed coefficients before absolute values, or a different hard instance.

The realistic-size conclusion is unchanged.  The exponent improvement from
\(1/18\) to \(1/16\) raises the bare \(N=1024\) scale only from about \(1.47\)
to \(1.54\), far below six and before theorem constants.  Track A and the
alternative-instance/passive-counterprotocol parts of Track C therefore
remain first-class objectives, not follow-up exposition work.

Reproduction:

- `searches/level_nine_tree_centered_repair.py` enumerates all 200 Type-A
  histories, all seven dangerous partitions, every retained Type-A and
  Type-B repair branch, and the reflected all-fresh cancellation class.
- `tests/level_nine_tree_centered_repair.py` protects the two forced
  derivative sites, exact case counts, all decay minima, reflection map, and
  the \(1/16\) exponent.
