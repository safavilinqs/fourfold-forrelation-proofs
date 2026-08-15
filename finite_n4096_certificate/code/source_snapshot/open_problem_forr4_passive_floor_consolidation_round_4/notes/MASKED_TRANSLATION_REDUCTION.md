# Masked translation reduction

Date: 2026-07-17

Status: proved exact covariance, the arbitrary-law twirling reduction, and the full-group projective Clifford block formula. This is a valid cancellation mechanism for the original 354-entry dependency set. The current actual-mask registry repairs all 354 entries; the final twelve use the separate joint shared-quintic theorem rather than this translation reduction. The translation reduction was not needed to prove a universal coefficient-one lemma; it remains a valid structural explanation and diagnostic tool.

## Structural inventory

The original 354 affected entries reduce under support complementation and path
reversal to 97 disjoint structural orbits:

| class | entries | orbits |
|---|---:|---:|
| universal septimic | 96 | 24 |
| universal multicubic | 14 | 5 |
| universal double-cubic | 24 | 8 |
| universal noncubic | 124 | 36 |
| recovered universal | 96 | 24 |
| total | 354 | 97 |

This is the theorem-template inventory. Individual translation shapes within
a template are parameters of its reduced Fourier blocks, not new theorem
classes.

## Exact one-link covariance

Write $W_{x,y}=(-1)^{x\cdot y}$ for the order-$q$ Walsh matrix. Translate
the left-link coordinates by $(a,b)$ and the right-link coordinates by
$(c,d)$, where the first component shifts rows and the second shifts
columns. For left support $L$ and right support $R$, define

$$
\epsilon_{a,b,c,d}(L,R)
=
\prod_{(r,u)\in L}W_{r,c}W_{a,c}W_{u,d}
\prod_{(s,v)\in R}W_{a,s}W_{b,v}W_{b,d}.
$$

Then the exact signed-permutation moment satisfies

$$
M(T_{a,b}L,T_{c,d}R)
=\epsilon_{a,b,c,d}(L,R)M(L,R).
$$

To see this, write the hidden signed permutation as
(P_{\pi(u),u}=s_u), then relabel it by

$$
\pi'(u)=\pi(u\mathbin\oplus b)\mathbin\oplus c,
\qquad
s'_u=W_{a,\pi'(u)}W_{u,d}s_{u\mathbin\oplus b}.
$$

Conditional on $\pi'$, the $s'_u$ remain independent uniform signs. The
translated $KP$ and $PK$ entries differ from their primed versions by the
displayed deterministic Walsh factors. Multiplying over the supports proves
the identity.

The regression checks 400 arbitrary exact link moments at each of $q=4$
and $q=8$.

## Chain and mask separation

Apply the one-link identity to each of the three independent links. For four
block translations $g=(g_1,g_2,g_3,g_4)$, the chain moment acquires a
character

$$
\epsilon_g(U_1,U_2,U_3,U_4)
=\prod_{i=1}^3\epsilon_{g_i,g_{i+1}}(U_i,U_{i+1}).
$$

Each factor is multiplicative in its supports. On a physical occurrence
entry, $U_i=A_i\mathbin\dot\cup B_i$. Therefore

$$
\epsilon_g(U)=\epsilon_g(A)\epsilon_g(B).
$$

The distinctness mask is invariant under translation. Hence simultaneous
translation of the row and column configurations transforms the physical
masked matrix by row and column permutations and diagonal signs. This
statement includes every omitted cross-cut mask; it does not pass through an
unmasked completed kernel.

The exact regression checks 400 chain identities and 400 occurrence-sign
factorizations at each of $q=4$ and $q=8$.

## Arbitrary-law twirling theorem

For row and column laws $p,r$, put

$$
\Phi(p,r)=
\left\|\operatorname{diag}(p)^{1/2}
K
\operatorname{diag}(r)^{1/2}\right\|_1.
$$

The covariance above gives

$$
\Phi(gp,gr)=\Phi(p,r)
$$

for every translation $g$. The weighted trace norm is the root fidelity
between two positive operators and is jointly concave in $p,r$. Thus

$$
\Phi\left(\frac1{|G|}\sum_g gp,
          \frac1{|G|}\sum_g gr\right)
\geq
\frac1{|G|}\sum_g\Phi(gp,gr)
=\Phi(p,r).
$$

Therefore a worst arbitrary pair of diagonal laws can be taken translation
invariant. This is a global reduction for the actual masked operator, not an
invariant-law assumption.

On a fixed four-block occurrence-size sector, the signed translation action is
a projective representation. Its cocycle depends only on the four block sizes
and the two translations, not on the support shape. This support independence
is checked exactly at $q=4$ and $q=8$. Translation-invariant laws are
mixtures of translation-orbit shapes. A projective Fourier/Clifford transform
block-diagonalizes the group-lifted kernel, as in the earlier adjacent
cubic--quintic twisted-spectrum formula. The remaining proof problem is to
bound the resulting small matrices uniformly over the two orbit-shape
probability simplexes for each of the 97 templates.

### Exact projective-type classification

The projective cocycle has an exact finite classification. Let
$m=\log_2q$ and let $e_i=|A_i|\bmod 2$ be the four split parities. Its
commutator rank over $\mathbb F_2$ is $m\rho(e)$, where

$$
\rho(e)=
\begin{cases}
0,&e\in\{0000,1111\},\\
4,&e\in\{0011,1100\},\\
8,&\text{for every other balanced two-odd pattern}.
\end{cases}
$$

This follows by inserting the coordinate basis translations into the exact
cocycle formula and row-reducing its alternating commutator form. The rank
formula is checked exactly against all 97 templates at $q=4,8,64$; the
normalized ranks are identical at all three orders. The template counts are

| normalized rank | templates |
|---:|---:|
| 0 | 21 |
| 4 | 26 |
| 8 | 50 |

For every template, the row sector and its complementary column sector have
the same rank. Thus the remaining projective Fourier analysis has three
Clifford types, not 97 unrelated projective representations. The orbit-shape
geometry within each type still has to be bounded.

## Exact full-group Clifford formula

The row and complementary-column cocycles are in fact identical, not merely
equal-rank. This is checked entrywise on their translation-generator tables
for all 97 templates at $q=4,8,64$. Consequently every mixed orbit-shape
kernel uses one common twisted group algebra. Substitution in the explicit
translation-sign formula shows that $c$ is a bicharacter, so its generator
table determines it on the full group. Another 4,800 exact evaluations at
$q=2,4,8,64$ check that basis-table extension against the original cocycle
formula.

Let $G=(\mathbb F_2)^{8\log_2q}$. Fix a row orbit shape $A_i$ and a column
orbit shape $B_j$, and write $\alpha_{B_j}(t)$ for the exact translation sign
of $B_j$. Define the physical masked symbol

$$
f_{ij}(t)=\alpha_{B_j}(t)K(A_i,T_tB_j).
$$

After diagonal row and column sign gauges, the complete lifted block is

$$
K_{ij}(g,h)=c(g,g\mathbin\oplus h)
f_{ij}(g\mathbin\oplus h),
$$

where $c$ is the common projective cocycle. This identity is for the actual
masked kernel. It is checked directly against a complete physical $q=2$
matrix.

Let $r$ be the commutator rank of $c$. Symplectic row reduction gives
$2^{8\log_2q-r}$ irreducible characters, each with matrix dimension
$n=2^{r/2}$. For translation-invariant shape laws $p_i,r_j$, the exact
objective is

$$
\Phi(p,r)=\frac{1}{|G|}
\sum_{\lambda}n
\left\|
\left(\sqrt{p_i r_j}\,\widehat f_{ij}(\lambda)\right)_{i,j}
\right\|_1.
$$

This formula retains every cross-shape block. It is therefore the required
full-group reduction, rather than a pure-orbit or sampled-subgroup surrogate.
The implementation builds a symplectic basis, performs the radical and Pauli
Walsh transforms, and checks the projective multiplication law. Eighty scalar
regular-matrix comparisons, six mixed matrix-valued comparisons with gradient
checks, six multiplication-law inventories, and one complete physical $q=2$
comparison give 92 algebraic or synthetic checks and one physical check, in
addition to the 4,800 bicharacter comparisons.

The symbol construction also exposes a useful path factorization. Each block
translation has only $q^2$ local values, the three signed link moments form
three $q^2\times q^2$ edge tables, and the full $|G|=q^8$ symbol is their path
product. At $q=4$, this replaces 65,536 independent three-link moment calls by
three 256-entry link tables.

## Exact screens

The deterministic screen chooses three physical seeds per structural
template and a six-dimensional translation subgroup. It builds the complete
$64\times64$ group-lifted matrix, imposes all masks, and evaluates every
entry with exact signed-permutation moments. Floating point is used only for
the final singular values.

At $q=4$, all 97 templates were screened:

$$
\max \Phi_{\rm tested}=0.0220970869121.
$$

The worst tested template is

$$
((1,1,9,1),(0,1,5,0)).
$$

For every tested orbit matrix, the optimized value and its concavity tangent
agree within $2\times10^{-10}$. Thus the uniform group law is the numerical
global maximizer on each tested finite orbit, although this does not optimize
mixtures of different orbit shapes.

The eight highest-risk $q=4$ structural representatives were then evaluated
with exact $q=8$ moments. The maximum was

$$
6.01796939448\times10^{-5}.
$$

No tested law approaches coefficient one. The large decrease from $q=4$
to $q=8$ makes a masked translation/Fourier contraction more plausible and
makes further generic mask-factor estimates less attractive.

### Mixed orbit shapes

The twirling theorem does not reduce a law to one pure orbit; it permits a
mixture of orbit shapes. To test this missing interaction, the mixed screen
uses three physical row shapes and three physical column shapes under one
common translation subgroup, evaluates all nine cross-shape blocks exactly,
and optimizes the two three-point probability simplexes.

On the eight highest-risk structural representatives, the maxima are

$$
\begin{aligned}
q=4:&\quad 0.0194045575507,\\
q=8:&\quad 0.000279017854455.
\end{aligned}
$$

The $q=8$ mixed value is larger than the corresponding pure-orbit maximum
$6.01797\times10^{-5}$, so orbit-shape mixing is a real effect and must be
retained in the theorem. It remains more than three orders of magnitude below
one in this selected screen. The final concavity-tangent gaps are
$2.86\times10^{-5}$ at $q=4$ and $2.69\times10^{-12}$ at $q=8$.

### Complete $q=4$ translation group

The full-group Clifford formula was evaluated on one complete pure orbit for
all 97 templates. A further 30 pure shapes were screened in each of seven
leading templates, for 210 focused trials. The results are

| screen | maximum coefficient |
|---|---:|
| one canonical full orbit for all 97 templates | $0.0836986669474$ |
| 210 focused complete pure orbits | $0.176776695297$ |

The focused maximum occurs in

$$
((1,1,9,1),(0,1,5,0))
$$

and equals $1/(4\sqrt2)$ to numerical precision.

Finite mixed-shape simplexes containing the empirically leading shapes were
then evaluated with the exact full-group formula:

| normalized rank | shapes | coefficient | concavity tangent upper |
|---:|---:|---:|---:|
| 0 | 3 | $0.0223749015526$ | $0.0223813195039$ |
| 4 | 5 | $0.166351733007$ | $0.166741363473$ |
| 8 | 3 | $0.176776693147$ | $0.176776695297$ |

These tangents upper-bound the indicated finite shape simplexes, subject to
floating singular-value certification. They do not bound unlisted orbit
shapes. No complete-group finite screen approaches coefficient one, but the
larger values show that full translation and shape selection matter: the
earlier subgroup maxima substantially understated the strongest tested
orbits.

### Rejected record-one anchor shortcut

A tempting follow-up was to isolate the 118 affected entries whose compatible
record triple is forced to $(1,1,1)$, anchor the squared three-link path at a
singleton block, and charge operator norm one on every record-one link. This
does not work. The accepted bound

$$
{q^r\over(q)_r}=1\qquad(r=1)
$$

controls the undecorated odd-record core, not the full support matrix after
even decorations are restored.

The existing exact $q=4$ full-sector regression gives

| degrees | record | full operator norm | maximum row energy | maximum column energy |
|---:|---:|---:|---:|---:|
| $1,1$ | 1 | $1$ | $1$ | $1$ |
| $1,3$ | 1 | $\sqrt3$ | $3$ | $1$ |
| $3,1$ | 1 | $\sqrt3$ | $1$ | $3$ |
| $3,3$ | 1 | $4$ | $35/9$ | $35/9$ |

Thus leaf elimination cannot insert a unit bound on a decorated record-one
link. This rejects the shortcut, not coefficient one: when a cubic is
simultaneously compatible with singleton links on both axes, the restricted
operator returns to one. Any record-one repair must retain that two-sided
compatibility or use the correct endpoint fixed-slice energies. The regression
is inherited from
`../open_problem_forr4_passive_floor_consolidation_round_2/searches/signed_permutation_full_sector_spectra.py`
and runs inside `run_round4_checks.sh`.

## Decision and next proof target

No counterexample was found. The route now has a concrete cancellation
mechanism:

1. twirl arbitrary laws without loss by exact masked covariance;
2. use the proved full-group Clifford formula for the exact rank-$0/4/8$
   blocks; and
3. prove a common upper bound of one for the remaining shape-indexed symbol
   matrices across the 18 residual complement/reversal orbits.

The projective transform itself is no longer the missing calculation. The
next theorem must bound the shape-indexed matrices $\widehat f_{ij}(\lambda)$
uniformly over all physical orbit shapes in those residual orbits. A useful proof should exploit their
three-link path factorization, retain two-sided decoration compatibility, and
treat the rank-four family first, because it is the only tested type with
visible mixed-shape amplification. Proving only the listed finite simplexes
would not close the remaining 52-entry quarantine.

Reproduce the exact identities and finite screens with:

    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/masked_translation_covariance.py
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/masked_translation_cocycle_inventory.py
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/masked_translation_full_group_screen.py
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/masked_translation_subspace_screen.py
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/masked_translation_mixture_screen.py
