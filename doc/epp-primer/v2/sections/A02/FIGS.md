# Section A.2 · figure specifications

Section file: `v2/sections/A02-coupling.html` (id `a-coupling`, prefix `a2`). Four figures, in reading order:

| # | File | Id | Kind | Tag |
|---|---|---|---|---|
| A2.1 | `fig-overlap.html` | `viz-a2-overlap` | new | interactive |
| A2.2 | `fig-vertex.html` | `viz-a2-vertex` | new (replaces v1's inline Feynman diagrams `#a-coupling-fd`) | interactive |
| A2.3 | `fig-timeline.html` | `viz-a2-timeline` | port of `viz-a-coupling-timeline` | diagram |
| A2.4 | `fig-bardeen.html` | `viz-a2-bardeen` | port of `viz-a-coupling-bardeen` | interactive |

v1 source of both ports: `doc/epp-primer/sections/A02-coupling.html` (timeline: the `<figure id="viz-a-coupling-timeline">` block; Bardeen: the `<figure id="viz-a-coupling-bardeen">` block and the second half of the `<script>` at the end of the file, after the comment `From Bloch to Bardeen`).

General rules for all four (from `v2/build/BUILDER_GUIDE.md` §5): v2 notation only (no $N_p$, $\varepsilon_{n\bk}$ for band energies, $\epsilon$, $\Omega$ without subscript, $M_0$, $l_{\bq\nu}$; band energies are $E_{n\bk}$, the Fermi energy $E_{\mathrm F}$, the cell volume $\Omega_{\mathrm V}$ (macro `\OV`), the electronic dielectric function $\varepsilon^{\mathrm e}$); vectors bold (`drawRich` with `{bold:true}`), labels upright (`\mathrm`), the imaginary unit and $\mathrm e$ upright; HTML readouts in MathJax `$…$`; colours read at draw time with `EPP.color(...)`; electrons in `electron`, phonons in `phonon`; caption in the four parts *What to notice / Marks / Try / Model*, about 200–260 words; equation references with `<a class="eqref" href="#a2-eq-…"></a>`; no source equation numbers anywhere in the figure. CSS classes prefixed `a2-<name>-`.

Equation ids of the section you will cite (tags as of this writing; always link by id):
`a2-eq-bloch` (A2.3), `a2-eq-rigid` (A2.8), `a2-eq-wave` (A2.10), `a2-eq-bvq` (A2.11), `a2-eq-select` (A2.12), `a2-eq-dV` (A2.14), `a2-eq-H` (A2.15), `a2-eq-g` (A2.16), `a2-eq-anatomy` (A2.17), `a2-eq-H0` (A2.19), `a2-eq-dw` (A2.20), `a2-eq-free` (A2.22), `a2-eq-bloch1929` (A2.23), `a2-eq-tf` (A2.24), `a2-eq-lindhard` (A2.25), `a2-eq-bardeen` (A2.26), `a2-eq-limit` (A2.27), `a2-eq-pines` (A2.28), `a2-eq-def` (A2.29). A.1: `a1-eq-u` (A1.27).

---

## A2.1 · `fig-overlap.html` · `viz-a2-overlap` · new

**Title.** "Why $\bk$ goes to $\bk+\bq$: a frozen phonon and the overlap that selects the final state" (or shorter; sentence case) `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. A frozen lattice wave of wavevector $q$ produces a potential change $\Delta V(x)=\xi\,\partial_qV(x)$ that is itself a wave: periodic up to the phase $\mathrm e^{\mathrm iqa}$ per cell, (A2.11).
2. The matrix element is an overlap: the integrand $\psi^*_{k'}(x)\,\partial_qV(x)\,\psi_k(x)$ repeats from cell to cell, turned by the phase $\mathrm e^{\mathrm i(k+q-k')a}$. Summed over the $N_\ell$ cells of the Born–von Kármán ring, the cell contributions cancel unless $k'=k+q$ (modulo $2\pi/a$): the selection rule (A2.12).
3. On the band, the allowed transition is a horizontal shift by $q$; when $k+q$ leaves the zone it folds back by $G=2\pi/a$ (umklapp), and the overlap is nonzero for the folded state, which is the same state.
4. The surviving matrix element, times the zero-point amplitude, gives $|g|$ of tens to about 150 meV for this model, matching the estimate "zero-point amplitude × potential gradient" (A2.17).

**Include position.** Section `a2-wave` ("A wave of displacements is a wave of potential"), right after the paragraph that ends: "In a chain with one atom per cell that moves only along the chain, the labels $\kappa$ and $\alpha$ drop out, and (A2.10) reads $\Delta V(x)=\xi\,\partial_qV(x)$, with $x$ the position along the chain. Figure A2.1 shows such a frozen wave, the potential change it creates, and what happens when that change is sandwiched between two Bloch states." The selection-rule derivation (A2.12) follows the figure.

**Physics (v2 notation, all in this section).**
- Chain of Figure A1.3 (A.1 `fig-waves`): one atom per cell, mass $M=12.01$ amu, spacing $a=1.42$ Å, nearest-neighbour spring $K=20$ eV/Å², motion along the chain. Phonon energy $\hbar\omega_q=2\cdot64.654\,\text{meV}\,\sqrt{K/M}\,|\sin(qa/2)|$, maximum 166.9 meV at $q=\pi/a$. $N_\ell=12$ cells on a Born–von Kármán ring; allowed wavevectors $k,q=2\pi s/(N_\ell a)$, integer $s$ from $-5$ to $6$.
- Ion potential (rigid ion, (A2.8)): a Gaussian well on each atom, $V_\kappa(x)=-V_{\mathrm a}\exp(-x^2/2w^2)$ with depth $V_{\mathrm a}=10$ eV and rms width $w=0.40$ Å (model values; in the caption write them as numbers, without new symbols). The potential change per unit amplitude of the wave is $\partial_qV(x)=\sum_\ell\mathrm e^{\mathrm iq\ell a}\,[-V_\kappa'(x-\ell a)]$, the 1D form of (A2.10) and (A2.14); use the periodic distance on the ring.
- Electron states: tight-binding Bloch sums on the ring, $\psi_k(x)\propto\sum_\ell\mathrm e^{\mathrm ik\ell a}\,\varphi(x-\ell a)$ with a Gaussian orbital $\varphi(x)=\exp(-x^2/2s^2)$, $s=0.50$ Å (do not name the orbital in the caption; $\varphi$ is only for this spec; never call it $\phi$, which is the phonon field), normalized numerically on the ring. One band. Band energies for the band panel: $E_k=-(5.0\ \text{eV})\cos ka$ (tight-binding band of width 10 eV, model; not computed from the orbitals, say so in *Model*).
- Matrix element $M(k',k)=\int_{\text{ring}}\psi^*_{k'}\,\partial_qV\,\psi_k\,\mathrm dx$ (eV/Å), per-cell contributions $c_\ell=\int_{\text{cell }\ell}\dots$; exactly $c_\ell=c_0\,\mathrm e^{\mathrm i(k+q-k')\ell a}$. Coupling $|g|=(\hbar/2M\omega_q)^{1/2}\,|M(k+q,k)|$, with $\hbar/2M\omega=2.0901\ \text{Å}^2/(M[\text{amu}]\,\hbar\omega[\text{meV}])$; result in meV.
- Grid: 80 points per cell (960 points) is enough; everything is a few ms. Prototype values (this spec's author checked them in Python): $s_k=2$, $s_q=3$: $\hbar\omega_q=118.0$ meV, zero-point amplitude 0.0384 Å, $|M|=1.52$ eV/Å, $|g|=58$ meV; $s_k=0$, $s_q=3$: $|g|=144$ meV; $s_k=3$, $s_q=6$: $|g|=132$ meV; $s_k=0$, $s_q=6$: $|g|=0$ by symmetry. Mismatched $k'$: $|M|\sim10^{-15}$ (show it as 0).
- $q=0$ ($s_q=0$) is a rigid translation: zero frequency, infinite zero-point amplitude, left out of all sums. Either exclude $s_q=0$ from the slider or show a note "rigid translation, left out" and blank the $g$ readout.

**What to draw.** Stack vertically (phone: one column, controls below the plots):
- (a) *Frozen wave and potential.* The chain over the 12 cells (unrolled ring, cell boundaries as faint ticks): rest positions (small grey dots), displaced atoms (filled, `phonon` colour), displacement $\mathrm{Re}\,u_\ell=\mathrm{Re}[\xi\,\mathrm e^{\mathrm iq\ell a}]$ with $\xi$ the zero-point amplitude, exaggerated (×12, state it). Behind: the static periodic potential $\sum_\ell V_\kappa(x-\ell a)$ (thin, `muted`), and the potential change $\mathrm{Re}\,\Delta V(x)$ (thick, `phonon` colour) on its own scale; label $\Delta V(x)=\xi\,\partial_qV(x)$. Axis: $x$ (Å), energy (eV).
- (b) *The three factors and their product.* $\mathrm{Re}\,\psi_k(x)$ (`electron`), $\mathrm{Re}\,\psi^*_{k'}(x)$ (`electron`, dashed), and the integrand $\mathrm{Re}[\psi^*_{k'}\,\partial_qV\,\psi_k]$ as a filled curve (`c5` or ink) over the same 12 cells; shade every cell with the sign of its integral, or mark each cell with a small arrow for $c_\ell$.
- (c) *Sum over cells.* A phasor dial in the style of A.1's `fig-bz` row-sum dials: the 12 complex numbers $c_\ell$ as arrows from the origin (or tip to tail), and their sum as a thick arrow. Matched: all 12 aligned, sum $=12\,c_0$. Mismatched: spread evenly, sum $=0$. Next to it, a small bar chart of $|M(k',k)|$ against all 12 values of $k'$: a single nonzero bar at $k'=k+q$ (folded), the selected $k'$ highlighted.
- (d) *Band.* $E_k=-(5.0\ \text{eV})\cos ka$ over $-\pi/a<k\le\pi/a$, the 12 allowed $k$ as hollow dots; initial $k$ filled (`electron`), the target $k+q$ as a ring, the selected $k'$ as a square. A curved arrow from $k$ to $k+q$ labelled $\mathbf q$ (bold in `drawRich`); if $k+q$ leaves the zone, continue the arrow into a faded repeat of the band and draw a dashed fold-back arrow labelled $-\mathbf G$ (here $G=\pm2\pi/a$), with the text "umklapp". Axis: $k$ in units of $\pi/a$, $E$ (eV).

**Controls (4).**
1. $q$: slider over integers $s_q\in\{-5,\dots,6\}$ ($q=2\pi s_q/12a$), shown as $q$ in units of $\pi/a$ and in Å⁻¹.
2. $k$: slider $s_k\in\{-5,\dots,6\}$.
3. $k'$: slider $s_{k'}\in\{-5,\dots,6\}$; clicking a dot in (d) or a bar in (c) also sets $k'$ (keyboard: the slider).
4. Check box "$k'$ follows $k+q$" (default on; moving the $k'$ slider turns it off).
Optional: Play/Pause (`EPP.anim`) that advances the phase of the frozen wave in (a) only, for intuition; not required.

**Readouts (HTML, MathJax).** $\hbar\omega_q$ (meV); zero-point amplitude $(\hbar/2M\omega_q)^{1/2}$ (Å); phase step $(k+q-k')a$ in units of $\pi$, reduced to $(-1,1]$; $|\sum_\ell c_\ell|/(12|c_0|)$ (1 or 0); $|\langle\psi_{k'}|\partial_qV|\psi_k\rangle|$ (eV/Å); $|g|$ (meV) when $k'=k+q$, else "0 (wavevectors do not match)"; "normal" or "umklapp, $G=\pm2\pi/a$".

**Initial state.** $s_k=2$, $s_q=3$, "$k'$ follows $k+q$" on, so $k'=k+q$ ($s=5$): aligned arrows, one bar, $|g|=58$ meV, normal process. The bar chart shows the selection rule even before any interaction.

**Caption key points.**
- *What to notice.* The integrand repeats from cell to cell, turned by $(k+q-k')a$; the 12 cell contributions add up only for $k'=k+q$ and cancel for every other $k'$: one bar. When $k+q$ leaves the zone, the nonzero bar sits at the folded wavevector.
- *Marks.* Colours (electron states, phonon displacement and $\Delta V$), exaggeration factor, ring = target $k+q$, square = selected $k'$, dashed fold-back arrow = $-\mathbf G$, dial arrows = $c_\ell$, thick arrow = their sum.
- *Try.* Untick "follows" and step $k'$ away from $k+q$: the arrows fan out and cancel. Set $s_k=4$, $s_q=4$: $k+q$ folds to $s=-4$ (umklapp) and the bar follows. Compare $s_q=1$ and $s_q=6$: the zero-point amplitude shrinks as $\omega_q$ grows, while the gradient changes with $k$.
- *Model.* Chain of Figure A1.3 ($M$, $a$, $K$ as above, carbon-like, illustrative); rigid-ion Gaussian wells (10 eV, 0.40 Å) and Gaussian orbitals (0.50 Å), model values; band $E_k=-(5.0\ \text{eV})\cos ka$ drawn for orientation only. Equations: (A2.10), (A2.11), (A2.12), (A2.14), (A2.16) and (A2.17), with $\xi$ the zero-point amplitude of (A1.27). The matrix element does not require energy conservation; real transitions also conserve energy (Sections A.7, A.9). Phonon energies (≤ 167 meV) are invisible on the eV scale of the band.

---

## A2.2 · `fig-vertex.html` · `viz-a2-vertex` · new

**Title.** "The vertices of the Hamiltonian and crystal momentum on the Brillouin zone" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. The first-order term of (A2.19) is a vertex: electron line $n\mathbf k$ in, $m\,\mathbf k+\mathbf q$ out, one phonon line; the same term absorbs a phonon $\mathbf q\nu$ or emits a phonon $-\mathbf q\nu$, and in both cases the electron gains $\mathbf q$.
2. On the two-dimensional Brillouin zone, $\mathbf k+\mathbf q$ is vector addition; if the sum leaves the first zone, it is brought back by a reciprocal lattice vector $\mathbf G$ (umklapp). Crystal momentum is conserved up to $\hbar\mathbf G$.
3. The second-order (Debye–Waller) term has two phonon lines at one point; contracting them ($\mathbf q'=-\mathbf q$, $\nu'=\nu$) gives a loop that shifts the band energy by an amount $\propto 2n_{\mathbf q\nu}+1$, nonzero at $T=0$, (A2.20).

**Include position.** Section `a2-hamiltonian` ("The whole Hamiltonian, term by term"), after the opening paragraph and before the display of (A2.19). Sentences before: "Collecting the independent electrons (A2.6), the phonons of Section A.1, the first-order coupling (A2.15) and a second-order coupling, explained below, gives the review's standard form of the electron–phonon Hamiltonian, the object the rest of Part A works with. Section A.1 called the phonon Hamiltonian $\hat H$; from here on it is $\hat H^{\mathrm{ph}}$, and $\hat H$ denotes the whole. Figure A2.2 draws the coupling terms as diagrams and shows the bookkeeping of crystal momentum at the first-order vertex on the Brillouin zone of the hexagonal sheet of Section A.1." It is cited again in Check yourself, question 2 (umklapp).

**Physics.**
- Hexagonal sheet of A.1 (hBN-like): $a=2.504$ Å, $\mathbf b_1=(2\pi/a)(1,-1/\sqrt3)$, $\mathbf b_2=(2\pi/a)(0,2/\sqrt3)$ (as in A.1 `fig-bz`), $|\mathbf b_j|=2.898$ Å⁻¹. First Brillouin zone: hexagon with corners $\mathrm K$ at $4\pi/3a=1.673$ Å⁻¹ on the $x$ axis and at multiples of 60°, edge midpoints $\mathrm M$ at $2\pi/(\sqrt3a)=1.449$ Å⁻¹.
- Kinematics only: $\mathbf k$ and $\mathbf q$ inside the first zone; $\mathbf k'=\mathbf k+\mathbf q$; if $\mathbf k'$ lies outside, find the $\mathbf G$ among the 6 shortest (and, for safety, the next shell) that brings $\mathbf k'-\mathbf G$ into the zone (closest to $\Gamma$). Normal process if $\mathbf G=0$.
- No band structure and no values of $g$ are needed. Labels follow (A2.15) and (A2.19): amplitude $g_{mn\nu}(\mathbf k,\mathbf q)/\sqrt{N_\ell}$; operator strings $c^\dagger_{m\mathbf k+\mathbf q}c_{n\mathbf k}\hat a_{\mathbf q\nu}$ (absorption) and $c^\dagger_{m\mathbf k+\mathbf q}c_{n\mathbf k}\hat a^\dagger_{-\mathbf q\nu}$ (emission), from $\hat\phi_{\mathbf q\nu}=\hat a_{\mathbf q\nu}+\hat a^\dagger_{-\mathbf q\nu}$.
- Debye–Waller panel: vertex with lines $\mathbf q\nu$ and $\mathbf q'\nu'$, outgoing electron $m\,\mathbf k+\mathbf q+\mathbf q'$, amplitude $g^{\mathrm{DW}}_{mn\nu\nu'}(\mathbf k,\mathbf q,\mathbf q')/N_\ell$; contracted loop with label $2n_{\mathbf q\nu}+1$ and the shift of (A2.20).

**What to draw.**
- (a) *Brillouin zone* (canvas, square, the larger panel): the hexagon, $\Gamma$ and the six nearest reciprocal lattice points as small dots with the neighbouring zones as faint hexagons (so a folded vector visibly lands in an equivalent zone), labels $\Gamma$, $\mathrm K$, $\mathrm M$ upright. Arrow $\mathbf k$ from $\Gamma$ (`electron`), arrow $\mathbf q$ from the tip of $\mathbf k$ (`phonon`), arrow $\mathbf k+\mathbf q$ from $\Gamma$ (`electron`, dashed). For umklapp: the vector $-\mathbf G$ from the tip of $\mathbf k+\mathbf q$ to its image (ink, dotted) and the folded arrow $\mathbf k+\mathbf q-\mathbf G$ from $\Gamma$ (`electron`, solid, lighter). Draggable tips (filled handles) for $\mathbf k$ and $\mathbf q$; `touch-action: none` on this canvas only. Bold vector labels via `drawRich` `{bold:true}`.
- (b) *First-order vertex* (SVG via `EPP.fd.draw`): electron line in ($n\mathbf k$), vertex dot, electron line out ($m\,\mathbf k+\mathbf q$; when umklapp, add "$=m\,\mathbf k+\mathbf q-\mathbf G$" below), phonon (wavy) line coming in labelled $\mathbf q\nu$ for absorption, or going out labelled $-\mathbf q\nu$ for emission. Under it, in HTML/MathJax: the amplitude $g_{mn\nu}(\mathbf k,\mathbf q)/\sqrt{N_\ell}$ and the operator string. If `EPP.fd` labels cannot be bold, put the labels in MathJax below or beside the SVG.
- (c) *Second order* (two small static SVG diagrams side by side): the Debye–Waller vertex with two phonon lines $\mathbf q\nu$, $\mathbf q'\nu'$ meeting the electron at one point (in $n\mathbf k$, out $m\,\mathbf k+\mathbf q+\mathbf q'$); and the contracted loop ($\mathbf q'=-\mathbf q$, electron in and out $n\mathbf k$) with the label $2n_{\mathbf q\nu}+1$.

**Controls (4).**
1. Drag the tip of $\mathbf k$ (pointer); keyboard: focus (a), arrow keys move the active tip by 0.05 Å⁻¹, Tab or space switches the active tip between $\mathbf k$ and $\mathbf q$.
2. Drag the tip of $\mathbf q$ (constrained to the first zone).
3. Segmented control: "absorb $\mathbf q\nu$" / "emit $-\mathbf q\nu$" (changes (b) and the operator string).
4. Button "Reset".
`aria-label`s on the canvas ("Brillouin zone with electron wavevector k, phonon wavevector q and their sum") and on each SVG.

**Readouts.** $\mathbf k$, $\mathbf q$, $\mathbf k+\mathbf q$, $\mathbf G$ as $(x,y)$ in Å⁻¹ and their lengths; "normal process ($\mathbf G=0$)" or "umklapp, $\mathbf G=\dots$"; the final-state label inside the zone.

**Initial state.** $\mathbf k=(1.339,0)$ Å⁻¹ ($0.8\,\Gamma\mathrm K$), $\mathbf q=(0.55,0.45)$ Å⁻¹, absorption: $\mathbf k+\mathbf q=(1.889,0.45)$ lies outside; $\mathbf G=\mathbf b_1+\mathbf b_2=(2.509,1.449)$ Å⁻¹ folds it to $(-0.620,-0.999)$ Å⁻¹: umklapp visible at once.

**Caption key points.**
- *What to notice.* The final state is always $\mathbf k+\mathbf q$; when that sum leaves the hexagon it is the same state as the folded vector $\mathbf k+\mathbf q-\mathbf G$ inside it (umklapp). Absorbing $\mathbf q$ and emitting $-\mathbf q$ give the electron the same crystal momentum.
- *Marks.* Colours and line styles; faint hexagons = equivalent zones; dotted arrow = $-\mathbf G$; wavy line = phonon; dot = vertex.
- *Try.* Drag $\mathbf q$ shorter until the process becomes normal; drag $\mathbf k$ toward $\Gamma$; switch absorb/emit and compare the phonon line and the operator string.
- *Model.* Zone of the hexagonal sheet of Section A.1 ($a=2.504$ Å); vertices of (A2.15) and (A2.19), with $\hat\phi_{\mathbf q\nu}=\hat a_{\mathbf q\nu}+\hat a^\dagger_{-\mathbf q\nu}$; the loop is the shift (A2.20), preview of the Debye–Waller self-energy (Section A.7). Energy conservation is not drawn: the Hamiltonian contains every vertex.

**v1 material.** v1 drew three inline diagrams with `EPP.fd.draw` (`#a-coupling-fd1`–`fd3`, script at the end of the v1 section) with labels "nk", "mk+q", "qν", "g^{DW}" in plain text. Reuse their geometry for (b) and (c), but with v2 labels (bold vectors), the emission reading added, and the kinematics panel (a) new.

---

## A2.3 · `fig-timeline.html` · `viz-a2-timeline` · port of `viz-a-coupling-timeline`

**Title.** "What replaces the potential gradient: from Bloch to first principles and Part B" `<span class="tag">diagram</span>`.

**What the reader must learn.** For half a century progress meant a better guess for one effective potential in Bloch's formula (A2.23); from 1985 the potential is computed self-consistently, and the open problems become numerical (dense grids) and electrostatic (the long-range polar part). Several steps (Fröhlich, Bardeen–Pines, the polar long-range vertex) are the ancestors of the EPP formalism of Part B.

**Include position.** Start of section `a2-history` ("A short history of $g$"), after the paragraph that ends: "The decades that followed replaced the potential gradient in Bloch's formula by better and better guesses. Each of those guesses is a special case of (A2.16): choose a potential, choose the electron states, and evaluate. Figure A2.3 gives an overview."

**Content: at most 12 cards, kept short** (year · names · one formula · at most 8 words of note · optional link). Formulas in v2 notation, as in the section:
1. 1929 · Bloch · $V_0\approx1.70$ eV, constant · deformable continuum; first theory of resistivity (model)
2. 1931 · Nordheim · $V_\kappa(\mathbf q)/\Omega_{\mathrm V}\propto q^{-2}$ · bare rigid ion; far too strong (model)
3. 1936 · Mott and Jones · ion cut at the Wigner–Seitz sphere · screening by truncation (model)
4. 1937 · Bardeen · $V_\kappa(\mathbf q)/\Omega_{\mathrm V}\varepsilon^{\mathrm e}(q)$ · self-consistent linear screening; link `#a2-eq-bardeen` (model)
5. 1950 · Bardeen and Shockley · $E_{\mathrm{def}}=\Omega_{\mathrm V}\,\partial E_{n\mathbf k}/\partial\Omega_{\mathrm V}$ · deformation potentials fitted to mobilities (model)
6. 1950 · Fröhlich, Pelzer and Zienau · $g\propto1/q$ · polar crystals; Section A.4 (model; chip "Part B precursor" → `#b-frohlich`)
7. 1955 · Bardeen and Pines · $\varepsilon^{\mathrm e}(q,\omega_{\mathbf q\nu})$ · first retarded vertex (model; chip → `#b-frohlich`)
8. 1959–61 · Phillips–Kleinman, Cohen–Heine, Sham · pseudopotential $g$ · improved rigid ion (pseudopotentials)
9. 1985 · Dacorogna, Cohen and Lam · $\langle\psi|\Delta_{\mathbf q\nu}V^{\mathrm{KS}}|\psi\rangle$ from frozen phonons · aluminium from first principles (first principles)
10. 1987–92 · density-functional perturbation theory · any $\mathbf q$ in one cell · Section A.3 (first principles)
11. 2007–15 · Wannier interpolation; polar long-range part · dense grids; Section A.4 (first principles; chip → `#b-foundations`)
12. draft · EPP manuscript · photon-exchange vertex, dynamic screening · Section B.8, `#b-frohlich` (electrodynamic)

Four categories with a colour bar on top of each card and a legend strip: model potentials (`c4`), pseudopotentials (`muted`), first principles (`c5`), electrodynamic, Part B (`photon`). Cards in a responsive grid (v1: `repeat(auto-fill, minmax(148px,1fr))`), one column of two or three cards per row on a phone.

**Controls.** None required (diagram). Optional: the legend items as four toggle buttons (keyboard accessible, `aria-pressed`) that dim the cards of the other categories.

**Readouts.** None.

**Initial state.** All cards shown.

**Caption key points.**
- *What to notice.* Until 1985 each step was a better guess for one effective potential; afterwards the potential is computed, and what remains hard is numerical (dense grids) and electrostatic (the long-range polar part). The chips mark the steps that the EPP formalism generalizes.
- *Marks.* Colour bar = category; chip = Part B precursor with a link.
- *Try.* (Only if the optional toggles exist.) Show only "first principles".
- *Model.* Dates and attributions as told by the review (no equation numbers of the source); formulas in the notation of this section with links to (A2.23), (A2.26), (A2.28), (A2.29), and to Sections A.3, A.4, B.8.

**What in v1 must change.**
- v1 uses source notation in the cards: $\epsilon(q)$, $\epsilon(q,\omega_{\mathbf q\nu})$, $E_1=\Omega\,\partial\varepsilon/\partial\Omega$, $\langle u|\Delta_{\mathbf q\nu}v^{\rm KS}|u\rangle$ "(G.11)", $\partial v^{\rm KS}=(\epsilon^{{\rm H}xc})^{-1}\partial v^{\rm en}$, "long-range $g^{\mathcal L}$ (G.180)", and for EPP "$g^{\rm dyn}=iJ^*_\sigma D^{\rm e}J_{nm}$ (cf. E.416)". Replace by the v2 forms above; no source equation numbers.
- v1's caption begins "Sources: Giustino §II.A–C, §VI.A." and has only *What to notice*; rewrite in the four-part structure.
- v1 styles everything inline; use classes `a2-timeline-…` and theme colours read from CSS variables.
- Shorten: v1 had 13 cards with longer notes; merge Giustino–Cohen–Louie (2007) and Verdi–Giustino/Sjakste (2015) into one card as above.
- The "EPP precursor" badges become "Part B precursor" chips; the last card is the EPP manuscript itself, linked to Section B.8.

---

## A2.4 · `fig-bardeen.html` · `viz-a2-bardeen` · port of `viz-a-coupling-bardeen`

**Title.** "From Bloch to Bardeen: screening the electron–ion vertex" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. The bare ion (Nordheim) diverges as $q^{-2}$ and gives $|g|\propto q^{-3/2}$ for acoustic phonons; any finite long-wavelength potential (Bloch's constant, Thomas–Fermi, Lindhard) gives $|g|\propto q^{1/2}$.
2. Thomas–Fermi and Lindhard screening both tend to $-\frac23E_{\mathrm F}$ as $q\to0$, (A2.27); the Mott–Jones cut stays five to six times deeper.
3. Lindhard screens less than Thomas–Fermi at short wavelengths, with a weak kink at $q=2k_{\mathrm F}$ (seed of the Kohn anomaly, Section A.6).
4. Bloch's constant 1.70 eV has the right size for the alkali metals; for Al the screened $|g|$ at $q=k_{\mathrm F}/2$ is about 270 meV, for Na about 90 meV.
5. Dynamic (Bardeen–Pines) screening changes the vertex only by a phase of about $3\times10^{-3}$ for acoustic phonons.

**Include position.** Section `a2-bardeen`, after the paragraph with the long-wavelength limit (A2.27) and the numbers, which ends: "… about 270 meV in aluminium and 90 meV in sodium at $q=k_{\mathrm F}/2$, with longitudinal sound velocities of 6400 and 3200 m/s (illustrative), against 1.7 and 1.0 eV for the bare ions. Figure A2.4 compares the screening models, from Bloch's constant to Lindhard's function, for these three metals." It is also cited earlier, in the Nordheim paragraph of `a2-free` ("about 1.7 eV in aluminium at $q=0.88$ Å⁻¹ (Figure A2.4 plots it)"). The intuition callout "why a metal removes the divergence" follows the figure.

**Physics (keep v1's code; it is correct). All in v2 notation:**
- Jellium with $Z_\kappa$ electrons per ion (label it $Z$ in the controls, "valence $Z$") at density parameter $r_s$: $k_{\mathrm F}a_0=1.91916/r_s$, $E_{\mathrm F}=(k_{\mathrm F}a_0)^2\times13.6057$ eV, $k_{\mathrm{TF}}^2=4k_{\mathrm F}/\pi a_0$, $\Omega_{\mathrm V}=Z\,\frac{4\pi}3(r_sa_0)^3$, $R_{\mathrm{WS}}=Z^{1/3}r_sa_0$, $q_{\mathrm D}=(6\pi^2/\Omega_{\mathrm V})^{1/3}$, $v_{\mathrm F}=(k_{\mathrm F}a_0)\times2.18769\times10^6$ m/s, $a_0=0.529177$ Å.
- Left panel, the potential per cell $V_0(q)$ in eV (the quantity that replaces $V_0$ in (A2.23)): Bloch $-1.70$ eV (drawn negative for comparison; only $|V_0|$ enters $|g|$; say so in *Marks*); Nordheim $V_\kappa(q)/\Omega_{\mathrm V}=-4\pi Z\cdot14.3996\ \text{eV Å}/(\Omega_{\mathrm V}q^2)$; Mott–Jones $\times[1-\cos(qR_{\mathrm{WS}})]$; Thomas–Fermi $/[1+(k_{\mathrm{TF}}/q)^2]$, (A2.24); Lindhard $/[1+(k_{\mathrm{TF}}/q)^2F_{\mathrm L}(q/2k_{\mathrm F})]$, (A2.25), with $F_{\mathrm L}(x)=\frac12+\frac{1-x^2}{4x}\ln|\frac{1+x}{1-x}|$ (guards: $F_{\mathrm L}\to1-x^2/3$ for small $x$, $F_{\mathrm L}(1)=\frac12$). Horizontal line at $-\frac23E_{\mathrm F}$, (A2.27); vertical line at $2k_{\mathrm F}$. The $y$ axis scales with $E_{\mathrm F}$ ($[-1.6E_{\mathrm F},0.14E_{\mathrm F}]$), $x=q/k_{\mathrm F}\in[0,3]$.
- Right panel, log–log, $|g|=(\hbar/2M\omega_q)^{1/2}\,q\,|V_0(q)|$ in meV for a longitudinal Debye phonon $\omega_q=v_{\mathrm s}q$, (A2.23) and (A2.26), with $\hbar\omega[\text{meV}]=6.5821\times10^{-3}\,v_{\mathrm s}[\text{m/s}]\,q[\text{Å}^{-1}]$ and $\hbar/2M\omega=2.0901\ \text{Å}^2/(M[\text{amu}]\,\hbar\omega[\text{meV}])$; no $N_\ell$ (the v2 convention: $g$ has none). Shade $q>q_{\mathrm D}$ (Debye model unreliable). Slope labels $\propto q^{-3/2}$ (bare) and $\propto q^{1/2}$ (finite $V_0$). $x\in[0.02,3]$, $y\in[1,10^5]$ meV.
- Readout of Bardeen–Pines, (A2.28): dynamic Lindhard at $q=k_{\mathrm F}$ ($z=q/2k_{\mathrm F}=\frac12$, $u=v_{\mathrm s}/v_{\mathrm F}$), $\mathrm{Re}\,\varepsilon^{\mathrm e}=1+(k_{\mathrm{TF}}/q)^2[\frac12+L(z-u)+L(z+u)]$ with $L(y)=\frac{1-y^2}{8z}\ln|\frac{y+1}{y-1}|$, $\mathrm{Im}\,\varepsilon^{\mathrm e}=(k_{\mathrm{TF}}/q)^2\frac{\pi}{2}u$ for $z+u<1$ (and $\frac{\pi}{8z}[1-(z-u)^2]$ for $|z-u|<1<z+u$); show $\mathrm{Im}/\mathrm{Re}$.
- Check values (this spec's author recomputed them): Al ($r_s=2.07$, $Z=3$): $k_{\mathrm F}=1.752$, $k_{\mathrm{TF}}=2.053$ Å⁻¹, $E_{\mathrm F}=11.70$ eV, $\Omega_{\mathrm V}=16.52$ Å³, $-\frac23E_{\mathrm F}=-7.80$ eV, Mott–Jones $-41.0$ eV, $q_{\mathrm D}/k_{\mathrm F}=0.874$, $v_{\mathrm s}/v_{\mathrm F}=3.2\times10^{-3}$, $\mathrm{Im}/\mathrm{Re}=3.0\times10^{-3}$, $|g|(k_{\mathrm F}/2)=270$ meV (Lindhard) and 1719 meV (bare). Na ($3.93$, 1): $E_{\mathrm F}=3.24$ eV, $-2.16$ eV, Mott–Jones $-10.4$ eV, 90 meV and 1007 meV. K ($4.86$, 1): $E_{\mathrm F}=2.12$ eV, $-1.41$ eV, Mott–Jones $-8.40$ eV, 52 meV.

**What to draw.** Two plot panels side by side on desktop, stacked on a phone (`viz-grid two`), controls below, readouts below the controls. Five series with v1's line styles: Bloch (`muted`, dotted), Nordheim (`c5`), Mott–Jones (ink, dashed), Thomas–Fermi (`c4`, dashed), Lindhard (`electron`, thick). Axis labels in v2 notation: $q/k_{\mathrm F}$; "$V_0(q)$ (eV)"; "$|g(q)|$ (meV)".

**Controls (5, as in v1).** Metal preset (Al, Na, K; segmented); $r_s$ 1.5–6 (step 0.01, units of $a_0$); valence $Z$ (1, 2, 3, 4); ion mass $M$ 5–210 amu; sound velocity $v_{\mathrm s}$ 1000–8000 m/s (step 50). Moving a slider clears the preset (v1 logic). Presets: Al ($r_s=2.07$, $Z=3$, $M=26.98$ amu, $v_{\mathrm s}=6400$ m/s), Na (3.93, 1, 22.99, 3200), K (4.86, 1, 39.10, 2000) (literature values, illustrative).

**Readouts (MathJax).** $k_{\mathrm F}$, $k_{\mathrm{TF}}$ (Å⁻¹); $E_{\mathrm F}$ (eV); $\Omega_{\mathrm V}$ (Å³); $-\frac23E_{\mathrm F}$ (eV); Mott–Jones limit $-\frac32Ze^2/R_{\mathrm{WS}}$ (eV); $q_{\mathrm D}/k_{\mathrm F}$; $v_{\mathrm s}/v_{\mathrm F}$; $\mathrm{Im}\,\varepsilon^{\mathrm e}/\mathrm{Re}\,\varepsilon^{\mathrm e}$ at $q=k_{\mathrm F}$ (Bardeen–Pines); optionally $|g|$ at $q=k_{\mathrm F}/2$ for Lindhard and for the bare ion.

**Initial state.** Al preset.

**Caption key points.**
- *What to notice.* The bare ion gives $|g|\propto q^{-3/2}$; every finite $V_0(0)$ gives $|g|\propto q^{1/2}$. Thomas–Fermi and Lindhard both reach $-\frac23E_{\mathrm F}$ at $q=0$, while the Mott–Jones cut stays several times deeper (readout).
- *Marks.* Series styles; the horizontal line at $-\frac23E_{\mathrm F}$, the vertical line at $2k_{\mathrm F}$; Bloch's constant drawn negative; shaded region beyond $q_{\mathrm D}$.
- *Try.* Switch to K: Bloch's 1.70 eV is close to $\frac23E_{\mathrm F}=1.41$ eV; look for the small kink of the Lindhard curve at $2k_{\mathrm F}$. Raise $v_{\mathrm s}$: all couplings drop as $v_{\mathrm s}^{-1/2}$ (zero-point amplitude). Read $\mathrm{Im}/\mathrm{Re}$: dynamic screening is a $10^{-3}$ effect for acoustic phonons.
- *Model.* Jellium with one point ion per cell, no umklapp, static screening without exchange and correlation; equations (A2.23)–(A2.28) and the Mott–Jones form of the text; presets as above; no $N_\ell$ in $g$.

**What in v1 must change.**
- Legend label "Bardeen: Lindhard (G.3–4)" uses source equation numbers: write "Bardeen: Lindhard".
- Readout keys "ε_F" → $E_{\mathrm F}$; "Ω" → $\Omega_{\mathrm V}$; "−⅔ε_F" → $-\frac23E_{\mathrm F}$ (also the plot label of the horizontal line); "V_{MJ}(0)" → "Mott–Jones, $q\to0$"; "Im ε/Re ε at k_F (G.5)" → $\mathrm{Im}\,\varepsilon^{\mathrm e}/\mathrm{Re}\,\varepsilon^{\mathrm e}$ at $q=k_{\mathrm F}$ (Bardeen–Pines). Use MathJax in HTML readouts.
- "Mott–Jones: WS cut" → "Mott–Jones: Wigner–Seitz cut".
- Caption: v1's caption is one block starting with *Model* and uses $\epsilon$, $\varepsilon_F$, $N_p$, $\Omega$ and source equation numbers "(G.5)"; rewrite in the four-part structure with v2 symbols and eqrefs to the section's equations.
- Rename the id and all classes (`viz-a2-bardeen`, `a2-bardeen-…`); keep the physics, the presets and the interaction logic.
