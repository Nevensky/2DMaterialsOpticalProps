# B.1 figure specifications

Section: `v2/sections/B01-foundations.html` (id `b-foundations`, prefix `b1`). Four figures, in reading order:

| # | File | Id | Kind | Tag |
|---|---|---|---|---|
| B1.1 | `B01/fig-umklapp.html` | `viz-b1-umklapp` | new | interactive |
| B1.2 | `B01/fig-vertices.html` | `viz-b1-vertices` | new | diagram |
| B1.3 | `B01/fig-kernel.html` | `viz-b1-kernel` | port of `viz-b-foundations-kernel` | interactive |
| B1.4 | `B01/fig-cell.html` | `viz-b1-cell` | port of `viz-b-foundations-cell` | interactive |

Rules that apply to all four (from `v2/build/BUILDER_GUIDE.md` §5 and `v2/NOTATION.md`):

- CSS classes prefixed with the figure: `b1-umklapp-…`, `b1-vertices-…`, `b1-kernel-…`, `b1-cell-…`.
- Caption: `p.cap-key` with **What to notice**, then `div.cap-more` with **Marks**, **Try**, **Model**; 200–260 words in total. Equation references in captions use `<a class="eqref" href="#b1-eq-…"></a>`; the ids used below exist in the section.
- Notation: v2 only. Gaussian units. $Z^*$ is dimensionless (the dipole is $eZ^*$). Write $\OV$ (not $\Omega_{\rm uc}$), $M_\kappa$, $M_{\mathrm r}$ (not $\mu$, which is the chemical potential), $\boldsymbol\varepsilon_\infty$ with components $\varepsilon_{\infty,xx}$, $\varepsilon_{\infty,zz}$ (never $\epsilon$, never $\perp$/$\parallel$ for the dielectric components: in v2, $\parallel$ means "in the plane of a sheet"), $g^{\mathrm{lr}}$ (not $g^{\mathcal L}$), $g^0$ for the bare vertex, $K_{\mathrm{nn}}$, $K_{\mathrm{nnn}}$ for springs (never $\alpha$, $\beta$, which are Cartesian indices), $\btau_\kappa$ for basis positions. No source equation numbers (no "E.131", "G.180") anywhere in a figure: the section's notes carry them. Upright labels: LO, LA, TO, $\Gamma$, $\mathrm H$, $\mathrm N$, "lr".
- Colours from `EPP.color(...)` at draw time; physics computed in closed form or tiny eigenproblems (`EPP.eigHerm2`), under 15 ms per redraw; fixed axis limits; 390 px phone layout with controls below plots; keyboard access and `aria-label`s; no console errors, no NaN.
- Test each with `bash v2/build/figtest.sh v2/sections/B01/fig-<name>.html` and inspect both screenshots.

Useful constants: $4\pi e^2=180.951$ eV Å (Gaussian; equals $e^2/\varepsilon_0$ in SI); $\hbar/(2M\omega)=2.0901$ Å² $/(M[\text{amu}]\cdot\hbar\omega[\text{meV}])$; $\hbar\omega=64.654$ meV $\times\sqrt{\omega^2[\text{eV Å}^{-2}\text{amu}^{-1}]}$.

---

## B1.1 · `fig-umklapp.html` · `viz-b1-umklapp` · new

**Title.** "Which Fourier components of the ions a lattice wave uses" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. A lattice wave $\bq$ couples to the electron only through the Fourier components $\bq+\bG$ of the ionic potentials: the marks sit on a copy of the reciprocal lattice shifted by $\bq$ (step 4 of the derivation of (B1.3)).
2. Each component has a weight, the projection $(\bq+\bG)\cdot\be$ times the form factor $V_\kappa(\bq+\bG)$, divided by $\OV$. A transverse wave has zero weight at $\bG=0$ and couples only through umklapp components.
3. For the bare potential the $\bG=0$ weight grows as $1/q$ (the Coulomb tail, (B1.5)); Thomas–Fermi screening makes it vanish linearly in $q$. The zero of the empty-core form factor suppresses components near $|\bq+\bG|=\pi/2r_{\mathrm c}$.

**Include position.** In `b1-bare`, after its second paragraph and before `<h4 id="b1-derive">`. Sentences just before the include line:
> "Fourier analysis gives the picture behind the result. The potential of one ion is a superposition of plane waves $\ee^{\ii\bq'\cdot\br}$ of every wavevector $\bq'$. A lattice wave of wavevector $\bq$ shakes the cells with the phases $\ee^{\ii\bq\cdot\bR_\ell}$, and only the plane waves that share this pattern add up over the whole crystal: those with $\bq'=\bq+\bG$, one for each reciprocal lattice vector $\bG$. Each of these components couples to the electron with a **weight**: the amplitude of the plane wave in the ion's potential times the component of the atomic motion along its wavevector, divided by the cell volume. Figure B1.1 draws the components and their weights for sodium."

The text also quotes this figure in `b1-reading` (a transverse wave keeps only umklapp weights) and in `b1-tail` (15 eV/Å bare, 0.6 eV/Å screened, at $q=0.3$ Å⁻¹ along [100]).

**Physics.**
- Weight of the component $\bq+\bG$: the bracket of (B1.4) (`b1-eq-bare`) without the zero-point amplitude, the basis phase and the overlap:
  $w(\bG)=\dfrac{1}{\OV}\,(\bq+\bG)\cdot\be\;V_\kappa(|\bq+\bG|)$ (eV/Å). One atom per cell at the origin, so the basis phase is 1.
- Bare empty-core (Ashcroft) pseudopotential: $V_\kappa(k)=-\dfrac{4\pi Z_\kappa e^2}{k^2}\cos(kr_{\mathrm c})$ (eV Å³).
- Screened option: divide by the Thomas–Fermi dielectric function $1+k_{\mathrm{TF}}^2/k^2$, i.e. $V_\kappa^{\rm scr}(k)=-\dfrac{4\pi Z_\kappa e^2\cos(kr_{\mathrm c})}{k^2+k_{\mathrm{TF}}^2}$ (Bardeen's vertex in its Thomas–Fermi limit, `a2-eq-bardeen`).
- Material: bcc sodium. $a=4.23$ Å, $\OV=a^3/2=37.84$ Å³, $Z_\kappa=1$, $r_{\mathrm c}=0.88$ Å, $k_{\mathrm{TF}}=1.489$ Å⁻¹ (from $k_{\mathrm F}=0.921$ Å⁻¹). Textbook values, illustrative; neither source gives them.
- Reciprocal-lattice plane $G_z=0$: $\bG=\frac{2\pi}{a}(h,k,0)$ with $h+k$ even; $2\pi/a=1.4854$ Å⁻¹. Cut of the first Brillouin zone at $q_z=0$: the square with corners $\mathrm H=\frac{2\pi}{a}(\pm1,0,0)$, $\frac{2\pi}{a}(0,\pm1,0)$, whose edges pass through $\mathrm N=\frac{\pi}{a}(\pm1,\pm1,0)$.
- $\bq=q(\cos\theta,\sin\theta,0)$. Polarization imposed: L, $\be=\hat\bq$; T, $\be=(-\sin\theta,\cos\theta,0)$. (Modes are exactly L or T only along $\Gamma$–H, $\theta=0$, and $\Gamma$–N, $\theta=45°$; elsewhere the figure imposes the polarization to show the geometry. Say so in **Model**.)
- Reference values (checked): $q=0.3$ Å⁻¹, $\theta=0$, L. Bare: $w(0)=-15.39$ eV/Å; largest $|w|$ over $\bG\neq0$ with $|h|,|k|\le4$: 1.41 eV/Å (at $\bG=\frac{2\pi}{a}(\pm2,0,0)$ and equivalents). Screened: $w(0)=-0.600$ eV/Å, largest umklapp 1.17 eV/Å. At $q=0.1$: bare $w(0)=-47.6$, screened $-0.214$ eV/Å. Small-$q$ limits: bare $q\,w(0)\to-4\pi Z_\kappa e^2/\OV=-4.78$ eV; screened $w(0)/q\to-4\pi Z_\kappa e^2/(\OV k_{\mathrm{TF}}^2)=-2.16$ eV Å². T mode: $w(0)=0$ exactly; largest umklapp 1.39 eV/Å (bare, $q=0.3$).

**What to draw.**
- Left panel (geometric, square, `viz-grid two` on desktop, stacked on phone): the plane $(q_x,q_y)$ in Å⁻¹, fixed limits $\pm4.6$ Å⁻¹, axes labelled "$q_x$ (Å⁻¹)", "$q_y$ (Å⁻¹)".
  - Light grey dots at every $\bG$ in range; label $\Gamma$ at the origin.
  - The zone cut: square outline (muted), with upright labels H (at a corner) and N (at an edge midpoint).
  - Arrow from $\Gamma$ to $\bq$ (ink), bold label **q**; a short double-headed arrow at its tip for $\be$ (phonon colour), bold label **e**.
  - At every point $\bq+\bG$ in range: a filled disc of area $\propto|w|$ (for example radius $r=7\,\text{px}\times\sqrt{|w|/(1\ \text{eV/Å})}$, capped at 26 px; a capped disc gets a thin outer ring and its value printed beside it), colour by the sign of $w$ (`c4` negative, `c5` positive, with a legend). Outline the $\bG=0$ disc in ink so it is always identifiable, also when its weight is zero (draw a small hollow ring then).
  - Dashed circle of radius $\pi/2r_{\mathrm c}=1.785$ Å⁻¹ about $\Gamma$, labelled "zero of $V_\kappa$" (the form factor depends on $|\bq+\bG|$ only, so the circle is centred on $\Gamma$).
- Right panel (plot, log–log, fixed limits): $x$ = $q$ from 0.02 to 1.0 Å⁻¹; $y$ = $|w(\bG=0)|$ from 0.01 to 300 eV/Å, for a longitudinal wave along the current direction. Curves: bare (solid, phonon colour) and screened (dashed, phonon colour). Thin muted line: the largest umklapp weight (bare) against $q$, to show that it hardly changes. Vertical cursor at the current $q$ with a point on each curve.

**Controls (4).**
1. $|\bq|$: range 0.02–1.00 Å⁻¹, step 0.01, initial 0.30.
2. Direction $\theta$ of $\bq$ from [100]: 0–90°, step 1°, initial 0°; two small preset buttons "$\Gamma$–H" ($\theta=0$) and "$\Gamma$–N" ($\theta=45°$) may sit in the same row.
3. Polarization: segmented L / T, initial L.
4. Ionic potential: segmented "bare" / "Thomas–Fermi screened", initial bare.

**Readouts.** $q$; $w(\bG=0)$ (eV/Å); largest $|w|$ over $\bG\neq0$ and its $\bG$ written as $(hk0)$; their ratio; for the bare potential $q\,|w(\bG=0)|$ (eV), which stays near 4.8 eV at small $q$ (the $1/q$ law).

**Initial state.** $q=0.30$ Å⁻¹, $\theta=0°$, L, bare: the disc at $\bq$ dominates, ten times the largest umklapp disc.

**Caption, key points.**
- *What to notice.* The discs sit on the reciprocal lattice shifted by $\bq$: these are the only Fourier components of the ionic potentials that a lattice wave $\bq$ can use. For the bare potential the disc at $\bq$ (the $\bG=0$ term) dominates and grows as $1/q$; a transverse wave has no weight there and couples only through umklapp terms.
- *Marks.* Grey dots: $\bG$. Square: cut of the first Brillouin zone at $q_z=0$, with $\Gamma$, H, N. Arrow **q**, double arrow **e**. Disc area: $|w|=|(\bq+\bG)\cdot\be\,V_\kappa(\bq+\bG)|/\OV$; colour: sign. Dashed circle: zero of the empty-core form factor. Right: $|w(\bG=0)|$ against $q$, bare and screened, log axes.
- *Try.* (1) Switch to T: the disc at $\bq$ disappears, the umklapp discs stay. (2) Shrink $q$ with the bare potential: $q\,|w(0)|$ stays near 4.8 eV while the umklapp discs barely change. (3) Switch on screening: at 0.3 Å⁻¹ the $\bG=0$ weight drops from 15 to 0.6 eV/Å and now vanishes linearly; umklapp weights shrink much less. (4) Turn $\bq$ to $\Gamma$–N.
- *Model.* Weight from (B1.4) (`b1-eq-bare`); Coulomb tail (B1.5) (`b1-eq-tail`); empty-core pseudopotential and Thomas–Fermi screening; Na parameters above (textbook, illustrative). The zero-point amplitude is common to all components and is left out; the overlap in (B1.4) decides which component a transition uses (for free electrons, exactly one). Only the plane $G_z=0$ is drawn. Polarizations are exactly L or T only along $\Gamma$–H and $\Gamma$–N.

---

## B1.2 · `fig-vertices.html` · `viz-b1-vertices` · new

**Title.** "Bare, screened and polar vertices" `<span class="tag">diagram</span>`.

**What the reader must learn.**
1. The bare vertex $g^0$ is the phonon shaking the potential of the unscreened ion.
2. The screened vertex $g$ of density-functional perturbation theory is the bare vertex plus every chain of electron–hole bubbles joined by Coulomb lines: $\Delta_{\bq\nu}V^{\mathrm{KS}}=(\varepsilon^{\mathrm e})^{-1}\Delta_{\bq\nu}V^{\mathrm{en}}$. A calculation that starts from $g$ must not add those bubbles again (the double-counting rule of `b1-dft`).
3. The polar vertex $g^{\mathrm{lr}}$ is the phonon's Born dipole $e\mathbf Z^*_\kappa$ talking to the electron through a Coulomb line already screened by the electrons, $4\pi e^2/(\bq\cdot\boldsymbol\varepsilon_\infty\cdot\bq)$. That screened line is what becomes the longitudinal photon in Section B.4.

**Include position.** In `b1-counting`, after its first paragraph and before `<h4 id="b1-dft">`. Sentences just before:
> "A perturbation expansion may then let every particle interact with every other, in any order, with no risk of including a process twice. The price is convergence. Neither interaction is small, the bare phonons are far from the real ones, and the bare vertex is too large by the very factor that screening removes. The two vertices of the two possible starting points are drawn in Figure B1.2."

The text later cites "(Figure B1.2, middle)" for the screened vertex and, in the closing `connect` callout, "the screened Coulomb line of the polar vertex in Figure B1.2". Keep the order bare, screened, polar (left to right, top to bottom on phones).

**Physics.**
- (a) Bare vertex, (B1.4) `b1-eq-bare`: one vertex dot where the phonon line meets the electron line; label $g^0$, sub-label $\Delta_{\bq\nu}V^{\mathrm{en}}$.
- (b) Screened vertex: an equation of diagrams, [blob vertex] $=$ [bare dot] $+$ [dot, bubble, Coulomb line $v$, electron] $+$ [dot, bubble, $v$, bubble, $v$, electron] $+\dots$ The bubble is a closed electron–hole loop ($\chi^0$); the Coulomb line is the bare $v$. The blob carries the label $(\varepsilon^{\mathrm e})^{-1}$; under the panel, the formula $g=\langle\psi_{m\bk+\bq}|(\varepsilon^{\mathrm e})^{-1}\Delta_{\bq\nu}V^{\mathrm{en}}|\psi_{n\bk}\rangle$. RPA; with exchange and correlation $v\to v+f^{\mathrm{xc}}$ (say in **Model** only).
- (c) Polar vertex, (B1.13) `b1-eq-polar` via the map (B1.12) `b1-eq-map`: the phonon line ends at a vertex labelled $e\mathbf Z^*_\kappa$ (bold Z); a screened Coulomb line (double or doubled-dashed, library type `screened`) labelled $4\pi e^2/\bq\cdot\boldsymbol\varepsilon_\infty\cdot\bq$ joins it to the electron line.
- External legs in every panel: incoming electron $n\bk$ (left), outgoing $m\,\bk+\bq$ (right), phonon $\bq\nu$ (from below). Same orientation in all panels.

**What to draw.** `EPP.fd.draw` (SVG): item types `fermion` (electron colour, arrows), `phonon` (wavy, phonon colour), `coulomb` (dashed), `screened`, `vertex`, `blob`, `text`. Panel (b) is wider than (a) and (c); on desktop use a two-row layout if needed: (a) and (c) side by side on the first row, (b) across the full width on the second row, or (a) (b) (c) in one row if it fits. Each panel has a small heading: "(a) bare vertex $g^0$", "(b) screened vertex $g$", "(c) polar vertex $g^{\mathrm{lr}}$". Labels in v2 notation (bold **q**, **k**, **Z**; upright "lr").

**Controls (2 toggles; keyboard accessible).**
1. "Show what a density-functional starting point already contains": shades the bubble chain of (b), and the screened line of (c), with a light band labelled "inside DFPT" and "inside $\boldsymbol\varepsilon_\infty$". Initial: off.
2. "Expand the screened line of (c)": replaces the screened line by its own series, $v+v\chi^0v+\dots$ drawn with the same bubbles as (b), showing that (c) contains the same screening as (b). Initial: off.

No readouts.

**Initial state.** The three diagrams, no shading, (c) unexpanded.

**Caption, key points.**
- *What to notice.* The screened vertex (b) is the bare vertex (a) plus every chain of electron–hole bubbles joined by Coulomb lines; a calculation that starts from (b) must not add those bubbles again. The polar vertex (c) is the long-range part of (b) in a polar insulator: the phonon's Born dipole reaches the electron through a Coulomb line already screened by the electrons.
- *Marks.* Arrowed lines: electrons ($n\bk$ in, $m\,\bk+\bq$ out). Wavy line: phonon $\bq\nu$. Dashed line: bare Coulomb interaction $v$. Doubled line: screened Coulomb interaction $4\pi e^2/\bq\cdot\boldsymbol\varepsilon_\infty\cdot\bq$. Loop: electron–hole bubble $\chi^0$. Blob: $(\varepsilon^{\mathrm e})^{-1}$. Dot in (c): Born dipole $e\mathbf Z^*_\kappa$.
- *Try.* (1) Shade the density-functional content: in a calculation with DFPT vertices everything shaded is already included. (2) Expand the line of (c): its bubbles are those of (b), which is why (c) may not be screened again either.
- *Model.* (a) is (B1.4), `b1-eq-bare`; (b) is the geometric series of the background box "Screening as a sum of bubbles", `a3-eq-scf`, in the random-phase approximation (exchange and correlation replace $v$ by $v+f^{\mathrm{xc}}$); (c) is (B1.13), `b1-eq-polar`, obtained from (a) by the replacement (B1.12), `b1-eq-map`. In Part B the doubled line of (c) becomes the longitudinal photon (Section B.4) and $\boldsymbol\varepsilon_\infty$ a photon self-energy (Section B.5).

---

## B1.3 · `fig-kernel.html` · `viz-b1-kernel` · port

**Port of** `viz-b-foundations-kernel` in `doc/epp-primer/sections/B01-foundations.html` (figure element at the paragraph "Two properties of (E.120) …", script block "Viz 2: anisotropic long-range kernel"). Keep its physics, the two-panel layout and the pointer interaction.

**Title.** "The long-range vertex depends on the direction of $\bq$" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. In the long-wave limit of an optical mode, (B1.16), $q\,|g^{\mathrm{lr}}|$ depends only on the direction $\hat\bq$: the vertex is nonanalytic at $\Gamma$.
2. It vanishes where $\hat\bq$ is perpendicular to the dipole of the mode: an in-plane mode with $\hat\bq$ along the $c$ axis is purely transverse and invisible to electrostatics.
3. An anisotropic $\boldsymbol\varepsilon_\infty$ reshapes the lobes through the denominator $\hat\bq\cdot\boldsymbol\varepsilon_\infty\cdot\hat\bq$; when $\varepsilon_{\infty,xx}>2\varepsilon_{\infty,zz}$ the in-plane lobes tilt away from the plane.

**Include position.** In `b1-dipole`, after the derivation box "Derivation: the potential of a dipole in an anisotropic medium, in real space", before `<h4 id="b1-wave">`. Sentences just before the derivation box (the box is closed by default):
> "It vanishes when $\bp\perp\bq$. Only the longitudinal component of a wave of dipoles produces a scalar potential. A transverse wave of dipoles produces a transverse electric field, which is a photon field and invisible to electrostatics. … When $\boldsymbol\varepsilon_\infty$ is anisotropic its denominator adds a second dependence on the direction. Figure B1.3 shows both for a layered crystal like hexagonal boron nitride."

Later text: "Figure B1.3 plots $q\,|g^{\mathrm{lr}}|$ against the direction" (`b1-optical`), and `b1-where` cites it as unaffected by the basis phase.

**Physics.**
- (B1.16) `b1-eq-opt`: $q\,|g^{\mathrm{lr}}|=\dfrac{4\pi e}{\OV}\,\dfrac{|\hat\bq\cdot\bp_\nu|}{\hat\bq\cdot\boldsymbol\varepsilon_\infty\cdot\hat\bq}$, with the cell dipole of one B–N pair moving against itself along $\hat\be$: $\bp_\nu=e\,\mathbf Z^*_{\mathrm B}\cdot\hat\be\,(\hbar/2M_{\mathrm r}\omega_\nu)^{1/2}$ (two atoms, $\mathbf Z^*_{\mathrm N}=-\mathbf Z^*_{\mathrm B}$, eigenvector $\be_{\mathrm B}=(M_{\mathrm N}/M_{\mathrm c})^{1/2}\hat\be$, $\be_{\mathrm N}=-(M_{\mathrm B}/M_{\mathrm c})^{1/2}\hat\be$). Hence
  $q\,|g^{\mathrm{lr}}|=\dfrac{4\pi e^2}{\OV}\,F(\hat\bq)\,\Big(\dfrac{\hbar}{2M_{\mathrm r}\omega_\nu}\Big)^{1/2}$, with the dimensionless angular factor $F(\hat\bq)=\dfrac{|\hat\bq\cdot\mathbf Z^*_{\mathrm B}\cdot\hat\be|}{\hat\bq\cdot\boldsymbol\varepsilon_\infty\cdot\hat\bq}$.
  (v1 code: `qg = E2EPS / OM * K * Math.sqrt(L2 / (MU * hw))` with `E2EPS = 180.95`, `L2 = 2.0901`: correct, keep, rename.)
- Uniaxial crystal, $c$ axis $=z$: $\boldsymbol\varepsilon_\infty=\mathrm{diag}(\varepsilon_{\infty,xx},\varepsilon_{\infty,xx},\varepsilon_{\infty,zz})$, $\mathbf Z^*_{\mathrm B}=\mathrm{diag}(Z^*_{xx},Z^*_{xx},Z^*_{zz})$; $\hat\bq=(\sin\theta,0,\cos\theta)$, $\theta$ from the $c$ axis. In-plane mode, $\hat\be=\hat x$: $F=Z^*_{xx}|\sin\theta|/(\varepsilon_{\infty,xx}\sin^2\theta+\varepsilon_{\infty,zz}\cos^2\theta)$. Out-of-plane mode, $\hat\be=\hat z$: $F=Z^*_{zz}|\cos\theta|/(\dots)$.
- Parameters (hBN-bulk-like, literature values, illustrative; neither source gives them): $\varepsilon_{\infty,xx}=4.87$, $\varepsilon_{\infty,zz}=2.95$, $Z^*_{xx}=2.71$, $Z^*_{zz}=0.81$; $\OV=18.08$ Å³ (one B–N pair); $M_{\mathrm B}=10.81$, $M_{\mathrm N}=14.007$ amu, $M_{\mathrm r}=6.101$ amu; $\hbar\omega_\nu$ held at the transverse values 169.9 meV (in-plane) and 96.7 meV (out-of-plane), as for the short-range phonons of the EPP manuscript.
- Reference values (checked): in-plane, $\theta=90°$: $F=0.5565$, $q|g^{\mathrm{lr}}|=0.2501$ eV Å⁻¹; in-plane, $\theta=50°$: $\hat\bq\cdot\boldsymbol\varepsilon_\infty\cdot\hat\bq=4.077$, $F=0.5092$, $0.2289$ eV Å⁻¹; out-of-plane, $\theta=0$: $F=0.2746$, $0.1636$ eV Å⁻¹. Tilt test $\varepsilon_{\infty,xx}=12$, $\varepsilon_{\infty,zz}=4$: in-plane maximum at $\theta=45°$ ($\sin^2\theta=\varepsilon_{\infty,zz}/(\varepsilon_{\infty,xx}-\varepsilon_{\infty,zz})$).
- Optional readout, the frequency of the mode with its own dipole field included (the single-mode form of `a4-eq-lot`): $(\hbar\omega)^2_{\rm field}=(\hbar\omega_\nu)^2+(64.654\ \text{meV})^2\times\dfrac{4\pi e^2}{\OV M_{\mathrm r}}\dfrac{(\hat\bq\cdot\mathbf Z^*_{\mathrm B}\cdot\hat\be)^2}{\hat\bq\cdot\boldsymbol\varepsilon_\infty\cdot\hat\bq}$ (with $4\pi e^2$ in eV Å, $\OV$ in Å³, $M_{\mathrm r}$ in amu). In-plane at $\theta=90°$: 198.0 meV (measured LO 199.6 meV); at $\theta=50°$: 190.0 meV. With this frequency $q|g^{\mathrm{lr}}|$ is lower by $(\omega_\nu/\omega_{\rm field})^{1/2}$: 7% at $\theta=90°$, the remark at the end of `b1-wave`.

**What to draw.**
- Left: polar plot in the $(x,z)$ plane, $c$ axis vertical: radius $=q|g^{\mathrm{lr}}|$ (eV Å⁻¹) along $\hat\bq$, for both modes (in-plane: solid, phonon colour; out-of-plane: dashed, `c5`); the selected mode bold, the other thin. Dotted rings with their values (fixed set of ring values; the plot range picks the smallest ring above the maximum, as in v1). Arrow from the origin along $\hat\bq$ (ink, bold label **q̂**) and a point where it meets the selected curve. Axis labels "$x$, in plane (eV Å⁻¹)" and "$z$, $c$ axis (eV Å⁻¹)".
- Small inset (new, geometric; top-left corner of the polar plot or beside it): two hexagonal sheets seen edge-on with the $c$ axis, a B and an N atom with arrows for the two polar modes (in-plane: antiparallel arrows along $x$; out-of-plane: antiparallel arrows along $z$). It tells the reader what "in-plane" and "out-of-plane" mean.
- Right: $q|g^{\mathrm{lr}}|$ against $\theta$ from 0 to 90° for both modes, same styles, vertical line at the current angle folded into 0–90°. Axis labels "$\theta$, angle between $\hat\bq$ and the $c$ axis (deg)" and "$q\,|g^{\mathrm{lr}}|$ (eV Å⁻¹)".

**Controls (5).**
1. $\varepsilon_{\infty,xx}$ (in plane): 1–16, step 0.01, initial 4.87.
2. $\varepsilon_{\infty,zz}$ (along $c$): 1–16, step 0.01, initial 2.95.
3. Direction $\theta$ of $\hat\bq$ from the $c$ axis: 0–360°, step 1°, initial 50° (keyboard access; clicking or dragging on either plot also sets it, as in v1).
4. Mode: segmented "in-plane" / "out-of-plane", initial in-plane.
5. Presets (one row of buttons): "hBN-like" (4.87, 2.95) and "isotropic" ($\varepsilon_{\infty,xx}=\varepsilon_{\infty,zz}=4.87$). The Born charges stay at the hBN values ($Z^*_{xx}=2.71$, $Z^*_{zz}=0.81$); drop v1's two $Z^*$ sliders.

**Readouts.** $\theta$; $\hat\bq\cdot\boldsymbol\varepsilon_\infty\cdot\hat\bq$; $F(\hat\bq)$; $q|g^{\mathrm{lr}}|$ (eV Å⁻¹) with $\hbar\omega_\nu$; optionally $\hbar\omega_{\rm field}$ (meV) and $q|g^{\mathrm{lr}}|$ evaluated with it.

**Initial state.** hBN-like preset, in-plane mode, $\theta=50°$.

**Caption, key points.**
- *What to notice.* At long wavelength the polar vertex of an optical mode is $1/q$ times a function of the direction of $\bq$ alone. It vanishes where $\hat\bq$ is perpendicular to the dipole of the mode: an in-plane mode with $\hat\bq$ along the $c$ axis is purely transverse and couples to no electrostatic potential.
- *Marks.* Polar plot: distance from the origin $=q|g^{\mathrm{lr}}|$ along $\hat\bq$; solid in-plane mode, dashed out-of-plane mode, bold the selected one; arrow **q̂**. Right: the same against $\theta$. Inset: the two polar modes of a B–N pair.
- *Try.* (1) Choose "isotropic": the lobes become the circles of a dipole pattern, $\propto|\cos|$ of the angle to the dipole. (2) Set $\varepsilon_{\infty,xx}=12$, $\varepsilon_{\infty,zz}=4$: the in-plane lobes tilt to 45°, because the denominator grows faster than the numerator near the plane. (3) Drag $\hat\bq$ around and watch $F$ and $q|g^{\mathrm{lr}}|$ change at fixed $q$.
- *Model.* (B1.16), `b1-eq-opt`, for one B–N pair; $\hbar\omega_\nu$ held at the transverse values, as for the short-range starting phonons of the EPP manuscript; Section A.4's vertex uses the full frequency, 198 meV for the in-plane mode with $\hat\bq$ along its dipole, which lowers $q|g^{\mathrm{lr}}|$ by 7%. At $q\to0$ the basis phases drop out, so the vertex with and without them agree here. Parameters as above (literature, illustrative).

**What in v1 is wrong or uses source notation, and must change.**
- Title and caption cite "(E.131′)", "(E.120)", "(E.131)", "(G.180)", "Verdi and Giustino's (G.180)", "Giustino's convention (no N)": remove all source equation numbers; cite (B1.16) and (B1.13) with `a.eqref`; say "per unit phonon field, with $N_\ell^{-1/2}$ outside" if normalization needs a word.
- $\hat\epsilon^\infty$, $\epsilon^\infty_\perp$ (in plane), $\epsilon^\infty_\parallel$ (along $c$): use $\boldsymbol\varepsilon_\infty$, $\varepsilon_{\infty,xx}$, $\varepsilon_{\infty,zz}$. The v1 meaning of $\parallel$ (along $c$) is the opposite of v2's 2D convention.
- $\hat Z^*$, $Z^*_\perp$, $Z^*_\parallel$ "in units of $|e|$": $\mathbf Z^*_{\mathrm B}$ with $Z^*_{xx}$, $Z^*_{zz}$, dimensionless.
- The symbol $K$ for the kernel clashes with the spring stiffness $K$ of NOTATION: call it the angular factor $F$ (or show it only as a readout labelled "angular factor").
- $\mu$ for the reduced mass: $M_{\mathrm r}$. $\Omega_{\rm uc}$: $\OV$. $g^{\mathcal L}$: $g^{\mathrm{lr}}$.
- v1 plots the dimensionless kernel; v2 plots $q|g^{\mathrm{lr}}|$ in eV Å⁻¹ on both panels (physical, and comparable between the two modes); keep $F$ as a readout.
- v1 has seven controls; reduce to five as listed, and add the $\theta$ slider for keyboard access.
- Caption restructured into What to notice / Marks / Try / Model.

---

## B1.4 · `fig-cell.html` · `viz-b1-cell` · port

**Port of** `viz-b-foundations-cell` in `doc/epp-primer/sections/B01-foundations.html` (figure element after the derivation "acoustic modes and the spurious coupling", script block "Viz 3: cell-choice test"). Keep its physics (checked, numbers below) and its plot; add the geometric panel.

**Title.** "Same crystal, two choices of cell" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. The same crystal can be grouped into cells in different ways; here the Ga plane at the origin is grouped either with the As plane above it or with the one below.
2. Switching cells turns the As eigenvector component by $\ee^{\ii\bq\cdot\bR_{\ell_0}}$ (step 2 before (B1.15)), but leaves the product $e_{\kappa\alpha,\nu}(\bq)\,\ee^{-\ii\bq\cdot\btau_\kappa}$ unchanged, (B1.15) (`b1-eq-cell`).
3. Hence the polar vertex with the basis phase, (B1.13), does not depend on the bookkeeping, while the phase-less form does: about 20% in $|g^{\mathrm{lr}}|^2$ of the LO mode at 0.5 Å⁻¹, and a spurious LA coupling growing as $\sqrt q$ ($|g|^2\propto1/q$) that survives even in a centrosymmetric toy.

**Include position.** In `b1-phase`, after its first paragraph and before the `ol.steps`. Sentences just before:
> "Nothing physical may depend on that choice. The phase convention of Section A.1, however, attaches the phase of a wave to the cell, so assigning an atom to a different cell changes the phase of its eigenvector component. Figure B1.4 draws a chain of GaAs planes along [111] booked in two ways and computes the polar vertex with and without the basis phase."

Later text quotes it: in `b1-where` (0.5–0.8% at 0.1 Å⁻¹, about 20% at 0.5 Å⁻¹, up for one cell and down for the other) and in the notes (GaAs piezoelectric numbers use this chain).

**Physics (v1 model, verified).**
- Longitudinal waves along [111] of zincblende GaAs, written as a chain of (111) planes with period $a/\sqrt3=3.2638$ Å ($a=5.653$ Å). Each As plane lies 2.4478 Å ($\sqrt3a/4$) above one Ga plane, joined to it by the single bond along [111] (plane force constant $K_{\mathrm{nn}}$), and 0.8159 Å ($a/4\sqrt3$) below the next Ga plane, joined to it by three inclined bonds (plane force constant $K_{\mathrm{nn}}/3$). Next-nearest neighbours give $2K_{\mathrm{nnn}}$ between adjacent like planes (Ga–Ga and As–As). $K_{\mathrm{nn}}=7.185$ eV/Å² (fitted to $\hbar\omega_{\mathrm{TO}}=33.3$ meV), $K_{\mathrm{nnn}}=0.451$ eV/Å² (fitted to the [100] sound velocity). $M_{\mathrm{Ga}}=69.72$, $M_{\mathrm{As}}=74.92$ amu; $Z^*_{\mathrm{Ga}}=+2.07$, $Z^*_{\mathrm{As}}=-2.07$; $\varepsilon_\infty=10.9$; $\OV=45.2$ Å³ (literature values, illustrative).
- Cells: Ga at the cell corner, $\tau_{\mathrm{Ga}}=0$. Cell A ("As above Ga"): $\tau_{\mathrm{As}}=+2.448$ Å. Cell B ("As below Ga"): $\tau_{\mathrm{As}}=-0.816$ Å. From A to B the As atom is reassigned by $\bR_{\ell_0}=-a/\sqrt3$ along [111], (B1.14).
- Dynamical matrix as in (A1.18) (`a1-eq-dyn`), phase on the partner's cell: on-site $\Phi_{\kappa\kappa}=K_{\rm in}+K_{\rm out}+2\,(2K_{\mathrm{nnn}})(1-\cos(qa/\sqrt3))$ for both atoms; $\Phi_{\mathrm{Ga,As}}(q)=-\big(K_{\rm in}+K_{\rm out}\,\ee^{\ii s\,qa/\sqrt3}\big)$ with, for cell A, $K_{\rm in}=K_{\mathrm{nn}}$, $K_{\rm out}=K_{\mathrm{nn}}/3$, $s=-1$, and for cell B, $K_{\rm in}=K_{\mathrm{nn}}/3$, $K_{\rm out}=K_{\mathrm{nn}}$, $s=+1$. Mass-scale and diagonalize with `EPP.eigHerm2` (v1 code does exactly this; keep).
- Centrosymmetric toy: equal springs $2K_{\mathrm{nn}}/3$ on both sides, As plane midway ($\tau_{\mathrm{As}}=\pm1.632$ Å), same masses and charges.
- Vertex: (B1.13) at $\bG=0$ for plane-wave electrons (overlap 1; only $\bG=0$ contributes): with the basis phase (solid) $q|g^{\mathrm{lr}}|=\dfrac{4\pi e^2}{\OV\varepsilon_\infty}\Big(\dfrac{\hbar}{2\omega}\Big)^{1/2}\Big|\sum_\kappa\dfrac{Z^*_\kappa e_\kappa}{\sqrt{M_\kappa}}\ee^{-\ii q\tau_\kappa}\Big|$; without it (dashed) the same with $\ee^{-\ii q\tau_\kappa}\to1$. Numerically `PRE = 180.95 / 45.2 / 10.9` and $(\hbar/2\omega)^{1/2}=\sqrt{2.0901/\hbar\omega[\text{meV}]}$ Å amu$^{1/2}$.
- Reference values (checked; eV Å⁻¹). GaAs-like, $q=0.1$ Å⁻¹: LO with phase 0.0316; without: cell A 0.0317 ($|g|^2$ ratio 1.005), cell B 0.0315 (0.992). $q=0.5$ Å⁻¹: LO with phase 0.0287; without: A 0.0314 (1.198), B 0.0253 (0.778). LA with phase 0.0082 ($q=0.1$), 0.0202 ($q=0.5$); LA without: A 0.0041 / 0.0076, B 0.0123 / 0.0283. Centrosymmetric toy: LA with phase 0.0000 ($q=0.1$), 0.0006 ($q=0.5$); without: 0.0079, 0.0177 (same for both cells, which are mirror images); LO ratio 0.993, 0.843. $\hbar\omega_{\mathrm{TO}}\to33.3$ meV as $q\to0$; at $q=0.5$: $\hbar\omega_{\mathrm{LO}}=33.10$, $\hbar\omega_{\mathrm{LA}}=15.32$ meV; LA sound velocity 5.06 km/s.
- Phasors at $q=0.5$ Å⁻¹, with the overall phase fixed so that the Ga component is real and positive (`EPP.fixPhase` or equivalent): LO, cell A: $e_{\mathrm{As}}=0.687\angle{-161.2°}$, product $e_{\mathrm{As}}\ee^{-\ii q\tau_{\mathrm{As}}}=0.687\angle128.6°$; cell B: $e_{\mathrm{As}}=0.687\angle105.3°$, product $0.687\angle128.6°$ (the same). LA, cell A: $0.726\angle18.8°$, product $0.726\angle{-51.4°}$; cell B: $0.726\angle{-74.7°}$, product $0.726\angle{-51.4°}$. The eigenvector phasor turns by $-qa/\sqrt3=-93.5°$ from A to B; the product does not move.

**What to draw.**
- Top panel (new, geometric): the chain of planes drawn along a horizontal [111] axis, about four periods, to scale: Ga planes (`c4`) and As planes (`c5`) as discs or short vertical bars, the short (0.816 Å) and long (2.448 Å) gaps visibly different; thin lines for the single bond (long gap) and a triple line or "×3" mark for the three inclined bonds (short gap). Shaded boxes (or brackets) for the cells of the current choice, one box per period, each containing one Ga and one As plane; tick marks at the cell corners labelled $\bR_\ell$ for a few cells. In cell 0, an arrow from the corner to its As plane labelled $\tau_{\mathrm{As}}$; when the choice is B, a dashed ghost of the A choice and an arrow between the two As planes labelled $\bR_{\ell_0}$ (bold R; length $a/\sqrt3$).
- Under the chain, two pairs of phasor dials for the selected mode (LO by default) at the cursor $q$ (unit circles, as in `A01/fig-honeycomb.html` or `A01/fig-diatomic.html`): pair 1, "eigenvector $e_\kappa$": the Ga and As components; pair 2, "with basis phase $e_\kappa\,\ee^{-\ii q\tau_\kappa}$". When the cell changes, the As hand of pair 1 jumps by $-93.5°$ (at 0.5 Å⁻¹) and the As hand of pair 2 stays: the visual content of (B1.15). Print modulus and angle under each hand.
- Bottom panel (as v1): plot of $q|g^{\mathrm{lr}}|$ (eV Å⁻¹, fixed 0–0.045) against $q$ along [111] (0–0.9 Å⁻¹): LO with phase (solid, phonon colour), LO without (dashed, phonon colour), LA with phase (solid, `c5`), LA without (dashed, `c5`); vertical cursor at $q$ with points on the four curves; legend "with basis phase, (B1.13)" and "without basis phase".

**Controls (4).**
1. Chain: segmented "GaAs [111]-like" / "centrosymmetric toy", initial GaAs-like.
2. Cell: segmented "As above Ga ($\tau_{\mathrm{As}}=+2.448$ Å)" / "As below Ga ($\tau_{\mathrm{As}}=-0.816$ Å)" (toy: $\pm1.632$ Å), initial "above".
3. Cursor $q$: range 0.01–0.90 Å⁻¹, step 0.01, initial 0.50; clicking the plot also sets it (v1 behaviour).
4. Mode shown on the dials: segmented LO / LA, initial LO.

**Readouts.** $\tau_{\mathrm{As}}$ (Å); $q\,\tau_{\mathrm{As}}$ (rad); $\hbar\omega_{\mathrm{LO}}$, $\hbar\omega_{\mathrm{LA}}$ (meV); $q|g_{\mathrm{LO}}|$ with / without phase (eV Å⁻¹); $|g_{\mathrm{LO}}|^2$ without / with; $q|g_{\mathrm{LA}}|$ with / without (eV Å⁻¹).

**Initial state.** GaAs-like, As above Ga, $q=0.50$ Å⁻¹, dials on LO: the dashed LO curve sits 20% (in $|g|^2$) above the solid one; switching the cell moves it 22% below.

**Caption, key points.**
- *What to notice.* The solid curves do not move when the cell changes; the dashed ones do, by about 20% in $|g_{\mathrm{LO}}|^2$ at 0.5 Å⁻¹, up for one cell and down for the other. On the dials, the As component of the eigenvector turns when the cell changes, while its product with the basis phase stays put.
- *Marks.* Chain: Ga and As planes along [111], shaded boxes = the chosen cells, $\tau_{\mathrm{As}}$ and $\bR_{\ell_0}$ arrows. Dials: eigenvector components and their products with $\ee^{-\ii q\tau_\kappa}$, Ga fixed real. Plot: LO (phonon colour) and LA (`c5`); solid with the basis phase, dashed without.
- *Try.* (1) Toggle the cell at 0.5 Å⁻¹: the dashed LO value jumps between 0.0314 and 0.0253 eV Å⁻¹, the solid one stays at 0.0287. (2) Move the cursor towards $q=0$: the dashed LA curve grows as $\sqrt q$, which means $|g|^2\propto1/q$. (3) Choose the centrosymmetric toy: the solid LA curve almost vanishes, the dashed one does not, and the two cells give the same dashed curves (mirror images).
- *Model.* (B1.13) at $\bG=0$ (`b1-eq-polar`) for plane-wave electrons, per unit phonon field with $N_\ell^{-1/2}$ outside; invariance (B1.15) (`b1-eq-cell`); eigenvectors from the dynamical matrix with the phase on the partner's cell, (A1.18) (`a1-eq-dyn`); the chain model and inputs above. In GaAs the solid LA curve is the internal-strain dipole of (B1.18); in the real crystal the response of the electrons at fixed ions, carried by dynamical quadrupoles that this vertex lacks, largely cancels it.

**What in v1 is wrong or uses source notation, and must change.**
- Legends, readouts and caption say "corrected (E.131′)", "printed (E.131)", "Giustino's sign of (E.15)", "§Z.3", "Giustino's normalization (no N)": remove; use "with basis phase, (B1.13)" and "without basis phase"; the sign convention is that of (A1.18).
- Springs called $\alpha$ and $\beta$ (Cartesian index letters): $K_{\mathrm{nn}}$, $K_{\mathrm{nnn}}$.
- Period called $c$ (the speed of light in v2): write $a/\sqrt3$.
- $m_{\rm Ga}$, $m_{\rm As}$: $M_{\mathrm{Ga}}$, $M_{\mathrm{As}}$. $\Omega_{\rm uc}$: $\OV$. $\epsilon^\infty$: $\varepsilon_\infty$. $g^{\mathcal L}$: $g^{\mathrm{lr}}$. Readout "$r_{\rm As}$" and "$q(r_{\rm As}-r_{\rm Ga})$": $\tau_{\mathrm{As}}$ and $q\,\tau_{\mathrm{As}}$. Born charges "$\pm2.07$" without "|e|" units (dimensionless).
- v1 has no drawing of the two cell choices; the brief requires one (top panel and dials above).
- Caption restructured into What to notice / Marks / Try / Model; v1's long Model paragraph (bond geometry and fits) becomes the Model part, shortened.
