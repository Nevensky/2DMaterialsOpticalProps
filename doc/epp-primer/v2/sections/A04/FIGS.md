# A.4 figure specifications (for the figure authors)

Section file: `v2/sections/A04-polar.html` (section id `a-polar`, prefix `a4`). Five figures, in reading order:

| # | file | id | kind |
|---|---|---|---|
| A4.1 | `A04/fig-born.html` | `viz-a4-born` | new |
| A4.2 | `A04/fig-dipoles.html` | `viz-a4-dipoles` | port of v1 `viz-a-polar-dipoles` |
| A4.3 | `A04/fig-lot.html` | `viz-a4-lot` | new |
| A4.4 | `A04/fig-frohlich.html` | `viz-a4-frohlich` | port of v1 `viz-a-polar-frohlich` |
| A4.5 | `A04/fig-wannier.html` | `viz-a4-wannier` | port of v1 `viz-a-polar-wannier` |

General rules for all five (BUILDER_GUIDE §5 and NOTATION):
- Gaussian units in every formula; numbers in eV, meV, Å, amu. $e^2=14.3996$ eV Å, $4\pi e^2=180.95$ eV Å, $\hbar c=1973.27$ eV Å,
  $\hbar^2/2m_{\mathrm e}=3.80998$ eV Å², $\hbar\omega[\text{meV}]=64.654\sqrt{\omega^2[\text{eV Å}^{-2}\,\text{amu}^{-1}]}$,
  $\hbar/(2M\omega)=2.0901$ Å² $/(M[\text{amu}]\,\hbar\omega[\text{meV}])$.
- v2 symbols only (no source symbols anywhere, including strings, aria-labels, legends): Born charge $\mathbf Z^*_\kappa$, $Z^*$ (dimensionless);
  $\varepsilon_\infty$ (canvas `ε_∞`), $\varepsilon_{\mathrm s}$ (`ε_s`), never ε^∞, ε^0, ϵ; cell volume $\OV$ (`Ω_V`); reduced mass $M_{\mathrm r}$ (`M_r`), never μ
  (μ is the chemical potential); cell mass $M_{\mathrm c}$ (`M_c`); $\omega_{\mathrm{LO}}$, $\omega_{\mathrm{TO}}$ (`ω_{LO}`, `ω_{TO}`, upright labels);
  long-range and short-range vertex $g^{\mathrm{lr}}$, $g^{\mathrm{sr}}$ (`g^{lr}`, `g^{sr}`), never $g^{\mathcal L}$, $g^{\mathcal S}$; Fröhlich constant $\alpha_{\mathrm F}$ (`α_F`);
  polaron length $r_{\mathrm p}$, inverse polaron length $1/r_{\mathrm p}$; light-cone wavevector $q_{\mathrm c}$; cutoff length $r_{\mathrm c}$; coarse grid $N_{\mathrm c}$;
  deformation potential $E_{\mathrm{def}}$ (as in Section A.2; not $E_1$); displacements $\bu_\kappa$, cells $\ell$, $\bR_\ell$ (never $p$, $\bR_p$, $\Delta\tau$, $N_p$);
  electron potential energy $V^{\mathrm{lr}}$; electric field $\bE$; polarization $\bP$; dipole $\mathbf p$. **σ is reserved for the conductivity**: never use it.
- Equation references in captions use `<a class="eqref" href="#a4-eq-…"></a>`; available ids: `a4-eq-born` (A4.1), `a4-eq-force` (A4.2),
  `a4-eq-split` (A4.3), `a4-eq-asr` (A4.4), `a4-eq-P` (A4.5), `a4-eq-E` (A4.6), `a4-eq-dd` (A4.7), `a4-eq-lot` (A4.8), `a4-eq-epss` (A4.9),
  `a4-eq-lst` (A4.10), `a4-eq-Vlr` (A4.11), `a4-eq-gmacro` (A4.12), `a4-eq-Vdip` (A4.13), `a4-eq-glr` (A4.14), `a4-eq-frohlich` (A4.15),
  `a4-eq-alphaF` (A4.16), `a4-eq-gR` (A4.17), `a4-eq-interp` (A4.18), `a4-eq-qc` (A4.19). Never write source equation numbers (G.180, E.131, …).
- Literature numbers are "literature values, illustrative" in the caption. Caption 200–260 words: What to notice / Marks / Try / Model.
- Colours from `EPP.color(...)` at draw time; displacements in the phonon colour, electric quantities (dipoles, field, bound charge, potential)
  in the photon colour family, electron quantities in the electron colour, as elsewhere in the primer.

---

## A4.1 `fig-born.html` · `viz-a4-born` · new

**Title.** "A Born charge is a tensor: the dipole created by a displaced atom" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. The dipole created by displacing an atom is linear in the displacement, $\mathbf p=e\,\mathbf Z^*_\kappa\cdot\bu_\kappa$, and is built column by
   column: $\mathbf p/e=u_{\kappa x}\,(\text{column }x)+u_{\kappa y}\,(\text{column }y)$, column $\beta$ being the dipole per unit displacement along $\beta$.
2. At a site of cubic symmetry the tensor is a number (dipole parallel to the displacement); at a site of lower symmetry (O in BaTiO₃) the
   dipole is not parallel to the displacement and is larger along the bond.
3. The Born charge differs from the nominal ionic charge: smaller in GaAs (2.07 against 3), much larger in BaTiO₃ (Ti 7.25 against 4,
   O −5.71 along the bond against −2): electrons follow and flow along bonds.

**Include position.** In subsection `a4-born-tensor` ("A tensor, not a number"), right after its only paragraph. The sentences just before:
"The oxygen of cubic BaTiO₃ sits between two titanium atoms, and only its bond axis is special: its tensor is diagonal, with $-5.71$ along
the bond and $-2.15$ across it. Displaced at an angle to the bond, it creates a dipole turned toward the bond direction (Figure A4.1). At
sites of still lower symmetry the tensor need not even be symmetric."

**Physics.** (A4.1) rearranged: dipole added to the cell $\mathbf p=e\,\mathbf Z^*_\kappa\cdot\bu_\kappa$, $p_\alpha=e\sum_\beta Z^*_{\kappa,\alpha\beta}u_{\kappa\beta}$;
nominal (rigid-ion) comparison $\mathbf p^{\mathrm{nom}}=e\,Z^{\mathrm{nom}}_\kappa\bu_\kappa$ (dashed). In-plane ($x,y$) components only. Values
(literature values, illustrative; cubic BaTiO₃ from Ghosez, Michenaud and Gonze, PRB 58, 6224 (1998)):
- GaAs, atom Ga, zincblende seen along [001]: $\mathbf Z^*=2.07\,\mathbb 1$, nominal +3 (As: −2.07, nominal −3).
- BaTiO₃ (cubic), atom Ti at a corner of the TiO₂ square net: $\mathbf Z^*=7.25\,\mathbb 1$, nominal +4.
- BaTiO₃, atom O on a Ti–O–Ti bond along $x$: $\mathbf Z^*=\mathrm{diag}(-5.71,-2.15)$ (in-plane block), nominal −2.
  (An O on a bond along $y$ has $\mathrm{diag}(-2.15,-5.71)$; Ba, +2.77, sits above the centre of the square.)
- Sum-rule readout for BaTiO₃: $xx$ components $2.77+7.25-5.71-2.15-2.15=0.01\approx0$ ((A4.4)); for GaAs $2.07-2.07=0$.

**What to draw.** One panel (canvas, square-ish, aspect ≈ 0.75 on desktop, full width on phone).
- A 2D net of 3×3 cells: for BaTiO₃, Ti at the square corners (filled circles), O at the edge midpoints (smaller filled circles, other colour),
  Ba as an open circle at each square centre (marked "above the plane" in Marks); for GaAs, a checkerboard of Ga and As columns.
- The selected atom (centre of the drawing) is drawn displaced by $\bu$ (exaggerated, e.g. 1 Å on screen per 0.05 Å physical, state the scale),
  with a ghost at its equilibrium site and the displacement arrow $\bu$ (phonon colour, label `u` bold).
- The dipole arrow $\mathbf p/e$ (photon colour, bold `p/e`) drawn from the atom, at a fixed scale stated in the legend (e.g. 100 px per e Å).
- The two columns as thin arrows: $u_{\kappa x}\,(\mathbf Z^*\hat{\mathbf x})$ and $u_{\kappa y}\,(\mathbf Z^*\hat{\mathbf y})$, drawn tip-to-tail
  (parallelogram construction) so that their sum is visibly $\mathbf p/e$; label them "column x × u_x", "column y × u_y".
- Optional dashed arrow: the nominal dipole $Z^{\mathrm{nom}}\bu$ (toggle).
- A small arc marking the angle between $\bu$ and $\mathbf p$ (or $-\mathbf p$ for an anion).
- Label the Ti–O–Ti bond axis when O is selected.

**Controls (4).**
- seg "atom": `Ga in GaAs` / `Ti in BaTiO₃` / `O in BaTiO₃`.
- range "direction of $\bu$" 0–360° (step 1°, aria-valuetext with degrees); the displacement can also be dragged on the canvas
  (`touch-action: none` only on the canvas; arrow keys rotate by 5°).
- range "$|\bu|$" 0–0.10 Å (step 0.005 Å).
- toggle "show nominal charge".

**Readouts.** $\bu$ (Å, components), $\mathbf p/e$ (e Å, components), $|\mathbf p|/(e|\bu|)$ (dimensionless), angle between $\bu$ and $\mathbf p$ (degrees,
measured from $\bu$ to $\mathbf p$ for a cation and to $-\mathbf p$ for an anion, so that it reads 0° at a cubic site), and the sum $\sum_\kappa Z^*_{\kappa,xx}$
of the selected crystal.

**Initial state.** O in BaTiO₃, $|\bu|=0.05$ Å at 30° from the bond ($x$) axis, nominal toggle off. Then $\mathbf p/e=(-0.247,-0.054)$ e Å, a dipole
pointing nearly opposite to $\bu$ but rotated 18° toward the bond: the key phenomenon is visible at once.

**Caption key points.**
- *What to notice.* The dipole is the sum of the two columns of $\mathbf Z^*$ weighted by the displacement components; at the oxygen of BaTiO₃ it
  is 2.7 times larger along the Ti–O bond than across it, so it is not parallel to the displacement. At the cubic Ga and Ti sites it is.
- *Marks.* Ghost = equilibrium site; phonon-coloured arrow = $\bu$; photon-coloured arrow = $\mathbf p/e$; thin arrows = the columns times
  $u_{\kappa x}$, $u_{\kappa y}$; dashed = nominal charge; open circles = Ba above the plane; scales of both arrows.
- *Try.* Rotate $\bu$ through 360° with O selected and watch the dipole lag behind or run ahead of it, aligning only along and across the bond;
  switch to Ti and see the angle drop to zero and the dipole grow to 7.25 e $|\bu|$ (nominal 4); switch on the nominal charge for Ga to see
  2.07 against 3.
- *Model.* $\mathbf p=e\,\mathbf Z^*_\kappa\cdot\bu_\kappa$ from (A4.1); sum rule (A4.4); values above, literature values, illustrative (cubic BaTiO₃ calculated
  by Ghosez, Michenaud and Gonze 1998); in-plane block only.

---

## A4.2 `fig-dipoles.html` · `viz-a4-dipoles` · port of v1 `viz-a-polar-dipoles` (v1 file `doc/epp-primer/sections/A04-polar.html`, static SVG)

**Title.** "A wave of Born-charge dipoles: LO makes a field, TO does not" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. In an LO wave the cell dipoles point along $\bq$; heads converge on planes half a wavelength apart, leaving sheets of positive and negative
   bound charge, with a longitudinal field between them, opposite to $\bP$. In a TO wave the dipoles lie side by side in each plane: no bound
   charge, no field.
2. The field amplitude does not depend on the wavelength; the potential energy of an electron, field times distance, grows in proportion to the
   wavelength ($\propto1/q$), a quarter wavelength out of phase with $\bP$.
3. The numbers are large: 0.76 MV/cm for a relative displacement of 0.01 Å in GaAs.

**Include position.** In subsection `a4-gauss` ("Gauss's law in a wave"), right after its first paragraph and before the paragraph that
introduces $\boldsymbol\varepsilon_\infty$ and the derivation of (A4.6). The sentences just before: "In an LO wave the dipoles point along $\bq$. Their heads
converge on some planes perpendicular to $\bq$ and their tails on others, so sheets of positive and negative bound charge appear half a
wavelength apart, like the plates of a stack of capacitors, with a field between them. In a TO wave the dipoles lie side by side within each
plane, the head of one next to the tail of its neighbour, and no net charge appears anywhere. Figure A4.2 shows both."

**Physics.** Two-atom cubic crystal, GaAs numbers: $Z^*=2.07$, $\OV=45.16$ Å³ (draw cells of side $\OV^{1/3}=3.56$ Å), $\varepsilon_\infty=10.9$; relative
displacement amplitude $|\bu_1-\bu_2|=0.01$ Å (fixed; state it). Wave along $x$, phase $\theta=qx-\omega t$ (animation advances $\omega t$; the
frequency is irrelevant for electrostatics).
- Dipole of cell $\ell$ (physical, real part): $\mathbf p_\ell=e Z^*|\bu_1-\bu_2|\cos(qX_\ell-\omega t)\,\hat{\mathbf n}$, with $\hat{\mathbf n}=\hat{\mathbf x}$ (LO) or $\hat{\mathbf y}$ (TO).
- Polarization (A4.5): $P(x)=P_0\cos(qx-\omega t)$, $P_0=eZ^*|\bu_1-\bu_2|/\OV=4.58\times10^{-4}$ e/Å² (0.73 μC/cm²).
- Bound charge $-\partial P_x/\partial x=qP_0\sin(qx-\omega t)$ for LO, zero for TO.
- Field (A4.6): LO $E_x(x)=-(4\pi/\varepsilon_\infty)P(x)$, amplitude $eE_0=4\pi e^2Z^*|\bu_1-\bu_2|/(\OV\varepsilon_\infty)=7.61\times10^{-3}$ eV/Å (0.761 MV/cm);
  TO zero.
- Electron potential energy (A4.11): $V^{\mathrm{lr}}(x)=-eA_0$ with $E_x=-\partial A_0/\partial x$, so $V^{\mathrm{lr}}(x)=-(eE_0/q)\sin(qx-\omega t)$; amplitude
  $eE_0/q=eE_0\lambda/2\pi$: 121 meV at $\lambda=100$ Å, 24 meV at 20 Å. TO: zero.

**What to draw.** Two stacked panels sharing the $x$ axis (stack naturally on phones).
- (a) Canvas: 5 rows × N cells of dipole arrows (one arrow per cell, at the cell centre; length ∝ $p_\ell$; photon colour); for LO, translucent
  bands (photon colour for +, electron colour or a muted second colour for −) at the bound-charge maxima, labelled "+" and "−"; field arrows
  $\bE$ (thin, between the sheets, pointing from + to −, i.e. against $\bP$) drawn in a strip above the lattice; the direction of $\hat\bq$ marked
  by a bold `q` arrow. For TO: arrows vertical, no bands, field strip empty with the text "no bound charge, $\bE=0$".
- (b) Plot vs $x$ (Å) over 2 wavelengths: $P_x/P_0$ (phonon or photon dashed), $E_x/E_0$, and $V^{\mathrm{lr}}(x)$ in meV on a right-hand scale (or
  normalized with the amplitude in the readout); an electron marker sitting at the minimum of $V^{\mathrm{lr}}$ (electron colour). Fixed y-limits.
- Axis labels: `x (Å)`, `V^{lr} (meV)`; legend entries `P_x`, `E_x`, `V^{lr}`.

**Controls (3).** seg "mode" `LO` / `TO`; range "wavelength $\lambda$" 20–200 Å (step 1 Å; show also in cells); Play/Pause of `EPP.anim`
(the wave travels; starts paused under reduced motion).

**Readouts.** $|\bP|$ (e/Å² and μC/cm²), $e|\bE|$ (eV/Å) and $|\bE|$ (MV/cm), amplitude of $V^{\mathrm{lr}}$ (meV), $q=2\pi/\lambda$ (Å⁻¹). For TO the field and
potential read 0.

**Initial state.** LO, $\lambda=60$ Å (17 cells), paused at a phase where a + sheet sits near the left third: bands, field and potential all
visible; readout 73 meV.

**Caption key points.**
- *What to notice.* The LO wave builds sheets of bound charge half a wavelength apart and a field between them, against $\bP$; the TO wave,
  with the same dipoles turned by 90°, has none. Longer waves leave the field unchanged but raise the potential swing in proportion.
- *Marks.* arrows = cell dipoles $\mathbf p_\ell$; bands = bound charge $-\nabla\cdot\bP$; field arrows; curves in (b); electron marker.
- *Try.* Switch LO/TO; double $\lambda$ and watch $|\bE|$ stay at 0.76 MV/cm while the $V^{\mathrm{lr}}$ amplitude doubles; play the wave and follow the
  electron's preferred position, a quarter wavelength from the dipole maximum.
- *Model.* (A4.5), (A4.6), (A4.11); GaAs $Z^*=2.07$, $\OV=45.16$ Å³, $\varepsilon_\infty=10.9$ (literature values, illustrative); relative displacement 0.01 Å.

**What changes from v1.** v1 was a static SVG with Giustino's notation: label `δp_p = e Z*·Δτ e^{iq·R_p}` (cell index $p$, $\Delta\tau$), $\rho_b$
labels, a potential named φ(x) (φ is the phonon field in v2: use $V^{\mathrm{lr}}$ for the electron energy, or $A_0$ for the potential), superscript
$\epsilon^\infty$, "U = −eφ", no TO comparison. v1's lower panel (one displaced ion, $R^{-2}$ against $R^{-3}$, labelled ∂V/∂τ) is dropped: the text
covers it. The physics of the top panel is kept (quarter-wavelength shift, $\propto1/q$ potential); add the TO mode, the animation and the
GaAs readouts.

---

## A4.3 `fig-lot.html` · `viz-a4-lot` · new

**Title.** "The modes at Γ depend on the direction of approach" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. $\Phi^{\mathrm{dd}}$ (A4.7) depends on the direction $\hat\bq$, not on $q$, so the optical frequencies at $\Gamma$ depend on the direction from which
   $\bq\to0$: the dynamical matrix is nonanalytic there.
2. In a cubic crystal (GaAs) the LO frequency is the same for every direction, but the LO pattern turns with $\hat\bq$.
3. In a uniaxial crystal (bulk hBN) the frequencies change with direction, and the two modes in the plane of $\hat\bq$ and the axis are mixtures,
   neither purely longitudinal nor transverse; along the axis the "LO" mode (103 meV) lies far below the in-plane TO modes (170 meV).

**Include position.** In subsection `a4-direction`, after its first paragraph and before the paragraph on bulk hBN. The sentences just before:
"In a crystal with one special axis, such as bulk hBN, the Born charges and $\boldsymbol\varepsilon_\infty$ differ along and across the axis, and the
frequencies themselves change with the direction of approach (Figure A4.3)."

**Physics.** Two-atom model, optical block at $\bq\to0$, $\hat\bq=(\sin\theta,0,\cos\theta)$, $\theta$ measured from the crystal axis $z$. With
$\mathbf Z^*_1=-\mathbf Z^*_2=\mathrm{diag}(Z_\perp,Z_\perp,Z_\parallel)$, $\boldsymbol\varepsilon_\infty=\mathrm{diag}(\varepsilon_\perp,\varepsilon_\perp,\varepsilon_\parallel)$, short-range TO frequencies
$\omega_\perp$ (in-plane) and $\omega_\parallel$ (axis), solve the $3\times3$ real symmetric eigenproblem for the direction $\bn$ of relative motion:
$$\sum_\beta\Big[\mathrm{diag}(\omega_\perp^2,\omega_\perp^2,\omega_\parallel^2)_{\alpha\beta}+\frac{4\pi e^2}{\OV M_{\mathrm r}}\,\frac{s_\alpha s_\beta}{\hat\bq\cdot\boldsymbol\varepsilon_\infty\cdot\hat\bq}\Big]n_\beta=\omega^2n_\alpha,\qquad s_\alpha=\sum_\gamma\hat q_\gamma Z^*_{1,\gamma\alpha}$$
(this is $\tPhi^{\mathrm{sr}}(0)+\Phi^{\mathrm{dd}}/\sqrt{M_\kappa M_{\kappa'}}$ of (A4.7) projected on the optical eigenvectors $\be_1=\sqrt{M_2/M_{\mathrm c}}\,\bn$, $\be_2=-\sqrt{M_1/M_{\mathrm c}}\,\bn$;
$s_\alpha$ is a code-internal helper, not a displayed symbol; in units eV, Å, amu the bracket is in eV Å⁻² amu⁻¹ and $\hbar\omega=64.654$ meV $\times\sqrt{\omega^2}$). The $y$ mode decouples ($\omega=\omega_\perp$); the $x$–$z$ block is $2\times2$
(`EPP.eigSym2`). Report the "longitudinal character" $(\hat\bq\cdot\bn)^2$ of each mode.
- GaAs (cubic): $\omega_\perp=\omega_\parallel=33.3$ meV, $Z_\perp=Z_\parallel=2.07$, $\varepsilon_\perp=\varepsilon_\parallel=10.9$, $\OV=45.16$ Å³, $M_{\mathrm r}=36.11$ amu → TO 33.3 meV
  (×2), LO 35.9 meV for every $\theta$.
- Bulk hBN (uniaxial; two-atom model per BN pair, literature values, illustrative, Caldwell et al., Nat. Commun. 5, 5221 (2014)):
  $\hbar\omega_\perp=169.9$ meV (1370 cm⁻¹), $\hbar\omega_\parallel=96.7$ meV (780 cm⁻¹), $\varepsilon_\perp=4.87$, $\varepsilon_\parallel=2.95$, $\OV=18.08$ Å³ (per BN pair,
  $a=2.504$ Å, $c=6.661$ Å), $M_{\mathrm r}=6.10$ amu; $Z_\perp=2.795$, $Z_\parallel=0.730$, chosen so that the model reproduces the measured LO
  frequencies 199.6 meV (1610 cm⁻¹, $\theta=90°$) and 102.9 meV (830 cm⁻¹, $\theta=0$). Check values (meV): θ=0: 102.9, 169.9, 169.9;
  45°: 98.5, 169.9, 189.3 with $(\hat\bq\cdot\bn)^2=0.43$, 0, 0.57; 90°: 96.7, 169.9, 199.6.
- Toggle "$\Phi^{\mathrm{dd}}$ off": drop the second term; all optical modes return to the TO values.

**What to draw.** Two panels side by side (`viz-grid side`, stacked on phones).
- (a) Geometry (canvas): the $x$–$z$ plane of reciprocal space near $\Gamma$: a circle of small $|\bq|$, the crystal axis $z$ vertical (labelled
  "axis" or `c`), a draggable bold arrow `q̂` at angle $\theta$; at the centre, the relative-motion directions $\bn$ of the three optical modes as
  double-headed arrows (each in the colour of its curve in (b)), the out-of-plane $y$ mode as ⊙. Optionally, beside the circle, a small sketch
  of the two sublattices moving along each $\bn$. Drag on the circle or use arrow keys to change $\theta$.
- (b) Plot: $\hbar\omega$ (meV) of the three optical modes against $\theta$ from 0° to 90° (lines; the mixed modes coloured by their longitudinal
  character, or two fixed colours plus a dashed $y$ mode), a vertical marker at the current $\theta$, dashed horizontal guides at the four TO/LO
  end values. Fixed y-limits per preset (GaAs 30–40 meV; hBN 90–210 meV).

**Controls (3).** select "crystal": `GaAs (cubic)` / `bulk hBN (uniaxial)`; range "direction θ of $\hat\bq$ from the axis" 0–90° (step 1°);
toggle "long-range force constants $\Phi^{\mathrm{dd}}$".

**Readouts.** For the current $\theta$: the three $\hbar\omega$ (meV) and their $(\hat\bq\cdot\bn)^2$; the angle between $\bn$ and $\hat\bq$ for the upper
mixed mode.

**Initial state.** bulk hBN, $\theta=45°$, $\Phi^{\mathrm{dd}}$ on: modes at 98.5, 169.9, 189.3 meV, the two mixed ones with longitudinal character
0.43 and 0.57, arrows visibly not along $\hat\bq$.

**Caption key points.**
- *What to notice.* In GaAs the LO frequency is the same whichever way $\bq$ approaches $\Gamma$, and only its pattern turns with $\hat\bq$; in hBN the
  frequencies themselves move, from 103 to 97 meV and from 170 to 200 meV, as $\hat\bq$ tilts from the axis into the layers.
- *Marks.* arrow = $\hat\bq$; double arrows = directions of relative motion $\bn$ of the modes; ⊙ = mode along $y$; curves = mode energies; guides
  = TO and LO values along and across the axis.
- *Try.* Sweep $\theta$ in hBN and watch the upper mode change from transverse to longitudinal; switch to GaAs and see flat lines; turn
  $\Phi^{\mathrm{dd}}$ off and see every splitting disappear.
- *Model.* (A4.7) projected on the optical modes (as in the derivation of (A4.8)); parameters above, literature values, illustrative; the Born
  charges of the hBN model are fitted to the measured LO frequencies; two-atom model of bulk hBN, electrostatic limit ($q\gg q_{\mathrm c}$).

---

## A4.4 `fig-frohlich.html` · `viz-a4-frohlich` · port of v1 `viz-a-polar-frohlich` (v1 file `doc/epp-primer/sections/A04-polar.html`, figure + the first IIFE of its script)

**Title.** "Fröhlich vertex explorer" `<span class="tag">interactive</span>` (keep).

**What the reader must learn.**
1. $|g^{\mathrm{lr}}|$ grows as $1/q$ and beats the acoustic deformation-potential coupling ($\propto q^{1/2}$) at long wavelengths (GaAs: below 0.37 Å⁻¹).
2. Two routes give the same LO coupling: dielectric (A4.15), from $\varepsilon_\infty$, $\varepsilon_{\mathrm s}$, $\omega_{\mathrm{LO}}$; and Born charge, (A4.12) reduced, with the LO frequency from (A4.8).
   $Z^*\to0$ kills the Born-route vertex and its splitting together; $\varepsilon_{\mathrm s}\to\varepsilon_\infty$ does the same for the dielectric route.
3. $\alpha_{\mathrm F}$ ranges from weak (GaAs 0.07) to intermediate (LiF, NaCl ≈ 5 with $m^*=m_{\mathrm e}$).
4. Scales: light cone $q_{\mathrm c}\sim10^{-5}$–$10^{-4}$ Å⁻¹ ≪ $1/r_{\mathrm p}\sim10^{-2}$–$10^{-1}$ Å⁻¹ ≪ zone (≈1 Å⁻¹).

**Include position.** End of subsection `a4-alpha` ("Fröhlich's coupling constant"), after its last paragraph. The sentences just before: "When
$\alpha_{\mathrm F}$ exceeds about 6, the electron digs a potential well in the lattice deep enough to trap itself. The review calls the $1/q$ coupling
'the physical basis for the phenomenon of electron self-trapping in polarons'. Figure A4.4 puts these numbers side by side."

**Physics (Gaussian; all in eV, Å, amu).**
- Dielectric route (A4.15): $|g^{\mathrm{lr}}|\,q=\big[\frac{4\pi e^2}{\OV}\frac{\hbar\omega_{\mathrm{LO}}}{2}(\frac1{\varepsilon_\infty}-\frac1{\varepsilon_{\mathrm s}})\big]^{1/2}$ (eV Å⁻¹), with $4\pi e^2=180.95$ eV Å.
- TO from (A4.10): $\omega_{\mathrm{TO}}=\omega_{\mathrm{LO}}(\varepsilon_\infty/\varepsilon_{\mathrm s})^{1/2}$.
- Born route: $\hbar^2\omega^{2}_{\mathrm{LO},Z}=\hbar^2\omega^2_{\mathrm{TO}}+64.654^2\,\frac{4\pi e^2Z^{*2}}{\OV M_{\mathrm r}\varepsilon_\infty}$ [meV²] from (A4.8); vertex (step 2 of the Fröhlich
  derivation) $|g|\,q=\frac{4\pi e^2Z^*}{\OV\varepsilon_\infty}\big(\frac{\hbar}{2M_{\mathrm r}\omega_{\mathrm{LO},Z}}\big)^{1/2}$ with $\hbar/(2M\omega)=2.0901/(M\,\hbar\omega[\text{meV}])$ Å².
- $Z^*$ implied by the dielectric data: $Z^{*2}=(\omega^2_{\mathrm{LO}}-\omega^2_{\mathrm{TO}})\,\OV M_{\mathrm r}\varepsilon_\infty/(4\pi e^2)$ (in consistent units).
- Acoustic comparison: $|g_{\mathrm{ac}}|=E_{\mathrm{def}}\,(\hbar q/2M_{\mathrm c}v_{\mathrm s})^{1/2}$ (for GaAs $\hbar/2M_{\mathrm c}v_{\mathrm s}=4.22\times10^{-4}$ Å).
- $\alpha_{\mathrm F}=(1/\varepsilon_\infty-1/\varepsilon_{\mathrm s})\sqrt{13.606\,(m^*/m_{\mathrm e})/\hbar\omega_{\mathrm{LO}}[\text{eV}]}$ (A4.16); $1/r_{\mathrm p}=\sqrt{(m^*/m_{\mathrm e})\,\hbar\omega_{\mathrm{LO}}/3.80998}$ Å⁻¹;
  $q_{\mathrm c}=\sqrt{\varepsilon_\infty}\,\hbar\omega_{\mathrm{LO}}/1973.27$ Å⁻¹ (A4.19).
- Presets (literature values, illustrative; keep v1's): ($\varepsilon_\infty$, $\varepsilon_{\mathrm s}$, $\hbar\omega_{\mathrm{LO}}$ meV, $\OV$ Å³, $Z^*$, $m^*/m_{\mathrm e}$, $M_{\mathrm r}$ amu, $M_{\mathrm c}$ amu,
  $v_{\mathrm s}$ m/s, $E_{\mathrm{def}}$ eV): GaAs (10.9, 12.9, 36.2, 45.16, 2.07, 0.067, 36.11, 144.64, 5200, 7); GaN, wurtzite treated as isotropic,
  one formula unit (5.35, 9.5, 92, 22.8, 2.65, 0.2, 11.66, 83.73, 8000, 8); MgO (2.96, 9.8, 89, 18.7, 1.98, 0.35, 9.65, 40.30, 9000, 5); LiF
  (1.93, 9.0, 82, 16.3, 1.05, 1.0, 5.08, 25.94, 6600, 5); NaCl (2.34, 5.9, 32.7, 44.9, 1.12, 1.0, 13.95, 58.44, 4500, 5).
- Check values: GaAs $|g|q=0.0321$ (dielectric) and 0.0305 (Born) eV Å⁻¹, TO 33.28 meV, LO from $Z^*$ 35.93 meV, implied $Z^*=2.19$,
  $\alpha_{\mathrm F}=0.071$, $1/r_{\mathrm p}=0.0252$ Å⁻¹, $q_{\mathrm c}=6.06\times10^{-5}$ Å⁻¹, LO/acoustic crossing 0.37 Å⁻¹. LiF: 0.430 and 0.429 eV Å⁻¹, $\alpha_{\mathrm F}=5.24$.
  MgO: 0.319 eV Å⁻¹, $\alpha_{\mathrm F}=1.72$. NaCl: 0.130 eV Å⁻¹, $\alpha_{\mathrm F}=5.26$. GaN: 0.173 eV Å⁻¹, $\alpha_{\mathrm F}=0.44$.

**What to draw.** Keep v1's layout: (a) log–log $|g|$ (meV) against $q$ (Å⁻¹) from $10^{-6}$ to 1, y from $10^{-1}$ to $10^{9}$ meV: solid LO
(dielectric route), dashed LO (Born route), dotted acoustic; shaded light cone $q<q_{\mathrm c}$; vertical line at $1/r_{\mathrm p}$; slope labels `∝ q^{−1}`,
`∝ q^{1/2}`. (b) bars: TO and LO from the dielectric data (LST) and from the Born charge. Legend labels in words plus v2 symbols, e.g.
"LO, Fröhlich (A4.15): ε_∞, ε_s", "LO, from Z^* and ε_∞", "acoustic, deformation potential E_{def}"; vline label `1/r_p`; band label "light cone".
Bar ticks "from ε_s/ε_∞" and "from Z^*, ε_∞".

**Controls (4; v1 had 7).** select "material" (5 presets); range $\varepsilon_\infty$ 1–15; range $\varepsilon_{\mathrm s}$ ($\ge\varepsilon_\infty$, up to 20); range $Z^*$ 0–3.
Fix $\hbar\omega_{\mathrm{LO}}$, $\OV$, $E_{\mathrm{def}}$, $m^*$, $M_{\mathrm r}$, $M_{\mathrm c}$, $v_{\mathrm s}$ per preset and show them in a small parameter line.

**Readouts.** $\hbar\omega_{\mathrm{TO}}$ (from (A4.10)), $\hbar\omega_{\mathrm{LO}}$ from $Z^*$, $Z^*$ implied by the dielectric data, $|g^{\mathrm{lr}}|\,q$ for both routes (eV Å⁻¹),
$\alpha_{\mathrm F}$ with regime (weak < 1 ≤ intermediate < 6 ≤ strong), $1/r_{\mathrm p}$, $q_{\mathrm c}$ (Å⁻¹, scientific notation with MathJax or rich text).

**Initial state.** GaAs preset.

**Caption key points.**
- *What to notice.* The LO coupling grows as $1/q$ and is larger than the acoustic one below 0.37 Å⁻¹ in GaAs; the two LO routes nearly coincide
  for every preset; the light cone lies three to four decades below $1/r_{\mathrm p}$.
- *Marks.* curves, bars, shading, vline as above.
- *Try.* Move $Z^*$ alone: only the Born route responds, and at $Z^*=0$ its vertex and its splitting vanish together; set $\varepsilon_{\mathrm s}=\varepsilon_\infty$: the
  dielectric route vanishes; switch to LiF and see $\alpha_{\mathrm F}$ jump to about 5.
- *Model.* (A4.15), (A4.8), (A4.10), (A4.12) reduced, (A4.16), (A4.19); acoustic $|g|=E_{\mathrm{def}}q(\hbar/2M_{\mathrm c}\omega)^{1/2}$ with $\omega=v_{\mathrm s}q$; presets
  literature values, illustrative; $\OV$ is the volume of a two-atom cell (one formula unit for wurtzite GaN); in 3D the $1/q$ divergence is
  integrable because $\dd^3q$ brings $q^2\dd q$.

**What changes from v1.** v1 used SI ($e^2/\varepsilon_0=180.95$ eV Å is numerically the same as $4\pi e^2$ Gaussian: keep the number, change the
formulas and words), $\epsilon^\infty$/$\epsilon^0$, μ for the reduced mass, $\Omega$, "$N_p\to1$", legend labels with source equation numbers ("(G.7)",
"(G.180)", "(G.6)"), $E_1$, `q_p`, the bare symbol α, seven controls. Replace all with the v2 symbols above, use (A4.x) references, rename the
classes `a4-frohlich-…`, rewrite the caption in the four-part structure. The physics and the numbers stay (the code's `E2 = 180.95`,
`HB2 = 4.1803e-3` eV² per amu Å⁻² usage is correct); recheck the Born-route formula against step 2 of the Fröhlich derivation.

---

## A4.5 `fig-wannier.html` · `viz-a4-wannier` · port of v1 `viz-a-polar-wannier` (v1 file `doc/epp-primer/sections/A04-polar.html`, figure + the second IIFE of its script)

**Title.** "Why polar crystals break Fourier interpolation: a one-dimensional lab" `<span class="tag">interactive · model</span>`.

**What the reader must learn.**
1. Fourier interpolation from a coarse grid works when the real-space coefficients decay inside the Wigner–Seitz cell of the grid.
2. A $1/q$ singularity has coefficients that never decay (the 1D analogue of the $1/R^2$ dipole tail): the interpolant rings near every $\Gamma$,
   and its error falls only slowly with $N_{\mathrm c}$.
3. Subtracting the analytic long-range part, interpolating the smooth remainder and adding the long-range part back fixes it: $N_{\mathrm c}=8$ is
   then exact to about $10^{-5}$.
4. The result does not depend on the cutoff $r_{\mathrm c}$ once the coarse grid resolves the remainder; a large $r_{\mathrm c}$ leaves more of the dipole
   field in $g^{\mathrm{sr}}$.

**Include position.** End of subsection `a4-wannier`, after its last paragraph and before subsection `a4-limits`. The sentences just before:
"For anatase TiO₂ the remainder $g^{\mathrm{sr}}$ stays below about 150 meV and varies smoothly, while the full $g$ diverges on both sides of $\Gamma$; a coarse
$4\times4\times4$ grid then reproduces direct calculations. Phonon dispersions of polar crystals are interpolated in the same way, with $\Phi^{\mathrm{dd}}$
playing the part of $g^{\mathrm{lr}}$. Figure A4.5 reproduces the whole logic in one dimension."

**Physics (keep v1's model; lattice constant $a=1$, arbitrary units, label it a model).**
- "Data": $g(q)=A\sin^2(qa/2)+B\sin^2(qa)+C\sum_n{}'(q+nb)^{-1}$ with $b=2\pi/a$, $A=1$, $B=0.35$; the lattice sum is $\frac{Ca}{2}\cot(qa/2)$; the term
  $q+nb=0$ is excluded (the rule $\bG\neq-\bq$ of (A4.14)), so the sample at $\Gamma$ is the smooth part only.
- Analytic long-range part with cutoff: $g^{\mathrm{lr}}(q)=C\sum_{|n|\le6}{}'\ee^{-r_{\mathrm c}^2(q+nb)^2}/(q+nb)$ (the factor of the text, $\ee^{-r_{\mathrm c}^2|\bq+\bG|^2}$).
- Coarse grid $q_j=2\pi j/(N_{\mathrm c}a)$; coefficients $g(R_\ell)=N_{\mathrm c}^{-1}\sum_j\ee^{-\ii q_jR_\ell}g(q_j)$, $R_\ell=\ell a$, kept inside the Wigner–Seitz cell of the
  grid, $|\ell|\le N_{\mathrm c}/2$, with half weights on its edge for even $N_{\mathrm c}$; interpolant $\sum_\ell\ee^{\ii qR_\ell}g(R_\ell)$: the 1D analogue of
  (A4.18) and its inverse. With subtraction, interpolate $g-g^{\mathrm{lr}}$ and add $g^{\mathrm{lr}}$ back analytically.
- True coefficients of $\frac C2\cot(q/2)$: magnitude $C/2$ for every $\ell\ne0$ (they never decay; a $1/q$ singularity in 1D is a step in $R$).
  The coarse grid aliases them into a ramp $\frac C2(1-2|\ell|/N_{\mathrm c})$ that reaches zero at the cell edge without real decay.

**What to draw.** Keep v1: (a) $g(q)$ over two zones, $q$ in units of $\pi/a$ from −2 to 2, y from −2.2 to 2.8: true $g$ (dense grid, ink),
coarse samples (dots), interpolant (dashed), analytic $g^{\mathrm{lr}}$ when subtraction is on (dotted). (b) $|g(R_\ell)|$ on a log axis ($10^{-9}$ to 3)
against $R_\ell$ (cells, 0–8): full data (filled) and after subtraction (hollow), the true tail $C/2$ (hline), the Wigner–Seitz edge (vline).
Axis labels `q (π/a)`, `g(q) (arb. units)`, `R_ℓ (cells)`, `|g(R_ℓ)|`; legend "interpolated g", "interpolated g^{sr} + g^{lr}", "analytic g^{lr}".

**Controls (4).** range "polarity $C$ (∝ $Z^*/\varepsilon_\infty$)" 0–0.5; range "coarse grid $N_{\mathrm c}$" 4–16; range "cutoff $r_{\mathrm c}$ (units of $a$)" 0.2–1.4
(v1's σ = $r_{\mathrm c}^2$ from 0.05 to 2 $a^2$); toggle "subtract $g^{\mathrm{lr}}$ before interpolating".

**Readouts.** rms error away from $\Gamma$: direct and with $g^{\mathrm{lr}}$ subtracted; the range of $R_\ell$ kept ($|R_\ell|\le\lfloor N_{\mathrm c}/2\rfloor$ cells).

**Initial state.** $C=0.3$, $N_{\mathrm c}=8$, $r_{\mathrm c}=0.71$ ($r_{\mathrm c}^2=0.5$, v1's default), subtraction off: ringing near both Γ points visible.

**Caption key points.**
- *What to notice.* Without subtraction the interpolant rings near Γ and its error falls only about as $1/N_{\mathrm c}$, because the true coefficients
  never decay; after subtraction the remainder decays like a Gaussian and $N_{\mathrm c}=8$ is exact to about $10^{-5}$.
- *Marks.* curves and points as above.
- *Try.* Toggle the subtraction at $N_{\mathrm c}=8$; set $C=0$ (a nonpolar crystal) and see fast convergence without any help; raise $r_{\mathrm c}$ with
  subtraction on and see the remainder's coefficients grow, so that a coarser grid no longer suffices: the independence of $r_{\mathrm c}$ holds once the
  grid resolves the remainder.
- *Model.* the 1D formulas above, analogues of (A4.14), (A4.17), (A4.18); $a=1$, $A=1$, $B=0.35$; arbitrary units.

**What changes from v1.** v1 uses σ for the cutoff (σ is reserved for the conductivity) → $r_{\mathrm c}$ with $\ee^{-r_{\mathrm c}^2(q+nb)^2}$; $g^{\mathcal L}_\sigma$,
$g^{\mathcal S}$ → $g^{\mathrm{lr}}$, $g^{\mathrm{sr}}$; $R_p$ → $R_\ell$; $N_c$ → $N_{\mathrm c}$; "(G.179)", "(G.178)", "footnote 7's cutoff, called $a$ there" → (A4.18) and the
text's $r_{\mathrm c}$; classes `a4-wannier-…`; caption rewritten in the four-part structure. The numerics stay as in v1.
