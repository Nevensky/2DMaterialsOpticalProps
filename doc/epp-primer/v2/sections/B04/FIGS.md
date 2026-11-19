# B.4 figures: specifications for the figure authors

Section file: `v2/sections/B04-fields.html` (section id `b-fields`, prefix `b4`). Read the section first: every
symbol below is defined there, and every equation is cited by its id (render it with `<a class="eqref" href="#…"></a>`).
Four figures, in reading order:

| # | File | Id | Kind |
|---|---|---|---|
| B4.1 | `B04/fig-helmholtz.html` | `viz-b4-helmholtz` | port of `viz-b-fields-helmholtz` |
| B4.2 | `B04/fig-drift.html` | `viz-b4-drift` | port of `viz-b-fields-drift` |
| B4.3 | `B04/fig-coulomb.html` | `viz-b4-coulomb` | new |
| B4.4 | `B04/fig-vertices.html` | `viz-b4-vertices` | new (diagram; reuses the v1 vertex glyphs) |

Common rules (BUILDER_GUIDE §5): CSS classes prefixed `b4-<name>-`; colours read at draw time with `EPP.color(...)`
(electron = fermion lines, phonon, photon, `c4`, `c5`, `muted`, `ink`); canvas text with `EPP.drawRich`, bold vectors
with `{bold:true}`; readouts in MathJax; 390 px without horizontal overflow (controls below plots); keyboard access and
`aria-label`s; no NaN, no console errors; under 15 ms per redraw. Notation: $\lambda$ is the photon **polarization**
index (1, 2 transverse, 3 longitudinal), so never write λ for a wavelength (write $2\pi/q$); $\sigma$ is reserved for
the conductivity; a tilde means mass-scaled, so no $\tilde A$; photon operators are $b$, never $a$; Gaussian units with
$e^2=14.3996$ eV Å, $\hbar c=1973.27$ eV Å, $\hbar=0.6582119569$ eV fs.

Equation ids of the section used below: `b4-eq-fourier` (B4.12, $\bA(\br)=\sum_\bq\bA_\bq\ee^{\ii\bq\cdot\br}$),
`b4-eq-triad` (B4.13, the triad), `b4-eq-Hfield` (B4.14, field energy $H^{\mathrm T}+H^{\mathrm L}$), `b4-eq-commq`
(B4.16), `b4-eq-ladder` (B4.17), `b4-eq-coulomb` (B4.20), `b4-eq-V` (B4.11), `b4-eq-Vp` (B4.24), `b4-eq-psi` (B4.21).

---

## B4.1 `fig-helmholtz.html` · `viz-b4-helmholtz` · port

- **v1 source.** `doc/epp-primer/sections/B04-fields.html`: the `<figure id="viz-b-fields-helmholtz">` element
  (lines 73–79) and the script block "Viz A: Helmholtz splitter" (lines 193–376), plus the helpers at the top of the
  section script (lines 166–191: `dot`, `cross`, `add`, `sub`, `scl`, `nrm`, `f2`, `arrow`). Keep the two-panel
  design (3D triad view with drag-to-rotate; real-space strips), the chip labels and the halo arrows.
- **What the reader must learn.** One Fourier mode of $\bA$ has a vector amplitude $\bA_\bq$ that splits into a
  transverse part $\bA^{\mathrm T}_\bq$ (across $\bq$) and a longitudinal part $\bA^{\mathrm L}_\bq=\hat\bq(\hat\bq\cdot\bA_\bq)$.
  Only the transverse part has a curl, hence a magnetic field, hence the "spring" of (B4.14); the longitudinal part has
  no curl, only a divergence, which for the electric field of the same pattern is a charge density (Gauss's law). The
  triad $\be_1,\be_2,\be_3=\hat\bq$ follows $\bq$.
- **Include position.** In "One Fourier mode, three polarizations" (`b4-modes`), after the Helmholtz paragraph and before
  the h4 "The polarization triad". Sentences just before: "Two consequences carry the whole section. The magnetic field
  comes only from the transverse part, $\bB=\nabla\times\bA^{\mathrm T}$. The charge in Gauss's law sees only the
  longitudinal part, $\nabla\cdot\bE=\nabla\cdot\bE^{\mathrm L}$. Figure B4.1 shows the two parts of one Fourier mode."
  The triad equation (B4.13) comes *after* the figure, so the caption must define $\be_\lambda$ and $A_{\bq\lambda}$ itself.
- **Physics** (arbitrary units, $|\bA_\bq|=1$).
  - One standing mode with a real amplitude: $\bA(\br)=\bA_\bq\cos(\bq\cdot\br)$ (the real field of the pair $\pm\bq$ in (B4.12)).
  - $\hat\bq=(\sin\theta_q\cos\varphi_q,\ \sin\theta_q\sin\varphi_q,\ \cos\theta_q)$.
  - Triad (B4.13): reference axis $\hat\bn=\hat x$ (use $\hat z$ when $|\hat x\cdot\hat\bq|>0.98$);
    $\be_1=\mathrm{normalize}\big(\hat\bn-(\hat\bn\cdot\hat\bq)\hat\bq\big)$, $\be_2=\hat\bq\times\be_1$ (right-handed),
    $\be_3=\hat\bq$.
  - Amplitude from two angles: $\bA_\bq=\cos\vartheta\,\be_3+\sin\vartheta(\cos\psi\,\be_1+\sin\psi\,\be_2)$, where
    $\vartheta$ is the angle between $\bA_\bq$ and $\bq$ and $\psi$ its azimuth about $\bq$ measured from $\be_1$.
  - Components $A_{\bq\lambda}=\be_\lambda\cdot\bA_\bq$; $\bA^{\mathrm T}_\bq=A_{\bq1}\be_1+A_{\bq2}\be_2$,
    $\bA^{\mathrm L}_\bq=A_{\bq3}\hat\bq$.
  - Real space along $\hat\bq$ (coordinate $\xi=\hat\bq\cdot\br$): $\bB=\nabla\times\bA=-\sin(q\xi)\,\bq\times\bA_\bq$,
    magnitude $q|\bA^{\mathrm T}_\bq|$, perpendicular to the plane of $\bq$ and $\bA^{\mathrm T}_\bq$;
    $\nabla\cdot\bA=-q\,A_{\bq3}\sin(q\xi)$. Colour scales in units of $q|\bA_\bq|$.
  - The identity behind the magnetic term of (B4.14): $|\hat\bq\times\bA_\bq|^2=A_{\bq1}^2+A_{\bq2}^2$.
- **What to draw.**
  - (a) 3D oblique projection, drag to rotate the view (yaw, pitch), also with arrow keys when the canvas has focus
    (`tabindex="0"`, `aria-label`). Cartesian axes $x,y,z$ thin and muted; $\bq$ long teal (`c4`) arrow labelled
    **q**; unit arrows **e**$_1$, **e**$_2$, **e**$_3=\hat{\mathbf q}$ (`fg2`, bold labels); the transverse plane as a
    faint photon-tinted disc with a dashed rim; $\bA_\bq$ thick ink arrow labelled **A**$_{\mathbf q}$; its parts
    $\bA^{\mathrm T}_\bq$ (photon, solid) and $\bA^{\mathrm L}_\bq$ (photon, dashed) with dotted construction lines of the
    parallelogram; labels **A**$^{\rm T}_{\mathbf q}$, **A**$^{\rm L}_{\mathbf q}$ (upright T, L). When one part vanishes,
    label the amplitude "**A**$_{\mathbf q}$ = **A**$^{\rm L}_{\mathbf q}$" or "= **A**$^{\rm T}_{\mathbf q}$".
  - (b) Three horizontal strips over two wavelengths along $\hat\bq$, in the plane of $\bq$ and $\bA_\bq$: (1) "total:
    **A** = **A**$^{\rm T}$ + **A**$^{\rm L}$" arrows on a plain background; (2) "transverse **A**$^{\rm T}$: shear wave,
    colour = $B$" with the background coloured by the out-of-plane $B$; (3) "longitudinal **A**$^{\rm L}$: compression
    wave, colour = ∇·**A**". Diverging ramp (photon = positive, electron = negative), arrows with a surface-coloured halo.
    Zero messages "**A**$^{\rm T}$ = 0: **B** = 0 everywhere" and "**A**$^{\rm L}$ = 0: ∇·**A** = 0 everywhere". Axis ticks
    $0,\ \pi/q,\ 2\pi/q,\ 3\pi/q,\ 4\pi/q$ (not λ), axis label "position along $\hat{\mathbf q}$".
- **Controls (5).** $\theta_q$ polar angle of $\bq$ (0–180°, step 1°); $\varphi_q$ azimuth of $\bq$ (0–360°); $\vartheta$
  angle between $\bA_\bq$ and $\bq$ (0–180°); $\psi$ azimuth of $\bA_\bq$ about $\bq$ from $\be_1$ (0–360°); a preset
  button pair "**A** ∥ **q**" ($\vartheta=0$) and "**A** ⊥ **q**" ($\vartheta=90°$), which set the slider exactly.
- **Readouts.** $A_{\bq1}$, $A_{\bq2}$, $A_{\bq3}$; $|\hat\bq\times\bA_\bq|^2$ next to $A_{\bq1}^2+A_{\bq2}^2$ (always equal);
  longitudinal share $A_{\bq3}^2/|\bA_\bq|^2$ in %; the components of $\be_1$, $\be_2$, $\be_3$.
- **Initial state.** $\theta_q=62°$, $\varphi_q=90°$ ($\bq$ in the $y$–$z$ plane, as in the source's triad figure),
  $\vartheta=40°$, $\psi=30°$: both parts visible, strip (2) and strip (3) both coloured.
- **Caption key points.**
  - *What to notice:* only the transverse part has a curl, so only it carries a magnetic field and costs magnetic (spring)
    energy; the longitudinal part has a divergence and nothing pulls it back.
  - *Marks:* colours and line styles above; the disc is the plane spanned by $\be_1,\be_2$; define
    $A_{\bq\lambda}=\be_\lambda\cdot\bA_\bq$ and $\lambda=1,2$ transverse, $3$ longitudinal.
  - *Try:* "**A** ∥ **q**": $B$ vanishes everywhere, only the compression strip remains; "**A** ⊥ **q**": a pure shear
    wave with no divergence, the shape of a photon; move $\theta_q,\varphi_q$: the triad follows $\bq$ and the readout
    identity holds at every angle.
  - *Model:* one mode $\bA_\bq\cos(\bq\cdot\br)$, arbitrary units; the triad of (B4.13) built from the reference axis
    $x$, right-handed (the handedness is a convention); $\bB$ and $\nabla\cdot\bA$ as above; the identity is the magnetic
    term of (B4.14).
- **What must change from v1.**
  - Notation: v1's $\tilde A_\alpha$, $\tilde{\mathbf A}$ (a tilde means mass-scaled in v2) become $A_{\bq\lambda}$ and
    $\bA_\bq$; "A_t", "A_l" become $\bA^{\mathrm T}$, $\bA^{\mathrm L}$; the index α of the polarizations becomes λ.
  - The strip axis used λ for the wavelength; λ is now the polarization index, so use multiples of $\pi/q$.
  - The triad: v1 followed the source's left-handed $\be_2=-\hat\bq\times\be_1$; use $\be_2=\hat\bq\times\be_1$.
  - Five sliders plus two buttons become four sliders plus a preset pair, with the angle parametrization above, so the
    presets are exact.
  - Caption: rewrite in What to notice / Marks / Try / Model; replace the (E.xxx) references by eqrefs to
    `b4-eq-triad` and `b4-eq-Hfield`.
  - Ids and classes: `viz-b4-helmholtz`, classes `b4-helmholtz-…`; readout labels in MathJax; canvas keyboard rotation.

---

## B4.2 `fig-drift.html` · `viz-b4-drift` · port

- **v1 source.** Same v1 file: `<figure id="viz-b-fields-drift">` (lines 93–101) and the script block "Viz B: oscillator
  vs drift" (lines 378–502). Keep the three plots (two phase-space panels and the energy-versus-time plot), the flow
  arrows, the quantum-state overlay and the animation through `EPP.anim`.
- **What the reader must learn.** A transverse mode is an oscillator: its phase point circles, and its energy sloshes
  between electric and magnetic parts, changing in quanta $\hbar\omega_q=\hbar cq$. The longitudinal mode has no
  spring: its momentum (minus the longitudinal electric field) is constant and its coordinate $A_{\bq3}$ drifts; with no
  charge it does not move at all. Gauss's law fixes its momentum by the charge, so its energy is purely electric, the
  Coulomb energy of the next subsection.
- **Include position.** End of "The field energy, mode by mode" (`b4-energy`), before the h3 "Transverse modes are
  photons". Sentences just before: "For $\lambda=3$ only the kinetic term is there. The longitudinal coordinate moves like
  a free particle: its momentum, the longitudinal electric field, stays constant unless charges move, and the coordinate
  drifts at a constant rate, $\partial_tA_{\bq3}=-cE_{\bq3}$. Figure B4.2 shows both motions in phase space."
  The source callout after (B4.19) also refers to this figure ("the drifting particle").
- **Physics.**
  - Dimensionless quadratures of one real mode, from (B4.17): $X=A_{\bq\lambda}\big/\sqrt{2\pi\hbar c^2/\OVSC\omega_q}$,
    $P=-E_{\bq\lambda}\big/\sqrt{2\pi\hbar\omega_q/\OVSC}$, so $b=(X+\ii P)/2$ and $[X,P]=2\ii$.
  - Energies from (B4.14): transverse $H_{\mathrm T}=\frac{\hbar\omega_q}{4}(X^2+P^2)$, electric part
    $\frac{\hbar\omega_q}4P^2$, magnetic part $\frac{\hbar\omega_q}4X^2$; longitudinal $H_{\mathrm L}=\frac{\hbar\omega_q}{4}P_3^2$.
    For $\lambda=3$ the same scale $\omega_q=cq$ is used only so that both panels share units (a convention; nothing in
    $H^{\mathrm L}$ fixes it).
  - Motion: transverse $X(t)=X_0\cos\omega_qt+P_0\sin\omega_qt$, $P(t)=P_0\cos\omega_qt-X_0\sin\omega_qt$ (clockwise);
    longitudinal $P_3(t)=P_0$, $X_3(t)=X_0+\omega_qP_0\,t$ (this is $\partial_tA_{\bq3}=-cE_{\bq3}$).
  - Gauss's law: $E_{\bq3}=-4\pi\ii\rho_\bq/q$, so $P_3$ is set by the charge; $P_0=0$ means no charge.
  - Quantum overlay: vacuum $\langle X^2\rangle=\langle P^2\rangle=1$ (unit circle) and a coherent state (disc of radius 1
    about the classical point) for the transverse mode; for the longitudinal mode a narrow horizontal band at $P_3=P_0$
    (Gauss's law makes $P_3$ sharp; $X_3$ is a gauge direction, undetermined).
  - Numbers: $\omega_q=cq$, $\hbar\omega_q=1973.27\ \text{eV Å}\times q$; period $2\pi/\omega_q$ in fs; wavelength $2\pi/q$
    in µm. At $q=10^{-4}$ Å⁻¹: 0.197 eV, 6.3 µm, 21 fs.
- **What to draw.** Two square phase-space plots, $X,P\in[-4.5,4.5]$, fixed ticks: (a) transverse, axes "$X\propto
  A_{\mathbf q1}$ (vector potential)" and "$P\propto-E_{\mathbf q1}$ (electric field)", flow arrows and constant-energy
  circles (muted), the orbit and the moving point (photon colour), title "transverse, λ = 1, 2" and
  "$H_{\rm T}=\hbar\omega_q(X^2+P^2)/4$"; (b) longitudinal, axes "$X_3\propto A_{\mathbf q3}$", "$P_3\propto-E_{\mathbf q3}$",
  horizontal flow arrows and level lines, the drifting point with a dashed trail, "runs off →" when it leaves the frame,
  title "longitudinal, λ = 3" and "$H_{\rm L}=\hbar\omega_qP_3^2/4$ (no spring)"; message "$P_3=0$ (no charge): $X_3$ is
  pure gauge" when $P_0=0$. Below: energies (in units of $\hbar\omega_q$) against time (fs) over a fixed window: transverse
  electric (photon), magnetic (`c5`), total (ink), longitudinal electric only (photon, dashed, legend "L: electric only
  (= Coulomb)"); a vertical marker at the current time.
- **Controls (4 + Play).** Photon wavevector $q$, log slider $10^{-0.5}$–$10$ in units of $10^{-4}$ Å⁻¹; initial $X_0$
  ($-2$ to $2$); initial $P_0$ ($-2$ to $2$; for $\lambda=3$ it is set by the charge); toggle "show quantum states";
  Play/Pause from `EPP.anim` (one transverse period every 3 s of real time).
- **Readouts.** $\hbar\omega_q$ (meV or eV), $2\pi/q$ (µm), period $2\pi/\omega_q$ (fs), mean photon number of the coherent
  state $\bar n=(X_0^2+P_0^2)/4$, $E_{\mathrm T}=\hbar\omega_q(\bar n+\frac12)$, $H_{\mathrm L}=\hbar\omega_qP_0^2/4$.
- **Initial state.** $q=1\times10^{-4}$ Å⁻¹ (mid-infrared, $\hbar\omega_q=197$ meV), $X_0=-1.5$, $P_0=1$, quantum states on,
  animation running (paused under reduced motion).
- **Caption key points.**
  - *What to notice:* the transverse point circles and its energy sloshes between electric and magnetic parts, an
    oscillator quantized in steps of $\hbar\omega_q$; the longitudinal point never turns back, and with $P_0=0$ it does
    not move. Its energy is purely electric, constant and set by the charge: the Coulomb energy.
  - *Marks:* panels, colours, the vacuum circle, the coherent-state disc, the Gauss band.
  - *Try:* set $P_0=0$ (no charge: the longitudinal point stops, the transverse one still circles); raise $q$ (faster
    rotation and larger quanta; the drift keeps its form); switch the quantum states on and off.
  - *Model:* quadratures from (B4.17), energies from (B4.14), Gauss's law $E_{\bq3}=-4\pi\ii\rho_\bq/q$; the λ = 3 panel
    uses the transverse scale only for plotting.
- **What must change from v1.**
  - Notation: v1's $\tilde A_\alpha$, $\tilde p_{A\alpha}$ and the index α become $A_{\bq\lambda}$, $E_{\bq\lambda}$ and λ
    (the source's $\tilde p_{A3}=\ii\rho_\bq/cq$ is our $E_{\bq3}=-4\pi\ii\rho_\bq/q$).
  - The readout "λ = 2π/q" must become "$2\pi/q$": λ is the polarization index.
  - References (E.204), (E.208)–(E.209), (E.235′) become eqrefs to `b4-eq-Hfield`, `b4-eq-ladder`, `b4-eq-coulomb`.
  - Caption in the four-part structure; ids and classes `viz-b4-drift`, `b4-drift-…`; the Play button host class renamed.

---

## B4.3 `fig-coulomb.html` · `viz-b4-coulomb` · new

- **What the reader must learn.** The energy stored in the field of static charges, which is purely longitudinal, is the
  Coulomb energy: two self-energies plus an interaction. The interaction is the cross term $\bE_1\cdot\bE_2/4\pi$ of the
  field energy density, spread through space. Between two opposite charges the cross term is positive, yet its integral,
  $-e^2/d$, is negative: the energy saving sits in the large region outside, where the two fields cancel. The
  self-energy grows as the charge shrinks (infinite for point charges), which is why the self terms are removed from
  (B4.20).
- **Include position.** In "The longitudinal mode is the Coulomb interaction" (`b4-coulomb`), after the paragraph on
  $E_{\bq3}=-4\pi\ii\rho_\bq/q$, before the h4 "From field energy to Coulomb energy". Sentences just before: "It is the field
  that a static charge distribution produces: a static field has no curl, so it is purely longitudinal, and the energy
  stored in it is the familiar electrostatic energy of the charges. Figure B4.3 shows where that energy sits for two
  charges." The section later says: "in Figure B4.3 the self-energy of a smeared charge grows as its width shrinks".
- **Physics** (Gaussian units, $e^2=14.3996$ eV Å; lengths in Å, energies in eV).
  - Two Gaussian charge clouds $Q_1$, $Q_2$ (in units of $e$) of width $w$, centred at $z=\mp d/2$ on the $z$ axis:
    $\rho_i(\br)=Q_i\,(2\pi w^2)^{-3/2}\exp\!\big(-|\br-\br_i|^2/2w^2\big)$. (Do not call the width σ; do not call it $s$,
    which labels electrons in B.2 and B.4.)
  - Field of one cloud: $\bE_i=Q_ie\,\frac{\br-\br_i}{|\br-\br_i|^3}\Big[\mathrm{erf}\Big(\frac{|\br-\br_i|}{\sqrt2\,w}\Big)-\sqrt{\tfrac2\pi}\,\frac{|\br-\br_i|}{w}\,\ee^{-|\br-\br_i|^2/2w^2}\Big]$;
    the total field $\bE=\bE_1+\bE_2$ is curl-free (longitudinal).
  - Energy density $\frac{E^2}{8\pi}=\frac{E_1^2}{8\pi}+\frac{E_2^2}{8\pi}+\frac{\bE_1\cdot\bE_2}{4\pi}$.
  - Geometry of the cross term: outside the clouds $\bE_1\cdot\bE_2\propto Q_1Q_2\,(\br-\br_1)\cdot(\br-\br_2)$, which
    vanishes on the sphere whose diameter joins the two centres (Thales' theorem). For opposite charges the cross term is
    positive inside that sphere and negative outside; for like charges the reverse.
  - Integrals (analytic): self-energy $U_i=Q_i^2e^2/(2\sqrt\pi\,w)$; interaction
    $U_{12}=\frac{Q_1Q_2e^2}{d}\,\mathrm{erf}\big(\frac{d}{2w}\big)\to\frac{Q_1Q_2e^2}{d}$ for $d\gg w$ (Coulomb's law). The
    total field energy is $U_1+U_2+U_{12}$, the classical form of (B4.20).
  - Check values: $Q_1=+1$, $Q_2=-1$, $d=4$, $w=0.5$: $U_{12}=-3.600$ eV, $U_i=8.12$ eV.
  - Optional numerical check (only if under 15 ms): integrate $\bE_1\cdot\bE_2/4\pi$ on an axisymmetric $(\rho,z)$ grid
    with weight $2\pi\rho$, using a mapped radial coordinate to reach large distances; it must reproduce $U_{12}$ within
    1 %. Otherwise show the analytic values only.
- **What to draw.**
  - (a) Map of the plane containing both charges, fixed limits $x\in[-8,8]$ Å, $z\in[-10,10]$ Å (aspect 4:5): heat map
    of the cross-term density $\bE_1\cdot\bE_2/4\pi$ (diverging: photon colour positive, electron colour negative,
    symmetric signed-log or clipped scale with a colour bar), or of the total $E^2/8\pi$ (sequential, log scale) when
    selected. Draw the Thales circle (diameter from one centre to the other) as a thin dashed line. Overlay 8–12 field lines of the total field per charge (RK integration from a ring around each positive
    charge, ending on negative charges or at the frame), muted ink. Charges as discs of radius $w$ labelled $Q_1$,
    $Q_2$ with their sign; a dimension line with the label $d$.
  - (b) Plot of energy against separation, $d\in[0,12]$ Å, fixed $y$ range $[-15,15]$ eV: $U_{12}(d)$ solid (photon
    colour), $Q_1Q_2e^2/d$ dashed (muted), a dot at the current $d$; axis labels "separation $d$ (Å)" and "interaction
    energy (eV)".
- **Controls (4).** Separation $d$ (1–12 Å, step 0.1); second charge $Q_2$ (segmented: $-2e$, $-e$, $+e$, $+2e$;
  $Q_1=+e$ fixed); width $w$ (0.2–1.5 Å); map: "cross term $\mathbf E_1\cdot\mathbf E_2/4\pi$" or "total $E^2/8\pi$".
- **Readouts.** $U_{12}$ (eV); Coulomb's law $Q_1Q_2e^2/d$ (eV); self-energy per charge $U_1$, $U_2$ (eV); total field
  energy (eV); if implemented, the numerical integral of the cross term.
- **Initial state.** $Q_1=+e$, $Q_2=-e$, $d=4$ Å, $w=0.5$ Å, cross-term map: positive (rose) between the charges, negative
  (blue) lobes around and beyond them, and $U_{12}=-3.60$ eV in the readout.
- **Caption key points.**
  - *What to notice:* the interaction energy of two charges is the cross term of their field energy; between opposite
    charges it is positive, yet its integral $-e^2/d$ is negative, because the saving sits outside, where the fields
    cancel.
  - *Marks:* colour scale and sign, field lines, discs of radius $w$, the dashed Thales circle on which the cross term
    changes sign, the energy curves.
  - *Try:* switch $Q_2$ to $+e$ (the map changes sign and $U_{12}$ becomes positive); shrink $w$ (the self-energies grow as
    $1/w$, while $U_{12}$ stays at Coulomb's value); bring the clouds closer than about $2w$ ($U_{12}$ falls below
    Coulomb's law because they overlap).
  - *Model:* Gaussian charges, analytic fields and energies above, Gaussian units; static fields are longitudinal, so the
    field energy is $H^{\mathrm L}$ of `b4-eq-coulomb` written for two classical charges; the self terms are what the
    section removes by normal ordering and by dropping on-site dipole terms.

---

## B4.4 `fig-vertices.html` · `viz-b4-vertices` · new (diagram)

- **Starting point.** v1 drew these glyphs inline in the section, not as a figure: the `div.fd` rows
  `#b-fields-fd-couplings` (v1 lines 58–63) and `#b-fields-fd-coulomb` (lines 131–135), drawn by the script block
  "Vertex glyphs" (lines 504–534) with `EPP.fd.draw`. Reuse those item lists, regrouped into one figure.
- **What the reader must learn.** Every term of $V'$ (B4.24) is a vertex: two electron–photon vertices (paramagnetic, one
  photon; diamagnetic, two photons at one point), two ion–photon vertices (paramagnetic, a phonon turns into a photon;
  diamagnetic, two photons at an ion site and no phonon line), and three Coulomb vertices from (B4.20)
  (electron–electron, the bare $v$; electron–ion, the bare Fröhlich vertex; ion–ion, the dipole–dipole interaction).
  No glyph joins an electron line directly to a phonon line: electrons and phonons interact only through the field.
- **Include position.** In "Electrons and ions as operators" (`b4-second`), after the ion momentum operator (B4.22) and
  its explanation, just before the h4 "Free part and interaction". Sentence just before: "Figure B4.4 collects the
  couplings that the next equation contains." The four-couplings paragraph after (B4.11) also cites it.
- **Physics.** The terms of (B4.24) (`b4-eq-Vp`), equivalently the four couplings of (B4.11) (`b4-eq-V`) and the Coulomb
  energy (B4.20), with the current (B4.21).
- **What to draw.** Two labelled rows, which wrap to one column at phone width:
  - Row "Photon vertices (minimal coupling)", 4 glyphs: (1) electron paramagnetic: fermion line in and out (electron
    colour, arrow) with one wavy photon line, label $-\frac1c\,\hat{\mathbf j}^{\rm para}\cdot\hat{\mathbf A}$; (2) electron
    diamagnetic: fermion line with two photon lines from one point, label
    $\frac{e^2}{2m_{\rm e}c^2}\hat n_{\rm e}\hat{\mathbf A}^2$; (3) ion paramagnetic: zigzag phonon line (phonon colour)
    turning into a photon line, label $-\frac ec\frac{\hat{\mathbf p}_{\ell\kappa}}{M_\kappa}\cdot\mathbf Z^{*\mathsf T}_\kappa\cdot\hat{\mathbf A}$;
    (4) ion diamagnetic: two photon lines meeting at a small phonon-coloured disc (the ion site), no phonon line, label
    $\frac{e^2}{2c^2M_\kappa}|\mathbf Z^{*\mathsf T}_\kappa\cdot\hat{\mathbf A}|^2$.
  - Row "Coulomb vertices (longitudinal field)", 3 glyphs with a dashed photon-coloured Coulomb line: (5)
    electron–electron, two fermion lines, label "bare $v$"; (6) electron–ion, a fermion line and a phonon line, label
    "bare Fröhlich vertex"; (7) ion–ion, two phonon lines, label "dipole–dipole".
  - Under each glyph a one-line "becomes" note in muted text: (1) $\Pi^{\rm e}$, current–current part (B.5); (2)
    diamagnetic part of $\Pi^{\rm e}$ (B.5); (3) $\Pi^{\rm ion}$ (B.5), phonon self-energy (B.7); (4) diamagnetic part of
    $\Pi^{\rm ion}$ (B.5); (5) $W^{\rm e}$ (B.9); (6) Fröhlich vertex (B.8); (7) LO–TO splitting (B.7, B.9).
  - Optional and light: hover or keyboard focus on a glyph highlights it and its label; no sliders, so the tag stays
    "diagram". Labels as MathJax in HTML under the SVG glyphs; each glyph has an `aria-label`.
- **Controls.** None (diagram). If the highlight is implemented, glyphs are focusable (`tabindex="0"`).
- **Initial state.** All seven glyphs shown.
- **Caption key points.**
  - *What to notice:* no vertex joins an electron line to a phonon line; electron–phonon coupling needs a photon or a
    Coulomb line in between.
  - *Marks:* straight arrowed lines are electrons, zigzag lines phonons, wavy lines transverse photons, dashed lines the
    instantaneous Coulomb interaction; dots are vertices, the small disc an ion site.
  - *Try:* if interactive, focus a glyph to see its term.
  - *Model:* the terms of `b4-eq-Vp`; the Coulomb pieces from squaring $\rho=-e\,n_{\rm e}+\rho^{\rm ion}$ in
    `b4-eq-coulomb`.
- **What must change from the v1 glyphs.** v1 labels used source notation ("π" for the ion momentum, "Q" for the
  displacement, "$Z^*Z^*/m_a$"); use $\hat{\mathbf p}_{\ell\kappa}$, $\hat{\mathbf u}_{\ell\kappa}$ and
  $\mathbf Z^*_\kappa\mathbf Z^{*\mathsf T}_\kappa/M_\kappa$. Put the glyphs inside one `figure.viz` with id
  `viz-b4-vertices`, classes `b4-vertices-…`, and a four-part caption.
