# B.8 figures: specifications for the figure authors

Section file: `v2/sections/B08-frohlich.html` (id `b-frohlich`, prefix `b8`). Read it before building any figure:
every symbol, equation tag and number below is defined there. Figures appear in this order, so the page numbers them
B8.1–B8.5: `bubble`, `pairing`, `continuity`, `dynamic`, `lightline`.

Common rules (from `v2/build/BUILDER_GUIDE.md` §5 and NOTATION):

- v2 notation only, in canvas text (`EPP.drawRich`), readouts (MathJax `$…$` in HTML), legends and aria-labels.
  Never write source symbols: no `Ω_uc`, `Ω^e`, `σ` as a branch index, `ε` written `\epsilon`, `Π_jj`, `χ^0` for the
  pair factor, `λ_c`, `m_a`, `e_{aνσ}`, `r^0_a`, `F(q,ω)` for the Lindhard function (F is the phonon propagator).
- Labels: upright for words (`ω_{LO}`, `ε_∞`, `k_{TF}`, `v_F`, `q_c` with upright c), bold vectors (`{bold:true}`),
  italic variables. Phonon branch index is ν; Cartesian indices α, β (γ, δ inner); bands n, m.
- Colours read at draw time with `EPP.color(...)`. Suggested roles: photon lines and TO curves `photon`; LO curves `c5`;
  ionic vertices and ionic currents `c4`; electronic vertices, electron lines and carriers `electron`; phonon lines
  `phonon`; references `muted`.
- Units on every axis; fixed axis limits; 2–5 controls with units; no NaN; < 15 ms per redraw; works at 390 px with
  controls below the plots; keyboard access and aria-labels; animations through `EPP.anim`.
- Captions: "What to notice" (1–2 sentences, visible), then `cap-more` with Marks, Try, Model; about 200–260 words in
  total. Reference equations with `<a class="eqref" href="#b8-eq-…"></a>` (ids listed per figure).
- Material data (literature values, illustrative; the same as v1 and Section A.4):

| | $\varepsilon_\infty$ | $\varepsilon_{\mathrm s}$ | $\hbar\omega_{\mathrm{LO}}$ (meV) | $\hbar\omega_{\mathrm{TO}}$ (meV) | $\OV$ (Å³) | $Z^*$ (dimensionless) | $M_{\mathrm r}$ (amu) |
|---|---|---|---|---|---|---|---|
| GaAs | 10.9 | 12.9 | 36.2 | 33.3 | 45.2 | 2.07 | 36.1 |
| MgO | 2.96 | 9.8 | 89.0 | 49.1 | 18.7 | 1.98 | 9.65 |
| NaCl | 2.34 | 5.9 | 32.7 | 20.3 | 44.9 | 1.12 | 13.9 |

- Constants: $e^2=14.3996$ eV Å, so $4\pi e^2=180.95$ eV Å; $\hbar c=1973.27$ eV Å; $\hbar^2/m_{\mathrm e}=7.61996$ eV Å²;
  $\hbar^2/(\mathrm{amu}\,\text{Å}^2)=4.1802$ meV; $\hbar/(2M\omega)=2.0901$ Å²$/(M[\mathrm{amu}]\,\hbar\omega[\mathrm{meV}])$.

---

## 1. `fig-bubble.html` · id `viz-b8-bubble` · NEW (diagram with a counting control)

**Title.** "One bubble between a screened and a bare photon line" `<span class="tag">interactive</span>`.

**What the reader must learn.** The electron-induced phonon self-energy is a single electron–hole bubble between a
dressed photon line and a bare one. Grouping each photon line with the two vertices at its ends defines the screened
vertex $g^{\mathrm{dyn}}$ and the bare vertex $g^{\mathrm{dyn},0}$. Dressing both lines counts a chain of $k$ bubbles
$k$ times; leaving both bare keeps only one bubble. Part A's nonadiabatic phonon self-energy has the same structure.

**Include position.** Section `b8-open` ("Opening the electron bubble"), after the background box "Reading a Feynman
diagram" and before the four steps. The main-text sentences just before it: "The vertex that couples a phonon to one
electronic transition appears when the chain is cut open at one bubble. Figure B8.1 draws the result: a single bubble
between two photon lines, one dressed and one bare."

**Physics.** Equations `b8-eq-open` (B8.3) $\Delta\Pi^{\mathrm{ph}}_{\nu\nu'}=-\sum J^*_{\alpha\nu}D^{\mathrm e}_{\alpha\gamma}\Pi^{\mathrm e}_{\gamma\delta}D^0_{\delta\beta}J_{\beta\nu'}$,
`b8-eq-vertex` (B8.6) for the two vertices, `b8-eq-vform` (B8.7)
$\Delta\Pi^{\mathrm{ph}}_{\nu\nu'}=\frac1{N_\ell}\sum_{\bk nm}g^{\mathrm{dyn}}_{mn\nu}\,\mathcal B_{n\bk,m\bk+\bq}\,g^{\mathrm{dyn},0}_{mn\nu'}$,
and the count `b8-eq-count` (B8.8): with $x=D^0\Pi^{\mathrm e}$ (one bubble on a bare line, a number),
- one bare + one screened: $D^{\mathrm e}\Pi^{\mathrm e}D^0=D^0(x+x^2+x^3+\dots)$, every chain once (correct);
- two screened: $D^{\mathrm e}\Pi^{\mathrm e}D^{\mathrm e}=D^0(x+2x^2+3x^3+\dots)$, chain of $k$ bubbles $k$ times (double counting);
- two bare: $D^0\Pi^{\mathrm e}D^0=D^0x$, only $k=1$ (screening to first order).

**What to draw.** Use `EPP.fd.draw` (SVG, theme-aware) inside the figure; no plot.
- Top panel, the main diagram (it redraws the EPP manuscript's Fig. 14 in v2 notation), left to right on one line:
  phonon line (zigzag, `phonon`) labelled ν; dot (`c4`) labelled $J^*_{\alpha\nu}$; photon line labelled
  $D^{\mathrm e}_{\alpha\gamma}$ (heavy wavy `photon2` when dressed, thin `photon` when bare); dot (`electron`)
  labelled $J_\gamma^{\,n\bk,m\bk+\bq}$; electron–hole bubble of two arrowed electron lines, upper labelled $n\bk$, lower
  $m\bk+\bq$, with $\mathcal B_{n\bk,m\bk+\bq}$ written inside; dot (`electron`) labelled $(J_\delta^{\,n\bk,m\bk+\bq})^*$;
  photon line labelled $D^0_{\delta\beta}$ (thin) or $D^{\mathrm e}_{\delta\beta}$ (heavy), depending on the control;
  dot (`c4`) labelled $J_{\beta\nu'}$; phonon line ν′. Under the left half a bracket labelled
  "$g^{\mathrm{dyn}}_{mn\nu}$ screened" (or "$g^{\mathrm{dyn},0}_{mn\nu}$ bare" when that line is bare); under the right half
  a bracket labelled "$g^{\mathrm{dyn},0}_{mn\nu'}$ bare" (or "screened"). Keep the labels upright for words ("bare",
  "screened"), italic for symbols.
- Middle panel, the expansion: three short rows for chains of $k=1,2,3$ bubbles (bare thin lines and small bubbles
  in a row between the two ionic dots). Beside each row a counter "×1", "×2", "×3" or "×0" showing how many times the
  selected pairing generates that chain. Colour a counter `muted` when it equals 1, `photon` when it is larger
  (double counting), `c4` when it is 0 (missing screening).
- Bottom panel, Part A for comparison: phonon line $\bq\nu$, small dot labelled "bare" (the matrix element of
  $\Delta_{\bq\nu}V^{\mathrm{en}}$), Kohn–Sham bubble, grey disk labelled "screened $g$", phonon line. One line of text:
  "Part A puts the screened end on the right; the order does not matter."
- At 390 px stack the panels; scale the SVG width to the container (no horizontal scroll).

**Controls.** One `seg` control "photon lines": [one bare + one screened] (default) · [two screened] · [two bare]. It
switches the line styles and brackets of the top panel and the counters of the middle panel. Optionally a toggle
"show Part A" (default on).

**Readout.** One line: "chains of 1, 2, 3 bubbles counted 1, 1, 1 times: correct" / "1, 2, 3 times: double counting" /
"1, 0, 0 times: screening missing". Give the readout `aria-live="polite"`.

**Initial state.** One bare + one screened: heavy line on the left, thin on the right, counters ×1 ×1 ×1.

**Caption key points.**
- What to notice: the electrons enter as one bubble between a dressed and a bare photon line; each line with its two
  end vertices is one electron–phonon vertex, screened on the left and bare on the right.
- Marks: zigzag = phonon; thin wavy = free photon line $D^0$; heavy wavy = dressed line $D^{\mathrm e}$; `c4` dots =
  ionic vertices $J_{\alpha\nu}$; `electron` dots = electronic vertices $J_\gamma^{\,n\bk,m\bk+\bq}$; the loop = pair
  factor $\mathcal B$; counters = how often each chain of bubbles is generated.
- Try: switch to "two screened" and watch the counters grow as 1, 2, 3; switch to "two bare" and see the chains with
  two or more bubbles disappear.
- Model: `b8-eq-open`, `b8-eq-vertex`, `b8-eq-vform`, `b8-eq-count`; the Part A diagram is the bare–screened
  structure of Section A.6 (`a6-eq-pina`, an export of A.6; link it with `a.eqref`).

**Notation warnings.** EPP Fig. 14 labels the bubble $\chi^0_{n\mathbf km\mathbf k+\mathbf q}$ and the vertices
$J^{\rm ion}$, $J^{\rm e}$, branch σ: do not reproduce these. EPP draws both photon lines in the same colour; here the
dressed line is heavy.

---

## 2. `fig-pairing.html` · id `viz-b8-pairing` · PORT of `viz-b-frohlich-pairing` (v1 `doc/epp-primer/sections/B08-frohlich.html`)

**Title.** "Counting bubbles: which vertices go into the self-energy" `<span class="tag">interactive</span>`.

**What the reader must learn.** Only one bare and one screened vertex give the measured LO frequency. Two screened
vertices barely screen; two bare vertices overscreen until the lattice becomes unstable.

**Include position.** End of `b8-counting` ("Why one line is dressed and one is bare"), after (B8.9). Sentences just
before: "…the three pairings give $\hbar\omega_{\mathrm{LO}}=35.9$ meV, close to the measured 36.2 meV; 54.1 meV with two
screened vertices; and an imaginary frequency with two bare vertices, a spurious lattice instability that sets in above
$\varepsilon_\infty\approx2.6$. Figure B8.2 plots the three cases against $\varepsilon_\infty$ for three materials."

**Physics.** `b8-eq-pairing` (B8.9): $\omega_{\mathrm{LO}}^2=\omega_{\mathrm{TO}}^2+(\omega_{\mathrm p}^{\mathrm{ion}})^2\,\zeta(\varepsilon)$,
$(\omega_{\mathrm p}^{\mathrm{ion}})^2=4\pi e^2Z^{*2}/(\OV M_{\mathrm r})$, so
$(\hbar\omega_{\mathrm p}^{\mathrm{ion}})^2[\mathrm{meV}^2]=180.95\,Z^{*2}/(\OV[\text{Å}^3]\,M_{\mathrm r}[\mathrm{amu}])\times4.1802\times10^{3}$.
$\zeta=1/\varepsilon$ (one bare + one screened), $(\varepsilon^2-\varepsilon+1)/\varepsilon^2$ (two screened),
$2-\varepsilon$ (two bare). Model: cubic diatomic crystal, LO mode, $c\to\infty$, $q\to0$, static electronic screening;
the bare frequency is $\omega_{\mathrm{TO}}$ (no local-field correction on this route). Check values: GaAs
$\hbar\omega_{\mathrm p}^{\mathrm{ion}}=44.57$ meV; at $\varepsilon=10.9$: 35.93, 54.13 meV, unstable; instability for
$\varepsilon>2+\omega_{\mathrm{TO}}^2/(\omega_{\mathrm p}^{\mathrm{ion}})^2=2.558$. MgO: 128.2 meV; 89.2, 123.2 meV,
unstable (threshold 2.15). NaCl: 39.0 meV; 32.6, 39.5 meV, unstable (threshold 2.27).

**What to draw.** One `EPP.Plot`: x = electronic screening $\varepsilon_\infty$ (dimensionless), fixed range 1–15;
y = $\hbar\omega_{\mathrm{LO}}$ (meV), range $[-0.6W, 1.22W]$ with $W=\sqrt{\omega_{\mathrm{TO}}^2+(\omega_{\mathrm p}^{\mathrm{ion}})^2}$
of the selected material (fixed per material). Imaginary frequencies are plotted as negative, $\mathrm{sign}(\omega^2)\sqrt{|\omega^2|}$,
inside a faint `photon` band below 0 labelled "ω² < 0: lattice unstable". Curves: one bare + one screened (`phonon`,
solid), two screened (`c5`, dashed), two bare (`photon`, dotted). Horizontal lines: bare $\omega_{\mathrm{TO}}$ (`muted`,
labelled "bare ω_{TO}"), measured $\omega_{\mathrm{LO}}$ (`phonon`, dotted, labelled "measured ω_{LO}"). Vertical dotted
line at the material's $\varepsilon_\infty$ (`muted`), vertical cursor line (`ink`) with three points on the curves.
Legend below the plot.

**Controls.** `select` "material": GaAs (default), MgO, NaCl; selecting a material moves the cursor to its
$\varepsilon_\infty$. `range` "cursor ε_∞": 1–15, step 0.01. Clicking or dragging on the plot moves the cursor
(`plot.onPointer`).

**Readouts.** $\varepsilon_\infty$ (cursor); $\hbar\omega_{\mathrm{LO}}$ for "bare × screened", "screened × screened",
"bare × bare" (or "unstable (ω² < 0)"); measured $\hbar\omega_{\mathrm{LO}}$; $\hbar\omega_{\mathrm p}^{\mathrm{ion}}$.

**Initial state.** GaAs, cursor at 10.9.

**Caption key points.**
- What to notice: at $\varepsilon_\infty=1$ (no electrons) the three pairings agree; only bare × screened lands on the
  measured LO; two screened vertices give 54 meV instead of 36 meV in GaAs, and two bare vertices make $\omega^2<0$
  above $\varepsilon_\infty\approx2.6$.
- Marks: the three curves, the two horizontal references, the material's $\varepsilon_\infty$, the cursor.
- Try: drag the cursor from 1 to 15 and watch the bare × bare curve cross zero; switch to MgO, whose small
  $\varepsilon_\infty=2.96$ keeps the miscounts smaller but still visible; compare the measured LO line with the solid
  curve at each material's $\varepsilon_\infty$.
- Model: `b8-eq-count`, `b8-eq-pairing`; parameters from the table above (literature values, illustrative); static,
  long-wavelength, electrostatic limit; one LO mode of a cubic diatomic crystal.

**What changes from v1.** v1 axis label "ε^e" and readout keys "ε^e", "ħω_p (unscreened)", legend "(E.412)", caption
symbols $\Omega^0$, $\Delta\Omega^{\rm e}$, $\omega_\sigma$, $\Omega_{\rm uc}$, $\mu$ (reduced mass), $\epsilon^\infty$, the SI
remark with $Z^*_{\rm G}$, and the source tag E.405 must all go. Use $\varepsilon_\infty$ (axis and cursor), $\omega_{\mathrm{TO}}$ for the
bare frequency, $\omega_{\mathrm p}^{\mathrm{ion}}$, $\OV$, $M_{\mathrm r}$, $Z^*$ (dimensionless) and the v2 equation ids. Rename id
and classes (`b8-pairing-…`). The physics and numbers of v1 are right; keep them. Rewrite the caption in the four-part
structure.

---

## 3. `fig-continuity.html` · id `viz-b8-continuity` · NEW (geometric, animated)

**Title.** "Charge and current of one electronic transition" `<span class="tag">interactive</span>`.

**What the reader must learn.** Averaged over a cell, an electronic transition $n\bk\to m\bk+\bq$ is a charge wave
$\propto\ee^{\ii(\bq\cdot\br-\omega_{\mathrm{tr}}t)}$. Its current along $\bq$ is fixed by continuity,
$\bq\cdot\mathbf j=\omega_{\mathrm{tr}}\rho$; its current across $\bq$ is free and moves no charge. A longitudinal (LO)
ionic current talks only to the charge, through the Coulomb line; a transverse (TO) ionic current talks only to the
transverse current, through light.

**Include position.** Section `b8-continuity`, right after its first paragraph and before the background box on
creation and annihilation operators. Sentences just before: "So the longitudinal part of the current vertex is fixed by
the charge of the transition, while the transverse part is free. Figure B8.3 shows the charge and the current of one
transition."

**Physics.** `b8-eq-identity` (B8.11): $\hbar\,\bq\cdot\mathbf j^{\,n\bk,m\bk+\bq}=(E_{m\bk+\bq}-E_{n\bk})\,\rho^{\,n\bk,m\bk+\bq}$,
i.e. $\bq\cdot\mathbf j=\omega_{\mathrm{tr}}\rho$. Draw the cell-averaged (macroscopic) transition densities as running
waves: charge $\rho(\br,t)=\rho_0\cos(\bq\cdot\br-\omega_{\mathrm{tr}}t)$; current
$\mathbf j(\br,t)=\big(j_\parallel\hat\bq+j_\perp\mathbf t\big)\cos(\bq\cdot\br-\omega_{\mathrm{tr}}t)$ with
$j_\parallel=\omega_{\mathrm{tr}}\rho_0/q$ fixed and $j_\perp$ free (slider). Check in the code that
$\partial_t\rho+\nabla\cdot\mathbf j=0$ holds for any $j_\perp$: $\nabla\cdot\mathbf j=-q\,j_\parallel\sin(\cdot)$ and
$\partial_t\rho=\omega_{\mathrm{tr}}\rho_0\sin(\cdot)$. Work in reduced units (wavelength and period of the wave = 1);
no material data are needed. Coupling readout (from `b8-eq-ratio`, `b8-eq-corollary`): an LO ionic current couples to
$\hat\bq\cdot\mathbf j=\omega_{\mathrm{tr}}\rho/q$ with weight $\mathcal K_{\mathrm{LO}}$, independent of $j_\perp$; a TO
ionic current couples to $\mathbf t\cdot\mathbf j=j_\perp$ with weight $\mathcal K_{\mathrm{TO}}$, which is negligible
unless $q\lesssim q_{\mathrm c}$.

**What to draw.** One `EPP.canvas` (aspect about 0.6, min height 240 px). Axes: $\hat\bq$ horizontal (label
"$\hat{\mathbf q}$", bold), $\mathbf t$ vertical (label "$\mathbf t$", bold); about two wavelengths across. Background:
the charge wave as a diverging colour field (`electron` for $\rho\gt0$, `photon` for $\rho\lt0$, transparent at 0),
stripes perpendicular to $\bq$; label "ρ" in a small legend ("charge of the transition"). Arrows on a coarse grid (about
9 × 5) for $\mathbf j(\br,t)$ in `ink`; their horizontal part is $j_\parallel\cos(\cdot)$, their vertical part
$j_\perp\cos(\cdot)$. A dashed vertical line moving with the crests at the phase velocity $\omega_{\mathrm{tr}}/q$. At
the right or below, two small fixed arrows: "LO ionic current" along $\hat\bq$ (`c5`) and "TO ionic current" along
$\mathbf t$ (`photon`), each with its coupling readout.

**Controls.** `range` "transverse current $j_\perp/j_\parallel$": −2 to 2, step 0.05, default 0. `seg` "transition":
[upward, $\omega_{\mathrm{tr}}\gt0$] (default) · [downward, $\omega_{\mathrm{tr}}\lt0$], which reverses the drift of the
stripes and the sign of $j_\parallel$. Play/Pause from `EPP.anim` (calm: one period in about 4 s; paused under reduced
motion, then show $t=0$).

**Readouts.** $j_\parallel=\omega_{\mathrm{tr}}\rho/q$ ("fixed by continuity"); $j_\perp$ ("free"); "charge piled up per
half period: set by $j_\parallel$ only"; "LO couples to: the charge ρ (Coulomb, every q)"; "TO couples to: $j_\perp$
(light, only for $q\lesssim q_{\mathrm c}$)".

**Initial state.** $j_\perp=0$, upward transition, animation running (paused under reduced motion): stripes drift to
the right with horizontal arrows that point forward on positive crests and backward on negative ones.

**Caption key points.**
- What to notice: the stripes of charge move because current flows along $\bq$ into each growing crest; that current is
  fixed by the charge, while current across $\bq$ flows along the crests and changes nothing.
- Marks: colour field = cell-averaged charge of the transition, $\rho\cos(\bq\cdot\br-\omega_{\mathrm{tr}}t)$; arrows =
  current; dashed line = a crest moving at $\omega_{\mathrm{tr}}/q$; the two side arrows = ionic currents of an LO and a
  TO phonon.
- Try: raise $j_\perp/j_\parallel$ to 2: the arrows tilt but the stripes and their motion do not change; reverse the
  transition and watch both the drift and $j_\parallel$ flip.
- Model: `b8-eq-identity` (continuity for one transition), `b8-eq-corollary` (the longitudinal electronic vertex is
  $\omega$ times the charge), `b8-eq-ratio` (LO and TO kernels); reduced units; the cell-periodic structure of the
  Bloch functions is averaged out ($\bG=0$).

---

## 4. `fig-dynamic.html` · id `viz-b8-dynamic` · PORT of `viz-b-frohlich-dynamic` (v1 `doc/epp-primer/sections/B08-frohlich.html`)

**Title.** "Screening the vertex at the phonon frequency, not at zero" `<span class="tag">interactive</span>`.

**What the reader must learn.** In a doped polar semiconductor, static (adiabatic) screening by the carriers kills the
polar vertex at small $q$, while screening at the phonon frequency leaves its $1/q$ form intact for
$q\lesssim\omega/v_{\mathrm F}$; near $n_{\mathrm e}\approx7\times10^{17}$ cm⁻³ the screened plasmon is resonant with the LO
phonon.

**Include position.** End of `b8-dynamic` ("Screened at the phonon frequency, not at zero"). Sentences just before:
"Near $n_{\mathrm e}\approx7\times10^{17}$ cm⁻³ the factor $\varepsilon_\infty-\omega_{\mathrm p}^2/\omega_{\mathrm{LO}}^2$
crosses zero: the screened plasmon of the carriers is resonant with the LO phonon, the coupled modes of Section B.9.
Figure B8.4 compares the two screenings."

**Physics.** Screening factor of the LO vertex in an isotropic crystal from `b8-eq-relation` (B8.17):
$|g^{\mathrm{dyn}}/g^{\mathrm{dyn},0}|=q^2/|\bq\cdot\boldsymbol\varepsilon^{\mathrm e}\cdot\bq|=1/|\varepsilon^{\mathrm e}(q,\omega)|$,
with $\varepsilon^{\mathrm e}(q,\omega)=\varepsilon_\infty-v(q)\chi^0(q,\omega)$: valence screening $\varepsilon_\infty=10.9$ plus
the $T=0$ Lindhard response of conduction electrons of effective mass $m^*=0.067\,m_{\mathrm e}$ (GaAs, literature values,
illustrative). Write it as $\varepsilon^{\mathrm e}=\varepsilon_\infty+(k_{\mathrm{TF}}^2/q^2)\,L(z,u)$ with $z=q/2k_{\mathrm F}$,
$u=(\hbar\omega+\ii\eta)/(\hbar v_{\mathrm F}q)$ and $L$ the Lindhard function normalized to $L\to1$ for $q\to0$ at
$\omega=0$ (do not call it F). $k_{\mathrm F}=(3\pi^2n_{\mathrm e})^{1/3}$, $E_{\mathrm F}=\hbar^2k_{\mathrm F}^2/2m^*$,
$\hbar v_{\mathrm F}=\hbar^2k_{\mathrm F}/m^*$, $k_{\mathrm{TF}}^2=6\pi e^2n_{\mathrm e}/E_{\mathrm F}$,
$\omega_{\mathrm p}^2=4\pi n_{\mathrm e}e^2/m^*$. Broadening $\eta=\hbar\omega/70$ so that $\omega\to0$ gives the static curve.
Static curve: $L(z,0)=\tfrac12+\frac{1-z^2}{4z}\ln\big|\frac{1+z}{1-z}\big|$. The v1 code (complex Lindhard with a
large-|u| series to avoid cancellation) is correct; keep it, with the variable renamed. Check values:
$n_{\mathrm e}=10^{17}$ cm⁻³: $k_{\mathrm F}=0.0144$ Å⁻¹, $E_{\mathrm F}=11.7$ meV, $\omega_{\mathrm{LO}}/v_{\mathrm F}=0.0222$ Å⁻¹,
$k_{\mathrm{TF}}/\sqrt{\varepsilon_\infty}=0.0146$ Å⁻¹, $2k_{\mathrm F}=0.0287$ Å⁻¹, $\hbar\omega_{\mathrm p}/\sqrt{\varepsilon_\infty}=13.7$ meV,
$\varepsilon^{\mathrm e}(q\to0,\omega_{\mathrm{LO}})=9.33$. Zero crossing of $\varepsilon_\infty-\omega_{\mathrm p}^2/\omega_{\mathrm{LO}}^2$ at
$n_{\mathrm e}=6.9\times10^{17}$ cm⁻³. At $10^{18}$ cm⁻³: $\hbar\omega_{\mathrm p}=143.5$ meV, $\varepsilon^{\mathrm e}(q\to0,\omega_{\mathrm{LO}})=-4.80$.

**What to draw.** One `EPP.Plot`, log–log: x = $q$ (Å⁻¹), fixed $10^{-4}$–1; y = screening factor
$|g^{\mathrm{dyn}}/g^{\mathrm{dyn},0}|$, fixed $10^{-6}$–10. Curves: dynamic $1/|\varepsilon^{\mathrm e}(q,\omega)|$
(`electron`, solid); adiabatic $1/\varepsilon^{\mathrm e}(q,0)$ (`electron`, dashed); undoped reference $1/\varepsilon_\infty$
(`muted`, dotted). Vertical markers: $\omega/v_{\mathrm F}$ (`phonon`), $k_{\mathrm{TF}}/\sqrt{\varepsilon_\infty}$ (`muted`, label
"k_{TF}/√ε_∞"), $2k_{\mathrm F}$ (`muted`), staggered labels so they never overlap.

**Controls.** `range` (log) "carrier density $n_{\mathrm e}$": $10^{15}$–$10^{19}$ cm⁻³, default $10^{17}$. `range` "vertex
frequency ħω": 0–80 meV, step 0.1, default 36.2 (= $\hbar\omega_{\mathrm{LO}}$ of GaAs). Optionally a button "ħω = ħω_{LO}".

**Readouts.** $\hbar\omega_{\mathrm p}/\sqrt{\varepsilon_\infty}$ (meV); $E_{\mathrm F}$ (meV);
$\varepsilon^{\mathrm e}(q\to0,\omega)=\varepsilon_\infty-\omega_{\mathrm p}^2/\omega^2$ ("∞ (static)" at ω = 0);
$k_{\mathrm{TF}}/\sqrt{\varepsilon_\infty}$ (Å⁻¹); $\omega/v_{\mathrm F}$ (Å⁻¹). Use MathJax in the readout keys.

**Initial state.** $n_{\mathrm e}=10^{17}$ cm⁻³, $\hbar\omega=36.2$ meV: the two curves merge above about 0.03 Å⁻¹ and
separate below, the dynamic one flattening at $\approx1/9.3$, the static one falling as $q^2$.

**Caption key points.**
- What to notice: for $q\gg\omega/v_{\mathrm F}$ the carriers follow the field and the curves merge; for smaller $q$
  static screening suppresses the vertex as $q^2$ (its $1/q$ becomes $q$), while at the phonon frequency it stays
  screened only by $\varepsilon_\infty-\omega_{\mathrm p}^2/\omega^2$.
- Marks: solid dynamic, dashed static, dotted undoped; the three vertical markers.
- Try: raise $n_{\mathrm e}$ towards $7\times10^{17}$ cm⁻³ and watch the dynamic curve shoot up at small $q$ as
  $\varepsilon^{\mathrm e}(q\to0,\omega_{\mathrm{LO}})$ crosses zero (plasmon–LO resonance, Section B.9); lower ħω to 0 and see
  the dynamic curve collapse onto the static one.
- Model: `b8-eq-static`, `b8-eq-relation` with the Lindhard model above; $T=0$, parabolic band, $\eta=\hbar\omega/70$;
  GaAs parameters, literature values, illustrative.

**What changes from v1.** v1 labels $|g^{\rm dyn}/g^{\rm dyn,0}|$ (fine, keep), "ε^e(q, 0)", "ε^e(q, ω)" (keep as
$\varepsilon^{\mathrm e}$ with upright e), "1/ε_∞" (keep), but its caption writes the Lindhard function as $F$ (now the phonon
propagator: rename), $\epsilon^\infty$, $\hat\epsilon^{\rm e}$, $k_{\rm TF}$, $\omega_{\rm pl}$ (use $\omega_{\mathrm p}$, the v2 symbol
for the electronic plasma frequency), $m_e$, $n$ (use $n_{\mathrm e}$; $n$ is a band index here), SI formulas with
$\varepsilon_0$ (put SI in the Model part only, or drop them), "q_{TF}" (label it $k_{\mathrm{TF}}/\sqrt{\varepsilon_\infty}$),
readout "ħω_{pl}/√ε_∞" (use ω_p). Rename id and classes (`b8-dynamic-…`). Rewrite the caption in the four-part structure.

---

## 5. `fig-lightline.html` · id `viz-b8-lightline` · PORT of `viz-b-frohlich-lightline` (v1 `doc/epp-primer/sections/B08-frohlich.html`)

**Title.** "LO and TO coupling across the light line" `<span class="tag">interactive</span>`.

**What the reader must learn.** The LO kernel is flat at every $q$ (the Coulomb line knows no $c$), so the Fröhlich
$1/q$ runs straight through the light cone. The TO kernel equals the LO one deep inside the cone, peaks at the light
line $q_{\mathrm c}=\sqrt{\varepsilon_\infty}\,\omega/c$ and falls as $(q_{\mathrm c}/q)^2$ outside; making $c$ larger pushes
$q_{\mathrm c}$ to zero and decouples TO phonons. The cone is 4–5 orders of magnitude below the zone boundary.

**Include position.** Section `b8-lightline`, after the bullet list that reads (B8.19), before the caveat "When the LO–TO
split is clean". Sentence just before: "Figure B8.5 plots both kernels, bare and screened, and the LO vertex itself
across the light line."

**Physics.** `b8-eq-ratio` (B8.19): $\mathcal K^0_{\mathrm{LO}}=1$, $\mathcal K^0_{\mathrm{TO}}=\omega^2/(\omega^2-c^2q^2)$;
with a constant electronic dielectric constant, $\mathcal K_{\mathrm{LO}}=1/\varepsilon$,
$\mathcal K_{\mathrm{TO}}=\varepsilon^{-1}[1-q^2/q_{\mathrm c}^2]^{-1}$, $q_{\mathrm c}=\sqrt\varepsilon\,\omega/c$. Regularize the light-line
pole with a small loss: $\varepsilon=\varepsilon_\infty+\ii\varepsilon_2$ ($\varepsilon_2$ is the absorptive part of
NOTATION §12.3), so $|\mathcal K_{\mathrm{TO}}|=1/|\varepsilon_\infty+\ii\varepsilon_2-(q/q_{\mathrm{vac}})^2|$ with
$q_{\mathrm{vac}}=\omega/c$, $|\mathcal K_{\mathrm{LO}}|=1/|\varepsilon_\infty+\ii\varepsilon_2|$; for the bare TO kernel use
$1/|1+\ii\varepsilon_2/\varepsilon_\infty-(q/q_{\mathrm{vac}})^2|$ as a numerical broadening. Evaluate the TO kernel at
$\omega=\omega_{\mathrm{TO}}$. Panel (b), the LO vertex to leading order in $q$ (eigenvectors at $\bq\to0$, phases 1,
intraband plane-wave matrix element 1, $\omega_{\mathrm{LO}}$ in the zero-point amplitude):
screened $|g^{\mathrm{dyn}}|=\frac{4\pi e^2}{\OV}\frac{Z^*}{\varepsilon_\infty q}\sqrt{\hbar/2M_{\mathrm r}\omega_{\mathrm{LO}}}$ (`b8-eq-static`),
bare $=\varepsilon_\infty$ times that (`b8-eq-bare`), and Fröhlich's
$\frac1q\big[\frac{4\pi e^2}{\OV}\frac{\hbar\omega_{\mathrm{LO}}}{2}(\frac1{\varepsilon_\infty}-\frac1{\varepsilon_{\mathrm s}})\big]^{1/2}$
(`a4-eq-frohlich`). v2 convention: no $\sqrt{N_\ell}$ anywhere (v1's "$\sqrt N|g|$" is simply $|g|$ here). Check values
(GaAs, $q=0.1$ Å⁻¹): screened 304.0 meV, bare 3314 meV, Fröhlich 321.0 meV; $q_{\mathrm c}=5.571\times10^{-5}$ Å⁻¹,
$\omega_{\mathrm{TO}}/c=1.688\times10^{-5}$ Å⁻¹. MgO: $q_{\mathrm c}=4.28\times10^{-5}$ Å⁻¹, screened 3193 meV at 0.1 Å⁻¹;
NaCl: $1.57\times10^{-5}$ Å⁻¹, 1308 meV.

**What to draw.** Two `EPP.Plot` panels in a `viz-grid two` (stacked on phones), both log–log.
- (a) kernels: x = $q$ (Å⁻¹), fixed $10^{-8}$–1 (or $q/q_{\mathrm c}$, $10^{-3}$–$10^3$, by the axis control); y =
  $|\mathcal K|$ (dimensionless, "normalized to the Coulomb value"), fixed $10^{-6}$–$10^4$. Curves: LO bare (`c5` dashed,
  flat at 1), LO screened (`c5` solid, flat at $1/\varepsilon_\infty$), TO bare (`photon` dashed), TO screened (`photon`
  solid). Vertical lines at $\omega/c$ (`muted`, label "ω/c") and $q_{\mathrm c}$ (`photon`, label "q_c" with upright c).
  A faint band from $10^{-3}$ to 1 Å⁻¹ labelled "wavevector grids" (only in the absolute-q mode).
- (b) LO vertex: x = $q$ (Å⁻¹), fixed $10^{-6}$–1; y = $|g|$ (eV), fixed $10^{-2}$–$10^6$. Curves: bare
  $|g^{\mathrm{dyn},0}|$ (`c5` dashed), screened $|g^{\mathrm{dyn}}|$ (`c5` solid, label "Verdi–Giustino"), Fröhlich (`muted`
  dotted). A faint `photon` band from $10^{-6}$ to $q_{\mathrm c}$ labelled "light cone" when it is wide enough.

**Controls.** `select` "material": GaAs (default), MgO, NaCl. `range` (log) "speed of light ×": 1–$10^6$, default 1
(label it in words, e.g. "speed of light (× real value)"; do not use λ, the photon polarization index). `range` (log)
"loss $\varepsilon_2/\varepsilon_\infty$": $10^{-3}$–0.1, default 0.01 (do not call it γ, a damping rate in v2). `seg` "panel
(a) axis": [q (Å⁻¹)] · [q/q_c].

**Readouts.** $q_{\mathrm c}=\sqrt{\varepsilon_\infty}\,\omega_{\mathrm{TO}}/c$ (Å⁻¹); $|\mathcal K_{\mathrm{TO}}/\mathcal K_{\mathrm{LO}}|$
at $q=10^{-3}$ Å⁻¹; $|g^{\mathrm{dyn}}|$ at $q=0.1$ Å⁻¹ (meV); Fröhlich at 0.1 Å⁻¹ (meV); $|g^{\mathrm{dyn},0}/g^{\mathrm{dyn}}|=\varepsilon_\infty$.

**Initial state.** GaAs, real speed of light, loss 0.01, absolute-$q$ axis: the TO peak sits at $5.6\times10^{-5}$ Å⁻¹,
far left of the grid band, and the LO lines are flat.

**Caption key points.**
- What to notice: the LO kernel is flat, so the Coulomb channel never feels the finite speed of light and the Fröhlich
  $1/q$ in (b) runs straight through the light cone; the TO kernel matches LO inside the cone, peaks at $q_{\mathrm c}$
  and falls as $(q_{\mathrm c}/q)^2$.
- Marks: dashed = free line $D^0$ (bare vertex), solid = dressed line $D^{\mathrm e}$ (screened vertex); `c5` LO, `photon`
  TO; the two vertical lines; the grid band; in (b) bare, screened and Fröhlich.
- Try: raise the speed of light by $10^4$ and watch $q_{\mathrm c}$ move left and the TO kernel vanish at every grid
  wavevector; switch the axis to $q/q_{\mathrm c}$ to see the universal resonance shape; lower the loss to sharpen the peak.
- Model: `b8-eq-ratio` with a constant $\varepsilon_\infty$ plus a small loss; (b) from `b8-eq-static` and `b8-eq-bare` at
  leading order in $q$, $\omega_{\mathrm{LO}}$ in the zero-point amplitude; Fröhlich from `a4-eq-frohlich`; the TO vertex has
  no universal form, since it needs transverse current matrix elements. Parameters from the table (literature values,
  illustrative).

**What changes from v1.** v1 uses source notation throughout: "λ_c" for the speed-of-light factor (λ is the photon
polarization index in v2), "γ" for the loss tangent, "ε^e = ε^∞(1+iγ)", "$\sqrt N|g|$" (v2's $g$ has no $N$ inside: plot
$|g|$), "G Eq. 7", "E p.95′/p.98′", "K" for the kernels (use $\mathcal K$; $K$ is a spring stiffness in v2), "q_c" with
italic c, $\epsilon^\infty$, $\Omega_{\rm uc}$, $\mu$ for the reduced mass, "e²/ε₀" in comments. The TO kernel formulas,
the material data and the numbers are right; keep them. Rename id and classes (`b8-lightline-…`). Rewrite the caption in
the four-part structure.
