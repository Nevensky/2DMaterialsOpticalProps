# Section X.2 · Figure specifications

Section file: `v2/sections/X02-dictionary.html` (id `x-dictionary`, prefix `x2`).
Three figures, numbered in reading order: X2.1 `continuity` (new), X2.2 `units` (port), X2.3 `glyphs` (port).

General rules for all three (from `v2/build/BUILDER_GUIDE.md` §5):
- v2 notation only. No source symbols or source equation numbers anywhere in a figure (labels, readouts,
  legends, aria-labels, inspector texts). Refer to equations only with `<a class="eqref" href="#…">` and to
  sections with `<a href="#…">Section …</a>`. Use only equation ids that exist in X.2 or are listed as exports
  in `v2/build/SECTION_BRIEFS.md`.
- `σ` is reserved for the conductivity and `μ` for the chemical potential: neither may appear in figure strings
  (the linter flags `σ`). The reduced mass is $M_{\mathrm r}$, the static dielectric constant $\varepsilon_{\mathrm s}$.
- Colours read at draw time with `EPP.color(...)`; animations through `EPP.anim`; 390 px phone layout with
  controls below the plot; keyboard access and `aria-label`s; no NaN; no console errors.
- Caption: `cap-key` (What to notice) plus `cap-more` (Marks, Try, Model), about 200–260 words.
- Test with `bash v2/build/figtest.sh v2/sections/X02/fig-<name>.html` and inspect both screenshots.

Equation ids of X.2 (tags in brackets) that the figures may cite:
`x2-eq-cont` (X2.1, $\omega\rho=\bq\cdot\bj$), `x2-eq-D0L` (X2.2, $D^{0\mathrm L}_{\alpha\beta}=-4\pi c\,\hat q_\alpha\hat q_\beta/(\omega+\ii\eta)^2$),
`x2-eq-coul` (X2.3, $-\frac1c\bj_2^*\cdot D^{0\mathrm L}\cdot\bj_1=\rho_2^*(4\pi/q^2)\rho_1$),
`x2-eq-T` (X2.4, $-\frac1c\bj_2^*\cdot D^{0\mathrm T}\cdot\bj_1=4\pi\,\bj_2^{\mathrm T*}\cdot\bj_1^{\mathrm T}/[(\omega+\ii\eta)^2-c^2q^2]$),
`x2-eq-De` (X2.5), `x2-eq-Pirho` (X2.6), `x2-eq-eps` (X2.7), `x2-eq-vertex` (X2.8), `x2-eq-subst` (X2.9),
`x2-eq-frohlich` (X2.10), `x2-eq-gB` (X2.11), `x2-eq-residue` (X2.12).

---

## 1. `fig-continuity.html` · id `viz-x2-continuity` · new

**Title.** "Only the current along q carries charge" (`<span class="tag">interactive</span>`; write **q** bold in the
title with MathJax, `$\bq$`).

**What the reader must learn.** By continuity, $\omega\rho(\bq,\omega)=\bq\cdot\bj(\bq,\omega)$ (X2.1): a current wave
builds a charge wave only through its component along the wavevector. The longitudinal part $\bj^{\mathrm L}$
carries the charge that the Coulomb interaction sees; the transverse part $\bj^{\mathrm T}$ flows along the
wavefronts, carries no charge, and couples only to transverse photons, which are weaker than Coulomb by
$(\omega/cq)^2$ outside the light cone. This is the geometric content of all three substitutions (X2.9).

**Include position.** In `#x2-charge`, right after the paragraph that ends
"…A current perpendicular to the wavevector flows along the wavefronts. Along its own direction it has the same
size everywhere, so inflow and outflow balance and no charge appears. [Figure X2.1] shows both cases and every
direction in between." The derivation of (X2.1) follows the figure.

**Physics (v2 notation).**
- One Fourier component (schematic wavelength on screen): $\bj(\br,t)=\mathrm{Re}\big[j_0(\cos\theta\,\hat\bq+\sin\theta\,\mathbf t)\,\ee^{\ii(\bq\cdot\br-\omega t)}\big]$,
  with $\mathbf t$ the in-plane unit vector perpendicular to $\hat\bq$ and $\theta$ the angle between current and $\bq$.
- Charge from continuity (X2.1): $\rho(\br,t)=\mathrm{Re}\big[(q j_0\cos\theta/\omega)\,\ee^{\ii(\bq\cdot\br-\omega t)}\big]$,
  in phase with the longitudinal current; the region where the current converges ($\partial_t\rho>0$) lies a
  quarter wavelength ahead of the charge maximum, so the charge wave travels with the current wave.
- Interaction energies for two equal current waves (amplitudes only; the common factor $4\pi|j_0|^2$ drops out of the ratio):
  Coulomb, from (X2.3): $U_{\mathrm L}\propto\cos^2\theta/\omega^2$; transverse, from (X2.4): $U_{\mathrm T}\propto\sin^2\theta/|\omega^2-c^2q^2|$.
  Ratio $R=U_{\mathrm T}/U_{\mathrm L}=\tan^2\theta\;\omega^2/|\omega^2-c^2q^2|$, which is $\approx\tan^2\theta\,(\omega/cq)^2$ for $cq\gg\omega$
  and diverges on the light line $q=\omega/c$ (display "on the light line" when $|1-c^2q^2/\omega^2|<10^{-3}$;
  display "no Coulomb part" at $\theta=90^\circ$ and "0" at $\theta=0$).
- Constants: $\hbar c=1973.27$ eV Å, so $\omega/c$ in Å⁻¹ is $\hbar\omega[\text{eV}]/1973.27$. Default
  $\hbar\omega=36.2$ meV (the LO phonon of GaAs, literature value, illustrative) gives $\omega/c=1.83\times10^{-5}$ Å⁻¹;
  at $q=0.1$ Å⁻¹, $(\omega/cq)^2=3.4\times10^{-8}$.

**What to draw.**
- Main canvas (aspect about 2:1 on desktop, about 1.2:1 on phone): a rectangle showing 2.5 wavelengths along $x$.
  The wavevector $\bq$ points along $+x$; draw it once as an arrow labelled **q** (`drawRich` bold), and a short
  dashed line marking one wavelength $2\pi/q$.
- Background: vertical stripes shaded by $\rho(x,t)$: positive charge in a tint of `'photon'` (rose), negative in a
  tint of `'electron'` (blue), opacity proportional to $|\rho|/\rho_{\max}$ with $\rho_{\max}$ the value at $\theta=0$
  (so the stripes fade as $\theta\to90^\circ$). Small "+" and "−" glyphs at the extrema help colour-blind readers.
- Current arrows on a grid (about 9 columns × 5 rows, fewer on phone): direction $\cos\theta\,\hat\bq+\sin\theta\,\mathbf t$,
  signed length proportional to $\cos(qx-\omega t)$ (arrows shrink, vanish and reverse), colour `'ink'`.
- Inset (top-right of the canvas on desktop, below the canvas on phone, or a second small canvas): the
  decomposition of the current amplitude. Arrow $\bj$ (`'ink'`), its projection $\bj^{\mathrm L}$ along $\hat\bq$
  (`'photon'`), its transverse part $\bj^{\mathrm T}$ (`'c4'` teal), dashed projection lines, the angle $\theta$ as an
  arc, and the unit vector $\hat\bq$. Labels with `drawRich` bold vectors: `j`, `j^{L}`, `j^{T}`, `q̂`.
- Animate the phase $\omega t$ with `EPP.anim` (calm: one wavelength every 4–5 s on screen; reduced motion honoured;
  Play/Pause label from `EPP.anim`). The drawn wavelength is fixed; $q$ and $\omega$ enter only the readouts.

**Controls (4).**
1. Current direction $\theta$: 0–180°, step 1°, default 30° (`aria-valuetext` "θ = 30°").
2. Wavevector $q$: logarithmic slider $10^{-6}$–1 Å⁻¹, default 0.1 Å⁻¹ (readouts only).
3. Frequency $\hbar\omega$: 1–200 meV, default 36.2 meV (readouts only).
4. Play/Pause of the animation.
Optional: two small buttons "along q (θ = 0°)" and "across q (θ = 90°)".

**Readouts (MathJax in HTML).** $|\bj^{\mathrm L}|/|\bj|=|\cos\theta|$; charge amplitude $|\rho|$ in units of
$q|\bj|/\omega$ (equals $|\cos\theta|$); light-line wavevector $\omega/c$ (Å⁻¹, scientific notation);
$(\omega/cq)^2$; ratio of transverse to Coulomb energy $R$ (scientific notation, with the special cases above).

**Initial state.** $\theta=30^\circ$, $q=0.1$ Å⁻¹, $\hbar\omega=36.2$ meV, animation running: the stripes are clearly
visible (87% of maximum) and the inset shows a sizeable transverse part; $R\approx1\times10^{-8}$.

**Caption.**
- *What to notice.* The charge stripes are set by the component of the current along $\bq$ alone; at $\theta=90^\circ$
  the current flows along the wavefronts and the stripes vanish. The transverse part, invisible to the Coulomb
  interaction, couples to photons with a strength suppressed by $(\omega/cq)^2$ unless $q$ is near the light line.
- *Marks.* Rose and blue stripes: positive and negative charge density; black arrows: current density; inset:
  $\bj$ split into $\bj^{\mathrm L}$ (rose) along $\hat\bq$ and $\bj^{\mathrm T}$ (teal) across it; arrow **q**: wavevector.
- *Try.* (1) $\theta=0^\circ$ and $\theta=90^\circ$: full stripes, then none. (2) Pause: the charge peaks where the current
  along $\bq$ is largest and grows a quarter wavelength ahead, where the arrows converge. (3) At $\theta=45^\circ$ lower
  $q$ towards $\omega/c\approx2\times10^{-5}$ Å⁻¹: $R$ grows to order 1 and diverges on the light line.
- *Model.* Plane wave with continuity (X2.1) `x2-eq-cont`; energies from (X2.3) `x2-eq-coul` and (X2.4) `x2-eq-T`;
  $\hbar c=1973.27$ eV Å; default frequency the GaAs LO phonon, 36.2 meV. A physical example is the polarization
  current $\partial_t\bP$ of an optical phonon or the electrons of a plasma wave. Lengths on screen are schematic.

---

## 2. `fig-units.html` · id `viz-x2-units` · port of `viz-x-dictionary-units`

**v1 source.** `doc/epp-primer/sections/X02-dictionary.html`: the `<figure id="viz-x-dictionary-units">` element
(lines 207–213) and its code, the block "2. one Fröhlich vertex, two notations" inside the section-level script
(lines 324–403). Move the code into the figure file's own `<script>` (IIFE, `document.getElementById('viz-x2-units')`).

**Title.** "One Fröhlich vertex in two unit systems and two bookkeepings" (`interactive`).

**What the reader must learn.** (i) The Gaussian and SI forms (X2.10) of the long-range vertex $g^{\mathrm{lr}}$ give the
same number. (ii) The Part B vertex $g^{\mathrm{dyn}}$, built on the bare phonon at $\omega_{\mathrm{TO}}$, is larger by
$(\omega_{\mathrm{LO}}/\omega_{\mathrm{TO}})^{1/2}$ (X2.11), but multiplied by the square root of the residue weight
$\omega_{\mathrm{TO}}/\omega_{\mathrm{LO}}$ (X2.12) it lies exactly on Part A's line. (iii) The unscreened vertex
$g^{\mathrm{dyn},0}$ is a factor $\varepsilon_\infty$ above. (iv) Fröhlich's form coincides with $g^{\mathrm{lr}}$ when $Z^*$ obeys the LO–TO relation.

**Include position.** In `#x2-bookkeeping`, after the paragraph "For GaAs, with $\hbar\omega_{\mathrm{TO}}=33.3$ meV, the ratio is
$(36.2/33.3)^{1/2}=1.043$, a vertex 4% larger. For MgO, with $\hbar\omega_{\mathrm{TO}}=49.1$ meV and $\hbar\omega_{\mathrm{LO}}=89$ meV,
it is 1.35. [Figure X2.2] plots the vertices of both bookkeepings, in both unit systems, for these crystals."
The residue equation (X2.12) and its derivation follow the figure.

**Physics (v2 notation; keep v1's algebra, which is correct).** Cubic crystal, two atoms per primitive cell (fcc),
longitudinal optical mode along $\bq$, $\bG=0$, unit Bloch overlap, static screening $\varepsilon^{\mathrm e}\to\varepsilon_\infty$,
leading order in $q$ (eigenvectors at $\bq\to0$, basis phase 1).
- $M_{\mathrm r}=M_1M_2/(M_1+M_2)$, $\OV=a^3/4$.
- Constants: $e^2=14.3996$ eV Å (Gaussian) $=e^2/4\pi\varepsilon_0$; $e^2/\varepsilon_0=180.95$ eV Å; $\hbar^2/(\text{amu Å}^2)=4.1802$ meV
  (NOTATION §7; v1 used 4.1803); $\hbar c=1973.27$ eV Å.
- Born charge from the LO–TO relation `a4-eq-lot`: $Z^{*2}=(\omega_{\mathrm{LO}}^2-\omega_{\mathrm{TO}}^2)\,\varepsilon_\infty\OV M_{\mathrm r}/4\pi e^2$
  (code: `Z = sqrt((ELO²−ETO²)·eps·Ω·Mr/(4π·E2·H2))`, energies in eV); presets give GaAs 2.18, NaCl 1.13, MgO 1.97.
- Ionic plasma frequency, both systems: $(\hbar\omega^{\mathrm{ion}}_{\mathrm p})^2=4\pi e^2Z^{*2}\hbar^2/(\OV M_{\mathrm r})$ (Gaussian)
  $=e^2Z^{*2}\hbar^2/(\varepsilon_0\OV M_{\mathrm r})$ (SI): the same number, computed both ways.
- $\hbar\omega_{\mathrm{LO}}=[(\hbar\omega_{\mathrm{TO}})^2+(\hbar\omega^{\mathrm{ion}}_{\mathrm p})^2/\varepsilon_\infty]^{1/2}$ (`a4-eq-lot`);
  $\varepsilon_{\mathrm s}=\varepsilon_\infty\omega_{\mathrm{LO}}^2/\omega_{\mathrm{TO}}^2$ (Lyddane–Sachs–Teller, `b3-eq-lst`).
- Part A, (X2.10) `x2-eq-frohlich`: $|g^{\mathrm{lr}}(q)|=\frac{4\pi e^2}{\OV}\frac{Z^*}{\varepsilon_\infty q}\big(\frac{\hbar}{2M_{\mathrm r}\omega_{\mathrm{LO}}}\big)^{1/2}$.
- Part B, (X2.11) `x2-eq-gB`: $|g^{\mathrm{dyn}}|=(\omega_{\mathrm{LO}}/\omega_{\mathrm{TO}})^{1/2}|g^{\mathrm{lr}}|$ (bare phonon at $\omega_{\mathrm{TO}}$);
  unscreened $|g^{\mathrm{dyn},0}|=\varepsilon_\infty|g^{\mathrm{dyn}}|$.
- Part B with residue, (X2.12) `x2-eq-residue`: $|g^{\mathrm{dyn}}|\,(\omega_{\mathrm{TO}}/\omega_{\mathrm{LO}})^{1/2}=|g^{\mathrm{lr}}|$.
- Fröhlich, `a4-eq-frohlich`: $|g_{\mathrm F}(q)|=\frac1q\big[\frac{4\pi e^2}{\OV}\,\frac{\hbar\omega_{\mathrm{LO}}}{2}\big(\frac1{\varepsilon_\infty}-\frac1{\varepsilon_{\mathrm s}}\big)\big]^{1/2}$.
- Upper polariton at $\bq\to0$ without electrons, `b7-eq-410`: $[(\hbar\omega_{\mathrm{TO}})^2+(\hbar\omega^{\mathrm{ion}}_{\mathrm p})^2]^{1/2}$.
- Polariton crossover $\sqrt{\varepsilon_\infty}\,\omega_{\mathrm{TO}}/c$ (Å⁻¹), below which TO phonons couple to light (`b8-eq-ratio`).
- Presets (literature values, illustrative; keep v1's): GaAs $M=69.723,\,74.922$ amu, $a=5.653$ Å, $\varepsilon_\infty=10.9$,
  $\hbar\omega_{\mathrm{TO}}=33.3$, $\hbar\omega_{\mathrm{LO}}=36.2$ meV; NaCl $22.990,\,35.453$ amu, $5.640$ Å, $2.34$, $20.3$, $32.7$ meV;
  MgO $24.305,\,15.999$ amu, $4.212$ Å, $2.96$, $49.1$, $89.0$ meV.
- Checks: GaAs with $Z^*=2.07$ (not the preset) gives $|g^{\mathrm{lr}}|=304$ meV at $q=0.1$ Å⁻¹ (quoted in the text);
  GaAs preset ($Z^*=2.18$): 320 meV; MgO preset at $q=0.1$ Å⁻¹: $|g^{\mathrm{lr}}|\approx3.18$ eV.

**What to draw.** `EPP.Plot`, log–log, $|g|$ (eV) against $q$ (Å⁻¹), fixed limits $x\in[10^{-3},1]$, $y\in[10^{-2},10^4]$, aspect about 0.56.
- Part A $g^{\mathrm{lr}}$ with $\omega_{\mathrm{LO}}$: `'ink'`, solid, width 2.5.
- Part B $g^{\mathrm{dyn}}$ with bare $\omega_{\mathrm{TO}}$: `'photon'`, dashed.
- Part B times the residue factor: `'c4'`, dotted, width 2.5 (lies on the ink line).
- Part B unscreened $g^{\mathrm{dyn},0}$: `'c5'`, dash-dot.
- Fröhlich: `'muted'` hollow circles at six $q$ values.
- Readout $q$: vertical `'muted'` line.
- Legend labels (`drawRich`): "Part A: g^{lr}, ω_{LO}", "Part B: g^{dyn}, bare ω_{TO}", "Part B × residue",
  "Part B unscreened: g^{dyn,0}", "○ Fröhlich". Add the circle to the legend properly (v1 injected HTML into `.legend`).

**Controls (5).** Material preset (GaAs, NaCl, MgO; sets $Z^*$, $\varepsilon_\infty$, $\hbar\omega_{\mathrm{TO}}$), default MgO;
Born charge $Z^*$ 0.2–4 (dimensionless); $\varepsilon_\infty$ 1–15; bare $\hbar\omega_{\mathrm{TO}}$ 5–200 meV; readout $q$,
logarithmic 0.001–1 Å⁻¹, default 0.05. A small "reset to preset" button is welcome.

**Readouts.** $4\pi e^2/\OV$ (Gaussian) and $e^2/\varepsilon_0\OV$ (SI), eV/Å², side by side and equal;
$\hbar\omega^{\mathrm{ion}}_{\mathrm p}$ (Gaussian | SI), meV; $\hbar\omega_{\mathrm{LO}}$ from the LO–TO relation; $\varepsilon_{\mathrm s}$;
$|g^{\mathrm{lr}}|(q)$ and $|g^{\mathrm{dyn}}|(q)$ (meV below 1 eV); ratio $(\omega_{\mathrm{LO}}/\omega_{\mathrm{TO}})^{1/2}$; $M_{\mathrm r}$, $\OV$;
upper polariton at $\bq\to0$ without electrons (meV); crossover $\sqrt{\varepsilon_\infty}\,\omega_{\mathrm{TO}}/c$ (Å⁻¹, scientific notation).

**Initial state.** MgO preset (largest effect: ratio 1.35), $q=0.05$ Å⁻¹.

**Caption.**
- *What to notice.* The dotted Part B line with the residue lies exactly on Part A's solid line, and Fröhlich's circles
  sit on both: paired with its own phonon propagator, each bookkeeping gives the same coupling. With the bare
  $\omega_{\mathrm{TO}}$ the Part B vertex is larger by $(\omega_{\mathrm{LO}}/\omega_{\mathrm{TO}})^{1/2}$ (4% for GaAs, 35% for MgO); the unscreened one by a further $\varepsilon_\infty$.
- *Marks.* The five line styles above; vertical line: readout $q$; the paired Gaussian and SI readouts agree to the last digit.
- *Try.* Switch presets; drag $\hbar\omega_{\mathrm{TO}}$ towards $\hbar\omega_{\mathrm{LO}}$ (small splitting: the Part B line approaches Part A);
  set $\varepsilon_\infty=1$ (screened and unscreened coincide); vary $Z^*$ (all lines scale, $\hbar\omega_{\mathrm{LO}}$ follows).
- *Model.* As in Physics above, with (X2.10)–(X2.12), `a4-eq-lot`, `a4-eq-frohlich`; $Z^*$ chosen so that the LO–TO
  relation reproduces the measured $\hbar\omega_{\mathrm{LO}}$ (GaAs 2.18, against a literature 2.07). For a longitudinal current the
  photon line does not depend on $c$ (Section B.8), so these vertices hold at every $q$; finite $c$ matters only for TO
  modes, below the crossover readout. At $q=0.5$ Å⁻¹ the omitted basis phase would change $|g|^2$ by up to about 20% in GaAs.

**What in v1 must change.**
- Source notation everywhere: $Z^*_{\mathrm{EPP}}=|e|Z^*$ readout (drop it, or show $eZ^*$ in $(\text{eV Å})^{1/2}$ with $e=3.7947$),
  $\omega_\sigma$ (→ bare $\omega_{\mathrm{TO}}$), $\varepsilon^\infty$ written `ε^∞` (→ $\varepsilon_\infty$), `ε^0` (→ $\varepsilon_{\mathrm s}$), `μ` for the reduced
  mass (→ $M_{\mathrm r}$), `Ω` (→ $\OV$), `√N|g^{dyn}|` (→ $|g^{\mathrm{dyn}}|$, with $N_\ell^{-1/2}$ outside as for $g$).
- Source equation labels in legend, readouts and caption ("Giustino (G.180)", "EPP p. 95′/98′", "E.454", "E.410′",
  "G.2", "G.7") must go; use "Part A / Part B" and v2 equation links.
- `σ` appears in v1 strings (`ω_σ`) and is forbidden by the linter.
- v1 caption is a single "Model" paragraph; rewrite in the four-part structure. Remove links to `#z-corrections`.
- Classes and ids: `x-dict-*` → `x2-units-*`; the code must live in the figure file.

---

## 3. `fig-glyphs.html` · id `viz-x2-glyphs` · port of `viz-x-dictionary-glyphs`

**v1 source.** `doc/epp-primer/sections/X02-dictionary.html`: the `<figure id="viz-x-dictionary-glyphs">` element
(lines 226–277, with the inspector texts) and its code, the block "3. diagram dictionary" of the section-level
script (lines 405–505, using `EPP.fd.draw`). Move the code into the figure file.

**Title.** "Diagram dictionary" (`interactive`).

**What the reader must learn.** Each diagram element of Part A has a counterpart in Part B, or none (the transverse
photon). A Part A electron–phonon vertex is not one Part B glyph but a composite of three: an electron current
vertex, a photon line and the ionic current vertex. Sending $c\to\infty$ collapses the Part B column onto Part A,
with the Born charges in place of the bare ionic charges and only the dipolar coupling kept.

**Include position.** In `#x2-glyphs`, after the background box "Reading a Feynman diagram", which follows the
paragraph ending "…A dashed rose line is longitudinal (the Coulomb interaction), a wavy one is transverse
(photons), a wavy line over a dashed one carries both polarizations, as $D^0$ does, and a heavy or doubled line is
dressed. Dots are charge vertices and squares are current vertices." The callout "Physical picture: reading a
diagram across" follows the figure.

**Layout.** Two columns (stacked on phone): "Part A · electrostatic" and "Part B · electrodynamic". Each glyph is a
`<button>` with an SVG glyph (`EPP.fd.draw`, 104×40, `aria-hidden`), a symbol (MathJax) and a short name. Below
the columns: an inspector panel (`role="status"`, `aria-live="polite"`) showing the translation of the hovered,
focused or pinned glyph; counterparts light up (rose border and tint), the selected glyph gets the accent border,
the others fade. Keep v1's interaction (hover, focus, tap to pin, Esc to release, scroll the inspector into view
on phone) and visual language: blue arrowed electron lines, amber zigzag phonons, rose electromagnetic lines
(dashed = longitudinal, wavy = transverse, wavy over dashed = $D^0$, heavy wavy over double dashed = $D^{\mathrm e}$),
dots = charge vertices, squares = current vertices (blue square for the electron current vertex, teal `'c4'`
square for the ionic one).

**Glyphs (v2 notation).**
- Part A: `gG` $G$, electron propagator; `gF` $F$, phonon propagator on adiabatic phonons; `gv` $v$, instantaneous
  Coulomb line; `gW` $W^{\mathrm e}$, $W$, screened interaction; `gn` charge vertex (a Coulomb line ends on an electron);
  `gZ` ionic charge vertex (a displaced ion acts through the screened gradient of its potential, bare charge $eZ_\kappa$);
  `gb` bare electron–phonon vertex $\langle\psi_{m\bk+\bq}|\Delta_{\bq\nu}V^{\mathrm{en}}|\psi_{n\bk}\rangle$ (small dot);
  `gg` screened vertex $g$ (grey disk).
- Part B: `eG` electron lines (only inside bubbles); `eF` $F$ on bare phonons; `e0L` $D^{0\mathrm L}$, Coulomb as the
  longitudinal photon; `e0T` $D^{0\mathrm T}$, transverse photon; `e0` $D^0=D^{0\mathrm L}+D^{0\mathrm T}$; `eDe` $D^{\mathrm e}$, photon line
  dressed by electrons; `erho` charge vertex $\rho^{\,n\bk,m\bk+\bq}$ (dot); `ej` electron current vertex
  $j_\alpha^{\,n\bk,m\bk+\bq}$ (blue square); `eJ` ionic current vertex $J_{\alpha\nu}(\bq,\omega)$ (teal square).
  (v1's separate `eJe` glyph, the dynamic electron vertex, merges into `ej`; the inspector text of `ej` may mention
  that Section B.8 uses its dynamic form.)

**Counterparts.** gG↔eG; gF↔eF; gv↔e0L, e0; gW↔eDe; gn↔erho, ej; gZ↔eJ; gb↔ej, e0, eJ (composite $g^{\mathrm{dyn},0}$);
gg↔ej, eDe, eJ (composite $g^{\mathrm{dyn}}$). Reverse: eG→gG; eF→gF; e0L→gv; e0T→none; e0→gv, gb; eDe→gW, gg;
erho→gn; ej→gn, gb, gg; eJ→gZ, gb, gg.

**Inspector texts (key content; v2 notation and v2 links only).**
- gG / eG: $G$ (Section A.5) ↔ electron lines inside the bubbles of $\Pi^{\mathrm e}$ (`b5-eq-Pie`); both use Kohn–Sham
  states in practice; Part B has no Dyson equation for electrons, hence no counterpart of the Fan–Migdal self-energy.
- gF / eF: $F$ on adiabatic phonons (`a5-eq-F`) ↔ $F^0$ on bare phonons (`b5-eq-F0`): same definition and units, Dyson
  equations of the same form (`a6-eq-dysonF`, `b7-eq-dysonF`); the reference phonons differ (adiabatic, with
  $\Phi^{\mathrm{dd}}$; bare, short-range only), and the dressed Part B phonon carries the weight of (X2.12).
- gv: $v(\bq)=4\pi e^2/q^2$ between charge vertices ↔ $D^{0\mathrm L}$ (X2.2) between current vertices, exactly Coulomb by (X2.3);
  $D^0$ adds the transverse part; $c\to\infty$ returns $D^{0\mathrm L}$.
- gW: $W^{\mathrm e}=v+v\chi^{\mathrm{irr}}W^{\mathrm e}$ (`x1-eq-dysonW`), $W=W^{\mathrm e}+W^{\mathrm{ph}}$ (`a5-eq-split`) ↔ $D^{\mathrm e}$
  (`b5-eq-dysonDe`) and $D$ (`b6-eq-dysonD`); $e^2D^{\mathrm e}_{\rho\rho}=W^{\mathrm e}$; for $c\to\infty$, (X2.5).
- gn / erho: the charge vertex is the same object in both parts, $\rho^{\,n\bk,m\bk+\bq}=-e\langle\psi_{n\bk}|\ee^{-\ii\bq\cdot\br}|\psi_{m\bk+\bq}\rangle$;
  (X2.8) rewrites it as the longitudinal part of the current vertex, so `ej` lights up too.
- gZ / eJ: Part A's ion acts through the screened gradient of its potential with the bare charge (`a5-eq-Wph`) ↔
  the ionic current vertex $J_{\alpha\nu}(\bq,\omega)$ (`b7-eq-J`): Born charges, basis phase $\ee^{-\ii(\bq+\bG)\cdot\btau_\kappa}$, a factor
  $\omega$ (a vibrating dipole is a current), long-range dipole only; for $c\to\infty$ only $\bq\cdot\bj_\nu$ couples, the
  numerator of $g^{\mathrm{lr}}$ (`a4-eq-glr`).
- gb: unscreened vertex ↔ composite $g^{\mathrm{dyn},0}$: electron current vertex, $D^0$, $J$ (`b8-eq-vertex`); for
  $c\to\infty$ it is the dipolar part of the bare vertex, $g^{\mathrm{lr}}$ with $\varepsilon_\infty\to1$ (Section B.8).
- gg: statically screened $g$ of DFPT (`a3-eq-scf`) ↔ composite $g^{\mathrm{dyn}}$ with $D^{\mathrm e}$ (`b8-eq-vertex`); for
  $c\to\infty$, $g^{\mathrm{lr}}$ with $\varepsilon_\infty\to\varepsilon^{\mathrm e}(\bq,\omega)$ (`b8-eq-static`), screened at the phonon frequency;
  both phonon self-energies keep one bare and one screened vertex.
- e0L: $D^{0\mathrm L}$ (X2.2) ↔ $v$; the quantized longitudinal photon has a vanishing propagator, while $D^{0\mathrm L}$ is the
  Coulomb energy of the longitudinal field (Section B.4): both hold.
- e0T: $D^{0\mathrm T}$ (`b5-eq-D0T`) ↔ no counterpart: Part A has light only as a classical external field (Section A.9);
  it vanishes as $1/c$ and matters for $q\lesssim\omega/c$.
- e0: vacuum Maxwell propagator in the temporal gauge ↔ $v$, its $c\to\infty$ limit; inside $g^{\mathrm{dyn},0}$ it plays the
  bare Coulomb potential of the displaced ions.
- eDe: $D^{\mathrm e}$ ↔ $W^{\mathrm e}$; at finite $c$ its transverse part describes light in the electronic medium, which $W^{\mathrm e}$ cannot.
- ej: current vertex $j_\alpha^{\,n\bk,m\bk+\bq}$ (`b5-eq-jel`) ↔ charge vertex, for its longitudinal part only (X2.8); the
  transverse part couples to $D^{0\mathrm T}$ and has no Part A counterpart; diamagnetic $A^2$ terms accompany current vertices.

**Electrostatic-limit notes** (hidden until the toggle is on; shown under the glyph): e0L "unchanged: the energy
$-\frac1c\bj\cdot D^{0\mathrm L}\cdot\bj$ does not depend on $c$"; e0T "vanishes as $1/c$" (also fade the button); e0 "→ $D^{0\mathrm L}$";
eDe "→ $-(c\,q_\alpha q_\beta/\omega^2)\,W^{\mathrm e}/e^2$, (X2.5)"; ej "only $\bq\cdot\bj$ survives, a charge vertex by (X2.8)";
eJ "only $\bq\cdot\bj_\nu$ couples, the numerator of $g^{\mathrm{lr}}$".

**Controls (2).** Glyph selection (hover, focus, tap to pin; Esc releases) and the toggle "electrostatic limit $c\to\infty$".

**Initial state.** `gb` (bare vertex) pinned, so that the three Part B glyphs of its composite are lit at load and the
inspector shows its translation; limit toggle off.

**Caption.**
- *What to notice.* A Part A vertex is a composite of three Part B glyphs: the bare vertex is an electron current
  vertex, the free line $D^0$ and the ionic vertex $J$. $D^{0\mathrm T}$ has no counterpart. With $c\to\infty$ the Part B column
  collapses onto Part A, with Born charges where Part A has bare ionic charges.
- *Marks.* Colours and line styles as above; dots: charge vertices; squares: current vertices (blue electron, teal
  ion); accent border: selected glyph; rose border: its counterparts.
- *Try.* Select $g$: $D^{\mathrm e}$ lights up instead of $D^0$. Select $D^{0\mathrm T}$: nothing lights in Part A. Switch on the limit.
- *Model.* The translations of Section X.2, (X2.2), (X2.5), (X2.8), and the definitions `b5-eq-dysonDe`, `b8-eq-vertex`,
  `b7-eq-J`, `a5-eq-split`, `a5-eq-Wph`.

**What in v1 must change.**
- Source symbols: $W_{\mathrm e}$, $W_{\mathrm{ph}}$, $\mathcal D$ (phonon), $g^{\mathrm b}$, $g^{\mathrm{cc}}$, $\Omega^{\mathrm e}$, $\Pi_{jj}$, $D^{0l}$, $D^{0t}$, $D^{\mathrm e}_{11}$,
  EPP's numeric $\rho_{n\bk,m\bk+\bq}$, $J^\mu$ with Cartesian $\mu,\nu$, $Z_\kappa\nabla\delta$, $\br^0_{\mathrm a}$, $m_{\mathrm a}$: convert to the v2 symbols above.
- Source equation numbers in every inspector text (G.70, G.109, E.251, E.374, E.428′, …) must go; cite v2 equations
  with `a.eqref` (ids above) and sections with links.
- Column titles "Giustino 2017 · electrostatic" and "EPP · electrodynamic" → "Part A · electrostatic", "Part B · electrodynamic".
- The remark that EPP colours its figures differently belongs to the section's footnote (already there), not to the figure.
- Links to `#z-corrections` removed; classes `x-dict-*` → `x2-glyphs-*`; v1 caption has only "What to notice":
  add Marks, Try and Model.
