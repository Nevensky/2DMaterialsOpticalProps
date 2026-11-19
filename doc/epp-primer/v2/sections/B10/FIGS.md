# B.10 figures: specifications for the figure authors

Section file: `v2/sections/B10-2d.html` (id `b-2d`, prefix `b10`). Three figures, in reading order:

| File | id | Kind |
|---|---|---|
| `B10/fig-field.html` | `viz-b10-field` | port of v1 `viz-b-2d-field` |
| `B10/fig-disp.html` | `viz-b10-disp` | port of v1 `viz-b-2d-disp` |
| `B10/fig-keldysh.html` | `viz-b10-keldysh` | new |

v1 source of both ports: `doc/epp-primer/sections/B10-2d.html` (figure elements at lines 17–26 and 113–120; one shared
`<script>` at the end, "Figure 1" block lines ~146–245, "Figure 2" block lines ~247–453; helper `sci()` and constants at
the top of the script). Each v2 figure is its own file with its own `<script>`.

Notation reminders (v2, NOTATION.md): in-plane wavevector $Q_\parallel$ (canvas text `Q_∥`), cell area $\Omega_{\mathrm S}$
(`Ω_S`), reduced mass $M_{\mathrm r}$ (`M_r`), Born charge $Z^*$ dimensionless, $\hbar\omega_{\mathrm{TO}}$ (`ħω_{TO}`),
screening length $r_0$, environment $\varepsilon_{\mathrm{env}}$ (`ε_{env}`), sheet polarization $\mathbf P^{\mathrm{2D}}$
(bold, `P^{2D}` with `{bold:true}`), sheet dielectric function $\varepsilon^{\mathrm{2D}}$. Forbidden in figure strings:
`σ` (reserved for the conductivity; the lint flags it), `λ` as a wavelength (λ is the photon polarization index), EPP
symbols ($\Omega_{\rm uc}$, $S_{\rm uc}$, $v_Q$, $\epsilon_{\rm 2D}$, $\Omega^{\rm e}$, $\beta$ for the decay constant),
EPP equation numbers ("E.472"). Equation references go in the caption with `<a class="eqref" href="#b10-eq-…">`.
Lengths in Å, energies in meV or eV.

Numbers used below (from NOTATION §7): $e^2=14.3996$ eV Å, $\hbar c=1973.27$ eV Å, $\hbar^2/2m_{\mathrm e}=3.80998$ eV Å²,
$\hbar^2/(\mathrm{amu}\cdot\text{Å}^2)=4.1802$ meV. Hence $\hbar^2\cdot2\pi e^2/(\mathrm{amu}\cdot\text{Å}^2)=0.378205$ eV² Å
and $(\hbar\omega_{\mathrm p}^{\mathrm{2D}})^2=689.42\;n^{\mathrm{2D}}[\text{Å}^{-2}]\,Q_\parallel[\text{Å}^{-1}]/(m^*/m_{\mathrm e})$ eV²
($10^{13}$ cm⁻² $=10^{-3}$ Å⁻²).

---

## 1. `fig-field.html` · `viz-b10-field` · port of v1 `viz-b-2d-field`

**Title.** "Bulk versus sheet: the depolarizing field of an LO wave" `<span class="tag">interactive</span>`.

**What the reader must learn.** The bound charge of an LO polarization wave forms infinite planes in a bulk crystal and
lines in a sheet. In bulk the field between the planes is a capacitor field, the same at every wavelength. In the sheet
the field lines leave the plane and arc out to a height of about $1/Q_\parallel$, so the in-plane field falls in
proportion to $Q_\parallel$. A TO wave has no bound charge and no field in either geometry.

**Include position.** In `<h3 id="b10-picture">`, after the third paragraph and before "The section follows this picture
in steps:". The sentences just before the include:
> "The longer the wavelength, the higher the arcs and the smaller the share of each line that runs inside the sheet, where
> it can push the ions. The restoring field in the plane therefore weakens in proportion to the wavevector, and at the zone
> centre only the springs remain: the LO and TO frequencies meet at $\Gamma$. Figure B10.1 compares the two geometries."

**Physics** (vacuum, no electronic screening unless the optional toggle is on; side view, $x$ along $\mathbf Q_\parallel$,
$z$ normal to the sheet).
- Sheet, LO: $\mathbf P^{\mathrm{2D}}=P_0\cos(Q_\parallel x)\,\hat x$, the real part of (B10.4) `b10-eq-P2d`. Bound charge
  per area $\rho^{\mathrm{2D}}=-\partial_xP^{\mathrm{2D}}=Q_\parallel P_0\sin(Q_\parallel x)$. Potential (B10.1) `b10-eq-A0`:
  $A_0=2\pi P_0\sin(Q_\parallel x)\,\mathrm e^{-Q_\parallel|z|}$. Field: $E_x=-2\pi Q_\parallel P_0\cos(Q_\parallel x)\,\mathrm e^{-Q_\parallel|z|}$
  (B10.5) `b10-eq-E2d`, $E_z=2\pi Q_\parallel P_0\sin(Q_\parallel x)\,\mathrm{sgn}(z)\,\mathrm e^{-Q_\parallel|z|}$.
  The field lines are exactly the contours of the stream function $\cos(Q_\parallel x)\,\mathrm e^{-Q_\parallel|z|}=c_k$; take
  $c_k=(k-\tfrac12)/N$, $k=1..N$ ($N=5$), which spaces them at equal flux (the v1 code does this correctly). A line of level
  $c$ spans $|x-x_m|\le\arccos(c)/Q_\parallel$ around a charge maximum $x_m$ and peaks at height $\ln(1/c)/Q_\parallel$.
- Bulk, LO: $\mathbf P=P_0\cos(qx)\,\hat x$ everywhere; bound charge planes; $E_x=-4\pi P_0\cos(qx)$, uniform in $z$
  (`a4-eq-E` with $\varepsilon_\infty=1$): horizontal field lines whose density does not change with the wavelength.
- Ratio of in-plane fields for a layer of thickness $d$ cut from the bulk (same dipoles, $\mathbf P=\mathbf P^{\mathrm{2D}}/d$):
  $E^{\mathrm{sheet}}_\parallel/E^{\mathrm{bulk}}=Q_\parallel d/2$, with $d=3.33$ Å (hBN). At a wavelength of 100 Å it is 0.105.
- TO (both panels): dipoles perpendicular to $\mathbf Q_\parallel$ (along $y$, out of the drawing plane), $\nabla\cdot\mathbf P=0$:
  no bound charge, no field.
- Optional screening toggle: bulk field divided by $\varepsilon_\infty=4.87$; sheet field divided by $1+r_0Q_\parallel$ with
  $r_0=10$ Å (B10.7) `b10-eq-eps2d`. Implement by keeping the flux per line fixed: draw only the levels
  $c_k=(k-\tfrac12)(1+r_0Q_\parallel)/N<1$ in the sheet and $\mathrm{round}(N_{\mathrm{rows}}/\varepsilon_\infty)$ rows in bulk
  (use $N_{\mathrm{rows}}=10$ unscreened so that two remain). Skip the toggle if it clutters; then add nothing else.

**What to draw.** Two canvases side by side (`viz-grid two`, stacked on phones), each with a small header in HTML:
"Bulk crystal: charge planes" and "Sheet at $z=0$: field leaks out to $\sim1/Q_\parallel$". Window 100 Å wide, $\pm35$ Å high.
- Dipoles: arrows along $x$ in the phonon colour `EPP.color('phonon')`, length $\propto|\cos|$, at spacing 5 Å; several rows
  in bulk, one row at $z=0$ in the sheet. In TO mode draw ⊙ and ⊗ circles of size $\propto|\cos|$ instead.
- Bound charge: a horizontal colour gradient, $+$ in `'photon'`, $-$ in `'c4'` (alpha $\propto|\sin|$), over the whole bulk
  panel and as a thin strip ($\pm1.6$ Å) in the sheet; "+" and "−" glyphs at the charge extrema.
- Field lines: `--fg-2`, width 1.3, arrowheads from $+$ to $-$ (in the plane the field opposes $\mathbf P$).
- Sheet: the plane $z=0$ as a thin `--fg` line; a dashed `--muted` line at $z=1/Q_\parallel$ labelled `z = 1/Q_∥` (hide it when
  $1/Q_\parallel$ leaves the window).
- In TO mode, a centred note in each panel: "TO: no bound charge, no field".

**Controls (2–3).**
1. Wavelength $2\pi/Q_\parallel$: range slider, log, 20–160 Å, initial 50 Å, format "50 Å". Label "wavelength 2π/Q_∥ (Å)".
2. Polarization: seg control LO / TO, initial LO.
3. Optional toggle "electronic screening" (off initially), as above.

**Readouts.** $Q_\parallel$ (Å⁻¹, 3 significant digits), $1/Q_\parallel$ (Å), "sheet/bulk in-plane field = Q_∥d/2" (with $d=3.33$ Å;
when screening is on, $Q_\parallel d\,\varepsilon_\infty/[2(1+r_0Q_\parallel)]$). In TO mode the ratio reads "no field".

**Initial state.** LO, wavelength 50 Å ($Q_\parallel=0.126$ Å⁻¹, $1/Q_\parallel=8$ Å, ratio 0.21), screening off.

**Caption parts.**
- *What to notice.* The bulk field is the same at every wavelength; in the sheet the field lines arc out to about
  $1/Q_\parallel$ and the in-plane field falls in proportion to $Q_\parallel$, so at long wavelength nothing is left to stiffen the
  LO mode.
- *Marks.* Arrows: Born-charge dipoles (phonon colour). Shading and glyphs: bound charge (+ and −). Grey lines: field lines
  at equal flux spacing, so their density is the field strength. Dashed line: $z=1/Q_\parallel$.
- *Try.* Stretch the wavelength from 20 to 160 Å: the bulk panel does not change, the sheet's arcs grow, and the ratio falls
  from 0.52 to 0.065. Switch to TO: no charge, no field in either panel. (With the toggle: screening weakens both, the sheet
  less at long wavelength.)
- *Model.* (B10.4), (B10.1) and (B10.5) for the sheet, `a4-eq-E` with $\varepsilon_\infty=1$ for the bulk; point dipoles in
  vacuum; $d=3.33$ Å as in hBN; the field lines are exact contours of $\cos(Q_\parallel x)\,\mathrm e^{-Q_\parallel|z|}$.

**What changes from v1.**
- id `viz-b-2d-field` → `viz-b10-field`; classes `b-2d-f3`, `b-2d-f2` → `b10-field-bulk`, `b10-field-sheet`.
- v1 works in nm (window 10 nm, slider 2–16 nm); v2 uses Å (100 Å, 20–160 Å).
- The readout key `λ` and the slider label "wavelength λ = 2π/Q" use λ, the photon polarization index in v2: write
  "wavelength 2π/Q_∥". `Q` → `Q_∥`, `1/Q` → `1/Q_∥`; the readout key `E_∥^{2D}/E^{3D} = Qd/2 (d = 3.33 Å)` → "sheet/bulk
  in-plane field = Q_∥d/2".
- v1's caption uses $P$ and $P_{\rm 2D}$ (scalar, source style) and has the order Model / What to notice: rewrite with
  $\mathbf P$, $\mathbf P^{\mathrm{2D}}$ and the four-part structure, with eqrefs to (B10.1), (B10.4), (B10.5).
- v1 has one control; add the LO/TO seg control (quality bar: 2–5 controls).
- aria-labels: update to v2 wording ("Bulk LO wave: bound-charge planes with uniform field lines"; "Sheet LO wave: field lines
  arching out of the plane to a height of order 1/Q_∥").
- The v1 physics (potential, stream function, flux spacing, arrow directions) is correct and stays.

---

## 2. `fig-disp.html` · `viz-b10-disp` · port of v1 `viz-b-2d-disp`

**Title.** "The LO branch of a polar sheet" `<span class="tag">interactive</span>`, with a subtitle span that changes with the
view: "in-plane modes, electrostatic" / "… doped sheet, dynamic screening" / "… finite speed of light".

**What the reader must learn.** LO equals TO at $\Gamma$ for every parameter choice; the LO branch rises linearly with a slope
divided by $\varepsilon_{\mathrm{env}}$ and saturates at a plateau set by $r_0$ alone, close to the bulk LO; the bulk shift is a
constant. In a doped sheet the $\sqrt{Q_\parallel}$ plasmon anticrosses the LO branch within one sample, with a splitting that
grows as the square root of the coupling. A finite speed of light matters only inside the light cone, $Q_\parallel\lesssim10^{-4}$ Å⁻¹,
where the modes radiate.

**Include position.** In `<h4 id="b10-lo-model">`, right after the paragraph that follows (B10.14) and before
`<h4 id="b10-slope">`. The sentences just before the include:
> "and a substrate replaces $1+r_0Q_\parallel$ by $\varepsilon_{\mathrm{env}}+r_0Q_\parallel$. TO is not shifted at any $Q_\parallel$.
> Figure B10.2 plots (B10.14) for hBN and for a weakly polar MoS₂-like sheet."

Later text refers to its "doped" and "finite c" views (subsections `b10-doped` and `b10-light`).

**Physics.** Energies in eV inside the code. $C\equiv\hbar^2\,2\pi e^2Z^{*2}/(\Omega_{\mathrm S}M_{\mathrm r})=0.378205\,Z^{*2}/(\Omega_{\mathrm S}[\text{Å}^2]\,M_{\mathrm r}[\mathrm{amu}])$ eV² Å,
$\varepsilon_{\mathrm b}(Q_\parallel)\equiv\varepsilon_{\mathrm{env}}+r_0Q_\parallel$.
- *Undoped view* — (B10.14) `b10-eq-lo` with the environment of (B10.8) `b10-eq-env`:
  $(\hbar\omega_{\mathrm{LO}})^2=(\hbar\omega_{\mathrm{TO}})^2+CQ_\parallel/\varepsilon_{\mathrm b}$. TO: flat at $\hbar\omega_{\mathrm{TO}}$.
  Bulk LO: $(\hbar\omega_{\mathrm{TO}})^2+2C/(d\,\varepsilon_\infty)$ (`a4-eq-lot` with $\Omega_{\mathrm V}=\Omega_{\mathrm S}d$).
  Plateau: $(\hbar\omega_{\mathrm{TO}})^2+C/r_0$. Tangent at $\Gamma$ with slope (B10.15) `b10-eq-slope`,
  $\hbar\,\mathrm d\omega/\mathrm dQ_\parallel=C/(2\hbar\omega_{\mathrm{TO}}\varepsilon_{\mathrm{env}})$, drawn over the lower half of the plot.
  Crossover marker at $Q_\parallel=\varepsilon_{\mathrm{env}}/r_0$ (if inside the axis).
  Right panel: the squared shift $(\hbar\omega_{\mathrm{LO}})^2-(\hbar\omega_{\mathrm{TO}})^2$ against $Q_\parallel$, with the bulk constant
  $2C/(d\varepsilon_\infty)$ and the unscreened line $CQ_\parallel/\varepsilon_{\mathrm{env}}$ ($r_0=0$).
- *Doped view* — (B10.25) `b10-eq-drude` and (B10.26) `b10-eq-coupled`: $(\hbar\omega_{\mathrm p}^{\mathrm{2D}})^2=689.42\,n^{\mathrm{2D}}Q_\parallel/(m^*/m_{\mathrm e})$
  eV², roots in $x=\omega^2$ of $\varepsilon_{\mathrm b}x^2-x\,[\varepsilon_{\mathrm b}\omega_{\mathrm{TO}}^2+\omega_{\mathrm p}^2+CQ_\parallel]+\omega_{\mathrm{TO}}^2\omega_{\mathrm p}^2=0$
  (upper root with $+\sqrt{\ }$, lower as $c/(a\,x_+)$ to avoid cancellation, as v1). Curves on $0\le Q_\parallel\le0.1$ Å⁻¹: TO, undoped
  LO (dashed), uncoupled plasmon $\omega_{\mathrm p}/\sqrt{\varepsilon_{\mathrm b}}$ (dashed), the two coupled modes. Right panel: the
  phonon spectral function $-\frac1\pi\mathrm{Im}F_{\nu\nu}$ at the selected $Q_\parallel$,
  $A(\omega)=2\omega_{\mathrm{TO}}D_i/[\pi(D_r^2+D_i^2)]$ (eV⁻¹, energies in eV), with
  $D=\omega^2-\omega_{\mathrm{TO}}^2-CQ_\parallel/\varepsilon^{\mathrm{2D}}(\omega)+\mathrm i\gamma_{\mathrm{ph}}\omega$ written as $D_r+\mathrm iD_i$,
  $\varepsilon^{\mathrm{2D}}(\omega)=\varepsilon_{\mathrm b}-\omega_{\mathrm p}^2/(\omega^2+\mathrm i\gamma_{\mathrm{el}}\omega)$,
  $\hbar\gamma_{\mathrm{ph}}=1$ meV, $\hbar\gamma_{\mathrm{el}}=5$ meV (retarded convention of NOTATION §8; v1's code already does this).
  Vertical marks at the undoped LO and at the uncoupled plasmon.
- *Finite-c view* — (B10.28) `b10-eq-loc` and (B10.29) `b10-eq-rad`, sheet in a homogeneous medium $\varepsilon_{\mathrm{env}}$:
  light-line wavevector $Q_{\mathrm L}=\sqrt{\varepsilon_{\mathrm{env}}}\,\omega_{\mathrm{TO}}/c$ (in code `QL = sqrt(env)*ħω_TO/ħc`).
  For $Q_\parallel>Q_{\mathrm L}$ solve $\zeta^2=Q_\parallel^2-Q_{\mathrm L}^2-k\,\zeta/(\varepsilon_{\mathrm{env}}+r_0\zeta)$, $k=\varepsilon_{\mathrm{env}}C/(\hbar c)^2$, by
  bisection on $0<\zeta<\sqrt{Q_\parallel^2-Q_{\mathrm L}^2}$ (60 steps); shift $(\hbar\Delta\omega)$ from $(\hbar\omega)^2=(\hbar\omega_{\mathrm{TO}})^2+C\zeta/(\varepsilon_{\mathrm{env}}+r_0\zeta)$.
  For $Q_\parallel<Q_{\mathrm L}$: $q_z=\sqrt{Q_{\mathrm L}^2-Q_\parallel^2}$,
  $\hbar\gamma^{\mathrm{rad}}_{\mathrm p}=C\,\varepsilon_{\mathrm{env}}q_z/[\hbar\omega_{\mathrm{TO}}(\varepsilon_{\mathrm{env}}^2+r_0^2q_z^2)]$,
  $\hbar\gamma^{\mathrm{rad}}_{\mathrm s}=C\,\hbar\omega_{\mathrm{TO}}/[(\hbar c)^2q_z]$. At $\Gamma$ both are $C/(\hbar c\sqrt{\varepsilon_{\mathrm{env}}})$.
  Electrostatic shift for comparison: from (B10.14). Log-log axes, $Q_\parallel$ from $10^{-6}$ to $10^{-2}$ Å⁻¹, energies in μeV.
  Right panel: $\zeta/Q_\parallel$ (the reduction of the sheet kernel by the finite speed of light; 1 in the electrostatic limit),
  linear $y$ from 0 to 1.1, log $x$.

**Presets** (literature values, illustrative; the section's table uses the same numbers).

| | $\Omega_{\mathrm S}$ (Å²) | $M_{\mathrm r}$ (amu) | $Z^*$ | $\hbar\omega_{\mathrm{TO}}$ (meV) | $r_0$ (Å) | $\varepsilon_{\mathrm{env}}$ | $d$ (Å) | $\varepsilon_\infty$ | $n^{\mathrm{2D}}$ (cm⁻²) | $m^*/m_{\mathrm e}$ | spectrum $Q_\parallel$ (Å⁻¹) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| hBN monolayer | 5.43 | 6.10 | 2.7 | 170 | 10 | 1 | 3.33 | 4.87 | $10^{13}$ | 0.5 | 0.03 |
| MoS₂-like (weakly polar) | 8.65 | 38.4 | 1.0 | 48 | 40 | 1 | 6.15 | 15.0 | $10^{12}$ | 0.5 | 0.05 |

**Values to reproduce** (check the code against them). hBN: $C=0.0832$ eV² Å; slope 0.245 eV Å; $\hbar\omega_{\mathrm{LO}}(0.1\ \text{Å}^{-1})=181.8$ meV;
plateau 192.9 meV; bulk LO 197.9 meV; doped: crossing near 0.029 Å⁻¹, minimum splitting 40.6 meV; light line
$8.6\times10^{-5}$ Å⁻¹ ($2\pi/Q_{\mathrm L}=7.3$ μm), $\hbar\gamma^{\mathrm{rad}}(\Gamma)=42.2$ μeV, electrostatic shift at the light line 21 μeV.
MoS₂-like: $C=1.14\times10^{-3}$ eV² Å; slope 11.9 meV Å; 48.24 meV at 0.1 Å⁻¹; plateau 48.30 meV; bulk 48.26 meV; doped: crossing
near 0.052 Å⁻¹, splitting 4.3 meV (undoped shift there 0.2 meV); $\hbar\gamma^{\mathrm{rad}}(\Gamma)=0.58$ μeV.

**What to draw.** Two plots (`viz-grid two`): left $\hbar\omega$ (meV) against $Q_\parallel$ (Å⁻¹), right as described per view.
Colours: LO and phonon-like curves `'phonon'`, plasmon `'electron'`, bulk LO `'c5'`, TO and reference lines `'muted'`, light-cone
band `'photon'` (alpha 0.07). Text labels on the plot: "LO = TO at Γ" next to the point at $(0,\hbar\omega_{\mathrm{TO}})$, "bulk LO", and a
small grey note in a corner with the subtitle of the view. Fixed axis limits per preset and view; change them only when the
data leave the range or fill less than 40% of it (v1's `sticky` helper does this; keep it).

**Controls (at most 5 visible).**
1. View: seg control "undoped" / "doped" / "finite c" (initial "undoped").
2. Material: select "hBN monolayer" / "MoS₂-like (weakly polar)" (initial hBN); setting it loads the preset values.
3. $\varepsilon_{\mathrm{env}}$: range 1–10, step 0.05, label "ε_{env} (substrate, capping)".
4. $r_0$: range 0–80 Å, step 0.5, label "r_0 (Å)".
5. Doped view only: $n^{\mathrm{2D}}$, log range $10^{11}$–$10^{14}$ cm⁻², label "carrier density n^{2D} (cm^{-2})".
The $Q_\parallel$ of the spectrum is chosen by clicking or dragging on the left plot in the doped view (a vertical marker with two
dots on the coupled modes) and with the left/right arrow keys when the plot has focus (step 0.001 Å⁻¹, range 0.002–0.1); say so in
its aria-label. $Z^*$, $\hbar\omega_{\mathrm{TO}}$ and $m^*$ are fixed by the preset (v1 had sliders for them; drop them to stay within 5).

**Readouts.** Undoped: "ħ²·2πe²Z*²/(Ω_S M_r)" (eV² Å), slope at Γ (meV Å), $\hbar\omega_{\mathrm{LO}}$ at 0.1 Å⁻¹, plateau, crossover
$\varepsilon_{\mathrm{env}}/r_0$ (Å⁻¹), bulk LO. Doped: uncoupled plasmon at the selected $Q_\parallel$, upper and lower coupled mode, crossing
$Q_\parallel$, minimum splitting. Finite c: light line $Q_{\mathrm L}$ (Å⁻¹), $2\pi/Q_{\mathrm L}$ (μm), $\hbar\gamma^{\mathrm{rad}}$ at $\Gamma$ (μeV),
electrostatic shift at $Q_{\mathrm L}$ (μeV). Readout keys in rich text (`ħω_{LO}`, `γ^{rad}`).

**Initial state.** Undoped view, hBN, $\varepsilon_{\mathrm{env}}=1$, $r_0=10$ Å: LO leaves 170 meV linearly and bends towards 193 meV,
below the dashed bulk line at 198 meV.

**Caption parts.**
- *What to notice.* LO equals TO at $\Gamma$ for every parameter choice. Raising $\varepsilon_{\mathrm{env}}$ lowers the slope but not the
  plateau; raising $r_0$ lowers the plateau and moves the crossover $\varepsilon_{\mathrm{env}}/r_0$ to smaller $Q_\parallel$.
- *Marks.* Solid phonon-coloured curve: 2D LO, (B10.14); grey dashed: TO; dashed `c5`: bulk LO; dotted: plateau and the tangent of
  slope (B10.15); vertical line: crossover. Doped view: electron-coloured dashed curve, the uncoupled plasmon; solid curves, the
  coupled modes; right, the spectral function. Finite-c view: shaded band, inside the light cone.
- *Try.* (1) Set $\varepsilon_{\mathrm{env}}=4.5$ (hBN encapsulation): the slope drops 4.5 times, the plateau stays. (2) Move $r_0$ from
  0 to 40 Å. (3) In the doped view, drag the spectrum marker through the crossing and watch the weight pass from one peak to the
  other; switch to MoS₂-like to see a 4 meV splitting from a 0.2 meV shift. (4) In the finite-c view, read the 42 μeV width at $\Gamma$.
- *Model.* (B10.14), (B10.15), (B10.25)–(B10.29), with `a4-eq-lot` for the bulk line; presets as in the table, literature values,
  illustrative; $\varepsilon_{\mathrm{env}}$ static (for a polar substrate the value at the phonon frequency is what enters); doped hBN is
  a thought experiment; the Drude sheet holds for $\omega\gg v_{\mathrm F}Q_\parallel$; in-plane modes only.

**What changes from v1.**
- id `viz-b-2d-disp` → `viz-b10-disp`; classes `b-2d-vt`, `pA`, `pB`, `c1`, `c2`, `ro` → `b10-disp-…`.
- Plot and legend strings in source notation must change: `'TO (bare ω_σ)'` (σ is reserved for the conductivity; the lint flags it)
  → `'TO (ω_{TO})'`; `'bulk LO, 3D (cf. E.454)'` → `'bulk LO'`; `'2D LO, E.472'` → `'2D LO'`; `'plateau √((ħω_T)^2 + C/r_0)'` →
  `'plateau'`; `'Q^* = ε_{env}/r_0'` → `'ε_{env}/r_0'`; `'A_{ph}(Q, ω)…'` and `'A_{ph} (eV^{-1})'` → `'−Im F/π'`;
  `'L^+ (coupled)'`, `'L^− (coupled)'` → `'upper coupled mode'`, `'lower coupled mode'`; `'electrostatic: v_{eff} = v_Q'`,
  `'finite c: v_{eff}/v_Q = β/Q'`, `'v_{eff}(Q, ω_T)/v_Q'` → `'electrostatic: 1'`, `'finite c: ζ/Q_∥'`, `'ζ/Q_∥'`; every `Q` in axis
  labels and readouts → `Q_∥`; `ħω_T` → `ħω_{TO}`; `n` → `n^{2D}`; readout `'C'` → `'ħ²·2πe²Z*²/(Ω_S M_r)'`.
- The caption cites EPP equation numbers and uses $S_{\rm uc}$, $\mu$ (reduced mass), $\epsilon_{\rm env}$, $\epsilon_{\rm 2D}$, $\mathbf j$,
  $F_\sigma$: rewrite in v2 notation ($\Omega_{\mathrm S}$, $M_{\mathrm r}$, $\varepsilon_{\mathrm{env}}$, $\varepsilon^{\mathrm{2D}}$, $F_{\nu\nu}$) with eqrefs,
  in the four-part structure.
- Reduce the controls from eight to at most five (drop the $Z^*$, $\hbar\omega_{\mathrm T}$ and $m^*$ sliders; replace the spectrum-$Q$ slider
  by click/drag plus arrow keys on the plot).
- The v1 formulas agree with (B10.14), (B10.25)–(B10.29): keep them (constants K = 0.378214 → 0.378205, PL = 689.42, HC = 1973.27).

---

## 3. `fig-keldysh.html` · `viz-b10-keldysh` · new

**Title.** "The Keldysh interaction and the exciton it binds" `<span class="tag">interactive</span>`.

**What the reader must learn.** In a sheet the interaction between two charges is the Coulomb law of the environment,
$e^2/\varepsilon_{\mathrm{env}}r$, only at large distance; inside $r_0/\varepsilon_{\mathrm{env}}$ it grows only logarithmically and lies far
below the bare Coulomb law. This sets the exciton binding (0.57 eV in freestanding MoS₂ instead of 13.6 eV for an unscreened
2D hydrogen problem). The same function, $1/\varepsilon^{\mathrm{2D}}(Q_\parallel)$, also sets the LO dispersion: the exciton samples it
near $Q_\parallel\sim1/a_{\mathrm X}$.

**Include position.** At the start of `<h3 id="b10-keldysh">`, after the background box "Wannier excitons in two dimensions" and
before "The interaction in real space is the inverse transform of (B10.9)." The sentences just before (above the box):
> "At large distances they run mostly through the environment, where only $\varepsilon_{\mathrm{env}}$ screens. The crossover lies near
> $r_0/\varepsilon_{\mathrm{env}}$, and it is the same $r_0$ that sets the LO dispersion. Figure B10.3 shows the interaction in real space
> and the exciton it binds."

**Physics.**
- Keldysh interaction (B10.22) `b10-eq-keldysh`, computed from its integral form:
  $W^{\mathrm{e,2D}}(r)=\frac{e^2}{r_0}\int_0^\infty\exp\!\big(-\frac{\varepsilon_{\mathrm{env}}r}{r_0}\sinh u\big)\,\mathrm du$ for $r_0>0$
  (Simpson or trapezoid, 200 points on $0\le u\le u_{\max}$, $u_{\max}=\operatorname{asinh}(40\,r_0/\varepsilon_{\mathrm{env}}r)+1$);
  $W=e^2/(\varepsilon_{\mathrm{env}}r)$ for $r_0=0$. Check: MoS₂ freestanding ($r_0=40$ Å) gives $W(10\ \text{Å})=0.616$ eV against
  $e^2/r=1.440$ eV.
- Limits (B10.23) `b10-eq-klimits`: $e^2/(\varepsilon_{\mathrm{env}}r)$ for $r\gg r_0/\varepsilon_{\mathrm{env}}$;
  $(e^2/r_0)[\ln(2r_0/\varepsilon_{\mathrm{env}}r)-0.5772]$ for $r\ll r_0/\varepsilon_{\mathrm{env}}$.
- Exciton (B10.24) `b10-eq-wannier`, variational estimate of the derivation box: trial $\psi_{\mathrm X}\propto\mathrm e^{-r/a_{\mathrm X}}$ with the trial radius $a_{\mathrm X}$,
  $E(a_{\mathrm X})=\frac{\hbar^2}{2m_{\mathrm X}a_{\mathrm X}^2}-\langle W\rangle$,
  $\langle W\rangle=\frac{2e^2}{a_{\mathrm X}}\int_0^{\pi/2}\frac{\cos\theta\,\mathrm d\theta}{\varepsilon_{\mathrm{env}}+(2r_0/a_{\mathrm X})\tan\theta}$
  (100-point midpoint rule); minimize over $\ln a_{\mathrm X}\in[\ln0.5,\ln300]$ by golden section (60 steps); $E_{\mathrm b}=-\min E$, and the
  exciton radius is the minimizing $a_{\mathrm X}$. With $\hbar^2/2m_{\mathrm e}=3.80998$ eV Å²,
  $\hbar^2/(2m_{\mathrm X}a_{\mathrm X}^2)=3.80998/(m_{\mathrm X}/m_{\mathrm e})/a_{\mathrm X}^2$ eV (lengths in Å). Do not call the trial radius $a$: in
  the section $a$ is the lattice constant.
  Reference without sheet screening ($r_0=0$): $E_{\mathrm b}=2m_{\mathrm X}e^4/(\varepsilon_{\mathrm{env}}^2\hbar^2)=54.42\ \mathrm{eV}\times(m_{\mathrm X}/m_{\mathrm e})/\varepsilon_{\mathrm{env}}^2$.
- Radial probability of the trial state: $2\pi r|\psi_{\mathrm X}|^2=(4r/a_{\mathrm X}^2)\,\mathrm e^{-2r/a_{\mathrm X}}$.
- Values to reproduce (variational; the exact values of the text are 2–3% higher): MoS₂ ($r_0=40$ Å, $m_{\mathrm X}=0.25$):
  $E_{\mathrm b}=0.558$ eV, $a_{\mathrm X}=10.1$ Å at $\varepsilon_{\mathrm{env}}=1$; 0.315 eV at 2.45; 0.187 eV at 4.5. hBN ($r_0=10$ Å,
  $m_{\mathrm X}=0.5$): 1.83 eV, 3.7 Å. $r_0=0$, $m_{\mathrm X}=0.25$, $\varepsilon_{\mathrm{env}}=1$: 13.6 eV, 1.06 Å.

**What to draw.** Two plots (`viz-grid two`).
- Left: log-log, $r$ from 1 to 1000 Å, $W$ from $10^{-3}$ to 20 eV (fixed). Curves: Keldysh $W^{\mathrm{e,2D}}(r)$ (solid, `'electron'`);
  environment Coulomb $e^2/(\varepsilon_{\mathrm{env}}r)$ (dashed, `'c5'`); bare $e^2/r$ (thin dotted, `'muted'`); the logarithmic asymptote of
  (B10.23), dotted, drawn only for $r<r_0/\varepsilon_{\mathrm{env}}$ and where positive. Vertical line at $r_0/\varepsilon_{\mathrm{env}}$ labelled
  `r_0/ε_{env}` (hidden for $r_0=0$). Vertical line at $a_{\mathrm X}$ labelled `a_X`, and the radial probability as a light filled
  curve (`'electron'`, alpha 0.15) scaled to the upper third of the panel, labelled "exciton (arb. units)". Axis labels
  `r (Å)` and `W^{e,2D}(r) (eV)`.
- Right: linear, $Q_\parallel$ from 0 to 0.5 Å⁻¹, $y$ from 0 to 1.05: the curve $1/\varepsilon^{\mathrm{2D}}(Q_\parallel)=1/(\varepsilon_{\mathrm{env}}+r_0Q_\parallel)$
  (solid, `'phonon'`), the value $1/\varepsilon_{\mathrm{env}}$ as a dashed horizontal line, a vertical line at the crossover
  $\varepsilon_{\mathrm{env}}/r_0$, and a shaded band from $0.5/a_{\mathrm X}$ to $2/a_{\mathrm X}$ labelled "exciton". A short note on the plot:
  "the LO shift at Q_∥ uses the same curve". Axis labels `Q_∥ (Å^{-1})` and `1/ε^{2D}(Q_∥)`.

**Controls (4).**
1. Material: select "MoS₂" ($r_0=40$ Å, $m_{\mathrm X}=0.25$) / "hBN" ($r_0=10$ Å, $m_{\mathrm X}=0.5$); loads $r_0$ and $m_{\mathrm X}$.
2. $r_0$: range 0–80 Å, step 0.5, label "r_0 (Å)".
3. $\varepsilon_{\mathrm{env}}$: range 1–10, step 0.05, label "ε_{env}" (show "1 freestanding, 2.45 SiO₂, 4.5 hBN" in the label's
   value format or as tick hints).
4. $m_{\mathrm X}$: range 0.1–1.0, step 0.01, label "m_X / m_e".

**Readouts.** $E_{\mathrm b}$ (variational, meV below 1 eV, else eV), $a_{\mathrm X}$ (Å), $r_0/\varepsilon_{\mathrm{env}}$ (Å),
"2D hydrogen (no sheet screening)" $E_{\mathrm b}$, and the ratio of the two.

**Initial state.** MoS₂, $r_0=40$ Å, $\varepsilon_{\mathrm{env}}=1$, $m_{\mathrm X}=0.25$: the Keldysh curve leaves the Coulomb line near
40 Å and flattens; $E_{\mathrm b}=0.56$ eV against 13.6 eV, $a_{\mathrm X}\approx10$ Å.

**Caption parts.**
- *What to notice.* Far apart, the two charges interact as in the environment; inside $r_0/\varepsilon_{\mathrm{env}}$ the interaction grows
  only logarithmically, far below the bare Coulomb law, and the exciton is bound 24 times more weakly than an unscreened 2D
  hydrogen atom.
- *Marks.* Solid: Keldysh interaction (B10.22); dashed: $e^2/\varepsilon_{\mathrm{env}}r$; dotted: bare $e^2/r$ and the logarithmic limit
  (B10.23); vertical lines: $r_0/\varepsilon_{\mathrm{env}}$ and $a_{\mathrm X}$; shaded: the exciton's radial probability. Right: the shared
  screening function and the range of $Q_\parallel$ the exciton samples.
- *Try.* (1) Raise $\varepsilon_{\mathrm{env}}$ to 2.45 (SiO₂) and 4.5 (hBN encapsulation): the binding falls to 0.32 and 0.19 eV and the
  crossover moves inwards. (2) Set $r_0=0$: the 2D hydrogen result, $2m_{\mathrm X}e^4/(\varepsilon_{\mathrm{env}}^2\hbar^2)$. (3) Choose hBN: a
  compact exciton of 3.7 Å bound by 1.8 eV, where the Wannier model is only qualitative.
- *Model.* (B10.22)–(B10.24) and the variational estimate of the derivation box (a lower bound on $E_{\mathrm b}$, 2–3% below the exact
  value); strictly 2D sheet, static screening; $r_0$, $m_{\mathrm X}$ and $\varepsilon_{\mathrm{env}}$ are literature values, illustrative.

Cost: about 150 log-spaced $r$ values × 200 quadrature points plus the minimization, well under 15 ms per redraw.
