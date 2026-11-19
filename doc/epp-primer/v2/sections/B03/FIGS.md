# B.3 figures · specifications for the figure authors

Section: `v2/sections/B03-dielectric.html` (id `b-dielectric`, prefix `b3`). Four figures, in reading order:

| # | File | Id | Kind | v1 source |
|---|---|---|---|---|
| B3.1 | `fig-cell.html` | `viz-b3-cell` | new | — |
| B3.2 | `fig-lt.html` | `viz-b3-lt` | port | `viz-b-dielectric-lt` in `doc/epp-primer/sections/B03-dielectric.html` (static SVG, lines 56–208) |
| B3.3 | `fig-local.html` | `viz-b3-local` | new | reuses the microscopic mode and the level diagram of `viz-b-dielectric-lab` (same v1 file, script lines 268–497) |
| B3.4 | `fig-lab.html` | `viz-b3-lab` | port | `viz-b-dielectric-lab` in the same v1 file (figure lines 235–245, script lines 268–497) |

Read `v2/build/BUILDER_GUIDE.md` §5 first. Common rules for all four:

- Notation exactly as in the section: $\omega_{\mathrm{TO}}$, $\omega_{\mathrm{LO}}$, $\omega_{\mathrm p}^{\mathrm{ion}}$, $\varepsilon(\omega)$, $\varepsilon_\infty$, $\varepsilon_{\mathrm s}$, $\bP^{\mathrm{ion}}$, $\bP^{\mathrm e}$, $\mathbf D$, $\bE$, $\bE^{\mathrm{loc}}$, $Z^*$ (dimensionless), $\OV$ (cell volume), $M_{\mathrm r}$ (reduced mass), $K$ and $\omega_0$ (the short-range spring of Section B.2 and its resonance), $\bu_{\mathrm{rel}}$ and $\xi_{\mathrm{rel}}$ (relative displacement, as in B.2), $\zeta$, $W$, $W^{\mathrm e}$, $W^{\mathrm{ph}}$, $v$, $\gamma$ (damping rate), $R$. Canvas text with `EPP.drawRich`: `ω_{TO}`, `ω_{LO}`, `ω_p^{ion}`, `ε_∞`, `ε_s`, `P^{ion}`, `P^{e}`, bold vectors with `{bold:true}`. Never: `ε` written as the source's ϵ, `ω̃`, `ω_O`, `Ω_uc`, `W_e`, `W_ph`, `|e|` after $Z^*$, `M (u)`, equation numbers of the form `E.159`, `σ` for anything.
- Units: energies in meV on axes; cm⁻¹ only as a secondary readout or hover value (1 meV = 8.06554 cm⁻¹). Lengths in Å or pm, masses in amu, fields in kV/cm.
- Physics constants: $e^2=14.3996$ eV Å; $\hbar\omega=64.654\ \text{meV}\times\sqrt{K[\text{eV/Å}^2]/M[\text{amu}]}$; hence $(\hbar\omega_{\mathrm p}^{\mathrm{ion}})^2=64.654^2\times4\pi\cdot14.3996\,Z^{*2}/(\OV M_{\mathrm r})=7.564\times10^5\ \text{meV}^2\times Z^{*2}/(\OV[\text{Å}^3]M_{\mathrm r}[\text{amu}])$.
- Material data (literature values, illustrative; the same numbers as the section text and its footnote `fn-b3-neutral`):

| | $\hbar\omega_{\mathrm{TO}}$ (meV) | $\hbar\omega_{\mathrm{LO}}$ (meV) | $\varepsilon_\infty$ | $\varepsilon_{\mathrm s}$ (meas.) | $Z^*$ | $\OV$ (Å³) | $M_1$, $M_2$ (amu) | $M_{\mathrm r}$ (amu) | $\hbar\gamma$ (meV) |
|---|---|---|---|---|---|---|---|---|---|
| GaAs | 33.3 (268.6 cm⁻¹) | 36.2 (292 cm⁻¹) | 10.9 | 12.9 | 2.07 | 45.17 | Ga 69.72, As 74.92 | 36.12 | 0.31 |
| NaCl | 20.33 (164 cm⁻¹) | 32.73 (264 cm⁻¹) | 2.34 | 5.9 | 1.12 | 44.85 | Na 22.99, Cl 35.45 | 13.95 | 0.99 |
| 3C-SiC | 98.69 (796 cm⁻¹) | 120.51 (972 cm⁻¹) | 6.52 | 9.72 | 2.70 | 20.72 | Si 28.09, C 12.01 | 8.41 | 0.62 |
| hBN, in plane | 169.86 (1370 cm⁻¹) | 199.61 (1610 cm⁻¹) | 4.87 | 6.70 | 2.71 | 18.08 | B 10.81, N 14.01 | 6.10 | 0.62 |

  Derived checks (the figures must reproduce them): GaAs $\hbar\omega_{\mathrm p}^{\mathrm{ion}}=44.57$ meV from $Z^*$; $\omega_{\mathrm{LO}}$ from $Z^*$ via (B3.10) = 35.93 meV; LST $\varepsilon_{\mathrm s}=12.88$; $K=M_{\mathrm r}\omega_{\mathrm{TO}}^2=9.58$ eV/Å². NaCl $\hbar\omega_{\mathrm p}^{\mathrm{ion}}=38.95$ meV, LST 6.06.
- Equation references in captions use `<a class="eqref" href="#b3-eq-…"></a>`; the ids are listed under each figure.
- CSS classes prefixed with the figure: `b3-cell-…`, `b3-lt-…`, `b3-local-…`, `b3-lab-…`.
- Colours read at draw time with `EPP.color(...)`; suggested roles: ions and displacements `phonon` (amber); polarization `c5` (violet); electric field $\bE$ `photon` (rose); electronic polarization `electron` (teal); $\mathbf D$ `muted`; text `ink`.
- Test with `bash v2/build/figtest.sh v2/sections/B03/fig-<name>.html`, then open both screenshots.

---

## B3.1 · `fig-cell.html` · `viz-b3-cell` · new

**Title.** "Two masses, one spring and one field" `<span class="tag">interactive</span>`

**What the reader must learn.**
1. At wavelengths much longer than the cell every cell moves in step: the crystal is one two-atom molecule per cell, repeated.
2. A uniform field pushes cation and anion with equal and opposite forces, so the centre of mass of a neutral cell stays at rest and only the bond stretch $\bu_{\mathrm{rel}}=\bu_1-\bu_2$ responds: one damped oscillator of mass $M_{\mathrm r}$ and stiffness $K$, resonant at $\omega_{\mathrm{TO}}$.
3. The lighter ion moves more: $\boldsymbol\xi_1/\boldsymbol\xi_2=-M_2/M_1$. The dipole of a cell is $eZ^*\bu_{\mathrm{rel}}$.
4. If the Born charges did not add up to zero (an unphysical setting), the centre of mass would be driven like a free charge, with an amplitude that grows as $1/\omega^2$: the reason for the sum rule.

**Include position.** In subsection "All cells in step" (`#b3-q0`), right after the paragraph that ends: "Inside one region the crystal behaves as a single two-atom molecule per cell, repeated in lockstep. Figure B3.1 shows such a region driven by an oscillating field." The reduction steps and (B3.1) follow it.

**Physics** (section equations `b3-eq-eom` (B3.1), `b3-eq-rel` (B3.2), `b3-eq-xi` (B3.3), `b3-eq-pion` (B3.4); derivation box "the driven two-mass problem with arbitrary charges").
- Charges $Z^*_1=Z^*$ (cation) and $Z^*_2=-Z^*+\delta$, where $\delta=Z^*_1+Z^*_2$ is the charge imbalance (0 is physical).
- Field $E(t)=E_0\cos\omega t$ along $x$, $E_0=100$ kV/cm $=10^{-3}$ V/Å (force $eZ^*E_0$ in eV/Å when $E_0$ is in V/Å).
- Relative and centre-of-mass amplitudes (complex, $\ee^{-\ii\omega t}$ convention; damping on the relative coordinate only, as in (B3.2)):
  - $\xi_{\mathrm{rel}}=e\,(Z^*_1/M_1-Z^*_2/M_2)\,E_0/(\omega_{\mathrm{TO}}^2-\omega^2-\ii\gamma\omega)$, with $\omega_{\mathrm{TO}}^2=K/M_{\mathrm r}$; for $\delta=0$ this is (B3.3).
  - $\xi_{\mathrm c}=-e\,\delta\,E_0/(M_{\mathrm c}\,\omega^2)$, $M_{\mathrm c}=M_1+M_2$ (free acceleration; plot only for $\omega\ge1$ meV).
  - $\xi_1=\xi_{\mathrm c}+(M_2/M_{\mathrm c})\,\xi_{\mathrm{rel}}$, $\xi_2=\xi_{\mathrm c}-(M_1/M_{\mathrm c})\,\xi_{\mathrm{rel}}$.
  - Dipole per cell $p=e(Z^*_1\xi_1+Z^*_2\xi_2)$.
  - Physical displacement $u_\kappa(t)=\mathrm{Re}[\xi_\kappa\ee^{-\ii\omega t}]$.
- Units: with $\omega$ in meV use $\omega^2$ in meV² and convert $e E_0/M$ through $\hbar^2/(\text{amu Å}^2)=4.1802$ meV: $\xi[\text{Å}]=\dfrac{(eZE_0)[\text{eV/Å}]\cdot 4180.2}{M[\text{amu}]\,(\hbar\omega_0^2-\hbar\omega^2-\ii\hbar\gamma\,\hbar\omega)[\text{meV}^2]}$, where $4180.2=4.1802\ \text{meV}\times1000$ converts eV to meV. Check: static GaAs $\xi_{\mathrm{rel}}=eZ^*E_0/K=2.07\times10^{-3}/9.58=2.16\times10^{-4}$ Å $=0.022$ pm; at resonance about 2.3 pm.
- Material: GaAs ($M_1=69.72$, $M_2=74.92$ amu, $Z^*=2.07$, $K=9.58$ eV/Å², the short-range spring of Section B.2, so that $\hbar\omega_0=\hbar\omega_{\mathrm{TO}}=33.3$ meV with only the macroscopic field acting, $\hbar\gamma=0.31$ meV). The mass-ratio control keeps $M_{\mathrm r}=36.12$ amu and $K$ fixed: $M_1=M_{\mathrm r}(1+\rho)$, $M_2=M_{\mathrm r}(1+1/\rho)$ with $\rho=M_1/M_2$, so $\omega_{\mathrm{TO}}$ does not move.

**What to draw.**
- Top panel (geometric, canvas, about 2:1): a row of 5 cells along $x$ (thin dashed cell boundaries). In each cell a cation (filled disc, "+", label `Ga` in the middle cell) and an anion (open disc, "−", `As`), joined by a zigzag spring labelled `K` in the middle cell. All five cells move identically (lockstep is the point).
- A large field arrow above the row, labelled `𝐄(t)` (bold), whose length and direction follow $E_0\cos\omega t$.
- Per ion, an amber displacement arrow from the equilibrium position (thin tick) to the current position, scaled by a visual gain so the largest displacement in the current state is about 0.3 of the ion spacing; state "displacements exaggerated" in the caption.
- Per cell, a small cross at the centre of mass: it stays fixed for $\delta=0$ and swings when $\delta\neq0$.
- Per cell, a violet dipole arrow `p` (below the ions), length $\propto$ the instantaneous dipole.
- Bottom panel (plot, `EPP.Plot`, log $y$): $|\xi_{\mathrm{rel}}|$ (ink, label `|ξ_{rel}|`) and $|\xi_{\mathrm c}|$ (rose, label `|ξ_c|`, drawn only when $\delta\neq0$) in pm against $\hbar\omega$ from 1 to 60 meV; fixed limits $y\in[10^{-3},10^{2}]$ pm; dashed vertical line at $\omega_{\mathrm{TO}}$ labelled `ω_{TO}`; a point marker at the current drive frequency.

**Controls (4).**
1. Drive energy $\hbar\omega$: 1–60 meV, step 0.1, default 25 meV (also settable by clicking the plot).
2. Mass ratio $M_1/M_2$: 0.2–5 (logarithmic), default 0.931 (GaAs).
3. Charge imbalance $Z^*_1+Z^*_2$: −0.3 to +0.3, step 0.01, default 0, labelled "(unphysical unless 0)".
4. Play/Pause through `EPP.anim` (reduced motion honoured; the static frame shows $t=0$ and a quarter period later as faint ghosts).

**Readouts** (`EPP.readout`, MathJax labels): $\hbar\omega_{\mathrm{TO}}$; $|\xi_1|$, $|\xi_2|$ (pm, for $E_0=100$ kV/cm); $\xi_1/\xi_2$ (real when $\gamma\to0$; show as a signed number, expected $-M_2/M_1$ for $\delta=0$); phase lag of $u_{\mathrm{rel}}$ behind $E$ (degrees); $|\xi_{\mathrm c}|$ (pm); dipole per cell $|p|$ in $e\,$pm.

**Initial state.** GaAs, $\hbar\omega=25$ meV (below resonance: ions in phase with the force, clearly visible), $\delta=0$, animation playing.

**Caption key points.**
- *What to notice.* All cells move together; cation and anion swing against each other while each cell's centre of mass stays still. The field drives only the bond stretch, a single oscillator of reduced mass with its resonance at $\hbar\omega_{\mathrm{TO}}=33.3$ meV.
- *Marks.* Filled/open discs: cations/anions; zigzag: the effective spring $K$; rose arrow: field; amber arrows: displacements (exaggerated); violet arrows: cell dipoles; crosses: centres of mass. Plot: $|\xi_{\mathrm{rel}}|$ ink, $|\xi_{\mathrm c}|$ rose, dashed line $\omega_{\mathrm{TO}}$, dot: current drive.
- *Try.* (1) Sweep $\hbar\omega$ through 33.3 meV: the stretch peaks, the lag passes 90°, and above resonance the ions move against the force. (2) Change the mass ratio: the lighter ion moves more, $\xi_1/\xi_2=-M_2/M_1$, while the resonance stays put. (3) Give the cell a net charge: the centre of mass starts to swing, with an amplitude rising as $1/\omega^2$ toward low frequency like a free charge, which an insulator cannot have.
- *Model.* <a class="eqref" href="#b3-eq-eom"></a>–<a class="eqref" href="#b3-eq-xi"></a> with damping on the relative coordinate; GaAs values above; $E_0=100$ kV/cm.

**Notes for the author.** Good model: `A01/fig-diatomic.html` (plot + animation + readouts). Keep the top panel within 390 px: 5 cells at phone width means about 70 px per cell; use 3 cells below 420 px if labels collide. `touch-action` only if you make anything draggable (nothing needs to be).

---

## B3.2 · `fig-lt.html` · `viz-b3-lt` · port of `viz-b-dielectric-lt`

**Title.** "Why the LO wave carries a field and the TO wave does not" `<span class="tag">interactive</span>`

**What the reader must learn.**
1. Cut the crystal into slabs perpendicular to $\bq$. In a longitudinal (LO) wave $\bP^{\mathrm{ion}}$ is normal to the slab faces, so the faces carry polarization charge ($-\nabla\cdot\bP\neq0$), and the field of these charges, $\bE=-4\pi\bP^{\mathrm{ion}}/\varepsilon_\infty$, points against $\bP^{\mathrm{ion}}$ everywhere: $\mathbf D=0$, an extra restoring force, a higher frequency.
2. In a transverse (TO) wave $\bP^{\mathrm{ion}}$ lies in the faces, the faces stay neutral, $\bE=0$ and $\mathbf D=4\pi\bP^{\mathrm{ion}}$: the springs act alone.
3. The electrons screen the LO field by $\varepsilon_\infty$, which sets the size of the LO–TO splitting.

**Include position.** Subsection "Zeros are LO, poles are TO" (`#b3-lt`), right after the slab paragraph, which ends: "In a transverse wave the ions move perpendicular to $\bq$, $\bP$ lies in the slab faces, the faces stay neutral, and no field arises: the springs act alone. Section A.4 drew the Born dipoles of the same two waves (Figure A4.x)." Maxwell's conditions (B3.8) follow it.

**Physics** (`b3-eq-max` (B3.8), `b3-eq-Elo` (B3.9), `b3-eq-lo` (B3.10), `b3-eq-pion` (B3.4)).
- One wavelength along $x$ ($\bq\parallel\hat x$): $P^{\mathrm{ion}}(x,t)=P_0\cos(qx-\omega t)$, along $\hat x$ (LO) or $\hat y$ (TO).
- LO: polarization charge density $-\partial P_x/\partial x$; with $N$ discrete slabs, face charge between slabs $j$ and $j+1$ is $\sigma_j=P_j-P_{j+1}$ (sign shown). $E_x=-4\pi P^{\mathrm{ion}}_x/\varepsilon_\infty$, $P^{\mathrm e}=(\varepsilon_\infty-1)E/4\pi$, $D=E+4\pi(P^{\mathrm{ion}}+P^{\mathrm e})=0$.
- TO: $E=0$, $P^{\mathrm e}=0$, $\mathbf D=4\pi\bP^{\mathrm{ion}}$.
- Frequencies for the readout (GaAs-like): $\hbar\omega_{\mathrm{TO}}=33.3$ meV fixed, $\hbar\omega_{\mathrm p}^{\mathrm{ion}}=44.57$ meV fixed, $\omega_{\mathrm{LO}}^2=\omega_{\mathrm{TO}}^2+(\omega_{\mathrm p}^{\mathrm{ion}})^2/\varepsilon_\infty$: 35.9 meV at $\varepsilon_\infty=10.9$, 55.6 meV at $\varepsilon_\infty=1$. Added stiffness $4\pi e^2Z^{*2}/(\OV\varepsilon_\infty)$ = 1.57 eV/Å² at 10.9, 17.2 eV/Å² at 1 ($Z^*=2.07$, $\OV=45.17$ Å³).
- Electrostatic limit ($cq\gg\omega$) throughout; say so in the Model line.

**What to draw** (canvas preferred, so that it scales to 390 px without a horizontal scrollbar).
- Two panels stacked: "LO" on top, "TO" below (side by side above 900 px is optional). Each panel: an arrow `𝐪` pointing right; a strip one wavelength long divided into $N$ thin vertical slabs (light outlines).
- Inside each slab two ions (cation filled "+", anion open "−") displaced in antiphase along $x$ (LO) or $y$ (TO), with amber displacement arrows.
- Below the strip, rows of arrows sampled at the slab centres, each row labelled at its left end: `𝐏^{ion}` (violet); `𝐏^{e}` (teal; LO only, otherwise "0"); `𝐄` (rose); `𝐃` (grey). For LO the $\mathbf D$ row shows "𝐃 = 0"; for TO the $\bE$ row shows "𝐄 = 0" and the $\mathbf D$ row arrows equal $4\pi\bP^{\mathrm{ion}}$.
- LO only: on each slab face a column of small "+" or "−" signs (rose for +, blue `electron` for −), the number of signs proportional to $|\sigma_j|$ (at most 4). In "continuum" mode replace the signs by a smooth colour band under the strip showing $-\partial P_x/\partial x$ (diverging colours, zero transparent).
- Arrow lengths share one scale per panel; the $\bE$ arrows in LO have length $4\pi|P^{\mathrm{ion}}|/\varepsilon_\infty$ on the same scale as $4\pi P$, so lowering $\varepsilon_\infty$ visibly lengthens them.

**Controls (3).**
1. Electronic screening $\varepsilon_\infty$: 1–15, step 0.1, default 10.9 (GaAs).
2. Slabs / continuum: segmented control, default "slabs" with $N=12$.
3. Play/Pause (`EPP.anim`): the wave travels to the right (phase $qx-\omega t$ advancing); the face charges travel with the nodes of $\bP^{\mathrm{ion}}$. Static frame at phase 0 under reduced motion.

**Readouts.** LO: $E/(4\pi P^{\mathrm{ion}})=-1/\varepsilon_\infty$; added stiffness (eV/Å²); $\hbar\omega_{\mathrm{LO}}$ (meV, from (B3.10) with the GaAs values above). TO: $E=0$; $\hbar\omega_{\mathrm{TO}}=33.3$ meV.

**Initial state.** $\varepsilon_\infty=10.9$, slabs, playing.

**Caption key points.**
- *What to notice.* In the LO wave the polarization meets the slab faces head-on, so every face carries polarization charge, and the field of these charges points against $\bP^{\mathrm{ion}}$ in every slab: $\mathbf D=0$ and the ions feel an extra restoring force. In the TO wave the faces stay neutral and there is no field.
- *Marks.* Discs: ions; amber: displacements; violet: $\bP^{\mathrm{ion}}$; teal: $\bP^{\mathrm e}$; rose: $\bE$; grey: $\mathbf D$; signs on faces: polarization charge.
- *Try.* (1) Lower $\varepsilon_\infty$ to 1: the LO field grows to $-4\pi\bP^{\mathrm{ion}}$ and $\hbar\omega_{\mathrm{LO}}$ rises to 55.6 meV. (2) Switch to continuum: the face charges merge into the density $-\nabla\cdot\bP$, largest where $\bP$ changes sign. (3) Play: the charges travel with the wave.
- *Model.* <a class="eqref" href="#b3-eq-max"></a>–<a class="eqref" href="#b3-eq-lo"></a>, electrostatic limit $cq\gg\omega$; GaAs-like $\hbar\omega_{\mathrm{TO}}=33.3$ meV, $\hbar\omega_{\mathrm p}^{\mathrm{ion}}=44.6$ meV.

**What must change from v1.**
- v1 is a static SVG with `min-width:420px` inside an `overflow-x:auto` wrapper, so phones scroll sideways: rebuild responsively (canvas or a scalable SVG without min-width).
- v1 uses the source notation $\epsilon(\tilde\omega_{\rm LO})$, $\tilde\omega_{\rm TO}$ and $\rho_{\rm b}$: write $\varepsilon$, $\omega_{\mathrm{LO}}$, $\omega_{\mathrm{TO}}$ and "polarization charge $-\nabla\cdot\bP$".
- v1 draws $\bE=-4\pi\bP$ for $\varepsilon_\infty=1$ only and has no $\mathbf D$ or $\bP^{\mathrm e}$: add the $\varepsilon_\infty$ control, the electronic polarization and the $\mathbf D$ row.
- v1 shows a chain with charge rectangles at the nodes; v2 asks for the slab geometry of the text (charges on slab faces).
- v1's caption is a single paragraph; rewrite as What to notice / Marks / Try / Model. Rename ids and classes (`b-dielectric-…` → `b3-lt-…`).

---

## B3.3 · `fig-local.html` · `viz-b3-local` · new

**Title.** "Local field: what depends on the model and what does not" `<span class="tag">interactive</span>`

**What the reader must learn.**
1. Lorentz's construction: the field at an ion is the macroscopic field plus $4\pi\bP/3$ from the polarization charge on the wall of a spherical cavity; the point dipoles inside the sphere cancel at the centre on cubic sites.
2. A local field proportional to $\bP$ only renormalizes the spring: $\omega_{\mathrm{TO}}^2=\omega_0^2-\zeta(\omega_{\mathrm p}^{\mathrm{ion}})^2$ (B3.14).
3. The measurable $\omega_{\mathrm{TO}}$, $\omega_{\mathrm{LO}}$, $\omega_{\mathrm{LO}}^2-\omega_{\mathrm{TO}}^2=(\omega_{\mathrm p}^{\mathrm{ion}})^2/\varepsilon_\infty$ and $\varepsilon_{\mathrm s}$ do not depend on $\zeta$; only the unmeasurable $\omega_0$ does. With the EPP choice $\zeta=1/3$ and $\varepsilon_\infty>3$, $\omega_{\mathrm{LO}}$ falls below $\omega_0$.
4. If $\zeta(\omega_{\mathrm p}^{\mathrm{ion}})^2>\omega_0^2$ the TO mode goes soft ($\omega_{\mathrm{TO}}^2<0$): Cochran's ferroelectric instability, with $\varepsilon_{\mathrm s}\to\infty$ on approach.

**Include position.** Subsection "Local field versus macroscopic field" (`#b3-local`), after the paragraph on Lorentz's construction, which ends: "Inside the sphere, keep the point dipoles of the lattice: on sites of cubic symmetry their fields at the centre cancel." The next line is "The local field at the ion is therefore the Lorentz field" with (B3.13).

**Physics** (`b3-eq-pion` (B3.4), `b3-eq-loc` (B3.13), `b3-eq-zeta` (B3.14), `b3-eq-lo` (B3.10), `b3-eq-lst` (B3.12); the Clausius–Mossotti derivation box `#b3-cm`).
- $\zeta=0$ ("local field inside $K$"), $\zeta=1/3$ ("Lorentz field, EPP"), $\zeta=1/(\varepsilon_\infty+2)$ ("Clausius–Mossotti").
- For a preset, $\omega_0$ is inferred so that the preset's measured $\omega_{\mathrm{TO}}$ is reproduced at the preset's $Z^*$: $\omega_0^2=\omega_{\mathrm{TO,meas}}^2+\zeta\,(\omega_{\mathrm p}^{\mathrm{ion}})^2$. Then, with $Z^*$ (hence $\omega_{\mathrm p}^{\mathrm{ion}}$) varied at fixed $\omega_0$: $\omega_{\mathrm{TO}}^2=\omega_0^2-\zeta x^2$, $\omega_{\mathrm{LO}}^2=\omega_{\mathrm{TO}}^2+x^2/\varepsilon_\infty$ with $x=\omega_{\mathrm p}^{\mathrm{ion}}$; $\varepsilon_{\mathrm s}=\varepsilon_\infty\omega_{\mathrm{LO}}^2/\omega_{\mathrm{TO}}^2$; unstable where $\zeta x^2>\omega_0^2$.
- GaAs checks: $\hbar\omega_0=33.30$ ($\zeta=0$), 35.54 (CM), 42.09 meV ($\zeta=1/3$); $\hbar\omega_{\mathrm{LO}}=35.93$ meV in all three at $Z^*=2.07$; Szigeti charge $3Z^*/(\varepsilon_\infty+2)=0.48$. NaCl: $\omega_0=20.33$, 27.62, 30.32 meV; $\omega_{\mathrm{LO}}=32.58$ meV (above $\omega_0$ even for $\zeta=1/3$, since $\varepsilon_\infty<3$).

**What to draw.**
- Left (or top on phones) — geometric panel, canvas: a square patch of a cubic lattice (about 9×9 sites) of point dipoles, all violet arrows along $+\hat z$ (up). A circle (the Lorentz sphere, cross-section) around the central site, radius about 3 lattice spacings. On the circle, polarization charge on the cavity wall: "−" signs on the upper arc (the wall toward which the dipoles point) and "+" signs on the lower arc, density $\propto|\cos\theta|$, none at the equator. Charge sign rule: the medium's surface charge is $\bP\cdot\hat{\mathbf n}$ with $\hat{\mathbf n}$ pointing out of the medium, into the cavity; at the top of the cavity this is $-P$. Inside the circle, the lattice dipoles are drawn lighter, with a note "fields cancel at the centre (cubic sites)". At the centre, three arrows from the central ion: $\bE$ (rose, label `𝐄`), $\tfrac{4\pi}{3}\bP$ (violet, label `4π𝐏/3`), and their sum $\bE^{\mathrm{loc}}$ (ink, label `𝐄^{loc}`).
  - The model control changes the panel: $\zeta=0$ greys out the cavity construction and labels the central ion "local field already inside $K$"; $\zeta=1/3$ shows the full construction; Clausius–Mossotti adds small teal ellipses (polarizable electron clouds) at the sites and labels the ion charge $eZ_{\mathrm S}$.
- Right (or below) — level diagram, `EPP.Plot`: $x=\hbar\omega_{\mathrm p}^{\mathrm{ion}}$ (meV), $y=\hbar\omega$ (meV). Curves for the selected model: $\omega_0$ (grey dashed horizontal, label `ω_0`), $\omega_{\mathrm{TO}}(x)$ (ink, `ω_{TO}`), $\omega_{\mathrm{LO}}(x)$ (amber, `ω_{LO}`). A vertical line at the current $x$ with dots at the computed $\omega_{\mathrm{TO}}$, $\omega_{\mathrm{LO}}$, and hollow markers at the preset's measured values. Arrows from $\omega_0$ to $\omega_{\mathrm{TO}}$ (down) and to $\omega_{\mathrm{LO}}$ at the current $x$. A rose-shaded band where $\zeta x^2>\omega_0^2$ labelled "unstable (soft TO mode)". Fixed axes per preset: GaAs $x\in[0,100]$, $y\in[0,60]$ meV; NaCl $x\in[0,80]$, $y\in[0,50]$; 3C-SiC $x\in[0,320]$, $y\in[0,180]$.

**Controls (4).**
1. Local-field model: segmented, "in $K$ (ζ = 0)" / "Lorentz, EPP (ζ = 1/3)" / "Clausius–Mossotti"; default Lorentz.
2. Material preset: GaAs (default), NaCl, 3C-SiC (cubic crystals only; the Lorentz factor needs cubic sites).
3. Born charge $Z^*$: 0–4, step 0.01, default the preset's (moves $x$ at fixed $\omega_0$; at large $Z^*$ the Lorentz model goes unstable).
4. Electronic screening $\varepsilon_\infty$: 1–15, default the preset's (changes $\omega_{\mathrm{LO}}$ and, in the CM model, $\zeta$).

**Readouts.** $\hbar\omega_{\mathrm p}^{\mathrm{ion}}$; $\hbar\omega_0$; $\hbar\omega_{\mathrm{TO}}$ ("unstable" if $\omega_{\mathrm{TO}}^2<0$); $\hbar\omega_{\mathrm{LO}}$; $(\hbar\omega_{\mathrm{LO}})^2-(\hbar\omega_{\mathrm{TO}})^2=(\hbar\omega_{\mathrm p}^{\mathrm{ion}})^2/\varepsilon_\infty$ in meV² (shown with the note "same for every ζ"); $\varepsilon_{\mathrm s}$ (LST); in CM mode the Szigeti charge $Z_{\mathrm S}$.

**Initial state.** GaAs, Lorentz model ($\zeta=1/3$), $Z^*=2.07$, $\varepsilon_\infty=10.9$: the diagram shows $\omega_0=42.1$ meV above both $\omega_{\mathrm{TO}}=33.3$ and $\omega_{\mathrm{LO}}=35.9$ meV.

**Caption key points.**
- *What to notice.* Switching the local-field model moves only $\omega_0$, the frequency of the springs alone; $\omega_{\mathrm{TO}}$, $\omega_{\mathrm{LO}}$, their difference of squares and $\varepsilon_{\mathrm s}$ stay where experiment puts them. With the Lorentz field and GaAs's $\varepsilon_\infty=10.9$, $\omega_{\mathrm{LO}}$ even lies below $\omega_0$.
- *Marks.* Left: violet dipoles; circle: the Lorentz sphere; ± signs: polarization charge on the cavity wall; rose $\bE$, violet $4\pi\bP/3$, ink $\bE^{\mathrm{loc}}$. Right: grey dashed $\omega_0$, ink $\omega_{\mathrm{TO}}$, amber $\omega_{\mathrm{LO}}$, hollow markers: measured values; shaded: unstable.
- *Try.* (1) Switch models: $\hbar\omega_0$ jumps from 33.3 to 35.5 to 42.1 meV in GaAs while the readouts of the observables do not change. (2) With the Lorentz model raise $Z^*$ until $\omega_{\mathrm{TO}}$ reaches zero: the soft mode of a ferroelectric, with $\varepsilon_{\mathrm s}$ diverging. (3) Choose NaCl ($\varepsilon_\infty=2.34<3$): $\omega_{\mathrm{LO}}$ stays above $\omega_0$ for every model.
- *Model.* <a class="eqref" href="#b3-eq-loc"></a>, <a class="eqref" href="#b3-eq-zeta"></a> with $\omega_{\mathrm p}^{\mathrm{ion}}$ from <a class="eqref" href="#b3-eq-pion"></a> and the LST relation <a class="eqref" href="#b3-eq-lst"></a>; presets in the table above (literature values, illustrative); Clausius–Mossotti with Szigeti charges is not in the sources.

**Notes for the author.** Start from v1's `viz-b-dielectric-lab` script: its `model()` microscopic branch with `sOf()` is exactly $\zeta$ ($s$ there), and its level diagram (`pV`) and warnings are the right-hand panel. Rename every symbol ($\omega_{\rm O}\to\omega_0$, $\tilde\omega\to\omega$, $\Omega_{\rm uc}\to\OV$, $s\to\zeta$, $\omega_p\to\omega_{\mathrm p}^{\mathrm{ion}}$) and switch the axes from cm⁻¹ to meV.

---

## B3.4 · `fig-lab.html` · `viz-b3-lab` · port of `viz-b-dielectric-lab`

**Title.** "Reststrahlen lab: $\varepsilon(\omega)$, reflectance, energy loss and $W$" `<span class="tag">interactive</span>`

**What the reader must learn.** One dielectric function gives four observables. $\mathrm{Im}\,\varepsilon$ peaks at $\omega_{\mathrm{TO}}$ (infrared absorption); $R\approx1$ in the Reststrahlen band between $\omega_{\mathrm{TO}}$ and $\omega_{\mathrm{LO}}$, where $\varepsilon<0$; the loss function $-\mathrm{Im}(1/\varepsilon)$ peaks at $\omega_{\mathrm{LO}}$ (EELS); $W/v=1/\varepsilon$ splits into the flat electronic $1/\varepsilon_\infty$ and a lattice part with its pole at $\omega_{\mathrm{LO}}$, going from $1/\varepsilon_{\mathrm s}$ at $\omega=0$ to $1/\varepsilon_\infty$ far above. The LST relation fixes $\varepsilon_{\mathrm s}$.

**Include position.** Subsection "The Reststrahlen band" (`#b3-rest`), after the first paragraph, which ends: "For NaCl the band runs from 20.3 to 32.7 meV, wavelengths of 61 to 38 µm; for GaAs from 33.3 to 36.2 meV, or 37.2 to 34.2 µm. Figure B3.4 computes the spectra from (B3.11)." The reflectance formula (B3.15) follows; the subsection on $W$ refers to "panel (d)".

**Physics** (`b3-eq-fact` (B3.11), `b3-eq-R` (B3.15), `b3-eq-loss` (B3.16), `b3-eq-W` (B3.17), `b3-eq-lst` (B3.12), `b3-eq-lo` (B3.10)).
- $\varepsilon(\omega)=\varepsilon_\infty(\omega_{\mathrm{LO}}^2-\omega^2-\ii\gamma\omega)/(\omega_{\mathrm{TO}}^2-\omega^2-\ii\gamma\omega)$ (retarded, $\ee^{-\ii\omega t}$, $\mathrm{Im}\,\varepsilon>0$ is absorption).
- $R=|(\sqrt\varepsilon-1)/(\sqrt\varepsilon+1)|^2$ with the branch $\mathrm{Im}\sqrt\varepsilon\ge0$ (`EPP.C.sqrt` returns the principal root, which has $\mathrm{Re}\ge0$; for $\mathrm{Im}\,\varepsilon\ge0$ this also has $\mathrm{Im}\ge0$, so it is the right branch here).
- Loss $-\mathrm{Im}(1/\varepsilon)$.
- $W/v=1/\varepsilon$; $W^{\mathrm e}/v=1/\varepsilon_\infty$; $W^{\mathrm{ph}}/v=1/\varepsilon-1/\varepsilon_\infty=\frac{1}{\varepsilon_\infty}\frac{\omega_{\mathrm{LO}}^2-\omega_{\mathrm{TO}}^2}{\omega^2-\omega_{\mathrm{LO}}^2+\ii\gamma\omega}$ (plot real parts).
- $\varepsilon_{\mathrm s}=\varepsilon_\infty\omega_{\mathrm{LO}}^2/\omega_{\mathrm{TO}}^2$; $\hbar\omega_{\mathrm p}^{\mathrm{ion}}=[\varepsilon_\infty(\omega_{\mathrm{LO}}^2-\omega_{\mathrm{TO}}^2)]^{1/2}$; implied $Z^*=[(\hbar\omega_{\mathrm p}^{\mathrm{ion}})^2\OV M_{\mathrm r}/7.564\times10^5\ \text{meV}^2]^{1/2}$ with the preset's $\OV$, $M_{\mathrm r}$.
- Materials: the table at the top; hBN in plane is uniaxial and is treated as isotropic (say so in the Model line).

**What to draw.** Four `EPP.Plot` panels in a 2×2 grid (one column below 700 px), all with $x=\hbar\omega$ in meV, fixed limits per preset ($[0,\ 1.6\,\hbar\omega_{\mathrm{LO}}]$ rounded up: GaAs 0–60, NaCl 0–55, 3C-SiC 0–200, hBN 0–320 meV), hover giving meV and cm⁻¹:
- (a) `ε(ω)`: Re ε (ink, solid), Im ε (rose, dashed); dotted horizontal line at $\varepsilon_\infty$ labelled `ε_∞`; dot at $\omega=0$ labelled `ε_s`. Fixed $y$ range per preset (about $[-2.5,3.5]\times\varepsilon_{\mathrm s}$), clipping the resonance peaks.
- (b) `reflectance R`: $R$ (rose), $y\in[0,1.04]$; the band shaded with legend "Reststrahlen band (ε < 0)".
- (c) `−Im 1/ε`: loss function (violet), label "−Im 1/ε (EELS)".
- (d) `W/v = 1/ε`: Re $W/v$ total (rose), $W^{\mathrm e}/v=1/\varepsilon_\infty$ (teal, dashed, flat, label `W^{e}/v = 1/ε_∞`), Re $W^{\mathrm{ph}}/v$ (amber, label `W^{ph}/v`); dotted line at $1/\varepsilon_{\mathrm s}$ labelled `1/ε_s`.
- In every panel: the band $\omega_{\mathrm{TO}}<\omega<\omega_{\mathrm{LO}}$ shaded amber (alpha about 0.1) and vertical lines `TO`, `LO`. Keep v1's dense sampling near $\omega_{\mathrm{TO}}$ and $\omega_{\mathrm{LO}}$ so the narrow peaks are resolved.

**Controls (5).**
1. Preset: GaAs (default), NaCl, 3C-SiC, hBN (in plane). Loading a preset sets the other four controls and the axis limits.
2. $\hbar\omega_{\mathrm{TO}}$: 0.5–1.6 × preset value (meV), readout also in cm⁻¹.
3. $\hbar\omega_{\mathrm{LO}}$: from 1.002 $\omega_{\mathrm{TO}}$ to 1.8 × preset value (clamp just above $\omega_{\mathrm{TO}}$, as v1 does, with a note that $\omega_{\mathrm{LO}}\le\omega_{\mathrm{TO}}$ would need $\varepsilon_{\mathrm s}<\varepsilon_\infty$).
4. $\varepsilon_\infty$: 1–15, step 0.01.
5. Damping $\hbar\gamma$: 0.05–3 meV.

**Readouts.** $\hbar\omega_{\mathrm{TO}}$ and $\hbar\omega_{\mathrm{LO}}$ (meV · cm⁻¹); $\varepsilon_{\mathrm s}$ (LST); $\hbar\omega_{\mathrm p}^{\mathrm{ion}}$ from (B3.10); implied $Z^*$ (with the preset's $\OV$, $M_{\mathrm r}$); $1/\varepsilon_\infty-1/\varepsilon_{\mathrm s}$; the band in µm ($\lambda[\mu\text{m}]=1239.84/\hbar\omega[\text{meV}]$ — use the word "wavelength", not the symbol λ, which v2 reserves for photon polarization).

**Initial state.** GaAs: 33.3 / 36.2 meV, $\varepsilon_\infty=10.9$, $\hbar\gamma=0.31$ meV (so $\varepsilon_{\mathrm s}=12.88$ and the band spans 34.2–37.2 µm).

**Caption key points.**
- *What to notice.* Infrared absorption ($\mathrm{Im}\,\varepsilon$) peaks at TO, while the loss function and the lattice part of $W$ peak at LO; between them $\varepsilon<0$ and the crystal reflects almost everything. $W$ starts at $v/\varepsilon_{\mathrm s}$ and ends at $v/\varepsilon_\infty$.
- *Marks.* Shaded: Reststrahlen band; ink Re ε, rose dashed Im ε; rose $R$; violet loss; teal dashed $W^{\mathrm e}/v$; amber $W^{\mathrm{ph}}/v$; rose $W/v$; dotted lines $\varepsilon_\infty$, $1/\varepsilon_{\mathrm s}$.
- *Try.* (1) Choose NaCl: the small $\varepsilon_\infty$ gives a band wider than half the TO frequency. (2) Drag $\omega_{\mathrm{LO}}$ down to $\omega_{\mathrm{TO}}$, a nonpolar crystal: the band, the loss peak and $W^{\mathrm{ph}}$ disappear. (3) Raise $\gamma$: the plateau drops below 1 and the peaks broaden, but absorption stays at TO and loss at LO. Also note the reflectance dip just above LO, where $\varepsilon=1$ (36.5 meV in GaAs).
- *Model.* <a class="eqref" href="#b3-eq-fact"></a> with one damping $\gamma$; <a class="eqref" href="#b3-eq-R"></a>, <a class="eqref" href="#b3-eq-loss"></a>, <a class="eqref" href="#b3-eq-W"></a>, <a class="eqref" href="#b3-eq-lst"></a>; presets from the table (literature values, illustrative; hBN uniaxial, in-plane values used as if isotropic).

**What must change from v1.**
- v1 has nine controls (parameter mode, preset, local-field model, $\varepsilon_\infty$, $\gamma$, plus $\omega_{\mathrm{TO}}$, $\omega_{\mathrm{LO}}$ or $Z^*$, $\Omega_{\rm uc}$, $M$, $\omega_{\rm O}$), a fifth plot (the level diagram) and a warnings callout. Keep only the spectroscopic mode with the five controls above; the microscopic mode, the local-field selector, the level diagram and the warnings move to `fig-local.html`.
- Source notation everywhere in v1 strings and labels: `ε(0)`→`ε_s`, `ω̃_{TO}`→`ω_{TO}`, `ω_p`→`ω_p^{ion}`, `Ω_{uc}`→`Ω_V`, `W_e/v`→`W^{e}/v`, `Re W_{ph}/v`→`W^{ph}/v`, `Z^* (|e|)`→`Z^*` (dimensionless), `reduced mass M (u)`→`M_r (amu)`, `Eq. (E.160)`, `§B.5–B.7` and other source/v1 references removed.
- Axis unit: v1 plots cm⁻¹; v2 plots meV with cm⁻¹ in hover and readouts.
- Default preset: NaCl in v1, GaAs in v2 (the section's numbers are GaAs).
- v1's caption puts Model first and has no Marks/Try; rewrite in the v2 order. Rename ids and classes (`b-dielectric-p-eps` etc. → `b3-lab-…`).
- Check v1's constants against the table: its `KWP = 4.92064e7` (cm⁻² units) is consistent; in meV use $7.564\times10^5$ meV² $\times Z^{*2}/(\OV M_{\mathrm r})$.
