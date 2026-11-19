# X.1 figure specs (section `x-light-cone`, prefix `x1`)

Three figures, in reading order: `polwave` (new), `dyson` (port), `cone` (port).
Notation is v2 (NOTATION.md): Gaussian units, bold vectors, upright labels (`\mathrm{L}`, `\mathrm{T}`, `\mathrm{e}`,
`\mathrm{ion}`, `\mathrm{light}`), no source symbols anywhere in the figure (labels, readouts, aria-labels, strings).
Equation references in captions use `<a class="eqref" href="#x1-eq-…"></a>`; the ids are listed per figure.
Constants: ħc = 1973.27 eV Å, ħ = 6.582120e-16 eV s. All three figures follow the quality bar of
`v2/build/BUILDER_GUIDE.md` §5 (390 px, keyboard, aria-labels, `EPP.color` at draw time, `EPP.anim` for motion).

Symbols the figures use (all defined in the section): $\bq$, $q$, $\omega$, $\hbar\omega$, $c$, $\bP$, $\bP^{\mathrm L}$,
$\bP^{\mathrm T}$, $\bE^{\mathrm L}$, $\bE^{\mathrm T}$, $\varepsilon_\infty$, $q_{\mathrm{light}}=\sqrt{\varepsilon_\infty}\,\omega/c$,
$D$, $D^0$, $D^{0\mathrm L}$, $D^{0\mathrm T}$, $D^{\mathrm L}$, $D^{\mathrm T}$, $D^{\mathrm e}$, $\Pi$, $\Pi^{\mathrm e}$, $\Pi^{\mathrm{ion}}$,
$v$, $\chi^{\mathrm{irr}}$, $W$, $W^{\mathrm e}$, $W^{\mathrm{ph}}$, $F$, $\omega_{\mathrm{LO}}$, $\omega_{\mathrm{TO}}$.
Do not use $\lambda$ for anything (photon polarization index), nor $\sigma$ (conductivity), nor $x$ for $cq/\omega$ (Cartesian).
The speed-of-light multiplier of the v1 figures, called $\lambda_c$ there, has **no symbol** in v2: call it
"speed of light: c ×" in labels and "the factor by which c is multiplied" in captions.

---

## 1. `fig-polwave.html` · id `viz-x1-polwave` · new

**Kind.** New figure (no v1 counterpart). Title: "The field of a polarization wave" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. Longitudinal and transverse are geometric: $\bP$ along $\hat\bq$ piles up bound charge in sheets, $\bP$ across $\hat\bq$ does not.
2. The longitudinal field is $-4\pi\bP^{\mathrm L}$ (divided by $\varepsilon_\infty$ with electrons) at **every** $q$: the instantaneous
   Coulomb field of the charge sheets, independent of $c$.
3. The transverse field depends on $q$ only through $q/q_{\mathrm{light}}$: equal to the longitudinal field deep inside the light
   cone, resonant at $q=q_{\mathrm{light}}$, and outside the cone tiny ($\propto(q_{\mathrm{light}}/q)^2$) and **parallel** to $\bP$.

**Include position.** In h4 `x1-polwave`, right after the paragraph that ends
"They divide the longitudinal field by the electronic dielectric constant $\varepsilon_\infty$ and slow light to $c/\sqrt{\varepsilon_\infty}$:
the transverse resonance moves to $cq=\sqrt{\varepsilon_\infty}\,\omega$, and deep inside the cone the transverse field is divided by
$\varepsilon_\infty$ too (shown in (X1.13) below). Figure X1.1 draws both waves and their fields as $q$ crosses the resonance." and before h4 `x1-free` ("The free photon propagator").
The sentences just before it are the three-regime list (pattern slower than light / at the speed of light / faster than light)
and the paragraph on the electrons.

**Physics.**
- Polarization wave $\bP\,\ee^{\ii(\bq\cdot\br-\omega t)}$ with $\bq=q\hat{\mathbf x}$; current $\bj=-\ii\omega\bP$; bound charge
  $\rho=-\ii\bq\cdot\bP$ (nonzero only for $\bP^{\mathrm L}$).
- Fields, ions only (vacuum), eq. (X1.4) `x1-eq-EP`:
  $\bE^{\mathrm L}=-4\pi\bP^{\mathrm L}$, $\bE^{\mathrm T}=-4\pi\bP^{\mathrm T}/(1-c^2q^2/\omega^2)$.
- Fields with electrons, from (X1.13) `x1-eq-iso` with $\varepsilon^{\mathrm L}=\varepsilon^{\mathrm T}=\varepsilon_\infty$ and
  $\bE=(\omega^2/c)D\bP$: $\bE^{\mathrm L}=-4\pi\bP^{\mathrm L}/\varepsilon_\infty$, $\bE^{\mathrm T}=-4\pi\bP^{\mathrm T}/(\varepsilon_\infty-c^2q^2/\omega^2)$.
- Use $r_{\mathrm L}=E^{\mathrm L}/(-4\pi P^{\mathrm L})=1/\varepsilon$ and $r_{\mathrm T}=E^{\mathrm T}/(-4\pi P^{\mathrm T})=1/(\varepsilon-(cq/\omega)^2)$,
  with $\varepsilon=1$ (ions only) or $\varepsilon_\infty$. Add a small damping only for display near the resonance:
  $1/(\varepsilon-(cq/\omega)^2+\ii\,0.02\,\varepsilon)$; arrows use the real part (the physical field is
  $\mathrm{Re}[r\,\cdot(-4\pi\bP)\,\ee^{\ii(qx-\omega t)}]$, so the imaginary part shows up as a phase shift; it is negligible
  except within a few percent of $q_{\mathrm{light}}$).
- Materials (literature values, illustrative): GaAs TO $\hbar\omega=33.3$ meV, $\varepsilon_\infty=10.9$
  ($\omega/c=1.69\times10^{-5}$ Å⁻¹, $q_{\mathrm{light}}=5.6\times10^{-5}$ Å⁻¹); hBN in-plane TO $\hbar\omega=169.9$ meV,
  $\varepsilon_\infty=4.87$ ($\omega/c=8.6\times10^{-5}$, $q_{\mathrm{light}}=1.9\times10^{-4}$ Å⁻¹).
- Check number for the caption: GaAs, ions only, $q=0.025$ Å⁻¹: $|r_{\mathrm T}|=(\omega/cq)^2/(1-(\omega/cq)^2)\approx4.6\times10^{-7}$.

**What to draw.**
- Two horizontal strips, one per row, sharing the $x$ axis, each showing **one wavelength** $2\pi/q$ (rescaled to the width):
  - row "longitudinal: $\bP$ along $\bq$" (LO): 13–17 sites with $\bP$ arrows along $x$ (amber, `phonon`), length
    $\propto\cos(qx-\omega t)$; between them, shaded bands for the bound charge $\rho\propto-\partial_xP_x$ (`c5` for +, `electron`
    for −, alpha 0.18) with small "+" / "−" glyphs; $\bE$ arrows (rose, `photon`) at the same sites, slightly offset vertically,
    length $=|r|\times$ (the $\bP$ arrow length), direction from the sign of $\mathrm{Re}\,r$ (antiparallel to $\bP$ for
    $\mathrm{Re}\,r>0$).
  - row "transverse: $\bP$ across $\bq$" (TO): the same with arrows along $y$; no charge bands; a muted note
    "no bound charge: $\nabla\cdot\bP=0$".
  - an arrow for $\bq$ (bold label **q**, `drawRich` with `{bold:true}`) at the left of each strip, pointing along $+x$.
  - a scale bar under the strips: "one wavelength $=2\pi/q=$ …" with units chosen automatically (nm, μm, mm).
  - cap the drawn $\bE$ length at 2.5× the $\bP$ length; when capped, draw the arrowhead hollow and print "resonance" next to it.
  - when $|r_{\mathrm T}|<0.02$ the transverse arrows vanish; print beside the row "$|\bE^{\mathrm T}|/4\pi|\bP|=$ 4.6×10⁻⁷ (too small to
    draw)" in rich text.
- Below the strips, a small plot (EPP.Plot, aspect ≈ 0.45): $|E/4\pi P|$ against $q$ on log–log axes, $q\in[10^{-6},10^{-1}]$
  Å⁻¹, $y\in[10^{-8},10^{2}]$, fixed limits. Lines: longitudinal $|r_{\mathrm L}|$ (dashed `photon`, flat), transverse $|r_{\mathrm T}|$
  (solid `photon`, peak at $q_{\mathrm{light}}$, slope −2 outside); a vertical line at $q_{\mathrm{light}}$ labelled `q_{light}`; a shaded
  band $q<q_{\mathrm{light}}$ labelled "light cone" (`photon`, alpha 0.08); a dot at the current $q$ on each line.
  Axis labels: `q (Å^{-1})`, `|E| / 4π|P|`.
- Optional animation (EPP.anim, Play/Pause): the phase $\omega t$ advances so that both patterns move to the right at the
  phase velocity, one period every 3 s regardless of $\omega$. Reduced motion: static frame at $\omega t=0$.

**Controls (4).**
1. wavevector $q$: log range, $10^{-6}$ to $10^{-1}$ Å⁻¹, `aria-valuetext` like "5.0 × 10⁻⁵ per ångström"; clicking in the
   plot sets $q$.
2. material: seg "GaAs TO (33.3 meV)" / "hBN TO (169.9 meV)" (sets $\hbar\omega$ and $\varepsilon_\infty$).
3. electrons: seg "ions only (vacuum)" / "with electrons ($\varepsilon_\infty$)".
4. Play/Pause (from EPP.anim).

**Readouts** (MathJax in HTML): $q$; $q_{\mathrm{light}}$ (with $\varepsilon=1$ or $\varepsilon_\infty$ as selected); $q/q_{\mathrm{light}}$;
wavelength $2\pi/q$; phase velocity $\omega/q$ in units of the speed of light in the medium, $c/\sqrt\varepsilon$;
$E^{\mathrm L}/(-4\pi P^{\mathrm L})$; $E^{\mathrm T}/(-4\pi P^{\mathrm T})$ with its sign and the words "antiparallel to $\bP$" or
"parallel to $\bP$".

**Initial state.** GaAs, ions only, $q=5.0\times10^{-5}$ Å⁻¹ (three times $\omega/c$): LO field $=4\pi P$ antiparallel, TO field
$1/8$ of it and **parallel** to $\bP$, both visible; animation paused.

**Caption key points.**
- *What to notice.* The longitudinal field stays at $-4\pi\bP$ (divided by $\varepsilon_\infty$ with electrons) whatever $q$: it is
  the Coulomb field of the charge sheets, set up instantly. The transverse field equals it only deep inside the light cone,
  diverges at $q=q_{\mathrm{light}}$, and outside the cone shrinks as $(q_{\mathrm{light}}/q)^2$ and turns parallel to $\bP$.
- *Marks.* Amber arrows: polarization $\bP$; rose arrows: field $\bE$, length in units of $4\pi|\bP|$, hollow head when capped at
  the resonance; shaded bands with + and −: bound charge $-\nabla\cdot\bP$; scale bar: the wavelength. Plot: dashed longitudinal,
  solid transverse, shaded light cone, dots at the current $q$.
- *Try.* Sweep $q$ from $10^{-6}$ to $10^{-1}$ Å⁻¹ and watch only the TO row change. Switch the electrons on: both fields shrink
  by $1/\varepsilon_\infty$ inside the cone, and the resonance moves from $\omega/c$ to $\sqrt{\varepsilon_\infty}\,\omega/c$. Set
  $q=0.025$ Å⁻¹, the scale of electron–phonon scattering in GaAs, and read the transverse field: about $5\times10^{-7}$ of the
  longitudinal one.
- *Model.* Maxwell's equations in the temporal gauge (X1.2) `x1-eq-wave` with the polarization current $\bj=-\ii\omega\bP$; fields
  (X1.4) `x1-eq-EP`, and with electrons (X1.13) `x1-eq-iso` through $\bE=(\omega^2/c)D\bP$; a damping of 2% of $\varepsilon$ added
  for display; GaAs and hBN values as above (literature values, illustrative).

---

## 2. `fig-dyson.html` · id `viz-x1-dyson` · port

**Kind.** Port of v1 `viz-x-light-cone-dyson`, in `doc/epp-primer/sections/X01-light-cone.html` (the `<figure>` at
lines 54–75 and the second IIFE of the `<script>`, "Figure 2: Giustino's W chain vs EPP's D chain"). Title: "Two Dyson
chains, one structure" `<span class="tag">interactive</span>`.

**What the reader must learn.** The screened Coulomb interaction of Part A and the photon propagator of Part B are built by
the same chain: bare line, bubble, line, bubble…, and the lattice enters both as a phonon line between two vertices. Taking
$c\to\infty$ removes every transverse (wavy) line, and the Part B column becomes the Part A column glyph for glyph.

**Include position.** In h3 `x1-dyson`, right after the paragraph that ends "…Iterating either equation gives the chain of its
diagrams: the bare line, then the line with one insertion of the response ($\chi^{\mathrm{irr}}$ or $\Pi$, drawn as a bubble), then
two insertions, and so on. Figure X1.2 draws the two chains side by side." (which follows (X1.11) `x1-eq-dysonW` and (X1.12) `x1-eq-dysonD`), and before "The dictionary between the two columns:"
and its table.

**Physics (rows; left = Part A, right = Part B).**
1. *The bare line.* Left: $v$, the instantaneous Coulomb line, $v(\br-\br')\delta(t-t')$, (X1.1) `x1-eq-v`.
   Right: $D^0=D^{0\mathrm L}+D^{0\mathrm T}$, (X1.5) `x1-eq-D0`: a dashed line ($D^{0\mathrm L}$, the Coulomb line rewritten for
   currents, (X1.6) `x1-eq-D0L`) plus a wavy line ($D^{0\mathrm T}$, free photons).
2. *Screening by the electrons.* Left: $W^{\mathrm e}=v+v\chi^{\mathrm{irr}}W^{\mathrm e}$ (nuclei clamped; (X1.11) `x1-eq-dysonW`
   with the electrons only). Right: $D^{\mathrm e}=D^0+D^0\Pi^{\mathrm e}D^{\mathrm e}$, first line of (X1.15) `x1-eq-twostep`.
3. *Adding the lattice.* Left: $W=W^{\mathrm e}+W^{\mathrm{ph}}$ with $W^{\mathrm{ph}}$ = screened vertex × $F$ × screened vertex
   (`a5-eq-split`, `a5-eq-Wph`): drawn as $W^{\mathrm e}$ – vertex – $F$ (phonon zigzag) – vertex – $W^{\mathrm e}$; the vertex is the
   charge vertex of a displaced nucleus, and $W^{\mathrm e}$ ending on it is the screened vertex.
   Right: $D=D^{\mathrm e}+D^{\mathrm e}\Pi^{\mathrm{ion}}D$, second line of (X1.15) `x1-eq-twostep`, drawn with the ionic bubbles
   summed: $D^{\mathrm e}$ – vertex – $F$ – vertex – $D^{\mathrm e}$, where the vertex is the Born-charge current vertex $J$ of
   Section B.7 (`gl-b7-J`) and the summed ionic part is $-JFJ^\dagger$ (the minus sign is not drawn; say so in the caption).
- $c\to\infty$: $D^{0\mathrm T}\to0$, so every wavy line disappears; $D^0\to D^{0\mathrm L}$ (dashed); the dressed lines keep only
  their longitudinal part, $D^{\mathrm e}\to-\frac{c}{\omega^2}q_\alpha q_\beta W^{\mathrm e}/e^2$ and likewise $D\to\propto W$, (X1.14)
  `x1-eq-DLW`; the electronic bubble's longitudinal part is $e^2\chi^{\mathrm{irr}}$, (X1.10) `x1-eq-conn`.

**What to draw.** Keep the v1 layout and the `EPP.fd.draw` glyphs, renamed:
- two columns (`viz-grid two`) with headers "Part A · Coulomb line" and "Part B · photon line"; three rows with headings
  "1 · The bare line", "2 · Screening by the electrons", "3 · Adding the lattice". Under each row heading, a one-line source note
  per column in v2 terms, with equation links: row 1 "(X1.1)" / "(X1.5)"; row 2 "(X1.11), nuclei clamped" / "(X1.15), first
  line"; row 3 "Section A.5" / "(X1.15), second line" (use `a.eqref` in the HTML).
- glyphs: `coulomb` (dashed, `photon` colour) for $v$ and $D^{0\mathrm L}$; `screened` (double dashed) for $W^{\mathrm e}$, $W$ and for
  the longitudinal parts of $D^{\mathrm e}$ and $D$; `photon` (wavy) for $D^{0\mathrm T}$; `photon2` (heavy wavy) for the dressed
  transverse parts; a wavy line drawn on top of a dashed one is the full field with both polarizations ($D^0$, $D^{\mathrm e}$, $D$);
  `blob` in `electron` colour for the bubbles $\chi^{\mathrm{irr}}$ (left) and $\Pi^{\mathrm e}$ (right); `phonon2` zigzag for $F$ in both
  columns; `vertex` dots: plain `ink` on the left (charge vertex of a displaced nucleus), `c4` teal on the right with labels
  $J$ and $J^\dagger$ (Born-charge current vertex).
- labels (rich text): `v`, `W^e`, `W`, `χ^{irr}`, `F`, `D^0`, `D^{0L}`, `D^{0T}`, `D^e`, `D`, `Π^e`, `J`, `J^†` with upright label
  letters (e, irr, L, T are labels).
- at $c\to\infty$ the right column shows the same chains without wavy lines; keep the labels $D^0$, $D^{\mathrm e}$, $D$ and add one
  muted line in row 1 "$D^{0\mathrm T}\to0$" and, under rows 2–3, "only the longitudinal part survives: $D^{\mathrm e}\propto W^{\mathrm e}$,
  $D\propto W$".

**Controls (2).**
1. seg "speed of light": "finite c (Part B)" / "c → ∞ (electrostatic)".
2. seg "highlight": "none" / "bare line" / "bubble" / "vertex" / "phonon line": thickens the matching glyphs in **both**
   columns (and dims the rest to 35% opacity), and writes the correspondence in the readout:
   - bare line: "$v$ ↔ $D^0=D^{0\mathrm L}+D^{0\mathrm T}$; $D^{0\mathrm L}$ is $v$ rewritten for currents, (X1.6)";
   - bubble: "$\chi^{\mathrm{irr}}$ (density response) ↔ $\Pi^{\mathrm e}$ (current response); its longitudinal part is $e^2\chi^{\mathrm{irr}}$, (X1.10)";
   - vertex: "charge vertex of a displaced nucleus ↔ Born-charge current vertex $J$ (Section B.7)";
   - phonon line: "$F$ in both: the lattice enters as a phonon between two vertices".

**Readouts.** "transverse (wavy) lines in the Part B column: 10" (0 at $c\to\infty$; recount after the port); the highlight line
above; a short status line: "finite $c$: every Part B line carries both polarizations" / "$c\to\infty$: the Part A chain, with
$v\to D^{0\mathrm L}$ and densities → currents".

**Initial state.** finite $c$, highlight "none".

**Caption key points.**
- *What to notice.* The two columns have the same structure row by row. Switch to $c\to\infty$: every wavy line disappears,
  because $D^{0\mathrm T}$ vanishes, and the Part B column becomes the Part A column glyph for glyph. Coulomb was the longitudinal
  photon all along.
- *Marks.* Dashed: Coulomb line ($v$, $D^{0\mathrm L}$); double dashed: screened Coulomb line; wavy: transverse photon (heavy when
  dressed); wavy on dashed: the full field, both polarizations. Blue disks: bubbles $\chi^{\mathrm{irr}}$ and $\Pi^{\mathrm e}$; amber
  zigzag: phonon propagator $F$; plain dots: charge vertex of a displaced nucleus, so that $W^{\mathrm e}$ ending on it is the
  screened vertex of Section A.5; teal dots: Born-charge current vertex $J$; the minus sign of the summed ionic part
  $-JFJ^\dagger$ is not drawn.
- *Try.* Toggle the speed of light; use "highlight" to find each object's counterpart; compare row 3 with (X1.15).
- *Model.* (X1.1), (X1.5), (X1.11), (X1.12), (X1.15) and (X1.14) for the $c\to\infty$ view; row 3 left after Section A.5
  (`a5-eq-split`, `a5-eq-Wph`), row 3 right after Section B.7.

**What in v1 is wrong or in source notation and must change.**
- Column headers "Giustino · …", "EPP · …" and the source equation numbers in every `.src` line (G Eqs. 74, 94, 104–105, 110;
  E Eqs. 375, 388, 379, 387–393) → v2 headers and v2 equation links as above.
- $P_{\rm e}$ → $\chi^{\mathrm{irr}}$; $W_{\rm e}$ → $W^{\mathrm e}$; the nuclear density correlator $\mathcal D$ (Giustino's $D$) → the
  v2 form of row 3 with $F$ between two vertices; $D^{0l}$, $D^{0t}$ → $D^{0\mathrm L}$, $D^{0\mathrm T}$; $\Pi^{\rm e,0}$ → $\Pi^{\mathrm e}$
  (in v2 $\Pi^{\mathrm e}$ is the irreducible electronic part; no superscript 0); $D^{e,ll}$, $D^{ll}$ at $c\to\infty$ → keep $D^{\mathrm e}$,
  $D$ and add the "only the longitudinal part survives" line; $J^*$ → $J^\dagger$.
- v1 caption: "Glyphs" paragraph first, then "What to notice" → v2 structure (`cap-key` What to notice, then Marks, Try, Model);
  remove all source symbols ($Z_\kappa\nabla\delta(\br-\boldsymbol\tau^0_{\kappa p})$, $\epsilon_{\rm e}^{-1}\nabla V_\kappa$, E Eq. 441, …).
- The note string "EPP column at c → ∞: Giustino's chain, with v → D⁰ˡ and ρ → j." and the aria-labels ("W_e", "P_e",
  "Pi e0", "D0l") → v2 wording ("Part B column", "W^e", "chi irr", "D zero L").
- Class names `x-lc-*` → `x1-dyson-*`.

---

## 3. `fig-cone.html` · id `viz-x1-cone` · port

**Kind.** Port of v1 `viz-x-light-cone-cone`, in `doc/epp-primer/sections/X01-light-cone.html` (the `<figure>` at lines
91–98 and the first IIFE of the `<script>`, "Figure 1: how big is the light cone?"). Title: "The light cone, the zone and the
phonon branches, to scale" `<span class="tag">interactive</span>`.

**What the reader must learn.** At phonon frequencies the light cone is a sphere of radius ~$10^{-4}$ Å⁻¹, four decades below
the zone edge (~1 Å⁻¹): the phonon branches of Part A cross the light line only in a sliver next to Γ, where the TO branch
turns into polaritons. Outside the cone the transverse field falls as $(q_{\mathrm{light}}/q)^2$, while the longitudinal (Coulomb)
line is flat and does not know $c$. Raising $c$ shrinks the cone to the point $q=0$: Part A's world.

**Include position.** In h3 `x1-cone`, right after the paragraph that ends "…Inside it the transverse and longitudinal fields
are comparable. Outside it the transverse field falls off as $(q_{\mathrm{light}}/q)^2$. Figure X1.3 draws the cone, the Brillouin
zone and the phonon branches of GaAs to scale on a logarithmic axis." (which follows (X1.16) `x1-eq-ratio`), and before
"The numbers are what make Part A work."

**Physics.** Let $s$ be the speed-of-light factor (the figure's control; no symbol in labels), $c\to s\,c$, and
$\varepsilon=\varepsilon_\infty$ (crystal) or 1 (vacuum).
- Light lines: $\hbar\omega=\hbar s c\,q/\sqrt{\varepsilon}$ (medium) and $\hbar\omega=\hbar s c\,q$ (vacuum).
- Light-cone radius at the selected $\hbar\omega$: $q_{\mathrm{light}}=\sqrt\varepsilon\,\omega/(s c)$, (X1.16) `x1-eq-ratio`.
- Field per unit current, normalized: $|\omega^2D^{\mathrm L}/4\pi c|=1/\varepsilon$ (flat) and
  $|\omega^2D^{\mathrm T}/4\pi c|=1/|\varepsilon-(scq/\omega)^2+\ii\,0.01\,\varepsilon|$, from (X1.13) `x1-eq-iso` with
  $\varepsilon^{\mathrm L}=\varepsilon^{\mathrm T}=\varepsilon$; their ratio is (X1.16). (Normalization note for the caption: this is also
  $|\bE|/4\pi|\bP|$ for a polarization wave, as in Figure X1.1.)
- Polaritons (transverse modes of the crystal, isotropic Lorentz model, no damping): solve $\omega^2\varepsilon(\omega)=s^2c^2q^2$
  with $\varepsilon(\omega)=\varepsilon_\infty(\omega_{\mathrm{LO}}^2-\omega^2)/(\omega_{\mathrm{TO}}^2-\omega^2)$:
  $\omega_\pm^2=\big[b\pm\sqrt{b^2-4\varepsilon_\infty s^2c^2q^2\omega_{\mathrm{TO}}^2}\big]/(2\varepsilon_\infty)$ with
  $b=\varepsilon_\infty\omega_{\mathrm{LO}}^2+s^2c^2q^2$ (closed form). Lower branch → $\omega_{\mathrm{TO}}$ and upper → light line for
  $q\gg q_{\mathrm{light}}$; lower → 0 and upper → $\omega_{\mathrm{LO}}$ for $q\to0$. The LO branch is not coupled (flat).
  Draw only near the anticrossing, $q<30\,q_{\mathrm{light}}(\omega_{\mathrm{TO}})$; beyond that they coincide with TO and the light line.
- Phonon branches across the zone (illustrative models, along Γ–X for GaAs, Γ–K for hBN; $q_{\mathrm{edge}}=2\pi/a$ for GaAs,
  $4\pi/3a$ for hBN):
  - optical: $\hbar\omega(q)=\hbar\omega_\Gamma+(\hbar\omega_{\mathrm{edge}}-\hbar\omega_\Gamma)\sin^2(\pi q/2q_{\mathrm{edge}})$;
  - acoustic: $\hbar\omega(q)=\hbar\omega_{\mathrm{edge}}\big[(1+B)\sin(\pi t/2)+B\sin(3\pi t/2)\big]$ with $t=q/q_{\mathrm{edge}}$ and
    $B=(2r/\pi-1)/4$, $r=\hbar v\,q_{\mathrm{edge}}/\hbar\omega_{\mathrm{edge}}$: slope $v$ (the sound velocity) at Γ, the value
    $\hbar\omega_{\mathrm{edge}}$ at the zone edge, and zero slope there (as implemented; a $\tanh$ form cannot meet both the slope
    and the edge value);
  - GaAs ($a=5.653$ Å, $q_{\mathrm{edge}}=1.111$ Å⁻¹, $\varepsilon_\infty=10.9$): LO 36.2 → 29.8 meV, TO 33.3 → 31.3 meV,
    LA $v=4.7$ km/s → 27.8 meV, TA $v=3.3$ km/s → 9.9 meV;
  - hBN in plane ($a=2.504$ Å, $q_{\mathrm{edge}}=1.673$ Å⁻¹, $\varepsilon_\infty=4.87$): LO 199.6 → 160 meV, TO 169.9 → 160 meV,
    LA $v=20$ km/s → 135 meV, TA $v=13$ km/s → 115 meV (illustrative; treated as isotropic).
  - $\hbar v$ in meV Å: $v[\mathrm{km/s}]\times6.582$ (so 4.7 km/s → 30.9 meV Å).
- Characteristic wavevectors (bands in panel b, illustrative): electron–phonon scattering in polar semiconductors
  0.02–0.07 Å⁻¹ (GaAs 0.025 Å⁻¹, where a conduction electron's kinetic energy equals $\hbar\omega_{\mathrm{LO}}$); spacing of
  practical $\bq$ grids 0.01–0.05 Å⁻¹; zone edge from $\pi/a$ with $a=6$ Å (0.52 Å⁻¹) to the GaAs X point (1.11 Å⁻¹).
- Check numbers (GaAs, crystal, $s=1$, 33.3 meV): $q_{\mathrm{light}}=5.57\times10^{-5}$ Å⁻¹, $2\pi/q_{\mathrm{light}}=11.3$ μm,
  $q_{\mathrm{light}}/q_{\mathrm{edge}}=5.0\times10^{-5}$, volume fraction $(q_{\mathrm{light}}/q_{\mathrm D})^3=1.3\times10^{-13}$ with the Debye
  radius $q_{\mathrm D}=1.09$ Å⁻¹, $|D^{\mathrm T}/D^{\mathrm L}|=5.0\times10^{-6}$ at 0.025 Å⁻¹ and $2.5\times10^{-9}$ at the X point.

**What to draw.**
- Panel (a), "frequency against wavevector": log $q$ axis $[10^{-7},3]$ Å⁻¹ (ticks every decade, every second decade below
  520 px), **linear** $\hbar\omega$ axis, fixed per material: $[0,50]$ meV for GaAs and vacuum, $[0,240]$ meV for hBN. Axis labels
  `q (Å^{-1})`, `ħω (meV)`.
  - light cone: shade the region to the left of the medium light line (above it), `photon`, alpha 0.08, label "light cone";
  - medium light line solid `photon` (label `ħcq/√ε_∞`), vacuum light line thin dashed `photon` (label `ħcq`);
  - LO dashed `phonon`, TO solid `phonon`, LA and TA thinner solid `phonon` (labels `LO`, `TO`, `LA`, `TA` at their right ends);
  - polaritons: two solid curves in `c4` near the crossing, labelled "polaritons";
  - zone-edge band (`muted`, alpha 0.12) from $\pi/a$ ($a=6$ Å) to $q_{\mathrm{edge}}$, label "zone edge";
  - horizontal dotted `muted` line at the selected $\hbar\omega$ and a vertical `photon` line at $q_{\mathrm{light}}$ labelled `q_{light}`;
    if $\hbar\omega$ is above the axis, an arrow at the top edge with "ħω = 2.00 eV, above the plot".
  - On log axes the cone looks large; add one rich-text note in the panel: "linear scale: the cone is 5 × 10⁻⁵ of the zone
    radius" (computed live).
- Panel (b), "field per unit current": same $q$ axis, log $y\in[10^{-12},10^{3}]$, label `|ω^2 D / 4πc|`; longitudinal dashed
  `photon` (flat at $1/\varepsilon$, label "longitudinal = Coulomb: flat, independent of c", short form "Coulomb: flat" below
  560 px); transverse solid `photon` (label "transverse ∝ (q_{light}/q)^2" on its falling side); when $s\neq1$, the transverse curve
  at the real $c$ as a muted dotted reference ("transverse at the real c"); bands for the scattering wavevectors
  (`electron`, alpha 0.10) and the zone edge (`muted`, alpha 0.12); vertical `photon` line at $q_{\mathrm{light}}$; the light-cone band
  as in (a).
- Keep v1's dense sampling around $q_{\mathrm{light}}$ so the resonance peak is resolved, and v1's handling of a cone that falls off
  either edge ("$q_{\mathrm{light}}$ = … Å⁻¹: cone off scale" / "zone inside the cone").
- The v1 ladder strip is dropped; its useful marks are the bands of panel (b).

**Controls (3).**
1. seg "medium": "GaAs" / "hBN" / "vacuum ($\varepsilon=1$)" (sets $\varepsilon$, the branches, the zone edge, the $\hbar\omega$
   axis, and resets $\hbar\omega$ to the TO energy; vacuum keeps GaAs's axis and draws no branches).
2. $\hbar\omega$: log slider, 1 meV to 5 eV (shown as meV below 1 eV).
3. speed of light: log slider, "c ×" from $10^{-2}$ to $10^{6}$, default 1, `aria-valuetext` "speed of light times 1".

**Readouts.** $q_{\mathrm{light}}=\sqrt{\varepsilon}\,\omega/c$ (Å⁻¹); $2\pi/q_{\mathrm{light}}$ (μm); "$q_{\mathrm{light}}$ / zone edge" (the ratio
of $q_{\mathrm{light}}$ to the zone-edge wavevector); "zone volume inside the cone", the cube of the ratio of $q_{\mathrm{light}}$ to the radius
of a sphere with the zone's volume, $(6\pi^2/\OV)^{1/3}$ (GaAs 1.09 Å⁻¹; for hBN use the area fraction of the 2D zone and say so);
$|D^{\mathrm T}/D^{\mathrm L}|$ at 0.025 Å⁻¹; the same at the zone edge. Write these labels in words; do not introduce symbols for the
zone-edge wavevector or the sphere radius (the variable names `qEdge`, `qD` in the code are fine). Readouts in HTML with MathJax labels; plain-number values with Unicode superscripts (as in v1's `sci()`).

**Initial state.** GaAs, $\hbar\omega=33.3$ meV, speed of light ×1: the cone is a wall at $5.6\times10^{-5}$ Å⁻¹, the polaritons
bend there, LO/TO run flat to the zone edge at ~1 Å⁻¹; panel (b) shows the peak at $q_{\mathrm{light}}$ and the $(q_{\mathrm{light}}/q)^2$ tail.

**Caption key points.**
- *What to notice.* At phonon frequencies the light cone is a sliver near $10^{-4}$ Å⁻¹, four decades below the zone edge; the
  TO branch turns into polaritons only there. The longitudinal (Coulomb) line in (b) is flat and independent of $c$; the
  transverse one peaks on the cone and falls as $(q_{\mathrm{light}}/q)^2$: about $5\times10^{-6}$ of Coulomb where electrons scatter
  in GaAs, $10^{-9}$ at the zone edge.
- *Marks.* Shaded: light cone; solid and dashed rose lines: light lines in the medium and in vacuum; amber: phonon branches
  (dashed LO); teal: polaritons; grey band: zone edge; blue band in (b): wavevectors of electron–phonon scattering; vertical line:
  $q_{\mathrm{light}}$; dotted horizontal line: the selected $\hbar\omega$.
- *Try.* Push the speed of light up to $10^6\,c$: the cone and the polaritons leave the left edge, and LO and TO run flat to the
  axis—Part A's world. Pull it down to $10^{-2}c$: even light a hundred times slower keeps the cone inside 1% of the zone.
  Choose "vacuum" and $\hbar\omega=2$ eV: the cone of a visible exciton, at $10^{-3}$ Å⁻¹.
- *Model.* Light lines and (X1.16) `x1-eq-ratio`; propagators from (X1.13) `x1-eq-iso` with $\varepsilon^{\mathrm L}=\varepsilon^{\mathrm T}=\varepsilon_\infty$
  (electrons of an insulator below the gap) and 1% damping; polaritons from $\omega^2\varepsilon(\omega)=c^2q^2$ with the Lorentz
  form of $\varepsilon(\omega)$ (derived in Section B.3, plotted in full in Section B.6); illustrative branch models and literature
  values as listed; the speed-of-light control multiplies $c$ everywhere.

**What in v1 is wrong or in source notation and must change.**
- The speed-of-light factor is called $\lambda_c$ (control label "speed of light factor λ_c (c → λ_c c)", caption, comments):
  $\lambda$ is the photon polarization index in v2 → "speed of light: c ×", no symbol.
- $q_c$ → $q_{\mathrm{light}}$ in labels, readouts and caption; readout keys "q_c = √ε_∞ ω/c", "q_c/q_D", "|D_T/D_L| at …" →
  "$q_{\mathrm{light}}=\sqrt{\varepsilon_\infty}\,\omega/c$", "light cone / zone edge", "$|D^{\mathrm T}/D^{\mathrm L}|$ at …".
- Caption cites "E Eqs. (372) and (374)", "E Eq. 380", "$4\pi c\,\Pi^{\rm e}/\omega^2=\epsilon_\infty-1$" → v2 equation links (X1.13),
  (X1.16) and $\varepsilon$ (not $\epsilon$).
- The figure is not geometric: v1 plots only $|\omega^2D/4\pi c|$ and a ladder of bars. v2 adds panel (a) (frequency against
  wavevector with the light lines, the shaded cone, the zone edge and the phonon branches) as the main panel, and keeps the
  propagator plot as panel (b).
- The ladder rows "$\omega_{\rm ph}/v_F$" and "$1/r_0$ (2D Keldysh length)" use symbols X.1 does not define ($v_{\mathrm F}$; $r_0$ is
  Section B.10's) → drop them; keep the scattering, grid and zone-edge marks as bands in (b) (and $k_{\mathrm{TF}}$ only if defined in
  the caption: the Thomas–Fermi screening wavevector of a metal, 1.5–2.1 Å⁻¹ for $r_s=2$–4).
- The four frequency presets (GaAs TO, 50 meV, hBN TO, 2 eV exciton) mixed materials and frequencies → one "medium" seg plus the
  $\hbar\omega$ slider.
- v1 colours `muted` for the zone band and `photon` for both light lines are fine; check every label in v2 notation
  (`ħω`, `q_{light}`, `ε_∞`, upright `L`, `T`).
