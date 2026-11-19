# B.9 figures · specifications for the figure authors

Section file: `v2/sections/B09-electrostatic.html` (id `b-electrostatic`, prefix `b9`). Read the section before you build:
the figures must use its notation and its equation numbers. Four figures, in reading order:

| File | Id | Kind | Placed in |
|---|---|---|---|
| `fig-limit.html` | `viz-b9-limit` | new (diagram with two toggles) | end of `#b9-picture` |
| `fig-sectors.html` | `viz-b9-sectors` | port of v1 `viz-b-electrostatic-sectors` | `#b9-model`, before `#b9-counting` |
| `fig-nonanalytic.html` | `viz-b9-nonanalytic` | new | end of `#b9-dd`, before `#b9-matrix` |
| `fig-plasmon.html` | `viz-b9-plasmon` | port of v1 `viz-b-electrostatic-plasmon` | end of `#b9-drude`, before `#b9-static` |

The two ports come from the v1 file `doc/epp-primer/sections/B09-electrostatic.html`: the `<figure class="viz" id="…">`
elements, and the matching blocks inside the single `<script>` at the end of that section (the block commented
`figure 1: sectors vs x = cq/omega` and the block commented `figure 2: LO-phonon–plasmon coupled modes`). The v1 section
script also draws two `div.fd` diagram rows; they are not ported as such, but `fig-limit` replaces them (see its entry).

Common rules for all four (from `v2/build/BUILDER_GUIDE.md` §5):

- Equation references in captions: `<a class="eqref" href="#b9-eq-…"></a>` (it renders "(B9.n)" with the parentheses; do not
  add your own). Ids used below: `b9-eq-free` (B9.1), `b9-eq-sectors` (B9.2), `b9-eq-closed` (B9.3), `b9-eq-iso` (B9.4),
  `b9-eq-limit` (B9.5), `b9-eq-We` (B9.6), `b9-eq-Wt` (B9.7), `b9-eq-piph` (B9.8), `b9-eq-pole` (B9.9), `b9-eq-static`
  (B9.10), `b9-eq-cubic` (B9.11), `b9-eq-gl` (B9.12), `b9-eq-first` (B9.13), `b9-eq-matrix` (B9.14), `b9-eq-ddG` (B9.15),
  `b9-eq-tot` (B9.16), `b9-eq-lopl` (B9.17), `b9-eq-Lpm` (B9.18), `b9-eq-tf` (B9.19).
- No source notation anywhere in a figure (labels, readouts, legends, aria-labels, strings): no $\epsilon$ (write
  $\varepsilon$), no $\Omega^{\rm e}$, $\Omega_{\rm uc}$, $j_\sigma$, $D_{11}$, $ll/lt/tt$, no "E.4xx" or "G Eq." labels, and no σ
  in any string (σ is reserved for the conductivity; the linter flags it). Branch index is $\nu$.
- Labels: upright labels as whole words (L, T, LT, TL, e, ph, dd, ion, LO, TO, p), bold vectors with `{bold:true}`
  ($\hat{\mathbf q}$, $\mathbf j_\nu$), italic variables. HTML readouts use MathJax `$…$`; canvas text uses `EPP.drawRich`.
- CSS classes prefixed with the figure name (`b9-limit-…`, `b9-sectors-…`, `b9-nonanalytic-…`, `b9-plasmon-…`).
- Colours read at draw time with `EPP.color(...)`: field lines `photon`, electronic bubbles and carriers `electron`,
  phonons `phonon`, ionic vertices and secondary curves `c4`, tertiary `c5`, faded items `muted`.
- Controls below the plots, 390 px with no horizontal overflow, keyboard access, `aria-label`s, no NaN, no console
  errors, fixed axis limits. Caption: "What to notice" (1–2 sentences) in `p.cap-key`, then Marks, Try, Model in
  `div.cap-more`, about 200–260 words in total.
- Test: `bash v2/build/figtest.sh v2/sections/B09/fig-<name>.html`, then look at both screenshots.

---

## 1. `fig-limit.html` · `viz-b9-limit` · new

**Title.** "From the dressed photon line to the Gonze–Lee term" `<span class="tag">interactive</span>` (two toggles).

**What the reader must learn.** The whole section in one picture. With a finite speed of light the propagator dressed by
electrons has four sectors, sorted by the line at each end (Coulomb or photon). As $c\to\infty$ the three sectors with a
photon end vanish, the purely longitudinal chain survives and sums to the screened interaction $W^{\mathrm e}$, and that
line, placed between two ionic vertices that see only $\bq\cdot\bj_\nu$, is the phonon self-energy. At $\omega=0$ it is
the Gonze–Lee dipole–dipole term; with $\omega$ kept it gives the plasmon–phonon modes.

**Include position.** Last element of `<h3 id="b9-picture">`, right after the plan list. Sentence just before:
"<figref viz-b9-limit> draws the whole chain as diagrams." (preceded by the numbered plan of the section).

**Physics** (no numerics; scalings at fixed $\bq,\omega$, from the section):

- free lines (B9.1): $D^{0\mathrm L}\propto c$, $D^{0\mathrm T}\propto 1/c$ for $cq\gg\omega$; bubble $\Pi^{\mathrm e}\propto 1/c$;
- sectors (B9.2): $D^{\mathrm{e,L}}\propto c$; $D^{\mathrm{e,LT}}$, $D^{\mathrm{e,TL}}$, $D^{\mathrm{e,T}}\propto 1/c$ (for the T sector this is
  its leading term, the free transverse line; the first-order term drawn falls as $c^{-3}$). On the physical scale $D/c$:
  L is of order $c^0$, the other three of order $c^{-2}$ (B9.5);
- survivor (B9.6): $\lim_{c\to\infty}D^{\mathrm e}_{\alpha\beta}=-\frac{c}{e^2\omega^2}q_\alpha q_\beta W^{\mathrm e}$,
  $W^{\mathrm e}=v+v\chi^0W^{\mathrm e}$; in the four-sector names of Section B.6, the survivor is $D_{\rho\rho}$;
- phonon self-energy (B9.8): two ionic vertices $\bq\cdot\bj^*_\nu$ and $\bq\cdot\bj_{\nu'}$ joined by $W^{\mathrm e}(\bq,\omega)$;
- at $\omega=0$: the dipole–dipole force constants $\Phi^{\mathrm{dd}}$ (B9.12); with $\omega$ kept: the pole condition (B9.9)
  and the coupled modes $\mathrm L_\pm$ (B9.18).

**What to draw.** SVG diagrams with `EPP.fd.draw` (item types `coulomb` = dashed, `photon` = wavy, `photon2` = heavy wavy,
`screened` = double dashed, `blob`, `phonon` = zigzag, `vertex`, `text`). Two rows; on phones each row wraps.

- Row (a), label "(a) four sectors, lowest order": four small diagrams, each "line – blob – line", 2×2 grid on phones,
  1×4 on desktop:
  - L–Π–L: dashed, blob, dashed. Caption under it: $D^{\mathrm{e,L}}\propto c$ and, smaller, "Coulomb–Coulomb, $D_{\rho\rho}$".
  - L–Π–T: dashed, blob, wavy. $D^{\mathrm{e,LT}}\propto 1/c$, "$D_{\rho\beta}$".
  - T–Π–L: wavy, blob, dashed. $D^{\mathrm{e,TL}}\propto 1/c$, "$D_{\alpha\rho}$".
  - T–Π–T: wavy, blob, wavy. $D^{\mathrm{e,T}}\propto 1/c$, "$D^{\mathrm T}_{\alpha\beta}$".
  Line labels $D^{0\mathrm L}$ (dashed) and $D^{0\mathrm T}$ (wavy) on the first diagram that uses each; blob label
  $\Pi^{\mathrm e}$ (colour `electron`). Field lines colour `photon`.
- Row (b), label "(b) the survivor and what it does": a chain of four panels with operators between them:
  1. heavy wavy line labelled $D^{\mathrm e}$;
  2. operator "$c\to\infty$" (an arrow with this label);
  3. dashed – blob – dashed – blob – dashed – "…" (labels $D^{0\mathrm L}$, $\Pi^{\mathrm e}$);
  4. operator "=";
  5. double-dashed line labelled $W^{\mathrm e}$, with a small caption "screened Coulomb, $W^{\mathrm e}=v+v\chi^0W^{\mathrm e}$";
  6. operator "between two ionic vertices";
  7. phonon zigzag – teal vertex – double-dashed $W^{\mathrm e}(\bq,\omega)$ – teal vertex – phonon zigzag; vertex labels
     $\hat{\bq}\cdot\mathbf j^*_\nu$ and $\hat{\bq}\cdot\mathbf j_{\nu'}$ written as `q·j_ν^*` and `q·j_{ν'}` with bold q, j
     (colour `c4-t`); caption "$\Pi^{\mathrm{ph}}_{\nu\nu'}$ (B9.8)";
  8. a final label that depends on the second toggle: "$\omega=0$: dipole–dipole term $\Phi^{\mathrm{dd}}$, Gonze–Lee
     (B9.12)" or "$\omega$ kept: dynamic screening, $\mathrm L_\pm$ modes (B9.18)"; in the second case the line label
     reads $W^{\mathrm e}(\bq,\omega)$, in the first $W^{\mathrm e}(\bq,0)$.

**Controls** (2): segmented "speed of light: finite | $c\to\infty$" (default finite) and segmented "screening:
$\omega$ kept | $\omega=0$" (default $\omega=0$). With $c$ finite all four diagrams of row (a) are in full colour and row
(b) is dimmed except panel 1. With $c\to\infty$ the LT, TL and T diagrams turn `muted` with a small "→ 0" badge, the L
diagram and row (b) are in full colour. No animation needed; if you add a fade, use `EPP.anim` (reduced motion).

**Readout** (one line, HTML): "on the physical scale $D/c$: L $\sim c^0$, LT, TL, T $\sim c^{-2}$" (static text) and the
current state: "surviving sectors: 4" or "surviving sectors: 1 ($D_{\rho\rho}$)".

**Initial state.** $c$ finite, $\omega=0$. Accessibility: each diagram SVG gets an `aria-label` that names its lines
(e.g. "longitudinal line, electronic bubble, transverse line: mixed sector LT, scales as 1 over c").

**Caption key points.**
- What to notice: at finite $c$ all four sectors exist; for $c\to\infty$ only the one with Coulomb lines at both ends
  survives, it sums to $W^{\mathrm e}$, and between two ionic vertices it becomes the phonon self-energy and, at
  $\omega=0$, the Gonze–Lee term.
- Marks: dashed rose = free longitudinal (Coulomb) line $D^{0\mathrm L}$; wavy rose = free transverse photon $D^{0\mathrm T}$;
  heavy wavy = $D^{\mathrm e}$; double dashed = $W^{\mathrm e}$; blue blob = electronic bubble $\Pi^{\mathrm e}$; amber zigzag =
  phonon; teal dots = ionic vertices, which in this limit see only $\bq\cdot\mathbf j_\nu$; grey = sectors that vanish.
  The scaling under each diagram is that of the whole sector at fixed $\bq$ and $\omega$; for T it comes from the free
  transverse line, while the first-order term drawn falls as $c^{-3}$.
- Try: switch $c$ to $\infty$ and see which sectors survive; switch the screening between $\omega=0$ and $\omega$ kept and
  see the last panel change from the Gonze–Lee term to the coupled LO-phonon–plasmon modes.
- Model: scalings from (B9.1) and the counting behind (B9.5); the survivor (B9.6); the self-energy (B9.8); the
  dipole–dipole term (B9.12); the coupled modes (B9.18). The sector names $D_{\rho\rho}$, $D_{\rho\beta}$, $D_{\alpha\rho}$,
  $D^{\mathrm T}_{\alpha\beta}$ are those of Section B.6.

---

## 2. `fig-sectors.html` · `viz-b9-sectors` · port of v1 `viz-b-electrostatic-sectors`

**Title.** "Which sector survives $c\to\infty$?" → in v2 a statement, not a question: "Four sectors, one survivor: the
dressed propagator of a uniaxial insulator" `<span class="tag">interactive</span>`.

**What the reader must learn.** In a concrete medium the longitudinal sector of $D^{\mathrm e}$ settles, outside the light
cone, on $1/(\hat{\bq}\cdot\boldsymbol\varepsilon\cdot\hat{\bq})$ and does not depend on $c$ (in units of $4\pi c/\omega^2$),
while the mixed and transverse sectors fall as $(cq/\omega)^{-2}$. At fixed $\bq$ and $\omega$, $cq/\omega\propto c$, so
moving right along the axis is sending $c\to\infty$. Phonons anywhere in the zone sit far to the right, deep in the
electrostatic regime. The mixed sector exists only in an anisotropic crystal, and vanishes for $\theta=0^\circ$, $90^\circ$
or $\varepsilon_{xx}=\varepsilon_{zz}$.

**Include position.** In `<h4 id="b9-model">`, after the derivation box of (B9.3) and before `<h4 id="b9-counting">`.
Sentences just before (the paragraph that ends before the derivation box): "… In an anisotropic crystal it does: the
polarization $\bP=\frac{1}{4\pi}(\boldsymbol\varepsilon-\mathbf 1)\bE$ of a longitudinal field has a component across
$\bq$. The mixed sectors are then nonzero, but still small outside the light cone. <figref viz-b9-sectors> computes all
four for hexagonal boron nitride."

**Physics.**
- Closed form (B9.3): $d\equiv\frac{\omega^2}{4\pi c}D^{\mathrm e}=\big[s^2(\mathbf 1-\hat{\bq}\hat{\bq}^{\mathsf T})-\boldsymbol\varepsilon\big]^{-1}$
  with $s=cq/\omega$ (use any internal name; on screen write "cq/ω"; do not call it $x$, which is a Cartesian direction).
- Uniaxial tensor $\boldsymbol\varepsilon=\mathrm{diag}(\varepsilon_{xx},\varepsilon_{xx},\varepsilon_{zz})\,(1+0.02\,\ii)$, optic axis $z$;
  the loss tangent 0.02 only regularizes the light-line poles.
- $\hat{\bq}=(\sin\theta,0,\cos\theta)$, $\mathbf t_1=(\cos\theta,0,-\sin\theta)$, $\mathbf t_2=(0,1,0)$.
- Plotted norms: $|{\rm L}|=|\hat{\bq}\cdot d\cdot\hat{\bq}|$; $|{\rm LT}|=\sqrt{|\hat{\bq}\cdot d\cdot\mathbf t_1|^2+|\hat{\bq}\cdot d\cdot\mathbf t_2|^2}$
  (equal to $|{\rm TL}|$ by symmetry of $d$); $|{\rm T}|$ = Frobenius norm of the $2\times2$ block on $\mathbf t_1,\mathbf t_2$.
- Limits to check: for $s\to\infty$, $|{\rm L}|\to1/(\hat{\bq}\cdot\boldsymbol\varepsilon\cdot\hat{\bq})$ (0.2558 at $\theta=45^\circ$
  with hBN values; 0.2557 with the 0.02 loss), $|{\rm LT}|$ and $|{\rm T}|\propto s^{-2}$ (at $\theta=45^\circ$: $|{\rm LT}|=2.46\times10^{-7}$, $|{\rm T}|=1.41\times10^{-6}$
  at $s=10^3$; $2.46\times10^{-11}$ and $1.41\times10^{-10}$ at $s=10^5$). Light lines of the medium at $s^2=\varepsilon$.
- The v1 code (`inv3`, `sectors`) does this correctly; keep it.
- Material: hexagonal boron nitride, $\varepsilon_{xx}=4.87$ (in plane), $\varepsilon_{zz}=2.95$ (along the c axis)
  (literature values, illustrative).
- Phonon band: phonons with energy $\hbar\omega$ and $10^{-3}\le q\le1$ Å⁻¹ have $cq/\omega=\hbar c\,q/\hbar\omega$,
  $\hbar c=1973.27$ eV Å. Default $\hbar\omega=170$ meV (the in-plane TO phonon of hBN): band from 11.6 to $1.16\times10^4$.

**What to draw.** One log–log plot (keep the v1 `EPP.Plot` set-up): $x$ axis "cq/ω" from 0.1 to $10^5$; $y$ axis
"norm of a sector of (ω²/4πc) D^e" from $10^{-11}$ to $100$, ticks every two decades. Curves: L solid (`photon`), LT = TL
dashed (`c4`), T dotted (`c5`). Horizontal muted line at $1/(\hat{\bq}\cdot\boldsymbol\varepsilon\cdot\hat{\bq})$ labelled
"1/(q̂·ε·q̂)" (bold q̂, `EPP.drawRich`). Muted dashed reference of slope −1 anchored on LT at $cq/\omega=10$ and drawn only
while LT ≠ 0 (the printed counting for T is slope −2 here), labelled "slope −1: printed power counting". Shaded band from 0.1 to $\sqrt{\max\varepsilon}$ labelled "light
cone"; shaded band for the phonons (amber, `phonon`, low alpha) labelled "phonons, 10⁻³–1 Å⁻¹". Text near the right end:
"LT, T ∝ (cq/ω)^{-2}" (or "T ∝ (cq/ω)^{-2}" when LT ≡ 0).

**Controls** (4): angle $\theta$ between $\hat{\bq}$ and the optic axis, 0–90°, step 1°, default 45°; $\varepsilon_{xx}$
(in plane) 1–15, step 0.01, default 4.87; $\varepsilon_{zz}$ (along the axis) 1–15, step 0.01, default 2.95; phonon energy
$\hbar\omega$ 10–200 meV, default 170 meV (moves only the phonon band).

**Readouts** (HTML, MathJax): $\hat{\bq}\cdot\boldsymbol\varepsilon\cdot\hat{\bq}$; $1/(\hat{\bq}\cdot\boldsymbol\varepsilon\cdot\hat{\bq})$;
$|{\rm L}|$ at $cq/\omega=10^5$; slope of LT for $cq/\omega>10^3$ (−2.00, or "none: LT ≡ 0"); slope of T for $cq/\omega>10^3$
(−2.00); phonon band in $cq/\omega$.

**Initial state.** $\theta=45^\circ$, hBN values, $\hbar\omega=170$ meV: all three curves visible, LT nonzero, slopes −2.

**Caption key points.**
- What to notice: outside the light cone the longitudinal sector settles on $1/(\hat{\bq}\cdot\boldsymbol\varepsilon\cdot\hat{\bq})$,
  independent of $c$, while the mixed and transverse sectors fall as $(cq/\omega)^{-2}$, faster than the slope −1 that the
  EPP manuscript prints for the mixed sectors. Every phonon of the zone lies far to the right.
- Marks: solid rose L, dashed teal LT (= TL), dotted T; grey line $1/(\hat{\bq}\cdot\boldsymbol\varepsilon\cdot\hat{\bq})$; dashed
  grey slope −1; rose band = light cone ($cq\lesssim\sqrt{\varepsilon}\,\omega$, where the transverse poles sit); amber band
  = phonons of energy $\hbar\omega$ with $q$ from $10^{-3}$ to 1 Å⁻¹. At fixed $\bq$ and $\omega$, $cq/\omega\propto c$, so moving
  right is sending $c\to\infty$.
- Try: set $\theta$ to 0° or 90°, or $\varepsilon_{xx}=\varepsilon_{zz}$: LT vanishes (in an isotropic medium a longitudinal
  field never feeds the transverse one). Make the crystal more anisotropic: LT grows but keeps its slope. Move $\hbar\omega$:
  the phonon band shifts but stays far outside the light cone.
- Model: (B9.3) for a frequency-independent dielectric tensor, i.e. $\Pi^{\mathrm e}=\frac{\omega^2}{4\pi c}(\boldsymbol\varepsilon-\mathbf 1)$;
  the general counting behind (B9.5); loss tangent 0.02 to keep the light-line poles finite; hBN values (literature,
  illustrative). The physical self-energy uses $D^{\mathrm e}/c$, which is exactly what is plotted.

**What in v1 is wrong or in source notation, and must change.**
- Id and classes: `viz-b-electrostatic-sectors` → `viz-b9-sectors`; the v1 figure has generic classes `.ps .ctrls .ro`
  inside a shared section script → own file, own IIFE, classes `b9-sectors-…`.
- Notation: $\epsilon$ → $\varepsilon$; $d^{\rm e}$ → "(ω²/4πc) D^e"; sector names "ll", "lt = tl", "tt" → "L", "LT = TL",
  "T"; $\epsilon_\perp$, $\epsilon_\parallel$ → $\varepsilon_{xx}$, $\varepsilon_{zz}$ (the ∥ symbol means "in plane" in the 2D
  notation of the primer); "x = cq/ω" → "cq/ω"; readout key "ε_L = q̂·ε·q̂" → "q̂·ε·q̂" (no symbol $\varepsilon_L$).
- Source labels to remove: "1/ε_L (E.446)", "slope −1, as E.440 implies for lt", "slope −1 (E.440)"; the caption's
  references to E.388, E.450, E.446, E.440, $\Omega^{\rm e}_{\sigma\sigma}$, $\bj_\sigma$, $\Omega_{\rm uc}$ and its "[derived]" marks.
- The amber band was fixed at $\hbar\omega\approx36$ meV, a GaAs number, while the medium is hBN: tie it to the new
  $\hbar\omega$ control, default 170 meV.
- Caption: rewrite in the four-part structure; "Model" first in v1 must move to the end.

---

## 3. `fig-nonanalytic.html` · `viz-b9-nonanalytic` · new

**Title.** "At Γ the frequencies depend on the direction of approach" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. The dipole–dipole term is nonanalytic: in a uniaxial crystal (wurtzite GaN) the frequencies at $\bq\to0$ change with
   the angle $\theta$ between $\hat{\bq}$ and the optic axis $z$. The LO-like mode moves from $\omega_{\mathrm{LO},z}$ at
   $\theta=0$ to $\omega_{\mathrm{LO},x}$ at $\theta=90^\circ$.
2. First-order perturbation theory on the bare modes, (B9.13), which is the diagonal approximation, fails when the
   dipole–dipole term mixes bare modes that are nearly degenerate. The full matrix condition (B9.14), which is what
   DFPT does, gives one mode polarized nearly along $\hat{\bq}$ (LO-like) and one across it (TO-like).
3. When the bare TO splitting exceeds the dipole–dipole coupling (weakly polar crystal), the modes stay along the
   crystal axes and first order works.

**Include position.** End of `<h4 id="b9-dd">`, before `<h4 id="b9-matrix">`. Sentences just before: "… In a uniaxial
crystal, such as wurtzite GaN, the LO and TO frequencies themselves change with the angle $\theta$ between $\hat\bq$ and
the optic axis $z$, between their values $\omega_{\mathrm{LO},z}$, $\omega_{\mathrm{TO},z}$ for modes polarized along $z$ and
$\omega_{\mathrm{LO},x}$, $\omega_{\mathrm{TO},x}$ for modes polarized across it (<figref viz-b9-nonanalytic>)." followed by
the paragraph "Acoustic modes are left alone at $\bq\to0$. … discussed in Section B.1." The next subsection (`b9-matrix`)
quotes the 45° numbers of this figure.

**Physics.** Keep only the three polar optical modes at Γ (Loudon's model of a uniaxial crystal). Bare modes, from the
short-range lattice: $\nu=x,y$ (the $\mathrm E_1$ pair, bare frequency $\omega_{\mathrm{TO},x}$) and $\nu=z$ (the $\mathrm A_1$ mode,
bare frequency $\omega_{\mathrm{TO},z}$). Mode currents $\mathbf j_x=j_x\hat{\mathbf x}$, $\mathbf j_y=j_x\hat{\mathbf y}$,
$\mathbf j_z=j_z\hat{\mathbf z}$, written through mode plasma frequencies
$(\omega^{\mathrm{ion}}_{\mathrm p,x})^2=4\pi j_x^2/\OV=\varepsilon_{xx}(\omega^2_{\mathrm{LO},x}-\omega^2_{\mathrm{TO},x})$ and
$(\omega^{\mathrm{ion}}_{\mathrm p,z})^2=\varepsilon_{zz}(\omega^2_{\mathrm{LO},z}-\omega^2_{\mathrm{TO},z})$ (the generalization of (B9.11)).
With (B9.8) at $\omega=0$ in (B9.14), at $\bq\to0$ along $\hat{\bq}=(\sin\theta,0,\cos\theta)$, the matrix
$\frac2\hbar\sqrt{\omega_\nu\omega_{\nu'}}\Pi^{\mathrm{ph}}_{\nu\nu'}=\frac{4\pi}{\OV}\frac{(\hat{\bq}\cdot\mathbf j_\nu)(\hat{\bq}\cdot\mathbf j_{\nu'})}{\hat{\bq}\cdot\boldsymbol\varepsilon_\infty\cdot\hat{\bq}}$
is the nonanalytic limit of (B9.12) in the mode basis. With $a(\theta)=1/(\varepsilon_{xx}\sin^2\theta+\varepsilon_{zz}\cos^2\theta)$
and a Born-charge scale $s$ (multiplies every $Z^*$, so the plasma frequencies by $s$):

$$M(\theta)=\begin{pmatrix}\omega^2_{\mathrm{TO},x}+a\,s^2(\omega^{\mathrm{ion}}_{\mathrm p,x})^2\sin^2\theta & a\,s^2\,\omega^{\mathrm{ion}}_{\mathrm p,x}\omega^{\mathrm{ion}}_{\mathrm p,z}\sin\theta\cos\theta\\ a\,s^2\,\omega^{\mathrm{ion}}_{\mathrm p,x}\omega^{\mathrm{ion}}_{\mathrm p,z}\sin\theta\cos\theta & \omega^2_{\mathrm{TO},z}+a\,s^2(\omega^{\mathrm{ion}}_{\mathrm p,z})^2\cos^2\theta\end{pmatrix}$$

in the basis of the bare $x$ and $z$ modes; the $y$ mode is transverse for every $\theta$ in this plane and stays at
$\omega_{\mathrm{TO},x}$ (the "ordinary" TO).
- Exact (matrix, B9.14): the two eigenvalues of $M$ (closed-form 2×2), $\hbar\omega=\sqrt{\cdot}$; the upper one is
  "LO-like", the lower one "TO-like". Eigenvector $(c_x,c_z)$ of each: displacement direction at angle
  $\phi=\mathrm{atan2}(c_x,c_z)$ from the $z$ axis, taken mod 180°.
- First order (diagonal, B9.13): $\hbar\omega_x=\hbar\sqrt{M_{xx}}$ and $\hbar\omega_z=\hbar\sqrt{M_{zz}}$, following the bare
  $x$ and $z$ modes.
- Mixing strength readout: $|M_{xz}|/|M_{xx}-M_{zz}|$ (large ⇒ first order fails).
- Material, wurtzite GaN (literature values, illustrative): $\hbar\omega_{\mathrm{TO},z}=65.9$ meV ($\mathrm A_1$, 532 cm⁻¹),
  $\hbar\omega_{\mathrm{TO},x}=69.3$ meV ($\mathrm E_1$, 559 cm⁻¹), $\hbar\omega_{\mathrm{LO},z}=91.0$ meV (734 cm⁻¹),
  $\hbar\omega_{\mathrm{LO},x}=91.9$ meV (741 cm⁻¹), $\varepsilon_{xx}=\varepsilon_{zz}=5.35$. Then
  $\hbar\omega^{\mathrm{ion}}_{\mathrm p,x}=139.6$ meV and $\hbar\omega^{\mathrm{ion}}_{\mathrm p,z}=145.1$ meV. (1 meV = 8.06554 cm⁻¹.)
- TO anisotropy control $\Delta=\hbar(\omega_{\mathrm{TO},x}-\omega_{\mathrm{TO},z})$ keeps $\hbar\omega_{\mathrm{TO},x}=69.3$ meV and the two
  plasma frequencies fixed; GaN has $\Delta=3.35$ meV.
- Check values (computed with this model): GaN ($s=1$, $\Delta=3.35$ meV) at $\theta=45^\circ$: exact 91.46 and 67.59 meV,
  first order 81.37 ($x$) and 79.46 ($z$) meV, LO-like displacement at 47.3° from $z$, mixing 6.2; at $\theta=15^\circ$: exact
  91.07 and 69.05 meV, first order 71.02 and 89.54 meV, LO-like at 16.2°. At $\theta=0$ and 90° first order is exact
  (91.00/69.28 and 91.87/65.93 meV). With $s=0.2$ at 45°: exact 69.90 and 66.43 meV, first order 69.81 and 66.53 meV
  (mixing 0.17). With $\Delta=0$, $s=1$ at 45°: exact 92.67 and 69.28 meV. Largest LO at $s=1.5$: 114.9 meV.

**What to draw.**
- Panel (a), plot: $x$ axis $\theta$ (deg) 0–90; $y$ axis $\hbar\omega$ (meV), fixed range 45–120. Exact branches solid:
  LO-like (`phonon`), TO-like (`c4`), ordinary TO along $y$ (`c5`, flat). First-order curves dashed thin (`muted` or
  `photon`), labelled "first order, $x$ mode" and "first order, $z$ mode". Faint horizontal reference ticks at the right edge for
  $\omega_{\mathrm{LO},z}$, $\omega_{\mathrm{LO},x}$, $\omega_{\mathrm{TO},z}$, $\omega_{\mathrm{TO},x}$. Vertical cursor at the current $\theta$
  with dots on the exact curves.
- Panel (b), geometry (beside (a) on desktop, below on phones): the quarter plane $x$–$z$ with the $z$ axis vertical,
  labelled "z, optic axis", and the $x$ axis horizontal. A bold arrow $\hat{\bq}$ at angle $\theta$ from $z$ (`ink`, label
  with `{bold:true}`). Through the origin, double-headed arrows for the displacement directions of the two extraordinary
  modes: LO-like (`phonon`), TO-like (`c4`), each labelled. The ordinary TO along $y$ as a ⊙ symbol at the origin (out of the
  plane), labelled "TO along y". Thin dashed arrows along $x$ and $z$ for the bare modes. Arc marks for $\theta$ and for the
  angle between the LO-like direction and $\hat{\bq}$.
- Small legend under (a).

**Controls** (4): $\theta$, 0–90°, step 1°, default 45°; Born-charge scale $s$, 0–1.5, step 0.01, default 1 (1 = GaN);
TO anisotropy $\Delta$, −20 to 20 meV, step 0.05, default 3.35 meV; toggle "show first order (diagonal)", default on.

**Readouts** (HTML): $\hbar\omega$ of the LO-like, TO-like and ordinary modes at $\theta$; first order $\hbar\omega_x$,
$\hbar\omega_z$; angle between the LO-like displacement and $\hat{\bq}$; mixing strength $|M_{xz}|/|M_{xx}-M_{zz}|$.

**Initial state.** GaN, $\theta=45^\circ$, first order shown: the two dashed first-order curves cross near 80 meV between
the exact branches at 91.5 and 67.6 meV, and the LO-like arrow in (b) lies close to $\hat{\bq}$.

**Caption key points.**
- What to notice: in GaN the LO-like frequency moves between 91.0 and 91.9 meV as $\hat{\bq}$ turns from the optic axis to
  the plane: at Γ the frequency depends on the direction of approach. First order, which follows the bare $x$ and $z$
  modes, is far off at intermediate angles; the true modes are polarized along and across $\hat{\bq}$.
- Marks: solid curves, exact modes from the matrix condition; dashed, first order; in (b) the bold arrow is $\hat{\bq}$, the
  coloured double arrows the displacement directions of the two modes in the $x$–$z$ plane, ⊙ the ordinary TO mode along $y$.
- Try: (1) $\theta=0^\circ$ or $90^\circ$: first order is exact, nothing mixes. (2) $s=0.2$: the dipole–dipole coupling falls
  below the bare TO splitting; the modes stay along $x$ and $z$ and first order agrees with the exact curves. (3) $\Delta=0$:
  the bare modes are degenerate, and any polarity mixes them completely.
- Model: the nonanalytic limit of (B9.12) in the basis of the three polar modes; matrix condition (B9.14) with (B9.8) at
  $\omega=0$; first order (B9.13); $\varepsilon_\infty$ isotropic 5.35; GaN parameters (literature, illustrative). Wurtzite has
  further modes ($\mathrm E_2$, $\mathrm B_1$) that carry no dipole and do not enter.

---

## 4. `fig-plasmon.html` · `viz-b9-plasmon` · port of v1 `viz-b-electrostatic-plasmon`

**Title.** "LO-phonon–plasmon modes: dynamic against static screening" `<span class="tag">interactive</span>`.

**What the reader must learn.** With the frequency dependence of $W^{\mathrm e}$ kept, a doped polar semiconductor has two
coupled longitudinal modes $\mathrm L_\pm$ that anticross where the plasmon meets $\omega_{\mathrm{LO}}$, and the lattice weight
passes from $\mathrm L_+$ to $\mathrm L_-$ as the density grows. Part A's adiabatic step (static screening) gives a single mode
that slides down to $\omega_{\mathrm{TO}}$ and misses the plasmon branch. At larger $q$ ($qv_{\mathrm F}\gg\omega_{\mathrm{LO}}$) the ridges
move onto the static curve and the plasmon dissolves in the particle–hole continuum, where the static step is right.

**Include position.** End of `<h4 id="b9-drude">`, right before `<h4 id="b9-static">`. Sentences just before: "… At
$n_{\mathrm c}=10^{18}$ cm⁻³, with $\hbar\omega_{\mathrm{TO}}=33.3$ and $\hbar\omega_{\mathrm{LO}}=36.2$ meV, the modes sit at
$\hbar\omega_{\mathrm L+}=47.8$ meV and $\hbar\omega_{\mathrm L-}=30.3$ meV. The lower mode always lies below both $\omega_{\mathrm{TO}}$
and $\omega_{\mathrm p}$, the upper one above both $\omega_{\mathrm{LO}}$ and $\omega_{\mathrm p}$. Raman scattering on n-type GaAs observed
this pair, after Varga had predicted it." The text after the figure (`b9-static`, `b9-window`) refers back to it.

**Physics** (Gaussian; meV, Å; $e^2=14.3996$ eV Å, $\hbar^2/m_{\mathrm e}=7.61996$ eV Å², $\hbar^2/2m_{\mathrm e}=3.80998$ eV Å²).
- Carriers (one parabolic band, $T=0$): $n_{\mathrm c}$ [cm⁻³] × $10^{-24}$ = Å⁻³; $k_{\mathrm F}=(3\pi^2n_{\mathrm c})^{1/3}$;
  $E_{\mathrm F}=\hbar^2k_{\mathrm F}^2/2m^*$; $\hbar v_{\mathrm F}=\hbar^2k_{\mathrm F}/m^*$; screened plasma frequency
  $(\hbar\omega_{\mathrm p})^2=4\pi e^2\hbar^2n_{\mathrm c}/(\varepsilon_\infty m^*)$; $k^2_{\mathrm{TF}}=6\pi n_{\mathrm c}e^2/E_{\mathrm F}$.
- Ionic strength from the undoped crystal: $(\omega^{\mathrm{ion}}_{\mathrm p})^2=\varepsilon_\infty(\omega^2_{\mathrm{LO}}-\omega^2_{\mathrm{TO}})$ (B9.11).
- Panel (a), dynamic, $q\to0$ (Drude): $\omega_{\mathrm L\pm}$ from (B9.18); plasmon $\omega_{\mathrm p}$; resonance density where
  $\omega_{\mathrm p}=\omega_{\mathrm{LO}}$.
- Panel (a), static at the chosen $q$: (B9.19) with the static Lindhard factor,
  $\varepsilon^{\mathrm e}(q,0)=\varepsilon_\infty+\frac{k^2_{\mathrm{TF}}}{q^2}F_0\!\big(\frac{q}{2k_{\mathrm F}}\big)$,
  $F_0(z)=\frac12+\frac{1-z^2}{4z}\ln\big|\frac{1+z}{1-z}\big|$ ($F_0\to1$ for $z\to0$, $F_0(1)=\frac12$), and
  $\bar\omega^2=\omega^2_{\mathrm{TO}}+(\omega^{\mathrm{ion}}_{\mathrm p})^2/\varepsilon^{\mathrm e}(q,0)$.
- Panel (b), heat map at the chosen $q$, with the full $T=0$ Lindhard function
  $\varepsilon^{\mathrm e}(q,\omega)=\varepsilon_\infty+\frac{k^2_{\mathrm{TF}}}{q^2}F(z,u)$, $z=q/2k_{\mathrm F}$,
  $u=(\omega+\ii\gamma_{\mathrm e})/qv_{\mathrm F}$,
  $F(z,u)=\frac12+\frac1{8z}\Big\{[1-(z-u)^2]\ln\frac{z-u+1}{z-u-1}+[1-(z+u)^2]\ln\frac{z+u+1}{z+u-1}\Big\}$ (principal complex
  logarithms), and for $|u|\gg1+z$ the expansion $F\simeq-\frac{1}{3u^2}-\big(\frac15+\frac{z^2}3\big)\frac1{u^4}$. The v1
  function `lind(z, ur, ui)` implements exactly this; keep it. Two quantities, chosen by a control:
  - loss function $-\mathrm{Im}\,1/\varepsilon(q,\omega)$ with the total dielectric function (B9.16)
    $\varepsilon=\varepsilon^{\mathrm e}(q,\omega)+\frac{(\omega^{\mathrm{ion}}_{\mathrm p})^2}{\omega^2_{\mathrm{TO}}-\omega^2-\ii\gamma_{\mathrm{ph}}\omega}$; it
    is the spectral weight of the total screened interaction $W=v/\varepsilon$ of Section A.5;
  - phonon spectral function $-\mathrm{Im}\,F_\nu(q,\omega)$ in meV⁻¹, with
    $F_\nu=\frac{2\hbar\omega_{\mathrm{TO}}}{(\hbar\omega)^2-(\hbar\omega_{\mathrm{TO}})^2-(\hbar\omega^{\mathrm{ion}}_{\mathrm p})^2/\varepsilon^{\mathrm e}(q,\omega)+\ii\,\hbar\gamma_{\mathrm{ph}}\,\hbar\omega}$
    (the diagonal propagator with the self-energy (B9.8); same poles as the loss function, different weights).
  Display broadening $\hbar\gamma_{\mathrm e}=2.5$ meV, $\hbar\gamma_{\mathrm{ph}}=0.8$ meV. Plot $\log_{10}$ of the quantity, colour
  range $10^{-3}$ to 1 (as in v1).
- Presets (literature values, illustrative): GaAs $\hbar\omega_{\mathrm{TO}}=33.3$, $\hbar\omega_{\mathrm{LO}}=36.2$ meV,
  $\varepsilon_\infty=10.9$, $m^*=0.067\,m_{\mathrm e}$, $y$ range 0–70 meV; GaN ($\mathrm E_1$ modes, the same values as
  `viz-b9-nonanalytic`) $\hbar\omega_{\mathrm{TO}}=69.3$, $\hbar\omega_{\mathrm{LO}}=91.9$ meV, $\varepsilon_\infty=5.35$, $m^*=0.20\,m_{\mathrm e}$,
  $y$ range 0–150 meV.
- Check values: GaAs at $10^{18}$ cm⁻³: $k_{\mathrm F}=0.0309$ Å⁻¹, $E_{\mathrm F}=54.4$ meV, $\hbar v_{\mathrm F}=3.52$ eV Å,
  $\hbar\omega_{\mathrm p}=43.45$ meV, $k_{\mathrm{TF}}=0.0706$ Å⁻¹; $\hbar\omega_{\mathrm L+}=47.75$, $\hbar\omega_{\mathrm L-}=30.30$ meV;
  resonance $6.9\times10^{17}$ cm⁻³; $\omega_{\mathrm{LO}}/qv_{\mathrm F}=3.4$ at $q=3\times10^{-3}$ Å⁻¹. GaN: resonance
  $6.6\times10^{18}$ cm⁻³; $\hbar\omega_{\mathrm p}=35.9$ meV at $10^{18}$ cm⁻³.

**What to draw.** Two panels stacked on phones, side by side on desktop (keep the v1 layout).
- (a) mode energies against carrier density: $x$ axis $n_{\mathrm c}$ (cm⁻³), log, $10^{16}$–$10^{20}$, ticks at each decade;
  $y$ axis $\hbar\omega$ (meV). Curves: $\mathrm L_+$ and $\mathrm L_-$ solid amber (`phonon`), labelled "L₊", "L₋" (upright L);
  plasmon $\omega_{\mathrm p}$ dotted blue (`electron`); static (B9.19) at the chosen $q$ dashed teal (`c4`); horizontal muted
  lines at $\omega_{\mathrm{TO}}$ and $\omega_{\mathrm{LO}}$ (labels "ω_{TO}", "ω_{LO}, undoped"); hollow circle at the resonance
  ($\omega_{\mathrm p}=\omega_{\mathrm{LO}}$); dashed vertical cursor at the chosen density.
- (b) heat map of $\log_{10}$ of the chosen spectral function over $n_{\mathrm c}$ and $\hbar\omega$ (same axes as (a)), with
  colour bar; overlays: $\mathrm L_\pm$ ($q\to0$) muted dotted, the static curve dashed teal, the density cursor. Colour
  `photon` for the loss function, `phonon` for the phonon spectral function.

**Controls** (5): material (segmented GaAs | GaN, default GaAs; switching resets $m^*$ and the cursor: GaAs $10^{18}$,
GaN $10^{19}$ cm⁻³); carrier density $n_{\mathrm c}$ (range slider on $\log_{10}n_{\mathrm c}$ from 16 to 20, default 18;
clicking or dragging on either panel also sets it, as in v1, but the slider gives keyboard access); wavevector $q$
(log, $10^{-3}$–$10^{-1}$ Å⁻¹, default $3\times10^{-3}$); effective mass $m^*/m_{\mathrm e}$ (0.02–1, default 0.067); panel (b)
shows (segmented: loss function | phonon spectral function, default loss function).

**Readouts** (HTML, MathJax): $n_{\mathrm c}$; $\hbar\omega_{\mathrm p}$; $\hbar\omega_{\mathrm L+}$, $\hbar\omega_{\mathrm L-}$; static
$\hbar\bar\omega$ at $q$; $\omega_{\mathrm{LO}}/qv_{\mathrm F}$; $k_{\mathrm{TF}}/\sqrt{\varepsilon_\infty}$; resonance density.

**Initial state.** GaAs, $n_{\mathrm c}=10^{18}$ cm⁻³, $q=3\times10^{-3}$ Å⁻¹, loss function: the two dynamic branches
straddle the resonance circle, the static curve lies at $\omega_{\mathrm{TO}}$ over most of the range, and the heat map shows
two ridges on the dotted $\mathrm L_\pm$ curves.

**Caption key points.**
- What to notice: static screening only slides the LO mode down to $\omega_{\mathrm{TO}}$; the dynamic solution has two branches
  that anticross where the plasmon meets $\omega_{\mathrm{LO}}$, and the phonon spectral function shows the lattice weight
  passing from $\mathrm L_+$ to $\mathrm L_-$.
- Marks: (a) amber $\mathrm L_\pm$ from (B9.18) (Drude carriers, $q\to0$); dotted blue the plasmon $\omega_{\mathrm p}$; dashed teal
  the adiabatic result (B9.19) at the chosen $q$, with the static Lindhard function; grey lines $\omega_{\mathrm{TO}}$, $\omega_{\mathrm{LO}}$;
  circle the resonance $\omega_{\mathrm p}=\omega_{\mathrm{LO}}$. (b) heat map on a log scale, broadened for display.
- Try: raise $q$ until $\omega_{\mathrm{LO}}/qv_{\mathrm F}\ll1$: the ridges move onto the static curve and the plasmon dissolves in
  the particle–hole continuum, where the static step is right. Switch (b) to the phonon spectral function and move the
  density through the resonance. Switch to GaN: the same physics at ten times the density.
- Model: carriers in one parabolic band at $T=0$, $\omega^2_{\mathrm p}=4\pi n_{\mathrm c}e^2/\varepsilon_\infty m^*$; LO mode of a
  cubic polar crystal (B9.16) with $(\omega^{\mathrm{ion}}_{\mathrm p})^2$ from the undoped crystal; full Lindhard function in (b);
  broadenings $\hbar\gamma_{\mathrm e}=2.5$, $\hbar\gamma_{\mathrm{ph}}=0.8$ meV; GaAs and GaN values (literature, illustrative).
  Raman scattering on n-type GaAs observed the $\mathrm L_\pm$ pair (Mooradian and Wright, 1966); no data are plotted.

**What in v1 is wrong or in source notation, and must change.**
- Id and classes: `viz-b-electrostatic-plasmon` → `viz-b9-plasmon`; own file and IIFE; classes `b9-plasmon-…`.
- Notation: v1 $\omega_p$ is the **ionic** plasma frequency → $\omega^{\mathrm{ion}}_{\mathrm p}$; v1 $\omega_{\rm pl}$ (unscreened carrier
  plasma frequency) and $\tilde\omega_{\rm pl}=\omega_{\rm pl}/\sqrt{\epsilon_\infty}$ → v2 $\omega_{\mathrm p}$, which is the **screened**
  one, $\omega^2_{\mathrm p}=4\pi n_{\mathrm c}e^2/(\varepsilon_\infty m^*)$ (the v1 variable `wt` = $\hbar\omega_{\mathrm p}$; `wpl2` is
  $\varepsilon_\infty(\hbar\omega_{\mathrm p})^2$). $\epsilon_{\rm tot}$ → $\varepsilon$ (no superscript = total); $\epsilon_\infty$ →
  $\varepsilon_\infty$; $F_\sigma$ → $F_\nu$; $n$ → $n_{\mathrm c}$; "L_±" → upright $\mathrm L_\pm$, $\omega_{\mathrm L\pm}$.
- Source labels to remove: "static E.454 at q" → "static, adiabatic (B9.19)" (canvas: "static at q"); "L_± dynamic E.453
  (q → 0)" → "L₊, L₋, dynamic (q → 0)"; readout keys "L_+, L_− (E.453)", "static (E.454)", "ħω_{pl}/√ε_∞" → the readouts
  listed above; caption references to E.453, E.454, E.405, "G Eq. 104", "EPP's digression on p. 88", "Giustino's
  $W_{\rm e}+W_{\rm ph}$" → v2 equation references and "the total screened interaction $W=v/\varepsilon$ of Section A.5".
- No string may contain σ (the v1 heat-map label "−Im F_σ" fails the linter) → "−Im F_ν (meV⁻¹)".
- GaN preset 69.2/92.0 meV → 69.3/91.9 meV, the values used in `viz-b9-nonanalytic` and in the section's note on GaN.
- Accessibility: v1 set the density only by clicking the panels; add the density slider.
- v1 caption: "Model" came first and the SI form of the plasma frequency was in the caption; rewrite in the four-part
  structure with the Gaussian form (the SI form is in the section's notes).
