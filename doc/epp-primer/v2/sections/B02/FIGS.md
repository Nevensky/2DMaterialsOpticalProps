# B.2 figure specs (section `b-lagrangian`, prefix `b2`)

Three figures, in reading order:

| # | file | id | kind | v1 source |
|---|---|---|---|---|
| B2.1 | `B02/fig-map.html` | `viz-b2-map` | port | `viz-b-lagrangian-map`, `doc/epp-primer/sections/B02-lagrangian.html` (figure markup lines 62–90, script lines 160–274) |
| B2.2 | `B02/fig-lorentz.html` | `viz-b2-lorentz` | new | none |
| B2.3 | `B02/fig-drive.html` | `viz-b2-drive` | port | `viz-b-lagrangian-drive`, same v1 file (figure markup lines 134–141, script lines 276–380) |

General rules for all three (BUILDER_GUIDE §5): v2 notation only (NOTATION.md), no source symbols and no source
equation numbers anywhere in the figure (labels, legends, readouts, aria-labels, info cards, caption). Refer to equations
only with `<a class="eqref" href="#…"></a>` (ids below) and to sections with `<a href="#b-fields">Section B.4</a>`.
CSS classes prefixed `b2-map-`, `b2-lorentz-`, `b2-drive-`. Colours read at draw time with `EPP.color(...)`.
Caption: `p.cap-key` "What to notice" + `div.cap-more` with Marks / Try / Model, about 200–260 words in total.

Symbols of B.2 that the figures use (all defined in the section, glossary ids in brackets):
- $s$, $\br_s$, $\mathbf v_s$: electron label, position, velocity [`gl-b2-s`]; electron charge $-e$ with $e>0$, mass
  $m_{\mathrm e}$ [`gl-b2-const`].
- $\bA$, $A_0$, $\bE$, $\bB$ [`gl-b2-A`]; $\bA^{\mathrm L}$, $\bA^{\mathrm T}$, $\bE^{\mathrm L}$, $\bE^{\mathrm T}$ longitudinal and
  transverse parts [`gl-b2-LT`].
- $V^{\mathrm{en}}_0$ frozen crystal potential [`gl-b2-V0`]; $\Phi^{\mathrm{sr}}$ short-range force constants (home A.4).
- $\mathbf Z^*_\kappa$ Born charge (dimensionless; the dipole of atom $\ell\kappa$ is $e\mathbf Z^*_\kappa\cdot\bu_{\ell\kappa}$).
- $\bP$ polarization, $\bj^{\mathrm e}$, $\bj^{\mathrm{ion}}=\partial_t\bP$, $\rho^{\mathrm{ion}}=-\nabla\cdot\bP$ [`gl-b2-j`, `gl-b2-P`].
- $\mathbf p_s=m_{\mathrm e}\mathbf v_s-\frac ec\bA$ canonical momentum [`gl-b2-ps`]. **$\mathbf p$ is a momentum in v2, never a
  dipole.** Write dipoles in words or as $e\mathbf Z^*\cdot\bu$; write the polarization as $\bP$ or, along the field, its component $P_x$ ($\bP$ is never a scalar).
- $\omega_{\mathrm c}=eB/m_{\mathrm e}c$ [`gl-b2-wc`]; $Z^*$, $K$, $u_{\mathrm{rel}}$, $\omega_0$, $E_0$, $\varphi$ of the driven
  diatomic crystal [`gl-b2-drive`]; $\gamma$ damping rate; $\xi$ complex amplitude; $M_{\mathrm r}$ reduced mass.
- $\chi$ is reserved for density responses and $\sigma$ for the conductivity: do not use them for a susceptibility or a width.

Equation ids of B.2 (tags in brackets): `b2-eq-pot` (B2.1), `b2-eq-gauge` (B2.2), `b2-eq-lt` (B2.3), `b2-eq-EL` (B2.4),
`b2-eq-constraint` (B2.5), `b2-eq-L` (B2.6), `b2-eq-V0` (B2.7), `b2-eq-field` (B2.8), `b2-eq-je` (B2.9), `b2-eq-P` (B2.10),
`b2-eq-Ldens` (B2.11), `b2-eq-dipole` (B2.12), `b2-eq-site` (B2.13), `b2-eq-euler` (B2.14), `b2-eq-ps` (B2.15),
`b2-eq-lorentz` (B2.16), `b2-eq-ions` (B2.17), `b2-eq-rel` (B2.18), `b2-eq-xi` (B2.19).
In the Lagrangian (B2.6) each line is one term, called term 1 to term 6 from the top (no labels are printed):
1. $\sum_s[\frac12m_{\mathrm e}\mathbf v_s^2-V^{\mathrm{en}}_0(\br_s)]$; 2. $\sum\frac12M_\kappa\dot u^2_{\ell\kappa\alpha}$;
3. $-\frac12\sum\Phi^{\mathrm{sr}}uu$; 4. $\frac1{8\pi}\int\dd^3r\,(E^2-B^2)$; 5. $-\frac ec\sum_s\mathbf v_s\cdot\bA(\br_s,t)$;
6. $\frac ec\sum_{\ell\kappa}(\mathbf Z^*_\kappa\cdot\dot\bu_{\ell\kappa})\cdot\bA(\btau_{\ell\kappa},t)$.

---

## B2.1 · `fig-map.html` · `viz-b2-map` · port of `viz-b-lagrangian-map`

**Title.** "Who talks to whom: the Lagrangian as a coupling map" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. In the EPP Lagrangian every line that leaves the electrons or the ions ends on the field; there is no direct
   electron–ion line (no $g$).
2. Each line is one term of (B2.6), split into its longitudinal part (Coulomb, instantaneous) and its transverse part
   (light, finite $c$).
3. Taking $c\to\infty$ removes the transverse lines and leaves only Coulomb lines: the electrostatic world of Part A.
4. In Part A's Hamiltonian the field is gone, the longitudinal field has collapsed into an instantaneous Coulomb line, and an
   explicit vertex $g$ joins electrons and phonons. The physics inside $g$ travels along the two field paths in EPP.

**Include position.** Start of `h3#b2-noeph`, right after:
> No term of (B2.6) contains both an electron coordinate and an ionic displacement. There is no $g$. Figure B2.1 draws the
> structure: every line that leaves the electrons or the ions ends on the field.

**Physics (v2 notation).** Purely structural: the terms 1–6 of `b2-eq-L`, the split (B2.3) `b2-eq-lt` into
$\bA^{\mathrm L}$ and $\bA^{\mathrm T}$, the Coulomb field (B2.4) `b2-eq-EL`, the dipole-energy form (B2.12) `b2-eq-dipole`.
Forward links for the cards (all exported ids): `b4-eq-H`, `b4-eq-coulomb`, `b4-eq-A`, `b5-eq-Pie`, `b5-eq-Piion`,
`b5-eq-jion`, `b7-eq-piph`, `b8-eq-vertex`, `b9-eq-We`, `b3-eq-eps`; Part A: `a2-eq-H0`, `a2-eq-H`, `a2-eq-g`, `a5-eq-W`,
`a4-eq-glr`, `a4-eq-dd`, `a4-eq-lot`, `a3-eq-scf`.

**What to draw.** Keep the v1 SVG layout (viewBox 380×316, nodes as rounded rectangles, edges as clickable groups with a
halo, FS font scaling on narrow screens).
- Nodes:
  - electrons (left, electron colour): title "electrons", sub-label EPP view "**r**_s, **v**_s", Part A view "ψ_{n**k**}";
  - ions (right, phonon colour): title "ions · phonons", sub-label EPP view "dipole eZ^*·**u**", Part A view "u_{ℓκα}";
  - transverse field (top, photon colour, solid border): title "transverse field **A**^T", sub "photons";
  - longitudinal field (bottom, photon colour, dashed border): title "longitudinal field **A**^L", sub "= Coulomb".
  Bold vectors with rich text `{bold:true}` (or `font-weight` on the vector tspans), upright labels T, L, en, sr.
- EPP edges (one per term part, as in v1; wavy photon line = transverse, dashed photon line = longitudinal):
  - `et` term 5 T, wavy, electrons→**A**^T, label "−(e/c) **v**·**A**^T";
  - `it` term 6 T, wavy, ions→**A**^T, label "(e/c)(Z^*·**u̇**)·**A**^T";
  - `el` term 5 L, dashed, electrons→**A**^L, label "−(e/c) **v**·**A**^L";
  - `il` term 6 L, dashed, ions→**A**^L, label "(e/c)(Z^*·**u̇**)·**A**^L";
  - `ei` ghost dotted line electrons–ions, italic label "no direct term";
  - self-loops: `e0` term 1 on electrons, label "V^{en}_0"; `i0` terms 2 and 3 on ions (zigzag loop), label "Φ^{sr}";
    `t0` term 4 T on **A**^T (wavy loop), label "ħcq"; `l0` term 4 L on **A**^L (plain loop), label "(E^L)^2/8π".
- Part A edges: `g_g` solid vertex line electrons–ions with a dot, label "g_{mnν}(**k**,**q**)" and small "screened, DFPT";
  `g_e0` loop on electrons "E_{n**k**}"; `g_v` dashed loop near the bottom "v δ(t−t′)"; `g_i0` zigzag loop on ions
  "ħω_{**q**ν}". Placeholders where the field nodes were: "no dynamical field / (light: external drive only)" at the top,
  "Coulomb: instantaneous v, / folded into E_{n**k**}, ω_{**q**ν} and g" at the bottom (italic, muted).
- Legend strip under the SVG (v1 HTML legend, converted): wavy = "couples to the transverse field $\bA^{\mathrm T}$",
  dashed = "longitudinal field $\bA^{\mathrm L}$ (Coulomb)", solid with dot = "explicit vertex $g$ (Part A)",
  circle = "free term of one subsystem".
- Info panel (`role="status"`, `aria-live="polite"`) with one card per edge plus two intro cards; cards in HTML with
  MathJax. Rewrite every card in v2 notation, with eqrefs instead of equation numbers. Content of the cards:
  - `intro-epp`: "Click, tap or Tab to a line. Each line is one term of the Lagrangian (B2.6); the panel gives the term and
    what it becomes after quantization (Section B.4). The dotted line in the middle: no term couples electrons and ions
    directly."
  - `intro-g`: "Part A's world. The transverse field is gone, the longitudinal field has collapsed into the instantaneous
    Coulomb interaction $v$, and an explicit vertex $g$ (eqref `a2-eq-g`) joins electrons and phonons."
  - `e0`: term 1, $\sum_s[\frac12m_{\mathrm e}\mathbf v_s^2-V^{\mathrm{en}}_0(\br_s)]$ with (B2.7): electrons in the frozen
    crystal potential; the ions' potentials do not follow the ions, so there is no short-range (deformation-potential)
    coupling; quantized: band electrons $E_{n\bk}$.
  - `i0`: terms 2 and 3: harmonic ions on short-range springs $\Phi^{\mathrm{sr}}$; the dipole–dipole part $\Phi^{\mathrm{dd}}$
    (`a4-eq-dd`) is left out because the field regenerates it (else counted twice); quantized: bare phonons.
  - `t0`: term 4, transverse part, $\frac1{8\pi}\int[(\bE^{\mathrm T})^2-B^2]$, written as
    $\frac1{8\pi}\int[(\partial_t\bA^{\mathrm T})^2/c^2-(\nabla\times\bA^{\mathrm T})^2]$ by (B2.8): the curl costs magnetic
    energy, a spring; each transverse wave is an oscillator of frequency $cq$, a photon (field expansion `b4-eq-A`).
  - `l0`: term 4, longitudinal part, $\frac1{8\pi}\int(\bE^{\mathrm L})^2$: no curl, no spring; Gauss's law fixes
    $\bE^{\mathrm L}$ at each instant (B2.4), and this energy becomes the instantaneous Coulomb energy
    $\frac12\iint\rho(\br)\rho(\br')/|\br-\br'|$ of electrons and ionic bound charges, $\rho=\rho^{\mathrm e}+\rho^{\mathrm{ion}}$,
    $\rho^{\mathrm{ion}}=-\nabla\cdot\bP$ (`b4-eq-coulomb`).
  - `et`: term 5, transverse part: in the Hamiltonian, electron–photon couplings linear and quadratic in $\bA$ (`b4-eq-H`);
    optical absorption and emission; the electronic photon self-energy $\Pi^{\mathrm e}$ (`b5-eq-Pie`).
  - `el`: term 5, longitudinal part: electrons are sources of $\bE^{\mathrm L}$; two such couplings make the
    electron–electron Coulomb interaction; screened by the electrons it is $W^{\mathrm e}$ (`b9-eq-We`).
  - `it`: term 6, transverse part: the polarization current of a TO phonon emits and absorbs photons; ion–photon couplings
    linear and quadratic in $\bA$ (`b4-eq-H`); infrared absorption at $\omega_{\mathrm{TO}}$, phonon polaritons, the ionic
    photon self-energy $\Pi^{\mathrm{ion}}$ (`b5-eq-Piion`) with the ionic current vertex (`b5-eq-jion`).
  - `il`: term 6, longitudinal part: a longitudinal polarization wave carries bound charge $-\nabla\cdot\bP$; its Coulomb
    energy with itself is the dipole–dipole interaction that lifts LO above TO (Sections B.3, B.7; `a4-eq-lot`); its Coulomb
    energy with the electrons is the unscreened polar (Fröhlich-type) coupling.
  - `ei`: no direct term: electron and ion meet only in $\bA$, by the Coulomb path (electron → $\bA^{\mathrm L}$ → bound
    charge, instantaneous) and the light path (electron current → $\bA^{\mathrm T}$ → polarization current, finite $c$); the
    vertex returns as a composite, electron current → photon dressed by the electrons → ionic current (`b8-eq-vertex`).
  - `g_e0`: Part A: band energies $E_{n\bk}$ of the working Hamiltonian (`a2-eq-H0`); they already contain the mean-field
    (Hartree and exchange–correlation) part of the electron–electron interaction. EPP counterpart: term 1 plus the
    Coulomb energy of $\bA^{\mathrm L}$.
  - `g_v`: Part A: the instantaneous Coulomb interaction $v$, the only interaction of the many-body Hamiltonian of
    Section A.5; screened it is $W$ (`a5-eq-W`). EPP counterpart: the energy of the longitudinal field.
  - `g_i0`: Part A: phonons $\hbar\omega_{\bq\nu}$ from density-functional perturbation theory; their force constants contain
    the dipole–dipole term, so the LO–TO splitting is already in $\omega_{\bq\nu}$ (`a4-eq-dd`). EPP removes that part from
    term 3 and regenerates it through $\bA^{\mathrm L}$.
  - `g_g`: Part A: $\hat H^{\mathrm{e\text{-}ph}}=N_\ell^{-1/2}\sum g_{mn\nu}(\bk,\bq)c^\dagger_{m\bk+\bq}c_{n\bk}\phi_{\bq\nu}$
    (`a2-eq-H`); $g$ is an input, the matrix element of the screened potential change (`a3-eq-scf`), whose polar long-range
    part is `a4-eq-glr`. It already contains the screened Coulomb field of the displaced ions, the physics EPP sends along
    the field paths.

**Controls (2).**
- `seg` "View": `EPP Lagrangian` / `Part A Hamiltonian` (initial: EPP).
- `toggle` "electrostatic limit, c → ∞" (EPP view only; disabled or hidden in the Part A view): fades the transverse node
  and the three transverse edges (`t0`, `et`, `it`) to 0.2 opacity and makes them non-interactive, leaving the Coulomb lines.
  Its card (`intro-static`): "With $c\to\infty$ only the longitudinal field remains: every interaction is instantaneous
  Coulomb, which, screened by the electrons, is the world of Part A. The transverse lines matter only near the light cone,
  $q\lesssim\sqrt{\varepsilon_\infty}\,\omega_{\mathrm{TO}}/c=5.6\times10^{-5}$ Å⁻¹ in GaAs."
- Plus the clickable edges themselves (Tab order, Enter/Space selects, Escape clears; `role="button"`, `aria-label` per
  edge in words, e.g. "Term 5, transverse part: electron current coupled to the transverse field").

**Readouts.** None besides the info panel.

**Initial state.** EPP view, nothing selected, card `intro-epp`; the dotted "no direct term" line visible in the middle.

**Caption key points.**
- *What to notice.* Every line from electrons or ions ends on the field; nothing joins them directly. In the Part A view the
  field disappears into an instantaneous Coulomb line and an explicit vertex $g$, whose physics EPP carries along the two
  field paths.
- *Marks.* Node colours: electrons (electron colour), ions (phonon colour), field (photon colour; dashed border = longitudinal).
  Wavy line: coupling to the transverse field; dashed: longitudinal (Coulomb); loops: free term of one subsystem; dotted
  grey: no direct term; solid line with a dot (Part A): explicit vertex. Labels are the terms of (B2.6) with their numbers.
- *Try.* (1) Click the dotted middle line, then the two paths around it. (2) Switch on "c → ∞": only Coulomb lines remain.
  (3) Switch to the Part A view and click $g$; compare with the EPP paths.
- *Model.* Structure of (B2.6) with the split (B2.3); Coulomb energy of the longitudinal field from Section B.4; no
  numbers except the light-cone estimate for GaAs ($\varepsilon_\infty=10.9$, $\hbar\omega_{\mathrm{TO}}=33.3$ meV,
  literature values).

**What in v1 must change.**
- Source notation in SVG labels: `r_e, v_e` → **r**_s, **v**_s; `p = Z^*Q` → dipole eZ^*·**u**; `Δτ_{καp}` (Giustino's
  displacement) → u_{ℓκα}; `V_{cr}` → V^{en}_0; `k` → Φ^{sr}; `ε_{nk}` → E_{n**k**}; `A_t`, `A_l` → **A**^T, **A**^L;
  `E_l^2/8π` → (E^L)^2/8π; edge labels `(e/c) v·A_t` → `−(e/c) v·A^T` (the electron charge is −e) and `(1/c) ṗ·A_t` →
  `(e/c)(Z^*·u̇)·A^T`.
- Every card cites source equations (E.134, E.136, E.214, E.235′, E.236, E.237′, E.238, G Eq. 1, G Eqs. 62–67, G Eq. 38,
  G Eq. 180) and uses source symbols ($\Psi^\dagger\Psi$, $Z^*_{a\alpha\mu}\pi_{ia\mu}$, $\epsilon$, $g^{\rm dyn}=i\sum J^*D^{\rm e}J^{\rm e}$):
  rewrite in v2 notation with eqrefs as listed above. The factor-$e$ remark about the printed Coulomb energy belongs to B.4,
  not here: drop it.
- Seg labels `EPP Lagrangian (E.134)` / `Giustino H (G Eq. 1)` → `EPP Lagrangian` / `Part A Hamiltonian`; "Giustino's
  world" → "Part A".
- v1 had one control; add the c → ∞ toggle.
- v1 caption is one paragraph with "What to notice" in it; rewrite in the v2 structure. Tag `interactive schematic` →
  `interactive`.
- Section links `§B.4` etc. → `<a href="#b-fields">Section B.4</a>`.

---

## B2.2 · `fig-lorentz.html` · `viz-b2-lorentz` · new

**Title.** "An electron moved by the vector potential alone" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. In the temporal gauge both forces on an electron come from $\bA$: the electric force from $\partial_t\bA$ at the electron,
   the magnetic force from how $\bA$ varies along its path. No scalar potential is needed.
2. The canonical momentum $\mathbf p_s=m_{\mathrm e}\mathbf v_s-\frac ec\bA(\br_s,t)$ (B2.15) differs from the kinetic momentum
   $m_{\mathrm e}\mathbf v_s$.
3. A component of $\mathbf p_s$ is conserved along any direction on which $\bA$ (and the Lagrangian) does not depend: with
   $B=0$ all of $\mathbf p_s$ is conserved, so $m_{\mathrm e}\mathbf v_s$ simply follows $\frac ec\bA(t)$; with $B\neq0$ in the
   gauge used here, $p_y$ stays conserved while $p_x$ changes.

**Include position.** In `h4#b2-lorentz`, after the first paragraph and before the derivation steps:
> … Figure B2.2 shows a free electron in a vector potential that produces a uniform oscillating electric field and a static
> magnetic field $B$, which bends the path into circles at the cyclotron frequency $\omega_{\mathrm c}=eB/m_{\mathrm e}c$
> ($\hbar\omega_{\mathrm c}=0.116$ meV per tesla). The arrows compare the kinetic momentum $m_{\mathrm e}\mathbf v_s$ with the
> canonical momentum of the Lagrangian, which differs from it by $-\frac ec\bA$.

**Physics.** Free electron ($V^{\mathrm{en}}_0=0$), mass $m_{\mathrm e}$, charge $-e$, starting at rest at the origin at $t=0$.
Vector potential (temporal gauge, Gaussian):
$$\bA(\br,t)=-\frac{cE_0}{\omega}\sin(\omega t)\,\hat{\mathbf x}+Bx\,\hat{\mathbf y},$$
so by (B2.2) $\bE=E_0\cos(\omega t)\,\hat{\mathbf x}$ (uniform) and $\bB=B\hat{\mathbf z}$ (uniform, out of the plane).
Equation of motion (B2.16) `b2-eq-lorentz` with $V^{\mathrm{en}}_0=0$. Numerics in SI ($\bA_{\rm SI}=\bA/c$):
- $\omega_{\mathrm c}=eB/m_{\mathrm e}$; $\dot v_x=-(e/m_{\mathrm e})E_0\cos\omega t-\omega_{\mathrm c}v_y$, $\dot v_y=+\omega_{\mathrm c}v_x$;
- $p_x=m_{\mathrm e}v_x+(eE_0/\omega)\sin\omega t$, $p_y=m_{\mathrm e}v_y-eBx$ (SI form of (B2.15) `b2-eq-ps`).
- Integrate with RK4 over a fixed window of 6 drive periods ($T=2\pi/\omega$), about 4000–20000 steps; recompute the whole
  trajectory only when a slider changes (well under 15 ms), then animate a marker along it. Optional closed form:
  $w=v_x+\ii v_y$ obeys $\dot w=-(e/m_{\mathrm e})E_0\cos\omega t+\ii\omega_{\mathrm c}w$.
- Momentum unit: $eE_0/\omega$ (the amplitude of $m_{\mathrm e}v_x$ at $B=0$). Length unit: nm.
- Constants: $e=1.602177\times10^{-19}$ C, $m_{\mathrm e}=9.109384\times10^{-31}$ kg, $\hbar=1.054572\times10^{-34}$ J s;
  $\hbar\omega_{\mathrm c}=0.11577$ meV per tesla.
- Fixed field amplitude $E_0=1$ kV/cm (the motion is linear in $E_0$; state it in the caption).
- Reference values (author's RK4 check, 6 periods, $E_0=1$ kV/cm):
  - $B=0$, $\hbar\omega=0.5$ meV: $x$ from $-61.0$ to 0 nm, $y=0$; $|m_{\mathrm e}\mathbf v|_{\max}=1.00$; $p_x\equiv0$.
  - $B=2$ T, $\hbar\omega=0.5$ meV (initial state): $\hbar\omega_{\mathrm c}=0.2315$ meV, $T=8.27$ ps; $x\in[-73.2,76.7]$ nm,
    $y\in[-56.5,42.1]$ nm; $|m_{\mathrm e}\mathbf v|_{\max}=1.86$; $p_x\in[-0.64,0.86]$; $p_y\equiv0$ to $10^{-14}$.
  - $B=4.32$ T, $\hbar\omega=0.5$ meV (resonance $\omega=\omega_{\mathrm c}$): orbit grows to about 550 nm; $|m_{\mathrm e}\mathbf v|_{\max}\approx19$.
  - $B=2$ T, $\hbar\omega=2$ meV: excursions about 4 nm.
  - Peak speeds stay below 1% of $c$: at most $v/c=7\times10^{-3}$ (2.1×10⁶ m/s), at cyclotron resonance with the lowest drive,
    $\hbar\omega=0.1$ meV and $B\approx0.85$ T. The nonrelativistic Lagrangian holds.

**What to draw.**
- Panel (a), a square canvas (`EPP.canvas`, aspect ≈ 1, max height about 360 px): the $x$–$y$ plane.
  - Faint grid; when $B\neq0$, a sparse pattern of small ⊙ marks (field out of the page) labelled "**B**" (bold).
  - The full trajectory over the window as a thin line (electron colour, 40% alpha) and the part already travelled as a
    full-strength line; the electron as a dot at the current position.
  - Three arrows from the electron, with one common momentum scale (scale bar "momentum $eE_0/\omega$" in a corner):
    kinetic momentum $m_{\mathrm e}$**v** (electron colour, solid), the field term $-\frac ec$**A**$(\br_s,t)$ (photon colour,
    solid, thinner), and their sum, the canonical momentum **p** (fg colour, dashed). Labels at the arrow tips with
    `drawRich`: "m_e **v**", "−(e/c)**A**", "**p**" (bold vectors via `{bold:true}`; "e" in m_e upright).
  - An inset arrow at the top left for the field $\bE(t)$ (photon colour), length ∝ $\cos\omega t$, label "**E**(t)".
  - Length scale bar in nm. The view is fitted to the whole trajectory (with 10% margin) each time a slider changes and stays
    fixed during the animation (no jumping).
- Panel (b), an `EPP.Plot` below (aspect ≈ 0.4): momentum components against time $t/T$ (0–6), in units of $eE_0/\omega$:
  $m_{\mathrm e}v_x$ (electron colour, solid), $p_x$ (electron colour, dashed), $m_{\mathrm e}v_y$ (`c4`, solid), $p_y$ (`c4`,
  dashed); a vertical line at the current time. Fixed $y$ limits per parameter set (symmetric, fitted to the window).
  Axis labels "t / T (drive periods)" and "momentum (units of eE_0/ω)". Legend with the four series.

**Controls (3).**
- range "magnetic field $B$": 0–5 T, step 0.05, initial 2 T, format "x.xx T".
- range "drive $\hbar\omega$": 0.1–2 meV, log, initial 0.5 meV, format "x.xx meV".
- Play/Pause from `EPP.anim` (calm: one drive period in about 2 s of animation; starts paused under reduced motion). The
  slider change restarts the marker at $t=0$.

**Readouts.** $\hbar\omega_{\mathrm c}$ (meV); $\omega/\omega_{\mathrm c}$ (or "∞" at $B=0$); $t$ (ps); at the marker:
$m_{\mathrm e}\mathbf v$ and $\mathbf p$ as $(x,y)$ components in units of $eE_0/\omega$; "conserved: $p_x$, $p_y$" when
$B=0$, "conserved: $p_y$" otherwise; peak $v/c$.

**Initial state.** $B=2$ T, $\hbar\omega=0.5$ meV, paused at $t\approx1.3\,T$ (so that the three arrows are clearly
different), with the trajectory drawn: a looping orbit, $p_y$ flat at zero, $m_{\mathrm e}v_y$ oscillating.

**Caption key points.**
- *What to notice.* The kinetic momentum $m_{\mathrm e}\mathbf v$ and the canonical momentum $\mathbf p=m_{\mathrm e}\mathbf v-\frac ec\bA$
  differ. $p_y$ never changes, because $\bA$ does not depend on $y$; with $B=0$, $p_x$ is constant too, and the velocity just
  follows $\frac ec\bA(t)$.
- *Marks.* Arrows: kinetic momentum (electron colour), the term $-\frac ec\bA$ at the electron (photon colour), canonical
  momentum (dashed). ⊙: magnetic field out of the page. Top-left arrow: the uniform field $\bE(t)=-\frac1c\partial_t\bA$.
  In (b) solid lines are kinetic and dashed lines canonical components.
- *Try.* (1) Set $B=0$: the path is a straight segment along $x$, $p_x$ and $p_y$ are flat, and $m_{\mathrm e}v_x$ oscillates as
  $-\sin\omega t$. (2) Set $B=4.3$ T at $\hbar\omega=0.5$ meV, where $\omega=\omega_{\mathrm c}$: the orbit spirals outward (cyclotron
  resonance), and $p_y$ still stays at zero. (3) Raise $\hbar\omega$ to 2 meV: the electron barely moves.
- *Model.* (B2.16) with $V^{\mathrm{en}}_0=0$ and the vector potential above, which by (B2.2) gives a uniform $\bE=E_0\cos\omega t\,\hat{\mathbf x}$
  and $\bB=B\hat{\mathbf z}$; canonical momentum (B2.15); free electron of mass $m_{\mathrm e}$, $E_0=1$ kV/cm, start at rest at
  the origin; RK4 over six periods. The magnetic part of $\bA$ is written as $Bx\,\hat{\mathbf y}$ (one of many gauge choices
  with $A_0=0$), so the Lagrangian does not depend on $y$.

---

## B2.3 · `fig-drive.html` · `viz-b2-drive` · port of `viz-b-lagrangian-drive`

**Title.** "Driven ions and electrons: amplitude and phase" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. The relative displacement of a polar ion pair is a driven damped oscillator (B2.18)–(B2.19): in phase below resonance,
   90° behind at resonance, in antiphase above.
2. At resonance the polarization current $\dot P_x$ is in phase with $E$, so the field does work on the lattice (infrared
   absorption).
3. At phonon frequencies the bound electrons follow the field with no delay (the constant $\varepsilon_\infty$), while free
   carriers lag by 90°–180°.

**Include position.** End of `h4#b2-drive`, after the paragraph:
> The field that drives the ions in (B2.18) is the field at the atoms … The figure sets $\omega_0$ to the measured transverse
> optical frequency of GaAs and reads $E$ as the macroscopic field, which is what the photon-exchange treatment of Sections
> B.5–B.7 amounts to.

and the paragraph before it, which introduces the three species ("Figure B2.3 compares the ion pair with two other charges
driven by the same field …").

**Physics (keep the v1 code; it is correct).** For each species, the complex displacement amplitude of a charge $q$ (in
units of $e$, sign irrelevant for the dipole) of mass $m$:
$$x=\frac{qeE_0/m}{\omega_0^2-\omega^2-\ii\gamma\omega},$$
which for the ion pair is (B2.19) `b2-eq-xi` with $q=Z^*$, $m=M_{\mathrm r}$ (from (B2.18) `b2-eq-rel`). The induced dipole
of one particle is $qe\,x$ for any sign of $q$, its response $\propto q^2$; the phase plotted is that of the dipole (and of
the polarization) behind $E$: $\varphi=\arg(\omega_0^2-\omega^2+\ii\gamma\omega)\in(0°,180°)$. For an electron (charge $-e$)
the displacement itself points opposite to its dipole; say so in the caption.
Polarization along the field $P_x=n\,qe\,x$ with $n$ the number density; current $j_x=\dot P_x=-\ii\omega P_x$; contribution to the
dielectric function $4\pi P_x/E$ (Gaussian; equal to $P_x/\varepsilon_0E$ in SI); time-averaged absorbed power $\frac12\omega|P_x|E_0\sin\varphi$.
(The field is $\bE=E\hat{\mathbf x}$; write $P_x$, $j_x$ for the components, never a scalar $P$.)
Parameters (literature values, illustrative):
- ion pair (GaAs-like): $Z^*=2.07$, $M_{\mathrm r}=36.11$ amu, $\OV=45.17$ Å³ ($n=1/\OV$), $\hbar\omega_0$ slider (initial
  33.3 meV $=\hbar\omega_{\mathrm{TO}}$), $\hbar\gamma$ slider (initial 0.31 meV); Ga–As bond 2.448 Å; zero-point amplitude of
  the relative motion $(\hbar/2M_{\mathrm r}\omega_0)^{1/2}=0.0417$ Å at 33.3 meV;
- bound electron (classical Lorentz model, not part of (B2.6)): $q=1$, $m=m_{\mathrm e}$, $\hbar\omega_0=4.97$ eV, $\hbar\gamma=0.5$ eV,
  $n=8/\OV$; chosen so that $1+4\pi ne^2/(m_{\mathrm e}\omega_0^2)=10.9=\varepsilon_\infty$ of GaAs (check:
  $\hbar^2\cdot4\pi ne^2/m_{\mathrm e}=244.2$ eV², so $\hbar\omega_0=\sqrt{244.2/9.9}=4.966$ eV);
- free carrier (Drude): $q=1$, $m=0.067\,m_{\mathrm e}$, $\omega_0=0$, $\hbar\gamma=3.5$ meV, $n=10^{17}$ cm⁻³.
- Units in code (v1, keep): $\hbar^2/(\mathrm{amu}\,\text{Å}^2)=4.1802$ meV, $\hbar^2/m_{\mathrm e}=7.61996$ eV Å²; field in kV/cm,
  $e\times1$ kV/cm $=10^{-5}$ eV/Å; displacement in Å.
Reference values at the initial state ($E_0=100$ kV/cm, drive at 33.3 meV): ion pair $|x|=0.0232$ Å, lag 90.0°;
static ion response $eZ^*E_0/K=2.16\times10^{-4}$ Å with $K=M_{\mathrm r}\omega_0^2=9.58$ eV/Å².

**What to draw (keep the v1 layout).**
- (a) `EPP.Plot` amplitude $|x|$ (Å) against $\hbar\omega$ (meV), log–log, $x$ 1 meV–20 eV, $y$ $10^{-7}$–$10^3$ Å, three curves
  (ion pair: phonon colour; bound electron: electron colour; free carrier: `c5`, dashed), hlines "Ga–As bond length" and
  "zero-point amplitude of the ion pair", vline at the drive frequency, direct curve labels.
- (b) `EPP.Plot` phase lag $\varphi$ (deg, 0–180, ticks every 45°) against $\hbar\omega$, same $x$ axis; hline at 90° with
  the text "$\dot P_x$ in phase with $E$"; texts "follows $E$" near 0° and "opposes $E$" near 180°; direct labels.
- (c) `EPP.Plot` time traces over two drive periods for the selected species: $E(t)/E_0$ (photon colour),
  $P_x(t)/|P_x|_{\max}$ (species colour), $j_x(t)/|j_x|_{\max}$ with $j_x=\dot P_x$ (`c4`, dashed).
- Clicking or dragging on (a) or (b) sets the drive frequency.

**Controls (5, as v1, relabelled).**
- range "short-range resonance $\hbar\omega_0$" (ion pair): 10–60 meV, step 0.1, initial 33.3 meV.
- range "ion damping $\hbar\gamma$": 0.05–5 meV, log, initial 0.31 meV.
- range "field amplitude $E_0$": 1–1000 kV/cm, log, initial 100 kV/cm.
- range "drive $\hbar\omega$": 1 meV–20 eV, log, initial 33.3 meV.
- seg "time trace for": ion pair / bound electron / free carrier (initial: ion pair).

**Readouts.** $\hbar\omega$; $|x|$ of the selected species (Å); "$P_x$ lags $E$ by" $\varphi$; "$j_x=\dot P_x$ leads $E$ by"
$90°-\varphi$; contribution to $\varepsilon$, $4\pi P_x/E$ (complex, Gaussian = SI $P_x/\varepsilon_0E$); $|j|$ (A/cm²); absorbed
power $\langle jE\rangle$ (W/cm³); a warning "|u_rel| > 10% of the bond: beyond the harmonic model" for the ion pair. State
in the caption that the practical readouts are in SI units.

**Initial state.** As v1: drive at 33.3 meV (on resonance), ion-pair trace, $E_0=100$ kV/cm. The 90° lag and the in-phase
current are visible immediately.

**Caption key points.**
- *What to notice.* At the ionic resonance the displacement lags the field by exactly 90°, so the polarization current is in
  phase with $E$ and power is absorbed. At phonon frequencies the bound electrons respond with no lag: that instantaneous
  response is $\varepsilon_\infty$.
- *Marks.* Colours per species; vertical dashed line: drive frequency; horizontal lines in (a): bond length and zero-point
  amplitude; in (c) the field, the polarization and its current.
- *Try.* (1) Move the drive across 33.3 meV and watch the phase jump from 0° to 180° within $\sim\hbar\gamma$. (2) Lower
  $\hbar\gamma$: the peak grows as $\omega_0/\gamma$. (3) Select the free carrier: its lag is 135° at $\hbar\omega=3.5$ meV, tends to
  90° (resistive) only for $\omega\ll\gamma$, below the plotted range, and to 180° (inductive) for $\omega\gg\gamma$.
- *Model.* (B2.18)–(B2.19) for the ion pair with GaAs-like values (above); the same driven-oscillator form for a harmonically
  bound valence electron (model) and a Drude carrier; damping added by hand (the harmonic Lagrangian has none). $E$ is the
  field acting on the charges; for the ions it is read as the macroscopic field, with $\omega_0$ set to the TO frequency (no
  local-field correction, as in the photon-exchange route of Sections B.5–B.7; Section B.3 treats the local field).

**What in v1 must change.**
- Labels: "bare optical ħω_O" → "short-range resonance $\hbar\omega_0$"; legend "ion pair (E.143′)" → "ion pair" (no source
  equation numbers anywhere); readout key "4πχ" → "$4\pi P_x/E$" ($\chi$ is reserved for density responses); "p lags E by" →
  "$P_x$ lags $E$ by" ($\mathbf p$ is a canonical momentum in v2); "dP/dt" → "$j_x=\dot P_x$"; time-trace labels "p(t)/|p|" →
  "P_x(t)/|P_x|_max".
- Caption: $\Omega_{\rm uc}$ (EPP's cell volume) → $\OV$; "$Z^*=2.07|e|$" → "$Z^*=2.07$ (charge $eZ^*$)"; "$M=36.1$ u" →
  "$M_{\mathrm r}=36.1$ amu"; $\epsilon_\infty$ → $\varepsilon_\infty$; "for an electron, $e<0$" → the electron has charge
  $-e$ with $e>0$; "4πχ_Gauss = χ_SI" → "$4\pi P_x/E$ (Gaussian) $=P_x/\varepsilon_0E$ (SI)"; section links §B.5–§B.7 → v2 links;
  "What to notice" moves from the end to `p.cap-key`; restructure into Marks / Try / Model.
- Rename the id and classes (`b-lag-amp` → `b2-drive-amp`, etc.).
- Physics: unchanged (formulas verified against (B2.19); bound-electron $\omega_0$ reproduces $\varepsilon_\infty=10.9$).
