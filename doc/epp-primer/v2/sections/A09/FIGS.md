# A.9 figure specifications

Section file: `v2/sections/A09-applications.html` (id `a-applications`, prefix `a9`).
Four figures, in reading order:

| # | file | id | kind |
|---|---|---|---|
| A9.1 | `A09/fig-indirect.html` | `viz-a9-indirect` | new (geometric band diagram) |
| A9.2 | `A09/fig-edge.html` | `viz-a9-edge` | port of v1 `viz-a-applications-edge` |
| A9.3 | `A09/fig-drift.html` | `viz-a9-drift` | new (geometric, k space) |
| A9.4 | `A09/fig-tc.html` | `viz-a9-tc` | port of v1 `viz-a-applications-tc` |

General rules for all four (BUILDER_GUIDE §5): v2 notation only (NOTATION.md); CSS classes prefixed
`a9-<name>-`; colours read at draw time with `EPP.color(...)` (`electron`, `phonon`, `photon`, `c5`,
`muted`, …); canvas text with `EPP.drawRich` (bold vectors with `{bold:true}`, upright labels:
`Γ`, `X`, `TA`, `TO`, `LO`, subscripts such as `E_{g}` rendered with an upright g where the label is a word);
HTML readouts with MathJax `$…$`; no `σ` anywhere in strings (it is reserved for the conductivity); fixed
axis limits; 2–5 controls with units; works at 390 px; keyboard access and `aria-label`s; caption of about
200–260 words in the structure *What to notice / Marks / Try / Model*, with `a.eqref` links to the section's
equations. Test with `bash v2/build/figtest.sh v2/sections/A09/fig-<name>.html`.

Equation ids of the section that the captions may link (tag in brackets):
`a9-eq-A` (A9.1, light wave), `a9-eq-hel` (A9.2, electron–light coupling), `a9-eq-alpha` (A9.3,
absorption coefficient), `a9-eq-amp` (A9.4, two-path amplitude), `a9-eq-indirect` (A9.5, phonon-assisted
absorption), `a9-eq-edge` (A9.6, parabolic edge), `a9-eq-wl` (A9.7, Williams–Lax), `a9-eq-drude` (A9.8),
`a9-eq-j` (A9.9, current), `a9-eq-bte` (A9.10, linearized Boltzmann equation), `a9-eq-mob` (A9.11, SERTA
mobility), `a9-eq-tautr` (A9.12, τ and τ_tr), `a9-eq-bcs` (A9.13), `a9-eq-fs` (A9.14, Fermi-surface average),
`a9-eq-a2fs` (A9.15, isotropic Eliashberg function), `a9-eq-lamfs` (A9.16, λ_ep), `a9-eq-wlog` (A9.17,
ω_log), `a9-eq-tc` (A9.18, McMillan–Allen–Dynes), `a9-eq-mustar` (A9.19, Morel–Anderson), `a9-eq-lamiw`
(A9.20), `a9-eq-screened` (A9.21). From other sections only: `a8-eq-varshni` (one-oscillator gap of A.8),
`a4-eq-frohlich`, `a7-eq-a2f`, `a7-eq-lambda`.

---

## 1. `fig-indirect.html` · `viz-a9-indirect` · new

**Title.** "Indirect absorption in silicon: two steps through a virtual state" `<span class="tag">interactive</span>`.

**What the reader must learn.**
- On a band diagram the photon arrow is vertical (it carries energy, ≈1.2 eV, and almost no wavevector) and the
  phonon arrow is almost horizontal (it carries the wavevector $\bq\approx\bk_{\mathrm c}-\bk_{\mathrm v}\approx1$ Å⁻¹
  and only 18–58 meV of energy).
- The transition can happen in two orders. In each the electron passes through a **virtual** state: a real band
  state reached with the wrong energy. The energy mismatch is exactly the denominator of (A9.4).
- Phonon emission ($s=+1$) needs $\hbar\omega\ge E_{\mathrm g}+\hbar\omega_{\bq\nu}$ and carries the weight
  $n_{\bq\nu}+1$; phonon absorption ($s=-1$) works already from $E_{\mathrm g}-\hbar\omega_{\bq\nu}$ but carries
  $n_{\bq\nu}$, which vanishes at low temperature.
- The photon's own wavevector, drawn to scale, is invisible on the scale of the zone.

**Include position.** In `#a9-light`, after the paragraph on the absorption coefficient and its `dl.where`,
right before `<h4 id="a9-steps">`. The sentences just before:
"To compute $\alpha_{\mathrm{abs}}$ one adds up $\hbar\omega$ times the rate of every transition, divides by the
volume of the crystal, and divides by $I$. [Figure A9.1] draws the bands of silicon from $\Gamma$ to X and the two
ways in which a phonon can complete a transition that light alone cannot make."
(The derivation of (A9.4)–(A9.5) follows the figure.)

**Physics (v2 notation).**
- Schematic silicon-like bands along $\Gamma$–X, $k\in[0,k_X]$ with $k_X=2\pi/a=1.157$ Å⁻¹ ($a=5.431$ Å),
  energies in eV measured from the valence-band maximum (literature values, illustrative):
  - valence band $E_{\mathrm v}(k)=-2.9\,\mathrm{eV}\times\sin^2(\pi k/2k_X)$: maximum 0 at $\Gamma$, $-2.9$ eV at X,
    $\approx-2.7$ eV at $k_{\mathrm c}$ (curvature near Γ corresponds to a hole mass of about 0.7 $m_{\mathrm e}$);
  - lowest conduction band: minimum 1.17 eV (the low-temperature gap of silicon; the figure has no temperature
    dependence, so this is the gap $E_{\mathrm g}$ that enters the δ function of (A9.5)) at $k_{\mathrm c}=0.85\,k_X=0.983$ Å⁻¹, value 3.4 eV at
    $\Gamma$ (the direct gap), $\approx1.30$ eV at X; near the minimum on the X side use the longitudinal mass
    0.92 $m_{\mathrm e}$, i.e. $E_{\mathrm c}=E_{\mathrm g}+3.81\,\mathrm{eV\,Å^2}\,(k-k_{\mathrm c})^2/0.92$; on the Γ side a smooth
    curve (for example $E_{\mathrm g}+(3.4-E_{\mathrm g})(1-k/k_{\mathrm c})^2$) joining 3.4 eV at Γ.
  The exact shapes are not important; the four reference values (0, 1.17 eV at $k_{\mathrm c}$, 3.4 eV at Γ, $-2.9$ eV
  at X) are.
- Assisting phonon: TO, $\hbar\omega_{\bq\nu}=57.8$ meV at $\bq\approx\bk_{\mathrm c}$ (mention TA, 18.4 meV, in the
  caption only). Bose occupation $n_{\bq\nu}=1/(\ee^{\hbar\omega_{\bq\nu}/\kB T}-1)$, $\kB=0.0861733$ meV/K.
- Initial state: valence state at $k_{\mathrm i}$ (default 0, the maximum). ($k_{\mathrm i}$, $k_{\mathrm f}$ are names used in this
  spec only; in the figure label the states $n\bk$ and $m\,\bk+\bq$ and the phonon wavevector $\bq$, as the section does.) Energy conservation, from the δ function of
  (A9.5): the final state lies on the conduction band at $k_{\mathrm f}$ with
  $E_{\mathrm c}(k_{\mathrm f})=E_{\mathrm v}(k_{\mathrm i})+\hbar\omega-s\hbar\omega_{\bq\nu}$; take the root of this equation on the Γ side of the
  minimum if it exists, else report "no final state". The phonon wavevector is $q=k_{\mathrm f}-k_{\mathrm i}$.
- Denominators of (A9.4) (`a9-eq-amp`), with the intermediate band chosen as the nearest band state at the
  intermediate wavevector:
  - photon first: intermediate state = conduction band at $k_{\mathrm i}$; denominator
    $D_1=E_{\mathrm c}(k_{\mathrm i})-E_{\mathrm v}(k_{\mathrm i})-\hbar\omega$ (≈ +2.2 eV at the default);
  - phonon first: intermediate state = valence band at $k_{\mathrm f}$; denominator
    $D_2=E_{\mathrm v}(k_{\mathrm f})-E_{\mathrm v}(k_{\mathrm i})+s\hbar\omega_{\bq\nu}$ (≈ −2.7 eV at the default).
  ($D_1$, $D_2$ are figure-local names; in the readouts write them out as "first denominator", "second
  denominator" with their formulas, not as new symbols.)
- Photon wavevector in the crystal $q_{\mathrm L}=n_{\mathrm r}\omega/c$ with $n_{\mathrm r}=3.5$ and
  $\hbar c=1973.27$ eV Å: $1.95\times10^{-3}$ Å⁻¹ at 1.1 eV, about $1/500$ of $k_{\mathrm c}$.

**What to draw.**
- Panel (a), main, a band diagram: horizontal axis wavevector $k$ along $\Gamma$–X (0 to 1.157 Å⁻¹, ticks every
  0.2 Å⁻¹, labels `Γ` and `X` upright at the ends, axis label "wavevector k along Γ–X (Å⁻¹)"); vertical axis energy
  (−3.2 to 4.0 eV, label "energy (eV)").
  - Valence band and conduction band as solid lines in `electron`; shade the gap lightly. A double arrow at
    $k_{\mathrm c}$ labelled $E_{\mathrm g}$ (indirect gap), another at Γ labelled "direct gap 3.4 eV".
  - Initial state: filled dot on the valence band at $k_{\mathrm i}$. Final state: open dot on the conduction band at
    $k_{\mathrm f}$.
  - Photon-first path: vertical arrow in `photon` from the initial state up by $\hbar\omega$ to a **dashed short
    horizontal level** (`muted`) = the virtual state (energy available, not on a band); a thin bracket from that level
    to the real conduction-band state at the same $k$, labelled "mismatch = first denominator"; then the phonon arrow
    in `phonon` (wavy or zigzag), almost horizontal, from the virtual level to the final state (it drops or rises by
    $\hbar\omega_{\bq\nu}$).
  - Phonon-first path: phonon arrow from the initial state, almost horizontal, to a dashed virtual level at
    $k_{\mathrm f}$ and energy $E_{\mathrm v}(k_{\mathrm i})-s\hbar\omega_{\bq\nu}$; bracket down to the real valence state at $k_{\mathrm f}$,
    labelled "mismatch = second denominator"; then a vertical photon arrow up to the final state.
  - The selected path is drawn fully; with "both", the second path is drawn at 40% opacity.
  - A bold vector label $\bq$ (drawRich bold) on the phonon arrow and $\hbar\omega$ on the photon arrow.
  - At the foot of the photon arrow, the photon wavevector to scale: a horizontal tick of length $q_{\mathrm L}$
    (less than a pixel) with the label "photon wavevector, to scale".
- Panel (b), small, below (a) on phones: the two time orders as diagrams, read left to right: an electron line
  $n\bk\to n'\to m\,\bk+\bq$ (`electron`, arrows), the light as a wavy `photon` line ending in a **cross** (a
  classical external field, not a propagator), the phonon as a zigzag `phonon` line leaving the second (photon
  first) or the first (phonon first) vertex. Labels on the intermediate line: $n'\bk$ (photon first) or
  $n'\,\bk+\bq$ (phonon first). Highlight the selected order. (This replaces the v1 Feynman-diagram block
  `fd-a-applications-hbb`, which v2 does not keep in the section text.)

**Controls (4).**
1. Photon energy $\hbar\omega$: range 1.00–1.40 eV, step 0.005, default 1.25 eV.
2. Process: segmented "phonon emission ($s=+1$)" / "phonon absorption ($s=-1$)", default emission.
3. Path shown: "photon first" / "phonon first" / "both", default both.
4. Temperature $T$: 4–600 K, default 296 K (affects only the weight readout).
Optional: drag the initial-state dot along the valence band, $k_{\mathrm i}\in[0,0.3]$ Å⁻¹ (keyboard: arrow keys).

**Readouts.** $\hbar\omega$ (eV); threshold $E_{\mathrm g}+s\hbar\omega_{\bq\nu}$ (eV); $|\bq|$ (Å⁻¹) or "no final
state below threshold"; $\hbar\omega_{\bq\nu}$ (meV); first and second denominators (eV); phonon weight
$n_{\bq\nu}+\tfrac12+\tfrac s2$ at $T$; ratio $q_{\mathrm L}/|\bq|$.

**Initial state.** $\hbar\omega=1.25$ eV, emission, both paths, $T=296$ K, $k_{\mathrm i}=0$. The final state then lies at
1.192 eV on the conduction band just on the Γ side of $k_{\mathrm c}$; denominators ≈ +2.15 eV and ≈ −2.7 eV;
weight $n+1=1.12$.

**Caption key points.**
- *What to notice.* Light provides the energy and the phonon the crystal momentum: the photon arrow is vertical,
  the phonon arrow nearly horizontal. Each order passes through a virtual state whose energy mismatch, 2 to 3 eV in
  silicon, is a denominator of (A9.4).
- *Marks.* Colours (electron bands, photon, phonon); dashed levels are virtual states, brackets the mismatches;
  the cross in (b) marks the classical light field; the tick at the foot of the photon arrow is the photon
  wavevector to scale.
- *Try.* Lower $\hbar\omega$ below $E_{\mathrm g}+\hbar\omega_{\bq\nu}$ with emission: no final state. Switch to absorption:
  final states return down to $E_{\mathrm g}-\hbar\omega_{\bq\nu}$, but cool to 4 K and watch the weight $n_{\bq\nu}$
  vanish. Switch paths and compare the denominators.
- *Model.* Schematic silicon bands (values above, literature, illustrative), TO phonon 57.8 meV (TA 18.4 meV also
  assists), $n_{\mathrm r}=3.5$; energy conservation and denominators from (A9.4)–(A9.5)
  (`a9-eq-amp`, `a9-eq-indirect`); the intermediate band is taken as the nearest one, whereas (A9.4) sums over all.

---

## 2. `fig-edge.html` · `viz-a9-edge` · port of v1 `viz-a-applications-edge`

v1 source: the `<figure class="viz" id="viz-a-applications-edge">` element and the "Viz A: indirect absorption
edge" block of the `<script>` at the end of `doc/epp-primer/sections/A09-applications.html`.

**Title.** "Phonon-assisted absorption edge of silicon: two slopes and a moving gap" `<span class="tag">interactive</span>`.

**What the reader must learn.**
- Each phonon branch gives two thresholds, $\hbar\omega_\nu$ below and above the gap; on a $\sqrt{\alpha_{\mathrm{abs}}}$ plot the
  absorption and emission parts are nearly straight lines with slopes $\propto\sqrt{n_\nu}$ and $\sqrt{n_\nu+1}$:
  the two slopes.
- Cooling freezes out the absorption part (weight $n_\nu$) and keeps the emission part ($n_\nu+1\to1$).
- With clamped-nuclei energies (the bare formula) the onset never moves; with the temperature-dependent gap of
  Section A.8 in the δ function it moves down on heating, the 296 K versus 78 K difference that the review's
  Fig. 14 had to impose by a rigid shift.

**Include position.** End of `#a9-parabolic`, after the "two slopes" paragraph, right before
`<h4 id="a9-moving">`. The sentence just before: "[Figure A9.2] evaluates (A9.6) with these two phonons and lets the
gap follow the temperature, as explained next."

**Physics (v2 notation).**
- (A9.6) (`a9-eq-edge`):
  $\alpha_{\mathrm{abs}}(\omega)=\frac{C}{\hbar\omega}\sum_\nu w_\nu\big[(n_\nu+1)(\Delta E_+)^2\Theta(\Delta E_+)+n_\nu(\Delta E_-)^2\Theta(\Delta E_-)\big]$,
  $\Delta E_\pm=\hbar\omega-E_{\mathrm g}\mp\hbar\omega_\nu$, $n_\nu=n_{\mathrm B}(\omega_\nu)$.
- Silicon-like parameters (literature values, illustrative): TO $\hbar\omega_\nu=57.8$ meV, $w=1$; TA 18.4 meV,
  $w=0.3$; **$C=6.0\times10^3$ cm⁻¹ eV⁻¹** (recalibrated, see below).
- Gap: with the toggle on, the one-oscillator form of A.8 (`a8-eq-varshni`):
  $E^{\mathrm{qp}}_{\mathrm g}(T)=E_{\mathrm g}-|\Delta E^{\mathrm{ZPR}}_{\mathrm g}|\,[1+2n_{\mathrm B}(\omega_0)]$ with
  the clamped-nuclei gap $E_{\mathrm g}=1.23$ eV, $|\Delta E^{\mathrm{ZPR}}_{\mathrm g}|=60$ meV (slider), $\hbar\omega_0=30$ meV, giving
  1.170 eV at 4 K, 1.169 eV at 78 K, 1.117 eV at 296 K, 1.017 eV at 600 K. With the toggle off, the gap stays at
  $E_{\mathrm g}=1.23$ eV at all $T$ (clamped nuclei, as in the bare (A9.5)), and the
  $|\Delta E^{\mathrm{ZPR}}_{\mathrm g}|$ slider has no effect (grey it out).
- Check values with $C=6.0\times10^3$ and the toggle on (computed): at 296 K, $\alpha_{\mathrm{abs}}=1.1$ cm⁻¹ at
  1.10 eV, 42 at 1.20 eV, $1.06\times10^3$ at 1.5 eV, $4.6\times10^3$ at 2.0 eV; at 78 K, 0.54 at 1.20 eV, 46 at
  1.30 eV, 435 at 1.5 eV, $2.5\times10^3$ at 2.0 eV. These match the review's Fig. 14 (theory and experiment)
  within about 20% at both temperatures. Occupations: $n_{\mathrm{TA}}=0.946$, $n_{\mathrm{TO}}=0.116$ at 296 K;
  0.069 and $1.8\times10^{-4}$ at 78 K; $n_{\mathrm{TA}}=7\times10^{-24}$ at 4 K.

**What to draw** (keep the v1 layout).
- Left plot: $\alpha_{\mathrm{abs}}$ on a log axis, $10^{-2}$ to $10^5$ cm⁻¹, against $\hbar\omega$ from 0.85 to 2.0 eV;
  the current temperature (solid, `photon`) and a fixed 78 K reference (dashed, `muted`); vertical dotted line at
  $E^{\mathrm{qp}}_{\mathrm g}(T)$ labelled `E^{qp}_{g}(T)` (upright qp and g). Axis labels "photon energy ħω (eV)", "α_{abs} (cm^{-1})".
- Right plot: $\sqrt{\alpha_{\mathrm{abs}}}$ (cm$^{-1/2}$), fixed limits 0–45, against $\hbar\omega$ from 0.88 to 1.36 eV; the
  total (solid `photon`), the phonon-absorption part ($s=-1$, dotted `c5`) and the phonon-emission part ($s=+1$,
  dashed `phonon`); thin dotted vertical lines at the gap minus $\hbar\omega_{\mathrm{TO,TA}}$ (`c5`, where absorption
  switches on) and the gap plus $\hbar\omega_{\mathrm{TA,TO}}$ (`phonon`, where emission switches on); dashed line at the
  gap. "The gap" is the one in the δ function: $E^{\mathrm{qp}}_{\mathrm g}(T)$ with the toggle on (label `E^{qp}_{g}(T)`),
  $E_{\mathrm g}$ with it off (label `E_{g}`); label the left plot's line the same way. (At 296 K the total reaches about 21 at 1.36 eV; at 600 K about 42.)

**Controls (4)** (as in v1).
1. Temperature $T$: 4–600 K, step 1, default 296 K.
2. Zero-point renormalization $|\Delta E^{\mathrm{ZPR}}_{\mathrm g}|$: 0–80 meV, step 1, default 60 meV.
3. Toggle "renormalized gap $E^{\mathrm{qp}}_{\mathrm g}(T)$ in the δ function", default on (off: the clamped gap $E_{\mathrm g}$).
4. Toggle "show emission and absorption parts", default on.

**Readouts.** $E^{\mathrm{qp}}_{\mathrm g}(T)$ (eV, 3 decimals); $n_{\mathrm{TO}}$, $n_{\mathrm{TA}}$ (show "&lt; 10⁻⁶" when tiny);
$\alpha_{\mathrm{abs}}$ at $\hbar\omega$ equal to the gap in the δ function (cm⁻¹); $\alpha_{\mathrm{abs}}$ at 1.5 eV (cm⁻¹).

**Initial state.** $T=296$ K, $|\Delta E^{\mathrm{ZPR}}_{\mathrm g}|=60$ meV, gap toggle on, parts on: the tail below the gap
and the kinks at the four thresholds are visible; the 78 K reference lies lower and starts later.

**Caption key points.**
- *What to notice.* Cooling removes the absorption below the gap, because that part is weighted by $n_\nu$;
  at 4 K it is negligible ($n_{\mathrm{TA}}\sim10^{-23}$). With the clamped gap the onset never moves; with $E^{\mathrm{qp}}_{\mathrm g}(T)$
  it moves down on heating.
- *Marks.* Colours and dash patterns of the two plots; violet (`c5`) and amber (`phonon`) dotted lines are the
  absorption and emission thresholds.
- *Try.* Lower $T$ to 78 K and 4 K and watch the absorption part vanish while the emission slope stays; switch
  the gap toggle off: the onset jumps to the clamped-nuclei value, about 60 meV too high; set
  $|\Delta E^{\mathrm{ZPR}}_{\mathrm g}|$ to 0.
- *Model.* (A9.6) (`a9-eq-edge`), parameters above, $E^{\mathrm{qp}}_{\mathrm g}(T)$ from `a8-eq-varshni`; $C$ chosen so that
  $\alpha_{\mathrm{abs}}(1.5\ \mathrm{eV},296\ \mathrm K)\approx10^3$ cm⁻¹ as in the review's Fig. 14; valid only near the
  edge (parabolic bands, constant matrix elements).

**What in the v1 figure is wrong or in source notation, and must change.**
- **Calibration:** v1 uses $C=1.2\times10^4$ cm⁻¹ eV⁻¹, which gives $2.1\times10^3$ cm⁻¹ at 1.5 eV and $9\times10^3$ at
  2.0 eV (296 K), twice the review's Fig. 14 (about $10^3$ and $4\times10^3$). Use $C=6.0\times10^3$; the v1 right-plot
  limit (0–40) then becomes 0–45 to fit 600 K.
- Notation, aligned with Section A.8: `α` → $\alpha_{\mathrm{abs}}$; `E_g^{cn}` (clamped) → $E_{\mathrm g}$; `E_g(T)`
  (renormalized) → $E^{\mathrm{qp}}_{\mathrm g}(T)$; `|ΔE_g^{ZP}|` → $|\Delta E^{\mathrm{ZPR}}_{\mathrm g}|$; `ħω_eff` → $\hbar\omega_0$;
  `n_{TO}`, `n_{TA}` keep upright labels. All labels upright where they abbreviate words (g, qp, ZPR, TO, TA).
- References: "(G.200)", "G Eq. 168", "Hall–Bardeen–Blatt" as an equation name, "§A.8" → `a.eqref` links to
  `a9-eq-edge`, `a8-eq-varshni` and the section link to `#a-temperature`; "Giustino's Fig. 14" → "the review's Fig. 14".
- Caption: v1 puts *Model* first and *What to notice* last in one paragraph; rewrite in the v2 structure.
- Ids and classes: `viz-a-applications-edge`, `.plotA`, `.plotB`, `.ctrls`, `.ro` → `viz-a9-edge` and
  `a9-edge-…` classes; readout labels in MathJax.

---

## 3. `fig-drift.html` · `viz-a9-drift` · new

**Title.** "Drift in k space: which collisions destroy the current" `<span class="tag">interactive</span>`.

**What the reader must learn.**
- In a steady state the field displaces the occupied region of $\bk$ space against $\bE$ (for electrons,
  $\hbar\dot\bk=-e\bE$); collisions pull it back.
- SERTA counts every collision as a full loss of the drift ($1/\tau$); the current actually decays with the
  momentum relaxation rate $1/\tau_{\mathrm{tr}}$, which weights each collision by $1-\cos\theta$ (A9.12).
- When small-angle collisions dominate (polar, Fröhlich-like coupling), $\tau_{\mathrm{tr}}\gg\tau$: SERTA
  underestimates the mobility and the iterative solution is needed. When scattering is isotropic (silicon-like),
  $\tau_{\mathrm{tr}}=\tau$ and SERTA is exact.

**Include position.** End of `#a9-drude`, right before `<h4 id="a9-boltzmann">`. The sentences just before:
"…And does every collision destroy the drift? A collision that turns the velocity by a small angle leaves most of
it along the field. The Boltzmann equation answers both questions. [Figure A9.3] shows the steady state in
$\bk$ space and how much of the drift a single collision removes."
(The figure is referred to again in `#a9-serta`, where $q_0$ is introduced: "[Figure A9.3] uses a squared
matrix element proportional to $1/(q^2+q_0^2)$ …".)

**Physics (v2 notation).** A two-dimensional parabolic band (drawn in 2D for clarity; the section says so),
elastic scattering on the energy shell $|\bk'|=|\bk|=k$.
- Rate per final angle $R(\theta)\propto1/(|\bk'-\bk|^2+q_0^2)$ with $|\bk'-\bk|^2=2k^2(1-\cos\theta)$, i.e. a
  squared matrix element $\propto1/(q^2+q_0^2)$ with $\bq=\bk'-\bk$: a Fröhlich-like $1/q^2$ (Section A.4) cut off at
  $q_0$. Normalize so that $1/\tau=\int_0^{2\pi}R\,\dd\theta$ is the same for all $q_0$ (only the angular shape
  changes).
- $1/\tau=\int_0^{2\pi}R\,\dd\theta$, $1/\tau_{\mathrm{tr}}=\int_0^{2\pi}R\,(1-\cos\theta)\,\dd\theta$ (2D form of
  (A9.12), `a9-eq-tautr`). Closed form with $x=q_0/2k$:
  $\tau_{\mathrm{tr}}/\tau=(\sqrt{1+x^2}+x)/(2x)$; values: $x=0.05\to10.5$, $0.1\to5.52$, $0.3\to2.24$, $1\to1.21$,
  $3\to1.03$, $10\to1.00$. (In 3D the ratio is smaller, e.g. 2.4 at $x=0.1$, same trend; mention in *Model*.)
- For elastic scattering in an isotropic band the iterative solution of (A9.10) is exactly
  $\delta f\propto\tau_{\mathrm{tr}}\,\mathbf v\cdot\bE$, so $\mu_{\mathrm{IBTE}}/\mu_{\mathrm{SERTA}}=\tau_{\mathrm{tr}}/\tau$
  (write these as "mobility, full solution / mobility, SERTA" in the readout, not as new symbols).
- Displacement of the occupied disk: $\delta\bk=-e\bE\tau_{\mathrm{tr}}/\hbar$ (full solution) versus
  $-e\bE\tau/\hbar$ (SERTA), drawn in arbitrary units; scale so that the SERTA shift at the largest field is
  $0.04\,k$ and cap the drawn full shift at $0.5\,k$ (show "×N" beside the arrow when capped).

**What to draw.**
- Panel (a), square, k space: axes $k_x$, $k_y$ (labels in italic with subscripts), origin marked. The
  equilibrium Fermi disk of radius $k$ (dashed outline, `muted`); the displaced disk of the full solution (filled,
  light `electron`, solid outline) and the SERTA displacement as a second dashed outline; small "+" and "−"
  marks on the crescents where occupation is gained and lost. The field arrow $\bE$ (bold label) at a corner,
  pointing in $-x$, and the acceleration arrow of the electrons in $+x$.
- One selected electron on the disk edge at $\bk$ (dot, bold label $\bk$; default on the $+x$ side). A fan of
  ~24 thin arrows to final states $\bk'$ around the circle, each with opacity and width $\propto R(\theta)$
  (`phonon`). One highlighted arrow at the angle set by the slider, labelled with the bold $\bq=\bk'-\bk$ and the
  angle $\theta$ (arc).
- Panel (b), below on phones: $R(\theta)$ (normalized, `phonon`, filled lightly) and $R(\theta)(1-\cos\theta)$
  (`electron`, filled lightly) against $\theta$ from −180° to 180°; legend "all collisions (SERTA, area ∝ 1/τ)" and
  "weighted by 1 − cos θ (area ∝ 1/τ_{tr})"; a vertical marker at the highlighted $\theta$. Fixed $y$ limits
  (normalize $R$ to its maximum over the slider range, or plot on a fixed 0–1 scale relative to $R(0)$).

**Controls (4).**
1. Small-angle preference $q_0/2k$: log slider 0.05–10, default 0.1 (aria-valuetext with the value).
2. Highlighted collision angle $\theta$: −180° to 180°, step 1°, default 30° (also by dragging the final state on
   the circle).
3. Field strength: 0–1 (arbitrary units), default 0.6.
4. Preset buttons: "polar, GaAs-like" ($q_0/2k=0.05$) and "isotropic, silicon-like" ($q_0/2k=10$).

**Readouts.** $\tau_{\mathrm{tr}}/\tau$ (= mobility of the full solution over SERTA); $1-\cos\theta$ of the
highlighted collision; $|\bq|/k$ of the highlighted collision; share of collisions with $|\theta|<30°$.

**Initial state.** $q_0/2k=0.1$ ($\tau_{\mathrm{tr}}/\tau=5.5$), $\theta=30°$ ($1-\cos\theta=0.13$), field 0.6: the fan
is strongly forward-peaked and the full-solution disk is visibly displaced much more than the SERTA outline.

**Caption key points.**
- *What to notice.* With small-angle scattering most collisions turn the velocity only a little: they shorten the
  lifetime but hardly relax the current, so $\tau_{\mathrm{tr}}\gg\tau$ and SERTA underestimates the mobility.
- *Marks.* Dashed circle: equilibrium; filled disk: steady state of the full solution; dashed outline: SERTA
  steady state; fan arrows: collisions from $\bk$, opacity ∝ rate; (b) areas ∝ $1/\tau$ and $1/\tau_{\mathrm{tr}}$.
- *Try.* Press "isotropic": the two areas become equal and the disks coincide. Drag the highlighted collision
  from 0° to 180°: the weight $1-\cos\theta$ grows from 0 to 2. Raise the field: both shifts grow in proportion.
- *Model.* 2D parabolic band, elastic scattering, squared matrix element $\propto1/(q^2+q_0^2)$ (a screened
  Coulomb or Fröhlich-like form, Section A.4); (A9.12) (`a9-eq-tautr`), (A9.10) (`a9-eq-bte`), (A9.11)
  (`a9-eq-mob`); not a materials calculation (real LO scattering in GaAs is also inelastic); in 3D the ratio
  is smaller for the same $q_0/2k$.

---

## 4. `fig-tc.html` · `viz-a9-tc` · port of v1 `viz-a-applications-tc`

v1 source: the `<figure class="viz" id="viz-a-applications-tc">` element and the "Viz B: alpha^2F -> lambda,
omega_log, Tc" block of the `<script>` at the end of `doc/epp-primer/sections/A09-applications.html`.

**Title.** "From the Eliashberg function to $T_{\mathrm c}$" `<span class="tag">interactive</span>` (write the title
with MathJax or as "T<sub>c</sub>" with an upright c).

**What the reader must learn.**
- $\lambda_{\mathrm{ep}}$ comes mostly from low-frequency weight (the $1/\omega$ in (A9.16)); $\omega_{\log}$ sets the
  prefactor of (A9.18).
- $T_{\mathrm c}$ depends exponentially on $\mu^*$ and vanishes below $\lambda_{\mathrm{ep}}=\mu^*/(1-0.62\mu^*)$.
- Scaling $|g|^2$ by 1.4, the size of the GW and hybrid-functional corrections of the last subsection, multiplies
  $T_{\mathrm c}$ by 4.5 in weakly coupled Al but only by 1.3 in Pb.
- The isotropic formula underestimates two-gap MgB₂ and strongly coupled H₃S.

**Include position.** In `#a9-mcmillan`, after the callout "Physical picture: why $\mu^*$ is small", right before
`<h4 id="a9-eliashberg">`. The sentence just before: "[Figure A9.4] computes $\lambda_{\mathrm{ep}}$, $\omega_{\log}$ and
$T_{\mathrm c}$ from model Eliashberg functions of four superconductors."
(The figure is also cited in the paragraph of numbers after (A9.18).)

**Physics (v2 notation).**
- Model Eliashberg function: $\alpha^2F(\omega)=s_g\sum_i\frac{\lambda_i\omega_i}{2}\,G(\omega-\omega_i;\Delta\omega_i)$,
  normalized Gaussians of centre $\omega_i$ and root-mean-square width $\Delta\omega_i$, each contributing about
  $\lambda_i$ to $\lambda_{\mathrm{ep}}$; $s_g$ is the $|g|^2$ scale. (Name the scale "$|g|^2$ scale" in the UI; avoid a new
  symbol in the caption if possible. Do **not** call the width σ.)
- $\lambda_{\mathrm{ep}}=2\int_0^\infty\alpha^2F/\omega\,\dd\omega$ (A9.16, `a9-eq-lamfs`), $\omega_{\log}$ (A9.17,
  `a9-eq-wlog`), the running integral $2\int_0^{\omega}\alpha^2F(\omega')\,\dd\omega'/\omega'$ (call it "running
  $\lambda_{\mathrm{ep}}$ up to ω"; **not** $\lambda(\omega)$, because $\lambda(\ii\omega)$ is the Matsubara coupling of (A9.20)),
  $T_{\mathrm c}$ from (A9.18) (`a9-eq-tc`) with $\kB=0.0861733$ meV/K. Trapezoid rule on 800 points from 0 to the
  preset maximum, as in v1.
- Optional strong-coupling factors of Allen and Dynes (1975), not part of (A9.18), as in v1:
  $\Lambda_1=2.46(1+3.8\mu^*)$, $\Lambda_2=1.82(1+6.3\mu^*)\,\bar\omega_2/\omega_{\log}$,
  factor $[1+(\lambda/\Lambda_1)^{3/2}]^{1/3}\,[1+(\bar\omega_2/\omega_{\log}-1)\lambda^2/(\lambda^2+\Lambda_2^2)]$ with
  $\bar\omega_2=\langle\omega^2\rangle^{1/2}$ (keep these as code; in strings call them "strong-coupling factors";
  internal variable names only, no symbols in the caption).
- Presets: peaks $(\hbar\omega_i\ \mathrm{meV},\lambda_i,\hbar\Delta\omega_i\ \mathrm{meV})$, literature-inspired, illustrative:
  Al (20, 0.25, 3), (37, 0.18, 2.5), $x$-range 0–50 meV, measured 1.18 K;
  Pb (4.4, 0.9, 0.6), (8.5, 0.65, 1.1), 0–12 meV, 7.2 K;
  MgB₂ isotropic (65, 0.6, 5), (35, 0.1, 6), (90, 0.03, 5), 0–110 meV, 39 K;
  H₃S (30, 0.5, 6), (120, 1.0, 12), (170, 0.6, 10), 0–210 meV, about 200 K (at about 155 GPa).
- Check values (computed with the v1 code): Al $\lambda_{\mathrm{ep}}=0.437$, $\hbar\omega_{\log}=25.2$ meV (292 K),
  $T_{\mathrm c}=1.96/1.27/0.99/0.56$ K at $\mu^*=0.10/0.12/0.13/0.15$; Pb 1.579, 5.6 meV (65 K), 7.81/7.38/7.16/6.73 K;
  MgB₂ 0.737, 59.5 meV (691 K), 27.1/23.2/21.3/17.7 K; H₃S 2.136, 91.8 meV (1066 K), 160/153/150/143 K, with the
  strong-coupling factors 197 K ($\mu^*=0.10$) and 180 K (0.13). $|g|^2\times1.4$ at $\mu^*=0.12$ multiplies
  $T_{\mathrm c}$ by 4.47 (Al), 1.30 (Pb), 1.98 (MgB₂), 1.21 (H₃S). Threshold $\mu^*/(1-0.62\mu^*)=0.130$ at 0.12.

**What to draw** (keep the v1 layout).
- Left: $\alpha^2F(\omega)$ (filled, `phonon`) and the running $\lambda_{\mathrm{ep}}$ (dashed, `electron`) against
  $\hbar\omega$ (meV), $x$ limits per preset, fixed $y$ limit per preset (room for $|g|^2$ up to ×2, as in v1); a
  vertical line at $\hbar\omega_{\log}$ labelled `ω_{log}`.
- Right: $T_{\mathrm c}$ (K, log axis $10^{-2}$–$10^3$) against $\mu^*$ (0–0.25): the current curve (`electron`), a
  point at the chosen $\mu^*$, the $|g|^2\times1$ reference (dashed `muted`) whenever the scale differs from 1 or the
  strong-coupling factors are on, a band for the usual range 0.10–0.20, a horizontal dashed line at the measured
  $T_{\mathrm c}$ (`c5`, "measured (lit.)").

**Controls (4)** (as in v1).
1. Material: Al / Pb / MgB₂ (isotropic) / H₃S, default Al.
2. $\mu^*$: 0.05–0.25, step 0.005, default 0.12.
3. $|g|^2$ scale (for example DFT → GW): ×0.5 to ×2, step 0.05, default ×1.
4. Toggle "strong-coupling factors (Allen–Dynes, not in (A9.18))", default off.

**Readouts.** $\lambda_{\mathrm{ep}}$; $\hbar\omega_{\log}$ in meV and in K; $T_{\mathrm c}$ (K); threshold
$\mu^*/(1-0.62\mu^*)$; $1+\lambda_{\mathrm{ep}}$.

**Initial state.** Al, $\mu^*=0.12$, ×1, factors off: $\lambda_{\mathrm{ep}}=0.44$, $\hbar\omega_{\log}=25$ meV,
$T_{\mathrm c}=1.3$ K against 1.18 K measured.

**Caption key points.**
- *What to notice.* On the log axis $T_{\mathrm c}$ falls steeply with $\mu^*$ and collapses at the threshold.
  Raising $|g|^2$ by 40% multiplies $T_{\mathrm c}$ by about 4.5 for Al and 2 for MgB₂ at $\mu^*=0.12$, but only by
  1.3 for Pb: weak coupling is the most sensitive.
- *Marks.* Filled curve $\alpha^2F$; dashed running $\lambda_{\mathrm{ep}}$ (its final value is $\lambda_{\mathrm{ep}}$);
  vertical line $\omega_{\log}$; right: curve, point, reference, measured line, usual band.
- *Try.* Switch to MgB₂: 18–27 K for $\mu^*=0.15$–0.10 against 39 K, the isotropic formula cannot describe its two
  gaps; switch to H₃S (about 150 K against 200 K) and turn on the strong-coupling factors (180–197 K); scale
  $|g|^2$.
- *Model.* Gaussian peaks (presets above, literature-inspired, illustrative); (A9.16)–(A9.18)
  (`a9-eq-lamfs`, `a9-eq-wlog`, `a9-eq-tc`); strong-coupling factors of Allen and Dynes (1975) optional.

**What in the v1 figure is wrong or in source notation, and must change.**
- The Gaussian width is called `σ_i` in the v1 caption: σ is reserved for the conductivity; use
  $\Delta\omega_i$ (or "width") and keep σ out of all strings (the linter flags it).
- `λ` → $\lambda_{\mathrm{ep}}$; `T_c` → $T_{\mathrm c}$; `λ(ω) = 2∫α^2F/ω` (cumulative) → "running
  $\lambda_{\mathrm{ep}}$ up to ω" (the symbol $\lambda(\ii\omega)$ means the Matsubara coupling in v2); `λ_{crit}` →
  "threshold $\mu^*/(1-0.62\mu^*)$"; `f_1f_2` and "(not in G.204)" → "strong-coupling factors (Allen–Dynes, not in
  (A9.18))" ($f$ denotes occupations in v2).
- References "(G.204)", "(G.206)", "Eq. (G.204)", "Sec. XII" → `a.eqref` links to `a9-eq-tc`, `a9-eq-lamfs`,
  `a9-eq-wlog` and the section link to `#a9-dft`.
- Caption: v1 puts *Model* first and *What to notice* last; rewrite in the v2 structure.
- Ids and classes: `viz-a-applications-tc`, `.plotA`, `.plotB`, `.ctrls`, `.ro` → `viz-a9-tc`, `a9-tc-…`;
  readout labels in MathJax; legend strings in drawRich (`α^2F(ω)`, `ω_{log}`, `T_{c}` with upright c, `μ^*`).
- The v1 physics and numbers are correct (checked above) and stay.
