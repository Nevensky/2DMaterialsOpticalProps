# B.7 figures: specifications for the figure authors

Section file: `v2/sections/B07-phonon-se.html` (id `b-phonon-se`, prefix `b7`). Read it before building a figure: every
symbol, equation tag and number below is defined there. Figures appear in this order, so the page numbers them
B7.1–B7.3: `exchange`, `resum`, `limits`.

| File | id | Kind | v1 source |
|---|---|---|---|
| `fig-exchange.html` | `viz-b7-exchange` | new (diagram with a stepper) | replaces the inline `div.fd` diagrams of v1 (`b-phonon-se-fd-ion`, `-fd-de`, `-fd-f`, `-fd-om`, drawn by the section script of `../sections/B07-phonon-se.html`) |
| `fig-resum.html` | `viz-b7-resum` | port | view "Why resum?" of `viz-b-phonon-se-limits` in `../sections/B07-phonon-se.html` (function `drawResum`) |
| `fig-limits.html` | `viz-b7-limits` | port | view "Order of limits" of `viz-b-phonon-se-limits` in `../sections/B07-phonon-se.html` (functions `model`, `trans`, `drawLimits`) |

## Common rules (BUILDER_GUIDE §5, NOTATION)

- v2 notation only, in canvas text (`EPP.drawRich`, `EPP.fd.draw` labels), HTML readouts (MathJax `$…$`), legends and
  aria-labels. No source symbols: no `σ` (branch index or anything else: `σ` is the conductivity and the linter flags it in
  script strings), no `Ω^e`, `Ω^0`, `Ω_uc`, `Π_jj`, `Π^{ion,0}`, `m_a`, `Z^*_a` with charge, `ε` written `\epsilon`, no
  `λ_c` (λ is the photon polarization index in v2), no "E Eq. …" labels.
- Branch index ν; Cartesian α, β; labels upright as whole words: `ω_{TO}`, `ω_{LO}`, `ω_{UP}`, `ω_{LP}`, `ε_∞`, `q_c`
  (upright c), `Π^{ph}`, `Π^{ph,0}`, `Π^{ion}`, `Π^{ion,scr}`, `D^e`, `Z^{ph}`, `ω_p^{ion}`. Bold vectors with
  `{bold:true}` (**q**, **j**).
- The bare optical frequency of the model equals $\omega_{\mathrm{TO}}$ (macroscopic fields only, no local field), so
  write `ω_{TO}` where v1 wrote `ω_σ`.
- Colours at draw time with `EPP.color(...)`. Roles in this section: phonon lines and the LO branch `phonon`; photon
  lines and the light line `photon`; electron bubbles `electron`; vertices $J$, $J^*$ `c4`; upper polariton `c5`; lower
  polariton `c4`; references `muted`; the electron-free value `ink` dotted.
- Units on every axis; fixed axis limits (see each figure); 2–5 controls with units; no NaN; < 15 ms per redraw; works
  at 390 px with controls below the plots and no horizontal overflow; keyboard access and aria-labels; animations only
  through `EPP.anim`.
- Captions: "What to notice" (1–2 sentences, visible), then `cap-more` with Marks, Try, Model; about 200–260 words in
  total. Reference equations with `<a class="eqref" href="#b7-eq-…"></a>`. Credit "the EPP manuscript" without figure
  or equation numbers (source numbers belong in the section's footnotes).
- Material data (literature values, illustrative; the same as the section text and Section B.8):

| Preset | $Z^*$ | $\OV$ (Å³) | $M_{\mathrm r}$ (amu) | $\varepsilon_\infty$ | $\hbar\omega_{\mathrm{TO}}$ (meV) | measured $\hbar\omega_{\mathrm{LO}}$ (meV) | $a$ (Å) |
|---|---|---|---|---|---|---|---|
| GaAs | 2.07 | 45.16 | 36.11 | 10.9 | 33.3 | 36.2 | 5.653 |
| NaCl | 1.12 | 44.85 | 13.95 | 2.34 | 20.3 | 32.7 | 5.64 |
| hBN, in plane (treated as isotropic) | 2.71 | 18.08 | 6.10 | 4.87 | 169.9 | 199.6 | 2.504 |

- Constants: $4\pi e^2=180.95$ eV Å; $\hbar^2/(\mathrm{amu}\,\text{Å}^2)=4.1802$ meV; $\hbar c=1973.27$ eV Å $=1.97327\times10^{6}$ meV Å.
  The squared ionic plasma frequency is $(\hbar\omega^{\ion}_{\mathrm p})^2=4\pi e^2Z^{*2}\hbar^2/(\OV M_{\mathrm r})$:
  in meV², `wp2 = 180.95 * Z*Z / (Om * Mr) * 4.1802e3`.
- Derived numbers the figures must reproduce (computed with `/tmp/claude-0/-home-user-2DMaterialsOpticalProps/e798d873-fe28-512b-a4f7-e6f2d29e481a/scratchpad/B07/b7_numbers.py`):

| Preset | $\hbar\omega^{\ion}_{\mathrm p}$ | $\hbar\omega_{\mathrm{LO}}$ (B7.17) | electron-free $q\to0$ (B7.13) | static $\Pi^{\ph,0}$ (B7.12) | $q_{\mathrm c}$ (B7.22) | $q_{\mathrm c}/(2\pi/a)$ |
|---|---|---|---|---|---|---|
| GaAs | 44.58 meV | 35.93 meV | 55.64 meV | 29.84 meV | 5.57×10⁻⁵ Å⁻¹ | 5.0×10⁻⁵ |
| NaCl | 38.95 meV | 32.56 meV | 43.92 meV | 37.36 meV | 1.57×10⁻⁵ Å⁻¹ | 1.4×10⁻⁵ |
| hBN ∥ | 224.4 meV | 198.0 meV | 281.5 meV | 148.2 meV | 1.90×10⁻⁴ Å⁻¹ | 7.6×10⁻⁵ |

  GaAs polaritons (electrons on, real $c$): at $q=q_{\mathrm c}$, $\hbar\omega_{\mathrm{UP}}=40.73$, $\hbar\omega_{\mathrm{LP}}=27.23$ meV,
  $Z^{\ph}_{\mathrm{UP}}=0.599$; at $3q_{\mathrm c}$, 100.92 and 32.96 meV, $Z^{\ph}_{\mathrm{UP}}=0.023$; at $10q_{\mathrm c}$, 333.3 and
  33.272 meV (TO 28 µeV below $\omega_{\mathrm{TO}}$), $Z^{\ph}_{\mathrm{LP}}=0.998$.

---

## 1. `fig-exchange.html` · id `viz-b7-exchange` · NEW (diagram with a stepper)

**Title.** "Photon exchange dresses the phonon" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. The phonon self-energy is one dressed photon exchanged by the lattice: phonon $\nu'$ makes a current ($J$), the current
   makes a field (longitudinal and transverse) that the electrons screen ($D^{\rme}$), and the field drives phonon $\nu$
   ($J^*$): $\Pi^{\ph}_{\nu\nu'}=-\sum J^*_{\alpha\nu}D^{\rme}_{\alpha\beta}J_{\beta\nu'}$.
2. The three steps of the section, as moves on diagrams: (1) chains of the full photon propagator regroup into ionic bubbles
   joined by electron-dressed lines, because every stretch of electron bubbles sums to $D^{\rme}$; (2) each ionic bubble is
   minus a phonon line between two current vertices; (3) the chain becomes vertex × (phonon series) × vertex, and the phonon
   series is the dressed propagator $F$.
3. The sign bookkeeping: the term with $N$ photon lines has $N+1$ bubbles and $N$ brackets, $(-1)^{N+1}(-1)^N=-1$.

**Include position.** In h3 `b7-picture`, directly before h3 `b7-step1`. The sentences just before it:
"A term of a series is a chain of lines and bubbles, read as a product of matrices. Figure B7.1 draws the three steps."
The figure therefore comes before all the formalism it draws.

**Physics** (all drawn, nothing computed): (B7.1) `b7-eq-regroup`, (B7.2) `b7-eq-scr`, (B7.4) `b7-eq-Piion`,
(B7.5) `b7-eq-piph`, (B7.6) `b7-eq-series`, (B7.7) `b7-eq-dysonF`; the inset Dyson equation of $D^{\rme}$ is
`b5-eq-dysonDe` (Section B.5).

**What to draw.** Two parts, stacked (they stack naturally on phones).
- (a) **The loop, read physically** (one SVG, about 340 × 120 px): a thin amber zigzag (phonon $F^0$, label `F^0`) enters
  from the left, a teal dot labelled `J` (below it, small text "current"), a heavy rose wavy line (label `D^e`, small text
  "field, longitudinal + transverse, screened by electrons"; draw three or four small blue ellipses along it, the
  polarized electrons), a teal dot `J^*` (small text "drives the lattice"), a thin amber zigzag leaving to the right. Under
  the middle, "−Π^{ph}_{νν′}" with a brace spanning the dots and the wavy line. Branch labels `ν′` (left line) and `ν`
  (right line). Static is fine; an optional calm `EPP.anim` pulse travelling left to right along the loop is welcome if it
  honours reduced motion.
- (b) **The three steps** as three rows of `EPP.fd.draw` diagrams (each row a flex-wrap `div` so terms wrap on phones; no
  term wider than 340 px). Use `phonon`/`phonon2` (thin/heavy zigzag), `photon`/`photon2` (thin/heavy wavy), `blob`,
  `vertex`, `text` items.
  - Row 1, "Step 1 · resum the electrons": `D` (heavy wavy, ink) `=` heavy `D^e` `+` `D^e`–amber blob `Π^{ion}`–`D^e` `+`
    `D^e`–`Π^{ion}`–`D^e`–`Π^{ion}`–`D^e` `+ …`. Beneath it, smaller, the inset `D^e = D^0 + D^0`–blue blob `Π^e`–`D^e`.
  - Row 2, "Step 2 · open the ionic bubble": amber blob `Π^{ion}` `= −` [teal dot `J`–thin zigzag `F^0`–teal dot `J^*`].
    Put the minus sign in front of the bracket, not on a line.
  - Row 3, "Step 3 · a Dyson equation for phonons": `Π^{ion,scr} = −J` (heavy zigzag `F`) `J^*`, then
    `F = F^0 + F^0`[dot `J^*`, heavy wavy `D^e`, dot `J`]`F^0 + F^0 … F^0 … F^0 + …`, and finally the boxed
    `F = F^0 + F^0 Π^{ph} F` with a heavy zigzag `F`. Mark each `J^*`–`D^e`–`J` piece with a light bracket labelled
    `−Π^{ph}`.

**Controls (3).**
- `step`: segmented `1 · electrons`, `2 · ionic bubble`, `3 · phonon series`, `all` (default `all`). The selected row is
  at full opacity with a thin accent outline; the others are dimmed to opacity 0.35. Keyboard-operable.
- `N`, number of dressed photon lines: range 0–3, step 1 (default 1). Highlights the term with $N$ lines $D^{\rme}$ in row 1
  and the corresponding term $F^0(\Pi^{\ph}F^0)^N$ in row 3 (accent halo), so the reader sees they are the same term.
- `electrons`: toggle "show the electron bubbles inside $D^{\rme}$" (default off). When on, the first heavy `D^e` of row 1
  is redrawn as `D^0 + D^0 Π^e D^0 + D^0 Π^e D^0 Π^e D^0 + …` in a small sub-row, to show that each stretch is a sum over
  electron bubbles.

**Readouts** (HTML, MathJax): the highlighted term, e.g. for $N=2$: "term $\Pi^{\ion}D^{\rme}\Pi^{\ion}D^{\rme}\Pi^{\ion}$
$=-J\,F^0\Pi^{\ph}F^0\Pi^{\ph}F^0\,J^*$"; the sign count "3 bubbles: $(-1)^3$; 2 brackets $J^*D^{\rme}J=-\Pi^{\ph}$:
$(-1)^2$; total $-1$". Use `aria-live="polite"`.

**Initial state.** `step = all`, `N = 1`, electrons off: the reader sees all three rows and one exchanged photon
highlighted in rows 1 and 3.

**Caption key points.**
- *What to notice.* The diagrams that dress light with matter (Section B.6) also dress the lattice with light: grouping the
  electron bubbles into $D^{\rme}$ and opening each ionic bubble turns every chain into a phonon line interrupted by photon
  exchanges, and each exchange $J^*D^{\rme}J$ is minus the phonon self-energy.
- *Marks.* Zigzag: phonon (thin $F^0$, heavy $F$). Wavy: field (thin $D^0$, heavy $D^{\rme}$). Amber blob: ionic bubble
  $\Pi^{\ion}$; blue blob: electronic bubble $\Pi^{\rme}$. Teal dots: vertices $J$ (the phonon makes a current) and $J^*$
  (the field drives the phonon). Brackets: $-\Pi^{\ph}$.
- *Try.* Step through 1, 2, 3. Move $N$ and watch the same term light up in rows 1 and 3, with the sign count. Show the
  electrons inside $D^{\rme}$.
- *Model.* <a class="eqref" href="#b7-eq-regroup"></a>, <a class="eqref" href="#b7-eq-scr"></a>,
  <a class="eqref" href="#b7-eq-Piion"></a>, <a class="eqref" href="#b7-eq-piph"></a>, <a class="eqref" href="#b7-eq-series"></a>,
  <a class="eqref" href="#b7-eq-dysonF"></a>; redrawn after the diagrams of the EPP manuscript, with the electron bubbles
  in the random-phase approximation and only macroscopic fields.

**v1 to v2.** v1 drew these diagrams inline (`div.fd`) from the section script, with EPP notation and labels
`Π^{ion,0}` (bare), `Π^{ion}` (screened), `Ω^e` and `J_{ασ}`, `J_{βσ'}`. In v2 the irreducible bubble is `Π^{ion}`, the
screened correlator `Π^{ion,scr}`, the self-energy `Π^{ph}`, branches `ν`, `ν′`; the section has no script, so all
drawing lives in this figure.

---

## 2. `fig-resum.html` · id `viz-b7-resum` · PORT (view "Why resum?" of v1 `viz-b-phonon-se-limits`)

**Title.** "Why the whole series is needed: truncated series against the Dyson sum" `<span class="tag">interactive</span>`.

**What the reader must learn.** No finite truncation of the phonon series moves the peak away from the bare frequency;
higher orders only pile up sharper structure there. The Dyson sum moves it, and the dressed pole sits exactly on the edge
of the frequency window where the series diverges, $|F^0\Pi^{\ph}|=1$. Without electrons the self-energy is comparable to
the phonon energy and the series diverges almost everywhere below the true pole.

**Include position.** At the end of h4 `b7-resum`, before h3 `b7-vacuum`. The sentences just before it: "Without electrons
the static self-energy of a GaAs optical mode at $\bq\to0$ is 29.8 meV (computed below), comparable to the energy of the
mode itself, … the window then extends from zero frequency up to the dressed mode. Figure B7.2 shows the truncated series
piling up structure at the bare frequency while the full sum moves the peak."

**Physics.** One branch at $\bq\to0$, energies in meV, propagators in meV⁻¹:
- free propagator with a phenomenological damping $\gamma$ (NOTATION: damped oscillators $\omega^2-\omega_0^2+\ii\gamma\omega$):
  $\hbar F^0(\omega)=2\hbar\omega_{\mathrm{TO}}/[(\hbar\omega)^2-(\hbar\omega_{\mathrm{TO}})^2+\ii\,\hbar\gamma\,\hbar\omega]$, the retarded form of
  `b5-eq-F0` with $\eta\to\gamma$;
- constant static self-energy $\Pi=(\hbar\omega^{\ion}_{\mathrm p})^2/(2\hbar\omega_{\mathrm{TO}}\varepsilon)$ with $\varepsilon=\varepsilon_\infty$
  (electrons on; (B7.20) `b7-eq-limq`) or $\varepsilon=1$ (electrons off; (B7.12) `b7-eq-pi0`);
- truncated series $F_N=F^0\sum_{k=0}^{N}(\Pi F^0)^k$ (the series (B7.6) `b7-eq-series` cut after $N$ exchanges);
- Dyson sum $F=F^0/(1-F^0\Pi)$ ((B7.9) `b7-eq-diag`), pole at $\hbar\bar\omega=[(\hbar\omega_{\mathrm{TO}})^2+2\hbar\omega_{\mathrm{TO}}\Pi]^{1/2}$ ((B7.10) `b7-eq-pole`);
- plotted quantity: $-\mathrm{Im}\,F/\pi$ in meV⁻¹ (do not call it $A$ or $\mathcal A$: $\bA$ is the vector potential and
  $\mathcal A$ the electron spectral function).
  Check values (GaAs): electrons on, $\Pi=2.74$ meV, Dyson pole 35.93 meV, divergence window 30.4–35.9 meV; electrons off,
  $\Pi=29.84$ meV, pole 55.64 meV, window from 0 to 55.6 meV.

**What to draw.** One plot, $-\mathrm{Im}\,F/\pi$ (meV⁻¹) against $\hbar\omega$ (meV):
- shaded band(s) where $|F^0\Pi|>1$ (`muted`, alpha 0.13), legend "series diverges, $|F^0\Pi^{\ph}|>1$";
- bare $F^0$ dashed `muted`; Dyson sum solid `phonon`; truncated $F_N$ solid `c5`, width 1.8, legend "truncated, $N=…$";
- vertical lines at the bare $\hbar\omega_{\mathrm{TO}}$ (`muted`, label "bare ω_{TO}") and at the Dyson pole (`phonon`, label
  "Dyson pole"); labels drawn last, on top of the spikes; values outside the frame clipped (`breakJump`).
- Fixed limits: x from $\max(0.3,\ \hbar\omega_{\mathrm{TO}}-1.6\Delta)$ to $\hbar\bar\omega+1.2\Delta$ with $\Delta=\hbar\bar\omega-\hbar\omega_{\mathrm{TO}}$,
  depending only on the preset and the electron toggle (v1 also moved it with $\gamma$; do not). y from $-0.7p$ to $1.3p$
  with $p=2\omega_{\mathrm{TO}}/(\pi\gamma\bar\omega)$, which changes only when $\gamma$ changes, never with $N$.

**Controls (4).** Crystal preset (seg: GaAs, NaCl, hBN ∥); electrons on/off (toggle, "electrons on ($\varepsilon_\infty$ in $D^{\rme}$)");
truncation order $N$ (range 0–20, step 1); broadening $\hbar\gamma$ (range 0.1–2 meV, step 0.05).

**Readouts.** bare pole $\hbar\omega_{\mathrm{TO}}$; Dyson pole $\hbar\bar\omega$ (B7.10); static $\Pi^{\ph}$ (meV); position of the maximum
of $F_N$ (meV); $|F^0\Pi^{\ph}|$ at the Dyson pole (shows 1 up to the broadening).

**Initial state.** GaAs, electrons on, $N=4$, $\hbar\gamma=0.5$ meV: the truncated curve peaks at 33.3 meV while the Dyson
peak sits at 35.9 meV, at the edge of the shaded window.

**Caption key points.**
- *What to notice.* However many exchanges are kept, the truncated series peaks at the bare frequency; only the sum moves
  the mode, to the edge of the window where the series diverges.
- *Marks.* Dashed: bare phonon. Amber: Dyson sum. Violet: truncated series $F_N$. Shading: $|F^0\Pi^{\ph}|>1$. Vertical
  lines: bare $\omega_{\mathrm{TO}}$ and dressed pole.
- *Try.* Raise $N$: the violet structure at $\omega_{\mathrm{TO}}$ sharpens and changes sign but never moves. Switch the
  electrons off: the window spans almost all frequencies below the pole (GaAs: 0 to 55.6 meV). Shrink $\hbar\gamma$.
- *Model.* <a class="eqref" href="#b7-eq-series"></a>, <a class="eqref" href="#b7-eq-diag"></a>,
  <a class="eqref" href="#b7-eq-pole"></a>; static self-energy <a class="eqref" href="#b7-eq-pi0"></a> (electrons off) or
  <a class="eqref" href="#b7-eq-limq"></a> (on); presets from the table; damping $\gamma$ added by hand (the theory has
  none in an insulator).

**v1 to v2.** Split out of the two-view v1 figure (keep its `drawResum` numerics, which are correct: `fr = 2T dr/d2`,
`fi = -2T di/d2`, geometric sum, Dyson). Change: id and classes (`b7-resum-…`); `Ω` → `Π^{ph}`; `ω_σ` → `ω_{TO}`;
"Dyson, E Eq. 404" → "Dyson sum"; y label `A = −Im F/π` → `−Im F/π (meV^{-1})`; readout labels in v2 symbols; the
preset data and $\omega^{\ion}_{\mathrm p}$ now come from the literature Born charges (v1 derived $\omega_{\mathrm p}$ from the measured LO,
which gave GaAs 33.0 meV instead of 29.8 meV for the static self-energy); fixed x-range as above; the shared view switch
and the `lam` control are removed.

---

## 3. `fig-limits.html` · id `viz-b7-limits` · PORT (view "Order of limits" of v1 `viz-b-phonon-se-limits`)

**Title.** "Order of limits: phonon polaritons, LO and TO" `<span class="tag">interactive</span>`.

**What the reader must learn.** On a logarithmic $q$ axis: for $q\ll q_{\mathrm c}$ the upper polariton meets the LO line and
carries all the transverse phonon weight, while the lower branch is a photon; for $q\gg q_{\mathrm c}$ the TO phonon is back at
$\omega_{\mathrm{TO}}$ and the LO stays at $\omega_{\mathrm{LO}}$. Switching the electrons off moves the $q\to0$ mode to the
electron-free value (B7.13), the result of the EPP manuscript. Raising the speed of light pushes $q_{\mathrm c}$ off the left
edge: what remains is the electrostatic pair of (B7.21), flat TO and flat LO.

**Include position.** At the end of h4 `b7-finiteq`, after the derivation box "several polar modes, and why the basis
matters" and directly before h3 `b7-limits`. The sentences just before the box: "… at $10\,q_{\mathrm c}$ the lower polariton
carries 0.998. The phonon weight is not the share of energy in the lattice: the upper polariton at $q\to0$ also stores
energy in its macroscopic field." The order-of-limits formulas (B7.20)–(B7.22) follow the figure.

**Physics.** Cubic crystal with two atoms per cell, electrons as a constant background, speed of light $c'=s_c\,c$ with a
factor $s_c$ (call it "speed of light factor" in the UI; not $\lambda$):
- $\varepsilon=\varepsilon_\infty$ (electrons on) or $1$ (off); $\omega^2_{\mathrm{LO}}=\omega^2_{\mathrm{TO}}+(\omega^{\ion}_{\mathrm p})^2/\varepsilon$ (flat in $q$);
- polaritons, roots of (B7.17) `b7-eq-branches` with $c\to c'$: $b=\varepsilon\omega^2_{\mathrm{LO}}+c'^2q^2$,
  $\omega^2_{\mathrm{UP}}=[b+(b^2-4\varepsilon c'^2q^2\omega^2_{\mathrm{TO}})^{1/2}]/(2\varepsilon)$, and for numerical stability
  $\omega^2_{\mathrm{LP}}=c'^2q^2\omega^2_{\mathrm{TO}}/(\varepsilon\,\omega^2_{\mathrm{UP}})$ (product of the roots);
- phonon weights (B7.18) `b7-eq-weight`: $Z^{\ph}=[1+(\omega^{\ion}_{\mathrm p})^2c'^2q^2/(\varepsilon\omega^2-c'^2q^2)^2]^{-1}$ at each root,
  with the unscreened $(\omega^{\ion}_{\mathrm p})^2$; $Z^{\ph}_{\mathrm{UP}}+Z^{\ph}_{\mathrm{LP}}=1$ (v1's
  `X = cq² (ω_UP² − ω_TO²)² / (wp2 ω_UP²)`, `Z_UP = 1/(1+X)`, is the same expression evaluated on the root);
- light line $c'q/\sqrt{\varepsilon}$; crossover $q_{\mathrm c}=\sqrt{\varepsilon}\,\omega_{\mathrm{TO}}/c'$ ((B7.22) `b7-eq-qc`);
- electron-free $q\to0$ value $[\omega^2_{\mathrm{TO}}+(\omega^{\ion}_{\mathrm p})^2]^{1/2}$ ((B7.13) `b7-eq-410`);
- static electron-free self-energy $(\hbar\omega^{\ion}_{\mathrm p})^2/(2\hbar\omega_{\mathrm{TO}})$ ((B7.12) `b7-eq-pi0`).
- $q$ grid: 361 points, log-spaced from $10^{-7}$ to $10^{-2}$ Å⁻¹. With $\hbar c'q$ in meV: `cq = s_c * 1.97327e6 * q`.

**What to draw.** Two panels side by side on wide screens (`viz-grid side`), stacked on phones.
- (a) $\hbar\omega$ (meV) against $q$ (Å⁻¹), log $x$ from $10^{-7}$ to $10^{-2}$; $y$ from 0 to a fixed value per preset,
  `ceil(1.16 × electron-free value / 5) × 5` (GaAs 65, NaCl 55, hBN 330 meV), independent of the toggles.
  Curves: LO flat (`phonon`, legend "LO, longitudinal"); upper polariton $\omega_{\mathrm{UP}}$ (`c5`, solid); lower polariton
  $\omega_{\mathrm{LP}}$ (`c4`, dashed 7/4); light line $c'q/\sqrt{\varepsilon}$ (`photon`, dotted 2/3, no hover); horizontal
  $\omega_{\mathrm{TO}}$ (`muted`, label "ω_{TO}"); horizontal electron-free value (`ink`, dotted, label "no electrons, q → 0"),
  shown in both toggle states; a short tick at the right edge at the measured $\hbar\omega_{\mathrm{LO}}$ (`muted`, label
  "measured LO"); vertical $q_{\mathrm c}$ (`muted`, label "q_c") when inside the frame; a vertical cursor line at the
  chosen $q$ (accent, thin). Near the left edge, label the upper branch "ω_{UP} → ω_{LO}" when it is phonon-like there
  ($Z^{\ph}_{\mathrm{UP}}>0.5$), else "LO".
- (b) phonon weight $Z^{\ph}$ (0–1.06, ticks 0, 0.25, 0.5, 0.75, 1) against the same log $q$ axis (ticks $10^{-7}$,
  $10^{-5}$, $10^{-3}$): $Z^{\ph}_{\mathrm{UP}}$ (`c5`), $Z^{\ph}_{\mathrm{LP}}$ (`c4` dashed), vertical $q_{\mathrm c}$ and the cursor;
  text labels "phonon-like" and "photon-like" on the corresponding branches at the left edge.

**Controls (4).**
- Crystal preset: segmented GaAs, NaCl, hBN ∥ (default GaAs).
- Electrons: toggle "electrons on ($\varepsilon_\infty$ in $D^{\rme}$)" (default on).
- Speed of light factor $s_c$: log range 1 to $10^6$, readout "×1.0", "×10", "×1.0·10³" (default 1). This is the
  thought experiment $c\to\infty$.
- $q$ cursor: log range $10^{-7}$ to $10^{-2}$ Å⁻¹ (default: the preset's $q_{\mathrm c}$ at $s_c=1$, GaAs $5.57\times10^{-5}$ Å⁻¹);
  the cursor can also be set by clicking or dragging in either panel.

**Readouts** (HTML, MathJax): $\hbar\omega^{\ion}_{\mathrm p}$; static $\Pi^{\ph,0}$ (B7.12); electron-free $q\to0$ value (B7.13);
$\hbar\omega_{\mathrm{LO}}$ (B7.20, with the current toggle); $q_{\mathrm c}$ (B7.22) and its fraction of $2\pi/a$ (scientific
notation as HTML, e.g. 5.6×10<sup>−5</sup>); at the cursor: $\hbar\omega_{\mathrm{UP}}$, $\hbar\omega_{\mathrm{LP}}$,
$Z^{\ph}_{\mathrm{UP}}$. GaAs at the default cursor: 40.73 meV, 27.23 meV, 0.599.

**Initial state.** GaAs, electrons on, $s_c=1$, cursor at $q_{\mathrm c}$: the anticrossing is centred, the upper branch reaches
the LO line on the left, the TO line on the right is met by the lower branch.

**Caption key points.**
- *What to notice.* Left of $q_{\mathrm c}$ the transverse phonon has become the upper polariton at the LO frequency; right of it
  TO and LO are the electrostatic pair. The order of the limits $q\to0$ and $c\to\infty$ decides which you get; the
  electrons decide how large the shift is.
- *Marks.* Colours and line styles as listed; dotted line: the electron-free value of the EPP manuscript; tick: measured
  LO; vertical lines: $q_{\mathrm c}$ and the cursor.
- *Try.* Switch the electrons off: the $q\to0$ modes climb to the dotted line (GaAs 55.6 meV). Raise the speed of light:
  $q_{\mathrm c}$ slides off the left edge and only flat TO and LO remain, the limit $c\to\infty$ of Part A. Drag the cursor across
  $q_{\mathrm c}$ and watch the phonon weight pass from one branch to the other.
- *Model.* <a class="eqref" href="#b7-eq-branches"></a>, <a class="eqref" href="#b7-eq-weight"></a>,
  <a class="eqref" href="#b7-eq-410"></a>, <a class="eqref" href="#b7-eq-limq"></a>, <a class="eqref" href="#b7-eq-limc"></a>,
  <a class="eqref" href="#b7-eq-qc"></a>, with the electrons as a constant $\varepsilon_\infty$ (insulator below its gap), macroscopic fields
  only (the bare optical frequency equals $\omega_{\mathrm{TO}}$), literature Born charges from the table (GaAs gives
  $\hbar\omega_{\mathrm{LO}}=35.9$ meV against the measured 36.2 meV), hBN treated as isotropic, no damping.

**v1 to v2.** Keep v1's `trans` and `drawLimits` numerics (correct). Change:
- id and classes (`viz-b7-limits`, `b7-limits-…`); drop the second view (now `fig-resum`) and its controls `N`, `G`, and
  the `view` switch; add the $q$ cursor with its readouts.
- Notation: `ω_+`, `ω_–` → `ω_{UP}`, `ω_{LP}`; `Z_+`, `Z_–` → `Z^{ph}_{UP}`, `Z^{ph}_{LP}`; `ω_σ` and `ω_– → ω_σ = TO` → `ω_{TO}`;
  `E Eq. 410′` (canvas label and readout) → "no electrons, q → 0" with the caption citing (B7.13); `E Eq. 410e`, `410h`,
  `409` → (B7.20), (B7.22), (B7.12); `ħω_p (unscreened)` → `ħω_p^{ion}`; `λ_c` → speed of light factor; no `σ` anywhere.
- Physics data: v1 fixed $\omega_{\mathrm p}$ from the measured splitting, $\omega^2_{\mathrm p}=\varepsilon_\infty(\omega^2_{\mathrm{LO}}-\omega^2_{\mathrm{TO}})$, which
  implies $Z^*=2.18$ for GaAs and gives 57.5 meV for the electron-free value against 55.6 meV in the text. v2 uses the
  literature Born charges of the table, so the figure reproduces the text (55.6, 35.9, 33.3 meV), and shows the measured
  LO as a reference tick. Remove the "implied Z*" readout.
- Caption rewritten in the v2 structure; the v1 caption mixed model, marks and notices in one block.
