# A.8 figure specs (section `a-temperature`, prefix `a8`)

Four figures: one new (`qp`) and three ports of v1 (`varshni`, `lo`, `catastrophe`). Figure numbers follow the
include order: A8.1 `qp`, A8.2 `varshni`, A8.3 `lo`, A8.4 `catastrophe`.

Common rules for all four (BUILDER_GUIDE §5): file `v2/sections/A08/fig-<name>.html`, root
`<figure class="viz" id="viz-a8-<name>">`, CSS classes prefixed `a8-<name>-`, script guarded by
`if (!fig || !window.EPP) return;`. Colours read at draw time with `EPP.color(...)`. Canvas text with
`EPP.drawRich` in v2 notation (upright labels: `ω_{LO}`, `E^{qp}_g`, `q_{pol}`; italic variables). No
equation numbers written into canvas text or readouts (they go stale): the caption cites equations with
`<a class="eqref" href="#a8-eq-…"></a>`. Readouts in HTML with MathJax `$…$`. Controls below the plots; 390 px
without horizontal overflow; keyboard access and `aria-label`s; no NaN (guard `T = 0`, `x = 0`, log singularities).
`σ` must not appear in any string (it is reserved for the conductivity). Constants: ħ = 658.212 meV fs =
6.582120×10⁻¹⁶ eV s, k_B = 0.0861733 meV/K, ħ²/2m_e = 3.80998 eV Å², ħc = 1973.27 eV Å, Ry = 13605.693 meV.

Material parameters shared by `qp`, `lo` and `catastrophe` (literature values, illustrative; the same as the
section's note on (A8.3)):

| material | m*/m_e | ħω_LO (meV) | ε∞ | ε_s | Ω_V (Å³) | α_F (computed) |
|---|---|---|---|---|---|---|
| GaAs | 0.067 | 36.2 | 10.9 | 12.9 | 45.2 | 0.0714 |
| GaN | 0.20 | 92 | 5.35 | 9.5 | 45.7 | 0.444 |
| MgO | 0.35 | 89 | 2.96 | 9.8 | 18.7 | 1.725 |

Compute α_F = (1/ε∞ − 1/ε_s) · (m*/m_e)^{1/2} · (Ry/ħω_LO)^{1/2}; q_pol = (m*/m_e · ħω_LO[eV] / 3.80998)^{1/2} Å⁻¹
(GaAs 0.02523, GaN 0.06949, MgO 0.09042 Å⁻¹); q_D = (6π²/Ω_V)^{1/3} (GaAs 1.0942, GaN 1.0902, MgO 1.4685 Å⁻¹);
q_lc = √ε∞ · ħω_LO[eV] / 1973.27 Å⁻¹ (GaAs 6.057×10⁻⁵, GaN 1.078×10⁻⁴, MgO 7.760×10⁻⁵ Å⁻¹).
All numbers quoted below were checked numerically by the section author.

---

## 1. `fig-qp.html` · id `viz-a8-qp` · new

**Title.** "Quasiparticle energy of a polar band edge: on shell and off shell" `<span class="tag">interactive</span>`

**What the reader must learn.**
- The pole condition (A8.2) as a geometric construction: the quasiparticle energy is where the straight line
  ħω − E_k meets the curve Re Σ_k(ω).
- On shell = the curve read at ħω = E_k; off shell (self-consistent) = the intersection; linearized = the
  intersection of the line with the tangent of Re Σ at E_k, whose slope gives Z^qp.
- At weak coupling (GaAs) the three agree; at α_F ≈ 1.7 (MgO) they spread, a sign that second-order theory fails.
- Above the emission threshold ħω = ħω_LO the self-energy has an imaginary part: a width and a lifetime.
- The renormalized (polaron) band flattens as it approaches the threshold.

**Include position.** In h4 `a8-polaron` ("A worked example: the band edge of a polar semiconductor"), right after
the paragraph that ends: "…It is the inverse of the polaron length r_p of Section A.4, the size of the cloud of
lattice polarization that dresses the electron: 1/q_pol = r_p = 40 Å in GaAs. Figure A8.1 plots the self-energy and
solves (A8.2) graphically: the solution is where the straight line ħω − E_k crosses the curve Re Σ." The closed form
(A8.4) follows the figure.

**Physics.** One parabolic band E_k = ħ²k²/2m* (measured from its minimum), empty (f = 0), T = 0 (n_B = 0), one
dispersionless LO phonon, Fröhlich coupling |g(q)|² of (A8.3) (`a8-eq-gfr`). Self-energy (A8.4) (`a8-eq-sigfr`) and
its continuation above the threshold. Work in dimensionless variables:
- x = ħω/ħω_LO (energy variable, from the band minimum), y = E_k/ħω_LO, s = √y = k/q_pol;
- S(x, y) = Σ_k(ω)/(α_F ħω_LO), complex, from these piecewise closed forms (verified against direct numerical
  integration):
  - s = 0: x < 1: Re S = −1/√(1−x), Im S = 0; x > 1: Re S = 0, Im S = −1/√(x−1).
  - s > 0, x < 1: Re S = −arctan(s/√(1−x))/s, Im S = 0.
  - s > 0, x > 1, r = √(x−1): if s < r: Re S = 0, Im S = −artanh(s/r)/s; if s > r: Re S = −π/(2s),
    Im S = −ln[(s+r)/(s−r)]/(2s).
  - At x = 1 + y (s = r) Im S has a logarithmic singularity and Re S jumps from −π/(2s) to 0: lift the pen
    (NaN) there and clip to the frame.
- On shell (x = y): Re S = −arcsin(s)/s for y < 1, −π/(2s) for y > 1 (−π/2 at y = 1); Im S = 0 for y ≤ 1,
  −arcosh(s)/s for y > 1.
- Off shell, self-consistent: solve x − y = α_F Re S(x, y) by bisection on [y − 3α_F − 1, 1 + y − 10⁻⁹]
  (for s = 0 the upper end is 1 − 10⁻⁹). f(x) = x − y − α_F Re S is continuous and increasing on that interval,
  negative at the lower end and positive at the upper end, so the root x* is unique.
- Linearized: ∂Re S/∂x at x = y is −1/(2√(1−y)) for y < 1 and 0 for y > 1, so
  Z^qp = 1/(1 + α_F/(2√(1−y))) below the threshold (→ 0 as y → 1⁻) and Z^qp = 1 above; the linearized energy is
  x_lin = y + Z^qp α_F Re S(y, y).
- Band-edge check values (k = 0, in meV): GaAs −2.58 (on shell), −2.50 (self-consistent), −2.49 (linearized),
  Z = 0.966; GaN −40.9, −34.8, −33.4, Z = 0.818; MgO −153.5, −104.2, −82.4, Z = 0.537.
- Above the threshold, GaAs at E_k = 2ħω_LO: half width on shell Γ^τ_k = −Im Σ_k(E_k/ħ) = α_F ħω_LO · arcosh(√2)/√2 = 1.61 meV,
  τ_k = ħ/(2Γ^τ_k) = 204 fs (equals the T = 0 value of the `lo` figure).

**What to draw.** Two panels in `viz-grid two`, stacked on phones.
- (a) "Graphical solution at the chosen k". x axis: ħω/ħω_LO from −2.5 to 3 (label `ħω / ħω_{LO}`); y axis: energy in
  units of ħω_LO from −2.5 to 1.5 (label `energy / ħω_{LO}`), fixed limits.
  - α_F Re S(x, y) as a solid curve (`electron`), label `Re Σ_k`, clipped, pen lifted at x = 1 + y and near the
    k = 0 divergence at x → 1⁻.
  - α_F Im S(x, y) as a dashed curve (`c5`), label `Im Σ_k` (≤ 0), clipped.
  - The straight line x − y (`ink`, thin), label `ħω − E_k`.
  - A vertical dashed line at x = y (`muted`), label `on shell`, with a hollow marker at (y, α_F Re S(y, y)).
  - A filled marker at the self-consistent root (x*, x* − y) (`electron`), label `off shell`.
  - When the toggle is on: the tangent α_F[Re S(y,y) + ∂Re S/∂x·(x − y)] (dotted, `c4`) and a square marker at its
    intersection with the line, label `linearized`.
  - A vertical line at x = 1 (`phonon`), label `ħω_{LO}`: emission threshold.
- (b) "Dispersion of the dressed band". x axis: k/q_pol from 0 to 1.6; y axis: energy/ħω_LO from −2 to 2.6, fixed.
  - Bare band y = s² (`muted`, dashed), label `E_k`.
  - On-shell dispersion y + α_F Re S(y, y) (`electron`, dashed), label `on shell`.
  - Self-consistent dispersion x*(y) (`electron`, solid), label `off shell`, with a shaded band of ± its linewidth
    Γ^τ_k = α_F|Im S(x*, y)| (A8.2) around it where that is nonzero (`electron-t`); legend `off shell ± Γ^τ_k`.
  - Linearized dispersion y + Z^qp α_F Re S(y, y) (`c4`, dotted) when the toggle is on.
  - Horizontal line at 1 (`phonon`), label `ħω_{LO}`.
  - Cursor: vertical line at the chosen k and markers on the three curves. Click or drag in this panel sets k.

**Controls.**
1. material, seg: GaAs | GaN | MgO (sets α_F and ħω_LO from the table).
2. cursor, range: `E_k/ħω_{LO}` from 0 to 2.5, step 0.01 (k = q_pol √(E_k/ħω_LO)).
3. toggle: "show the linearized solution".

**Readouts.** $\alpha_{\mathrm F}$; $E_{\bk}$ (meV); on shell $\mathrm{Re}\,\Sigma_{\bk}(E_{\bk}/\hbar)$ (meV); off shell
$E^{\mathrm{qp}}-E_{\bk}$ (meV); linearized $Z^{\mathrm{qp}}\,\mathrm{Re}\,\Sigma_{\bk}$ (meV); $Z^{\mathrm{qp}}$;
on-shell half width $\Gamma^{\tau}_{\bk}=-\mathrm{Im}\,\Sigma_{\bk}(E_{\bk}/\hbar)$ (meV); lifetime $\tau_{\bk}=\hbar/2\Gamma^{\tau}_{\bk}$ (fs, "∞ below
the threshold"). Energies in meV are the dimensionless values times ħω_LO of the material.

**Initial state.** GaN, E_k = 0 (band edge), toggle on: the three markers are visibly distinct (−0.44, −0.38,
−0.36 in units of ħω_LO).

**Caption.**
- *What to notice.* At the band edge the on-shell value is the polaron shift −α_Fħω_LO. The self-consistent and
  linearized solutions differ from it by terms of order α_F²: 3% for GaAs, 15 to 18% for GaN, 30 to 50% for MgO,
  where second-order theory is no longer enough. Above E_k = ħω_LO the self-energy has an imaginary part: the
  electron can emit a real LO phonon.
- *Marks.* Solid and dashed curves in (a): Re Σ_k and Im Σ_k; thin line: ħω − E_k; hollow marker: on shell;
  filled marker: self-consistent solution; square and dotted tangent: linearized solution; vertical line: the
  emission threshold. In (b): bare band, on-shell and self-consistent dispersions, the shaded linewidth ±Γ^τ_k, the cursor.
- *Try.* Switch to MgO and watch the three markers part. Move k towards E_k = ħω_LO: the curve steepens near the
  threshold, Z^qp drops towards 0 and the dressed band flattens below the threshold. Move above it in GaAs:
  at E_k = 2ħω_LO the on-shell Γ^τ_k is 1.6 meV and τ_k = 204 fs.
- *Model.* (A8.3) and (A8.4) with their continuation above the threshold, solved with (A8.2); empty band, T = 0,
  one dispersionless LO phonon; parameters of the table (literature values, illustrative). Cite with
  `a8-eq-gfr`, `a8-eq-sigfr`, `a8-eq-bw`.

---

## 2. `fig-varshni.html` · id `viz-a8-varshni` · port of v1 `viz-a-temperature-varshni` (v1 file `sections/A08-temperature.html`)

**Title.** "The Varshni curve: zero-point plateau, linear regime and the ZPR from the asymptote" `interactive`

**What the reader must learn.**
- The shape of E^qp_g(T) from (A8.12): a plateau set by zero-point motion, then a linear fall.
- The straight high-temperature asymptote extrapolates to the clamped-nuclei gap E_g, not to the measured
  E^qp_g(0); the difference is the ZPR. This is how ZPRs are read off experiments.
- A softer phonon starts the linear regime earlier and makes it steeper (slope −2k_B|ΔE^ZPR_g|/ħω_0).
- In the one-phonon model the half width of a state away from the band edge follows the same Bose factor.

**Include position.** In h3 `a8-varshni`, after the bullet list of model assumptions and the sentence
"Figure A8.2 shows the result, which the following steps derive." The derivation steps ending in (A8.12) follow.

**Physics.**
- Gap: E^qp_g(T) = E_g − |ΔE^ZPR_g| coth(ħω_0/2k_BT), with coth → 1 at T = 0 (`a8-eq-varshni`).
- Asymptote: E_g − |ΔE^ZPR_g| · 2k_BT/ħω_0. Plateau edge: T = ħω_0/(3k_B), where coth(3/2) = 1.105.
- Half width (one-phonon model, lifetime subsection): Γ^τ(T) = Γ^τ(T = 0) coth(ħω_0/2k_BT); lifetime
  τ = ħ/(2Γ^τ) = 658.212 meV fs / (2Γ^τ[meV]) (`a8-eq-taudef`).
- Presets (all ZPR as positive magnitudes):

| preset | E_g (eV) | ZPR (meV) | ħω_0 (meV) | Γ^τ(T = 0) (meV) | T_max (K) | y range (eV) |
|---|---|---|---|---|---|---|
| Toy, ħω_0 = 50 meV (matches the printed curves) | 1.00 | 100 | 50 | 50 | 500 | 0.75–1.05 |
| Toy, ħω_0 = 100 meV (as captioned) | 1.00 | 100 | 100 | 50 | 500 | 0.75–1.05 |
| Diamond direct gap, GW coupling | 7.73 | 628 | 125 | unchanged | 800 | 6.70–7.80 |
| Diamond direct gap, LDA coupling | 7.73 | 400 | 125 | unchanged | 800 | 6.70–7.80 |
| Si indirect gap | 1.23 | 60 | 30 | unchanged | 600 | 1.00–1.25 |
| custom | unchanged | slider | slider | slider | unchanged | expand if needed |

  Check values: toy 50 meV: E^qp_g(500 K) = 0.809 eV, asymptote 0.828 eV, slope −0.345 meV/K, plateau edge 193 K;
  toy 100 meV: 0.878 eV and 0.914 eV at 500 K; diamond GW: 7.102 eV at 0 K, 7.092 at 300 K, 7.029 at 500 K, 6.951 at
  650 K, slope −0.866 meV/K, plateau edge 484 K; diamond LDA: 7.33 eV at 0 K, 7.234 at 650 K; Si: 1.170 eV at 0 K,
  1.115 eV at 300 K, slope −0.345 meV/K, plateau edge 116 K; Γ^τ = 50 meV → τ = 6.58 fs.

**What to draw.** Two panels in `viz-grid two`.
- (a) gap: x axis T (K) from 0 to T_max of the preset; y axis E^qp_g (eV), limits fixed per preset (table), expanded
  only if the curve leaves the frame. Layers: plateau band from 0 to ħω_0/3k_B (`muted`, alpha ≈ 0.09, text
  `plateau`); horizontal line at E_g (`muted`, text `E_g`); asymptote (`ink`, dashed, label `high-T asymptote`);
  the curve E^qp_g(T) (`electron`, label `E^{qp}_g(T)`); a double arrow at T ≈ 0.08 T_max from E_g down to
  E^qp_g(0), text `|ΔE^{ZPR}_g|`; temperature cursor (vertical dashed line and point).
- (b) seg-selectable: half width Γ^τ(T) (meV) with a horizontal line at Γ^τ(T = 0) (text `Γ^τ(T = 0)` with its value), or lifetime
  τ(T) (fs) with a horizontal line at ħ/2Γ^τ(T = 0) (text with its value); cursor as in (a). y limits fixed for the
  current preset and width.
- Clicking either panel moves the temperature cursor T*.

**Controls (5).**
1. preset, select (table; label "preset"), wide.
2. `|ΔE^{ZPR}_g|`, range 0–700 meV, step 1.
3. `ħω_0`, range 10–200 meV, step 1.
4. `Γ^{τ}(T = 0)` (half width at T = 0), range 1–100 meV, step 1.
5. right panel, seg: "half width" | "lifetime".
Moving a slider away from a preset's values switches the select to "custom" (as in v1).

**Readouts.** $E_{\mathrm g}$; $E^{\mathrm{qp}}_{\mathrm g}(0)$; $E^{\mathrm{qp}}_{\mathrm g}(T^*)$ with T* in K;
high-T slope (meV/K); plateau edge $\hbar\omega_0/3\kB$ (K); $\Gamma^{\tau}(T^*)$ (meV); $\tau(T^*)$ (fs).

**Initial state.** Preset "Toy, ħω_0 = 50 meV", T* = 300 K, right panel "half width".

**Caption.**
- *What to notice.* The plateau is zero-point motion. The dashed asymptote meets E_g at T = 0, and the arrow between
  E_g and the low-temperature gap is the ZPR. The half width follows the same Bose factor as the shift.
- *Marks.* Blue curve: E^qp_g(T); dashed line: high-temperature asymptote; grey line: E_g; arrow: |ΔE^ZPR_g|; shaded
  band: k_BT < ħω_0/3; right: half width or lifetime; dashed vertical lines: the temperature cursor.
- *Try.* Lower ħω_0 at fixed ZPR: the plateau shortens and the slope −2k_B|ΔE^ZPR_g|/ħω_0 steepens. Pick "Toy,
  ħω_0 = 100 meV" and compare with the review's Fig. 3, whose curves follow 50 meV. Pick diamond with GW coupling:
  7.10 eV at 0 K and 6.95 eV at 650 K, as in the review's Fig. 13. Switch to lifetime: 50 meV gives 6.6 fs.
- *Model.* (A8.12) (`a8-eq-varshni`), the one-phonon half width of the lifetime subsection and τ from (A8.13)
  (`a8-eq-taudef`). Toy parameters from the review's Fig. 3. Diamond: E_g = 7.73 eV and ħω_0 = 125 meV from a
  one-oscillator fit to the review's Fig. 13, ZPR 628 meV with GW coupling or 0.40 eV read off its LDA curve. Si:
  literature values, illustrative. Γ^τ(T = 0) outside the toy presets is illustrative.

**What must change from v1.**
- id `viz-a-temperature-varshni` → `viz-a8-varshni`; classes `pv-*` → `a8-varshni-*`; script wrapped as in the
  common rules (v1 used one shared script for three figures; split out this figure's code).
- Notation in labels, readouts, aria-labels and option strings: `E_g^{cn}` → `E_g`; `E_g(T)` → `E^{qp}_g(T)`;
  `ZPR`, `|ΔE_g^{ZP}|` → `|ΔE^{ZPR}_g|`; `Γ`, `Γ^{ZP}`, `Γ(T)` → `Γ^{τ}`, `Γ^{τ}(T = 0)`, `Γ^{τ}(T)`
  (the superscript τ keeps the linewidth apart from the zone centre Γ); `τ = ħ/2Γ` → `τ = ħ/2Γ^{τ}`; `n(ħω_0/k_BT)` → `n_B(ω_0)`; `ħω_0/3k_B` keeps its form with `k_B` upright.
- Remove the curve label "E_g(T), Eq. (G.168)" and every source equation number in strings.
- Remove the `T_max` slider (fixed per preset) and the "inverse Varshni" sign toggle (the section discusses the
  inverse effect in words): 7 controls → 5.
- v1 recomputed the y limits on every slider move (axis jumps); use fixed limits per preset.
- Rewrite the caption in the What to notice / Marks / Try / Model structure; v1's caption cites "(G.168)" and Fig. 3
  in source notation.

---

## 3. `fig-lo.html` · id `viz-a8-lo` · port of v1 `viz-a-temperature-lo` (v1 file `sections/A08-temperature.html`), with a new geometric panel

**Title.** "LO-phonon scattering in a polar semiconductor: the emission threshold" `interactive`

**What the reader must learn.**
- Geometrically: a carrier can emit an LO phonon only if its band has a state one phonon energy lower; below
  E_k = ħω_LO that state would lie below the band bottom. Absorption is always possible but needs a phonon.
- In the rate (A8.15): absorption ∝ n_B(ω_LO), emission ∝ n_B + 1 with a square-root onset at the threshold.
- Lifetime against T: below the threshold it diverges as T → 0 (no phonons to absorb); above it, it levels off at
  the spontaneous-emission value (GaAs about 200 fs, GaN about 13 fs at E_k = 2ħω_LO).

**Include position.** In h4 `a8-lo`, after the paragraph that ends "…Below this threshold it can only absorb phonons,
at a rate proportional to n_B(ω_LO), which vanishes as T → 0. Figure A8.3 shows the threshold in the band diagram and
in the rates." The closed form (A8.15) follows the figure. The section's derivation box refers to the band diagram
as **panel (a)**: keep that panel first.

**Physics.**
- (A8.15) (`a8-eq-lorate`), with x = E_k/ħω_LO, n = n_B(ω_LO) = 1/(exp(ħω_LO/k_BT) − 1) (0 at T = 0):
  rate = R0 · [n · arsinh(√x)/√x + (n + 1) · arsinh(√(x − 1))/√x · Θ(x − 1)], R0 = 2α_Fω_LO in ps⁻¹ =
  2α_F · (ħω_LO[eV]/6.582120×10⁻¹⁶) × 10⁻¹²; τ = 1000/rate fs.
- GaAs: α_F = 0.0714, R0 = 7.85 ps⁻¹; GaN: α_F = 0.444, R0 = 124.1 ps⁻¹. Check values: GaAs E_k = 2ħω_LO:
  τ = 204 fs (T = 0), 117 fs (300 K); E_k = ħω_LO/2: 418 fs (300 K), ∞ at T = 0; E_k = 1.5ħω_LO: 237 fs (0 K),
  129 fs (300 K); GaN E_k = 2ħω_LO: 12.9 fs (0 K), 12.1 fs (300 K).
- Band diagram (geometry of the delta functions in (A8.14)): band E_k/ħω_LO = (k/q_pol)². A carrier at k* = √x*
  has emission final states at |k'| = √(x* − 1) (exist only for x* > 1) and absorption final states at
  |k'| = √(x* + 1). In a 1D cut through k*, the momentum transfers run from k* − k' to k* + k', the integration limits
  of the derivation box (emission: k ± √(k² − q_pol²); absorption: √(k² + q_pol²) ∓ k, in units of q_pol).
- Mean free path: v = ħk/m* = (2 · 3.80998 eV Å² · k[Å⁻¹]/(m*/m_e))/(6.582120×10⁻¹⁶ eV s), in Å/fs after × 10⁻¹⁵;
  ℓ = vτ in nm.

**What to draw.** Three panels: (a) on top, full width, compact (about 200–230 px tall, `EPP.canvas`); (b) and (c)
below in `viz-grid two` (stacked on phones); controls after all panels.
- (a) "Band diagram": k/q_pol from −2.4 to 2.4 (horizontal), energy/ħω_LO from −0.6 to 6.2 (vertical; axis labels
  `k / q_{pol}` and `E / ħω_{LO}`). Parabola (`electron`). The carrier as a filled dot at (k*, x*) on the right
  branch. Emission level x* − 1: horizontal dashed line (`phonon`) with hollow dots at ±√(x* − 1) and a downward
  arrow of length 1 from the carrier, text `emit ħω_{LO}` and the weight `n_B + 1 = …`. Absorption level x* + 1:
  dashed line (`c5`) with hollow dots at ±√(x* + 1), upward arrow, text `absorb ħω_{LO}` and `n_B = …`. Thin
  horizontal arrows on each final level from k* to +k' and to −k' (smallest and largest momentum transfer). When
  x* < 1, draw the emission level greyed below the band bottom with the text `no final state: E_k < ħω_{LO}`.
- (b) rate against carrier energy: x axis `E_k / ħω_{LO}` from 0 to 5; y axis `1/τ (ps^{-1})` from 0 to a fixed top
  per material (the 600 K maximum × 1.08, rounded). Lines: total at T (`electron`, solid, label `total`), emission
  (`phonon`, dashed, `emission ∝ n_B + 1`), absorption (`c5`, dotted, `absorption ∝ n_B`), total at T = 0 (`muted`,
  `total at T = 0`); vertical line at 1 (`phonon`, text `emission threshold`); cursor line and point. Click sets
  the cursor.
- (c) lifetime against temperature: x axis T (K) from 0 to 600; y axis τ (fs), log, 1 to 10⁶. Curves at
  E_k = ħω_LO/2 (`muted`, dotted, `below threshold`), E_k = 2ħω_LO (`c4`, dashed, `above threshold`) and at the
  cursor energy (`electron`); vertical line at T and a point. Return NaN where the rate is 0 (pen lifted).

**Controls (3).**
1. material, seg: GaAs | GaN.
2. `T`, range 0–600 K, step 5.
3. cursor `E_k / ħω_{LO}`, range 0.05–5, step 0.01.

**Readouts.** $\alpha_{\mathrm F}$; $2\alpha_{\mathrm F}\omega_{\mathrm{LO}}$ (ps⁻¹); $n_{\mathrm B}(\omega_{\mathrm{LO}})$
at T; $E_{\bk}$ (meV); $1/\tau$ (ps⁻¹); $\tau$ (fs); $\hbar/\tau=2\Gamma^{\tau}$ (meV); mean free path $v\tau$ (nm).

**Initial state.** GaAs, T = 300 K, cursor 1.5 (above the threshold; τ = 129 fs).

**Caption.**
- *What to notice.* The rate jumps on at E_k = ħω_LO, rising as (E_k − ħω_LO)^{1/2}, where spontaneous emission
  becomes possible; in the band diagram the emission level reaches the band bottom exactly there. As T → 0 the
  lifetime below the threshold diverges, while above it levels off at the spontaneous-emission value (about 200 fs
  in GaAs, 13 fs in GaN at E_k = 2ħω_LO).
- *Marks.* (a): the band, the carrier, the final levels one LO energy below and above with their final states,
  and the extreme momentum transfers; (b): total, emission and absorption rates at T and the total at T = 0;
  (c): lifetimes at three energies.
- *Try.* Set T = 0 and move the cursor below 1: the rate vanishes and the emission level drops below the band.
  Raise T at E_k = ħω_LO/2: absorption switches on as n_B grows. Switch to GaN: α_F is six times larger and the
  rates about 16 times faster.
- *Model.* (A8.15), derived from (A8.14) with the Fröhlich coupling (A8.3): one parabolic band, one dispersionless
  LO phonon, empty band (f = 0), no nonparabolicity, no other valleys, no Z^qp factor. Parameters of the table
  (literature values, illustrative). Cite `a8-eq-lorate`, `a8-eq-tau`, `a8-eq-gfr`.

**What must change from v1.**
- id `viz-a-temperature-lo` → `viz-a8-lo`; classes `pl-*` → `a8-lo-*`; own script (v1 shared one script).
- Notation: the carrier energy `ε`, `ε^*` → `E_k` (ε is the dielectric function in v2); `α` → `α_F`;
  `n(ħω_{LO}, T)` → `n_B(ω_{LO})`; remove "Eq. (G.170)" from the curve label; axis label `carrier energy E_k / ħω_{LO}`.
- Add panel (a), the geometric band diagram (new); v1 had only the two plots.
- Rewrite the caption in the four-part structure. v1's physics (functions `rAbs`, `rEm`, R0, the mean free path) is
  correct and can be kept.

---

## 4. `fig-catastrophe.html` · id `viz-a8-catastrophe` · port of v1 `viz-a-temperature-frohlich` (v1 file `sections/A08-temperature.html`)

**Title.** "The Fröhlich catastrophe: convergence of the band-edge shift with the q grid" `interactive`

**What the reader must learn.**
- With the phonon energy kept, the shift of the conduction-band minimum converges to −α_Fħω_LO once the grid
  resolves q ≲ q_pol (N_1 ≳ q_D/q_pol); the adiabatic shift grows linearly with N_1 without limit.
- Where the contributions come from: per logarithmic interval of q, the adiabatic integrand grows as 1/q, the
  nonadiabatic one peaks at q_pol and vanishes below it.
- The light cone q < q_lc lies 400 to 1200 times below q_pol: finite c is not what cures the divergence.

**Include position.** In h3 `a8-catastrophe`, after the paragraph that ends "…A uniform grid of N_1 points along
each direction resolves wavevectors down to about q_min = q_D/N_1. Figure A8.4 follows both prescriptions as the grid
is refined." The closed forms (A8.17) follow the figure.

**Physics.**
- (A8.17) (`a8-eq-catastrophe`) with q_min = q_D/N_1, T = 0, spherical zone:
  ΔE^NA_c(N_1) = −α_Fħω_LO (2/π)[arctan(q_D/q_pol) − arctan(q_D/(N_1 q_pol))];
  ΔE^A_c(N_1) = −α_Fħω_LO (2/π)(q_pol/q_D)(N_1 − 1).
- Contribution per logarithmic interval, in units of α_Fħω_LO, zero for q > q_D:
  nonadiabatic (2/π) q q_pol/(q² + q_pol²); adiabatic (2/π) q_pol/q. Both filled between q_min and q_D.
- Derived readouts: converged nonadiabatic shift −α_Fħω_LO (2/π) arctan(q_D/q_pol); N_1 for 90% of it
  q_D/(q_pol tan(0.1 arctan(q_D/q_pol))); light cone resolved at N_1 = q_D/q_lc; adiabatic shift cut at q_lc
  −α_Fħω_LO (2/π) q_pol (1/q_lc − 1/q_D); share of the nonadiabatic shift from q < q_lc:
  arctan(q_lc/q_pol)/arctan(q_D/q_pol).
- Check values: GaAs α_Fħω_LO = 2.58 meV, converged −2.55 meV, N_1(90%) = 278, ΔE^A_c = −3.76 meV at N_1 = 100 and
  −37.9 meV at 1000, light cone resolved at N_1 = 1.8×10⁴, adiabatic cut at q_lc −0.685 eV, share 1.6×10⁻³;
  GaN 40.9 meV, −39.2 meV, 103, cut −16.8 eV, share 1.0×10⁻³; MgO 153.5 meV, −147.5 meV, 107, cut −114 eV,
  share 5.7×10⁻⁴.

**What to draw.** Two panels in `viz-grid two`.
- (a) |ΔE_c| (meV) against N_1, log–log, N_1 from 4 to 1000 (x label `grid size N_1 (points per direction)`); y
  limits fixed per material (v1: 10^floor(log10(0.02 α_Fħω_LO)) to 10^ceil(log10(1.2 |ΔE^A_c(1000)|))).
  Nonadiabatic curve (`electron`, solid, label `nonadiabatic (ħω_{LO} kept)`); adiabatic curve (`c5`, dashed,
  label `adiabatic (Allen–Heine)`); horizontal line at α_Fħω_LO (`muted`, text with its value); vertical line at
  N_1 = q_D/q_pol (`phonon`, text `q_{min} = q_{pol}`); cursor line at N_1 with points on both curves.
- (b) contribution per ln q against q: x log from 10⁻⁶ to 3 Å⁻¹ (label `q (Å^{-1})`), y linear 0 to 0.5 (label
  `|dΔE/d ln q| / α_Fħω_{LO}`); the adiabatic curve leaves the frame at small q (that is the divergence; clip it).
  Fills between q_min and q_D (adiabatic `c5`, alpha 0.10; nonadiabatic `electron`, alpha 0.22); both curves; vertical
  lines at q_pol (`phonon`), q_min (`muted`, dashed), q_D (`muted`); when the toggle is on, a band 10⁻⁶ < q < q_lc
  (`photon`, alpha 0.16, text `light cone`) and a line at q_lc.

**Controls (3).**
1. material, seg: GaAs | GaN | MgO.
2. `N_1`, range 4–1000, logarithmic.
3. toggle: "mark the light cone q < q_lc".

**Readouts.** $\alpha_{\mathrm F}$ (with "outside weak coupling" when > 1); $q_{\mathrm{pol}}$ (Å⁻¹) and
$1/q_{\mathrm{pol}}$ (Å); $q_{\mathrm{lc}}$ (Å⁻¹); $\Delta E^{\mathrm{NA}}_{\mathrm c}(N_1)$ and
$\Delta E^{\mathrm A}_{\mathrm c}(N_1)$ (meV); converged nonadiabatic shift (meV); $N_1$ for 90%; $N_1$ that resolves the
light cone; adiabatic shift cut at $q_{\mathrm{lc}}$ (eV); share of the nonadiabatic shift from $q\lt q_{\mathrm{lc}}$ (%).

**Initial state.** GaAs, N_1 = 200, light cone on: the nonadiabatic curve is near its plateau (86%), the adiabatic
one already three times larger and still rising.

**Caption.**
- *What to notice.* The adiabatic shift grows linearly with N_1; the nonadiabatic one levels off near −α_Fħω_LO (at 0.96–0.99 of it,
  because the zone is finite) once q_min is well below q_pol (half of it is reached at q_min = q_pol). The light cone sits 400 to 1200 times below
  q_pol and holds about 10⁻³ of the converged shift: the phonon energy ħω_LO cures the divergence, not the speed of
  light.
- *Marks.* (a): nonadiabatic and adiabatic |ΔE_c| against the grid size, the polaron shift α_Fħω_LO, the grid where
  q_min = q_pol, the cursor. (b): contribution per logarithmic interval of q, so that equal areas mean equal
  contributions; filled areas are what the grid includes; markers q_pol, q_min, q_D; rose band: the light cone.
- *Try.* Increase N_1 from 4 to 1000 in GaAs: the nonadiabatic shift settles near −2.5 meV, the adiabatic one reaches
  −37.9 meV. Switch to MgO: α_F ≈ 1.7, outside weak coupling, and convergence needs only N_1 ≈ 100 because q_pol is
  larger. Turn the light cone on and compare its width with the region where the integrands live.
- *Model.* (A8.17) for the conduction-band minimum of one parabolic band coupled to one dispersionless LO phonon by
  (A8.3), T = 0, Rayleigh–Schrödinger, spherical zone of radius q_D = (6π²/Ω_V)^{1/3}, grid cutoff q_min = q_D/N_1.
  Parameters of the table (literature values, illustrative). Cite `a8-eq-catastrophe`, `a8-eq-power`, `a8-eq-gfr`.

**What must change from v1.**
- id `viz-a-temperature-frohlich` → `viz-a8-catastrophe`; classes `pf-*` → `a8-catastrophe-*`; own script.
- Notation: `N_q` → `N_1` (in v2 N_q is the number of wavevectors of the grid); `q_p` → `q_{pol}`; `q_c` → `q_{lc}`;
  `α` → `α_F`; readout keys `ΔE^{NA}(N_q)`, `ΔE^{A}(N_q)` → `ΔE^{NA}_c(N_1)`, `ΔE^{A}_c(N_1)`; curve label
  "adiabatic AH (G.167)" → "adiabatic (Allen–Heine)"; y label `|ΔE_{CBM}|` → `|ΔE_c|`.
- Rewrite the caption in the four-part structure (v1 cites G.166/G.167 and EPP §7.2 in source form).
- v1's physics (functions `NA`, `AD`, `matParams`, the per-ln-q integrands and the readouts) is correct and can be
  kept.
