# Builder guide: Electron–Phonon–Photon Primer

You are writing ONE section fragment of a long, single-page pedagogical web primer.
Part A explains Giustino, Rev. Mod. Phys. 89, 015003 (2017). The Bridge connects it to the EPP manuscript.
Part B develops the EPP formalism ("Theory of the Electron-Phonon-Photon coupling").
Readers are physics grad students / researchers (computational condensed matter, optics of 2D materials: RPA, BSE, W).
The page must be rigorous AND pedagogical: the physical picture first, then the equations, then the interactive figure that makes the equations tangible.

The scaffold (CSS, fonts, MathJax, the `EPP` JS library) already exists. Your fragment is concatenated into
`<main>` in filename order. Do NOT write `<html>`, `<head>`, `<body>`, `<style>` or `<link>` tags, and do not load external scripts.

## Files
- Scaffold CSS: `build/head.html` (read it to see the available classes). Library: `build/lib.js` (read it, the API is below).
- Study notes (your primary source): `notes/G1..G6.md` (Giustino), `notes/E1..E9.md` (EPP), `notes/SYNTHESIS.md` (bridge, notation dictionary, outline).
- PDFs (authoritative for equations: check any equation you show against the page image with Read(file_path, pages="a-b")):
  - Giustino: `<giustino2017.pdf (not in the repo)>` (journal page 015003-N = PDF page N)
  - EPP: `<EPPv2.pdf (not in the repo)>`
- Example fragment showing every component: `test/sections/A01-demo.html`.

## Fragment skeleton
```html
<section class="chapter" id="a-phonons" data-num="A.2" data-title="Short TOC title">
  <header class="ch-head">
    <p class="eyebrow"><span class="num">A.2</span><span>Giustino §III.A · EPP §1</span></p>
    <h2>Full chapter title</h2>
    <p class="lede">One or two sentences: what the reader will understand by the end.</p>
  </header>

  <p>Running text… inline math $\omega_{\bq\nu}$ …</p>
  <div class="eq">$$ ... \tag{G.38} $$</div>

  <figure class="viz" id="viz-a-phonons-dispersion">
    <p class="viz-title">Title of the interactive <span class="tag">interactive</span></p>
    <div class="viz-grid two"><div class="plotA"></div><div class="plotB"></div></div>
    <div class="ctrls"></div>
    <div class="ro"></div>
    <figcaption><strong>What to notice.</strong> …</figcaption>
  </figure>

  <details class="derivation"><summary>Derivation: from (G.36) to (G.38)</summary><div class="d-body">…</div></details>
  <div class="callout intuition"><p class="callout-title">Physical picture</p><p>…</p></div>
  <div class="callout connect"><p class="callout-title">Bridge to EPP</p><p>…</p></div>
  <div class="callout caveat"><p class="callout-title">Subtlety</p><p>…</p></div>
  <div class="callout source"><p class="callout-title">Note on the source</p><p>…</p></div>
</section>
<script>
(() => {
  // all JS for this section, scoped in an IIFE; query only inside this section's figures
})();
</script>
```
Part headers (only where the outline says so) are a `<div class="part-head" id="part-b" data-part-title="Part B · EPP formalism">` with `<p class="part-label">`, `<h2>`, `<p>`.

## Content rules
- Section ids, numbers and titles come from your assignment. Every element id you create must start with your section id (e.g. `viz-a-phonons-...`) so nothing collides.
- Length: 900–2000 words of prose per section plus equations. Dense but readable. Short paragraphs. Use `<h3>` sub-headings to structure.
- Equations: transcribe faithfully from the PDF and keep the paper's notation. Tag equations with their source number: `\tag{G.38}` for Giustino Eq. (38), `\tag{E.112}` for EPP Eq. (112). If you simplify or combine, don't tag it with a source number; say "cf. (E.112)" in the text instead. Never invent a source equation number.
- Where the papers differ in units (Giustino mostly SI-like with e²/4πε₀ or atomic units, EPP probably Gaussian), say which convention a displayed equation uses.
- MathJax: standard LaTeX + AMS. Defined macros: `\bk \bq \br \bR \bG \bA \bE \bB \bj \bJ \bP \bQ \bu \be \bp \bZ \dd \ii \Tr \ket{} \bra{} \braket{}{} \mel{}{}{} \hc`. Do not define new macros. Break long display equations with `\begin{aligned}…\end{aligned}` so they fit ~760px (phone users can scroll inside `.eq`). Use `$$…$$` only inside `<div class="eq">` (or inside `.d-body`/callouts). Never put a literal `$` in prose.
- Use the particle color classes in prose sparingly (they use high-contrast text variants of the hues): `<span class="t-e">electron</span>`, `<span class="t-p">phonon</span>`, `<span class="t-g">photon</span>`.
- Callouts: `intuition` = physical picture, `connect` = Giustino↔EPP links (at least one per section), `caveat` = approximations/subtle points, `source` = typos/inconsistencies in the source (EPP is a draft; flag, don't propagate).
- Don't promise things the page doesn't do. Don't write meta text like "In this section we will…" more than once. No emoji. Plain, direct sentences. No em-dash asides.
- Cite precisely: "Giustino §IV.B.3, Eq. (G.79)" or "EPP §5.2.1, Eq. (E.212)".

## Visual language (fixed across the page)
- electron = blue (`'electron'`), solid line with arrow
- phonon = amber (`'phonon'`), zigzag line
- photon (transverse) = rose (`'photon'`), wavy line; Coulomb / longitudinal photon = rose dashed (`'coulomb'`); screened Coulomb W = rose double dashed (`'screened'`)
- extra series colors: `'c4'` (teal), `'c5'` (violet), `'muted'`, `'ink'`. This five-hue order (electron, phonon, photon, c4, c5) has been run through a colour-vision-deficiency validator. `'c4'` next to `'photon'` is only borderline-distinguishable for deutan viewers in dark mode, so any `'c4'` series MUST also carry a dash pattern or a direct text label.
- Use these semantic colors in plots: e.g. phonon branches amber, photon/light line rose, electron bands blue. Use dashes to distinguish further series. Never hard-code hex colors; always pass color names (they follow the light/dark theme).
- Chart rules: at most 5 coloured series per plot (fold the rest into small multiples or a second panel); ONE y-axis per plot (never dual axes: use two panels); lines 2px; a legend whenever there are ≥2 labelled series, plus direct `text` labels for ≤4 series where space allows; recessive grid. Annotation text is automatically drawn in the text variant of its hue. Heat maps use a single-hue sequential ramp (or diverging with `color2` for signed quantities); always give a colorbar label with units.

## EPP library API (window.EPP)
- `new EPP.Plot(el, {xlabel, ylabel, xlim:[a,b], ylim:[a,b], aspect:0.58, height, xname, hoverX:fn, hoverY:fn, xfmt, yfmt, xticks:[v|{v,label}], yticks, logx, logy, legend:true, hover:true, margin:{l,r,t,b}, ariaLabel})`
  - `plot.draw(layers, {xlim, ylim, xlabel, ylabel}?)` redraws (cheap; call it on every slider change). Layers:
    - `{type:'line', x:[], y:[], color, width, dash:[6,4], label, hover:false, breakJump}` (breaks at NaN/Infinity and jumps > 3× y-range: return NaN to lift the pen at poles)
    - `{type:'fill', x, y1, y2?, base?, color, alpha, label}`
    - `{type:'points', x, y, color, r, hollow}`
    - `{type:'hline'|'vline', at, color, dash, width, text, textAt:'left', textDy}`
    - `{type:'band', x0, x1, color, alpha, text}` vertical shaded region; `{type:'hband', y0, y1, ...}`
    - `{type:'text', x, y, text, color, size, align:'left'|'center'|'right', dx, dy, baseline}`
    - `{type:'arrow', x0, y0, x1, y1, color, width, head}`; `{type:'segments', segs:[[x0,y0,x1,y1],...], color}`
    - `{type:'heatmap', x0, x1, y0, y1, nx, ny, z, vmin, vmax, color:'electron', color2?:'photon' (diverging: color2 = low end, neutral at the midpoint, color = high end), gamma, label:'A(k,E)', colorbar:'A (eV^{-1})', cfmt}`: `z` is a Float64Array of length nx*ny, row-major with j running over y from y0 to y1 and i over x from x0 to x1 (`z[j*nx+i]`). Give the plot `margin:{r:60}` so the colorbar fits. Pass a NEW array object when values change (the rendered image is cached per array). Keep nx*ny ≲ 60k for fast redraws. Hover shows the value (set `yname` in the Plot options).
  - Labels/axis titles accept rich text: `ω_{LO}`, `q^2`, `ħω_0`, `ε_∞` (sub `_x`/`_{xy}`, super `^x`/`^{xy}`).
  - `plot.onPointer = (x, y, type) => {}` (type 'down'|'move'|'up') for click/drag selection in data coordinates.
- `EPP.controls(el, specs, onChange)` → `{state, set(key, v, silent), refresh()}`. Specs: `{key, type:'range', label, min, max, step, value, fmt:v=>string, log:true?}`, `{key, type:'select', label, options:[[value,label],...], value}`, `{key, type:'toggle', label, value}`, `{key, type:'seg', label, options:[[value,label],...], value}`, `{key, type:'button', label, onClick}`. onChange(state) fires once at start and after every change (rAF-batched). Labels accept rich text. Values from 'seg'/'select' are strings.
- `EPP.readout(el).set({'ω_{LO}': '36.2 meV', ...})` small live key/value line.
- `EPP.canvas(el, {aspect, height, minHeight, maxHeight})` → `{canvas, ctx, w, h, clear(), onResize(cb)}` for custom drawings (lattice animations, field arrows). Redraw in `onResize`, which also fires on theme change. Read colors with `EPP.color('phonon')`, `EPP.alpha('photon', 0.3)`, `EPP.css('--fg')` at draw time (never cache colors across frames, the theme can change).
- `EPP.anim(figureEl, (t, dt) => draw(t), {buttonHost: el, autoplay})` → Play/Pause button, pauses offscreen, starts paused under reduced motion. Use for mode animations. Keep animations calm.
- `EPP.fd.draw(el, {w, h, label, items})` Feynman diagrams (SVG, theme-aware). Put diagrams inside `<div class="fd">`, optionally each in `<figure>` with `<figcaption>`; use `<span class="op">=</span>` / `+` between diagrams. Items:
  - lines `{t:'fermion'|'fermion2'|'phonon'|'phonon2'|'photon'|'photon2'|'coulomb'|'screened'|'plain', a:[x,y], b:[x,y], bend:0, arrow:true, label:'k+q', lpos:'right', ldist, color, amp, wl, width}`. `*2` = dressed propagator (double line for fermions, heavy stroke for zigzag/wavy). `bend` > 0 bulges to the screen-left of a→b (upwards for a left-to-right line); labels sit on the screen-left side by default, `lpos:'right'` for the other side.
  - `{t:'vertex', at:[x,y], r}`, `{t:'blob', at, r, label:'Π', color}`, `{t:'text', at, text, color, size, anchor, upright}`.
- Numerics: `EPP.linspace(a,b,n)`, `EPP.sample(f,a,b,n)`, `EPP.clamp`, `EPP.bose(E,kT)`, `EPP.boseT(E_eV,T_K)`, `EPP.fermi(E,mu,kT)`, `EPP.const` (hbar_eVs, kB_eV, c_ms, e_C, me_kg, amu_kg, eps0, a0_A, Ry_eV, Ha_eV, meV_per_cm1, THz_per_meV), complex `EPP.C` (c, add, sub, mul, scale, div, abs, arg, conj, inv, sqrt, exp, log; numbers are [re, im]), `EPP.eigSym2(a,b,d)`, `EPP.eigHerm2(a,d,[bre,bim])`, `EPP.jacobi(A)`, `EPP.bisect(f,a,b)`, `EPP.trapz(x,y)`, `EPP.fmt(v,digits)`, `EPP.fix(v,decimals)`, `EPP.richHTML(s)`.

## Interactive figure quality bar
- Every figure computes real physics from the equations in the section (closed form or light numerics, < 15 ms per redraw). State the model and the parameter values in the caption. Prefer realistic material parameters with units (and say which material they mimic). If you use a toy model, label it as a model.
- Axes have labels with units. Set sensible fixed limits so curves don't jump around as sliders move (unless autoscaling is the point).
- Controls: 2–5 meaningful sliders/toggles with physical units in the label and value format. Initial state should already show the key phenomenon.
- The caption ends with "What to notice" guidance in one to three sentences.
- Add a readout of key derived numbers when useful (e.g. ω_LO, λ, T_c, ZPR).
- Figures must look right at 390px width (stacked) and 820px. Use `viz-grid two` / `side` for side-by-side panels (stack on phones).
- No console errors. Guard against NaN. No external requests.

## Self-check before you finish
1. Assemble and check your fragment in isolation:
   `python3 build/assemble.py test/<your-id>.html <dir containing only your fragment>` then
   `node build/check.js test/<your-id>.html test/shots-<your-id> --shots`
   (copy your fragment into an empty temp dir like `test/frag-<your-id>/` first). The JSON report must show no errors, no TeX errors, no overflow, no blank canvases.
2. Look at the desktop and phone screenshots (Read the PNGs) once and fix anything visibly broken (overlapping labels, clipped text, empty plot).
3. Re-verify every tagged equation against the PDF page image.
Run commands from the scratchpad directory: `doc/epp-primer`.


## Page-wide conventions (decided; follow them)
- Units: Part A uses Giustino's SI conventions and drops $N_p$ from $g$ (his Eq. 38 convention; note the $N_p$ inconsistency of G Eqs. 2, 57, 58, 61, 180 in a `source` callout only where relevant). Part B uses the EPP manuscript's Gaussian units, with Born charges $Z^*$ carrying the dimension of charge ($Z^*_{\rm EPP}=|e|Z^*_{\rm G}$). For numbers use $e^2/4\pi\varepsilon_0 = 14.3996$ eV Å, $e^2/\varepsilon_0=180.95$ eV Å, $\hbar c = 1973.27$ eV Å, $\hbar^2/2m_e = 3.80998$ eV Å², $k_B = 0.0861733$ meV/K, 1 eV = 8065.54 cm⁻¹, $\hbar^2/(\mathrm{amu}\,\text{Å}^2) = 4.1803$ meV.
- Renamed symbols, to avoid clashes between the papers. Giustino's phonon Green function is written $\mathcal D$ and his phonon self-energy $\Pi^{\rm ph}$. EPP's photon propagator is $D$ (with $D^0$, $D^{0l}$, $D^{0t}$, $D^{\rm e}$), its photon self-energy $\Pi_{jj}$ (with $\Pi^{\rm e}$, $\Pi^{\rm ion}$), and its phonon self-energy $\Omega^{\rm e}$. EPP's phonon propagator is $F$. Say once in your section if you rename relative to the PDF.
- Use retarded forms in visualizations: denominators $\omega^2-\omega_0^2+i\gamma\omega$.
- Literature numbers not printed in either paper must be flagged as illustrative (e.g. "(literature value, illustrative)").
- "Retardation" has two meanings (Giustino: frequency-dependent electronic screening; EPP: finite speed of light). Always say which one you mean.
- The list of known source pitfalls is in `notes/SYNTHESIS.md` §6. Read it in full before writing; never propagate those errors.
- Cross-link other sections with `<a href="#section-id">§B.6</a>`. Do not re-derive another section's main result; link to it.

## Section map (the whole page)
Part A · Electron–phonon theory (Giustino 2017). The Part A header goes at the top of A.1's fragment.
- A.1 `a-lattice` Lattice dynamics and phonons. Giustino §III.A, App. B; EPP §1.
- A.2 `a-coupling` The electron–phonon Hamiltonian and a century of $g$. Giustino §I–II, Eqs. 1–11, 38.
- A.3 `a-dfpt` DFPT: the vertex is screened. Giustino §III.B.
- A.4 `a-polar` Polar materials: Fröhlich, Verdi–Giustino and the LO–TO splitting. Giustino §II.A.3, §VI.A.3.
- A.5 `a-hedin-baym` Hedin–Baym: the lattice inside $W$. Giustino §IV.
- A.6 `a-phonon-se` Phonon self-energy: adiabatic vs nonadiabatic. Giustino §V.A, §VII.
- A.7 `a-electron-se` Electron self-energy: Fan–Migdal, Debye–Waller, kinks and satellites. Giustino §V.B.1–2, §V.B.5–6, §VIII.
- A.8 `a-temperature` Band gaps vs temperature, lifetimes and the Fröhlich catastrophe. Giustino §V.B.3–4, §IX.A, §X.A.
- A.9 `a-applications` Light as an external field, transport and superconductivity. Giustino §IX.B, §X.B, §XI, §XII.
Bridge · From electron–phonon to electron–phonon–photon. The Bridge header goes at the top of X.1's fragment.
- X.1 `x-light-cone` From instantaneous Coulomb to the full electromagnetic field.
- X.2 `x-dictionary` A dictionary between the two formalisms (concept table, notation, diagram glyphs).
Part B · The EPP formalism. The Part B header goes at the top of B.1's fragment.
- B.1 `b-foundations` Same lattice, bare vertex, electrostatic polar coupling. EPP §§1–3.
- B.2 `b-lagrangian` One Lagrangian for electrons, ions and fields. EPP §4.1.
- B.3 `b-dielectric` The polar crystal as a dielectric. EPP §4.2.
- B.4 `b-fields` Quantized fields: transverse photons, and Coulomb as the longitudinal photon. EPP §§4.3–4.4.
- B.5 `b-photon-se` Propagators and the RPA photon self-energy. EPP §§4.5–4.7, §5.
- B.6 `b-propagator` Charge conservation, the total propagator and polaritons. EPP §6.
- B.7 `b-phonon-se` The electrodynamic phonon self-energy and the LO–TO splitting. EPP §7.1.
- B.8 `b-frohlich` The electrodynamic Fröhlich vertex. EPP §7.2.
- B.9 `b-electrostatic` The electrostatic limit: recovering Giustino. EPP §7.3.
- B.10 `b-2d` Two-dimensional crystals: no LO–TO splitting at Γ. EPP §8.
- Z.1 `z-outlook` What EPP adds and what is still open (written separately).
