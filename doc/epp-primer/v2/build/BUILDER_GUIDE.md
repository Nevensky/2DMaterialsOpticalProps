# Primer v2 · Builder guide for section and figure authors

Read this file, `v2/NOTATION.md` (binding) and your entry in `v2/build/SECTION_BRIEFS.md` before you write anything.

## 0. What v2 is

The primer explains, from the ground up, the first-principles theory of electron–phonon coupling, following
F. Giustino, *Rev. Mod. Phys.* **89**, 015003 (2017) ("the review", Part A). It then explains the EPP manuscript,
"Theory of the Electron-Phonon-Photon coupling" (draft v2), which adds the full electromagnetic field (Part B). A bridge (X.1, X.2) connects the two.

- **Reader.** A physics or chemistry graduate student. They know quantum mechanics (including the harmonic
  oscillator and first-order perturbation theory), classical electrodynamics, Bloch's theorem and band structures, and
  Fourier series. They have **not** read the two sources. They know second quantization, Green's functions, Dyson
  equations, RPA and GW only vaguely. Teach those in background boxes.
- **v2 is a rewrite of v1** (`doc/epp-primer/sections/*.html`). v1 keeps the source notations and is the content
  reference: its physics, coverage, numbers and figures survive unless your brief says otherwise.
  v2 changes three things:
  1. one notation (`NOTATION.md`), with the sources' symbols only in footnotes;
  2. self-contained, step-by-step derivations;
  3. geometric figures that show what the indices, vectors and matrix elements mean.
- **The model to imitate** is `v2/sections/A01-lattice.html` with `v2/sections/A01/fig-*.html`, reviewed twice by
  physics, naive-reader, notation and visual reviewers. Read A.1 in full before you write. Match its structure, density,
  tone and markup.
- **Errors in the sources.** v1's Z sections (`../sections/Z02-errata.html`, `../sections/Z03-corrections.html`) and
  NOTATION §12.4 list known errors. v2 presents the correct result in the main text and states the correction once, in a
  footnote.

## 1. Files and ids

| Item | Pattern | Example |
|---|---|---|
| section file | `v2/sections/<ID>-<slug>.html` (slug from v1) | `v2/sections/A02-coupling.html` |
| figure files | `v2/sections/<ID>/fig-<name>.html`, included by `<!-- include: <ID>/fig-<name>.html -->` | `A02/fig-overlap.html` |
| figure specs | `v2/sections/<ID>/FIGS.md` (written by the section author, read by figure authors) | |
| section element | `<section class="chapter" id="<v1 id>" data-num="A.2" data-title="…">` | `id="a-coupling"` |
| prefix | `a2`, `a9`, `x1`, `x2`, `b1`, …, `b10` | |
| heading ids | `<prefix>-<slug>` | `a2-anatomy` |
| equation ids | `<prefix>-eq-<slug>`; exported ids are fixed by the brief | `a2-eq-g` |
| equation tags | `\tag{A2.1}`, `\tag{X1.1}`, `\tag{B10.4}`, numbered in reading order | |
| figure ids | `viz-<prefix>-<name>` (name = file name without `fig-`) | `viz-a2-overlap` |
| footnote ids | `fn-<prefix>-<slug>` | `fn-a2-g` |
| glossary ids | `gl-<prefix>-<slug>`; the ids in NOTATION §12 are fixed | `gl-a2-g` |
| CSS classes in figures | prefixed with the figure: `a2-overlap-…` | |

X.1 starts with the bridge header and B.1 with the Part B header, before the `<section>`. Copy them from v1
(`<div class="part-head" id="part-x" data-part-title="Bridge · …">` and `id="part-b"`) and rewrite their text
in v2 style, as A.1 does with `part-a`.

Never edit another section's files, the build tools, `NOTATION.md` or this guide. If you need a symbol that NOTATION
does not fix, choose one that follows its rules, define it, and report it in your final answer.
Do not commit. The orchestrator commits.

## 2. Section skeleton, in order

1. `<header class="ch-head">` with:
   - an eyebrow: `<p class="eyebrow"><span class="num">A.2</span><span>topic · topic · topic</span></p>`;
   - the `<h2>` title;
   - a `<p class="lede">` of 2–4 sentences: what the section does and why.
2. A **"What you need"** paragraph: earlier sections, with links, and the background the reader should have, naming the
   background boxes that teach it.
3. The body: `<h3 id>` subsections with `<h4 id>` subsubsections. In each topic, follow this order:
   1. the physical picture, in words;
   2. the simplest concrete model, with a figure;
   3. the general formalism;
   4. what it means and where it goes next.
4. Callouts where they help:
   - `intuition` for a physical picture;
   - `caveat` for what the approximation leaves out;
   - `connect` for "Where these objects go next" (one at the end of the section);
   - `source` for a short note on how the sources present it, the only callout allowed source notation.
5. **Check yourself**: a `details.derivation` with 2–3 short questions and worked answers (see the end of A.1).
6. **Notes**: `<div class="notes"><h4>Notes</h4><ol class="footnotes">…</ol></div>`, in order of first reference.
7. **Glossary**: `<details class="glossary-box"><summary>Symbols defined in this section</summary><dl class="glossary">…</dl></details>`.

## 3. Markup, exactly as in A.1

```html
<div class="eq" id="a2-eq-g">$$ g_{mn\nu}(\bk,\bq)=\langle\psi_{m\bk+\bq}|\Delta_{\bq\nu}V^{\mathrm{KS}}|\psi_{n\bk}\rangle \tag{A2.5}$$</div>
<dl class="where">
  <dt>$g_{mn\nu}(\bk,\bq)$</dt><dd>electron–phonon matrix element … (eV)</dd>
</dl>

<ol class="steps">
  <li><p>What is done, with the formula.</p><span class="why">Why this step is allowed or why it is taken.</span></li>
</ol>

<details class="derivation"><summary>Derivation: …</summary><div class="d-body"> … </div></details>
<details class="background"><summary>Green's functions <span class="skip">skip if familiar</span></summary><div class="b-body"> … </div></details>
<div class="callout intuition"><p class="callout-title">Physical picture: …</p><p>…</p></div>

<sup class="fn"><a href="#fn-a2-g">3</a></sup>                       footnote marker (numbers are reassigned by position; write them in order)
<li id="fn-a2-g">Giustino: Eq. (38) … EPP: Eq. (…) …</li>         footnote, inside ol.footnotes

<a class="sym" href="#gl-a1-e">$e_{\kappa\alpha,\nu}$</a>          symbol link: opens the glossary entry as a popover
<a class="eqref" href="#a4-eq-glr"></a>                              equation of any section: shows "(A4.7)"
<a class="figref" href="#viz-a1-diatomic"></a>                       figure of any section: shows "Figure A1.4"
<a href="#a-polar">Section A.4</a>                                   section link

<dt id="gl-a2-g" data-first="a2-anatomy">$g_{mn\nu}(\bk,\bq)$</dt>
<dd>meaning, units, where defined. <span class="gl-src">Giustino: $g_{mn\nu}(\bk,\bq)$ (Eq. 38). EPP: $g_{nm\alpha}(\bk,\bq)$, with $N^{-1/2}$ inside and shifted labels.</span></dd>
```

- `data-first` is the id of the heading or element where the symbol is defined.
- Tables go in `<div class="tbl"><table>…</table></div>`.
- In X.2 only, a source-symbol dictionary table may be wrapped in `<!-- srcmap -->…<!-- /srcmap -->`, which the linter skips.
- **Phones.** The page must work at 390 px. Split long display equations with `\begin{gathered}…\\…\end{gathered}` or
  `\begin{aligned}`. Keep each line to about 40 rendered characters. A.1 shows the pattern.
- Use the v2 macros: `\ee \ii \dd \bk \bq \br \bR \bG \bu \be \bA \bE \bj \bP \btau \OV \OS \OVSC \OSSC \OBZ \rme \ion \ph \kB \tPhi`.
  Labels are `\mathrm{…}`, never `{\rm …}`.

## 4. Writing rules (A.1 lessons included)

1. **Define every symbol at first use**, in words, with its range and units. A display equation that introduces symbols
   is followed by a `dl.where`. A symbol whose home is an earlier section (NOTATION §12) is reminded in a few words at
   its first use in your section and linked with `a.sym` to its home glossary entry. Do not define it again.
2. **Define every term before using it**: crystal momentum, adiabatic, screened, vertex, retarded, on-shell,
   quasiparticle, longitudinal, transverse, and so on. Bold the term at its definition.
3. **Derivations motivate every step.** Use `ol.steps`, one step per idea, each with a `span.why`. Never write "it can be
   shown", "it follows that", "simply", "obviously", "clearly". If a derivation is long, put it, complete, in a
   `details.derivation` and state the result and its meaning in the main text.
4. **Picture first.** Each topic starts with the physical picture and the simplest model, with a figure, before the
   general formula.
5. **Numbers.** Give realistic magnitudes with units: eV, meV, Å, amu, K, cm⁻¹ where customary. Name the material.
   Use NOTATION §7 for conversions. Check the dimensions of every key formula, and give units in the `where` list.
6. **Source notation only in footnotes** and in `span.gl-src`. The same goes for source equation numbers. For every key
   result, a footnote gives the source equation numbers and the source symbols, so a reader can cross-check. Credit
   ideas in the main text as "the review" or "the EPP manuscript", and authors by name (Fröhlich, Bardeen, Pines, Allen,
   Heine, Hedin, Baym, …).
7. **SI footnote** on each key equation with $e$, $c$ or a Coulomb kernel (NOTATION §7).
8. **Corrections.** The main text shows the correct form. One footnote says what the source prints and why it is changed.
9. **Plain, exact English.** Short sentences, one idea per paragraph, active voice. Bold only the term being defined. No
   hype and no rhetorical questions. Use "we" sparingly.
10. **Cross-references.** Link sections freely. Link equations and figures of other sections only through the ids in
    NOTATION §12.6 and the briefs' **exports** lists, with `a.eqref` and `a.figref`. Use `a.eqref` for your own
    equations too, so numbers never go stale.
11. **Length.** Write what self-containedness needs, typically 45–90 KB of HTML for the section file, figures excluded.
    Do not pad.
12. **Background boxes** teach material the reader may know, closed by default. The main text never depends on them.
    Assumed and taught in boxes:
    - second quantization;
    - Green's functions and their poles;
    - Dyson equations;
    - Matsubara and retarded functions;
    - RPA, GW and Hedin's equations;
    - gauge freedom;
    - Lagrangian and Hamiltonian field theory.

## 5. Figures

Each figure is one file:

```html
<figure class="viz" id="viz-a2-overlap">
  <p class="viz-title">Title in sentence case <span class="tag">interactive</span></p>
  … divs for plots, canvases, controls, readouts …
  <figcaption><p class="cap-key"><strong>What to notice.</strong> 1–2 sentences.</p><div class="cap-more">
    <p><strong>Marks.</strong> what each colour, line, arrow and label means.</p>
    <p><strong>Try.</strong> 2–3 concrete things to do with the controls, and what happens.</p>
    <p><strong>Model.</strong> the equations used (with <a class="eqref" …>), parameter values and their source, the material it mimics.</p>
  </div></figcaption>
</figure>
<script>
(() => {
  const fig = document.getElementById('viz-a2-overlap');
  if (!fig || !window.EPP) return;
  …
})();
</script>
```

- **Caption.** About 200–260 words in total. The figure number is added automatically. A title that is purely
  geometric (no controls) uses `<span class="tag">diagram</span>`.
- **Geometric figures.** Draw the objects the equations talk about wherever an index, a vector or a matrix element has
  a spatial meaning. Examples:
  - cells, atoms and displacement patterns;
  - wavevectors on a Brillouin zone, Fermi circles shifted by $\bq$;
  - band diagrams with transition arrows;
  - field lines and polarization charges;
  - the polarization triad $\be_1,\be_2,\be_3=\hat\bq$;
  - a dipole sheet;
  - the light cone drawn to scale.

  A figure comes before the formalism that uses it. Label in the primer's notation: bold vectors (`drawRich` with
  `{bold:true}`), upright labels, italic variables.
- **Ported v1 figures.** Start from the v1 code: the `<figure … id="viz-…">` element and the `<script>` after it, in the
  v1 section file named in your brief. Keep its physics and interactions unless they are wrong. Then:
  - rename the id and the classes;
  - convert every symbol to v2 notation (labels, readouts, legends, aria-labels and strings);
  - rewrite the caption in the structure above;
  - check every formula in the code against the section's equations.
- **Notation in figures.** As in the prose:
  - no source symbols;
  - $\sigma$ is reserved for the conductivity;
  - readouts in HTML use MathJax `$…$`;
  - canvas text uses `EPP.drawRich` rich text (`ω_{LO}`, `q^2`).
- **Physics.** Compute from the section's equations, closed form or light numerics, under 15 ms per redraw. State the
  model and parameters. The initial state already shows the key phenomenon.
- **Quality bar.**
  - Axes are labelled, with units.
  - Fixed limits, so curves do not jump.
  - 2–5 controls with units.
  - No console errors and no NaN.
  - No external requests.
  - Works at 390 px with no horizontal overflow: controls below plots, `touch-action: none` only on canvases that are dragged.
  - Keyboard access and `aria-label`s.
  - Colours read at draw time with `EPP.color(...)` (themes change).
  - Calm animations through `EPP.anim`, which honours reduced motion.
- **Library.** Read `v2/build/lib.js` (window.EPP) and the API notes in `doc/epp-primer/build/BUILDER_GUIDE.md`
  ("EPP library API", "Interactive figure quality bar"). v2 adds:
  - `EPP.eigHerm(re, im)`, a complex Hermitian eigensolver returning ascending `{vals, vecs:[{re,im}]}`;
  - `EPP.fixPhase(v)`;
  - the colour `'c5-t'`;
  - `drawRich`/`measureRich` option `{bold:true}`, and the Plot text-layer `bold`;
  - range `aria-valuetext`;
  - the Play/Pause label of `EPP.anim`.

  Good models:
  - `A01/fig-diatomic.html`: plot, animation and readouts;
  - `A01/fig-bz.html`: a geometric construction;
  - `A01/fig-forces.html`: draggable geometry;
  - `A01/fig-zeropoint.html`.
- **Test every figure:** `bash v2/build/figtest.sh v2/sections/<ID>/fig-<name>.html`.
  - `errors` must be `[]`. The MathJax warning about `lazyAlwaysTypeset` is expected.
  - Then open both screenshots it prints (desktop light, phone dark) with the Read tool. Fix overlaps, clipped or
    unreadable labels, empty panels and anything that looks wrong. Repeat until clean.

## 6. Section test and self-check

`bash v2/build/sectest.sh <ID> [--shots]` does the following:
- assembles the conventions, A.1, your section and stubs for the other sections;
- lints your files (0 problems required);
- runs the page checker. Errors must be empty except **EXTERNAL** links into other sections, and each of those must be
  an id from NOTATION §12 or a brief's exports list.

Then:
- Run `python3 v2/build/lint.py --inventory <scratch>/inv.md` and read the warnings for your section: symbols with no
  glossary entry.
- Read your section top to bottom as the reader in §0. Is every symbol and term defined before use? Does every step
  say why? Is every number given with units?

## 7. Final answer of an agent

Return a short report:
- the files written;
- the equation tags and exported ids;
- new symbols not in NOTATION;
- anything you could not resolve.
