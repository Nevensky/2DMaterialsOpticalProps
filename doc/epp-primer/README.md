# Electron–Phonon–Photon Primer

Interactive lecture notes that explain the electron–phonon theory reviewed in
F. Giustino, *Electron-phonon interactions from first principles*,
Rev. Mod. Phys. **89**, 015003 (2017), and then develop the EPP formalism of the
manuscript *Theory of the Electron-Phonon-Photon coupling*, connecting the two.

Open `out/epp-primer.html` for the assembled page. It is a fragment meant for the
claude.ai artifact viewer, which wraps it in a document skeleton. It loads
MathJax 3.2.2 from cdnjs and fonts from Google Fonts.

## Contents

| Part | Sections |
|---|---|
| Part A · Giustino 2017 | A.1 lattice dynamics · A.2 the coupling *g* · A.3 DFPT · A.4 polar materials · A.5 Hedin–Baym · A.6 phonon self-energy · A.7 electron self-energy · A.8 temperature effects · A.9 optics, transport, superconductivity |
| Bridge | X.1 Coulomb as the longitudinal photon · X.2 dictionary between the formalisms |
| Part B · EPP | B.1 foundations · B.2 Lagrangian · B.3 dielectric · B.4 quantized fields · B.5 photon self-energy · B.6 total propagator · B.7 electrodynamic phonon self-energy · B.8 electrodynamic Fröhlich vertex · B.9 electrostatic limit · B.10 2D crystals |
| Closing | Z.1 what EPP adds and what is open · Z.2 notes on the sources (errata for the EPP draft and Giustino 2017) |

Equation tags `(G.n)` and `(E.n)` refer to the equation numbers of the two PDFs.
The PDFs are not in this repository.

## Layout

- `sections/` — one HTML fragment per section (prose, equations, figure markup and a scoped script).
  The assembler orders them Part A (`A*`), Bridge (`X*`), Part B (`B*`), closing (`Z*`).
- `build/head.html` — title, design tokens, CSS and the MathJax configuration.
- `build/lib.js` — the `window.EPP` helper library: theme-aware canvas plots (lines, fills, heat maps),
  controls, animation loop, Feynman-diagram SVG builder and small numerics.
- `build/hero.html`, `build/foot.html` — opening figure with reading paths; references and contents-rail script.
- `build/assemble.py` — concatenates everything into `out/epp-primer.html` and builds the contents rail.
- `build/BUILDER_GUIDE.md` — conventions every section follows (markup, notation, units, colours, chart rules).
- `notes/` — page-by-page study notes of both papers (G1–G6: Giustino, E1–E9: EPP) and `SYNTHESIS.md`
  (bridge table, notation dictionary, known pitfalls).

## Build and check

```sh
cd doc/epp-primer
python3 build/assemble.py                       # writes out/epp-primer.html

# optional checks (need Node and Playwright with Chromium, and a local MathJax copy):
mkdir -p vendor && (cd vendor && npm pack mathjax@3.2.2 && mkdir mathjax && tar -xzf mathjax-3.2.2.tgz -C mathjax --strip-components=1)
node build/check.js out/epp-primer.html shots --eager   # console/TeX errors, overflow, blank canvases
build/shoot_all.sh "$PWD/figs"                          # screenshots of every figure, per section
python3 build/textdump.py                               # text-only dump in reading order (out/primer-text.md)
```

`check.js` serves the local MathJax copy in place of the CDN and blocks web fonts, so screenshots use
fallback fonts.
