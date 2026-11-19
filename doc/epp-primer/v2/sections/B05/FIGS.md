# Section B.5 · figure specifications

Section file: `v2/sections/B05-photon-se.html` (id `b-photon-se`, prefix `b5`). Five figures, in reading order:

| # | File | Id | Kind | Tag |
|---|---|---|---|---|
| B5.1 | `fig-poles.html` | `viz-b5-poles` | new | interactive |
| B5.2 | `fig-diagrams.html` | `viz-b5-diagrams` | new (replaces v1's inline Feynman diagrams `#b-photon-se-fd-dyson` and `#b-photon-se-fd-pi`) | interactive |
| B5.3 | `fig-electrons.html` | `viz-b5-electrons` | port of `viz-b-photon-se-electrons` | interactive |
| B5.4 | `fig-sumrule.html` | `viz-b5-sumrule` | port of `viz-b-photon-se-sumrule` | interactive |
| B5.5 | `fig-gauge.html` | `viz-b5-gauge` | port of `viz-b-photon-se-gauge` | interactive |

v1 source of the three ports: `doc/epp-primer/sections/B05-photon-se.html`. Each port has its `<figure id="viz-b-photon-se-…">` block in the body and its code in the single `<script>` at the end of that file: `// ---------- Viz B: electrons, truncated band sums ----------` (electrons), `// ---------- Viz C: oscillator-strength sum rule ----------` (sumrule) and `// ---------- Viz A: para + dia = gauge invariant ----------` (gauge). The v1 inline diagrams are drawn at the top of the same script (`// ---------- Feynman diagrams ----------`); reuse their `EPP.fd.draw` item lists for B5.2.

**General rules for all five** (from `v2/build/BUILDER_GUIDE.md` §5 and `v2/NOTATION.md`):
- v2 notation only. In particular: phonon branch $\nu$ (never $\sigma$, which is reserved for the conductivity; the linter flags `σ` in script strings); atom $\kappa$ with mass $M_\kappa$ (never $m_a$); reduced mass $M_{\mathrm r}$ (never $\mu$, which is the chemical potential); cell volume $\Omega_{\mathrm V}$ (macro `\OV`, never $\Omega_{\mathrm{uc}}$); Born charges $Z^*_\kappa$ dimensionless, so the dipole is $eZ^*$; electronic dielectric constant $\varepsilon_\infty$ (never $\epsilon^\infty$); static constant $\varepsilon_{\mathrm s}$; ionic plasma frequency $\omega^{\mathrm{ion}}_{\mathrm p}$ (never $\Omega_p$); electronic plasma frequency $\omega_{\mathrm p}$; transverse and longitudinal optical frequencies $\omega_{\mathrm{TO}}$, $\omega_{\mathrm{LO}}$; photon self-energy $\Pi$ with the parts $\Pi^{\mathrm{dia}}$ (contact) and $\Pi^{\mathrm{para}}$ (bubble); time-ordered phonon line $F^{\mathcal T}_\nu$, retarded line $F^0_\nu$.
- Labels upright where they abbreviate words (`\mathrm`), $\mathrm i$ and $\mathrm e$ upright, vectors bold (`drawRich` with `{bold:true}`); HTML readouts in MathJax `$…$`; canvas text with `EPP.drawRich` (`ω_{TO}`, `ω_p^{ion}`).
- No source equation numbers (E.281, E.297′, …) anywhere in a figure: cite this section's equations by id with `<a class="eqref" href="#b5-eq-…"></a>` in captions; canvas legends use words.
- Colours read at draw time with `EPP.color(...)`: electrons `electron`, phonons and ions `phonon`, photons `photon`, extra series `c4`, `c5`, neutral `ink`/`muted`.
- Caption in the four parts *What to notice* (in `p.cap-key`) / *Marks* / *Try* / *Model*, about 200–260 words in total. CSS classes prefixed `b5-<name>-`. Controls below the plots; works at 390 px; keyboard access and `aria-label`s; no NaN, no console errors; retarded denominators with a damping rate, $\omega^2-\omega_0^2+\mathrm i\gamma\omega$ (NOTATION §8), never the time-ordered $+\mathrm i\eta$ for plotted spectra (except where B5.1 shows the time-ordered line on purpose).
- Test with `bash v2/build/figtest.sh v2/sections/B05/fig-<name>.html` and look at both screenshots.

**Equation ids of the section** (tags as of this writing; always link by id):
`b5-eq-Fdef` (B5.1), `b5-eq-FT` (B5.2), `b5-eq-Fret` (B5.3), `b5-eq-F0` (B5.4), `b5-eq-D` (B5.5), `b5-eq-D0T` (B5.6), `b5-eq-AL` (B5.7), `b5-eq-D0L` (B5.8), `b5-eq-V` (B5.9), `b5-eq-jop` (B5.10), `b5-eq-dyson` (B5.11), `b5-eq-Pidef` (B5.12), `b5-eq-dysonDe` (B5.13), `b5-eq-eps` (B5.14), `b5-eq-dia` (B5.15), `b5-eq-crystal` (B5.16), `b5-eq-jel` (B5.17), `b5-eq-Pipara` (B5.18), `b5-eq-Pie` (B5.19), `b5-eq-lorentz` (B5.20), `b5-eq-pion` (B5.21), `b5-eq-jion` (B5.22), `b5-eq-Piionpara` (B5.23), `b5-eq-Piiondia` (B5.24), `b5-eq-sum` (B5.25), `b5-eq-Piion` (B5.26), `b5-eq-epsir` (B5.27). Other sections: `a1-eq-dyn` (A1.18, the dynamical matrix with the phase on the partner atom's cell), `a4-eq-asr` (Born-charge sum rule), `b3-eq-eps` (the one-mode dielectric function of B.3), `b3-eq-lst`, `a4-eq-lot`.

---

## B5.1 · `fig-poles.html` · `viz-b5-poles` · new

**Title.** "Free boson lines: where the poles sit, and why temperature drops out" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. The retarded free line $F^0_\nu$ (B5.4) has both poles just below the real axis, at $\omega=\pm\omega_0-\mathrm i\eta$. Its imaginary part is odd in $\omega$: negative at $+\omega_0$, positive at $-\omega_0$. It does not change with temperature.
2. The time-ordered line $F^{\mathcal T}_\nu$ (B5.2) at $T=0$ has its pole at $+\omega_0$ below the axis and its pole at $-\omega_0$ above it. For $\omega>0$ it is the same function as the retarded line; for $\omega<0$ its imaginary part has the opposite sign.
3. At $T>0$ the time-ordered line picks up the occupation $n_{\mathrm B}$: its imaginary part becomes negative at both $\pm\omega_0$, with weight $1+2n_{\mathrm B}$, while its real part stays the same. This is why EPP's rule (use $T=0$ boson lines, then make the $\mathrm i\eta$ retarded) loses nothing for harmonic bosons.

**Include position.** Section `b5-phonon` ("The free phonon line"), right after (B5.4) and the paragraph that ends: "They differ only at negative frequencies, where the time-ordered line has its pole above the real axis and the retarded line below. Figure B5.1 shows both lines, their poles, and how the time-ordered line changes with temperature while the retarded one does not." The subsection "Temperature: the rule $n_{\mathrm B}=0$" follows the figure.

**Physics (v2 notation).** One phonon mode of energy $\hbar\omega_0$ (the label $\bq\nu$ dropped). Work with energies $E=\hbar\omega$, $E_0=\hbar\omega_0$, and a finite broadening $\Gamma=\hbar\eta$ (the infinitesimal kept finite for display; say so in *Model*).
- Retarded (B5.4): $F^0(E)=\dfrac{2E_0}{(E+\mathrm i\Gamma)^2-E_0^2}$ (1/meV).
- Time-ordered at temperature $T$ (B5.2): $F^{\mathcal T}(E)=(1+n)\dfrac{2E_0}{E^2-E_0^2+2\mathrm iE_0\Gamma}-n\dfrac{2E_0}{E^2-E_0^2-2\mathrm iE_0\Gamma}$, with $n=1/(\mathrm e^{E_0/k_{\mathrm B}T}-1)$, $k_{\mathrm B}=0.0861733$ meV/K. (The $2E_0\Gamma$ gives both lines peaks of half width $\approx\Gamma$ near $\pm E_0$.)
- Plot both normalized by the peak size of the retarded line, $1/\Gamma$: $y=\Gamma\,\mathrm{Im}\,F$ and $y=\Gamma\,\mathrm{Re}\,F$. Then the retarded imaginary part peaks at $\mp1$ at $E=\pm E_0$, the time-ordered one at $-(1+2n)$ at both. Abscissa $x=E/E_0\in[-2.5,2.5]$, fixed. Ordinate fixed at $[-6.5,\,1.8]$ (with the control ranges below, $1+2n\le5.8$).
- Check values: $\hbar\omega_0=33.3$ meV (GaAs TO phonon), $T=300$ K: $n=0.381$, $1+2n=1.76$; $T=1000$ K: $1+2n=5.24$; $T=600$ K: $1+2n=3.21$.

**What to draw.**
- (a) *The two lines on the real axis.* $y=\Gamma\,\mathrm{Im}\,F^0$ solid (`ink`), $y=\Gamma\,\mathrm{Im}\,F^{\mathcal T}$ dashed (`phonon`), optional real parts as thin lines (identical for both lines and for every $T$; draw one thin `muted` curve labelled "Re, both"). Vertical guide lines at $x=\pm1$ labelled $\pm\omega_0$ (canvas `±ω_0`). Shade the half $x<0$ very lightly with the label "negative frequencies".
- (b) *The complex frequency plane.* Axes $\mathrm{Re}\,\omega/\omega_0\in[-2,2]$ and $\mathrm{Im}\,\omega$ (schematic, the distance $\eta$ from the axis exaggerated and fixed on screen). Retarded poles: two filled dots (`ink`) at $\pm1-\mathrm i\eta$. Time-ordered poles at $T=0$: open circles (`phonon`) at $+1-\mathrm i\eta$ and $-1+\mathrm i\eta$. For $T>0$ add the two "thermal" poles of the second term of (B5.2), open circles at $+1+\mathrm i\eta$ and $-1-\mathrm i\eta$, with radius growing as $\sqrt{n}$ (area $\propto n$) and the first pair with area $\propto1+n$. Label the upper half plane "no poles of a retarded function" and the real axis "real ω".
- Phone: (a) above (b), controls below.

**Controls (4).**
1. Temperature $T$: 0–1000 K, step 10, default 300 K.
2. Mode energy $\hbar\omega_0$: 30–200 meV, step 0.1, default 33.3 meV (GaAs TO). (Lower values would push $1+2n$ off the fixed scale.)
3. Broadening $\hbar\eta$ as a fraction of $\hbar\omega_0$: 0.02–0.15, default 0.05 (readout also in meV).
4. Check box "show real part" (default off).

**Readouts.** $n_{\mathrm B}(\omega_0)$; $1+2n_{\mathrm B}$ (equal to $\coth(\hbar\omega_0/2k_{\mathrm B}T)$, show both); $k_{\mathrm B}T/\hbar\omega_0$; "peak of Im at $+\omega_0$: retarded $-1$, time-ordered $-(1+2n_{\mathrm B})$" in units of $1/\hbar\eta$.

**Initial state.** $T=300$ K, $\hbar\omega_0=33.3$ meV, $\hbar\eta/\hbar\omega_0=0.05$, real part off: the dashed peaks are already 1.76 times the solid ones and both dashed peaks point down, while the solid curve is odd.

**Caption key points.**
- *What to notice.* The solid retarded line never moves when $T$ changes, and its imaginary part is odd in $\omega$. The dashed time-ordered line agrees with it at positive frequencies when $T=0$; heating multiplies its peaks by $1+2n_{\mathrm B}$ and turns the peak at $-\omega_0$ negative. Only the imaginary part knows the temperature.
- *Marks.* Solid: $F^0$ <a class="eqref" href="#b5-eq-F0"></a>; dashed: $F^{\mathcal T}$ <a class="eqref" href="#b5-eq-FT"></a>; thin: real part, the same for both. In (b) filled dots are the poles of $F^0$, open circles those of $F^{\mathcal T}$, their areas proportional to the weights $1+n_{\mathrm B}$ and $n_{\mathrm B}$.
- *Try.* Set $T=0$: the curves coincide for $\omega>0$ and are mirror images in sign for $\omega<0$; only one time-ordered pole is above the axis. Raise $T$ to 600 K: the dashed peaks grow by 3.2, the solid curve and both real parts do not change. Raise $\hbar\omega_0$ at fixed $T$: $n_{\mathrm B}$ falls and the dashed line returns to its $T=0$ form.
- *Model.* One harmonic mode, (B5.2) and (B5.4), with the infinitesimal $\eta$ kept finite as a broadening; default mode the TO phonon of GaAs, 33.3 meV (literature value, illustrative); ordinates in units of $1/\hbar\eta$.

---

## B5.2 · `fig-diagrams.html` · `viz-b5-diagrams` · new

**Title.** "Dressing the photon: the Dyson series and the four pieces of the self-energy" `<span class="tag">interactive</span>` (a diagram with one highlight control; use `diagram` if no control is built).

**What the reader must learn.**
1. The dressed photon line $D$ is the free line $D^0$ interrupted any number of times by the self-energy $\Pi$: the Dyson equation (B5.11) is the sum of the chain $D^0+D^0\Pi D^0+D^0\Pi D^0\Pi D^0+\dots$
2. In the random-phase approximation $\Pi$ has four pieces, two for each kind of charge (B5.12): the electronic bubble (two current vertices joined by an electron–hole pair), the electronic contact term (two photon lines at one point, weighted by the electron density), the ionic bubble (current vertex, free phonon line $F^0_\nu$, current vertex: (B5.23)) and the ionic contact term (two photon lines at an atom).
3. The electrons can be resummed first into $D^{\mathrm e}$, and the ions then dress $D^{\mathrm e}$ (B5.13).

**Include position.** Section `b5-couplings` ("Four couplings, two for each kind of charge"), right after the paragraph: "Feynman diagrams draw these couplings as vertices: a paramagnetic vertex joins one photon line to a current of matter, and a contact vertex joins two photon lines at one point. Figure B5.2 shows how they dress the photon line." The self-consistent-field derivation of the Dyson equation (subsection `b5-dyson`) follows the figure.

**Physics.** No numerics. Glyphs must match the primer's diagram dictionary (X.2 `fig-glyphs`): photon line wavy in `photon` colour, dressed (full) propagator as a double wavy line (`photon2`), free phonon line as a zigzag in `phonon` colour, electron lines solid with arrows in `electron` colour, self-energy as a blob (`blob`). Current vertices as small filled squares (electron vertex in `electron`, ionic vertex in `c4`), the electronic contact vertex as a filled dot, the ionic contact vertex as a ring in `phonon` colour at an atom. Draw with `EPP.fd.draw` (SVG; see the item list in `v2/build/lib.js`), as v1 did.

**What to draw** (four rows of diagrams, each row an SVG that scales to the width; on a phone the rows wrap so that every equation sign stays with its term):
- (a) *Dyson equation.* $D$ (double wavy) $=$ $D^0$ (single wavy) $+$ $D^0$ — blob $\Pi$ — $D$. Labels $D$, $D^0$, $\Pi$.
- (b) *The chain.* $D^0$ $+$ $D^0\,\Pi\,D^0$ $+$ $D^0\,\Pi\,D^0\,\Pi\,D^0$ $+\dots$, with small blobs.
- (c) *The four pieces of $\Pi$*, each in a box with a short caption under it:
  1. electronic bubble: two photon stubs, two `electron` squares, an electron loop with labels $n\mathbf k$ and $m\,\mathbf k+\mathbf q$ (write $\mathbf k$, $\mathbf q$ bold if the SVG label helper allows; otherwise plain and say so in *Marks*); caption "$\Pi^{\mathrm{e,para}}$: current–current bubble";
  2. electronic contact: two photon lines meeting at a dot, label $e^2n_{\mathrm e}/m_{\mathrm e}$; caption "$\Pi^{\mathrm{e,dia}}$: contact term";
  3. ionic bubble: two photon stubs, two `c4` squares labelled $j_{\alpha\nu}$ and $j^*_{\beta\nu}$, joined by the phonon zigzag labelled $F^0_\nu$; caption "$\Pi^{\mathrm{ion,para}}$: vertex, phonon line, vertex";
  4. ionic contact: two photon lines meeting at a `phonon` ring, label $e^2\mathbf Z^*_\kappa\mathbf Z^{*\mathsf T}_\kappa/M_\kappa$ (canvas/SVG text `e^2 Z^* Z^{*T}/M_κ`); caption "$\Pi^{\mathrm{ion,dia}}$: contact term at the atoms".
- (d) *Electrons first* (shown when the toggle is on): $D^{\mathrm e}=D^0+D^0\,\Pi^{\mathrm e}D^{\mathrm e}$ and $D=D^{\mathrm e}+D^{\mathrm e}\,\Pi^{\mathrm{ion}}D$, with $D^{\mathrm e}$ drawn as a double wavy line marked e, $\Pi^{\mathrm e}$ as a blob in `electron` colour and $\Pi^{\mathrm{ion}}$ in `phonon` colour.

**Controls (2).**
1. Radio group "highlight": none / electronic bubble / electronic contact / ionic bubble / ionic contact (keyboard accessible; clicking a box in row (c) selects it too). The selected box is outlined, the others are dimmed, and a readout under row (c) shows, in MathJax, the formula and one sentence:
   - electronic bubble: "$\frac{\mathrm i}{\hbar c}\int\mathrm dt\,\mathrm e^{\mathrm i\omega t}\langle\mathcal T\hat j^{\mathrm{para}}_\alpha\hat j^{\mathrm{para}}_\beta\rangle$, a sum over pairs of states <a class="eqref" href="#b5-eq-Pipara"></a>; it carries the restoring force of the electrons, the band gaps";
   - electronic contact: "$-\frac{e^2}{m_{\mathrm e}c}n_{\mathrm e}(\mathbf r)\,\delta_{\alpha\beta}\,\delta(\mathbf r-\mathbf r')$ <a class="eqref" href="#b5-eq-Pidef"></a>: the inertia of free electrons";
   - ionic bubble: "$-\frac{\hbar}{c\Omega_{\mathrm V}}\sum_\nu\frac{\omega_{\mathbf q\nu}}{2}F^0_\nu\,j_{\alpha\nu}j^*_{\beta\nu}$ <a class="eqref" href="#b5-eq-Piionpara"></a>: exact for harmonic phonons";
   - ionic contact: "$-\frac{e^2}{c\Omega_{\mathrm V}}\sum_{\kappa\gamma}Z^*_{\kappa,\alpha\gamma}Z^*_{\kappa,\beta\gamma}/M_\kappa$ on the diagonal <a class="eqref" href="#b5-eq-Piiondia"></a>: the inertia of ions with their springs cut".
2. Check box "electrons first" (default off): shows row (d).

**Initial state.** Highlight "ionic bubble" (the piece that later sections use most), row (d) hidden.

**Caption key points.**
- *What to notice.* Every term of the series is the free photon line interrupted by the same self-energy, and the self-energy itself has only four pieces: each kind of charge contributes a bubble with two current vertices and a contact term with two photon lines at one point.
- *Marks.* Wavy single: free line $D^0$ (here transverse, $D^0=D^{0\mathrm T}$ <a class="eqref" href="#b5-eq-D0T"></a>); wavy double: dressed line $D$; zigzag: free phonon line $F^0_\nu$ <a class="eqref" href="#b5-eq-F0"></a>; solid with arrows: electrons; squares: current vertices ($-\frac1c\hat\mathbf j\cdot\hat\mathbf A$); dot and ring: contact vertices of electrons and ions.
- *Try.* Select each piece and read its formula. Turn on "electrons first": the electronic insertions collapse into $D^{\mathrm e}$, and only the ionic pieces remain between dressed lines; this is the starting point of Section B.7.
- *Model.* The Dyson equation <a class="eqref" href="#b5-eq-dyson"></a>, the four pieces <a class="eqref" href="#b5-eq-Pidef"></a>, the regrouping <a class="eqref" href="#b5-eq-dysonDe"></a>; the Coulomb term is set aside, as in the section. The ionic pair is exact within the harmonic approximation; the electronic bubble neglects the electron–hole attraction (no excitons). After the EPP manuscript's diagrams (its Fig. 7 and Fig. 11 draw the full longitudinal-plus-transverse versions; source reference in the caption only as "the EPP manuscript").

---

## B5.3 · `fig-electrons.html` · `viz-b5-electrons` · port of `viz-b-photon-se-electrons`

**Title.** "Electrons: truncated band sums and the fake Drude metal" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. With a complete set of transitions, "contact term plus bubble" and the subtracted form (B5.19) give the same $\varepsilon(\omega)$; with a truncated band sum they do not.
2. The naive form "contact term + truncated bubble" leaves $-(1-S)\,\omega_{\mathrm p}^2/\omega^2$ in $\varepsilon$: the insulator looks like a metal ($\mathrm{Re}\,\varepsilon\to-\infty$ as $\omega\to0$), with the missing f-sum weight $1-S$ sitting in a Drude peak at zero frequency.
3. The subtracted form misses the same weight at high frequency instead, where it matters little, and its $\varepsilon(0)$ converges quickly with the number of transitions, because high transitions contribute $W_s/\omega_s^2$.

**Include position.** Section `b5-truncation` ("Truncated band sums"), right after the caveat callout "What the subtracted form buys in practice" and the sentence "Figure B5.3 shows the two forms side by side for a model insulator with the plasma frequency of silicon." The section on the ions (`b5-ions`) follows.

**Physics (keep v1's model; v2 notation).**
- Twenty transitions $s=1,\dots,20$ at $\hbar\omega_s=\hbar\omega_1+(s-1)\times1.2$ eV, with weights $w_s\propto\mathrm e^{-(s-1)/3}$ normalized to $\sum_{s=1}^{20}w_s=1$, and strengths $W_s=w_s\,\omega_{\mathrm p}^2$, $\hbar\omega_{\mathrm p}=16.6$ eV (valence electrons of Si, illustrative). All twenty exhaust the f-sum.
- Subtracted (gauge-invariant) form, (B5.20) with damping: $\varepsilon(\omega)=1+\sum_{s\le N_{\mathrm t}}W_s/(\omega_s^2-\omega^2-\mathrm i\gamma\omega)$.
- Naive form: the same plus $-(1-S)\,\omega_{\mathrm p}^2/\omega^2$, with $S=\sum_{s\le N_{\mathrm t}}w_s$ (the fraction of the f-sum captured).
- Running f-sum $\frac{2}{\pi\omega_{\mathrm p}^2}\int_0^\omega\omega'\,\mathrm{Im}\,\varepsilon(\omega')\,\mathrm d\omega'$: it rises to $S$ for the subtracted form (the remaining $1-S$ lies in transitions not kept); for the naive form add the zero-frequency Drude weight $1-S$ at $\omega\to0^+$ (v1 does this; keep it, with the break at 0).
- Check values with $\hbar\omega_1=3.6$ eV: $N_{\mathrm t}=3$: $S=0.633$, $\varepsilon(0)=10.58$ (subtracted), fake plasma energy $\hbar\omega_{\mathrm p}\sqrt{1-S}=10.06$ eV; $N_{\mathrm t}=1$: $S=0.284$, $\varepsilon(0)=7.03$; $N_{\mathrm t}=10$: $S=0.966$, $\varepsilon(0)=11.78$; all 20: $\varepsilon(0)=11.81$ (Si measures 11.7).
- v1's code computes exactly these (function `eps(x, N, d1, g)` returns $1+\sum W/(D^2-x^2-\mathrm igx)$; the naive curve adds `-(1-S)*WP2/(x*x)`); keep the numerics, rename the variables.

**What to draw** (as v1). (a) $\mathrm{Re}\,\varepsilon$ (solid, `electron`) and $\mathrm{Im}\,\varepsilon$ (filled, `electron`, alpha 0.14) of the subtracted form against $\hbar\omega$ from 0 to 20 eV, fixed ordinate $[-25,55]$; the naive $\mathrm{Re}\,\varepsilon$ dashed in `c5`; a horizontal line at $\varepsilon(0)$ with all 20 transitions; an arrow and the label "fake Drude $-(1-S)\,\omega_{\mathrm p}^2/\omega^2$" pointing to where the naive branch enters the frame. (b) The running f-sum from 0 to 32 eV, ordinate $[0,1.15]$: a line at 1 labelled "f-sum: $\omega_{\mathrm p}^2$", a line at $S$, the subtracted curve (solid) and the naive curve (dashed) with the label "zero-frequency weight $1-S$". Axis labels `ħω (eV)`, `ε(ω)`, `running f-sum / ω_p^2`.

**Controls (4).** Number of transitions kept $N_{\mathrm t}$: 1–20, step 1, default 3. First transition $\hbar\omega_1$: 2–5 eV, default 3.6 eV. Broadening $\hbar\gamma$: 0.1–1 eV, default 0.3 eV. Check box "show contact term + truncated bubble" (default on).

**Readouts.** $N_{\mathrm t}$; $S$; $\varepsilon(0)$ of the subtracted form; $\varepsilon(0)$ with all 20 transitions; fake plasma energy $\hbar\omega_{\mathrm p}\sqrt{1-S}$ (eV).

**Initial state.** $N_{\mathrm t}=3$, $\hbar\omega_1=3.6$ eV, $\hbar\gamma=0.3$ eV, naive curve on: the dashed curve dives to $-\infty$ below about 3 eV while the solid one stays at $\varepsilon(0)=10.6$.

**Caption key points.**
- *What to notice.* With three transitions the naive curve turns the insulator into a metal, yet it satisfies the f-sum: the contact term insists on the full $\omega_{\mathrm p}^2$ and puts what the bubble misses at zero frequency. The subtracted form misses the same weight at high energy instead, and its $\varepsilon(0)$ is already within 10 % of the converged value.
- *Marks.* Solid and filled: real and imaginary parts of the subtracted form <a class="eqref" href="#b5-eq-Pie"></a>; dashed: contact term plus truncated bubble; horizontal line: $\varepsilon(0)$ with all twenty transitions; in (b) the lines at 1 and $S$.
- *Try.* Raise $N_{\mathrm t}$ to 20: the two curves merge and $S=1$. Set $N_{\mathrm t}=1$: the fake plasma energy rises to 14 eV. Move $\hbar\omega_1$: $\varepsilon(0)$ follows $1+\sum_sW_s/\omega_s^2$.
- *Model.* Twenty oscillators of the form (B5.20) <a class="eqref" href="#b5-eq-lorentz"></a>, spacing 1.2 eV, weights $\propto\mathrm e^{-(s-1)/3}$ that exhaust the f-sum with $\hbar\omega_{\mathrm p}=16.6$ eV (silicon valence electrons; a model, not a band structure); damping as in a driven oscillator.

**What in v1 must change.**
- Index $i$ and the symbols $\Delta_i$, $f_i$ (NOTATION: $i$ is never an index; $\Delta$ is reserved; $f$ is the Fermi–Dirac occupation): use $s$, $\hbar\omega_s$, $w_s$, $W_s$; the first gap control is "first transition $\hbar\omega_1$".
- Legends "Im ε (E.281)", "Re ε, gauge invariant (E.281)" and similar: drop source equation numbers; use "subtracted form" and "contact term + truncated bubble". "TRK: ω_p^2" may stay as "f-sum: ω_p^2".
- Caption: v1's single paragraph becomes the four-part caption above; mention the Thomas–Reiche–Kuhn name once.
- Ids and classes: `viz-b5-electrons`, classes `b5-electrons-…`; read colours at draw time.

---

## B5.4 · `fig-sumrule.html` · `viz-b5-sumrule` · port of `viz-b-photon-se-sumrule`

**Title.** "Ionic sum rule: springs move infrared strength, never the total" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. The strengths $|j_\nu(q)|^2$ of the three modes change with the wavevector, the springs and the masses, but their sum is always $e^2\sum_\kappa Z^{*2}_\kappa/M_\kappa$: the completeness identity (B5.25) read as a sum rule.
2. The total is set by charges and masses only (it is the contact term (B5.24)), which is why the contact term can be computed without knowing any force constant.
3. At $q=0$ the acoustic mode carries no strength (Born-charge sum rule); at small $q$ it carries only a small internal-strain dipole. The basis phases change how the strength is shared at finite $q$, never the total.

**Include position.** Section `b5-sumrule` ("Completeness: the contact term is minus the static bubble"), right after the sum-rule display and the sentence "The springs move infrared strength between the modes, differently at each $\bq$, but they cannot change the total. Figure B5.4 shows this for a chain with three atoms per cell." The total ionic self-energy (B5.26) follows.

**Physics (v2 notation; keep v1's model).**
- A one-dimensional chain with three atoms per cell, motion along the chain only (three branches $\nu=1,2,3$). Masses $M_1=30$ amu, $M_2$ (control), $M_3=60$ amu. Springs $K_{12}=6$ eV/Å² and $K_{23}=4$ eV/Å² inside the cell (fixed), $K_{31}$ (control) between atom 3 of cell $\ell$ and atom 1 of cell $\ell+1$. Positions in the cell $\tau_\kappa=0,\ a/3,\ 2a/3$ (model).
- Force constants transformed with the phase on the partner atom's cell, (A1.18) <a class="eqref" href="#a1-eq-dyn"></a>: $\Phi_{11}(q)=K_{12}+K_{31}$, $\Phi_{22}=K_{12}+K_{23}$, $\Phi_{33}=K_{23}+K_{31}$, $\Phi_{12}=\Phi_{21}=-K_{12}$, $\Phi_{23}=\Phi_{32}=-K_{23}$, $\Phi_{13}=-K_{31}\,\mathrm e^{-\mathrm iqa}$, $\Phi_{31}=-K_{31}\,\mathrm e^{+\mathrm iqa}$. Dynamical matrix $\tilde\Phi_{\kappa\kappa'}=\Phi_{\kappa\kappa'}/\sqrt{M_\kappa M_{\kappa'}}$; diagonalize with `EPP.eigHerm` (replaces v1's 6×6 real embedding with `EPP.jacobi`); $\hbar\omega=64.654\ \text{meV}\times\sqrt{\text{eigenvalue}}$.
- Born charges $Z^*_1$, $Z^*_2$ (controls) and $Z^*_3=-Z^*_1-Z^*_2$, so that the Born-charge sum rule <a class="eqref" href="#a4-eq-asr"></a> holds.
- Ionic vertex (B5.22) at $G=0$ in one dimension: $j_\nu(q)=\sum_\kappa eZ^*_\kappa\,e_{\kappa,\nu}(q)\,\mathrm e^{-\mathrm iq\tau_\kappa}/\sqrt{M_\kappa}$, with normalized eigenvectors (the overall phase of each eigenvector drops out of $|j_\nu|^2$). Plot $|j_\nu|^2$ in units of $e^2/\text{amu}$. Target: $e^2\sum_\kappa Z^{*2}_\kappa/M_\kappa$, which equals $-c\,\Omega\,\Pi^{\mathrm{ion,dia}}_{00}$ per cell (B5.24).
- Check values at the defaults ($M_2=16$ amu, $K_{31}=3$ eV/Å², $Z^*_1=2$, $Z^*_2=-1.2$, so $Z^*_3=-0.8$; target $0.2340$): at $qa/\pi=0.35$: $\hbar\omega_\nu=7.9,\ 31.0,\ 57.7$ meV and $|j_\nu|^2=0.0005,\ 0.0617,\ 0.1718$ (sum 0.2340); at $q=0$: $0,\ 32.4,\ 57.5$ meV and $0,\ 0.0558,\ 0.1782$; at $qa/\pi=1$: $18.9,\ 24.5,\ 58.3$ meV and $0.0194,\ 0.0745,\ 0.1401$. The relative difference between sum and target is $\sim10^{-16}$.

**What to draw** (as v1). (a) The three branches $\hbar\omega_\nu(q)$ for $qa/\pi\in[0,1]$ (colours `muted`, `phonon`, `c5`; the third dashed), with dots at every sixth sample whose radius grows with $|j_\nu|$ (the infrared strength along each branch), and a vertical line at the chosen $q$; clicking or dragging in (a) sets $q$. Axis labels `q a / π`, `ħω_ν (meV)`. (b) A bar chart at the chosen $q$: one bar per mode ($\nu=1,2,3$, each labelled with its energy), and a stacked bar of the three labelled $\sum_\nu$; a dashed horizontal line at the target labelled `e^2 Σ_κ Z_κ^{*2}/M_κ (contact term)`. Axis label `|j_ν|^2 (e^2/amu)`. Fixed ordinate per target (v1's `niceUp` logic is fine).

**Controls (5).** Wavevector $qa/\pi$: 0–1, step 0.01, default 0.35. Inter-cell spring $K_{31}$: 0.5–15 eV/Å², default 3. Mass $M_2$: 4–100 amu, default 16. Born charge $Z^*_1$: $-3$ to 3, default 2. Born charge $Z^*_2$: $-3$ to 3, default $-1.2$.

**Readouts.** $\hbar\omega_\nu$ for $\nu=1,2,3$ (meV); $\sum_\nu|j_\nu|^2$; $e^2\sum_\kappa Z^{*2}_\kappa/M_\kappa$; $Z^*_3$; relative difference.

**Initial state.** The defaults above: at $qa/\pi=0.35$ the acoustic bar is nearly empty, the two optical bars share the strength unequally, and the stacked bar ends exactly on the dashed line.

**Caption key points.**
- *What to notice.* Moving $q$ or $K_{31}$ reshuffles strength between the three modes, but the stacked bar always ends on the dashed line: that is the completeness identity <a class="eqref" href="#b5-eq-sum"></a>. At $q=0$ the acoustic mode carries nothing.
- *Marks.* Branch colours, dot radius $\propto|j_\nu|$, the vertical line at the chosen $q$, bars and stacked bar, dashed line = contact term <a class="eqref" href="#b5-eq-Piiondia"></a>.
- *Try.* Drag $q$ from 0 to $\pi/a$: the acoustic strength grows from zero, the total stays. Change $K_{31}$: the branches move, the total stays. Change $M_2$ or a Born charge: the total itself changes, because the contact term contains masses and charges.
- *Model.* Chain of three atoms per cell with nearest-neighbour springs, dynamical matrix with the phase on the partner atom's cell <a class="eqref" href="#a1-eq-dyn"></a>, vertex <a class="eqref" href="#b5-eq-jion"></a> with the basis phases $\mathrm e^{-\mathrm iq\tau_\kappa}$, Born charges that add up to zero. All parameters are model values.

**What in v1 must change.**
- $\sigma$ as the branch index everywhere (axis label `|j_σ|^2`, ticks `σ=1`, `Σ_σ`, readout `ħω_σ`, legend `σ = 1 (acoustic)`): use $\nu$. The linter rejects `σ` in script strings.
- Masses `m_1, m_2, m_3` → $M_1, M_2, M_3$; springs `k_{12}, k_{23}, k_{31}` → $K_{12}, K_{23}, K_{31}$; Born charges `Z_1, Z_2` → $Z^*_1, Z^*_2$ (dimensionless; the dipole is $eZ^*$; units of $|j|^2$ are $e^2$/amu); readout `Σ_a Z_a^2/m_a` → $e^2\sum_\kappa Z^{*2}_\kappa/M_\kappa$.
- v1's caption phrase "one phase $e^{-iqa}$ on the inter-cell bond, Giustino's sign" → "phase on the partner atom's cell, (A1.18)"; the source equation numbers (E.295′), (E.304′), (E.309′) → eqrefs `b5-eq-jion`, `b5-eq-Piiondia`, `b5-eq-sum`.
- Diagonalize with `EPP.eigHerm` and fix phases with `EPP.fixPhase` (the physics is unchanged; v1's matrix already has the v2 phase convention, $\Phi_{13}=-K_{31}\mathrm e^{-\mathrm iqa}$).
- Caption into the four parts; ids and classes `viz-b5-sumrule`, `b5-sumrule-…`.

---

## B5.5 · `fig-gauge.html` · `viz-b5-gauge` · port of `viz-b-photon-se-gauge`

**Title.** "Contact term plus bubble: the gauge-invariant ionic self-energy and the infrared response" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. The ionic bubble at $\omega=0$ is exactly $+1$ and the contact term $-1$ (in units of $(\omega^{\mathrm{ion}}_{\mathrm p})^2/4\pi c$), so the total self-energy starts at zero: gauge invariance (B5.15), proved for the ions by completeness (B5.25).
2. Far above $\omega_{\mathrm{TO}}$ the total tends to the contact term: the lattice responds as free ions.
3. Through (B5.14) the total gives the infrared dielectric function (B5.27), with $\varepsilon(0)=\varepsilon_{\mathrm s}$, a Reststrahlen band between $\omega_{\mathrm{TO}}$ and $\omega_{\mathrm{LO}}$, and a predicted $\omega_{\mathrm{LO}}$ within 1 % of the measured value for four crystals. Dropping the contact term adds $+(\omega^{\mathrm{ion}}_{\mathrm p})^2/\omega^2$, a Drude divergence of the wrong sign.
4. The optical sum rule $\frac{2}{\pi}\int_0^\infty\omega\,\mathrm{Im}\,\varepsilon\,\mathrm d\omega=(\omega^{\mathrm{ion}}_{\mathrm p})^2$ holds whatever the damping.

**Include position.** Section `b5-ir` ("The infrared dielectric tensor"), right after the paragraph that ends "… it is the sum rule (B5.25) in optical form. Figure B5.5 assembles the pieces for four crystals." The closing callout "Where these objects go next" follows the figure.

**Physics (v2 notation; keep v1's model).**
- One polar mode at $\mathbf q\to0$ (cubic two-atom crystal; for hBN the in-plane mode of the uniaxial crystal) carries the whole sum rule: $(\hbar\omega^{\mathrm{ion}}_{\mathrm p})^2=0.75641\ \text{eV}^2\times Z^{*2}/(\Omega_{\mathrm V}[\text{Å}^3]\,M_{\mathrm r}[\text{amu}])$, i.e. $4\pi e^2Z^{*2}/(\Omega_{\mathrm V}M_{\mathrm r})$ (B5.24 for this crystal).
- Total (B5.26) with damping: $4\pi c\,\Pi^{\mathrm{ion}}=(\omega^{\mathrm{ion}}_{\mathrm p})^2\,\omega^2/(\omega^2_{\mathrm{TO}}-\omega^2-\mathrm i\gamma\omega)$. Contact (B5.24): $4\pi c\,\Pi^{\mathrm{ion,dia}}=-(\omega^{\mathrm{ion}}_{\mathrm p})^2$. Bubble $=$ total $-$ contact $=(\omega^{\mathrm{ion}}_{\mathrm p})^2\,\omega^2_{\mathrm{TO}}/(\omega^2_{\mathrm{TO}}-\omega^2-\mathrm i\gamma\omega)$, (B5.23); its value at $\omega=0$ is $+(\omega^{\mathrm{ion}}_{\mathrm p})^2$.
- Dielectric function (B5.27), (B3.6): $\varepsilon=\varepsilon_\infty+4\pi c\,\Pi^{\mathrm{ion}}/\omega^2=\varepsilon_\infty+(\omega^{\mathrm{ion}}_{\mathrm p})^2/(\omega^2_{\mathrm{TO}}-\omega^2-\mathrm i\gamma\omega)$. Without the contact term: $\varepsilon+(\omega^{\mathrm{ion}}_{\mathrm p})^2/\omega^2$. LO frequency from the zero: $\omega^2_{\mathrm{LO}}=\omega^2_{\mathrm{TO}}+(\omega^{\mathrm{ion}}_{\mathrm p})^2/\varepsilon_\infty$; $\varepsilon_{\mathrm s}=\varepsilon(0)=\varepsilon_\infty+(\omega^{\mathrm{ion}}_{\mathrm p}/\omega_{\mathrm{TO}})^2$.
- Presets (literature values, illustrative) — $Z^*$, $M_{\mathrm r}$ (amu), $\Omega_{\mathrm V}$ (Å³), $\varepsilon_\infty$, $\hbar\omega_{\mathrm{TO}}$ (meV), measured $\hbar\omega_{\mathrm{LO}}$ (meV) → predicted $\hbar\omega^{\mathrm{ion}}_{\mathrm p}$, $\hbar\omega_{\mathrm{LO}}$, $\varepsilon_{\mathrm s}$:
  - GaAs: 2.07, 36.1, 45.2, 10.9, 33.3, 36.2 → 44.6 meV, 35.9 meV, 12.69;
  - NaCl: 1.12, 13.9, 44.9, 2.34, 20.3, 32.7 → 39.0 meV, 32.6 meV, 6.03;
  - MgO: 1.98, 9.65, 18.7, 2.96, 49.1, 89.0 → 128.2 meV, 89.2 meV, 9.78;
  - hBN (in-plane): 2.71, 6.10, 18.1, 4.87, 169.9, 199.6 → 224.3 meV, 198.0 meV, 6.61.
  All predicted LO frequencies lie within 1 % of the measured ones.
- f-sum readout: integrate $\frac{2}{\pi(\omega^{\mathrm{ion}}_{\mathrm p})^2}\int_0^\infty\omega\,\mathrm{Im}\,\varepsilon\,\mathrm d\omega$ numerically (v1's substitution $x=\omega_{\mathrm{TO}}+\frac{\gamma}{2}\tan\theta$ works); it returns 1.000.

**What to draw** (as v1). (a) The three pieces in units of $(\omega^{\mathrm{ion}}_{\mathrm p})^2/4\pi c$ against $\hbar\omega$: contact term $-1$ (dashed, `c5`), bubble (real part, `phonon`), total (real part, thick `ink`), imaginary part of the total (dotted `ink`); a vertical line at $\omega_{\mathrm{TO}}$; three points at $\omega=0$ labelled "bubble(0) = +1", "total(0) = 0", "contact = −1". Ordinate fixed at $[-3,3]$. (b) $\mathrm{Re}\,\varepsilon$ (`ink`), $\mathrm{Im}\,\varepsilon$ (filled `phonon`), $\mathrm{Re}\,\varepsilon$ without the contact term (dashed `photon`), the band $\omega_{\mathrm{TO}}<\omega<\omega_{\mathrm{LO}}$ shaded with the label "Re ε < 0", a horizontal line labelled `ε_∞`, a dot at $\varepsilon(0)=\varepsilon_{\mathrm s}$. Abscissa per preset from 0 to about $2.2\,\omega_{\mathrm{LO}}$ (v1's fixed-per-preset limits; they grow only when $Z^*$ pushes the curves off scale). Axis labels `ħω (meV)`, `4πcΠ/(ω_p^{ion})^2`, `ε(ω)`.

**Controls (4).** Material (select: GaAs, NaCl, MgO, hBN in-plane; resets $Z^*$ to the preset). Born charge $Z^*$: 0.2–3.5, step 0.01 (default the preset value, 2.07). Damping $\hbar\gamma$: 0.1–10 meV, logarithmic, default 1 meV. Check box "show ε without the contact term" (default on).

**Readouts.** $\hbar\omega^{\mathrm{ion}}_{\mathrm p}$; $\hbar\omega_{\mathrm{TO}}$; predicted $\hbar\omega_{\mathrm{LO}}$ with the measured value in brackets; $\varepsilon_{\mathrm s}=\varepsilon(0)$; f-sum $\frac{2}{\pi(\omega^{\mathrm{ion}}_{\mathrm p})^2}\int\omega\,\mathrm{Im}\,\varepsilon\,\mathrm d\omega$.

**Initial state.** GaAs, $Z^*=2.07$, $\hbar\gamma=1$ meV, contact-free curve on: $\hbar\omega^{\mathrm{ion}}_{\mathrm p}=44.6$ meV, predicted LO 35.9 meV (measured 36.2), $\varepsilon_{\mathrm s}=12.7$, the dashed curve diverging upward at low frequency.

**Caption key points.**
- *What to notice.* The bubble at $\omega=0$ is exactly $+1$ and the contact term $-1$, so the total starts at zero; far above $\omega_{\mathrm{TO}}$ it tends to the contact term. Dropping the contact term adds $+(\omega^{\mathrm{ion}}_{\mathrm p})^2/\omega^2$, a Drude divergence with the wrong sign. The predicted LO frequency lands within 1 % of experiment for all four crystals.
- *Marks.* Dashed: contact term <a class="eqref" href="#b5-eq-Piiondia"></a>; thin: bubble <a class="eqref" href="#b5-eq-Piionpara"></a>; thick: total <a class="eqref" href="#b5-eq-Piion"></a>; dotted: its imaginary part; shaded band: the Reststrahlen band; dot: $\varepsilon_{\mathrm s}$.
- *Try.* Switch materials: the band moves from 20 meV (NaCl) to 170 meV (hBN). Raise $Z^*$: the band widens as $Z^{*2}$ and $\varepsilon_{\mathrm s}$ grows. Change the damping by a factor 100: the peaks broaden, the f-sum readout stays at 1.
- *Model.* One polar mode carrying the whole sum rule, (B5.26) and (B5.27) <a class="eqref" href="#b5-eq-epsir"></a> with a damping rate in place of $\eta$; electrons as the background $\varepsilon_\infty$; Born charges, masses, volumes and frequencies from the literature (illustrative); only the $\mathbf G=\mathbf G'=0$ element, with no Lorentz local field. The same function as Section B.3 <a class="eqref" href="#b3-eq-eps"></a>.

**What in v1 must change.**
- $\Omega_p$ for the ionic plasma frequency → $\omega^{\mathrm{ion}}_{\mathrm p}$ (canvas `ω_p^{ion}`); `ε^∞` → `ε_∞`; the reduced mass `μ` → $M_{\mathrm r}$; `Ω_uc` → $\Omega_{\mathrm V}$.
- Legends "dia (E.304′)", "para (E.297′)", "para + dia (E.311)" → "contact", "bubble", "total" (eqrefs only in the caption).
- Axis label `4πcΠ / Ω_p^2` → `4πcΠ/(ω_p^{ion})^2`; point labels "para(0) = +1", "dia = −1" → "bubble(0) = +1", "contact = −1".
- Caption into the four parts; ids and classes `viz-b5-gauge`, `b5-gauge-…`.
