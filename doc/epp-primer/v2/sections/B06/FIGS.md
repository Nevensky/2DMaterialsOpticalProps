# B.6 figures: specifications for the figure authors

Section file: `v2/sections/B06-propagator.html` (id `b-propagator`, prefix `b6`). Three figures, in reading order:

| # | File | Id | Kind |
|---|---|---|---|
| B6.1 | `B06/fig-sectors.html` | `viz-b6-sectors` | new (may reuse the diagram code of the v1 "Diagrams" tab) |
| B6.2 | `B06/fig-polariton.html` | `viz-b6-polariton` | new |
| B6.3 | `B06/fig-map.html` | `viz-b6-map` | port of v1 `viz-b-propagator-map` (`doc/epp-primer/sections/B06-propagator.html`) |

Conventions for all three (from the section, binding):

- Gaussian units. $\hbar c=1973.27$ eV Å. Wavevectors in Å⁻¹ (log axes) or μm⁻¹ (linear axes; 1 μm⁻¹ = 10⁻⁴ Å⁻¹). Energies $\hbar\omega$ in meV (phonons) or eV (electrons).
- Retarded functions: $\omega\to\omega+\ii\eta$, or $+\ii\gamma\omega$ for a damped oscillator. Time dependence $\ee^{-\ii\omega t}$.
- Symbols (v2): $D$ total propagator; $D^0=D^{0\mathrm L}+D^{0\mathrm T}$; $D^{0\mathrm L}=-\frac{4\pi c}{\omega^2}\hat q_\alpha\hat q_\beta$ (B6.10, `b6-eq-D0L`); $D^{0\mathrm T}=-\frac{4\pi c(\delta_{\alpha\beta}-\hat q_\alpha\hat q_\beta)}{(\omega+\ii\eta)^2-c^2q^2}$ (B6.2, `b6-eq-lines`); Coulomb line between charges $v/e^2=4\pi/q^2$; four sectors $D_{\rho\rho}$, $D_{\rho\beta}$, $D_{\alpha\rho}$, $D^{\mathrm T}_{\alpha\beta}$ (B6.3, `b6-eq-four`); correlators $\Pi_{\rho\rho}$, $\Pi_{\rho\beta}$, $\Pi_{\alpha\rho}$, $\Pi_{\alpha\beta}$ (B6.4, `b6-eq-induced`); $\Pi=\Pi^{\mathrm e}+\Pi^{\mathrm{ion}}$; L and T parts $X^{\mathrm L}$, $X^{\mathrm{LT}}$, $X^{\mathrm{TL}}$, $X^{\mathrm T}$ (B6.8, `b6-eq-LT`); $\varepsilon^{\mathrm L}$, $\varepsilon^{\mathrm T}$ (B6.17, `b6-eq-DLT`); polariton branches $\omega_{\mathrm{UP}}$, $\omega_{\mathrm{LP}}$ (B6.20, `b6-eq-polariton`); $\varepsilon_\infty$, $\varepsilon_{\mathrm s}$, $\omega_{\mathrm{TO}}$, $\omega_{\mathrm{LO}}$, $\omega^{\mathrm{ion}}_{\mathrm p}$ (ionic plasma frequency), $\omega_{\mathrm p}$ (electronic), $\OV$ (cell volume), branch index $\nu$.
- Forbidden in labels, readouts, aria-labels and strings: $\Pi_{jj}$, $\Pi_{\rho j}$, $\Pi_{j\rho}$, $D_{11}$, $D_{1\nu}$, $D_{\mu1}$, $D^{0l}$, $D^{0t}$, $\hat D^{0t}$, $\Omega_{\rm uc}$, $\Omega_p$, $\epsilon$ (use $\varepsilon$), $\epsilon_0$ for the static constant (use $\varepsilon_{\mathrm s}$), $\mu$ as a Cartesian index, $\sigma$ for a branch (σ is reserved for the conductivity), "E.350"-style source equation numbers. Refer to equations with `<a class="eqref" href="#b6-eq-…"></a>`.
- Colours: the whole electromagnetic field in the photon colour (`'photon'`); phonons `'phonon'`; electrons `'electron'`; totals and neutral marks `'muted'`/`'ink'`. Read colours at draw time with `EPP.color`.
- Material numbers (literature values, illustrative; the section's footnote states them): GaAs $\hbar\omega_{\mathrm{TO}}=33.3$ meV, $\hbar\omega_{\mathrm{LO}}=36.2$ meV, $\varepsilon_\infty=10.9$, so $\varepsilon_{\mathrm s}=\varepsilon_\infty\omega_{\mathrm{LO}}^2/\omega_{\mathrm{TO}}^2=12.88$ and $\hbar\omega^{\mathrm{ion}}_{\mathrm p}=\sqrt{\varepsilon_\infty(\omega_{\mathrm{LO}}^2-\omega_{\mathrm{TO}}^2)}=46.9$ meV; lattice constant $a=5.653$ Å. NaCl 20.3 / 32.7 meV, $\varepsilon_\infty=2.34$. SiC 98.7 / 120.5 meV, $\varepsilon_\infty=6.52$. hBN in-plane 169.9 / 199.6 meV, $\varepsilon_\infty=4.87$ (treated as isotropic; the section's caveat explains why that is only a model). Sodium plasmon $\hbar\omega_{\mathrm p}=5.7$ eV.

---

## B6.1 · `fig-sectors.html` · `viz-b6-sectors` · new

**Title.** "From four sectors to one tensor" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. The $4\times4$ bookkeeping of Coulomb lines and photon lines has four slots, $\rho$ and the three polarizations $\lambda=1,2,3$ of the basis $\be_1,\be_2,\be_3=\hat\bq$. Its $\lambda=3$ (longitudinal) row and column are empty: the vector potential of the Coulomb-gauge bookkeeping is transverse. The scalar slot $\rho$ does the longitudinal job (Coulomb).
2. The conversion of (B6.9) moves the $\rho$ row and column into the empty $\lambda=3$ row and column, multiplied by fixed factors. The result is one $3\times3$ tensor whose blocks are $D^{\mathrm L}$ (entry 33), $D^{\mathrm{LT}}$ (31, 32), $D^{\mathrm{TL}}$ (13, 23) and $D^{\mathrm T}$ (block 1–2 × 1–2).
3. In diagrams: at lowest order the four sectors are a scalar, a row, a column and a $3\times3$ block, which cannot be added; after the rule "charge vertex → current vertex, $v/e^2\to D^{0\mathrm L}$" they are all $3\times3$ and add up to $D^0+D^0\Pi D^0$; the same works at every order and gives the geometric series, the Dyson equation (B6.14).

**Include position.** At the start of `<h3 id="b6-sectors">From four sectors to one tensor</h3>`, right after its first paragraph, which ends: "Rewriting both ends of every sector in vector form turns the four sectors into four $3\times3$ tensors, which can be added. Figure B6.1 draws the result: the scalar slot of the $4\times4$ matrix moves into the empty longitudinal slot of a $3\times3$ one." The figure is also cited earlier (subsection "Lines and bubbles together", for the dot and square vertices) and later (subsection "One Dyson equation", for step 3).

**Physics** (all in $(\bq,\omega)$ space, homogeneous medium):
- Sectors, (B6.3) `b6-eq-four`: $A_0=D_{\rho\rho}\rho^{\mathrm{ext}}+\sum_\beta D_{\rho\beta}j^{\mathrm{ext}}_\beta$, $A_\alpha=D_{\alpha\rho}\rho^{\mathrm{ext}}+\sum_\beta D^{\mathrm T}_{\alpha\beta}j^{\mathrm{ext}}_\beta$. Free: $D_{\rho\rho}=v/e^2=4\pi/q^2$, $D^{\mathrm T}=D^{0\mathrm T}$, mixed sectors zero.
- Lowest order, (B6.5) `b6-eq-lowest`: $D_{\rho\rho}=\frac{v}{e^2}+\frac{v}{e^2}\Pi_{\rho\rho}\frac{v}{e^2}$; $D_{\rho\beta}=\frac{v}{e^2}\sum_\gamma\Pi_{\rho\gamma}D^{0\mathrm T}_{\gamma\beta}$; $D_{\alpha\rho}=\sum_\gamma D^{0\mathrm T}_{\alpha\gamma}\Pi_{\gamma\rho}\frac{v}{e^2}$; $D^{\mathrm T}=D^{0\mathrm T}+D^{0\mathrm T}\Pi D^{0\mathrm T}$.
- Conversion, (B6.9) `b6-eq-sectors`: $D^{\mathrm L}=-\frac{cq^2}{\omega^2}D_{\rho\rho}$, $D^{\mathrm{LT}}_\beta=-\frac{cq}{\omega}D_{\rho\beta}$, $D^{\mathrm{TL}}_\alpha=\frac{q}{\omega}D_{\alpha\rho}$, $D^{\mathrm T}$ unchanged; equivalently the tensor factors $-\frac{c}{\omega^2}q_\alpha q_\beta$, $-\frac{c}{\omega}q_\alpha$, $\frac{1}{\omega}q_\beta$, $1$.
- Longitudinal line, (B6.10) `b6-eq-D0L`: $D^{0\mathrm L}=-\frac{c}{\omega^2}q_\alpha q_\beta\frac{v}{e^2}=-\frac{4\pi c}{\omega^2}\hat q_\alpha\hat q_\beta$.
- Lowest order of the tensor, (B6.11) `b6-eq-lowest3`; Dyson equation, (B6.14) `b6-eq-dysonD`.
- The figure computes nothing numerically; it is a diagram figure with a stepper.

**What to draw.**
- *Top panel, geometric (SVG):* two matrix grids, side by side on wide screens, stacked on phones.
  - Left, the $4\times4$ matrix with rows and columns labelled $\rho$, $1$, $2$, $3$ (the last three are $\lambda$; a small triad icon with bold upright labels $\be_1$, $\be_2$, $\be_3=\hat\bq$ explains them; use `drawRich` with `{bold:true}` if drawn on canvas, or `<tspan font-weight>` in SVG). Cell fills: $(\rho,\rho)$ Coulomb corner, labelled $D_{\rho\rho}$; $(\rho,1)$, $(\rho,2)$ the row $D_{\rho\beta}$; $(1,\rho)$, $(2,\rho)$ the column $D_{\alpha\rho}$; block $(1,1)$–$(2,2)$ the transverse sector $D^{\mathrm T}$; every cell with a 3 (row 3, column 3, and $(\rho,3)$, $(3,\rho)$) hatched grey and marked "empty: no free longitudinal photon". Use the photon colour at two strengths for Coulomb corner and transverse block, a mixed tint for row and column.
  - Right, the $3\times3$ tensor with rows and columns $1,2,3$: block $(1,1)$–$(2,2)$ = $D^{\mathrm T}$; $(3,1)$, $(3,2)$ = $D^{\mathrm{LT}}$; $(1,3)$, $(2,3)$ = $D^{\mathrm{TL}}$; $(3,3)$ = $D^{\mathrm L}$.
  - Arrows from the left $\rho$ row, $\rho$ column and corner to the right row 3, column 3 and entry 33, labelled with the factors $\times(-cq/\omega)$, $\times(q/\omega)$, $\times(-cq^2/\omega^2)$. Arrows labelled in the photon text colour.
- *Bottom panel, Feynman diagrams* (EPP.fd; the `chain()` helper of the v1 tab can be reused): one block per step.
  - Step 1 (B6.5): four rows. (a) $D_{\rho\rho}=$ dashed Coulomb line $+$ [dashed – $\Pi^{\mathrm{ion}}$ bubble with two dots – dashed] $+$ [same with $\Pi^{\mathrm e}$] $+\dots$; (b) $D_{\rho\beta}=$ [dashed – bubble (dot left, square right) – wavy] for ion and e; (c) $D_{\alpha\rho}=$ [wavy – bubble (square left, dot right) – dashed]; (d) $D^{\mathrm T}_{\alpha\beta}=$ wavy $+$ [wavy – bubble (two squares) – wavy]. At the right of each row a small tag: "scalar · length²", "row of 3 · length × time", "column of 3 · length²", "3 × 3 · length × time".
  - Step 2 (B6.9–B6.11): the same four rows with every vertex a square and every dashed Coulomb line now the dashed line $D^{0\mathrm L}$ between squares; at the right of each row a boxed factor $-\frac{c}{\omega^2}q_\alpha q_\beta$, $-\frac{c}{\omega}q_\alpha$, $\frac{1}{\omega}q_\beta$, $1$; a fifth row (e), separated by a rule: $D=$ [dashed over wavy $D^0$] $+$ [$D^0$ – $\Pi^{\mathrm{ion}}$ – $D^0$] $+$ [$D^0$ – $\Pi^{\mathrm e}$ – $D^0$] $+\dots$
  - Step 3 (B6.14): (a) three second-order diagrams of different sectors, separated by commas: dashed–$\Pi^{\mathrm e}$(dots)–dashed–$\Pi^{\mathrm{ion}}$(dots)–dashed; dashed–$\Pi^{\mathrm e}$(dot, square)–wavy–$\Pi^{\mathrm{ion}}$(square, dot)–dashed; dashed–$\Pi^{\mathrm e}$(dots)–dashed–$\Pi^{\mathrm{ion}}$(dot, square)–wavy; (b) the same after the rule, now added with "+"; (c) summed into $D^0$–$\Pi$–$D^0$–$\Pi$–$D^0$; (d) the Dyson equation drawn: heavy double line $D$ $=$ $D^0$ $+$ $D^0$–$\Pi$–$D$.
- *Glyphs:* dashed = Coulomb line $v/e^2$ (between dots) or $D^{0\mathrm L}$ (between squares); wavy = $D^{0\mathrm T}$; wavy over dashed = $D^0$; heavy wavy over double dashes = the dressed $D$; all in the photon colour. Dots = charge vertices, squares = current vertices (ink). Bubbles: electron colour $\Pi^{\mathrm e}$, phonon colour $\Pi^{\mathrm{ion}}$, grey $\Pi$. Labels on lines and bubbles in v2 notation (`D^{0L}` as rich text `D^{0L}` is acceptable on canvas; prefer `D^{0\mathrm L}` in HTML/MathJax).

**Controls** (3):
1. Step (seg): "1 · four sectors", "2 · convert", "3 · Dyson". Default 1. Changing the step updates both panels: in step 1 the right grid is faded and the arrows hidden; in step 2 the $\rho$ row and column slide into slot 3 (one calm animation via `EPP.anim`, skipped under reduced motion) and the arrows with factors appear; in step 3 the right grid is labelled $D=D^0+D^0\Pi D$.
2. Toggle "show kinds and units" (default on): the row tags of step 1 and, in steps 2–3, "all blocks: 3 × 3 · length × time".
3. Toggle "electron and ion bubbles separately" (default on); off merges them into the grey $\Pi$ in every step.

**Readouts** (one line, `aria-live="polite"`): step 1 "objects: scalar + row + column + 3 × 3 (cannot be added)"; step 2 "four 3 × 3 blocks: $D^{\mathrm L}$, $D^{\mathrm{LT}}$, $D^{\mathrm{TL}}$, $D^{\mathrm T}$ (they add up to $D$)"; step 3 "order $n$: $D^0(\Pi D^0)^n$, a geometric series".

**Initial state.** Step 1, kinds shown, bubbles separate: the reader sees four rows of different kinds and, above them, the $4\times4$ grid with its hatched empty slot 3.

**Caption key points** (200–260 words in total).
- *What to notice.* The $4\times4$ matrix has an empty longitudinal slot, and its scalar slot does the longitudinal job. The conversion moves one into the other; only then do the four sectors add, line by line and order by order, into one Dyson equation.
- *Marks.* Grid fills (Coulomb corner, mixed row and column, transverse block, hatched empty cells); arrows and factors; line glyphs; dots and squares; bubble colours.
- *Try.* Step 1 → 2: watch the $\rho$ row and column move into slot 3 and the dots turn into squares. Step 3: see three diagrams from different sectors become terms of one series. Untick "separately" to merge $\Pi^{\mathrm e}+\Pi^{\mathrm{ion}}$ into $\Pi$.
- *Model.* (B6.3), (B6.5), (B6.7) the connection formulas, (B6.9), (B6.10), (B6.11), (B6.14); homogeneous medium in $(\bq,\omega)$ space, random-phase approximation (lowest-order bubbles, no mixed electron–ion bubbles); after the diagrams of the EPP manuscript, redrawn in this primer's glyphs.

**If reusing v1 code** (v1 `sections/B06-propagator.html`, the `chain()` builder, the "Diagrams" tab `b06-v-dia` and the inline block figure `b-propagator-blockfig`): rename classes to `b6-sectors-…`; replace every source label: $D_{11}$, $D_{1\nu}$, $D_{\mu1}$, $D_{\mu\nu}$ → $D_{\rho\rho}$, $D_{\rho\beta}$, $D_{\alpha\rho}$, $D^{\mathrm T}_{\alpha\beta}$; $D^{ll}$, $D^{lt}$, $D^{tl}$, $D^{tt}$ → $D^{\mathrm L}$, $D^{\mathrm{LT}}$, $D^{\mathrm{TL}}$, $D^{\mathrm T}$; $\Pi_{jj}$ → $\Pi$; $\hat D^{0t}$, $D^{0t}$, $D^{0l}$ → $D^{0\mathrm T}$, $D^{0\mathrm L}$; $v$ → $v/e^2$ (the charges sit in the vertices); "E Eq. …" and "EPP Fig. 9" → eqrefs to (B6.5), (B6.9), (B6.14). Replace the real-space operators $-\frac{c}{\omega^2}\partial_\mu\partial'_\nu$, $i\frac c\omega\partial_\mu$, $\frac{i}{\omega}\partial'_\nu$ (μ is not allowed as an index) by the $\bq$-space factors above. The geometric grids are new.

---

## B6.2 · `fig-polariton.html` · `viz-b6-polariton` · new

**Title.** "Phonon polaritons and the LO mode in GaAs" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. In the transverse channel light and the TO phonon anticross into an upper and a lower polariton. Between $\omega_{\mathrm{TO}}$ and $\omega_{\mathrm{LO}}$ no transverse wave exists (Reststrahlen band). The upper branch starts at $\omega_{\mathrm{LO}}$ at $q=0$; the lower branch is $cq/\sqrt{\varepsilon_{\mathrm s}}$ at small $q$ and tends to $\omega_{\mathrm{TO}}$ at large $q$.
2. The longitudinal mode is flat at $\omega_{\mathrm{LO}}$ for all $q$ and does not depend on $c$.
3. Scale: all of this happens at $q\sim\omega_{\mathrm{TO}}/c\approx10^{-5}$–$10^{-4}$ Å⁻¹, about $5\times10^{-5}$ of the zone; on a log axis out to the zone the branches are the flat TO and LO lines of Part A for $q\gtrsim10^{-3}$ Å⁻¹.
4. Character: the phonon fraction changes continuously along each branch; at $q=\sqrt{\varepsilon_\infty}\,\omega_{\mathrm{TO}}/c$ both branches are exactly half phonon and half light.
5. Order of limits: raising $c$ squeezes the polariton region toward $q=0$; the LO line does not move.

**Include position.** In `<h4 id="b6-polar">The simplest model: a polar insulator</h4>`, after the paragraph that ends: "For GaAs, $Z^*=2.07$, $\OV=45.2$ Å³ and $M_{\mathrm r}=36.1$ amu give $\hbar\omega^{\mathrm{ion}}_{\mathrm p}=44.6$ meV, against $\hbar\omega_{\mathrm{TO}}=33.3$ meV and $\varepsilon_\infty=10.9$. Figure B6.2 shows the modes this crystal supports; the formalism that produces them follows the figure." It precedes `<h4 id="b6-modes">`. Later text cites it for the limits of (B6.20) and for the phonon fraction.

**Physics.**
- Dielectric function (B6.16) `b6-eq-polar`, written with the measured LO frequency: $\varepsilon(\omega)=\varepsilon_\infty\frac{\omega_{\mathrm{LO}}^2-\omega^2}{\omega_{\mathrm{TO}}^2-\omega^2}$, i.e. $(\omega^{\mathrm{ion}}_{\mathrm p})^2=\varepsilon_\infty(\omega_{\mathrm{LO}}^2-\omega_{\mathrm{TO}}^2)$; $\varepsilon_{\mathrm s}=\varepsilon_\infty\omega_{\mathrm{LO}}^2/\omega_{\mathrm{TO}}^2$.
- Transverse modes (B6.19) `b6-eq-modesT` → branches (B6.20) `b6-eq-polariton`:
  $\omega^2_{\mathrm{UP},\mathrm{LP}}=\frac{\omega^2_{\mathrm{LO}}}{2}+\frac{c^2q^2}{2\varepsilon_\infty}\pm\sqrt{\Big(\frac{\omega^2_{\mathrm{LO}}}{2}+\frac{c^2q^2}{2\varepsilon_\infty}\Big)^2-\frac{c^2q^2\omega^2_{\mathrm{TO}}}{\varepsilon_\infty}}$.
  Use `cq = λ·1973.27e3·q` (meV, q in Å⁻¹) with the speed-of-light factor λ. For tiny $q$ use the stable form $\omega^2_{\mathrm{LP}}=P/\omega^2_{\mathrm{UP}}$ with $P=c^2q^2\omega_{\mathrm{TO}}^2/\varepsilon_\infty$, to avoid cancellation.
- Longitudinal mode (B6.18) `b6-eq-modesL`: $\omega=\omega_{\mathrm{LO}}$ for every $q$.
- Light lines: $\omega=cq$, $cq/\sqrt{\varepsilon_\infty}$, $cq/\sqrt{\varepsilon_{\mathrm s}}$ (with the same λ).
- Phonon fraction (derivation box at the end of the section, step 5): $f_{\mathrm{ph}}=\dfrac{M}{M+\varepsilon_\infty+\varepsilon(\omega)}$ with $M=\dfrac{(\omega^{\mathrm{ion}}_{\mathrm p})^2(\omega^2+\omega^2_{\mathrm{TO}})}{(\omega^2_{\mathrm{TO}}-\omega^2)^2}$, evaluated on each transverse branch. On the LO branch ($\bB=0$): $f_{\mathrm{ph}}=M/(M+\varepsilon_\infty)=(\omega^2_{\mathrm{TO}}+\omega^2_{\mathrm{LO}})/(2\omega^2_{\mathrm{LO}})$, constant.
- Characteristic wavevector: $q_{\times}=\sqrt{\varepsilon_\infty}\,\omega_{\mathrm{TO}}/c$ (write it as this expression in labels; the section gives it no symbol).
- Check values (GaAs, λ = 1): $\varepsilon_{\mathrm s}=12.88$; $\sqrt{\varepsilon_\infty}\omega_{\mathrm{TO}}/c=5.571\times10^{-5}$ Å⁻¹ $=0.557$ μm⁻¹, wavelength 11.3 μm; there $\hbar\omega_{\mathrm{UP}}=41.15$ meV, $\hbar\omega_{\mathrm{LP}}=26.95$ meV, $f_{\mathrm{ph}}=0.500$ on both (exact for every preset). At $10\times$ that wavevector $\hbar\omega_{\mathrm{LP}}=33.270$ meV (0.03 meV below TO). At $q=10^{-5}$ Å⁻¹, $\hbar\omega_{\mathrm{LP}}=5.49$ meV, and 10.90 meV with λ = 2. $f_{\mathrm{ph}}$ at $q\to0$: 0.923 (UP), 0.077 (LP). Distance from $\Gamma$ to the zone edge along [100] (the X point of the fcc zone): $2\pi/a=1.111$ Å⁻¹; use it as the "zone edge" mark.
- Cost: closed form, a few hundred points per redraw; far below 15 ms.

**What to draw.**
- One plot, $\hbar\omega$ (meV) against $q$. Fixed $y$ range $[0,\ 1.6\,\hbar\omega_{\mathrm{LO}}]$ per preset (GaAs 0–58 meV).
- $x$ axis: linear mode $[0,\ 4\sqrt{\varepsilon_\infty}\omega_{\mathrm{TO}}/c]$ in μm⁻¹ (GaAs 0–2.23 μm⁻¹), fixed for each preset; when λ changes the axis stays and the curves move; log mode $[10^{-7},\ 1]$ Å⁻¹, fixed.
- Reststrahlen band: horizontal shaded strip between $\omega_{\mathrm{TO}}$ and $\omega_{\mathrm{LO}}$ (phonon colour, low alpha), labelled inside "Reststrahlen: no transverse wave".
- Dashed horizontals at $\omega_{\mathrm{TO}}$ and $\omega_{\mathrm{LO}}$ (phonon colour), labelled $\omega_{\mathrm{TO}}$, $\omega_{\mathrm{LO}}$ (rich text `ω_{TO}`, `ω_{LO}`), placed so they never overlap (one label below its line at the left, the other above its line at the right).
- Light lines (photon colour, thin): $\omega=cq$ dashed, labelled "$cq$"; $cq/\sqrt{\varepsilon_\infty}$ and $cq/\sqrt{\varepsilon_{\mathrm s}}$ dotted, labelled `cq/√ε_∞` and `cq/√ε_s`.
- Polariton branches UP and LP: thick (≈3 px), drawn as short segments whose colour mixes the photon colour ($f_{\mathrm{ph}}=0$) and the phonon colour ($f_{\mathrm{ph}}=1$) (mix `EPP.rgb('photon')` and `EPP.rgb('phonon')` at draw time). Labels `ω_{UP}` and `ω_{LP}` near the right end.
- LO branch: thick phonon-colour line at $\omega_{\mathrm{LO}}$ in a distinct style (for example long dashes, or a solid line with a small "L" marker every 80 px), labelled "LO, longitudinal ($\varepsilon^{\mathrm L}=0$)". It sits on top of the UP start point at $q=0$.
- Vertical marker at $\sqrt{\varepsilon_\infty}\omega_{\mathrm{TO}}/c$ (muted, dotted), labelled `√ε_∞ ω_{TO}/c`. In log mode also a marker at $2\pi/a$ labelled "zone edge (X)" and a bracket "polariton region" spanning roughly $0.1$–$10$ times the marker.
- Cursor: a vertical line at the chosen $q$ with dots on UP, LP and LO; draggable on the canvas (`touch-action: none` on that canvas only), also set by the slider below.
- Small colour key (HTML legend or a bar): "phonon fraction 0 (light) → 1 (phonon)".

**Controls** (4):
1. Crystal (select): GaAs (default), NaCl, SiC, hBN in-plane.
2. $q$ axis (seg): "linear (μm⁻¹)" / "logarithmic (Å⁻¹)". Default linear.
3. Speed of light $c\times\lambda$ (log range, $\lambda=1$ to $1000$, default 1); `aria-valuetext` "c times 10" etc. All $c$ in the formulas become $\lambda c$; the LO line and the Reststrahlen band stay put.
4. Cursor $q$ (range): linear mode in units of $\sqrt{\varepsilon_\infty}\omega_{\mathrm{TO}}/c$, 0–4, default 1; log mode $10^{-7}$–$1$ Å⁻¹, default at the marker. The canvas drag sets the same value. Arrow keys move it when the slider has focus.

**Readouts** (MathJax in HTML, `aria-live="polite"`): $\hbar\omega_{\mathrm{TO}}$, $\hbar\omega_{\mathrm{LO}}$, $\varepsilon_\infty$, $\varepsilon_{\mathrm s}$; $\sqrt{\varepsilon_\infty}\,\omega_{\mathrm{TO}}/c$ in Å⁻¹ and μm⁻¹ and its wavelength in μm; at the cursor: $q$, $\hbar\omega_{\mathrm{UP}}$ with $f_{\mathrm{ph}}$, $\hbar\omega_{\mathrm{LP}}$ with $f_{\mathrm{ph}}$, and the ratio of $q$ to $2\pi/a$ (lattice constants, illustrative: GaAs 5.653 Å, NaCl 5.64 Å, 3C-SiC 4.36 Å; for hBN use the in-plane $\Gamma$–M distance $2\pi/(\sqrt3\,a)$ with $a=2.50$ Å).

**Initial state.** GaAs, linear axis, λ = 1, cursor at $\sqrt{\varepsilon_\infty}\omega_{\mathrm{TO}}/c$: readouts show 41.1 and 27.0 meV, both 50 % phonon; the anticrossing, the dark band and the flat LO line are all visible.

**Caption key points.**
- *What to notice.* Light and the TO phonon anticross into two polaritons around a band where no transverse wave exists, while the LO mode stays flat at $\omega_{\mathrm{LO}}$; everything happens at $q$ of order $10^{-5}$–$10^{-4}$ Å⁻¹, a few hundred-thousandths of the Brillouin zone.
- *Marks.* Branch colour = phonon fraction; shaded band = Reststrahlen; photon-coloured light lines; phonon-coloured TO/LO lines; LO drawn in its own style because it is longitudinal; marker and cursor.
- *Try.* (1) Switch to the logarithmic axis: beyond about $10^{-3}$ Å⁻¹ only flat TO and LO remain, the electrostatic picture of Part A. (2) Raise $c\times\lambda$ to 1000: the polariton region slides toward $q=0$ and the LO line does not move (the order of limits of Section B.7). (3) Put the cursor on the marker: both branches are exactly half phonon. (4) Choose NaCl: a much wider Reststrahlen band.
- *Model.* (B6.16) with $\varepsilon_\infty$ for the electrons, (B6.18)–(B6.20); Lyddane–Sachs–Teller for $\varepsilon_{\mathrm s}$; energy split from the derivation box; no damping; local isotropic response; presets are literature values (illustrative), hBN treated as isotropic.

---

## B6.3 · `fig-map.html` · `viz-b6-map` · port of v1 `viz-b-propagator-map`

**Title.** "Spectral function of the field, and charge from current" `<span class="tag">interactive</span>`.

**What the reader must learn.**
- *Tab 1, Spectral map (default).* The spectral function of the field, $\frac{\omega^2}{4\pi c}\mathrm{Im}\,D$, in the transverse and the longitudinal channel for four media. Vacuum: weight only on the light line, nothing in the L channel. Free-electron (Drude) metal: no transverse wave below $\omega_{\mathrm p}$, transverse branch $\omega^2=\omega_{\mathrm p}^2+c^2q^2$; in L a flat plasmon at $\omega_{\mathrm p}$. Two-level insulator: exciton-like polaritons; flat L line at $\omega_{\mathrm L}=\sqrt{\varepsilon^{\mathrm e}(0)}\Delta$. Polar crystal: polaritons around a dark Reststrahlen band; in L one flat LO line, the loss function. The longitudinal channel never shows the light line: it is Coulomb's law dressed by matter.
- *Tab 2, Charge from current.* The connection formula (B6.7) holds only with the diamagnetic term. For one polar mode whose current makes an angle $\theta$ with $\bq$: the charge response computed directly from the bound charge and the one obtained from the full $\Pi^{\mathrm{ion}}$ through (B6.7) coincide for every $\theta$; the paramagnetic part alone adds $-\cos^2\theta\,\omega^2_{\mathrm{TO}}/\omega^2$, which diverges as $\omega\to0$ and adds a wrong-sign Drude term to $\varepsilon^{\mathrm L}$. At $\theta=90°$ the mode current is transverse and carries no charge.
- The v1 "Diagrams" tab is **not** ported here: it becomes `viz-b6-sectors`.

**Include position.** In `<h4 id="b6-modes">`, after the derivation box on the polariton branches, following the paragraph: "The modes are poles of $D$, so they appear as peaks in its imaginary part, the spectral function of the field. Figure B6.3 maps it in both channels for four media: vacuum, a free-electron metal, a model insulator and the polar crystal." It is followed by the caveat "What the isotropic model leaves out". Tab 2 is cited earlier, from the derivation box "the ionic charge response, directly and through the connection formula" in `#b6-continuity`, as "the second tab of Figure B6.3".

**Physics** (keep v1's numerics, rename the symbols).
- *Tab 1.* Local isotropic response, $\varepsilon^{\mathrm L}=\varepsilon^{\mathrm T}=\varepsilon(\omega)$, in (B6.17) `b6-eq-DLT`: $D^{\mathrm T}=-\frac{4\pi c}{\omega^2\varepsilon-c^2q^2}$, $D^{\mathrm L}=-\frac{4\pi c}{\omega^2\varepsilon}$. Plotted: $\log_{10}$ of $\frac{\omega^2}{4\pi c}\mathrm{Im}\,D$, i.e. $\mathrm{Im}[-\omega^2/(\omega^2\varepsilon-c^2q^2)]$ (T) and the loss function $\mathrm{Im}(-1/\varepsilon)$ (L). For $\omega\gt0$ the spectral weight is $+\mathrm{Im}\,D$ with this section's sign of $D$ ($\bA=D\bj^{\mathrm{ext}}$). Models of $\varepsilon-1=\frac{4\pi c}{\omega^2}\Pi$:
  - vacuum: $0$, with a small $\eta$ for visibility;
  - Drude: $-\omega_{\mathrm p}^2/[\omega(\omega+\ii\gamma)]$, default $\hbar\omega_{\mathrm p}=5.7$ eV (sodium);
  - two-level insulator: $(\varepsilon^{\mathrm e}(0)-1)\Delta^2/(\Delta^2-\omega^2-\ii\gamma\omega)$, defaults $\Delta=3.4$ eV, $\varepsilon^{\mathrm e}(0)=4$ (one oscillator standing in for all interband transitions, so its L–T splitting $(\sqrt{\varepsilon^{\mathrm e}(0)}-1)\Delta$ is far larger than a real exciton's);
  - polar crystal: $\varepsilon_\infty-1+(\omega^{\mathrm{ion}}_{\mathrm p})^2/(\omega^2_{\mathrm{TO}}-\omega^2-\ii\gamma\omega)$ with $(\omega^{\mathrm{ion}}_{\mathrm p})^2=\varepsilon_\infty(\omega^2_{\mathrm{LO}}-\omega^2_{\mathrm{TO}})$, presets as in B6.2.
  - $\gamma$ is a fraction of $\omega_0=\omega_{\mathrm{TO}}$, $\omega_{\mathrm p}$, $\Delta$ or 1 eV. Keep v1's pixel-integrated Lorentzian (poles narrower than a pixel stay visible) and its 200 × 200 grid.
  - Pole lines at $\gamma=0$: T channel, polar: (B6.20); Drude: $\omega^2=\omega_{\mathrm p}^2+c^2q^2$; two-level: roots of $\omega^4-(\omega_{\mathrm L}^2+c^2q^2)\omega^2+c^2q^2\Delta^2=0$ with $\omega_{\mathrm L}^2=\varepsilon^{\mathrm e}(0)\Delta^2$. L channel: the zeros of $\varepsilon$, $\omega_{\mathrm{LO}}$, $\omega_{\mathrm p}$ or $\omega_{\mathrm L}$ (B6.18).
- *Tab 2.* One polar mode, current at angle $\theta$ to $\bq$, $\gamma=0.01\,\omega_{\mathrm{TO}}$ (retarded). Unit of the charge response: $|\Pi_{\rho\rho}(\omega\to0)|$ at $\theta=0$, i.e. $q^2|\bj_\nu|^2/(\OV\omega^2_{\mathrm{TO}})$. With $\mathrm{den}=\omega^2-\omega^2_{\mathrm{TO}}+\ii\gamma\omega$:
  - direct, from the bound charge (derivation box in `#b6-continuity`): $\cos^2\theta\,\omega^2_{\mathrm{TO}}/\mathrm{den}$;
  - (B6.7) with the full $\Pi^{\mathrm{ion}}$ (`b5-eq-Piion`): $-\frac{\omega^2_{\mathrm{TO}}}{\omega^2}\cdot\big(-\omega^2\cos^2\theta/\mathrm{den}\big)$, identical;
  - (B6.7) with the paramagnetic part alone: $-\frac{\omega^2_{\mathrm{TO}}}{\omega^2}\cdot\big(-\omega^2\cos^2\theta/\mathrm{den}+\cos^2\theta\big)$, since $\Pi^{\mathrm{ion,dia}}$ (`b5-eq-dia`) contributes $-\cos^2\theta$ in these units;
  - right panel: $\mathrm{Re}\,\varepsilon^{\mathrm L}=\varepsilon_\infty+(\omega^{\mathrm{ion}}_{\mathrm p})^2\,\Pi^{\mathrm L}/\omega^2$ in the same units, for the full and the paramagnetic-only $\Pi$; LO point $\omega_{\mathrm{LO}}(\theta)=\sqrt{\omega^2_{\mathrm{TO}}+\cos^2\theta\,(\omega^{\mathrm{ion}}_{\mathrm p})^2/\varepsilon_\infty}$.
  - v1's formulas for these curves are correct; only names change.

**What to draw.** As v1: tab 1, a heat map ($q$ in μm⁻¹ on $x$, $\hbar\omega$ on $y$, colour bar labelled "$\omega^2\,\mathrm{Im}\,D/4\pi c$", sequential ramp `'c5'`), the light line $\omega=cq$ (photon colour, dashed), the $\gamma=0$ pole lines (ink, short dash), reference lines $\omega_{\mathrm{TO}}$, $\omega_{\mathrm{LO}}$ (phonon colour) or $\omega_{\mathrm p}$, $\Delta$, $\omega_{\mathrm L}$ (electron colour), labels placed so they never collide. Tab 2, two plots side by side (stacked on phones): left the three charge-response curves (direct: phonon colour, thick; full route: ink, dashed on top of it; paramagnetic only: `'c5'`), with $\mathrm{Im}$ of the direct curve as a light fill and a vertical $\omega_{\mathrm{TO}}$ marker; right $\mathrm{Re}\,\varepsilon^{\mathrm L}$ for both routes, the band between $\omega_{\mathrm{TO}}$ and $\omega_{\mathrm{LO}}(\theta)$ shaded, $\varepsilon_\infty$ as a dashed reference, and the LO zero marked.

**Controls.**
- Tab selector (seg): "Spectral map" / "Charge from current" (navigation, not counted).
- Tab 1 (3 to 5 visible): model (seg: vacuum, Drude, insulator, polar; default polar); channel (seg: "transverse, $\mathrm{Im}\,D^{\mathrm T}$" / "longitudinal, $\mathrm{Im}\,D^{\mathrm L}$"; default T); damping $\gamma/\omega_0$ (log range 0.3 %–10 %, default 2 %); one model-specific group: polar → crystal (select); Drude → $\hbar\omega_{\mathrm p}$ (1–10 eV, step 0.1); insulator → $\Delta$ (1–5 eV) and $\varepsilon^{\mathrm e}(0)$ (1.5–10). Draw the $\gamma=0$ pole lines always (drop v1's toggle, to stay within five controls).
- Tab 2 (3): crystal (select); angle $\theta$ between $\bq$ and the mode current (0–90°, default 0); toggle "show the paramagnetic-only route" (default on).

**Readouts.** Tab 1 polar: $\hbar\omega_{\mathrm{TO}}$, $\hbar\omega_{\mathrm{LO}}$, $\varepsilon_{\mathrm s}=\varepsilon_\infty\omega^2_{\mathrm{LO}}/\omega^2_{\mathrm{TO}}$, $\hbar\omega^{\mathrm{ion}}_{\mathrm p}$, the mixing wavevector $\sqrt{\varepsilon_\infty}\,\omega_{\mathrm{TO}}/c$ (μm⁻¹) and its wavelength, $\hbar\gamma$. Drude: $\hbar\omega_{\mathrm p}$, plasma wavelength $2\pi c/\omega_{\mathrm p}$ (nm), $\hbar\gamma$. Insulator: $\Delta$, $\hbar\omega_{\mathrm L}=\sqrt{\varepsilon^{\mathrm e}(0)}\,\Delta$, L–T splitting, $\hbar\gamma$. Vacuum: "$\Pi=0$", $\hbar\eta$, "$D^{\mathrm L}=-4\pi c/\omega^2$ is real: no spectral weight". Tab 2: $\Pi_{\rho\rho}(\omega\to0)$ in the unit above ($-\cos^2\theta$; paramagnetic only: $\to-\infty$), $\varepsilon^{\mathrm L}(0)=\varepsilon_\infty+\cos^2\theta\,(\omega^{\mathrm{ion}}_{\mathrm p})^2/\omega^2_{\mathrm{TO}}$ (paramagnetic only: $\to+\infty$), $\hbar\omega_{\mathrm{LO}}(\theta)$ ("none" at 90°).

**Initial state.** Tab 1, polar crystal, GaAs, transverse channel, $\gamma=2$ %: the anticrossing, the dark band and the light line are visible at once.

**Caption key points** (one caption with a part per tab, shown with the tab, as in v1; each part 200–260 words with the four-part structure).
- *Tab 1. What to notice.* In the transverse channel of the polar crystal, light and the TO phonon anticross into two polaritons around a dark band, and the upper branch starts at $\omega_{\mathrm{LO}}$; the longitudinal channel has one flat line at $\omega_{\mathrm{LO}}$ (or $\omega_{\mathrm p}$) for every $q$ and nothing at the light line. *Marks.* Colour scale, light line, pole lines, reference lines. *Try.* Switch channel; choose Drude and see the gap below $\omega_{\mathrm p}$ in T and the flat plasmon in L; lower $\gamma$ and see the poles sharpen. *Model.* (B6.17), (B6.18)–(B6.20); the four models of $\varepsilon$; local response; presets (literature, illustrative); the factor $\omega^2$ removes the temporal-gauge $1/\omega^2$ of $D^{\mathrm L}$, so the L panel is the loss function, which vanishes where $\varepsilon\to\infty$ (at $\omega_{\mathrm{TO}}$ or $\Delta$); at $q=0$ the channels coincide.
- *Tab 2. What to notice.* The direct charge response and the one from the full current correlator through (B6.7) lie on top of each other for every angle; without the diamagnetic term the result diverges at $\omega\to0$. *Marks.* The three curves and the $\varepsilon^{\mathrm L}$ panel. *Try.* Set $\theta=90°$: no charge response at all (a transverse current carries no charge); toggle the paramagnetic route. *Model.* One mode of $\Pi^{\mathrm{ion}}$ (`b5-eq-Piion`) with $\Pi^{\mathrm{dia}}$ (`b5-eq-dia`); the direct route from the derivation box of `#b6-continuity`; three routes computed by separate formulas, so their agreement is a real check.

**What in v1 is wrong or in source notation, and must change.**
- Id `viz-b-propagator-map` → `viz-b6-map`; classes `b06-*` → `b6-map-*`. The v1 script also builds the inline block figure `b-propagator-blockfig` and the Diagrams tab from the same `chain()` helper: neither belongs in this port (see B6.1).
- Labels and strings: `'Π_{jj} model'` → "Π model"; `'transverse Im D_T'`, `'longitudinal Im D_L'` → "transverse, Im D^T", "longitudinal, Im D^L"; `'ε_0 = ε_∞ω_{LO}^2/ω_{TO}^2'` → `ε_s = ε_∞ω_{LO}^2/ω_{TO}^2`; `'ħΩ_p'` → `ħω_p^{ion}`; `'static ε_e(0)'` → `static ε^e(0)`; `'Π_{ρρ} from Π_{jj}'` → "Charge from current"; `'E.350 with Π_{jj} = para + dia'` → "(B6.7) with the full Π"; `'E.350 with Π^{para} only'` → "(B6.7) with the paramagnetic part only"; `'direct ⟨T PP⟩ (p. 77)'` → "direct, from the bound charge −∇·P"; `'Re ε_L'` → `Re ε^L`; `'ε_L, gauge-invariant Π_{jj}'` → `ε^L, full Π`; `'ε_L = 0, γ = 0'` → `ε^L = 0, γ = 0`; `'ω_{LO}(θ)'` stays.
- Caption: replace $\Omega_p$, $\epsilon$, $\epsilon_\infty$, $\epsilon_{\rm e}(0)$, $\Pi_{jj}$, $\Pi_L$, $D_T$, $D_L$, $\Omega_{\rm uc}$, $\bj_\sigma$ and every "E Eq. …", "E.352′", "E.297′", "EPP p. 77" by the v2 symbols and eqrefs above. The branch index is $\nu$, never $\sigma$.
- The polar readout `'polariton crossover √ε_0 ω_{TO}/c'` used the static constant; the section uses $\sqrt{\varepsilon_\infty}\,\omega_{\mathrm{TO}}/c$ (where the two branches are exactly half phonon and half light; 0.557 μm⁻¹ in GaAs). Change label and value.
- Drop the source remark on EPP's time-ordered $D$ from the caption; keep only "for $\omega\gt0$ the spectral weight is $+\mathrm{Im}\,D$".
- v1's "poles at γ = 0" toggle: drop it and always draw the lines (control budget).
- Physics checks done for this spec: the T-channel quantity, the L-channel loss function, the two-level and polar pole formulas and all tab-2 formulas in the v1 code agree with (B6.17)–(B6.20) and the derivation box; keep them.
