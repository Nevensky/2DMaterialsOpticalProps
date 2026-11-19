# Section A.6 · Figure specifications

Five figures for `v2/sections/A06-phonon-se.html` (section id `a-phonon-se`, prefix `a6`). Read
`v2/build/BUILDER_GUIDE.md` §5 first, and the section itself: every symbol below is defined there and must be used
exactly as written (the glossary of the section lists them). Figure numbers are assigned in order of appearance:
`dyson` = A6.1, `bubble` = A6.2, `fermi` = A6.3, `window` = A6.4, `kohn` = A6.5.

Common rules for all five:
- File `v2/sections/A06/fig-<name>.html`, root `<figure class="viz" id="viz-a6-<name>">`, CSS classes prefixed
  `a6-<name>-`, script guarded by `if (!fig || !window.EPP) return;`.
- v2 notation only (no Giustino symbols: no $\mathcal D$, $D^{\rm A}$, $\Pi^{\rm NA}$, $g^{\rm b}$ as a function of $\br$,
  $\Gamma$ as a vertex, $\varepsilon_{n\bk}$, $\Omega_p$, $\Omega$ for a renormalized frequency, $G$ for a coupling
  constant, $W$ for a reduced frequency, $Q$ for $q/k_{\mathrm F}$, $A_{\rm cell}$). $\sigma$ is reserved for the
  conductivity. Gaussian units.
- Equations are cited with `<a class="eqref" href="#a6-eq-…"></a>`; the tags below are the current ones.
- Canvas text through `EPP.drawRich` (`ω_{LA}`, `k_F`, `Π^{ph,NA}`), bold vectors with `{bold:true}` (**q**, **k**).
- Caption 200–260 words: What to notice / Marks / Try / Model.
- Test with `bash v2/build/figtest.sh v2/sections/A06/fig-<name>.html` and look at both screenshots.

Equation ids and tags of the section (for the captions):

| tag | id | content |
|---|---|---|
| A6.1 | `a6-eq-B` | spectral function $\mathcal B_{\bq\nu}(\omega)=-\frac1\pi\mathrm{Im}F_{\nu\nu}$, free line |
| A6.2 | `a6-eq-resp` | oscillator response $\langle\hat\phi\rangle=-F^0h$ |
| A6.3 | `a6-eq-osc` | $F^{-1}=(F^0)^{-1}-\Pi$, $\Pi=\hbar\,\delta K/2M\omega_0$ |
| A6.4 | `a6-eq-dysonF` | Dyson equation $F=F^0+F^0\Pi^{\ph}F$ (export) |
| A6.5 | `a6-eq-diag` | diagonal $F_{\nu\nu}=\hbar^{-1}2\omega_{\bq\nu}/(\omega^2-\omega_{\bq\nu}^2-2\omega_{\bq\nu}\Pi^{\ph}_{\nu\nu}/\hbar)$ |
| A6.6 | `a6-eq-pole` | pole equation $\omega^2=\omega_{\bq\nu}^2+2\omega_{\bq\nu}\Pi^{\ph}_{\nu\nu}(\bq,\omega)/\hbar$ |
| A6.7 | `a6-eq-dphi` | $\Pi^{\ph}$ from a force-constant correction $\delta\Phi$ |
| A6.8 | `a6-eq-shift` | first order: $\Delta\omega_{\bq\nu}=\mathrm{Re}\Pi^{\ph}_{\nu\nu}(\omega_{\bq\nu})/\hbar$, $\gamma_{\bq\nu}=-\mathrm{Im}\Pi^{\ph}_{\nu\nu}(\omega_{\bq\nu})/\hbar$ (export) |
| A6.9 | `a6-eq-lor` | Lorentzian line of half width $\hbar\gamma_{\bq\nu}$ |
| A6.10 | `a6-eq-tau` | lifetime $\tau^{\ph}_{\bq\nu}=1/2\gamma_{\bq\nu}$ |
| A6.11 | `a6-eq-pi` | split $\Pi^{\ph}=\Pi^{\mathrm{ph,A}}+\Pi^{\mathrm{ph,NA}}$ (export) |
| A6.12 | `a6-eq-hf` | Hellmann–Feynman force constants |
| A6.13 | `a6-eq-chi` | $\Pi^{\mathrm{ph,NA}}$ with $\chi(\omega)-\chi(0)$ and two bare vertices |
| A6.14 | `a6-eq-dysonNA` | $F=F^0+F^0\Pi^{\mathrm{ph,NA}}F$ (DFPT reference) |
| A6.15 | `a6-eq-bs` | one bare vertex, $\chi^{\mathrm{irr}}$, one screened vertex $\mathcal V_{\bq\nu}$, minus $\omega\to0$ |
| A6.16 | `a6-eq-count` | $\sum_{a,b}\chi^{\mathrm{irr}}(v\chi^{\mathrm{irr}})^{a+b}=\sum_n(n+1)\chi^{\mathrm{irr}}(v\chi^{\mathrm{irr}})^n$ |
| A6.17 | `a6-eq-pina` | working formula with occupations (export) |
| A6.18 | `a6-eq-gold` | golden-rule linewidth |
| A6.19 | `a6-eq-allen` | Allen's double-delta formula |
| A6.20 | `a6-eq-pair` | pair energy $E_{\bk+\bq}-E_\bk=(\hbar^2/m^*)(\bk\cdot\bq+q^2/2)$ |
| A6.21 | `a6-eq-edges` | continuum edges $\hbar\omega^{\max}(q)$, $\hbar\omega^{\min}(q)$ |
| A6.22 | `a6-eq-chord` | $\mathrm{Im}\chi^0=-(N_{\mathrm F}/\OS)L_{\mathrm{res}}/q$ (2D) |
| A6.23 | `a6-eq-imchi` | small-$\omega$ $\mathrm{Im}\chi^0$ |
| A6.24 | `a6-eq-2d` | doped-sheet model $\Pi^{\mathrm{ph,NA}}=|g|^2\OS[\chi^0(q,\omega)-\chi^0(q,0)]$ |
| A6.25 | `a6-eq-stiff` | stiffening $\hbar\Delta\omega_\bq\simeq2N_{\mathrm F}|g|^2$ above the continuum |
| A6.26 | `a6-eq-landau` | Landau-damping width inside the continuum |
| A6.27 | `a6-eq-saitta` | $\hbar\Delta\omega_{\bq\to0,\nu}\simeq2N_{\mathrm F}\langle|g_{nn\nu}|^2\rangle_{\mathrm{FS}}$ |
| A6.28 | `a6-eq-wp` | $(\omega^{\ion}_{\mathrm p})^2=4\pi Z^2e^2/M\OV$ |
| A6.29 | `a6-eq-jpi` | jellium self-energy $\Pi^{\ph}=\frac{\hbar\omega^{\ion}_{\mathrm p}}2[1/\varepsilon^{\rme}(q,\omega)-1]$ |
| A6.30 | `a6-eq-jdisp` | $\omega_{\mathrm{LA}}(q)^2=(\omega^{\ion}_{\mathrm p})^2/\varepsilon^{\rme}(q,0)$ |
| A6.31 | `a6-eq-bohm` | Bohm–Staver $v_{\mathrm s}=\omega^{\ion}_{\mathrm p}/k_{\mathrm{TF}}=v_{\mathrm F}(Zm_{\mathrm e}/3M)^{1/2}$ |
| A6.32 | `a6-eq-lindhard` | static Lindhard $\varepsilon^{\rme}(q,0)$, $x=q/2k_{\mathrm F}$ |
| A6.33 | `a6-eq-kohn` | slope singularity at $2k_{\mathrm F}$ (export) |
| A6.34 | `a6-eq-jdamp` | $\gamma_\bq/\omega_{\mathrm{LA}}\simeq(\pi/4)v_{\mathrm s}/v_{\mathrm F}$ |

Conventions fixed by the section: $\gamma_{\bq\nu}$ is the **half width at half maximum** (HWHM, full width $2\gamma$);
$\bar\omega_{\bq\nu}=\omega_{\bq\nu}+\Delta\omega_{\bq\nu}$ is the renormalized frequency; $N_{\mathrm F}$ is the density of
states at $E_{\mathrm F}$ **per spin and per cell**; spectral functions are in eV⁻¹ (meV⁻¹ on axes); $\hbar=0.6582119569$
meV ps, $\hbar^2/2m_{\mathrm e}=3.80998$ eV Å², 1 meV = 8.06554 cm⁻¹.

---

## 1. `fig-dyson.html` · `viz-a6-dyson` · Dyson equation and the moving pole

- **Kind.** Port of `viz-a-phonon-se-fd-dyson` (v1 file `doc/epp-primer/sections/A06-phonon-se.html`, the
  `<figure>` near line 31 and the first `<script>` after the section, items `fd-a-phonon-se-d1…d3`), extended with two
  interactive panels. Title: "The Dyson equation moves the phonon pole" `<span class="tag">interactive</span>`.
- **What the reader must learn.**
  1. The Dyson equation as a picture: $F=F^0+F^0\Pi^{\ph}F$ and its series $F^0+F^0\Pi^{\ph}F^0+\dots$
  2. Any finite number of terms leaves the line at the old frequency $\omega_{\bq\nu}$ (the partial sums have poles of
     higher order there and blow up); only the resummed series moves the pole, to $\bar\omega_{\bq\nu}$, and gives it
     the half width $\gamma_{\bq\nu}=-\mathrm{Im}\,\Pi^{\ph}/\hbar$.
  3. Retarded poles lie in the lower half of the complex frequency plane; $\mathrm{Re}\,\Pi^{\ph}$ moves them along the
     real direction, $\mathrm{Im}\,\Pi^{\ph}$ moves them down.
- **Include position.** End of the subsection `a6-oscillator`, after the paragraph that ends "… <figref dyson> shows
  what a self-energy of this kind, a static correction plus friction, does to the spectral line." and before
  `<h3 id="a6-dyson">`. The sentences just before: "An imaginary correction describes friction. … A negative imaginary
  part of the self-energy means that the oscillator loses energy to whatever produces the correction."
- **Physics.** One mode of energy $E_0=\hbar\omega_{\bq\nu}$; work in energy units $E=\hbar\omega$ (meV).
  - Free propagator (A6.5 with $\Pi=0$, display broadening $\hbar\eta=0.4$ meV fixed):
    $F^0(E)=2E_0/[(E+\ii\hbar\eta)^2-E_0^2]$ (meV⁻¹).
  - Self-energy model, causal (real part even, imaginary part odd in $\omega$): $\Pi(E)=R+\ii I\,E/E_0$, with
    $R=\mathrm{Re}\,\Pi^{\ph}(\omega_{\bq\nu})$ and $I=\mathrm{Im}\,\Pi^{\ph}(\omega_{\bq\nu})\le0$ the controls. It is a
    static spring correction plus friction, (A6.3).
  - Dressed: $F(E)=1/[1/F^0(E)-\Pi(E)]$ (A6.4); partial sums $S_N(E)=F^0\sum_{n=0}^{N-1}(\Pi F^0)^n$ (scalar).
  - Spectral functions $\mathcal B=-\frac1\pi\mathrm{Im}F$ and $-\frac1\pi\mathrm{Im}S_N$ (A6.1).
  - Exact poles (with $\eta\to0$): $E=\ii I\pm\sqrt{E_0^2+2E_0R-I^2}$, so $\hbar\gamma=-I$ exactly and
    $\hbar\bar\omega=\sqrt{E_0^2+2E_0R-I^2}$ (A6.6). First order (A6.8): $\hbar\bar\omega\simeq E_0+R$, $\hbar\gamma\simeq-I$.
  - Lifetime (A6.10): $\tau^{\ph}=\hbar/(2\hbar\gamma)=0.3291\ \text{ps}/(\hbar\gamma\,[\text{meV}])$.
  - Convergence of the series at energy $E$: $|F^0(E)\Pi(E)|<1$.
  - Illustrative values; $E_0=50$ meV is a typical optical phonon.
- **What to draw.**
  - (a) Diagrams with `EPP.fd`, two rows. Row 1: thick double zigzag ($F$) = thin zigzag ($F^0$) + thin zigzag – grey
    blob labelled `Π^{ph}` – thick double zigzag. Row 2: the series `F^0 + F^0Π^{ph}F^0 + F^0Π^{ph}F^0Π^{ph}F^0 + …`
    drawn with the same glyphs; the first $N$ terms at full opacity, the rest faint. Labels `F`, `F^0`, `Π^{ph}` (upright
    "ph").
  - (b) Plot of the spectral function versus $\hbar\omega$, fixed $x$ range 0–100 meV, fixed $y$ range −0.4 to 1.0
    meV⁻¹ (clip). Thin grey: free line. Dashed (second colour): partial sum $S_N$, clipped where it leaves the frame,
    with a small label "first N terms". Bold phonon colour: Dyson result. Vertical dashed lines at $\hbar\omega_{\bq\nu}$
    (muted) and $\hbar\bar\omega_{\bq\nu}$ (phonon colour). Light shaded band where $|F^0\Pi|\ge1$, labelled "series
    diverges".
  - (c) Complex frequency plane: horizontal Re $\hbar\omega$ from −100 to 100 meV, vertical Im $\hbar\omega$ from −12 to
    4 meV (axes labelled, not to scale). Real axis drawn; upper half plane lightly shaded with the label "no poles
    (causality)". Open circles at $\pm E_0-\ii\hbar\eta$ (free poles), filled dots at the exact poles
    $\pm\hbar\bar\omega-\ii\hbar\gamma$, small crosses at the first-order estimate, arrows from free to dressed poles.
    Labels `ω_{qν}`, `ω̄_{qν} − iγ_{qν}`.
- **Controls (4).** $\mathrm{Re}\,\Pi^{\ph}(\omega_{\bq\nu})$: −10 to +10 meV, step 0.1, default +4.
  $\mathrm{Im}\,\Pi^{\ph}(\omega_{\bq\nu})$: −10 to 0 meV, step 0.1, default −3. Terms of the series $N$: 1–8 plus "all",
  default 2. $\hbar\omega_{\bq\nu}$: 20–100 meV, default 50.
- **Readouts.** $\hbar\omega_{\bq\nu}$; $\hbar\bar\omega_{\bq\nu}$ exact and first order; $\hbar\Delta\omega_{\bq\nu}$;
  $\hbar\gamma_{\bq\nu}$ (HWHM) and full width $2\hbar\gamma_{\bq\nu}$ (meV and cm⁻¹); $\tau^{\ph}_{\bq\nu}$ (ps);
  $|F^0\Pi^{\ph}|$ at $\hbar\bar\omega$.
- **Initial state.** $R=+4$ meV, $I=-3$ meV, $N=2$, $E_0=50$ meV: the Dyson line sits at 53.8 meV with half width 3 meV
  ($\tau=0.11$ ps), while the two-term sum shows a dispersive spike at 50 meV; the shaded divergence band covers the
  free pole.
- **Caption key points.**
  - *What to notice:* the summed propagator has one Lorentzian line at the shifted frequency with half width
    $\gamma=-\mathrm{Im}\Pi^{\ph}/\hbar$; any finite number of terms keeps the line at $\omega_{\bq\nu}$ and grows without
    bound there.
  - *Marks:* glyphs of (a); grey free line, dashed partial sum (clipped), bold Dyson result; shaded band
    $|F^0\Pi^{\ph}|\ge1$; open circles free poles, dots dressed poles, crosses first order; shaded upper half plane.
  - *Try:* raise $N$ (the peak never moves); set Im to 0 (sharp but shifted line); make Re negative (softening);
    compare exact and first-order readouts at $|\mathrm{Re}\,\Pi^{\ph}|=10$ meV.
  - *Model:* (A6.1), (A6.4)–(A6.6), (A6.8)–(A6.10); the self-energy model is a static spring correction plus friction,
    (A6.3); values illustrative.
- **What must change from v1.** v1 is a static diagram in Giustino's notation ($\mathcal D$, $\mathcal D^{\rm A}$,
  $\Pi^{\rm NA}$, title "(G Fig. 1a)") with a one-paragraph caption. Use $F$, $F^0$, $\Pi^{\ph}$ (the reference is a
  general $F^0$, as in (A6.4)), sentence-case title without source references, the four-part caption, and add panels
  (b) and (c). Keep the v1 glyph conventions (double zigzag for the dressed propagator, grey blob for the self-energy).

---

## 2. `fig-bubble.html` · `viz-a6-bubble` · One bare vertex, one screened vertex

- **Kind.** Port of `viz-a-phonon-se-fd-bubble` (v1 file `doc/epp-primer/sections/A06-phonon-se.html`, the
  `<figure>` near line 60 and items `fd-a-phonon-se-b1…b4` of the first script), extended with a chain row and a small
  counting panel. Title: "Inside the phonon self-energy: one bare vertex, one screened vertex"
  `<span class="tag">interactive</span>`.
- **What the reader must learn.**
  1. The lowest-order self-energy is an electron bubble with two bare vertices.
  2. The screening is a chain of bubbles joined by Coulomb lines.
  3. Summing the chain into one end gives one bare vertex, the irreducible bubble and one screened vertex, minus the
     same at $\omega=0$ (A6.15); the screened vertex obeys its own Dyson equation.
  4. Screening both ends counts a chain with $n$ Coulomb lines $n+1$ times (A6.16); leaving both ends bare keeps only
     $n=0$.
- **Include position.** In `a6-vertices`, after the paragraph that ends "… <figref bubble> draws the chain and what
  follows from it." and before "Summing the chain into one end of the bubble takes three steps." The sentences just
  before: "Equation (A6.13) has the bare vertex at both ends and the full response $\chi$ between them. The response is a
  chain: the electrons respond to a potential with the irreducible response $\chi^{\mathrm{irr}}$, the induced density
  creates a Hartree potential through $v$, …"
- **Physics.** (A6.13), the chain $\chi=\chi^{\mathrm{irr}}+\chi^{\mathrm{irr}}v\chi^{\mathrm{irr}}+\dots$, the screened
  vertex $\mathcal V_{\bq\nu}=(\varepsilon^{\rme})^{-1}\Delta_{\bq\nu}V^{\mathrm{en}}$, which satisfies
  $\mathcal V=\Delta V^{\mathrm{en}}+v\chi^{\mathrm{irr}}\mathcal V$, (A6.15), (A6.16). Counting panel, uniform gas,
  static: with $x=v\chi^{\mathrm{irr}}=1-\varepsilon^{\rme}$, the three choices sum to
  - both bare: $\chi^{\mathrm{irr}}$ (only $n=0$);
  - one bare, one screened (exact): $\chi^{\mathrm{irr}}\sum_nx^n=\chi^{\mathrm{irr}}/\varepsilon^{\rme}$;
  - both screened: $\chi^{\mathrm{irr}}\sum_n(n+1)x^n=\chi^{\mathrm{irr}}/(\varepsilon^{\rme})^2$;
  so, relative to the exact result, both bare is $\times\varepsilon^{\rme}$ and both screened $\times1/\varepsilon^{\rme}$
  (closed forms of the geometric series). Metal at small $q$: $\varepsilon^{\rme}\approx k_{\mathrm{TF}}^2/q^2$; at
  $q=k_{\mathrm{TF}}/3$, $\varepsilon^{\rme}=10$.
- **What to draw** (`EPP.fd`, rows stacked, each with a short left label):
  - (a) "lowest order": zigzag – small black dot – bubble of two thin fermion lines `G^0` – small black dot – zigzag;
    label under each dot `Δ_{qν}V^{en}`.
  - (b) "screening chain": zigzag – dot – bubble – dashed Coulomb line `v` – bubble – … – dot – zigzag, with $n$ Coulomb
    lines set by a control (0–4); bubbles labelled `χ^{irr}`.
  - (c) "exact, (A6.15)": zigzag – small black dot – bubble of thick lines `G` with an open circle `Λ` at its right
    vertex – grey disk `𝒱_{qν}` – zigzag, followed by `− [ same ]_{ω=0}`.
  - (d) "screened vertex": grey disk = small dot + small dot – bubble `χ^{irr}` – dashed `v` – grey disk.
  - (e) Bar chart: number of times a chain with $n$ Coulomb lines is counted, $n=0\dots6$, for the three choices (both
    bare: 1,0,0,…; one screened: 1,1,1,…; both screened: 1,2,3,…). The chosen series in full colour, the others faint.
- **Controls (3).** Chain length $n$ in row (b): 0–4, default 2. Choice: "both bare / one screened / both screened",
  default "one screened". $\varepsilon^{\rme}$ (static): 1–21, default 10.
- **Readouts.** Count of the chain drawn in (b) under the chosen rule; static self-energy relative to the exact one for
  the chosen rule ($\varepsilon^{\rme}$, 1 or $1/\varepsilon^{\rme}$).
- **Initial state.** $n=2$, "one screened", $\varepsilon^{\rme}=10$: readouts "counted once", "×1"; switching shows "×10"
  and "×0.1".
- **Caption key points.**
  - *What to notice:* the phonon self-energy needs exactly one screened end; screening both counts a chain with $n$
    Coulomb lines $n+1$ times and makes the result too small by $\varepsilon^{\rme}$.
  - *Marks:* zigzag external phonon line; thin lines Kohn–Sham propagators $G^0$, thick lines interacting $G$; small
    black dot bare vertex $\Delta_{\bq\nu}V^{\mathrm{en}}$; grey disk screened vertex $\mathcal V_{\bq\nu}$; open circle
    vertex $\Lambda$; dashed line Coulomb interaction $v$.
  - *Try:* lengthen the chain; switch to "both screened" and watch the bars grow as $n+1$; set $\varepsilon^{\rme}$ to 20,
    a metal at $q\approx0.2\,k_{\mathrm{TF}}$.
  - *Model:* (A6.13), (A6.15), (A6.16) and the screened vertex of (A5.28) (`a5-eq-Wph`); uniform-gas closed forms.
    Contrast in one sentence: the electron self-energy (Section A.7) has two screened vertices because its phonon line is
    internal.
- **What must change from v1.** v1 uses Giustino's labels $g^{\rm b}$, $g$, $\Gamma$, $G_0$, $\Pi^{\rm NA}$ and the title
  "(G Fig. 1b, 1c)". In v2: bare vertex $\Delta_{\bq\nu}V^{\mathrm{en}}$, screened vertex $\mathcal V_{\bq\nu}$, Hedin
  vertex $\Lambda$ ($\Gamma$ is the zone centre), $G^0$, $\Pi^{\mathrm{ph,NA}}$. Keep the explicit subtraction
  $-[\,\cdot\,]_{\omega=0}$ of v1. Add rows (b), (d) and the counting panel (e); four-part caption.

---

## 3. `fig-fermi.html` · `viz-a6-fermi` · Which electron–hole pairs a phonon can create (new)

- **Kind.** New, geometric. Title: "Electron–hole pairs: the Fermi disc and its copy shifted by **q**"
  `<span class="tag">interactive</span>`.
- **What the reader must learn.**
  1. The pairs $\bk\to\bk+\bq$ allowed by the Pauli principle end in the crescent between the Fermi disc and its copy
     shifted by $\bq$.
  2. Pairs of one energy lie on a straight line perpendicular to $\bq$; the phonon of energy $\hbar\omega$ can decay only
     through the part of that line inside the crescent, and $\mathrm{Im}\,\chi^0$ is its length (A6.22).
  3. Above $\hbar\omega^{\max}(q)$ the line misses the crescent: no decay, and the electrons cannot follow the phonon.
  4. Pairs of zero energy are the crossing points of the two circles; they exist for $q<2k_{\mathrm F}$, merge when the
     circles touch at $q=2k_{\mathrm F}$ and then disappear: the geometric origin of the Kohn anomaly.
- **Include position.** In `a6-pairs`, after the paragraph that ends "Translate the Fermi disc by $\bq$: the translated
  disc holds the final states $\bk+\bq$ of all occupied states $\bk$. The part of it outside the Fermi disc is allowed by
  the Pauli principle; the overlap is blocked." and before "Three results follow from the picture." The equation just
  before it is (A6.20).
- **Physics.** Two-dimensional parabolic band, $T=0$, intraband pairs. Reduced units: lengths in $k_{\mathrm F}$,
  energies in $E_{\mathrm F}=\hbar^2k_{\mathrm F}^2/2m^*$; write $\tilde q=q/k_{\mathrm F}$, $\tilde w=\hbar\omega/E_{\mathrm F}$;
  $\bq$ along $+x$.
  - Pair energy (A6.20): $(E_{\bk+\bq}-E_\bk)/E_{\mathrm F}=2\tilde k_x\tilde q+\tilde q^2$.
  - Resonant line, in final-state space ($\bk'=\bk+\bq$): $\tilde k'_x=c\equiv\tilde q/2+\tilde w/(2\tilde q)$ (vertical
    line); in initial-state space $\tilde k_x=k^\star/k_{\mathrm F}=\tilde w/(2\tilde q)-\tilde q/2$.
  - Allowed part of the line (final states outside the Fermi disc and inside the shifted disc):
    $\tilde k'^2_y\in\big(\max[0,1-c^2],\ 1-(c-\tilde q)^2\big)$; length
    $\tilde L=2\big[\sqrt{1-(c-\tilde q)^2}-\sqrt{\max(0,1-c^2)}\big]$ when $1-(c-\tilde q)^2>\max(0,1-c^2)$, else 0. Note
    $1-(c-\tilde q)^2=1-c^2+\tilde w$.
  - (A6.22): $\mathrm{Im}\,\chi^0/(N_{\mathrm F}/\OS)=-\tilde L/\tilde q$. Check against (A6.23) for small $\tilde w$:
    $-\tilde w/(\tilde q\sqrt{1-\tilde q^2/4})$.
  - Edges (A6.21): $\tilde w^{\max}=\tilde q^2+2\tilde q$, $\tilde w^{\min}=\max(0,\tilde q^2-2\tilde q)$.
  - Zero-energy pairs: final states $\tilde k'=(\tilde q/2,\pm\sqrt{1-\tilde q^2/4})$ for $\tilde q<2$, one point at
    $\tilde q=2$, none beyond.
- **What to draw.**
  - (a) $\bk$ plane, square, $k_x/k_{\mathrm F}$ from −1.6 to 4.4 and $k_y/k_{\mathrm F}$ from −2.2 to 2.2 (fixed, equal
    scales). Fermi disc $|\bk|<k_{\mathrm F}$ filled with the electron colour (light), its radius labelled `k_F`. The
    Fermi circle shifted by $\bq$ (centre at $\bq$), dashed. Crescent (shifted disc outside the Fermi disc) shaded
    "allowed"; lens (overlap) hatched grey "blocked". Bold arrow **q** from the origin. Six to eight parallel thin
    arrows $\bk\to\bk+\bq$ from sample occupied $\bk$: tips in the crescent coloured by pair energy (sequential scale,
    legend in units of $E_{\mathrm F}$), tips in the lens grey. The resonant line $k'_x=c\,k_{\mathrm F}$ dashed across
    the frame; its segment(s) inside the crescent drawn bold in the accent colour, labelled "pairs with energy ħω".
    Small dots at the zero-energy points; at $\tilde q=2$ a label "circles touch: q = 2k_F".
  - (b) Small $(q,\hbar\omega)$ map: $x=q/k_{\mathrm F}\in[0,3]$, $y=\hbar\omega/E_{\mathrm F}\in[0,3]$; continuum between
    $\tilde w^{\min}$ and $\tilde w^{\max}$ shaded; dashed line $\hbar\omega=\hbar v_{\mathrm F}q$ ($\tilde w=2\tilde q$,
    the small-$q$ edge); marker at the current $(\tilde q,\tilde w)$, filled when inside the continuum, hollow outside;
    region labels "no pairs: phonon cannot decay" above and "no pairs" below the lower edge.
- **Controls (3).** $q/k_{\mathrm F}$: 0.05–3, default 0.8. $\hbar\omega/E_{\mathrm F}$: 0–2.5, default 0.4. Toggle "sample
  arrows" (default on). Optionally a button "q = 2k_F" that sets $\tilde q=2$, $\tilde w=0.02$.
- **Readouts.** $q/k_{\mathrm F}$; $\hbar\omega/E_{\mathrm F}$; "inside / above / below the continuum";
  $\hbar\omega^{\min}/E_{\mathrm F}$ and $\hbar\omega^{\max}/E_{\mathrm F}$ at this $q$; $L_{\mathrm{res}}/k_{\mathrm F}$;
  $\mathrm{Im}\,\chi^0$ in units of $N_{\mathrm F}/\OS$; number of zero-energy pairs (2, 1, 0).
- **Initial state.** $\tilde q=0.8$, $\tilde w=0.4$: inside the continuum ($\tilde w^{\max}=2.24$), $c=0.65$,
  $\tilde L=0.458$, $\mathrm{Im}\,\chi^0=-0.57\,N_{\mathrm F}/\OS$; two short bold segments near the crossing points.
- **Caption key points.**
  - *What to notice:* a phonon can decay only through pairs on the resonant line inside the crescent; above
    $\hbar\omega^{\max}$ the line misses it; pairs of zero energy exist only while the circles cross and end at
    $q=2k_{\mathrm F}$.
  - *Marks:* filled disc occupied states; dashed circle the Fermi circle shifted by **q**; shaded crescent allowed final
    states; hatched lens blocked; arrows $\bk\to\bk+\bq$ coloured by pair energy; bold segment resonant pairs; dots
    zero-energy pairs; right panel the continuum (A6.21), the line $\hbar v_{\mathrm F}q$ and the current point.
  - *Try:* lower $q$ at fixed $\hbar\omega$ until the bold segment vanishes (the phonon leaves the continuum); set
    $\hbar\omega$ near 0 and move $q$ through $2k_{\mathrm F}$ (the dots merge and vanish); watch
    $\mathrm{Im}\,\chi^0$ grow as the dots approach tangency.
  - *Model:* (A6.20)–(A6.23), two-dimensional parabolic band at $T=0$, reduced units $k_{\mathrm F}$, $E_{\mathrm F}$;
    the doped sheet of <figref window> and the Kohn anomaly of <figref kohn> use the same geometry.

---

## 4. `fig-window.html` · `viz-a6-window` · Nonadiabatic phonon in a doped two-dimensional semiconductor

- **Kind.** Port of `viz-a-phonon-se-2d` (v1 file `doc/epp-primer/sections/A06-phonon-se.html`, `<figure>` near line 96
  and the second `<script>`). Title: "A doped sheet: hardening above the continuum, Landau damping inside"
  `<span class="tag">interactive</span>`.
- **What the reader must learn.** Where the phonon is faster than the electrons (above the continuum, small $q$) it is
  hardened by $2N_{\mathrm F}|g|^2$ and undamped; inside the continuum it is Landau damped and returns to the adiabatic
  frequency; beyond $2k_{\mathrm F}$ it is nearly free; the window shrinks toward $q\lesssim\omega_0/v_{\mathrm F}$ as the
  doping grows.
- **Include position.** In `a6-window`, after the paragraph with the numbers that ends "… <figref window> computes the
  full self-energy, the shift and width, and the spectral function." and before "Near the upper edge of the continuum,
  first order and the full pole equation (A6.6) differ strongly."
- **Physics** (keep the v1 code, which was checked against direct $\bk$ integration):
  - $k_{\mathrm F}=(2\pi n_{\mathrm e})^{1/2}$ ($n_{\mathrm e}$ in cm⁻² × 10⁻¹⁶ → Å⁻²); $E_{\mathrm F}=3809.98\ \text{meV Å}^2\times k_{\mathrm F}^2/(m^*/m_{\mathrm e})$; $m^*=0.5\,m_{\mathrm e}$ fixed.
  - (A6.24) with $\chi^0=(2N_{\mathrm F}/\OS)\,\tilde\chi$:
    $\Pi^{\mathrm{ph,NA}}(q,\omega)=C\,[\tilde\chi(x,\zeta)-\tilde\chi(x,0)]$, coupling energy $C=2N_{\mathrm F}|g|^2$ (meV).
    In code only: $x=q/k_{\mathrm F}$, $\zeta=(\hbar\omega+\ii\hbar\eta)/E_{\mathrm F}$,
    $\tilde\chi(x,\zeta)=-\big[1+\big(s(\zeta-x^2)+s(-\zeta-x^2)\big)/2x^2\big]$, $s(y)=y\sqrt{1-4x^2/y^2}$ (principal
    branch); static value $-1$ for $x\le2$ and $-(1-\sqrt{1-4/x^2})$ for $x>2$. (Do not show $\tilde\chi$, $\zeta$ or
    $s$ in the visible text: the tilde means mass-scaled in v2; describe it in words as "the Lindhard function in units
    of $2N_{\mathrm F}/\OS$".)
  - Spectral function from (A6.5) with a display broadening $\eta$:
    $\mathcal B_\bq(\omega)=-\frac1\pi\mathrm{Im}\,\dfrac{2\hbar\omega_0}{(\hbar\omega)^2-(\hbar\omega_0)^2+2\ii\hbar\eta\,\hbar\omega-2\hbar\omega_0\Pi^{\mathrm{ph,NA}}}$ (meV⁻¹).
  - Pole: "first order", (A6.8): $\hbar\Delta\omega_\bq=\mathrm{Re}\,\Pi^{\mathrm{ph,NA}}(\omega_0)$,
    $\hbar\gamma_\bq=-\mathrm{Im}\,\Pi^{\mathrm{ph,NA}}(\omega_0)$ (replace v1's square-root form); "pole equation", (A6.6)
    iterated on the real axis as in v1. $\eta\to0^+$ (1e-4 meV) for $\Delta\omega$, $\gamma$, $\tau$; the control $\eta$
    only broadens the map. $\tau^{\ph}_\bq=\hbar/2\hbar\gamma_\bq$ (A6.10).
  - Edges (A6.21): $\hbar\omega^{\max}=E_{\mathrm F}(x^2+2x)$, $\hbar\omega^{\min}=\max[0,E_{\mathrm F}(x^2-2x)]$; $\omega_0$ inside
    the continuum for $x\in[\sqrt{1+\hbar\omega_0/E_{\mathrm F}}-1,\ \sqrt{1+\hbar\omega_0/E_{\mathrm F}}+1]$.
  - Limits to check: (A6.25) $\hbar\Delta\omega_\bq\to C$ at small $q$; (A6.26) $\hbar\gamma_\bq\simeq C\,\omega_0/(v_{\mathrm F}q\sqrt{1-(q/2k_{\mathrm F})^2})$ for $\omega_0\ll v_{\mathrm F}q$.
  - Values (illustrative, cell area of MoS₂ $\OS=8.65$ Å², one valley): $n_{\mathrm e}=10^{13}$ cm⁻² gives
    $k_{\mathrm F}=0.0793$ Å⁻¹, $E_{\mathrm F}=47.9$ meV; $N_{\mathrm F}=0.090$ eV⁻¹ per spin; $C=5$ meV ↔ $|g|=0.17$ eV;
    window $0.43<q/k_{\mathrm F}<2.43$ for $\hbar\omega_0=50$ meV; at $q=k_{\mathrm F}$, $\hbar\gamma\approx5$ meV,
    $\tau\approx0.07$ ps (first order; (A6.26) gives 3 meV but needs $\omega_0\ll v_{\mathrm F}q$).
- **What to draw** (as v1): (a) heat map of $\mathcal B_\bq(\omega)$, $x$ = $q/k_{\mathrm F}\in[0,3]$, $y=\hbar\omega\in[0,2\hbar\omega_0]$
  meV, square-root colour scale, colour bar "𝓑 (meV⁻¹)", continuum edges `ħω^{max}`, `ħω^{min}` in the electron colour,
  line `ħω_0`, dashed $q$ cursor. (b) $\hbar\Delta\omega_\bq$ and $\hbar\gamma_\bq$ (HWHM) versus $q/k_{\mathrm F}$, band
  shaded where $\omega_0$ is inside the continuum. (c) $\tau^{\ph}_\bq$ (ps, log axis). Click or drag in any panel moves
  the cursor. Controls below the plots.
- **Controls (5).** Electron density $n_{\mathrm e}$: $10^{12}$–$10^{14}$ cm⁻², log, default $10^{13}$. $\hbar\omega_0$:
  10–200 meV, default 50. Coupling $2N_{\mathrm F}|g|^2$: 0–15 meV, default 5. Map broadening $\hbar\eta$: 0.5–5 meV,
  default 1. Pole: "first order (A6.8)" / "pole equation (A6.6)", default first order.
- **Readouts.** $k_{\mathrm F}$ (Å⁻¹), $E_{\mathrm F}$ (meV), $\hbar\omega_0/E_{\mathrm F}$, window in $q/k_{\mathrm F}$; at the
  cursor: $q/k_{\mathrm F}$, $\Pi^{\mathrm{ph,NA}}$ at $\omega_0$ (or at the pole), $\hbar\bar\omega_\bq$, $\hbar\gamma_\bq$,
  $\tau^{\ph}_\bq$ (∞ when undamped).
- **Initial state.** v1 defaults ($n_{\mathrm e}=10^{13}$ cm⁻², 50 meV, 5 meV, 1 meV, first order, cursor at
  $q=k_{\mathrm F}$): hardening of about 5 meV at small $q$, width of about 5 meV at the cursor.
- **Caption key points.**
  - *What to notice:* above the continuum the phonon is hardened by about $2N_{\mathrm F}|g|^2$ and undamped (the
    electrons cannot follow, so the static softening built into $\omega_0$ is undone); inside it is Landau damped and
    returns toward $\omega_0$; beyond $q\approx2.4\,k_{\mathrm F}$ it is nearly free.
  - *Marks:* colour map, edges, $\hbar\omega_0$, cursor; shift and width curves; shaded continuum band; lifetime.
  - *Try:* raise $n_{\mathrm e}$ until $E_{\mathrm F}\gg\hbar\omega_0$ (the window shrinks toward $q\to0$); switch to the
    pole equation near the upper edge (the hardened pole stays sharp just above it); set the coupling to zero.
  - *Model:* (A6.24)–(A6.26), (A6.5), (A6.6), (A6.8), (A6.10), (A6.21); parameters above; no plasmon (the response has
    no Coulomb chain) and one mode, so no mode mixing.
- **What must change from v1.** Notation: $Q$ → axis label "q/k_F"; $W$ → internal only; $G$ → $2N_{\mathrm F}|g|^2$
  ($G$ is the electron propagator in v2); $\hbar\Pi^{\rm NA}$ → $\Pi^{\mathrm{ph,NA}}$; $B(Q,\omega)$ → $\mathcal B_\bq(\omega)$;
  $\hbar\omega_\pm$ → $\hbar\omega^{\max}$, $\hbar\omega^{\min}$; $\hbar(\Omega-\omega_0)$ → $\hbar\Delta\omega_\bq$; $\Omega$ →
  $\bar\omega$; $\varepsilon_F$ → $E_{\mathrm F}$; $n$ → $n_{\mathrm e}$; $A_{\rm cell}$ → $\OS$; "(G.137)", "(G.138–139)" →
  (A6.6), (A6.8). First order uses (A6.8) instead of the square-root form. Four-part caption; ids and classes
  `a6-window-…`.

---

## 5. `fig-kohn.html` · `viz-a6-kohn` · Jellium phonons: Bohm–Staver sound and the Kohn anomaly

- **Kind.** Port of `viz-a-phonon-se-jellium` (v1 file `doc/epp-primer/sections/A06-phonon-se.html`, `<figure>` near
  line 121 and the third `<script>`). Title: "Phonons of the electron gas: sound and the Kohn anomaly"
  `<span class="tag">interactive</span>`.
- **What the reader must learn.** Static screening turns the ionic plasma oscillation into sound with the Bohm–Staver
  velocity; with the Lindhard function the dispersion looks smooth but its slope has a logarithmic spike at
  $2k_{\mathrm F}$, which Thomas–Fermi screening (no Fermi surface) lacks; in jellium $2k_{\mathrm F}$ lies beyond $q_{\mathrm D}$ for
  $Z=1$; nonadiabatic corrections are of order $v_{\mathrm s}/v_{\mathrm F}\sim10^{-3}$.
- **Include position.** In `a6-kohn`, after the paragraph that ends "This is the **Kohn anomaly** (Kohn, 1959).
  <figref kohn> shows it for sodium, potassium and aluminium." and before "The geometric reason is the one of
  <figref fermi>." The equations just before are (A6.32) and (A6.33).
- **Physics** (keep the v1 code; rename):
  - (A6.30) $\omega_{\mathrm{LA}}(q)^2=(\omega^{\ion}_{\mathrm p})^2/\varepsilon^{\rme}(q,0)$; (A6.28)
    $(\hbar\omega^{\ion}_{\mathrm p})^2=0.75642\,Z^2/(\OV[\text{Å}^3]\,M[\text{amu}])$ eV² (Gaussian $4\pi Z^2e^2\hbar^2/M\OV$).
  - $\OV=Z\,\frac{4\pi}3(r_{\mathrm s}a_0)^3$, $a_0=0.529177$ Å; $k_{\mathrm F}a_0=1.91916/r_{\mathrm s}$;
    $E_{\mathrm F}=(k_{\mathrm F}a_0)^2\times13.6057$ eV; $k_{\mathrm{TF}}^2=4k_{\mathrm F}/\pi a_0$; $\hbar v_{\mathrm F}=2E_{\mathrm F}/k_{\mathrm F}$;
    $q_{\mathrm D}=(2/Z)^{1/3}k_{\mathrm F}$.
  - Lindhard (A6.32): bracket $L(x)=\frac12+\frac{1-x^2}{4x}\ln|\frac{1+x}{1-x}|$, $x=q/2k_{\mathrm F}$, derivative
    $L'(x)=\frac1{2x}-\frac{1+x^2}{4x^2}\ln|\frac{1+x}{1-x}|$ (use these names in code only; do not call it $F$, which is the
    phonon propagator in v2). Thomas–Fermi: $\varepsilon^{\rme}=1+k_{\mathrm{TF}}^2/q^2$.
  - Group velocity (A6.33): $\dd\omega_{\mathrm{LA}}/\dd q=-(\omega_{\mathrm{LA}}/2\varepsilon^{\rme})\,\dd\varepsilon^{\rme}/\dd q$,
    $\dd\varepsilon^{\rme}/\dd q=-2k_{\mathrm{TF}}^2L/q^3+(k_{\mathrm{TF}}^2/q^2)L'(x)/2k_{\mathrm F}$; km/s = (eV Å) × 151.927.
  - Bohm–Staver (A6.31): $v_{\mathrm s}=\omega^{\ion}_{\mathrm p}/k_{\mathrm{TF}}=v_{\mathrm F}(Zm_{\mathrm e}/3M)^{1/2}$; damping
    (A6.34): $\gamma_\bq/\omega_{\mathrm{LA}}\simeq(\pi/4)v_{\mathrm s}/v_{\mathrm F}$.
  - Keep v1's grid refinement around $2k_{\mathrm F}$ so the spike is resolved.
  - Presets (free-electron $r_{\mathrm s}$, literature values, illustrative): Na $r_{\mathrm s}=3.93$, $Z=1$, $M=22.99$ amu;
    K 4.86, 1, 39.10; Al 2.07, 3, 26.98. Sodium gives $\hbar\omega^{\ion}_{\mathrm p}=29.6$ meV, $k_{\mathrm F}=0.923$ Å⁻¹,
    $E_{\mathrm F}=3.24$ eV, $k_{\mathrm{TF}}/k_{\mathrm F}=1.61$, $v_{\mathrm s}=3.01$ km/s, $v_{\mathrm s}/v_{\mathrm F}=2.8\times10^{-3}$,
    $\gamma/\omega=2.2\times10^{-3}$, $q_{\mathrm D}/k_{\mathrm F}=1.26$; K: 16.5 meV; Al: 124 meV.
- **What to draw** (as v1): (a) $\hbar\omega_{\mathrm{LA}}(q)$ (meV) versus $q/k_{\mathrm F}\in[0,3]$: Lindhard solid, Thomas–Fermi
  dashed, Bohm–Staver line $v_{\mathrm s}q$ dotted, horizontal line `ħω_p^{ion}` (bare ions), vertical `2k_F` (electron
  colour) and `q_D` (muted). (b) Group velocity $\dd\omega_{\mathrm{LA}}/\dd q$ (km/s): Lindhard solid with the spike,
  Thomas–Fermi dashed, horizontal `v_s`, vertical `2k_F` labelled "Kohn anomaly", `q_D`. Fixed axes scaled to the preset
  ($y$ up to 1.3 $\hbar\omega^{\ion}_{\mathrm p}$ and 1.35 $v_{\mathrm s}$, as v1).
- **Controls (4).** Metal preset (Na / K / Al / custom); $r_{\mathrm s}$ 1.5–6; $Z$ 1, 2, 3; $M$ 5–210 amu.
- **Readouts.** $k_{\mathrm F}$, $E_{\mathrm F}$, $k_{\mathrm{TF}}/k_{\mathrm F}$, $\hbar\omega^{\ion}_{\mathrm p}$, $v_{\mathrm s}$,
  $v_{\mathrm s}/v_{\mathrm F}=(Zm_{\mathrm e}/3M)^{1/2}$, $\gamma_\bq/\omega_{\mathrm{LA}}\simeq(\pi/4)v_{\mathrm s}/v_{\mathrm F}$, $q_{\mathrm D}/k_{\mathrm F}$.
- **Initial state.** Sodium.
- **Caption key points.**
  - *What to notice:* sound at small $q$ with $v_{\mathrm s}=\omega^{\ion}_{\mathrm p}/k_{\mathrm{TF}}$ (3 km/s for Na); the Lindhard
    dispersion looks smooth but its slope spikes logarithmically at $2k_{\mathrm F}$; Thomas–Fermi has no spike.
  - *Marks:* curves and lines above.
  - *Try:* switch Na → Al (denser gas and $Z=3$: $\hbar\omega^{\ion}_{\mathrm p}=124$ meV and faster sound, while
    $q_{\mathrm D}/k_{\mathrm F}$ drops to 0.87, so $2k_{\mathrm F}$ lies even further beyond the Debye sphere; in a crystal the
    spike folds back by reciprocal lattice vectors); raise $M$ (sound slows as $M^{-1/2}$ while the anomaly stays at
    $2k_{\mathrm F}$); lower $r_{\mathrm s}$ (denser gas, larger $k_{\mathrm F}$ and $E_{\mathrm F}$).
  - *Model:* (A6.28)–(A6.34), ions as a uniform fluid in a free-electron gas, static Lindhard or Thomas–Fermi
    screening; presets free-electron values; the geometric reason is drawn in <figref fermi>.
- **What must change from v1.** $\Omega_p$ → $\omega^{\ion}_{\mathrm p}$; $\epsilon_{\rm e}(q,0)$ → $\varepsilon^{\rme}(q,0)$;
  $\Omega$ (volume per ion) → $\OV$; $\omega_q$ → $\omega_{\mathrm{LA}}(q)$; the Lindhard function is not called $F$; units
  Gaussian (v1 said SI); readout "Landau damping γ/ω" → "$\gamma_\bq/\omega_{\mathrm{LA}}$ (A6.34)"; four-part caption; ids
  and classes `a6-kohn-…`.
