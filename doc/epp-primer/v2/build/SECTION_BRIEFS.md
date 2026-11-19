# Section briefs for primer v2, A.2 to B.10

Each brief lists:
- **Source**: the v1 file to rewrite (in `doc/epp-primer/sections/`), the research notes to consult
  (`doc/epp-primer/notes/`), and the source pages.
- **Cover**: what the v2 section must teach. The derivations named here are shown step by step.
- **Figures**:
  - *port* = rebuild a v1 figure (its v1 id is given) in v2;
  - *new* = a suggested geometric figure.

  The author may change the new ones, but adds at most two. Total figures per section: 2–5.
- **Exports**: equation ids that must exist, on the equations described. Other sections link to them with `a.eqref`.
- **May link**: exports of other sections that this section is expected to cite (any id in NOTATION §12 is also fine).
- **Care**: physics pitfalls and corrections.

Source texts:
- PDFs: `/root/.claude/uploads/e798d873-fe28-512b-a4f7-e6f2d29e481a/9d9c9111-giustino2017.pdf` and `…/608c6413-EPPv2.pdf`.
- Text extracts: `/tmp/claude-0/-home-user-2DMaterialsOpticalProps/e798d873-fe28-512b-a4f7-e6f2d29e481a/scratchpad/pdf/g.txt` and `epp.txt`, with one file per EPP page, `p001.txt`…
- The notes `G1–G6`, `E1–E9` and `SYNTHESIS.md` summarize the sources page by page, with equation numbers.

---

## A.2 · `A02-coupling.html` · id `a-coupling` · "The electron–phonon Hamiltonian and its matrix element"
- **Source.** v1 `A02-coupling.html`. Notes G1 (history) and G2 (§III.B). Giustino §§II, III.B, Eqs. (1), (31)–(38).
- **Cover.**
  - Electrons from scratch: Bloch states $\psi_{n\bk}$ normalized on the BvK supercell, bands $E_{n\bk}$, occupations
    $f_{n\bk}$, $\mu$. Second quantization goes in a background box: $c$, $c^\dagger$, one-body operators as
    $\sum c^\dagger c$.
  - Expand the crystal potential felt by an electron to first order in the A.1 displacements, insert the normal-mode
    expansion (A1.27), and derive $\hat H^{\mathrm{e\text{-}ph}}$ and $g$ step by step. Show:
    - where the factor $N_\ell^{-1/2}$ comes from;
    - why $\bk\to\bk+\bq$, using crystal-momentum conservation from lattice periodicity: the lattice sum, then the
      Bloch theorem;
    - why $g$ is a zero-point amplitude times a potential gradient projected on the mode.
  - Units and the typical size of $g$: tens to hundreds of meV.
  - Giustino's Eq. (1) read term by term: electron band energy, phonon energy, linear coupling, second-order
    (Debye–Waller) coupling, preview only.
  - What Eq. (1) leaves out.
  - A compact history: rigid ions, Bloch, Nordheim, Fröhlich, Bardeen's screened vertex in a metal, Pines, the
    deformation potentials of semiconductors, then DFPT.
  - Bardeen's result as the first screened vertex. A point ion in jellium screened by Thomas–Fermi or Lindhard gives
    $4\pi Ze^2/(q^2\varepsilon(q))$. This motivates A.3.
- **Figures.**
  - *new* `overlap`: a frozen phonon in a 1D chain or 2D sheet. It shows the potential change it creates, and the
    matrix element as the overlap of $\psi^*_{m\bk+\bq}$, $\Delta V$ and $\psi_{n\bk}$. A cosine band beside it shows
    the transition $\bk\to\bk+\bq$, so the selection rule is visible: the overlap vanishes unless the wavevectors match.
  - *port* `viz-a-coupling-bardeen` → `bardeen`.
  - *port* `viz-a-coupling-timeline` → `timeline`, kept short.
- **Exports.**
  - `a2-eq-dV`: definition of $\partial_{\kappa\alpha,\bq}V$ and $\Delta_{\bq\nu}V$;
  - `a2-eq-g`: definition of $g_{mn\nu}(\bk,\bq)$;
  - `a2-eq-H`: $\hat H^{\mathrm{e\text{-}ph}}$;
  - `a2-eq-H0`: the full Eq. (1)-type Hamiltonian;
  - `a2-eq-bardeen`: Bardeen's screened vertex.
- **May link.** A.1 exports.
- **Care.** In the summed form the phase $\ee^{\ii\bq\cdot\bR_\ell}$ sits on the displaced atom's cell. Keep $N_\ell^{-1/2}$ outside $g$.

## A.3 · `A03-dfpt.html` · id `a-dfpt` · "Density-functional perturbation theory: the vertex is screened"
- **Source.** v1 `A03-dfpt.html`. Notes G2 and G4 (DW). Giustino §III.C–E, Eqs. (39)–(60) and (61)–(74) as relevant.
- **Cover.**
  - The self-consistent change of $V^{\mathrm{KS}}$ when atoms move in a wave $\bq$:
    $\Delta V^{\mathrm{KS}}=\Delta V^{\mathrm{en}}+\Delta V^{\mathrm H}+\Delta V^{\mathrm{xc}}$, each term derived.
  - $\Delta n$ from first-order perturbation theory, then the Sternheimer equation, derived from it and from the
    projector onto empty states, so that no empty states are needed.
  - The loop as a geometric series: $\Delta V^{\mathrm{KS}}=\Delta V^{\mathrm{en}}+(v+f^{\mathrm{xc}})\chi^0\Delta V^{\mathrm{KS}}$,
    so $\Delta V^{\mathrm{KS}}=(\varepsilon^{\mathrm e})^{-1}\Delta V^{\mathrm{en}}$, with the meaning of each factor.
  - Test-electron and test-charge dielectric functions, and when each applies.
  - Jellium, where Bardeen's vertex is recovered (link to `a2-eq-bardeen`).
  - Static means adiabatic: why DFPT's frequency-independent screening assumes $\hbar\omega_{\bq\nu}\ll$ electronic
    excitation energies, and when that fails (metals, doped semiconductors; preview of A.6).
  - The second-order (Debye–Waller) matrix element and the rigid-ion approximation.
  - Typical numbers: the Thomas–Fermi wavevector of a metal, the screening of a silicon vertex.
- **Figures.**
  - *port* `viz-a-dfpt-scf` → `scf`: the loop as a geometric series.
  - *port* `viz-a-dfpt-cloud` → `cloud`: a displaced ion and its screening cloud. Make it geometric, with radial
    profiles, the bare and screened potential, and the induced charge.
  - Optionally, *new*: Sternheimer as a projection, a band diagram showing that $\Delta\psi$ lives in the empty bands.
- **Exports.**
  - `a3-eq-dvks`: the three parts of $\Delta V^{\mathrm{KS}}$;
  - `a3-eq-sternheimer`;
  - `a3-eq-chi0`: $\chi^0$ as a sum over states;
  - `a3-eq-scf`: $\Delta V^{\mathrm{KS}}=(\varepsilon^{\mathrm e})^{-1}\Delta V^{\mathrm{en}}$;
  - `a3-eq-eps`: $\varepsilon^{\mathrm e}$, test-electron and RPA;
  - `a3-eq-gdw`: $g^{\mathrm{DW}}$.
- **May link.** `a2-eq-g`, `a2-eq-dV`, `a2-eq-bardeen`.
- **Care.** Spin factor 2 in $\chi^0$. State the sign convention of $\chi^0$ (negative static response).

## A.4 · `A04-polar.html` · id `a-polar` · "Polar crystals: Born charges, the Fröhlich vertex and LO–TO splitting"
- **Source.** v1 `A04-polar.html`. Notes G5 (§VI polar term) and G2. Giustino §VI.B, Eqs. (176)–(180), the Fröhlich
  model in §II.A.3, and the LO–TO discussion. Gonze–Lee is cited in the review.
- **Cover.**
  - **Born effective charges** from scratch: the dipole created per unit displacement; why $Z^*\neq Z_\kappa$, because
    electrons follow; the tensor meaning; the acoustic sum rule $\sum_\kappa\mathbf Z^*_\kappa=0$, by the same argument
    as (A1.6).
  - The polarization of a lattice wave, $\bP=(e/\OV)\sum_\kappa\mathbf Z^*_\kappa\cdot\bu_\kappa$.
  - Gauss's law in $\bq$ space gives the macroscopic longitudinal field $\bE=-4\pi\hat\bq(\hat\bq\cdot\bP)/(\hat\bq\cdot\boldsymbol\varepsilon_\infty\cdot\hat\bq)$,
    for LO but not TO.
  - The electron's energy in that field gives $g^{\mathrm{lr}}$, step by step:
    - the $1/q$ divergence;
    - the sum over $\bG$, with $(\bq+\bG)$ and $\boldsymbol\varepsilon_\infty$;
    - the overlap of Bloch periodic parts;
    - the phase $\ee^{-\ii(\bq+\bG)\cdot\btau_\kappa}$.
  - Why polar crystals break Wannier interpolation, and the fix of subtracting and adding $g^{\mathrm{lr}}$.
  - The Verdi–Giustino vertex and its Fröhlich limit (one LO mode, isotropic).
  - $\alpha_{\mathrm F}$, with GaAs and LiF numbers.
  - The same kernel in the force constants, $\Phi^{\mathrm{dd}}$: the LO–TO splitting
    $\omega_{\mathrm{LO}}^2-\omega_{\mathrm{TO}}^2=4\pi e^2Z^{*2}/(\OV M_{\mathrm r}\varepsilon_\infty)$ for a two-atom
    cubic crystal, derived.
  - Where electrostatics stops: $q\lesssim\omega_{\mathrm{TO}}/c$, a preview of Part B.
- **Figures.**
  - *new* `born`: a displaced ion in a polar crystal and the dipole it creates, with the tensor columns
    $eZ^*_{\kappa,\cdot\beta}$ as arrows for displacements along $x$ and $y$.
  - *port* `viz-a-polar-dipoles` → `dipoles`: Born-charge dipoles of an LO wave and the field, with the TO comparison.
  - *port* `viz-a-polar-wannier` → `wannier`.
  - *port* `viz-a-polar-frohlich` → `frohlich`.
- **Exports.**
  - `a4-eq-born`: definition of $Z^*$;
  - `a4-eq-asr`: Born-charge sum rule;
  - `a4-eq-P`: polarization of a wave;
  - `a4-eq-E`: longitudinal field;
  - `a4-eq-glr`: the Verdi–Giustino vertex;
  - `a4-eq-frohlich`: Fröhlich's $g$;
  - `a4-eq-alphaF`;
  - `a4-eq-dd`: $\Phi^{\mathrm{dd}}$;
  - `a4-eq-lot`: LO–TO splitting.
- **May link.** A.1 exports, `a2-eq-g`, `a3-eq-scf`.
- **Care.** The phase on the atom, $\ee^{-\ii(\bq+\bG)\cdot\btau_\kappa}$, as in the review's Eq. (180). The EPP
  manuscript's polar vertex lacks it (Part B fixes this). Gaussian units with an SI footnote. $Z^*$ is dimensionless.

## A.5 · `A05-hedin-baym.html` · id `a-hedin-baym` · "Hedin–Baym: electrons and nuclei with one Coulomb interaction"
- **Source.** v1 `A05-hedin-baym.html`. Note G3. Giustino §IV, Eqs. (75)–(122), and Table I.
- **Cover.**
  - The Hamiltonian of electrons and nuclei with the instantaneous Coulomb interaction only.
  - Background boxes: the electron propagator $G$ and its poles; the Dyson equation $G=G^0+G^0\Sigma G$; Hedin's
    equations.
  - Screening: $W=v+v\chi^{\mathrm{irr}}W$ and $\varepsilon=1-v\chi^{\mathrm{irr}}$, in the χ notation of NOTATION §12.3.
  - Letting the nuclei move: derive $W=W^{\mathrm e}+W^{\mathrm{ph}}$, with
    $W^{\mathrm{ph}}$ = screened vertex × phonon propagator $F$ × screened vertex. Define $F$ here (NOTATION §12.3).
  - The vertex is the screened gradient of the nuclear potential. That is DFPT's $g$ (link A.3).
  - A polar insulator: $W^{\mathrm{ph}}$ is the Fröhlich interaction, attractive below $\omega_{\mathrm{LO}}$.
  - Jellium with mobile ions: $W=v/\varepsilon$ with $\varepsilon=1+k_{\mathrm{TF}}^2/q^2-(\omega^{\mathrm{ion}}_{\mathrm p})^2/\omega^2$;
    the attraction below the ionic plasma frequency, and Bohm–Staver sound. These are Bardeen and Pines' results.
  - Phonons are nuclei interacting through $W^{\mathrm e}$ (Baym).
  - The full Hedin–Baym set, read as a map: the clickable Table I.
- **Figures.**
  - *port* `viz-a-hedin-baym-fd` → `fd`: Feynman diagrams of the lattice inside $W$.
  - *port* `viz-a-hedin-baym-polar` → `polar`.
  - *port* `viz-a-hedin-baym-jellium` → `jellium`.
  - *port* `viz-a-hedin-baym-map` → `map`.
  - Optionally, *new*: the geometric picture of $W^{\mathrm{ph}}$ as two electrons coupled by a lattice wave.
- **Exports.**
  - `a5-eq-dysonG`;
  - `a5-eq-W`: $W$ and $\varepsilon$;
  - `a5-eq-split`: $W=W^{\mathrm e}+W^{\mathrm{ph}}$;
  - `a5-eq-Wph`: $W^{\mathrm{ph}}$ as $gFg$;
  - `a5-eq-F`: definition of $F$, with $F^0$;
  - `a5-eq-hedin`: the Hedin–Baym set.
- **May link.** `a2-eq-g`, `a3-eq-scf`, `a3-eq-eps`, `a4-eq-frohlich`.
- **Care.**
  - $F$ in 1/energy, as in NOTATION.
  - The review's $D$ for the phonon propagator and $P$ for the polarization belong in footnotes only.
  - The ionic plasma frequency is $\omega^{\mathrm{ion}}_{\mathrm p}$ (home A.5, `gl-a5-wp`); the electronic one is $\omega_{\mathrm p}$.

## A.6 · `A06-phonon-se.html` · id `a-phonon-se` · "The phonon self-energy: shifts, widths and the Kohn anomaly"
- **Source.** v1 `A06-phonon-se.html`. Notes G4 (§V.A–B) and G5 (§VII). Giustino Eqs. (123)–(144) and (181)–(186).
- **Cover.**
  - Phonons as poles of $F$.
  - The Dyson equation $F=F^0+F^0\Pi^{\mathrm{ph}}F$. Derive, to first order:
    - the shift $\Delta\omega\simeq\mathrm{Re}\,\Pi^{\mathrm{ph}}/\hbar$;
    - the width $\gamma\simeq-\mathrm{Im}\,\Pi^{\mathrm{ph}}/\hbar$, stating the convention.
  - The adiabatic split: $\Pi^{\mathrm{ph,A}}$ is what DFPT computes.
  - Inside $\Pi^{\mathrm{ph,NA}}$: one bare and one screened vertex. Explain the double-counting logic with diagrams.
  - The working formula with occupations, and the static subtraction.
  - A doped 2D semiconductor and its nonadiabatic window, $\hbar\omega$ versus $\hbar v_{\mathrm F}q$.
  - Where nonadiabatic effects have been seen: graphene, MgB₂, doped diamond.
  - Jellium: adiabatic phonons, Bohm–Staver sound and the **Kohn anomaly**. Give the geometric reason: $q=2k_{\mathrm F}$
    spans the Fermi surface, with the log singularity of the Lindhard function.
- **Figures.**
  - *port* `viz-a-phonon-se-fd-dyson` → `dyson`.
  - *port* `viz-a-phonon-se-fd-bubble` → `bubble`.
  - *port* `viz-a-phonon-se-2d` → `window`.
  - *port* `viz-a-phonon-se-jellium` → `kohn`.
  - *new* `fermi`: a Fermi circle and the same circle shifted by $\bq$. Their overlap and tangency at $q=2k_{\mathrm F}$
    explain the anomaly; the intraband transitions $\bk\to\bk+\bq$ are shown as arrows.
- **Exports.**
  - `a6-eq-dysonF`;
  - `a6-eq-shift`: shift and width;
  - `a6-eq-pi`: $\Pi^{\mathrm{ph}}$ split into A and NA;
  - `a6-eq-pina`: the working formula;
  - `a6-eq-kohn`.
- **May link.** `a5-eq-F`, `a5-eq-Wph`, `a2-eq-g`, `a3-eq-chi0`.
- **Care.** Units: $\Pi^{\mathrm{ph}}$ in energy. Retarded with $+\ii\eta$. The review's $\Pi$ goes in a footnote only.

## A.7 · `A07-electron-se.html` · id `a-electron-se` · "The electron self-energy: Fan–Migdal, Debye–Waller and kinks"
- **Source.** v1 `A07-electron-se.html`. Note G4 (§V.C–D). Giustino Eqs. (145)–(174).
- **Cover.**
  - Clamped nuclei as the reference.
  - The phonon part of $\Sigma=\ii GW$: derive Fan–Migdal from $W^{\mathrm{ph}}$ at finite temperature, retarded. Show:
    - emission and absorption terms with $n_{\mathrm B}$ and $f$;
    - their physical reading as processes.
  - Debye–Waller: the crystal potential blurred by vibrations.
  - The spectral function $\mathcal A_{n\bk}(\omega)$, $Z^{\mathrm{qp}}$, and the mass enhancement $1+\lambda_{\mathrm{ep}}$.
  - The ARPES kink at $\hbar\omega_{\mathrm{ph}}$.
  - The Eliashberg function $\alpha^2F(\omega)$ and $\lambda_{\mathrm{ep}}=2\int\dd\omega\,\alpha^2F/\omega$.
  - Beyond Dyson: the cumulant and polaron satellites.
- **Figures.**
  - *new* `processes`: a band diagram. An electron at $\bk$ emits or absorbs a phonon $\bq\nu$ and lands at $\bk+\bq$.
    Show the energy windows, Pauli blocking by the Fermi sea, and the factors $n_{\mathrm B}+1-f$ and $n_{\mathrm B}+f$.
  - *port* `viz-a-electron-se-kink` → `kink`.
  - *port* `viz-a-electron-se-cumulant` → `cumulant`.
- **Exports.**
  - `a7-eq-fm`: Fan–Migdal;
  - `a7-eq-dw`: Debye–Waller self-energy;
  - `a7-eq-spectral`: $\mathcal A$;
  - `a7-eq-Z`: $Z^{\mathrm{qp}}$;
  - `a7-eq-a2f`;
  - `a7-eq-lambda`.
- **May link.** `a5-eq-Wph`, `a5-eq-dysonG`, `a2-eq-g`, `a3-eq-gdw`.
- **Care.** Give the sign conventions of $\mathrm{Im}\,\Sigma$. Avoid $A$ for the spectral function.

## A.8 · `A08-temperature.html` · id `a-temperature` · "Band gaps versus temperature, lifetimes and the Fröhlich catastrophe"
- **Source.** v1 `A08-temperature.html`. Note G5 (§VIII). Giustino §VIII, Eqs. (187)–(199).
- **Cover.**
  - Quasiparticle energies from Fan–Migdal plus Debye–Waller, on shell versus off shell.
  - Allen–Heine: derive, with steps, how the rigid-ion approximation and the acoustic sum rule write Debye–Waller with
    first-order matrix elements.
  - The Varshni curve from Bose factors: the zero-point plateau, the linear high-$T$ regime, and the ZPR from the
    high-$T$ asymptote. Diamond and Si numbers.
  - Lifetimes from $\mathrm{Im}\,\Sigma$, as a Fermi golden rule.
  - The LO emission threshold in polar semiconductors.
  - Adiabatic versus nonadiabatic: why the adiabatic ZPR diverges in polar materials, the Fröhlich catastrophe. Show
    the $\int\dd^3q/q^2\times1/(\text{energy denominator}\to0)$ argument.
  - Nonperturbative adiabatic methods: frozen phonons, the special displacement method, molecular dynamics.
- **Figures.**
  - *port* `viz-a-temperature-varshni` → `varshni`.
  - *port* `viz-a-temperature-lo` → `lo`.
  - *port* `viz-a-temperature-frohlich` → `catastrophe`.
  - Optionally, *new* `zg`: a supercell with a special displacement pattern, atoms moved by $\pm$ zero-point
    amplitudes per mode. This is geometric.
- **Exports.**
  - `a8-eq-qp`: on-shell energy shift;
  - `a8-eq-ah`: Allen–Heine;
  - `a8-eq-zpr`;
  - `a8-eq-varshni`;
  - `a8-eq-tau`: lifetime;
  - `a8-eq-catastrophe`.
- **May link.** `a7-eq-fm`, `a7-eq-dw`, `a3-eq-gdw`, `a4-eq-frohlich`, `a1-eq-asr`.
- **Care.** Factors of 2 in widths and lifetimes: state $1/\tau=2\,\mathrm{Im}\Sigma/\hbar$ with sign.

## A.9 · `A09-applications.html` · id `a-applications` · "Light absorption, transport and superconductivity"
- **Source.** v1 `A09-applications.html`. Notes G5 and G6. Giustino §§IX–XI.
- **Cover.**
  - Phonon-assisted (indirect) absorption from second-order perturbation theory. Light supplies the energy, the phonon
    the momentum. Derive the structure of the rate and the role of $n_{\mathrm B}$ and $n_{\mathrm B}+1$.
  - Temperature: two slopes, and a moving edge.
  - Transport: from lifetimes to mobility through the Boltzmann equation, in a sketch. State SERTA versus the iterative
    solution. Give Si and GaAs numbers.
  - Superconductivity: from $\alpha^2F$ to $T_{\mathrm c}$. Present McMillan–Allen–Dynes with $\mu^*$, $\omega_{\log}$,
    $\lambda_{\mathrm{ep}}$, and the Eliashberg idea.
  - Is the DFT vertex too weak? GW-corrected couplings.
  - What Part A assumed: light external, the Coulomb interaction instantaneous. This sets up the bridge.
- **Figures.**
  - *new* `indirect`: a band diagram of an indirect-gap semiconductor. The photon arrow is vertical and the phonon arrow
    horizontal, with $\bq=\bk_{\mathrm c}-\bk_{\mathrm v}$. Show both orders of the two steps and the virtual
    intermediate state.
  - *port* `viz-a-applications-edge` → `edge`.
  - *port* `viz-a-applications-tc` → `tc`.
- **Exports.**
  - `a9-eq-indirect`;
  - `a9-eq-mob`;
  - `a9-eq-tc`.
- **May link.** `a7-eq-a2f`, `a7-eq-lambda`, `a8-eq-tau`.
- **Care.** $\mu_{\mathrm{mob}}$, not $\mu$. Write $\varepsilon_2$ and $\alpha_{\mathrm{abs}}$ for absorption.

## X.1 · `X01-light-cone.html` · id `x-light-cone` · "From the Coulomb line to the full field"
- **Source.** v1 `X01-light-cone.html`. Note SYNTHESIS. EPP §§4–6 overview.
- **Cover.**
  - The bridge header (`part-x`) first.
  - Every interaction line in Part A is a Coulomb line.
  - Two meanings of "retardation": the frequency dependence of screening, and the finite speed of light.
  - What EPP puts in its place: the photon propagator $D$ and the photon self-energy $\Pi$. This section is their home;
    define them in words and with the formula of NOTATION §12.4.
  - Two Dyson equations with one structure: $W=v+v\chi^{\mathrm{irr}}W$ and $D=D^0+D^0\Pi D$.
  - How big the light cone is: $\omega_{\mathrm{LO}}/c\approx10^{-4}$ Å⁻¹ against a Brillouin zone of about 1 Å⁻¹.
  - The order of limits, $q\to0$ against $c\to\infty$.
- **Figures.**
  - *port* `viz-x-light-cone-dyson` → `dyson`.
  - *port* `viz-x-light-cone-cone` → `cone`. Make it geometric and to scale, on a log axis, with the zone, the light
    cone and the phonon branches.
- **Exports.**
  - `x1-eq-dysonW`;
  - `x1-eq-dysonD`;
  - `x1-eq-D`: definition of $D$.
- **May link.** `a5-eq-W`, `a5-eq-split`, `a4-eq-glr`, `a4-eq-lot`.

## X.2 · `X02-dictionary.html` · id `x-dictionary` · "Dictionary between the two formalisms"
- **Source.** v1 `X02-dictionary.html`. Note SYNTHESIS.
- **Cover.**
  - In v2 the notation is already unified, so this section maps **concepts and objects**. Three substitutions in v2
    notation:
    - Coulomb $v$ → photon propagator $D$;
    - density response $\chi$ → photon self-energy $\Pi$;
    - charge vertex → current vertex, with $\Pi^{\mathrm e}_{\rho\rho}=e^2\chi^{\mathrm{irr}}$.
  - A concept map.
  - A units dictionary, Gaussian against SI. Write one Fröhlich vertex both ways.
  - A diagram-glyph dictionary.
  - What EPP does not (yet) contain: Debye–Waller, the nonadiabatic electron self-energy, temperature effects,
    anharmonicity, beyond-RPA vertices.
  - Here, and only here, a table of Giustino, EPP and v2 symbols is allowed, inside `<!-- srcmap -->…<!-- /srcmap -->`.
    Derive it from NOTATION §9 and §12 and check every entry against the sources.
- **Figures.**
  - *port* `viz-x-dictionary-units` → `units`.
  - *port* `viz-x-dictionary-glyphs` → `glyphs`.
- **Exports.** `x2-eq-subst`: the three substitutions.
- **May link.** Any exports.

## B.1 · `B01-foundations.html` · id `b-foundations` · "The EPP starting point: lattice, bare vertex and polar coupling"
- **Source.** v1 `B01-foundations.html`. Notes E1–E3. EPP §§1–3. v1 `Z03-corrections.html` §1.
- **Cover.**
  - The Part B header (`part-b`) first.
  - Recap, in a few paragraphs with links, rather than re-deriving:
    - EPP §1's three-step diagonalization is A.1's, so link (A1.18) and `figref` to `viz-a1-diatomic`;
    - EPP's normal coordinates and phonons are A.1's.
  - The bare vertex (EPP §2): the unscreened electron–ion vertex.
  - Two starting points and the double-counting rule: start from bare vertices and let diagrams screen, or start from
    DFPT's screened vertex and never screen twice.
  - Polar coupling from Born charges (EPP §3) in v2 notation, consistent with A.4.
  - Why the basis phase is needed: derive the cell-choice invariance. Only $e_{\kappa\beta,\nu}(\bq)\,\ee^{-\ii(\bq+\bG)\cdot\btau_\kappa}$
    is unchanged when an atom is reassigned to a neighbouring cell.
  - The long-wave expansion of the vertex.
  - What carries forward into B.2–B.10.
- **Figures.**
  - *drop* `viz-b-foundations-modes`, which duplicates A.1; link it instead.
  - *port* `viz-b-foundations-kernel` → `kernel`: the anisotropic long-range kernel.
  - *port* `viz-b-foundations-cell` → `cell`: printed versus corrected vertex under a change of cell. Geometric: show
    the two cell choices.
- **Exports.**
  - `b1-eq-bare`;
  - `b1-eq-polar`: the corrected polar vertex;
  - `b1-eq-cell`: the invariant product.
- **May link.** A.1 exports, `a2-eq-g`, `a4-eq-glr`, `a4-eq-born`.
- **Care.**
  - EPP §3 is in SI units and §4 on is Gaussian.
  - From §4 on, EPP's $e$ is the signed electron charge and its $Z^*$ carries the charge.
  - State the correction once, in a footnote.

## B.2 · `B02-lagrangian.html` · id `b-lagrangian` · "One Lagrangian for electrons, ions and the field"
- **Source.** v1 `B02-lagrangian.html`. Note E4. EPP §4.1–4.2 start, Eqs. (~200)–(~250).
- **Cover.**
  - Gauge freedom, in a background box.
  - The temporal gauge $A_0=0$: what it means, and why it is allowed. Then $\bE=-\frac1c\partial_t\bA$ and
    $\bB=\nabla\times\bA$.
  - The Lagrangian term by term, each explained:
    - the electron field with minimal coupling;
    - the harmonic ions;
    - the ions coupled through their Born dipoles, $-\sum e\mathbf Z^*\bu\cdot\bE$ at the atom positions (basis
      phase rule);
    - the field energy $(E^2-B^2)/8\pi$.
  - Why the ions couple through their Born dipoles.
  - Why there is no electron–phonon term: it is carried by the field.
  - Euler–Lagrange equations, derived:
    - electrons, the Lorentz force;
    - ions, a lattice driven by the local field.
- **Figures.**
  - *port* `viz-b-lagrangian-map` → `map`: who talks to whom.
  - *port* `viz-b-lagrangian-drive` → `drive`: driven ions and electrons.
- **Exports.**
  - `b2-eq-L`: the Lagrangian;
  - `b2-eq-gauge`;
  - `b2-eq-lorentz`;
  - `b2-eq-ions`;
  - `b2-eq-P`: the Born-dipole polarization, with the atom positions.
- **May link.** `a4-eq-born`, `a4-eq-P`, `b1-eq-polar`.
- **Care.** Gaussian units throughout, with SI footnotes. The basis-phase rule: dipoles at $\btau_{\ell\kappa}$, not at $\bR_\ell$.

## B.3 · `B03-dielectric.html` · id `b-dielectric` · "The polar crystal as a dielectric"
- **Source.** v1 `B03-dielectric.html`. Notes E4 and E5 start. EPP §4.2.
- **Cover.**
  - Two masses, one spring and one field: the driven diatomic crystal of B.2 at $\bq\to0$.
  - Derive the Lorentz oscillator with the ionic plasma frequency $\omega_{\mathrm p}^{\mathrm{ion}}$.
  - Adding the electrons through $\varepsilon_\infty$ gives
    $\varepsilon(\omega)=\varepsilon_\infty+(\varepsilon_{\mathrm s}-\varepsilon_\infty)\omega_{\mathrm{TO}}^2/(\omega_{\mathrm{TO}}^2-\omega^2-\ii\gamma\omega)$.
  - Zeros are LO and poles are TO: why, in terms of fields.
  - The Lyddane–Sachs–Teller relation, derived, and what is model independent.
  - The Reststrahlen band, with reflectance and the loss function.
  - The lattice inside $W=v/\varepsilon(\omega)$: its lattice part is the Fröhlich interaction (link A.5).
- **Figures.**
  - *port* `viz-b-dielectric-lt` → `lt`: why the LO mode carries a field and the TO mode does not. Geometric:
    polarization charges on slabs.
  - *port* `viz-b-dielectric-lab` → `lab`: the Reststrahlen lab.
- **Exports.**
  - `b3-eq-eom`;
  - `b3-eq-eps`;
  - `b3-eq-lst`;
  - `b3-eq-W`.
- **May link.** `b2-eq-ions`, `a4-eq-lot`, `a5-eq-Wph`.
- **Care.** Local field versus macroscopic field, stated. Give GaAs numbers.

## B.4 · `B04-fields.html` · id `b-fields` · "Quantized fields: transverse photons and the Coulomb interaction"
- **Source.** v1 `B04-fields.html`. Note E5. EPP §§4.3–4.4.
- **Cover.**
  - Maxwell's equations from the field Lagrangian in the temporal gauge, derived.
  - The Legendre transform to the Hamiltonian, with the canonical momentum of $\bA$ ($\propto-\bE$). It gives four
    minimal couplings:
    - the paramagnetic $\bj\cdot\bA$;
    - the diamagnetic $A^2$;
    - the ion–field coupling;
    - field energy.
  - Quantizing the field in the polarization triad $\be_1,\be_2,\be_3=\hat\bq$.
  - The transverse modes are oscillators, with the field expansion of NOTATION §12.4 and photon operators $b$.
  - The longitudinal mode has no restoring force; it drifts. Gauss's law as a constraint turns it into the
    instantaneous Coulomb interaction. Derive this.
  - The final split $H=H_0+V'$.
- **Figures.**
  - *port* `viz-b-fields-helmholtz` → `helmholtz`: one Fourier mode and the triad. Geometric, in 3D projection.
  - *port* `viz-b-fields-drift` → `drift`: oscillator against drift in phase space.
- **Exports.**
  - `b4-eq-maxwell`;
  - `b4-eq-H`: the Hamiltonian;
  - `b4-eq-A`: the field expansion;
  - `b4-eq-coulomb`: the longitudinal photon as Coulomb;
  - `b4-eq-H0`: the split.
- **May link.** `b2-eq-L`, `b2-eq-gauge`.
- **Care.**
  - Gaussian normalization $(2\pi\hbar c^2/\OVSC\omega_q)^{1/2}$.
  - Photon operators are $b$, not $a$.
  - Use the polarization-vector convention needed for a Hermitian $\bA$: state $\be^*_\lambda(-\bq)=\be_\lambda(\bq)$, or the convention used.

## B.5 · `B05-photon-se.html` · id `b-photon-se` · "Free propagators and the photon self-energy"
- **Source.** v1 `B05-photon-se.html`. Notes E6 and E7. EPP §5, Eqs. (~246)–(~320), and Eqs. (379)–(385).
  v1 `Z03-corrections.html` §1.
- **Cover.**
  - Temperature and the "$n_{\mathrm B}=0$" rule: why the bosonic propagators used here do not depend on $T$.
  - Derive the free propagators, giving each the home of NOTATION §12.4:
    - the free phonon line $F^0$;
    - the free transverse photon line $D^{0\mathrm T}$;
    - the longitudinal photon line, which vanishes while Coulomb does not.
  - Dressing the transverse photon: $D^{\mathrm e}=D^0+D^0\Pi^{\mathrm e}D^{\mathrm e}$.
  - Gauge invariance fixes the diamagnetic term: $\Pi^{\mathrm{dia}}=-\Pi^{\mathrm{para}}(\omega=0)$ for an insulator.
  - The electronic current vertex $j_\alpha^{\,n\bk,m\bk+\bq}$.
  - The ionic self-energy, exactly, with the ionic current vertex $j_{\alpha\nu}(\bq)$ **with the basis phase**.
  - Completeness: the diamagnetic term is minus the static bubble.
  - The infrared dielectric tensor from $\Pi$; connect it to B.3.
- **Figures.**
  - *port* `viz-b-photon-se-electrons` → `electrons`: truncated band sums make a fake Drude metal.
  - *port* `viz-b-photon-se-sumrule` → `sumrule`.
  - *port* `viz-b-photon-se-gauge` → `gauge`.
- **Exports.**
  - `b5-eq-F0`;
  - `b5-eq-D0T`;
  - `b5-eq-D0L`;
  - `b5-eq-dysonDe`;
  - `b5-eq-Pie`: the electronic current–current self-energy;
  - `b5-eq-jel`;
  - `b5-eq-jion`;
  - `b5-eq-Piion`;
  - `b5-eq-dia`;
  - `b5-eq-epsir`.
- **May link.** `b4-eq-A`, `b4-eq-coulomb`, `a5-eq-F`, `b3-eq-eps`.
- **Care.**
  - The basis-phase correction goes in a footnote.
  - EPP's signed $e$ and its charge-carrying $Z^*$.
  - The factor $\omega^2/c$ and the units of $\Pi$ must be stated.

## B.6 · `B06-propagator.html` · id `b-propagator` · "The total propagator: connection formulas, plasmons and polaritons"
- **Source.** v1 `B06-propagator.html`. Note E8. EPP §6.0.1–6.0.4 and §6.1 start. v1 `Z03-corrections.html` §3.
- **Cover.**
  - One interaction, two kinds of lines.
  - Induced charges and currents from the four correlators $\Pi_{\rho\rho}$, $\Pi_{\rho\beta}$, $\Pi_{\alpha\rho}$,
    $\Pi_{\alpha\beta}$; define them here.
  - Charge conservation: derive the connection formulas step by step from the continuity equation **with the factor
    $e$**. Use charge density and charge current consistently, and give the correction in a footnote.
  - From four sectors to one tensor: the sectors $D_{\rho\rho}$, $D_{\rho\beta}$, $D_{\alpha\rho}$, $D_{\alpha\beta}$
    and the L and T projections.
  - What $D^0$ is: the vacuum Maxwell propagator.
  - One Dyson equation, $D=D^0+D^0\Pi D$.
  - An isotropic medium:
    - $\varepsilon^{\mathrm L}(\omega)=0$ gives plasmons and LO phonons;
    - $\varepsilon^{\mathrm T}(\omega)=c^2q^2/\omega^2$ gives polaritons.
- **Figures.**
  - *port* `viz-b-propagator-map` → `map`.
  - *new* `polariton`: an interactive dispersion $\omega(q)$ with the light lines $cq$ and $cq/\sqrt{\varepsilon_\infty}$,
    $\omega_{\mathrm{TO}}$, $\omega_{\mathrm{LO}}$, the upper and lower polariton branches and the LO branch. A log-$q$
    option shows that the polariton region is $q\sim\omega_{\mathrm{TO}}/c$. Use GaAs.
  - Optionally, *new* `sectors`: the 4×4 → 3×3 structure, drawn.
- **Exports.**
  - `b6-eq-induced`;
  - `b6-eq-continuity`;
  - `b6-eq-connection`;
  - `b6-eq-sectors`;
  - `b6-eq-dysonD`;
  - `b6-eq-modesL`: longitudinal modes;
  - `b6-eq-modesT`: polaritons.
- **May link.** `b5-eq-D0T`, `b5-eq-D0L`, `b5-eq-Pie`, `b5-eq-Piion`, `b3-eq-eps`, `x1-eq-dysonD`.
- **Care.** Charge and number densities must not mix. State the units of each sector.

## B.7 · `B07-phonon-se.html` · id `b-phonon-se` · "The phonon self-energy from photon exchange; phonon polaritons"
- **Source.** v1 `B07-phonon-se.html`. Note E9. EPP §6.1, Eqs. (386)–(410′). v1 `Z03-corrections.html` §2.
- **Cover.**
  - Step 1: resum the electrons into the photon line ($D^{\mathrm e}$).
  - Step 2: open the ionic bubble. Derive $\Pi^{\mathrm{ion}}$ written with $F^0$, and the dynamic vertex
    $J_{\alpha\nu}(\bq,\omega)$.
  - Step 3: a Dyson equation for phonons, with $\Pi^{\mathrm{ph}}_{\nu\nu'}=-\sum J^*_{\alpha\nu}D^{\mathrm e}_{\alpha\beta}J_{\beta\nu'}$,
    derived from the series.
  - EPP §7.1: switch the electrons off and set $q=0$.
  - What EPP Eq. (410′) is: the upper phonon polariton at $q\to0$, whose frequency equals $\omega_{\mathrm{LO}}$ in a
    cubic crystal. It is a transverse mode, not the LO–TO splitting. Show both, and why the numbers agree.
  - The order of limits.
- **Figures.**
  - *port* `viz-b-phonon-se-limits` → `limits`: the order of limits.
  - *new* `exchange`: a diagram of a phonon emitting and reabsorbing a photon dressed by electrons, with the series
    $F=F^0+F^0\Pi^{\mathrm{ph}}F^0+\dots$ drawn.
  - Optionally, *new*: the frequency of the coupled mode against $q$ near the light line, solved from the Dyson
    equation. It may reuse the B.6 polariton physics with the phonon viewpoint.
- **Exports.**
  - `b7-eq-Piion`;
  - `b7-eq-J`;
  - `b7-eq-piph`;
  - `b7-eq-dysonF`;
  - `b7-eq-410`: the corrected reading.
- **May link.** `b5-eq-jion`, `b5-eq-F0`, `b5-eq-dysonDe`, `b6-eq-modesT`, `a6-eq-dysonF`.
- **Care.** This is the second correction; give it in a footnote and in the text.

## B.8 · `B08-frohlich.html` · id `b-frohlich` · "The electrodynamic Fröhlich vertex"
- **Source.** v1 `B08-frohlich.html`. Note E9. EPP §7, pp. 86–99. v1 `Z03-corrections.html`, §3 and "Where the corrections meet".
- **Cover.**
  - A vertex made of a photon: the electron couples to a phonon by exchanging a photon with the ionic current.
  - Opening the electron bubble.
  - Two vertices and one bubble: which vertices go into the self-energy (counting).
  - From currents to charges: the continuity identity, with the factor $e$ restored.
  - The electrostatic limit gives Verdi–Giustino with $\varepsilon^{\mathrm e}(\bq,\omega)$, screened at the phonon
    frequency, not at zero frequency.
  - Across the light line: LO is Coulomb, while TO needs photons.
- **Figures.**
  - *port* `viz-b-frohlich-pairing` → `pairing`.
  - *port* `viz-b-frohlich-dynamic` → `dynamic`.
  - *port* `viz-b-frohlich-lightline` → `lightline`.
- **Exports.**
  - `b8-eq-vertex`: the electrodynamic vertex;
  - `b8-eq-identity`: the continuity identity;
  - `b8-eq-static`: the electrostatic limit;
  - `b8-eq-ratio`: LO and TO coupling across the light line.
- **May link.** `b7-eq-J`, `b5-eq-jel`, `b6-eq-connection`, `a4-eq-glr`, `a4-eq-frohlich`.
- **Care.** The factor-$e$ and basis-phase corrections meet here. Use the merged vertex of Z03.

## B.9 · `B09-electrostatic.html` · id `b-electrostatic` · "The electrostatic limit: recovering Part A"
- **Source.** v1 `B09-electrostatic.html`. Note E9. EPP §7.2 (or the relevant pages). Z03 E.454a.
- **Cover.**
  - Four sectors, one survivor: for $c\to\infty$ only $D_{\rho\rho}$ survives.
  - The survivor is $W^{\mathrm e}=D^{\mathrm e}_{\rho\rho}$.
  - Phonons dressed by $W^{\mathrm e}$ give the Gonze–Lee nonanalytic term: derive $\Phi^{\mathrm{dd}}$ of A.4 from the
    phonon self-energy.
  - Keeping $\omega$: LO-phonon–plasmon modes in doped polar semiconductors, from
    $\varepsilon(\omega)=\varepsilon_\infty-\omega_{\mathrm p}^2\varepsilon_\infty/\omega^2+\dots=0$.
  - What Part A gets exactly, and what it misses.
- **Figures.**
  - *port* `viz-b-electrostatic-sectors` → `sectors`.
  - *port* `viz-b-electrostatic-plasmon` → `plasmon`: the coupled modes against doping, for GaAs.
- **Exports.**
  - `b9-eq-limit`;
  - `b9-eq-We`;
  - `b9-eq-gl`: the Gonze–Lee term;
  - `b9-eq-lopl`.
- **May link.** `b6-eq-sectors`, `b6-eq-dysonD`, `b7-eq-piph`, `a4-eq-dd`, `a4-eq-lot`, `a5-eq-split`.

## B.10 · `B10-2d.html` · id `b-2d` · "Two-dimensional crystals: no LO–TO splitting"
- **Source.** v1 `B10-2d.html`. Notes E7 and E9 (2D replacements). The EPP 2D remarks; the literature cited in v1
  (Sohier, Keldysh, Cudazzo).
- **Cover.**
  - A sheet of point dipoles: the 2D Fourier transform of Coulomb, $v^{\mathrm{2D}}(Q_\parallel)=2\pi e^2/Q_\parallel$,
    derived.
  - The electrostatic limit: a strictly 2D screened interaction, $\varepsilon^{\mathrm{2D}}=1+r_0Q_\parallel$, with
    $r_0=2\pi\chi^{\mathrm{2D}}$ derived from the sheet polarizability.
  - Why LO meets TO at Γ: $\omega_{\mathrm{LO}}^2(Q)=\omega_{\mathrm{TO}}^2+\dfrac{2\pi e^2Z^{*2}Q}{\OS M_{\mathrm r}(1+r_0Q)}$
    and its slope. Derive it.
  - The same $r_0$ sets the exciton binding, through the Keldysh potential.
  - Doped sheets: plasmon–phonon modes, with the plasmon $\propto\sqrt{Q}$.
  - Finite $c$: when the sheet radiates; the light cone in 2D.
  - Give hBN and MoS₂ numbers.
- **Figures.**
  - *port* `viz-b-2d-field` → `field`: the depolarization field of an LO wave, bulk versus sheet. Geometric field lines.
  - *port* `viz-b-2d-disp` → `disp`.
  - Optionally, *new* `keldysh`: the Keldysh against Coulomb potential, with the $r_0$ crossover.
- **Exports.**
  - `b10-eq-v2d`;
  - `b10-eq-eps2d`;
  - `b10-eq-lo`;
  - `b10-eq-keldysh`.
- **May link.** `a4-eq-lot`, `b9-eq-We`, `b6-eq-modesT`.
- **Care.** Per-area quantities ($\OS$) and the $\parallel$ notation of NOTATION §3.
