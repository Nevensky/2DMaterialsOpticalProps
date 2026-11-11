# A.7 figure specifications

Section file: `v2/sections/A07-electron-se.html` (id `a-electron-se`, prefix `a7`). Four figures, in reading order:

| # | File | Id | Kind | v1 source |
|---|---|---|---|---|
| A7.1 | `A07/fig-processes.html` | `viz-a7-processes` | new | none |
| A7.2 | `A07/fig-diagrams.html` | `viz-a7-diagrams` | new | v1 inline diagrams `fd-a-electron-se-*` (not a figure in v1) |
| A7.3 | `A07/fig-kink.html` | `viz-a7-kink` | port | `viz-a-electron-se-kink` |
| A7.4 | `A07/fig-cumulant.html` | `viz-a7-cumulant` | port | `viz-a-electron-se-cumulant` |

v1 file for all ports: `doc/epp-primer/sections/A07-electron-se.html`. The kink figure is its lines 73–79 with the script
at lines 163–258; the cumulant figure is lines 106–112 with the script at lines 260–292. The v1 Feynman diagrams are
lines 21–41 with their `EPP.fd.draw` specs at lines 123–161.

General rules for all four (BUILDER_GUIDE §5): v2 notation only (no source symbols in labels, readouts, legends,
aria-labels); canvas text through `EPP.drawRich` (`E_{F}`, `ħω_0`, bold vectors with `{bold:true}`); colours read at
draw time with `EPP.color(...)`; fixed axis limits; controls below plots; 390 px without horizontal overflow; keyboard
access and `aria-label`s; no NaN, no console errors; caption of about 200–260 words in four parts (What to notice,
Marks, Try, Model) with `a.eqref` links to the section's equations. Energies of electron states are written $E$ (an
energy variable, glossary `gl-a7-Evar`), with $E=\hbar\omega$ where a self-energy is evaluated. Colours: electrons
`'electron'`, phonons `'phonon'`, secondary curves `'c4'`, `'c5'`, guides `'muted'`.

Notation used by the section (all defined there): $E_{n\bk}$ clamped-nuclei band energy; $E_{\mathrm F}$ Fermi energy;
$\mu$ chemical potential; $f(E)$ Fermi–Dirac function; $n_{\mathrm B}(\omega)$ Bose–Einstein function;
$\omega_0$ Einstein phonon frequency; $\Sigma^{\mathrm{FM}}$ Fan–Migdal self-energy; $\mathcal A$ spectral function
(calligraphic; canvas: the Unicode letter 𝒜); $Z^{\mathrm{qp}}$ quasiparticle weight; $E^{\mathrm{qp}}$ quasiparticle
energy; $\lambda_{\mathrm{ep}}$ coupling strength; $\lambda_{n\bk}$ mass-enhancement parameter of one state;
$m^*$, $m^{\mathrm{qp}}$ effective and quasiparticle mass; $k_{\mathrm F}$ Fermi wavevector; $N_{\mathrm F}$ density
of states per spin and per cell at the Fermi level; $E_{\mathrm c}$ cutoff; $\eta$ the retarded infinitesimal, kept
finite as a broadening $\hbar\eta$; $S$, $E_0$, $g_0$ the single-level model; $s$ the replica index.

---

## A7.1 `fig-processes.html` · `viz-a7-processes` · new

**Title.** "Emission and absorption of a phonon, with Pauli blocking" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. An electron in state $\bk$ at energy $E$ lands at $\bk+\bq$ on the band by emitting a phonon $-\bq\nu$ (final
   energy $E-\hbar\omega_0$) or absorbing a phonon $\bq\nu$ (final energy $E+\hbar\omega_0$); the horizontal component
   of the arrow is the crystal momentum $\hbar\bq$ given to or taken from the lattice.
2. The two weights of (A7.8): $n_{\mathrm B}+1-f$ for the final state one phonon energy below, $n_{\mathrm B}+f$ for the
   one above. Each is "out" plus "in": $(1-f)(n_{\mathrm B}+1)+f\,n_{\mathrm B}$ and $(1-f)\,n_{\mathrm B}+f(n_{\mathrm B}+1)$.
3. At $T=0$ there is a window of width $2\hbar\omega_0$ around $E_{\mathrm F}$ in which both weights vanish: no real
   process, infinite phonon-limited lifetime. Below $E_{\mathrm F}$ the second term is the decay of a hole (an
   electron from $E+\hbar\omega_0$ drops in and emits).
4. Temperature adds $n_{\mathrm B}$ to both weights and smears the window.

**Include position.** End of h4 `#a7-sea` ("A Fermi sea and a warm lattice"), before h3 `#a7-reference`. The sentences
just before it: "The derivation below produces exactly these weights. At zero temperature they leave a window of width
$2\hbar\omega_{\bq\nu}$ around the chemical potential in which an electron, or a hole, can neither emit nor absorb a real
phonon. [Figure A7.1] shows the processes on a band, for one dispersionless (Einstein) phonon of frequency $\omega_0$."

**Physics.** Eqs. (A7.8) `a7-eq-fm` and (A7.9) `a7-eq-im` for one band, one Einstein phonon $\hbar\omega_0$ and a
constant $|g|^2$. Model band: two-dimensional parabolic band $E_k=\hbar^2k^2/2m^*$ measured from its bottom, with
$m^*=0.5\,m_{\mathrm e}$ ($\hbar^2/2m_{\mathrm e}=3.80998$ eV Å²) and Fermi energy $E_{\mathrm F}=0.40$ eV above the
bottom, so $k_{\mathrm F}=0.229$ Å⁻¹; a 2D parabolic band has a constant density of states $N_{\mathrm F}$. Illustrative
parameters of a doped two-dimensional semiconductor; not a specific material. Measure energies from $E_{\mathrm F}$
($\mu=E_{\mathrm F}$); the band bottom is then at $E_{\mathrm b}=-0.40$ eV.
- Occupation $f(E)=1/(\ee^{E/\kB T}+1)$ (use `EPP.fermi`), $n_{\mathrm B}=1/(\ee^{\hbar\omega_0/\kB T}-1)$ (`EPP.bose`),
  $\kB=0.0861733$ meV/K; at $T=0$ use the step functions.
- Weights at probe energy $E$ (on shell, $\hbar\omega=E$):
  $w_{\mathrm{em}}(E)=\big[n_{\mathrm B}+1-f(E-\hbar\omega_0)\big]\,\theta(E-\hbar\omega_0-E_{\mathrm b})$ (a final state must
  exist above the band bottom); $w_{\mathrm{abs}}(E)=n_{\mathrm B}+f(E+\hbar\omega_0)$.
- With $|g|^2N_{\mathrm F}=\lambda_{\mathrm{ep}}\hbar\omega_0/2$ and $\lambda_{\mathrm{ep}}=0.5$ fixed (the convention of
  (A7.20)), $|\mathrm{Im}\,\Sigma^{\mathrm{FM}}(E)|=\frac{\pi\lambda_{\mathrm{ep}}\hbar\omega_0}{2}\,[w_{\mathrm{em}}+w_{\mathrm{abs}}]$
  and $\tau=\hbar/(2|\mathrm{Im}\,\Sigma|)$ with $\hbar=658.212$ meV fs (show "∞" when $|\mathrm{Im}\,\Sigma|=0$).
- Check values at the initial state ($E=+0.15$ eV, $\hbar\omega_0=60$ meV, $T=300$ K): $n_{\mathrm B}=0.109$,
  $f(E-\hbar\omega_0)=0.030$, $f(E+\hbar\omega_0)=3\times10^{-4}$, $w_{\mathrm{em}}=1.079$, $w_{\mathrm{abs}}=0.109$,
  $|\mathrm{Im}\,\Sigma|=56$ meV, $\tau=5.9$ fs.

**What to draw.** Two panels side by side on desktop, stacked on phones, sharing the vertical energy axis
$E-E_{\mathrm F}$ from −0.45 to +0.55 eV (fixed).
- (a) Band diagram, horizontal axis $k$ from −0.36 to +0.36 Å⁻¹ (label `k (Å^{-1})`). The parabola in `'electron'`;
  a thicker overlay of the same curve with opacity $f(E_k)$ shows the Fermi sea (fully opaque below $E_{\mathrm F}$ at
  $T=0$, smeared at finite $T$). Dashed `'muted'` line at $E_{\mathrm F}$, labelled `E_{F}`. A light `'phonon'`-tinted
  band between $E_{\mathrm F}-\hbar\omega_0$ and $E_{\mathrm F}+\hbar\omega_0$ across both panels, labelled
  "window at T = 0". The probed state: a filled dot on the branch $k\gt0$ at energy $E$, labelled `n` and bold `k`
  (rich text "n**k**"). Two horizontal dashed `'phonon'` guides at $E-\hbar\omega_0$ and $E+\hbar\omega_0$, labelled
  `E − ħω_0` and `E + ħω_0`. Final states: circles where each guide meets the band, at $\pm k'$ with
  $k'=\sqrt{2m^*(E\mp\hbar\omega_0-E_{\mathrm b})}/\hbar$ (none if below the band bottom); fill each circle with
  opacity equal to its occupation $f$. Arrows from the probed state to each final state, solid `'phonon'` for emission
  (downward, labelled "emit −**q**" on the forward one), dashed `'phonon'` for absorption (upward, "absorb **q**"),
  with line width and opacity proportional to the weight ($w_{\mathrm{em}}$ or $w_{\mathrm{abs}}$, clipped at 2) so
  that a blocked process fades out. On the forward emission arrow draw a small horizontal bracket labelled bold
  `q` = $k'-k$ to show the crystal momentum.
- (b) Weights against energy: horizontal axis "weight" from 0 to 4 (fixed). Curves $w_{\mathrm{em}}(E)$ (`'phonon'`
  solid, legend "emission: n_B + 1 − f"), $w_{\mathrm{abs}}(E)$ (`'phonon'` dashed, "absorption: n_B + f") and their
  sum (`'ink'`, "sum = |Im Σ| / (πλ_{ep}ħω_0/2)"). A horizontal `'muted'` line at the probe energy with dots on the
  three curves. The same window band and $E_{\mathrm F}$ line as in (a).
- `touch-action: none` only on panel (a), where the probe can be dragged along the band.

**Controls** (below the panels).
1. Probe energy $E-E_{\mathrm F}$: −0.38 to +0.50 eV, step 0.005 eV, default +0.15 eV (also set by dragging the dot
   in (a) along the band).
2. Phonon energy $\hbar\omega_0$: 30 to 150 meV, step 5 meV, default 60 meV.
3. Temperature $T$: 0 to 600 K, step 10 K, default 300 K.

**Readouts** (HTML with MathJax): $n_{\mathrm B}(\omega_0)$; $f(E-\hbar\omega_0)$; $f(E+\hbar\omega_0)$; emission weight
$n_{\mathrm B}+1-f$; absorption weight $n_{\mathrm B}+f$; $|\mathrm{Im}\,\Sigma^{\mathrm{FM}}|$ for
$\lambda_{\mathrm{ep}}=0.5$ (meV); $\tau=\hbar/2|\mathrm{Im}\,\Sigma|$ (fs).

**Initial state.** $E-E_{\mathrm F}=+0.15$ eV, $\hbar\omega_0=60$ meV, $T=300$ K: emission open, absorption weak, the
window visible in (b) as the region where the sum nearly vanishes.

**Caption key points.**
- *What to notice.* At low temperature both weights vanish within $\hbar\omega_0$ of the Fermi level, so states there
  cannot decay by phonons. Above the window an electron emits into empty states; below it the hole decays through the
  second term, as an electron from one phonon energy higher drops in.
- *Marks.* Band, Fermi sea opacity $f$, window band, guides at $E\pm\hbar\omega_0$, final-state circles filled by $f$,
  solid emission and dashed absorption arrows with width ∝ weight, the bracket $\bq$; in (b) the two weights and their
  sum.
- *Try.* (1) Set $T=0$ and drag the state into the window: both arrows fade, $|\mathrm{Im}\,\Sigma|=0$, $\tau=\infty$.
  (2) Move the state to $-0.2$ eV: the absorption-type term carries weight $f=1$, the hole decays. (3) Raise $T$ to
  600 K with $\hbar\omega_0=30$ meV: both weights grow by $n_{\mathrm B}\approx1.3$ and the window fills in.
- *Model.* (A7.8), (A7.9) with one Einstein phonon and constant $|g|^2$; 2D parabolic band, $m^*=0.5\,m_{\mathrm e}$,
  $E_{\mathrm F}=0.40$ eV; $|g|^2N_{\mathrm F}=\lambda_{\mathrm{ep}}\hbar\omega_0/2$ with $\lambda_{\mathrm{ep}}=0.5$;
  illustrative, not a specific material. Emission creates the phonon $-\bq\nu$, absorption removes $\bq\nu$, as for
  $g_{mn\nu}(\bk,\bq)$ of Section A.2.

---

## A7.2 `fig-diagrams.html` · `viz-a7-diagrams` · new

**Title.** "The electron–phonon self-energy as diagrams" `<span class="tag">diagram</span>`.

**What the reader must learn.**
1. The Dyson equation (A7.3) starts from the clamped-nuclei propagator $G^{\mathrm{cn}}$ and adds $\Sigma^{\mathrm{ep}}$.
2. $\Sigma^{\mathrm{ep}}=\Sigma^{\mathrm{FM}}+\Sigma^{\mathrm{DW}}+\Sigma^{\mathrm{dGW}}$ (A7.4): Fan–Migdal is an electron
   line with a phonon arc, a dynamic process; Debye–Waller is a phonon loop attached at one point, an equal-time
   (static) term; dGW is the difference of two electronic GW diagrams.
3. Vertex counting: the Fan–Migdal diagram has a screened vertex at both ends of the phonon line, while the phonon
   self-energy of Section A.6 (an electron bubble between two phonon lines) has one bare and one screened vertex.

**Include position.** In h3 `#a7-reference`, after its second paragraph and before the Dyson equation (A7.3). The
sentences just before it: "Letting the nuclei move changes three things. The electron can exchange phonons with the
lattice; the average potential of the vibrating nuclei differs from the potential of nuclei at rest; and the
propagator, density and screening of the electrons change slightly. [Figure A7.2] draws the three corresponding terms."

**Physics.** Eqs. (A7.3) `a7-eq-dyson`, (A7.4) `a7-eq-split`, (A7.5) `a7-eq-wph`; after the review's Fig. 2(a)–(e)
and Fig. 1(c). The vertex function $\Lambda$ is set to one (approximation 1 of `#a7-working`) and not drawn.

**What to draw.** Three rows of `EPP.fd.draw` diagrams (the helper exists in `v2/build/lib.js`; v1 lines 123–161 give
working item lists to start from), each row a `div.fd` with `span.op` signs and short `figcaption`s.
- (a) Dyson: thick double electron line labelled `G` = thin line labelled `G^{cn}` + thin `G^{cn}` – grey disk
  labelled `Σ^{ep}` – thick `G`.
- (b) $\Sigma^{\mathrm{ep}}$ = [Fan–Migdal: thick electron line `G` with a thick wavy phonon arc labelled `F`, grey disks
  (`'muted'` blobs) at both ends, caption "Σ^{FM}: two screened vertices"] + [Debye–Waller: a wavy phonon loop `F`
  attached at one point to the electron line through a hatched (or outlined) disk labelled `g^{DW}`, caption
  "Σ^{DW}: equal-time loop, static"] + [dGW: electron line `G` with a double dashed screened line `W^{e}` minus a thin
  `G^{cn}` with `W^{e,cn}`, caption "Σ^{dGW}: change of the electronic GW term (dropped)"].
- (c) Counting: [phonon self-energy: phonon line in, electron bubble, phonon line out; left vertex a small black dot
  labelled "bare", right vertex a grey disk labelled `g`, caption "Π^{ph} (Section A.6): one bare, one screened"]
  "versus" [Fan–Migdal again with both disks labelled `g`, caption "Σ^{FM}: both ends screened"].
- Labels in v2 notation only: no `D` for the phonon line (that is the review's symbol), no `W_e`, no `Γ` for the
  vertex, no `g^b`. Use `'electron'` for electron lines, `'phonon'` for phonon lines, `'photon'`/default for the
  screened Coulomb line (as `EPP.fd` does for `'screened'`).
- Each diagram gets an `aria-label` describing it in words.

**Controls.** None (diagram).

**Caption key points.**
- *What to notice.* Everything the moving nuclei add sits in $\Sigma^{\mathrm{ep}}$ on top of $G^{\mathrm{cn}}$; the
  Fan–Migdal arc carries frequency (dynamic), the Debye–Waller loop starts and ends at one time (static); the
  Fan–Migdal vertices are both screened, unlike the bubble of the phonon self-energy.
- *Marks.* Thin line: clamped-nuclei propagator; thick double line: full propagator; wavy line: phonon propagator
  $F$; double dashed line: interaction screened by electrons; grey disk: screened vertex $g$ (DFPT); black dot: bare
  vertex (the unscreened potential change $\Delta_{\bq\nu}V^{\mathrm{en}}$); hatched disk: second-order vertex
  $g^{\mathrm{DW}}$.
- *Try* (reading exercise, no controls). Follow the Fan–Migdal diagram from right to left: the electron emits a phonon
  at one screened vertex, propagates, and reabsorbs it at the other; in the bubble, one end of the phonon line meets
  the bare potential and the chain of further bubbles that screens it is already inside the other vertex.
- *Model.* (A7.3)–(A7.5) with $\Lambda=1$; after the review's Fig. 2; dGW is never computed and is dropped.

---

## A7.3 `fig-kink.html` · `viz-a7-kink` · port of `viz-a-electron-se-kink`

**Title.** "The ARPES kink: a band coupled to one Einstein phonon" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. In a metal ($E_{\mathrm F}\gg\hbar\omega_0$) the quasiparticle band is flattened by $Z^{\mathrm{qp}}\approx1/(1+\lambda)$
   within $\hbar\omega_0$ of $E_{\mathrm F}$ and sharp there; below $E_{\mathrm F}-\hbar\omega_0$ it returns to the bare
   band, broadened: the kink.
2. Where the slope of $\mathrm{Re}\,\Sigma$ exceeds one, the quasiparticle equation (A7.16) has three solutions at one
   $k$ (graphically: three crossings in the side panel), giving the S-shaped Brillouin–Wigner band.
3. When $E_{\mathrm F}\sim\hbar\omega_0$ the band is heavier and a polaron satellite appears below it that no solution of
   (A7.16) follows: the full spectral function (A7.15) carries new information.

**Include position.** In h3 `#a7-kink`, after its first paragraph and before the paragraph "The model of
[Figure A7.3], after the review, …". The sentences just before it: "In the real part, $\mathrm{Re}\,\Sigma$ falls
linearly inside the window, with the slope $-\lambda_{n\bk}$ of (A7.19), so the band flattens by $1/(1+\lambda_{n\bk})$;
it peaks at the window edge, and beyond it decays, so the band returns to its bare dispersion, broadened. The boundary at
$E_{\mathrm F}-\hbar\omega_0$ shows up in ARPES as a kink in the measured band. [Figure A7.3] computes it for a model
band."

**Physics** (keep v1's numerics, which are correct; verified again here).
- Band $E_k=\hbar^2k^2/2m^*$ from its bottom, $m^*=0.1\,m_{\mathrm e}$, $\hbar\omega_0=100$ meV, cutoff
  $E_{\mathrm c}=5$ eV above $E_{\mathrm F}$, $\hbar\eta=20$ meV, $T=0$, retarded Fan–Migdal only (no Debye–Waller),
  $|g|^2N_{\mathrm F}=\lambda_{\mathrm{ep}}\hbar\omega_0/2$ with $N_{\mathrm F}$ per spin.
- Energies measured from $E_{\mathrm F}$, $z=E+\ii\hbar\eta$ with $E=\hbar\omega-E_{\mathrm F}$.
- "2D, constant" density of states: Eq. (A7.20) `a7-eq-model`,
  $\Sigma=\frac{\lambda\hbar\omega_0}{2}\big[\ln\frac{z-\hbar\omega_0}{z-\hbar\omega_0-E_{\mathrm c}}+\ln\frac{z+\hbar\omega_0+E_{\mathrm F}}{z+\hbar\omega_0}\big]$
  (v1 `sigmaRaw`, branch `dim==='2'`).
- "3D, √E" density of states $N(E)=N_{\mathrm F}\sqrt{(E+E_{\mathrm F})/E_{\mathrm F}}$:
  $\Sigma=\frac{\lambda\hbar\omega_0}{2\sqrt{E_{\mathrm F}}}\big[I(z+\hbar\omega_0+E_{\mathrm F};E_{\mathrm F})+I(z-\hbar\omega_0+E_{\mathrm F};E_{\mathrm F}+E_{\mathrm c})-I(z-\hbar\omega_0+E_{\mathrm F};E_{\mathrm F})\big]$
  with $I(a;U)=\int_0^U\frac{\sqrt u}{a-u}\dd u=-2\sqrt U+\sqrt a\,\ln\frac{\sqrt a+\sqrt U}{\sqrt a-\sqrt U}$ (principal
  branches; v1 `I3`).
- Luttinger shift $\Sigma\to\Sigma-\mathrm{Re}\,\Sigma(0)$.
- Spectral function (A7.15) `a7-eq-spectral`: $\mathcal A(k,E)=\frac1\pi\frac{\hbar\eta-\mathrm{Im}\,\Sigma}{(E-E_k+E_{\mathrm F}-\mathrm{Re}\,\Sigma)^2+(\hbar\eta-\mathrm{Im}\,\Sigma)^2}$,
  clipped at 3 eV⁻¹ (as the review).
- Brillouin–Wigner band: $k(E)$ from $E_k-E_{\mathrm F}=E-\mathrm{Re}\,\Sigma(E)$ (A7.16); on-shell (Rayleigh–Schrödinger)
  band: $E=E_k-E_{\mathrm F}+\mathrm{Re}\,\Sigma(E_k-E_{\mathrm F})$; roots at the cursor $k$ by sign change and bisection.
- Verified numbers (3D, $\lambda=1.2$): $E_{\mathrm F}=10\hbar\omega_0$: $k_{\mathrm F}=0.162$ Å⁻¹,
  $\lambda_{n\bk}=-\partial\mathrm{Re}\,\Sigma/\partial E|_{E_{\mathrm F}}=1.11$ (the $\lambda_{n\bk}$ of (A7.19) at $k_{\mathrm F}$), $Z^{\mathrm{qp}}=0.47$, maximum slope
  of $\mathrm{Re}\,\Sigma$ 1.34 (S-shape present). $E_{\mathrm F}=2\hbar\omega_0$: at $k=0$ the spectral function peaks at
  $E-E_{\mathrm F}=-0.081$ eV (quasiparticle band bottom, $-0.81\hbar\omega_0$) and at $-0.247$ eV (satellite,
  $1.66\hbar\omega_0$ below the quasiparticle bottom). Constant DOS: the step of $|\mathrm{Im}\,\Sigma|$ outside the
  window is $\pi\lambda\hbar\omega_0/2=188$ meV. BW is multivalued where the slope of $\mathrm{Re}\,\Sigma$ exceeds one; its maximum is about $\lambda(\omega_0/4\eta-1/4)$ plus about 10% from the finite band, $1.115\lambda$ (3D) at the default, so the S needs $\lambda\gtrsim0.9$.

**What to draw** (as v1, renamed and relabelled).
- Left: heatmap of $\mathcal A(k,E)$, $k$ in Å⁻¹, vertical axis $E-E_{\mathrm F}$ (eV); colour `'ink'` ramp, colourbar
  "𝒜 (eV^{-1}), clipped at 3". Overlays: bare band $E_k$ (`'electron'` dashed, legend "bare band E_k"), Brillouin–Wigner
  quasiparticle band (`'electron'` solid, "quasiparticle band (A7.16)"), on-shell band (`'c4'` dotted, "on-shell
  (Rayleigh–Schrödinger)"), dashed guides at $0$ ("E_F"), $-\hbar\omega_0$ and $+\hbar\omega_0$ ("−ħω_0", "+ħω_0",
  `'phonon'`), the cursor $k$ as a vertical `'muted'` line, and hollow dots at the BW roots.
- Right: $\mathrm{Re}\,\Sigma$ (`'ink'`), $\mathrm{Im}\,\Sigma$ (`'c5'` dashed) in meV against the same energy axis, the line
  $E-E_k$ for the cursor $k$ (`'muted'` dotted, "E − E_k (cursor k)"), dots at its crossings with $\mathrm{Re}\,\Sigma$;
  fixed limits −450 to 350 meV.
- Window "zoom on E_F": $E-E_{\mathrm F}$ from $\max(-E_{\mathrm F}-3\hbar\omega_0,-5\hbar\omega_0)$ to $+2\hbar\omega_0$;
  "full band": from $-E_{\mathrm F}-3\hbar\omega_0$ to $+2\hbar\omega_0$, $k$ from 0 to $1.5k_{\mathrm F}$ (v1 logic).
- Dragging on the map sets the cursor (keep `touch-action: none` on the map only).

**Controls** (five; v1 had six).
1. Fermi energy $E_{\mathrm F}/\hbar\omega_0$ (measured from the band bottom): 0.5 to 15, step 0.5, default 10; the
   label also shows meV.
2. Coupling $\lambda_{\mathrm{ep}}$: 0 to 2, step 0.05, default 1.2.
3. Density of states (segmented): "3D, √E" (default) / "2D, constant".
4. Energy window (segmented): "zoom on E_F" (default) / "full band".
5. Cursor $k/k_{\mathrm F}$: 0 to 1.5, step 0.005, default 0.875.
The v1 broadening slider is dropped: $\hbar\eta=20$ meV fixed, the review's value.

**Readouts.** $k_{\mathrm F}$ (Å⁻¹); $\lambda_{n\bk}=-\partial\,\mathrm{Re}\,\Sigma/\partial E$ at $E_{\mathrm F}$ (no new symbol: the $\lambda_{n\bk}$ of (A7.19) at $k_{\mathrm F}$);
$Z^{\mathrm{qp}}$ at $E_{\mathrm F}$; $m^{\mathrm{qp}}/m^*=1+\lambda_{n\bk}$; maximum $\partial\,\mathrm{Re}\,\Sigma/\partial E$
(with "(three solutions)" when > 1); number of BW roots at the cursor.

**Initial state.** $E_{\mathrm F}=10\hbar\omega_0$, $\lambda_{\mathrm{ep}}=1.2$, 3D, zoom, cursor 0.875: the kink, the
sharp flattened band above $-\hbar\omega_0$, the broad band below, and the S-shape with three roots near the cursor.

**Caption key points.**
- *What to notice.* For $E_{\mathrm F}=10\hbar\omega_0$ the band is flattened by $Z^{\mathrm{qp}}\approx0.47$ within
  $\hbar\omega_0$ of $E_{\mathrm F}$ and sharp, then snaps back to the bare band, broad, below $-\hbar\omega_0$; where the
  slope of $\mathrm{Re}\,\Sigma$ exceeds one the quasiparticle band is S-shaped. For $E_{\mathrm F}=2\hbar\omega_0$ the band
  is heavier and a broad satellite about $1.7\hbar\omega_0$ below its bottom has no quasiparticle solution.
- *Marks.* Heatmap $\mathcal A(k,E)$ clipped at 3 eV⁻¹; bare, Brillouin–Wigner and on-shell bands; guides at $E_{\mathrm F}$
  and $\pm\hbar\omega_0$; side panel $\mathrm{Re}\,\Sigma$, $\mathrm{Im}\,\Sigma$ and the line $E-E_k$ whose crossings are
  the solutions of (A7.16).
- *Try.* (1) Set $E_{\mathrm F}/\hbar\omega_0=2$ and "full band": find the satellite below the quasiparticle band.
  (2) Lower $\lambda_{\mathrm{ep}}$ below about 0.9: the S-shape disappears (maximum slope below one). (3) Drag the cursor
  across the kink and watch the crossings in the side panel go from one to three and back.
- *Model.* (A7.20) and its √E variant, (A7.15), (A7.16), (A7.18); $m^*=0.1\,m_{\mathrm e}$, $\hbar\omega_0=100$ meV,
  $\hbar\eta=20$ meV, $E_{\mathrm c}=5$ eV, $|g|^2N_{\mathrm F}=\lambda_{\mathrm{ep}}\hbar\omega_0/2$ ($N_{\mathrm F}$ per spin),
  $T=0$, Luttinger shift, no Debye–Waller; after the review's Fig. 4. The readout $\lambda_{n\bk}$ differs slightly from the
  slider because the band is finite and not particle–hole symmetric. $E_{\mathrm F}=10\hbar\omega_0$ mimics a metal,
  $2\hbar\omega_0$ a degenerately doped semiconductor.

**What in v1 must change.**
- Id `viz-a-electron-se-kink` → `viz-a7-kink`; classes `pk-*` → `a7-kink-*`.
- Source notation everywhere: `ε_k` → $E_k$ (and the aria-label "E minus epsilon_k" → "E minus E_k"); `ε_F` → $E_{\mathrm F}$
  measured from the band bottom; `A(k,E)` → $\mathcal A$; "Brillouin–Wigner (G.164)" → "quasiparticle band (A7.16)";
  readouts `Z(k_F)` → $Z^{\mathrm{qp}}$ at $E_{\mathrm F}$, `m_{ep}/m^*` → $m^{\mathrm{qp}}/m^*$; all "(G.xxx)" references
  → `a.eqref` to the v2 equations.
- v1's caption is one block that mixes model and observations; rewrite in the four-part structure.
- The remark on the review's caption ($\varepsilon_{\max}$ vs $\varepsilon_c$, $N_F$ per spin, $\lambda=2$ or 1) lives in the
  section's footnote `fn-a7-fig4`; the caption only states "$N_{\mathrm F}$ per spin".
- Drop the broadening slider (six controls → five).
- v1 describes the 3D density of states as ∝ √(ε+ε_F) with energies from $E_{\mathrm F}$; write it as ∝ √E from the
  band bottom.

---

## A7.4 `fig-cumulant.html` · `viz-a7-cumulant` · port of `viz-a-electron-se-cumulant`

**Title.** "Dyson against cumulant for one level and one phonon" `<span class="tag">interactive</span>`.

**What the reader must learn.**
1. Both approaches give the same quasiparticle shift $S\hbar\omega_0$ and weight $1-S$ at small $S$.
2. At larger $S$ the Dyson equation with the Fan–Migdal self-energy gives exactly two poles: one satellite,
   $\hbar\omega_0\sqrt{1+4S}$ from the quasiparticle (too far), and too much quasiparticle weight ($Z_+$ instead of
   $\ee^{-S}$).
3. The cumulant gives a Poisson (Franck–Condon) ladder of replicas spaced by exactly $\hbar\omega_0$ (A7.27), the
   polaron satellites of ARPES; for this model it is exact.

**Include position.** In h3 `#a7-cumulant`, after the paragraph that introduces the single-level model and before the
paragraph "For this model both approaches can be solved exactly …". The sentences just before it: "A level of energy
$E_0$ is coupled with a constant matrix element $g_0$ to a dispersionless phonon of frequency $\omega_0$; the strength
of the coupling is the Huang–Rhys factor $S=|g_0|^2/(\hbar\omega_0)^2$, the mean number of phonons in the distortion
cloud. Removing the electron leaves the lattice distorted, and the distortion relaxes by emitting phonons.
[Figure A7.4] compares the two approaches."

**Physics** (keep v1's numerics, which are correct). Reduced variable $x=(E-E_0)/\hbar\omega_0$ with $E=\hbar\omega$;
$r=\hbar\eta/\hbar\omega_0$.
- Dyson with the Fan–Migdal self-energy of the occupied level (derivation box of `#a7-cumulant`):
  $\Sigma/\hbar\omega_0=S/(x+1+\ii r)$, $G=1/(x+\ii r-\Sigma/\hbar\omega_0)$, plotted as $-r\,\mathrm{Im}\,G$
  ($=\pi\hbar\eta\,\mathcal A$). Poles $x_\pm=\frac12(-1\pm\sqrt{1+4S})$, weights
  $Z_\pm=[1+S/(x_\pm+1)^2]^{-1}$.
- Cumulant (A7.27) `a7-eq-poisson`: $\pi\hbar\eta\,\mathcal A=\sum_s\ee^{-S}\frac{S^s}{s!}\frac{r^2}{(x-S+s)^2+r^2}$,
  summed to $s_{\max}=\max(8,\lceil S+8\sqrt{S+0.5}+6\rceil)$.
- Check values at the initial state ($S=1$): Dyson $x_+=0.618$, $Z_+=0.724$, satellite at $x_-=-1.618$ (gap
  $\sqrt5\,\hbar\omega_0=2.24\hbar\omega_0$); cumulant quasiparticle at $x=1$ with weight $\ee^{-1}=0.368$, replicas at
  $x=0,-1,-2,\dots$
- Model: independent-boson (Langreth core-hole) model; the hole branch measured by ARPES.

**What to draw** (as v1). One plot, horizontal axis $(E-E_0)/\hbar\omega_0$ from −9 to 5 (fixed), vertical axis
$\pi\hbar\eta\,\mathcal A(E)$ from 0 to 1.08 (fixed), so that an isolated peak has height equal to its weight. Curves:
cumulant (`'electron'` solid, "cumulant (exact for this model)"), Dyson (`'c4'` dashed, "Dyson with Fan–Migdal Σ").
Markers at the exact weights: hollow `'electron'` circles for the replicas $\ee^{-S}S^s/s!$ at $x=S-s$, filled `'c4'`
dots for $Z_\pm$ at $x_\pm$. Vertical `'muted'` line at $x=0$ labelled "bare level E_0". Text "Dyson satellite" next to
$x_-$ when $Z_->0.03$.

**Controls.**
1. Huang–Rhys factor $S$: 0 to 4, step 0.05, default 1.
2. Phonon energy $\hbar\omega_0$: 20 to 200 meV, step 5, default 100 meV.
3. Broadening $\hbar\eta$: 5 to 40 meV, step 1, default 12 meV.

**Readouts.** $Z^{\mathrm{qp}}$ Dyson ($Z_+$); $Z^{\mathrm{qp}}$ cumulant ($\ee^{-S}$); quasiparticle shift Dyson
($x_+\hbar\omega_0$, meV); shift cumulant ($S\hbar\omega_0$, meV); quasiparticle–satellite gap Dyson
($\hbar\omega_0\sqrt{1+4S}$) and cumulant ($\hbar\omega_0$), "no satellite" at $S=0$.

**Initial state.** $S=1$, $\hbar\omega_0=100$ meV, $\hbar\eta=12$ meV: the Dyson satellite at $-2.24\hbar\omega_0$ from
the quasiparticle against the cumulant ladder.

**Caption key points.**
- *What to notice.* At small $S$ both give the shift $S\hbar\omega_0$ toward the Fermi level. As $S$ grows, Dyson keeps
  one satellite, too far from the quasiparticle, and too much weight in the quasiparticle; the cumulant moves weight
  into replicas spaced by exactly $\hbar\omega_0$.
- *Marks.* Solid: cumulant; dashed: Dyson; hollow circles and filled dots: exact weights; vertical line: the bare level;
  each peak broadened by a Lorentzian of half-width $\hbar\eta$, scaled so that its height is its weight.
- *Try.* (1) $S=0.25$: shifts 20.7 and 25 meV, weights 0.85 and 0.78 (the "Check yourself" numbers). (2) $S=3$: the
  ladder's maximum moves to the second and third replicas while Dyson shows one far satellite. (3) Change
  $\hbar\omega_0$: the spacing of the replicas follows it exactly.
- *Model.* Derivation box of `#a7-cumulant`; (A7.26), (A7.27); one occupied level coupled with constant $g_0$ to one
  dispersionless phonon; $S=|g_0|^2/(\hbar\omega_0)^2$; the second-order cumulant is exact (Langreth 1970). The same
  failure of Dyson and the same cure appear for plasmon satellites in GW (GW versus GW plus cumulant). TiO₂ ARPES shows
  such replicas.

**What in v1 must change.**
- Id `viz-a-electron-se-cumulant` → `viz-a7-cumulant`; classes `pc-*` → `a7-cumulant-*`.
- `ε_0` → $E_0$ in axis label, legend ("bare level E_0"), aria-label and caption; axis `(E − ε_0)/ħω_0` →
  `(E − E_0)/ħω_0`; vertical axis `πη A(E)` → `πħη 𝒜(E)`.
- Readout labels `Z_{QP}` → $Z^{\mathrm{qp}}$; "QP" spelled "quasiparticle".
- Caption rewritten in the four-part structure, with links to (A7.26), (A7.27) and the section's derivation box.
