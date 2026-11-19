# Section A.5 · figure specifications

Five figures, in reading order. Read `v2/build/BUILDER_GUIDE.md` §5 (figures), `v2/NOTATION.md`, and the section file
`v2/sections/A05-hedin-baym.html` first; the equations named below are there, with their ids. Every symbol in labels,
readouts, legends, aria-labels and strings is v2 notation: no source symbols (no $\mathcal D$, $\Gamma$ for the vertex,
$P$, $P_{\rm e}$, $W_{\rm e}$, $W_{\rm ph}$, $\epsilon$, $\epsilon^0$, $\epsilon^\infty$, $\Omega_p$, $\Pi^{\rm ph}$ for the
atom-basis force constants, no "G Eq." references). Captions refer to equations only with
`<a class="eqref" href="#a5-eq-…"></a>`.

Colour convention, shared by all five figures and matching the prose:
- `electron` (blue): electronic quantities, $W^{\rme}$, electrons;
- `phonon` (amber): nuclear quantities, $W^{\ph}$, ions, the phonon propagator $F$;
- `photon` (rose): the total screened interaction $W$ and the bare Coulomb $v$ (Part B turns $v$ into the photon line);
- `ink` / `muted` for guides, `c4`/`c5` for secondary curves.

Equation ids used below: `a5-eq-toy` (A5.4), `a5-eq-osc` (A5.5), `a5-eq-epsw` (A5.6), `a5-eq-Wpolar` (A5.7),
`a5-eq-wp` (A5.8), `a5-eq-epsjel` (A5.9), `a5-eq-Wjel` (A5.10), `a5-eq-bs` (A5.11), `a5-eq-dysonG` (A5.13),
`a5-eq-Vtot` (A5.14), `a5-eq-Sigma` (A5.17), `a5-eq-vertex` (A5.18), `a5-eq-W` (A5.19), `a5-eq-We` (A5.20),
`a5-eq-chiion` (A5.22), `a5-eq-split` (A5.23), `a5-eq-dn` (A5.25), `a5-eq-F` (A5.26), `a5-eq-Wph` (A5.28),
`a5-eq-frohlich` (A5.30), `a5-eq-Phiw` (A5.31), `a5-eq-Finv` (A5.33), `a5-eq-nion` (A5.34), `a5-eq-hedin` (A5.35),
`a5-eq-v` (A5.3). Other sections: `a4-eq-frohlich`, `a4-eq-lot`, `b5-eq-D0L`.

---

## 1. `fig-exchange.html` · `viz-a5-exchange` · new

**Title.** "Two electrons coupled by a row of ions" `<span class="tag">interactive</span>`.

**What the reader must learn.** The lattice part of the screened interaction has the anatomy coupling × oscillator ×
coupling, (A5.4): electron 1 pushes each ion, each ion responds as a driven oscillator, and electron 2 feels the
displaced ions. Below the ions' frequency they move with the force and gather positive charge near electron 1, so
$W^{\ph}<0$ (attraction); above it they move against the force and $W^{\ph}>0$; far above, they cannot follow. Whether
the total $W=W^{\rme}+W^{\ph}$ is attractive depends on how $W^{\ph}$ compares with the screened repulsion $W^{\rme}$.

**Include position.** In `a5-toy`, right after the paragraph that ends: "With the parameters of Figure A5.1 (sodium mass,
$\hbar\omega_0=15$ meV, illustrative), the ions reduce the static repulsion at 2 Å from 0.97 to 0.54 eV, and at
$\omega=0.8\,\omega_0$ they reverse it, to $-0.24$ eV." The prose quotes these numbers, so the figure must reproduce them.

**Physics** (all lengths in Å, energies in eV, masses in amu).
- Ions at $\mathbf R_\ell=(\ell a,0)$, $\ell=-10,\dots,10$ (21 ions; draw about 9 of them, the ones with $|x|\le12$ Å),
  $a=3.0$ Å, charge $Ze$ with $Z=1$, mass $M=22.99$ amu (sodium), each an isotropic two-dimensional Einstein oscillator
  (motion in the plane of the drawing, $x$ and $y$) of frequency $\hbar\omega_0=15$ meV, damping $\hbar\gamma=0.3$ meV.
- Electrons on the line $y=h=2.0$ Å: electron 1 fixed at $(x_1,h)$ with $x_1=0$ (directly above ion 0); electron 2 at
  $(x_2,h)$.
- Screened interaction between two unit charges: $W^{\rme}(r)=(e^2/r)\,\ee^{-k_{\mathrm{TF}}r}$, $e^2=14.3996$ eV Å.
  Electron–ion energy $-Z\,W^{\rme}(d)$ with $d$ the distance.
- Force on ion $\ell$ exerted by an electron at $(x,h)$: $\bF_\ell(x)=Z\,W^{\rme\prime}(d)\,(\mathbf R_\ell-\mathbf r_{\rm e})/d$,
  with $W^{\rme\prime}(d)=-W^{\rme}(d)\,(1/d+k_{\mathrm{TF}})$; it points from the ion toward the electron.
- Displacement of ion $\ell$ (complex amplitude, driven at $\omega$):
  $\bu_\ell=\bF_\ell(x_1)/[M(\omega_0^2-\omega^2-\ii\gamma\omega)]$.
- Lattice part, (A5.4) with damping:
  $W^{\ph}(x_2,x_1;\omega)=\sum_\ell\bF_\ell(x_2)\cdot\bF_\ell(x_1)\,/\,[M(\omega^2-\omega_0^2+\ii\gamma\omega)]$; plot its real
  part. Units: $M\omega^2$ [eV/Å²] $=(\hbar\omega\,[\text{eV}])^2\,M[\text{amu}]/4.1802\times10^{-3}$, so $M\omega_0^2=1.237$ eV/Å².
- $W^{\rme}(x_2,x_1)=W^{\rme}(|x_2-x_1|)$ (both electrons on the same line); total $W=W^{\rme}+\mathrm{Re}\,W^{\ph}$.
- Check values ($k_{\mathrm{TF}}=1.0$ Å⁻¹, $x_2=2.0$ Å): $W^{\rme}=0.974$ eV; $\mathrm{Re}\,W^{\ph}=-0.44$ eV at $\omega=0$,
  $-1.22$ eV at $\omega=0.8\,\omega_0$, $+2.07$ eV at $\omega=1.1\,\omega_0$; total $W=0.54$, $-0.24$, $+3.04$ eV. At
  $\omega=0.8\,\omega_0$ the total is negative from about 1.65 to 4.4 Å.

**What to draw.**
- Top panel (canvas, geometric, about 2.2:1): the row of ions as amber circles with "+" (radius about 0.3 Å on the drawing
  scale), their equilibrium sites as faint crosses; the electron line $y=h$ as a faint dashed line; electron 1 as a filled
  blue dot labelled "1" and electron 2 as a blue ring labelled "2", draggable along the line (`touch-action: none` only
  on this canvas). On each ion: an amber arrow for the force $\bF_\ell(x_1)$ of electron 1 and a dark arrow for the
  displacement $\mathrm{Re}\,\bu_\ell$ (same exaggeration factor for all ions, fixed, chosen so that the largest arrow
  at the initial state is about 0.8 $a$). The displacement arrows flip when $\omega$ crosses $\omega_0$. Optionally, faint
  lines from electron 2 to the three nearest ions.
- Bottom panel (`EPP.Plot`): $x$ axis "$x_2-x_1$ (Å)" from $-12$ to $12$; $y$ axis "potential energy (eV)" fixed at
  $[-3,3]$ (clip; $W^{\rme}$ diverges at 0). Curves: $W^{\rme}$ (blue), $\mathrm{Re}\,W^{\ph}$ (amber), $\mathrm{Re}\,W$
  (rose, thicker); zero line; a vertical marker at the current $x_2$; light rose shading where $\mathrm{Re}\,W<0$.
- Labels in v2 notation via `drawRich`: "W^{e}", "Re W^{ph}", "Re W", "x_1", "x_2", "**F**_ℓ" (bold), "**u**_ℓ" (bold).

**Controls** (4).
1. Exchanged energy $\hbar\omega$: range 0–30 meV (0 to $2\omega_0$), step 0.1 meV, readout also as $\omega/\omega_0$.
2. Screening wavevector $k_{\mathrm{TF}}$: 0.5–2.0 Å⁻¹, step 0.05.
3. Position of electron 2, $x_2$: −12 to 12 Å, step 0.1 (also by dragging electron 2; keyboard: arrow keys on the slider).
4. Play/Pause (`EPP.anim`, paused by default): the ions oscillate as $\mathrm{Re}[\bu_\ell\ee^{-\ii\omega t}]$ (slowed
   down) while a halo on electron 1 pulses as $\cos\omega t$; at $\omega=0$ the animation shows the static displacement.

**Readouts.** At the current $x_2$: $W^{\rme}$, $\mathrm{Re}\,W^{\ph}$, $\mathrm{Re}\,W$ (eV, 2 decimals), with
"attractive"/"repulsive" for the total; $\omega/\omega_0$; "ions move with the force" ($\omega<\omega_0$) or "against the
force" ($\omega>\omega_0$). MathJax in HTML readouts.

**Initial state.** $\hbar\omega=12$ meV ($0.8\,\omega_0$), $k_{\mathrm{TF}}=1.0$ Å⁻¹, $x_2=2.0$ Å, paused. Readouts
$0.97$, $-1.22$, $-0.24$ eV: the attractive window is visible at once.

**Caption.**
- *What to notice.* The ions add a sum of (force from electron 1) × (displacement per unit force) × (force on electron 2),
  (A5.4). Below $\omega_0$ the ions move toward electron 1 and $W^{\ph}<0$; above $\omega_0$ they move against the force and
  the sign flips.
- *Marks.* Blue: electrons (1 is the source, 2 the probe) and $W^{\rme}$; amber: ions, the forces $\bF_\ell(x_1)$ and
  $\mathrm{Re}\,W^{\ph}$; dark arrows: ion displacements; rose: the total $W$, shaded where it is negative.
- *Try.* Sweep $\hbar\omega$ through 15 meV and watch the displacement arrows flip and $W^{\ph}$ change sign; set
  $\hbar\omega=0$ to see the ions halve the repulsion without reversing it; lower $k_{\mathrm{TF}}$ (weaker screening)
  and see both $W^{\rme}$ and $W^{\ph}$ grow; drag electron 2 to find where the total is attractive.
- *Model.* (A5.4) with damping; Einstein oscillators, a flat optical branch, so each ion is its own phonon mode and
  (A5.28) reduces to (A5.4) with the vertex equal to the zero-point amplitude times the force. Screened Coulomb
  $W^{\rme}(r)=(e^2/r)\ee^{-k_{\mathrm{TF}}r}$; $a=3.0$ Å, $h=2.0$ Å, $Z=1$, $M=22.99$ amu, $\hbar\omega_0=15$ meV,
  $\hbar\gamma=0.3$ meV. A toy with illustrative numbers, not a model of a real metal.

---

## 2. `fig-polar.html` · `viz-a5-polar` · port of v1 `viz-a-hedin-baym-polar` (v1 `sections/A05-hedin-baym.html`)

**Title.** "Lattice screening in a polar insulator" `<span class="tag">interactive</span>`.

**What the reader must learn.** In a polar insulator at long wavelength the electrons screen with a constant
$\varepsilon_\infty$, while the ions add a resonance: $W=v/\varepsilon(\omega)$ has a zero at $\omega_{\mathrm{TO}}$ and a pole at
$\omega_{\mathrm{LO}}$; statically the lattice screens more ($W=v/\varepsilon_{\mathrm s}$); $W^{\ph}<0$ below $\omega_{\mathrm{LO}}$;
the loss peaks at LO, the infrared absorption at TO; and $W^{\ph}$ is exactly Fröhlich's coupling squared times the LO
propagator, (A5.30).

**Include position.** End of `a5-polar`, after the paragraph: "The residue of the pole, $v(q)(1/\varepsilon_\infty-1/\varepsilon_{\mathrm s})$,
turns out to be the square of Fröhlich's matrix element: a later subsection shows that $W^{\ph}$ is Fröhlich's coupling
squared times the propagator of an LO phonon."

**Physics** (keep the v1 code; it is correct). With $E=\hbar\omega$ and damping $\hbar\gamma$:
- $\varepsilon(\omega)=\varepsilon_\infty(E_{\mathrm{LO}}^2-E^2-\ii\gamma E)/(E_{\mathrm{TO}}^2-E^2-\ii\gamma E)$, (A5.6) with damping
  in the convention $\omega^2-\omega_0^2+\ii\gamma\omega$; $\varepsilon_{\mathrm s}=\varepsilon_\infty E_{\mathrm{LO}}^2/E_{\mathrm{TO}}^2$.
- $W/v=1/\varepsilon$; $W^{\rme}/v=1/\varepsilon_\infty$; $W^{\ph}/v=1/\varepsilon-1/\varepsilon_\infty$, (A5.7).
- Fröhlich overlay, undamped: $\OV|g^{\mathrm{lr}}|^2F^0_{\mathrm{LO}}/v=(1/\varepsilon_\infty-1/\varepsilon_{\mathrm s})\,E_{\mathrm{LO}}^2/(E^2-E_{\mathrm{LO}}^2)$, (A5.30).
- Loss $-\mathrm{Im}(1/\varepsilon)$ and infrared absorption $\mathrm{Im}\,\varepsilon$ (both positive for $E>0$).
- Presets (literature values, illustrative): GaAs $E_{\mathrm{TO}}=33.3$, $E_{\mathrm{LO}}=36.2$ meV, $\varepsilon_\infty=10.9$,
  $\OV=45.2$ Å³; MgO 49.1, 89.0 meV, 2.96, 18.7 Å³; bulk hBN, in-plane modes, 169.9, 199.6 meV, 4.87, 18.1 Å³ per BN
  pair.
- Readouts: $\varepsilon_{\mathrm s}$ (LST); $\hbar\omega_{\mathrm p}^{\ion}=\sqrt{\varepsilon_\infty(E_{\mathrm{LO}}^2-E_{\mathrm{TO}}^2)}$ (new; GaAs
  46.9 meV, MgO 127.7, hBN 231.2); $W^{\rme}/v=1/\varepsilon_\infty$; $W(0)/v$ ($1/\varepsilon_{\mathrm s}$, or $1/\varepsilon_\infty$ when
  clamped); $\OV|g^{\mathrm{lr}}|^2/v(q)=(E_{\mathrm{LO}}/2)(1/\varepsilon_\infty-1/\varepsilon_{\mathrm s})$ in meV (GaAs 0.255);
  $|g^{\mathrm{lr}}|\,q=[2\pi e^2E_{\mathrm{LO}}(1/\varepsilon_\infty-1/\varepsilon_{\mathrm s})/\OV]^{1/2}$ in meV Å⁻¹ (GaAs 32.0; v1's
  `Math.sqrt(180.95/om*(elo/1000)/2*(ie-i0))*1000` is the same number).

**What to draw.** As v1: (a) real parts of $W/v$ (rose, total), $W^{\rme}/v$ (blue), $\mathrm{Re}\,W^{\ph}/v$ (amber), and the
dashed Fröhlich overlay (`c4`), versus $\hbar\omega$ from 0 to $2E_{\mathrm{LO}}$, fixed $y$ range $[-2.5,3.5]/\varepsilon_\infty$;
grey horizontal guide at $1/\varepsilon_{\mathrm s}$ labelled "1/ε_s"; vertical lines at TO and LO, labelled "ω_{TO}",
"ω_{LO}"; the Reststrahlen band TO–LO lightly shaded and labelled "Reststrahlen". (b) Loss $-\mathrm{Im}(1/\varepsilon)$
(rose) and $\mathrm{Im}\,\varepsilon$ (`c5`, dashed), log scale $10^{-4}$–$10^4$, same marks.

**Controls** (4, as v1). Material (GaAs, MgO, hBN in-plane); damping $\hbar\gamma$ 0.1–5 meV (log); toggle "clamp the
nuclei ($F=0$)", which leaves the flat $W^{\rme}/v$ and removes the lattice curves; toggle "show Fröhlich form
$\OV|g^{\mathrm{lr}}|^2F^0_{\mathrm{LO}}$".

**Initial state.** GaAs, $\hbar\gamma=0.8$ meV, nuclei free, Fröhlich overlay on.

**Caption.**
- *What to notice.* Statically the lattice screens more, $W(0)=v/\varepsilon_{\mathrm s}<W^{\rme}$; $W$ vanishes at TO and has
  its pole at LO, so the loss peaks at LO while the infrared absorption peaks at TO. The dashed Fröhlich curve lies on
  $W^{\ph}$: Fröhlich's coupling is the residue of $W$ at the LO pole.
- *Marks.* Rose total $W$, blue $W^{\rme}$, amber $W^{\ph}$, dashed Fröhlich form; guides at TO, LO and $1/\varepsilon_{\mathrm s}$;
  shaded Reststrahlen band, where even the total $W$ is negative.
- *Try.* Switch to MgO, whose large LO–TO splitting gives $\varepsilon_{\mathrm s}/\varepsilon_\infty=3.3$; clamp the nuclei to see
  the flat $W^{\rme}$ of a GW calculation; raise the damping and watch the pole broaden while the Fröhlich (undamped)
  curve stays sharp.
- *Model.* (A5.6), (A5.7) with damping; (A5.30) for the overlay; materials as above; $\varepsilon_{\mathrm s}$ from the
  Lyddane–Sachs–Teller relation. Long wavelength, homogeneous medium, one polar mode.

**Change from v1.**
- Ids and classes: `viz-a5-polar`, `a5-polar-…` (v1 used `hb-pa`, `hb-pb`).
- Strings: "W_e/v (clamped nuclei)" → "W^{e}/v = 1/ε_∞ (nuclei clamped)"; "Re W_{ph}/v" → "Re W^{ph}/v";
  "Fröhlich |g|^2𝒟^A_{LO}/ħv" → "Fröhlich Ω_V|g^{lr}|^2F^0_{LO}/v"; "1/ε^0" → "1/ε_s"; "𝒟 = 0: W = W_e" → "F = 0: W = W^{e}";
  readout keys "ε^0 (LST)" → "ε_s (LST)", "W_e/v = 1/ε^∞" → "W^e/v = 1/ε_∞", "|g|^2/v" → "Ω_V|g^{lr}|²/v(q)",
  "|g|·q" → "|g^{lr}| q"; toggle labels as above. Aria-labels likewise.
- Caption: v1's single paragraph mixed "Model" first, SI units and "(G Eq. 132)"; rewrite in the four-part structure,
  Gaussian units (the plotted ratios are unit-free), v2 symbols, eqrefs only.
- Add the $\hbar\omega_{\mathrm p}^{\ion}$ readout; add the "Reststrahlen" label.

---

## 3. `fig-jellium.html` · `viz-a5-jellium` · port of v1 `viz-a-hedin-baym-jellium` (v1 `sections/A05-hedin-baym.html`)

**Title.** "Jellium with mobile ions" `<span class="tag">interactive</span>`.

**What the reader must learn.** Below the phonon frequency $\omega_{\mathrm{LA}}(q)$ the total interaction between electrons
is negative at every $q$: the slow ions overscreen. Electronic screening turns the ionic plasma oscillation (flat at
$\omega_{\mathrm p}^{\ion}$) into an acoustic branch with the Bohm–Staver velocity $\omega_{\mathrm p}^{\ion}/k_{\mathrm{TF}}$; the
electronic plasma frequency $\omega_{\mathrm p}$ is two orders of magnitude higher.

**Include position.** End of `a5-jellium`, after the paragraph that ends: "With the full Lindhard dielectric function in
place of Thomas–Fermi, the dispersion acquires the Kohn anomaly at $q=2k_{\mathrm F}$ (Section A.6)."

**Physics** (keep the v1 code; it is correct).
- Free-electron parameters from $r_{\mathrm s}$: $k_{\mathrm F}=1.91916/(r_{\mathrm s}a_0)$, $a_0=0.529177$ Å; $n^{\rme}_0=k_{\mathrm F}^3/3\pi^2$;
  $\OV=Z/n^{\rme}_0$; $k_{\mathrm{TF}}^2=4k_{\mathrm F}/(\pi a_0)$.
- $\hbar\omega_{\mathrm p}^{\ion}=[4\pi Z^2e^2\hbar^2/(M\OV)]^{1/2}$, numerically $\sqrt{0.75642\,Z^2/(\OV[\text{Å}^3]\,M[\text{amu}])}$ eV,
  (A5.8). New readout $\hbar\omega_{\mathrm p}=[4\pi n^{\rme}_0e^2\hbar^2/m_{\mathrm e}]^{1/2}=\sqrt{4\pi n^{\rme}_0\cdot14.3996\cdot7.61996}$ eV.
- Electronic screening, static: Lindhard $\varepsilon^{\rme}(q)=1+(k_{\mathrm{TF}}/q)^2[\frac12+\frac{1-x^2}{4x}\ln|\frac{1+x}{1-x}|]$,
  $x=q/2k_{\mathrm F}$ (value $\frac12$ for the bracket at $x=1$); Thomas–Fermi $1+(k_{\mathrm{TF}}/q)^2$; or none, $\varepsilon^{\rme}=1$.
- $\omega_{\mathrm{LA}}^2(q)=(\omega_{\mathrm p}^{\ion})^2/\varepsilon^{\rme}(q)$; with damped free ions,
  $W=\big(v/\varepsilon^{\rme}\big)(\omega^2+\ii\gamma\omega)/(\omega^2+\ii\gamma\omega-\omega_{\mathrm{LA}}^2)$, (A5.9)–(A5.10);
  $W^{\rme}=v/\varepsilon^{\rme}$; $W^{\ph}=W-W^{\rme}$. Plot in units of $4\pi e^2/k_{\mathrm{TF}}^2$, in which
  $v/\varepsilon^{\rme}=k_{\mathrm{TF}}^2/[q^2\varepsilon^{\rme}(q)]$.
- $v_{\mathrm s}=\omega_{\mathrm p}^{\ion}/k_{\mathrm{TF}}$ in km/s, (A5.11).
- Presets (free-electron values, illustrative): Na $r_{\mathrm s}=3.93$, $Z=1$, $M=22.99$ amu; K 4.86, 1, 39.10; Al 2.07, 3,
  26.98; Pb 2.30, 4, 207.2. Na: $\hbar\omega_{\mathrm p}^{\ion}=29.6$ meV, $\hbar\omega_{\mathrm p}=6.05$ eV, $k_{\mathrm F}=0.923$ Å⁻¹,
  $k_{\mathrm{TF}}=1.49$ Å⁻¹, $v_{\mathrm s}=3.0$ km/s; at $q=0.6k_{\mathrm F}$ with Lindhard $\varepsilon^{\rme}=8.02$,
  $\hbar\omega_{\mathrm{LA}}=10.4$ meV. Al: 123.6 meV, 15.8 eV, 9.1 km/s.

**What to draw.** As v1: (a) heat map of $\mathrm{Re}\,W$ over $(q/k_{\mathrm F},\hbar\omega)$, $q/k_{\mathrm F}\in[0.02,2.5]$,
$\hbar\omega\in[0,1.3\,\hbar\omega_{\mathrm p}^{\ion}]$, rose for repulsion and teal (`c4`) for attraction, clipped at $\pm2$, with a
colour bar "Re W (units of 4πe²/k_{TF}²)"; the dashed phonon line $\hbar\omega_{\mathrm{LA}}(q)$; a vertical "cut" line at
the chosen $q$; click or drag on the map to choose $q$. (b) Cut at that $q$ versus $\hbar\omega$: $W^{\rme}$ (blue),
$\mathrm{Re}\,W^{\ph}$ (amber), $\mathrm{Re}\,W$ (rose); band below $\hbar\omega_{\mathrm{LA}}$ shaded teal and labelled
"W < 0"; vertical guide at $\hbar\omega_{\mathrm{LA}}$ labelled "ħω_{LA}".

**Controls** (4, as v1). Metal (Na, K, Al, Pb); wavevector $q/k_{\mathrm F}$ 0.05–2.5; electronic screening (Lindhard,
Thomas–Fermi, none); damping $\hbar\gamma$ 0.1–3 meV (log).

**Readouts.** $\hbar\omega_{\mathrm p}^{\ion}$, $\hbar\omega_{\mathrm p}$ (new), $k_{\mathrm F}$, $\varepsilon^{\rme}(q)$,
$\hbar\omega_{\mathrm{LA}}(q)$, $v_{\mathrm s}=\omega_{\mathrm p}^{\ion}/k_{\mathrm{TF}}$ ("not defined" without screening).

**Initial state.** Na, $q/k_{\mathrm F}=0.6$, Lindhard, $\hbar\gamma=0.4$ meV.

**Caption.**
- *What to notice.* Above the phonon line the interaction is the electron-screened repulsion $W^{\rme}$; below it the
  ions overscreen and $W$ turns attractive, at every $q$. Screening makes the ions acoustic.
- *Marks.* Heat map: rose repulsion, teal attraction; dashed line $\hbar\omega_{\mathrm{LA}}(q)$; cut panel: blue $W^{\rme}$,
  amber $W^{\ph}$, rose $W$, teal band where $W<0$.
- *Try.* Choose "none" for the screening: the ions oscillate at $\omega_{\mathrm p}^{\ion}$ for every $q$; switch to Thomas–Fermi
  and Lindhard to see sound appear; compare Al and Pb, whose heavy ions give a much lower phonon energy; note how far
  below $\hbar\omega_{\mathrm p}$ (eV) all of this happens.
- *Model.* (A5.9)–(A5.11): jellium with free, damped ions and static electronic screening, valid for
  $\omega\ll v_{\mathrm F}q$; free-electron parameters as above. Bardeen–Pines interaction; Bohm–Staver sound.

**Change from v1.**
- Ids and classes: `viz-a5-jellium`, `a5-jellium-…` (v1 used `hb-ja`, `hb-jb`).
- Strings: "ħω_q (phonon)" → "ħω_{LA}(q) (phonon)"; "W_e (screened Coulomb)" → "W^{e} (screened by electrons)";
  "Re W_{ph}" → "Re W^{ph}"; "W / W_{TF}" and "Re W/W_{TF}" → "W (units of 4πe²/k_{TF}²)" (no new symbol $W_{\rm TF}$);
  readout keys "ħΩ_p" → "ħω_p^{ion}", "ε_e(q)" → "ε^e(q)", "ħω_q" → "ħω_{LA}(q)", "Ω_p/k_{TF}" → "v_s = ω_p^{ion}/k_{TF}";
  label "metal (free-electron r_s)" can stay.
- Caption: remove "SI" and "(G Eq. 60)"; four-part structure; eqrefs only.
- Add the $\hbar\omega_{\mathrm p}$ readout.

---

## 4. `fig-fd.html` · `viz-a5-fd` · port of v1 `viz-a-hedin-baym-fd` (v1 `sections/A05-hedin-baym.html`)

**Title.** "The lattice inside W" `<span class="tag">interactive</span>` (two toggles; use "diagram" if the toggles are
dropped).

**What the reader must learn.** $W=W^{\rme}+W^{\ph}$, (A5.23): the lattice part is the electron-screened interaction to a
nucleus, the nuclear charge fluctuation $\chi^{\ion}$, and the electron-screened interaction from the nucleus. Written with
phonons, (A5.28), each end is a screened vertex $\mathcal V=W^{\rme}$ acting on the bare nuclear-charge vertex
$\Delta_{\bq\nu}n^{\ion}$, and the middle is the phonon propagator $F$: a phonon exchanged between two screened vertices,
crystal momentum $\bq$ in at one end and out at the other.

**Include position.** Start of `a5-split`, right after: "The result is an extra, slow channel in the screened
interaction: $W^{\rme}$ to a nucleus, propagation of the nuclear charge fluctuation, and $W^{\rme}$ again from the nucleus to
the other electron. Figure A5.4 draws the result before the derivation, which follows the review and needs five steps."

**Physics.** No numerics. Row 1 is (A5.23); row 2 is (A5.28) with (A5.27), and $\mathcal V_{\bq\nu}=(\varepsilon^{\rme})^{-1}\Delta_{\bq\nu}V^{\mathrm{en}}$.

**What to draw** (`EPP.fd.draw`, as v1, two `.fd` rows).
- Row 1: [screened line 1–2, label $W$] = [screened line, label $W^{\rme}$] + [screened line 1–3, label $W^{\rme}$] • [heavy
  phonon line 3–4, label $\chi^{\ion}$] • [screened line 4–2, label $W^{\rme}$]; small labels 1, 2, 3, 4 at the points;
  figcaptions under the three glyphs: "total W", "nuclei clamped: W^e (the W of GW)", "W^{ph} = W^e χ^{ion} W^e".
- Row 2: "W^{ph} =" then the text "(1/N_ℓ) Σ_{qνν'}", then [grey disk labelled $\mathcal V_{\bq\nu}$] [phonon line labelled
  $F_{\nu\nu'}(\bq,\omega)$ with a small arrow and "q" above it] [grey disk labelled $\mathcal V_{-\bq\nu'}$]; figcaption
  "screened vertex 𝒱 = (ε^e)^{-1}Δ_{qν}V^{en} at both ends".
- Optionally a third, small glyph after row 2: "𝒱 = W^e •" (a screened line ending on a dot), captioned "a screened vertex is
  W^e acting on the bare nuclear-charge vertex Δ_{qν}n^{ion}".
- Colours: screened lines as the library draws them (double dashed), $W$ labels rose, $W^{\rme}$ labels blue; phonon line
  amber; dots and disks as v1.

**Controls** (2, optional but recommended).
1. Toggle "clamp the nuclei ($F=0$)": greys out the second term of row 1 and all of row 2, and shows "W = W^e".
2. Toggle "static screening (DFPT)": relabels the disks $\mathcal V(\omega)\to\mathcal V(0)\simeq\Delta_{\bq\nu}V^{\mathrm{KS}}$ and adds the
   note "matrix elements = g_{mnν}(k,q)".

**Initial state.** Both toggles off.

**Caption.**
- *What to notice.* The lattice part of the screened interaction is a nuclear charge fluctuation between two
  electron-screened lines; in the phonon language, a phonon exchanged between two screened vertices.
- *Marks.* Double-dashed lines: screened interactions ($W$ total, $W^{\rme}$ electrons only); heavy amber line: the nuclear
  charge fluctuation $\chi^{\ion}$ in row 1, the phonon propagator $F$ in row 2; dots: the bare nuclear-charge vertex
  $\Delta_{\bq\nu}n^{\ion}$; grey disks: screened vertices $\mathcal V$; 3 and 4 are internal points, integrated over.
- *Try.* Clamp the nuclei: only $W^{\rme}$ survives, the GW case. Switch on static screening: the vertex becomes DFPT's
  $\Delta_{\bq\nu}V^{\mathrm{KS}}$, whose matrix elements are the $g_{mn\nu}(\bk,\bq)$ of Section A.2.
- *Model.* (A5.23), (A5.27), (A5.28); exact for Coulomb-coupled electrons and nuclei, one-phonon in row 2.

**Change from v1.**
- Ids: `viz-a5-fd` and span ids `fd-a5-fd-…` (v1: `fd-a-hedin-baym-w`, `-we`, `-wph`, `-g`).
- v1 title "The lattice inside W: Eqs. (G.104), (G.105), (G.110)" carries source equation numbers: drop them.
- v1 labels: "W_e" → "W^{e}" (via `label: 'W^{e}'`), "𝒟" → "χ^{ion}" (row 1) and "F_{νν'}(q,ω)" (row 2), "g" →
  "𝒱_{qν}" and "𝒱_{−qν'}", caption "W_ph = W_e 𝒟 W_e" → "W^{ph} = W^e χ^{ion} W^e", "screened vertex g = ε_e^{-1}∇V_κ" →
  "screened vertex 𝒱 = (ε^e)^{-1}Δ_{qν}V^{en}". Aria-labels likewise.
- v1 used the same letter for the nuclear density correlator and the displacement correlator; v2 shows $\chi^{\ion}$ in
  row 1 and $F$ (mode basis) in row 2.
- Caption: rewrite in the four-part structure.

---

## 5. `fig-map.html` · `viz-a5-map` · port of v1 `viz-a-hedin-baym-map` (v1 `sections/A05-hedin-baym.html`)

**Title.** "The Hedin–Baym equations as a map" `<span class="tag">interactive</span>`.

**What the reader must learn.** How the closed set (A5.35) feeds itself: an electronic loop (Hedin's pentagon
$G\to\Lambda\to\chi^{\mathrm{irr,e}}\to W^{\rme}\to\Sigma\to G$), a nuclear loop ($W^{\rme}\to\Phi(\omega)\to F\to W^{\ph}\to\Sigma$),
and the densities and $V^{\mathrm{tot}}$ at the bottom. Electrons and nuclei are joined only through $\varepsilon^{\rme}$ and
$W^{\rme}$ (into $W^{\ph}$ and $\Phi$), through $W^{\ph}$ in $\Sigma$, and through $\langle\hat n^{\ion}\rangle$ in $V^{\mathrm{tot}}$.
Clamping the nuclei leaves Hedin's pentagon. Every bare line is the instantaneous Coulomb $v$, which Part B replaces.

**Include position.** In `a5-set`, right after the paragraph that reads (A5.35) line by line and ends: "The electronic
equations are Hedin's (1965), the nuclear ones Baym's (1961), and the review collects them as the Hedin–Baym equations.
Figure A5.5 shows how they feed each other."

**Physics.** No numerics: a directed graph. Nodes (14) and their equations:
- $G$ (`a5-eq-dysonG`, with `a5-eq-G`), $\Sigma$ (`a5-eq-Sigma`), $\Lambda$ (`a5-eq-vertex`), $\chi^{\mathrm{irr,e}}$, $W^{\rme}$,
  $\varepsilon^{\rme}$ and $\langle\hat n^{\rme}\rangle$ (all `a5-eq-We`), $W^{\ph}$ and the screened vertex (`a5-eq-Wph`),
  $F$ (`a5-eq-Finv`, defined in `a5-eq-F`), $\Phi(\omega)$ (`a5-eq-Phiw`), $\langle\hat n^{\ion}\rangle$ (`a5-eq-nion`),
  $V^{\mathrm{tot}}$ (`a5-eq-Vtot`), $v$ (`a5-eq-v`), and the bare potential change $\Delta_{\bq\nu}V^{\mathrm{en}}$ (defined under
  `a5-eq-dn` and `a5-eq-Wph`).
- Directed edges "ingredient → quantity" (v1 edges renamed): $G\to\langle\hat n^{\rme}\rangle$; $F\to\langle\hat n^{\ion}\rangle$;
  $\langle\hat n^{\rme}\rangle,\langle\hat n^{\ion}\rangle,v\to V^{\mathrm{tot}}$; $V^{\mathrm{tot}},\Sigma\to G$;
  $G,\Lambda,W^{\rme},W^{\ph}\to\Sigma$; $G,\Lambda\to\chi^{\mathrm{irr,e}}$; $v,\chi^{\mathrm{irr,e}}\to W^{\rme}$;
  $v,\chi^{\mathrm{irr,e}}\to\varepsilon^{\rme}$; $\Sigma,G\to\Lambda$; $\varepsilon^{\rme},F,\Delta_{\bq\nu}V^{\mathrm{en}}\to W^{\ph}$;
  $W^{\rme}\to\Phi(\omega)$; $\Phi(\omega)\to F$.

**What to draw.** As v1: buttons positioned in a 680×440 frame (phones: 680/640), electronic pentagon left, nuclear loop
right, densities and $V^{\mathrm{tot}}$ at the bottom; SVG arrows; selected node and its edges dark; outline colours blue
(electronic), amber (nuclear), rose ($v$), ink ($V^{\mathrm{tot}}$). Node labels as HTML: $\Lambda$, $\Sigma$, $G$,
$W^{\ph}$, $W^{\rme}$, $F$, $\chi^{\mathrm{irr,e}}$, $\varepsilon^{\rme}$, $\Phi(\omega)$, $\langle n^{\rme}\rangle$, $v$, $V^{\mathrm{tot}}$,
$\langle n^{\ion}\rangle$, $\Delta_{\bq\nu}V^{\mathrm{en}}$ (upright superscripts). Panel below: the selected node's name,
its equation in v2 notation (copy from the section) with a link `<a class="eqref" href="#…"></a>`, one sentence of
meaning, and the "Built from: … · Feeds: …" line.

**Controls** (4).
1. Toggle "clamp the nuclei ($F=0$): Hedin": greys $F$, $\Phi(\omega)$, $W^{\ph}$, $\Delta_{\bq\nu}V^{\mathrm{en}}$; note: "With
   $F=0$, $W^{\ph}=0$ and $\langle\hat n^{\ion}\rangle=n^{\ion}_0$, a rigid lattice of point charges; what remains is Hedin's
   pentagon."
2. Toggle "GW vertex, $\Lambda=\delta\delta$": cuts the arrows into $\Lambda$; note: "$\chi^{\mathrm{irr,e}}=-2\ii\hbar GG$ becomes the
   $\chi^0$ of Section A.3 with Kohn–Sham propagators, and $\Sigma=\ii\hbar G(W^{\rme}+W^{\ph})$: GW plus the Fan–Migdal
   self-energy of Section A.7."
3. Toggle "adiabatic (static screening)" (new): marks $\Phi(\omega)$ and $W^{\ph}$; note: "With $W^{\rme}(\omega)\to W^{\rme}(0)$,
   $\Phi(0)$ is the Born–Oppenheimer force-constant matrix of Section A.1, which DFPT computes, and $\mathcal V(0)$ is DFPT's
   $\Delta_{\bq\nu}V^{\mathrm{KS}}$."
4. Toggle "show the instantaneous Coulomb $v(12)$": rose rings on nodes that contain $v$ explicitly ($V^{\mathrm{tot}}$,
   $W^{\rme}$, $\varepsilon^{\rme}$, $\Delta_{\bq\nu}V^{\mathrm{en}}$), rose fill on nodes that contain it through $W^{\rme}$ ($\Sigma$,
   $W^{\ph}$, $\Phi$), dashed rose arrows from $v$; note: "Part B keeps $v$ as the longitudinal photon propagator
   $D^{0\mathrm L}$ (`b5-eq-D0L`) and adds the transverse photons (Sections B.5–B.6)."

**Initial state.** $W^{\ph}$ selected; all toggles off.

**Caption.**
- *What to notice.* Electrons and nuclei meet only through $\varepsilon^{\rme}$ and $W^{\rme}$ (building $W^{\ph}$ and $\Phi$), through
  $W^{\ph}$ in $\Sigma$, and through $\langle n^{\ion}\rangle$ in $V^{\mathrm{tot}}$; clamping the nuclei leaves Hedin's pentagon.
- *Marks.* Node outlines: blue electronic, amber nuclear, rose the bare $v$, ink $V^{\mathrm{tot}}$; arrows from ingredients
  to the quantity built from them; the selected node's arrows dark.
- *Try.* Click $\Sigma$ and follow its inputs; clamp the nuclei; switch on the Coulomb view to see every slot that Part B
  fills with the photon propagator.
- *Model.* The closed set (A5.35): each node is one line of it plus the two inputs $v$ and $\Delta_{\bq\nu}V^{\mathrm{en}}$;
  harmonic, zero temperature.

**Change from v1.**
- Ids and classes: `viz-a5-map`, `a5-map-…` (v1: `hb-map`, `hb-node`, `hb-edges`, `hb-panel`, `hb-io`).
- Node renames: Γ → Λ (Γ is the zone centre); 𝒟 → F; $P_{\rm e}$ → $\chi^{\mathrm{irr,e}}$; $\epsilon_{\rm e}$ → $\varepsilon^{\rme}$;
  $W_{\rm e}$, $W_{\rm ph}$ → $W^{\rme}$, $W^{\ph}$; $\Pi^{\rm ph}$ (the review's force-constant-unit self-energy) → $\Phi(\omega)$ (in
  v2, $\Pi^{\ph}$ is A.6's mode-basis self-energy); $\langle n_{\rm n}\rangle$ → $\langle n^{\ion}\rangle$; $V_{\rm tot}$ →
  $V^{\mathrm{tot}}$; $\nabla V_\kappa$ → $\Delta_{\bq\nu}V^{\mathrm{en}}$. Update the `NAME` map (HTML labels) and the aria-labels.
- Panel texts: replace every "G Eqs. …" by eqrefs to the section's equations, rewrite every equation in v2 notation
  (spin factor 2 explicit, $\Lambda$, $V^{\mathrm{ext}}$ not $\varphi$, $\chi^{\mathrm{irr,e}}$), and drop "E.374", "E.375", "E.437"
  and the $D^{0l}_{\mu\nu}$ formula with Cartesian $\mu,\nu$.
- Add the "adiabatic" toggle and its note.
- Caption: v1's Marks-first paragraph → the four-part structure.
