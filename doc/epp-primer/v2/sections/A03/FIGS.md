# A.3 figure specs (section `a-dfpt`, prefix `a3`)

Three figures: two ports from v1 and one new geometric figure. Read `v2/build/BUILDER_GUIDE.md` §5 and
`v2/NOTATION.md` first. The section text that surrounds each figure is in `v2/sections/A03-dfpt.html`.

Common to all three:

- **Notation (v2, never the sources').** $\Delta V^{\mathrm{KS}}$, $\Delta V^{\mathrm{en}}$, $\Delta n$, $\chi^0$,
  $\varepsilon^{\mathrm e}$ (test electron), $\varepsilon^{\mathrm e}_{\mathrm{RPA}}$, $v(\bq)=4\pi e^2/q^2$, $f^{\mathrm{xc}}$,
  $k_{\mathrm F}$, $k_{\mathrm{TF}}$, $E_{\mathrm F}$ (not $\varepsilon_{\mathrm F}$), $N_{\mathrm F}=\partial n_0/\partial\mu$, $n_0$, $r_s$ (dimensionless, as in A.2), $a_0$, $R_{\mathrm{WS}}$ (A.2's Wigner–Seitz radius), $n^{\mathrm{ind}}$, $V^{\mathrm{scr}}_\kappa$,
  Lindhard function $F_{\mathrm L}(x)$ (as in Section A.2; $F$ alone is the phonon propagator), mixing weight $w$ (not $\beta$: $\beta$
  is a Cartesian index), error factor $\zeta$, cycle count $N_{\mathrm{cyc}}$, projector $\mathcal P^{\mathrm{occ}}_{\bk+\bq}$,
  $\Delta\psi^{\perp}_{n\bk}$, $E_{\mathrm{shift}}$. Canvas text in `EPP.drawRich` form: `ΔV^{KS}`, `ε^e`, `k_{TF}`,
  `F_L`.
  Dielectric functions are `ε` (varepsilon glyph), never `ϵ`.
- **Units.** Gaussian, $e^2=14.3996$ eV Å, $a_0=0.529177$ Å, Ry $=13.6057$ eV, $\hbar^2/2m_{\mathrm e}=3.80998$ eV Å²,
  $\hbar/m_{\mathrm e}=1.15768\times10^{-4}$ m²/s.
- **Electron gas (jellium) formulas** used by figures 1 and 3 ($r_s$ dimensionless, as in Section A.2: $\frac{4\pi}3(r_sa_0)^3=1/n_0$):
  - $k_{\mathrm F}=1.91916/(r_sa_0)$ Å⁻¹; $k_{\mathrm{TF}}^2=4k_{\mathrm F}/(\pi a_0)$, i.e. $(k_{\mathrm{TF}}/k_{\mathrm F})^2=0.66341\,r_s$;
    $E_{\mathrm F}=(1.91916/r_s)^2\times13.6057$ eV; $n_0=3/(4\pi r_s^3a_0^3)$; $R_{\mathrm{WS}}=Z^{1/3}r_sa_0$.
  - Lindhard function (A3.18): $F_{\mathrm L}(x)=\tfrac12+\frac{1-x^2}{4x}\ln\left|\frac{1+x}{1-x}\right|$, with guards
    $F_{\mathrm L}\simeq1-x^2/3$ for $x<10^{-4}$ and $F_{\mathrm L}(1)=\tfrac12$.
  - $v\chi^0(q)=-(k_{\mathrm{TF}}/q)^2\,F_{\mathrm L}(q/2k_{\mathrm F})$ (A3.18, A3.5).
  - Exchange-only kernel of the section (text after A3.6): $f^{\mathrm{xc}}=-\pi e^2/k_{\mathrm F}^2$, so
    $(v+f^{\mathrm{xc}})\chi^0=v\chi^0\,(1-q^2/4k_{\mathrm F}^2)$. Call it "exchange kernel (toy)"; it is the test-electron case.
  - Dielectric functions: Thomas–Fermi $1+k_{\mathrm{TF}}^2/q^2$ (A3.1); RPA $1+(k_{\mathrm{TF}}/q)^2F_{\mathrm L}$ (A3.19);
    test electron $1+(k_{\mathrm{TF}}/q)^2F_{\mathrm L}\,(1-q^2/4k_{\mathrm F}^2)$ (A3.13 with the kernel above).
- Link equations in captions with `<a class="eqref" href="#a3-eq-…"></a>`; ids used below: `a3-eq-tf` (A3.1), `a3-eq-den` (A3.4), `a3-eq-dh` (A3.5), `a3-eq-dxc` (A3.6), `a3-eq-chi0` (A3.8),
  `a3-eq-dn` (A3.10), `a3-eq-sternheimer` (A3.11), `a3-eq-loop` (A3.12), `a3-eq-scf` (A3.13), `a3-eq-series` (A3.14),
  `a3-eq-mix` (A3.15), `a3-eq-tc` (A3.16), `a3-eq-eps` (A3.17), `a3-eq-lindhard` (A3.18), `a3-eq-epsj` (A3.19).
- Literature numbers ($r_s$ of Al, Na, K) are flagged "(literature values, illustrative)".

---

## 1. `fig-cloud.html` · id `viz-a3-cloud`

**Kind.** Port of v1 `viz-a-dfpt-cloud`, file `doc/epp-primer/sections/A03-dfpt.html`: the `<figure … id="viz-a-dfpt-cloud">`
element (v1 lines 121–127), the derivation box above it (v1 lines 113–119, the formula for $S(r)$), and the block
"Figure 2" of the section-wide `<script>` (v1 lines 281–380, functions `F`, `TF2`, `kchi`, `Sprof`, `interp`, `draw2`).

**What the reader must learn.**
1. The bare perturbation of a displaced ion is a long-range dipole field ($\propto\cos\theta/r^2$); the electrons'
   induced charge cancels most of it within a few screening lengths $1/k_{\mathrm{TF}}$ (0.5–0.75 Å; 26 % of the ion's charge
   lies within $1/k_{\mathrm{TF}}$, 90 % within $3.9/k_{\mathrm{TF}}$), so only 8–14 % of the bare field survives at the
   Wigner–Seitz radius (RPA).
2. The cloud moves with the ion and holds exactly the ion's charge (perfect screening); by Gauss's theorem the surviving
   fraction $S(r)$ equals the fraction of the ion's charge not yet neutralized within $r$ (test-charge screening).
3. The full free-electron response (Lindhard, no exchange–correlation) differs from Thomas–Fermi at short distances and leaves Friedel ripples of
   period $\pi/k_{\mathrm F}$; the test-electron (exchange) screening differs again at $q\gtrsim k_{\mathrm F}$.

**Include position.** End of subsection `a3-tf` ("The simplest model: one ion in an electron gas"), right after:
"By Gauss's theorem $S(r)$ is also the fraction of the ion's charge that is not yet neutralized within the radius $r$.
[Figure] draws the displaced ion, its potential and its cloud, and compares the Thomas–Fermi model with the full response
of free electrons, the Lindhard function derived later in this section, without and with exchange."

**Physics.**
- A point ion of valence $Z$ (1, 2 or 3) at the origin of an electron gas of density parameter $r_s$, displaced
  by $\bu=u\hat{\mathbf z}$. Everything is per unit displacement $u$ (eV/Å for potentials, Å⁻⁴ for $\Delta n$).
- Bare ion: $V_\kappa(r)=-Ze^2/r$, $V_\kappa(q)=-4\pi Ze^2/q^2$. Screened: $V^{\mathrm{scr}}_\kappa(q)=V_\kappa(q)/\varepsilon(q)$
  with one of the three $\varepsilon$ above (A3.1, A3.13, A3.19).
- Surviving fraction of the field along the displacement axis (section text after A3.1; v1 derivation box):
  $S(r)=\partial_rV^{\mathrm{scr}}_\kappa/\partial_rV_\kappa=1+\frac2\pi\int_0^\infty\dd q\,\big[\tfrac1{\varepsilon(q)}-1\big]\big[\tfrac{\sin qr}{q}-r\cos qr\big]$.
  Thomas–Fermi in closed form: $S=(1+k_{\mathrm{TF}}r)\,\ee^{-k_{\mathrm{TF}}r}$. Lindhard: numerical (v1: trapezoid,
  $q\le12k_{\mathrm F}$, 3000 points, rotation recurrence for sin/cos). Exchange case: $1/\varepsilon-1\to+k_{\mathrm{TF}}^2/3q^2$
  at large $q$; subtract $c/(q^2+\kappa_0^2)$ with $c=k_{\mathrm{TF}}^2/3$, $\kappa_0=2k_{\mathrm F}$ and add its integral
  analytically, $c\,\frac{\pi}{2\kappa_0^2}\big[1-\ee^{-\kappa_0r}(1+\kappa_0r)\big]$ (as v1). Name this constant
  `kap0` or similar in code, not `kappa` alone (κ is the atom index).
- Screened perturbation on the map: $\Delta V^{\mathrm{KS}}(\br)=-\partial_zV^{\mathrm{scr}}_\kappa=-(Ze^2/r^2)\,S(r)\cos\theta$;
  bare: $\Delta V^{\mathrm{en}}(\br)=-(Ze^2/r^2)\cos\theta$ ($\theta$ measured from $+z$; the potential energy is lowered
  ahead of the moving ion).
- Induced electron density of the static ion: in $q$ space $n^{\mathrm{ind}}(q)=Z\,[1-1/\varepsilon(q)]$ for Thomas–Fermi
  and RPA, and $n^{\mathrm{ind}}(q)=Z\,a/\big(1+a\,(1-q^2/4k_{\mathrm F}^2)\big)$ with $a=(k_{\mathrm{TF}}/q)^2F_{\mathrm L}$ for
  the exchange case ($\Delta n=\chi^0\Delta V^{\mathrm{KS}}$, A3.8 and A3.13). Real space:
  $n^{\mathrm{ind}}(r)=\frac{1}{2\pi^2r}\int_0^\infty\dd q\,q\sin(qr)\,n^{\mathrm{ind}}(q)$; Thomas–Fermi closed form
  $Zk_{\mathrm{TF}}^2\ee^{-k_{\mathrm{TF}}r}/4\pi r$; $\int4\pi r^2n^{\mathrm{ind}}\,\dd r=Z$. For the displaced ion
  $\Delta n(\br)=-\partial_zn^{\mathrm{ind}}=-n^{\mathrm{ind}}'(r)\cos\theta$. Precompute radial arrays when $r_s$
  or the model changes (cache as v1 does); the map is then an interpolation, well under 15 ms.
- Checked numbers (RPA = Lindhard without xc; "exch" = test electron with the toy kernel):

  | metal | $r_s$, $Z$ | $k_{\mathrm F}$ (Å⁻¹) | $k_{\mathrm{TF}}$ (Å⁻¹) | $1/k_{\mathrm{TF}}$ (Å) | $E_{\mathrm F}$ (eV) | $\pi/k_{\mathrm F}$ (Å) | $R_{\mathrm{WS}}$ (Å) | $S(R_{\mathrm{WS}})$ TF / RPA / exch | bare field at $R_{\mathrm{WS}}$ → RPA, exch (eV/Å) |
  |---|---|---|---|---|---|---|---|---|---|
  | Al | 2.07, 3 | 1.752 | 2.053 | 0.487 | 11.70 | 1.79 | 1.58 | 0.166 / 0.121 / 0.116 | 17.3 → 2.10, 2.02 |
  | Na | 3.93, 1 | 0.923 | 1.490 | 0.671 | 3.24 | 3.40 | 2.08 | 0.185 / 0.139 / 0.284 | 3.33 → 0.46, 0.95 |
  | K  | 4.86, 1 | 0.746 | 1.340 | 0.746 | 2.12 | 4.21 | 2.57 | 0.142 / 0.076 / 0.215 | 2.18 → 0.17, 0.47 |

  Friedel ripples of $S(r)$ for Na (RPA) reach about +0.06 near 4 Å and −0.035 near 5.8 Å.

**What to draw.**
- **(a) "The displaced ion and its cloud"** (geometric, square, `EPP.canvas` or a Plot heatmap): the $x$–$z$ plane,
  $-6\le x,z\le6$ Å, $z$ up. The ion as a filled disc at the origin (phonon colour) with an arrow along $+z$ labelled
  **u** (`drawRich` with `{bold:true}`). A diverging colour map of the field chosen by the "map" control: screened
  $\Delta V^{\mathrm{KS}}$, bare $\Delta V^{\mathrm{en}}$, or induced $\Delta n$. Use one fixed compressed signed scale for
  the two potentials, e.g. $\mathrm{sign}(V)\,\mathrm{asinh}(|V|/V_0)/\mathrm{asinh}(V_{\max}/V_0)$ with $V_0=0.05$ eV/Å,
  $V_{\max}=15$ eV/Å, so that switching bare ↔ screened shows the contrast on the same colours; a colour bar with units
  (eV/Å; Å⁻⁴ for $\Delta n$). Clip $r<0.15$ Å. Dashed circles at $r=1/k_{\mathrm{TF}}$ and $r=R_{\mathrm{WS}}$, labelled
  `1/k_{TF}` and `R_{WS}`. Axis labels `x (Å)`, `z (Å)`.
- **(b) "Field that survives along the displacement axis"**: $S(r)$ for $0\le r\le10$ Å, fixed limits $[-0.1,1.05]$ (or
  a zoom toggle to $[-0.08,0.08]$ if a fifth control is free). Curves: Thomas–Fermi (muted, dotted), Lindhard RPA
  (electron colour), Lindhard + exchange (`c5`, dashed); the selected model thicker. Dashed line $S=1$ "bare ion";
  vertical lines at $R_{\mathrm{WS}}$ and $1/k_{\mathrm{TF}}$. Axis labels `r (Å)`, `S(r)`.
- **(c) "Where the screening electrons sit"**: the radial distribution $4\pi r^2n^{\mathrm{ind}}(r)/Z$ (Å⁻¹) of the static
  ion's cloud for the three models, $0\le r\le6$ Å, with a text mark "area = 1: the cloud holds the ion's charge".
- Panels (a) on top (full width on phones), (b) and (c) side by side below (`viz-grid two`, stacked on phones).

**Controls (5).**
1. `metal`: seg Al / Na / K / custom (sets $r_s$ and $Z$ as in the table; moving a slider switches to custom).
2. `rs`: range 1.5–6 (dimensionless), step 0.01.
3. `Z`: seg 1 / 2 / 3.
4. `model`: seg "Thomas–Fermi" / "Lindhard (RPA)" / "Lindhard + exchange".
5. `map`: seg "screened ΔV^{KS}" / "bare ΔV^{en}" / "induced Δn".

**Readouts.** $k_{\mathrm F}$, $k_{\mathrm{TF}}$ and $1/k_{\mathrm{TF}}$, $E_{\mathrm F}$, $\pi/k_{\mathrm F}$, $R_{\mathrm{WS}}$,
and "field at $R_{\mathrm{WS}}$: bare … eV/Å, screened … eV/Å (… %)" for the selected model.

**Initial state.** Na ($r_s=3.93$, $Z=1$), model Lindhard (RPA), map "screened". The map shows a compact
dipole confined within about 1 Å; the readout says 0.46 of 3.33 eV/Å (14 %) at $R_{\mathrm{WS}}=2.08$ Å.

**Caption (≈200–260 words).**
- *What to notice.* The bare field of the displaced ion (map → bare) fills the whole panel; the screened one fades within a
  few screening lengths $1/k_{\mathrm{TF}}$, and at $R_{\mathrm{WS}}$ only 8–14 % of it survives in the three presets (RPA). Beyond, only
  weak Friedel ripples of period $\pi/k_{\mathrm F}$ remain.
- *Marks.* Colours (sign of $\Delta V$ or $\Delta n$, compressed scale), the arrow **u**, the two circles, the three curve
  styles, the $S=1$ line.
- *Try.* (1) Switch the map between bare and screened. (2) Choose Thomas–Fermi: the ripples vanish, because the
  local rule has no $2k_{\mathrm F}$ kink. (3) Choose "Lindhard + exchange": in Na and K the field that survives at
  $R_{\mathrm{WS}}$ doubles to triples, in Al it barely changes; the difference comes from $q\gtrsim k_{\mathrm F}$, where the
  toy kernel is least reliable. (4) Raise $r_s$: the screening length $1/k_{\mathrm{TF}}\propto r_s^{1/2}$ grows.
- *Model.* Point ion in jellium, static screening diagonal in $q$: (A3.1), (A3.18), (A3.19), the test-electron
  function (A3.13) with the exchange kernel of the text after (A3.6); $S(r)$ and Gauss's theorem as in the text; $r_s$
  of Al, Na, K (literature values, illustrative).

**What in v1 must change.**
- v1 labels and readouts use source notation: `∂v/∂τ`, `ε_F`, "test charge: Lindhard RPA (G.55)", "test electron: toy xc
  (G.54)", `ε(q)` with the `(G.xx)` tags, `e²/4πε0` comments. Use $\Delta V^{\mathrm{KS}}$, $E_{\mathrm F}$, the v2 model names
  and eqrefs; Gaussian $e^2$ (same number, 14.3996 eV Å).
- v1's left panel ($1/\varepsilon(q)$ against $q/k_{\mathrm F}$) is replaced by the geometric map (a) and the induced-charge
  profile (c), as the brief asks ("geometric, with radial profiles, the bare and screened potential, and the induced
  charge").
- The v1 code lived in the section-wide script; move it into the figure file, rename classes `a3-cloud-…`.
- v1 caption order (Model first) → v2 structure (What to notice, Marks, Try, Model).
- The v1 control "S(r) axis: full / zoom on tail" may stay only if it fits in 5 controls; otherwise use the fixed range
  above.

---

## 2. `fig-sternheimer.html` · id `viz-a3-sternheimer`

**Kind.** New.

**What the reader must learn.**
1. The first-order change $\Delta\psi_{n\bk}$ of an occupied state has components on every state at $\bk+\bq$ (the
   wavevector of the perturbation is added: a band diagram with arrows from $\bk$ to $\bk+\bq$).
2. The components on occupied states can be huge (they diverge when an occupied state at $\bk+\bq$ has the energy
   $E_{n\bk}$), but they cancel in the density; only the empty-band part $\Delta\psi^{\perp}_{n\bk}$ matters (A3.10), and it
   is bounded by the gap.
3. The Sternheimer operator $H^{\mathrm{KS}}-E_{n\bk}$ is singular or badly conditioned exactly along those occupied
   directions; adding $E_{\mathrm{shift}}\mathcal P^{\mathrm{occ}}_{\bk+\bq}$ lifts them above zero without changing
   $\Delta\psi^{\perp}_{n\bk}$ (A3.11).

**Include position.** Subsection `a3-stern` ("The Sternheimer equation: occupied states only"), after the opening
paragraph, which ends: "[Figure] shows, for a model insulator, the components of $\Delta\psi_{n\bk}$ on every state at
$\bk+\bq$, which of them the density needs, and why the equation that determines them can be made well behaved."
The derivation steps (A3.10)–(A3.11) follow the figure.

**Physics (toy model, label it as such).**
- One-dimensional chain, lattice constant $a=3$ Å, three sites per cell, A, B, C at $x=0$, $a/3$, $2a/3$; on-site energies
  $E_{\mathrm A}$ (control, default −1.5 eV), $E_{\mathrm B}=0.5$ eV, $E_{\mathrm C}=1.5$ eV; nearest-neighbour hopping
  $t=2$ eV (A–B and B–C in the cell, C to A of the next cell). Bloch Hamiltonian with the phase at the cell corner:
  `H(k) = [[E_A, −t, −t e^{−ika}], [−t, E_B, −t], [−t e^{ika}, −t, E_C]]`; bands $E_1<E_2<E_3$ and eigenvectors with
  `EPP.eigHerm`. Band 1 is occupied (two electrons per cell with spin), bands 2 and 3 are empty. With the defaults the
  occupied band is 1.22 eV wide and the gap is 1.88 eV ($E_{\mathrm A}=-3$ eV: 0.88 and 3.10 eV; $E_{\mathrm A}=-0.5$ eV:
  1.48 and 1.15 eV).
- Perturbation: the displacement wave changes the on-site energy of atom A of cell $\ell$ by $\delta_{\mathrm A}\ee^{\ii qR_\ell}$,
  $\delta_{\mathrm A}=1$ eV/Å per unit amplitude (a deformation-potential toy for $\Delta V^{\mathrm{KS}}$). Its matrix element
  between Bloch states is $\langle\psi_{m,k+q}|\Delta V^{\mathrm{KS}}|\psi_{1,k}\rangle=\delta_{\mathrm A}\,c^*_{\mathrm A,m}(k+q)\,c_{\mathrm A,1}(k)$,
  with $c_{\mathrm A,m}$ the A component of the eigenvector of band $m$ (only moduli are shown, so the gauge does not matter).
- Components of the first-order change of the occupied state (step 1 of the derivation of A3.8, and A3.10):
  $C_m=\langle\psi_{m,k+q}|\Delta V^{\mathrm{KS}}|\psi_{1,k}\rangle/\big(E_1(k)-E_m(k+q)\big)$, $m=1,2,3$. $C_1$ is the occupied
  component (cancels in $\Delta n$; diverges when $E_1(k+q)=E_1(k)$, e.g. $q=0$ or $k+q=-k$); $|\Delta\psi^{\perp}|^2=|C_2|^2+|C_3|^2$.
- Spectrum of the shifted Sternheimer operator at $k+q$ (A3.11): $\lambda_1=E_1(k+q)-E_1(k)+E_{\mathrm{shift}}$ (occupied
  subspace), $\lambda_{2,3}=E_{2,3}(k+q)-E_1(k)$ (empty subspace); condition number $=\max|\lambda|/\min|\lambda|$ (say
  "condition number" in words; no symbol, κ is the atom index).
- Check that the figure may display: solving $(H(k+q)-E_1(k)+E_{\mathrm{shift}}\mathcal P^{\mathrm{occ}})x=-(1-\mathcal P^{\mathrm{occ}})\Delta V\psi_1$
  in the 3-dimensional space gives $|x|^2=|C_2|^2+|C_3|^2$ for every $E_{\mathrm{shift}}$ with $\lambda_1\ne0$.
- Verified numbers (defaults, $k=-0.35\,\pi/a$, $q=0.6\,\pi/a$): $E_1(k)=-3.873$ eV, $E(k+q)=(-3.985,1.057,3.428)$ eV,
  $|C|^2=(24.5,\ 0.0082,\ 0.0010)$ Å⁻², $|\Delta\psi^\perp|^2=0.0092$ Å⁻²; $\lambda=(-0.113,4.93,7.30)$ eV and condition number 65 at
  $E_{\mathrm{shift}}=0$; 5.3 at $E_{\mathrm{shift}}=1.5$ eV; $|x|^2=0.00922$ for $E_{\mathrm{shift}}=0.3$, 1 and 2 eV. At
  $k=-0.3\,\pi/a$ the occupied component is infinite and the operator singular.
  (Units of $|C|^2$: per Å² of displacement amplitude, since $\delta_{\mathrm A}$ is per Å.)

**What to draw.**
- **(a) Band diagram**: $E$ (eV, fixed range about −5 to 5) against $k$ in units of $\pi/a$ on $[-1,1]$. Band 1 thick in the
  electron colour with a light fill below it ("occupied"); bands 2 and 3 thin ("empty"). A filled dot for the occupied state
  at $k$ on band 1; a dashed vertical line at $k+q$ (folded into the zone); hollow dots on the three bands at $k+q$. Arrows
  from the state at $k$ to each state at $k+q$, width or opacity growing with $|C_m|$ (log-compressed); the arrow to band 1
  grey and dashed with the label "occupied: cancels in Δn", the arrows to bands 2 and 3 solid with the label "Δψ^⊥". A
  horizontal double-headed arrow under the axis from $k$ to $k+q$ labelled **q** (bold). Optionally the values $|C_m|^2$
  next to the arrowheads.
- **(b) Spectrum of the Sternheimer operator** at fixed $q$, against the occupied state's $k$ over the zone: $\lambda_1(k)$
  (depends on $E_{\mathrm{shift}}$; electron colour), $\lambda_2(k)$, $\lambda_3(k)$ (muted), the zero line, a shaded band
  $|\lambda|<0.1$ eV labelled "nearly singular", and a marker at the current $k$. With $E_{\mathrm{shift}}=0$ the curve
  $\lambda_1$ crosses zero at the $k$ where $E_1(k+q)=E_1(k)$; with $E_{\mathrm{shift}}$ above the occupied bandwidth it lies
  above zero everywhere. Axis labels `k (π/a)`, `eigenvalue (eV)`.

**Controls (4).**
1. `k`: occupied state's wavevector, range −1 to 1 ($\pi/a$), step 0.005.
2. `q`: perturbation wavevector, range −1 to 1 ($\pi/a$), step 0.005.
3. `Eshift`: $E_{\mathrm{shift}}$, range 0–3 eV, step 0.01.
4. `EA`: on-site energy of atom A, range −3 to −0.5 eV (sets the gap and the occupied bandwidth).

**Readouts.** $E_1(k)$; $E_m(k+q)$ for $m=1,2,3$; $|C_1|^2$ ("occupied, cancels"); $|\Delta\psi^{\perp}|^2$ ("needed");
smallest $|\lambda|$; condition number; occupied bandwidth and gap; a note "E_shift exceeds the occupied bandwidth" when true.

**Initial state.** $E_{\mathrm A}=-1.5$ eV, $k=-0.35\,\pi/a$, $q=0.6\,\pi/a$, $E_{\mathrm{shift}}=0$: the occupied arrow is
thick (|C₁|² ≈ 25), the empty ones thin (≈ 0.01), $\lambda_1$ is close to zero and the condition number is about 65.

**Caption (≈200–260 words).**
- *What to notice.* The occupied component can be arbitrarily large, and it diverges when the occupied band at $k+q$ has
  the energy of the state at $k$; yet it never enters the density. The empty components, which the density needs, stay
  small because they are divided by at least the gap.
- *Marks.* Bands, fill, dots, the arrow styles, **q**, the eigenvalue curves and the "nearly singular" band.
- *Try.* (1) Move $k$ until $k+q=-k$ (or set $q=0$): $|C_1|$ and the condition number blow up, $|\Delta\psi^{\perp}|^2$ does
  not change. (2) Raise $E_{\mathrm{shift}}$ above the occupied bandwidth: the $\lambda_1$ curve lifts above zero, the condition
  number drops to a few, and $|\Delta\psi^{\perp}|^2$ is unchanged. (3) Make $E_{\mathrm A}$ more negative: the gap widens and the
  empty components shrink, while $|C_1|^2$ grows as the occupied band flattens (as built).
- *Model.* Toy three-band chain (parameters above, not a real material), on-site deformation of atom A as the
  perturbation; equations (A3.10) and (A3.11).

---

## 3. `fig-scf.html` · id `viz-a3-scf`

**Kind.** Port of v1 `viz-a-dfpt-scf`, file `doc/epp-primer/sections/A03-dfpt.html`: the `<figure … id="viz-a-dfpt-scf">`
element (v1 lines 80–86), the bubble-chain diagrams in the section body (v1 lines 62–72, `div.fd` with
`a-dfpt-fd-dressed`, `…-bare`, `…-one`, `…-two`), and the blocks "diagrams" and "Figure 1" of the section-wide `<script>`
(v1 lines 145–279).

**What the reader must learn.**
1. The DFPT loop sums the geometric series (A3.14): the terms alternate in sign and grow when $|(v+f^{\mathrm{xc}})\chi^0|>1$,
   which in a metal is every long wavelength, although the sum $1/\varepsilon^{\mathrm e}$ is finite.
2. Linear mixing (A3.15) multiplies the error by $\zeta=1-w\varepsilon^{\mathrm e}$ per cycle: it converges only for
   $0<w\varepsilon^{\mathrm e}<2$, so any fixed $w$ fails at small enough $q$ (charge sloshing).
3. Kerker's weight $w(q)=w_0q^2/(q^2+q_0^2)$ keeps $w\varepsilon^{\mathrm e}$ within a factor 2 of $w_0$ at every $q$.

**Include position.** End of subsection `a3-mixing` ("Mixing and charge sloshing"), after: "Production codes use more
elaborate schemes that learn the screening from previous cycles, but the arithmetic of (A3.15) is the reason they are
needed: even in silicon, with $\varepsilon^{\mathrm e}\approx12$ for slowly varying potentials, the plain iteration diverges.
[Figure] runs the loop for one Fourier component of the electron gas."

**Physics.**
- One Fourier component of the electron gas, no local fields: $x(q)\equiv(v+f^{\mathrm{xc}})\chi^0=-(k_{\mathrm{TF}}/q)^2F_{\mathrm L}(q/2k_{\mathrm F})$
  (RPA), multiplied by $(1-q^2/4k_{\mathrm F}^2)$ when the exchange toggle is on. $\varepsilon^{\mathrm e}=1-x$ (A3.13, A3.19).
- Weight: $w(q)=w$ (plain) or $w_0q^2/(q^2+q_0^2)$ with $q_0=k_{\mathrm{TF}}$ and $w_0=w$ (Kerker). $\zeta=1-w(q)\,\varepsilon^{\mathrm e}$.
- Iterates from $\Delta V_0=\Delta V^{\mathrm{en}}$ (step 3 of the mixing derivation):
  $\Delta V/\Delta V^{\mathrm{en}}=1/\varepsilon^{\mathrm e}+\zeta^{N_{\mathrm{cyc}}}(1-1/\varepsilon^{\mathrm e})$, $N_{\mathrm{cyc}}=0\dots40$.
- Series terms for the diagram strip (A3.14): $x^0=1$, $x^1$, $x^2$ at the current $q$ (signed).
- Cycles to $10^{-6}$: smallest $N_{\mathrm{cyc}}$ with $|\zeta|^{N}\,|1-1/\varepsilon^{\mathrm e}|<10^{-6}/\varepsilon^{\mathrm e}$, or
  "never (diverges)" for $|\zeta|\ge1$ (v1 counted relative to $|\varepsilon-1|$; either is fine if the caption says which).
  Stability limit: $w<2/\varepsilon^{\mathrm e}$ (plain) or $w_0<2(q^2+q_0^2)/(q^2\varepsilon^{\mathrm e})$ (Kerker).
- Checked default ($r_s=3.93$, $q=0.5k_{\mathrm F}$, RPA): $(k_{\mathrm{TF}}/q)^2=10.43$, $F_{\mathrm L}(0.25)=0.979$,
  $\varepsilon^{\mathrm e}=11.21$, $1/\varepsilon^{\mathrm e}=0.0892$; with $w=0.15$, $\zeta=-0.681$ (converges, alternating);
  stable for $w<0.178$. With the exchange kernel $\varepsilon^{\mathrm e}=10.57$.

**What to draw.**
- **(a) The series as diagrams** (top strip, `EPP.fd.draw` inside `div.fd`, as in v1): screened vertex = bare vertex +
  bare vertex with one bubble and one interaction line + two bubbles + ⋯. Labels in v2: the dressed vertex
  `ΔV^{KS}`, the bare vertex `ΔV^{en}`, bubbles `χ^0`, the interaction line `v+f^{xc}` (dashed "coulomb" glyph), the phonon
  line `qν`. Under each diagram its value for the current $q$: `1`, `x = …`, `x^2 = …` (signed, 3 significant digits), so
  the reader sees the alternating, growing terms.
- **(b) The iterates**: $|\Delta V/\Delta V^{\mathrm{en}}|$ against $N_{\mathrm{cyc}}$ (0–40) on a log axis ($10^{-4}$ to $10^{2}$);
  filled points for positive values, hollow for negative; a dashed horizontal line at $1/\varepsilon^{\mathrm e}$ labelled
  `1/ε^e = …` (A3.13); a text "converges: |ζ| = …" or "diverges: |ζ| = …". Axis labels `cycle N_{cyc}`,
  `|ΔV / ΔV^{en}|`.
- **(c) The error factor**: $\zeta(q)$ against $q/k_{\mathrm F}$ on $[0,3]$, fixed $y$ range $[-3,1.2]$, for plain mixing
  (electron colour) and Kerker (`c5`), the active scheme solid and the other dashed; the band $|\zeta|<1$ shaded and
  labelled "converges"; text "ζ < −1: charge sloshing"; a marker at the current $q$, or a downward arrow with the value when
  off scale. Axis labels `q / k_F`, `ζ = 1 − w(q) ε^e(q)`.

**Controls (5).**
1. `q`: $q/k_{\mathrm F}$, 0.05–3, step 0.01.
2. `rs`: $r_s$, 1–6 (dimensionless), step 0.01; tag "(Na)" at 3.93.
3. `w`: mixing weight $w$ (Kerker: $w_0$), 0.05–1.5, step 0.01.
4. `mix`: seg "plain" / "Kerker".
5. `xc`: toggle "exchange kernel (test electron)".

**Readouts.** $\varepsilon^{\mathrm e}(q)$, $1/\varepsilon^{\mathrm e}$, $w(q)$, $\zeta$, cycles to $10^{-6}$, "stable for $w<$ …"
(or "$w_0<$ …").

**Initial state.** $q=0.5\,k_{\mathrm F}$, $r_s=3.93$ (Na), $w=0.15$, plain, exchange off: the iterates alternate in sign
(filled and hollow points) and settle on $1/\varepsilon^{\mathrm e}=0.089$; the diagram values read $1$, $-10.2$, $104$.

**Caption (≈200–260 words).**
- *What to notice.* At the default the loop converges with $\zeta<0$: the iterate overshoots alternately and settles on
  $1/\varepsilon^{\mathrm e}$. The plain series ($w=1$) diverges wherever $\varepsilon^{\mathrm e}>2$, which in a metal means
  every long wavelength; Kerker keeps $w\varepsilon^{\mathrm e}$ within a factor 2 of $w_0$.
- *Marks.* Diagram glyphs and their values; filled/hollow points; the $1/\varepsilon^{\mathrm e}$ line; the convergence band;
  the curve styles.
- *Try.* (1) Raise $w$ above the stable value in the readout: the iterates grow (charge sloshing). (2) Lower $q$ at fixed
  $w$: the failure moves in, since $\varepsilon^{\mathrm e}\simeq(k_{\mathrm{TF}}/q)^2$. (3) Switch to Kerker: the curve in (c)
  flattens near $1-w_0$. (4) Toggle exchange: $1/\varepsilon^{\mathrm e}$ rises slightly, a test electron is screened less.
- *Model.* (A3.12)–(A3.15) for one Fourier component, (A3.18)–(A3.19), the exchange kernel after (A3.6); Kerker with
  $q_0=k_{\mathrm{TF}}$; $r_s=3.93$ is Na-like (literature value, illustrative).

**What in v1 must change.**
- Source and clashing symbols: mixing parameter `β` → `w` ($\beta$ is a Cartesian index); contraction factor `r` → `ζ`
  ($r$ is a distance); cycle `n` → `N_{cyc}` ($n$ is a band index and the density); `∂v_n/∂v^{en}` → `ΔV/ΔV^{en}`;
  `ε(q)` → `ε^e(q)`; "toy xc kernel f_x" → "exchange kernel (toy)"; `(G.53)` and other source tags → eqrefs.
- Diagram labels: v1's `g`, `g^b`, `K`, and the caption sentence on Giustino's Fig. 1(d), are source notation; use
  `ΔV^{KS}`, `ΔV^{en}`, `v+f^{xc}`. The diagrams move from the section body into this figure, with the per-term values
  added.
- Readout keys `β(q)`, `r`, `cycles to 10^{−6}`, `stable for β <` → `w(q)`, `ζ`, `cycles to 10^{−6}`, `stable for w <`.
- v1's caption put Model first; restructure to What to notice, Marks, Try, Model; drop the SI comments (numbers unchanged).
- Classes `a3-scf-…`; the code moves from the section-wide script into the figure file.
