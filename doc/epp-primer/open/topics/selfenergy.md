# Electron self-energies with photon and phonon exchange

Labels used below: **[G]**, Giustino RMP 2017; **[EPP]**, the manuscript (draft v2, with the three corrections that primer v2 builds into Part B and states in footnotes of §B.1, §B.2 and §B.5–§B.8; primed and lettered E tags such as E.428′ and E.448a are the corrected or added equations of [corrections/EPP_corrections.tex](../../corrections/EPP_corrections.tex)); **[lit]**, published work in §6; **[derived here]**, a short derivation made for this note that still needs an independent check; **[proposal]**, a suggestion for the manuscript that appears in neither source; **[speculation]**. Notation and units are those of primer v2 ([§C.0](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#c-conventions)): Gaussian units with $e>0$, so the Coulomb interaction $v=4\pi e^2/q^2$ and the screened $W$ carry $e^2$; §A.$n$, §B.$n$ and §X.$n$ are sections of primer v2, and (A$n$.$m$), (B$n$.$m$), (X$n$.$m$) its equations. Overall volume normalizations of the self-energy formulas still have to be fixed against (G.157). References marked "unverified" could not be confirmed in the last check (§6).

## 1. What is missing in EPP

- **No electron Dyson equation.** EPP solves three bosonic Dyson equations: the total field $D=D^0+D^0\Pi D$ (E.379–381, §6), the electron-dressed field $D^{\mathrm{e}}$ (E.388) and the phonons $F=F^0+F^0\Pi^{\mathrm{ph}}F$ with $\Pi^{\mathrm{ph}}=-J^*D^{\mathrm{e}}J$ (E.398–399, §7; primer (B5.11), (B5.13), (B7.5), (B7.7)). Electrons appear only inside bubbles: $\Pi^{\mathrm{e}}$ (E.383) and $\chi^0$ (p. 92). There is no $\Sigma$ and no dressed $G$ (primer [§X.2](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#x2-missing)).
- **No Fan–Migdal (FM), no Debye–Waller (DW).** Nothing corresponds to (G.151)–(G.161). The same goes for everything built on them: the zero-point renormalization (ZPR) and Varshni shift (G.166–168), lifetimes (G.170), spectral functions (G.171–174) and the cumulant (G §V.B.6, p. 015003-32; primer §A.7, §A.8).
- **The derived vertex is used only once.** The electrodynamic Fröhlich vertices $g^{\mathrm{dyn},0}$ and $g^{\mathrm{dyn}}$ (E.415–416, (B8.6); their $c\to\infty$ limits are E p.95′, E p.98′ and E.448a, primer (B8.15)–(B8.17)) appear only inside the phonon self-energy (E.413, (B8.7)), never on an electron line (§B.8).
- **No electron–two-phonon vertex.** The electron–ion coupling is linear in the displacements: ions enter through Born dipoles, $\mathbf{P}$ linear in $\mathbf{u}$ (E.237, E.241; (B2.10)). The only $u^2$ term, $\rho^{\mathrm{ion}}\rho^{\mathrm{ion}}$ in the Coulomb part of E.237 (with $\rho^{\mathrm{ion}}=-\nabla\cdot\mathbf{P}$), is ion–ion and dresses the phonons, not the electrons. The crystal potential seen by the electrons is frozen, so there is no $\partial^2V^{\mathrm{en}}/\partial u\,\partial u'$ and hence no DW. The $A^2$ terms of E.237 are second order in the field, not in the displacements. There is no short-range channel either.
- **Temperature.** The "$n_{\mathrm{B}}=0$" rule (EPP §4.5, E.246–250) is exact for the retarded RPA boson propagators EPP computes, where $T$ enters only through $f$ in the bubbles (E.248). It fails for $\Sigma$, whose numerators carry $n_{\mathbf{q}\nu}$ (G.158; [§B.5](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b5-rule)).

## 2. Why it matters

**Magnitudes** [lit]:

- **Diamond.** The direct-gap ZPR is 628 meV with $GW$-corrected couplings, more than 40% above DFT (Antonius 2014; G Fig. 13).
- **Survey of many materials.** ZPRs often exceed 0.3 eV, and reach 0.7 eV, in light-element materials. In infrared-active crystals a large share, for many materials the majority, comes from the long-range Fröhlich term, which must be treated nonadiabatically (Miglio 2020).
- **Weak-coupling Fröhlich shift.** For one band edge the weak-coupling estimate is $-\alpha_{\mathrm{F}}\hbar\omega_{\mathrm{LO}}$: about −2.6 meV for the GaAs CBM ($\alpha_{\mathrm{F}}\approx0.071$) and −41 meV for the GaN CBM ($\alpha_{\mathrm{F}}\approx0.44$), both from model parameters (primer §A.8). The gap term is $-(\alpha^{\mathrm{c}}_{\mathrm{F}}+\alpha^{\mathrm{v}}_{\mathrm{F}})\hbar\omega_{\mathrm{LO}}$, with c and v the conduction and valence band edges, and $\alpha^{\mathrm{v}}_{\mathrm{F}}$ is usually the larger, since hole masses are heavier. For zinc-blende GaN, Nery & Allen find a Fröhlich gap correction of −0.05 eV out of a total −0.15 eV (G p. 44: 45 of 150 meV), split as −0.02 eV (CBM) and +0.03 eV (VBM). These are analytic small-$q$ corrections added to an adiabatic grid calculation, so they are not the same quantity as the model value $-\alpha_{\mathrm{F}}\hbar\omega_{\mathrm{LO}}$.
- **Strong coupling.** MgO ($\alpha_{\mathrm{F}}\approx1.7$ from the model parameters of primer §A.8) lies outside weak coupling. Nery et al. (2018) use MgO and LiF to compare cumulant and Dyson satellites.
- **Lattice screening in $W$.** It can lower $GW$ gaps by more than 15% in highly polar compounds (Botti & Marques 2013). Lambrecht et al. (2017) object that the polar term was not confined in $q$, a sign of how strongly the result depends on dynamics.
- **Light cone.** For GaAs, $q_{\mathrm{light}}(\omega_{\mathrm{LO}})=\sqrt{\varepsilon_\infty}\,\omega_{\mathrm{LO}}/c\approx6\times10^{-5}$ Å$^{-1}$, against the polaron wave vector $q_{\mathrm{pol}}=(2m^*\omega_{\mathrm{LO}}/\hbar)^{1/2}\approx0.025$ Å$^{-1}$. The cone holds a fraction $\approx(2/\pi)\,q_{\mathrm{light}}(\omega_{\mathrm{LO}})/q_{\mathrm{pol}}\approx1.5\times10^{-3}$ of the Fröhlich-shift integrand, so finite $c$ can affect at most about $10^{-3}$ of the shift. Even inside the cone it changes only the transverse (TO-polariton) channel, because the LO kernel is $c$-independent (primer §A.8, §B.8, §X.1).

**For an RPA/$GW$/BSE practitioner:**

1. **A $GW$ gap is a clamped-nuclei gap.** Without the ZPR, agreement with experiment at the 0.1–0.5 eV level can be fortuitous, especially in light-element and polar crystals such as hBN.
2. **Temperature-dependent BSE needs input.** It requires $\Sigma^{\mathrm{ep}}$ energies and widths (BSE+Allen–Heine, G p. 42). Phonon screening of excitons, its two-particle counterpart, lowers binding energies (Filip 2021; unverified, §6).
3. **In 2D the polar channel behaves differently.** The 2D Fröhlich vertex tends to a constant as $Q_\parallel\to0$ (Sohier 2016; unverified, §6; cf. E.469–472); primer v2 now writes this sheet vertex ([§B.10](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b10-coupling)). The adiabatic integrand $\propto|g|^2/Q_\parallel^2\times Q_\parallel\,\mathrm{d}Q_\parallel$ then diverges logarithmically, while the nonadiabatic integral converges [derived here]. The sheet's static screening, $\varepsilon^{\mathrm{2D}}(Q_\parallel)\approx\varepsilon_{\mathrm{env}}+r_0Q_\parallel$ with $\varepsilon_{\mathrm{env}}=(\varepsilon_{\mathrm{sub}}+\varepsilon_{\mathrm{cap}})/2$ ((B10.8)), sets both the Fröhlich vertex and the Rytova–Keldysh exciton potential. Encapsulation changes $\varepsilon_{\mathrm{env}}$, which enters both at once; $r_0=2\pi\chi^{\mathrm{2D}}$ ((B10.7)) is a property of the sheet.
4. **The add-on is cheap.** Your code already has $W^{\mathrm{e}}$. The lattice adds one pole per polar mode to its head, and that pole follows from EPP's equations (E.388, E.393, E.398; derived in §4.4 of this note).

## 3. What standard theory does

- **Exact regrouping** [G]. $G=G^{\mathrm{cn}}+G^{\mathrm{cn}}\Sigma^{\mathrm{ep}}G$ with $\Sigma^{\mathrm{ep}}=\Sigma^{\mathrm{FM}}+\Sigma^{\mathrm{DW}}+\Sigma^{\mathrm{dGW}}$ (G.151–155; (A7.3)–(A7.4)). FM is $\mathrm{i}\hbar G\Lambda W^{\mathrm{ph}}$ with $W^{\mathrm{ph}}=W^{\mathrm{e}}\chi^{\mathrm{ion}}W^{\mathrm{e}}$ (G.105, G.110; (A5.23)), screened at both ends (G.156; (A7.5)). DW is the static potential of the smeared nuclei (G.154, G.159). dGW has never been evaluated (G p. 27).
- **Working formulas** [G]. Approximations (i)–(vii) lead to (G.157)–(G.161) and to Allen–Heine (G.164–167). The rigid-ion DW (G.195–198) follows from the translational sum rules (G.193–194), and FM and DW partly cancel.
- **Polar materials** [lit]. The vertex is split as $g=g^{\mathrm{sr}}+g^{\mathrm{lr}}$ (Verdi & Giustino 2015) and Wannier-interpolated. The adiabatic formula diverges (G p. 43). The cures are nonadiabatic denominators, analytic small-$q$ integration (Nery & Allen 2016) and spectral functions.
- **Combining with $GW$.** Route (a) takes $G^{\mathrm{cn}}$ from $GW$ and adds FM+DW, optionally with $GW$-corrected $g$ (Antonius 2014). Route (b) puts lattice screening into $W$ (Botti & Marques; Lambrecht), which is G.104 used literally.
- **Beyond Dyson.** Options are the cumulant (Nery 2018; Kandolf 2022 for the doped Fröhlich model) and polaron self-trapping (Lafuente-Bartolome 2022; unverified, §6). For the 2D long-range vertex and 2D polarons, see the 2D topic note. Mobilities are reviewed by Poncé (2020), and free-carrier screening of the vertex is treated by Macheda (2022) (both unverified, §6; see also [plasmon.md](plasmon.md)).

## 4. How it could enter the EPP formalism

**4.1 Electron Dyson equation [proposal].** The interaction $V'$ (E.237; (B5.9)) gives two electron self-energies: a paramagnetic one from two $-\frac1c\mathbf{j}\cdot\mathbf{A}$ vertices and a diamagnetic tadpole from $\frac{e^2}{2m_{\mathrm{e}}c^2}n^{\mathrm{e}}A^2$:

$$G=G^{\mathrm{KS}}+G^{\mathrm{KS}}\big(\Sigma^{\mathrm{EPP}}-V^{\mathrm{xc}}\big)G,\qquad \Sigma^{\mathrm{EPP}}=\Sigma^{\mathrm{para}}+\Sigma^{\mathrm{dia}},$$

$$\Sigma^{\mathrm{para}}_{nn'\mathbf{k}}(\omega)=\mathrm{i}\hbar\sum_m\!\int\!\frac{\mathrm{d}^3q}{(2\pi)^3}\!\int\!\frac{\mathrm{d}\omega'}{2\pi}\mathrm{e}^{\mathrm{i}\omega'\eta}\,G_{m\mathbf{k}+\mathbf{q}}(\omega-\omega')\Big[-\frac1c\,j_\alpha^{\,n\mathbf{k},m\mathbf{k}+\mathbf{q}}D_{\alpha\beta}(\mathbf{q},\omega')\big(j_\beta^{\,n'\mathbf{k},m\mathbf{k}+\mathbf{q}}\big)^*\Big],$$

$$\Sigma^{\mathrm{dia}}(\mathbf{r})=\frac{e^2}{2m_{\mathrm{e}}c^2}\langle A^2(\mathbf{r})\rangle,\qquad \langle A_\alpha A_\alpha\rangle=-\mathrm{i}\hbar c\!\int\!\frac{\mathrm{d}\omega}{2\pi}\mathrm{e}^{\mathrm{i}\omega\eta}D_{\alpha\alpha}(\mathbf{r},\mathbf{r};\omega).$$

$\Sigma^{\mathrm{EPP}}$ contains the full electronic exchange (through $D^{0\mathrm{L}}$) and the RPA correlation (through the longitudinal $D^{\mathrm{e}}$, i.e. $\Sigma^{\mathrm{L}}[W^{\mathrm{e}}]$ of §4.3), so the $V^{\mathrm{xc}}$ already in $G^{\mathrm{KS}}$ must be subtracted. The kernel is $K=-\frac1c\,\mathbf{j}^*\cdot D\cdot\mathbf{j}$. With $D^{0\mathrm{L}}$ (E.374; (B5.8)) and continuity (E.428′; (X2.8)) it becomes $(\Delta E/\hbar\omega')^2\,v\,|\rho|^2$, which equals the Coulomb kernel $v\,|\rho|^2$ only at $\omega'=\Delta E/\hbar$. Here $\rho=\rho_{n\mathbf{k},m\mathbf{k}+\mathbf{q}}=\langle\psi_{n\mathbf{k}}|\mathrm{e}^{-\mathrm{i}\mathbf{q}\cdot\mathbf{r}}|\psi_{m\mathbf{k}+\mathbf{q}}\rangle$ is the number-density matrix element of E.106 (the charge vertex of (B8.10) is $-e\rho$). Off that shell the temporal-gauge para term differs from Coulomb exchange, and only para + dia reproduces it on shell (§4.3).

**4.2 Vector form of $W=W^{\mathrm{e}}+W^{\mathrm{ph}}$ [derived here from E.379, E.387–388, E.393, E.398].** Resumming the ions gives exactly, within EPP's RPA,

$$D=D^{\mathrm{e}}+D^{\mathrm{ph}},\qquad D^{\mathrm{ph}}=-D^{\mathrm{e}}JFJ^{\dagger}D^{\mathrm{e}},$$

with $J$ the dynamic ionic vertex (B7.3) and $F$ the dressed phonon propagator (B7.7). Primer v2 names only the Coulomb sector of $D^{\mathrm{ph}}$, $W^{\mathrm{ph}}=e^2D^{\mathrm{ph}}_{\rho\rho}$. The phonon part of $K$ is therefore $\frac1c(\mathbf{j}^*D^{\mathrm{e}}J)\,F\,(J^\dagger D^{\mathrm{e}}\mathbf{j})$, with screened vertices at both ends, as in (G.156). The phonon self-energy (E.412) instead pairs one bare vertex with one screened vertex. $\Sigma^{\mathrm{para}}[D^{\mathrm{e}}]$ is "$GW$ with the total, electron-dressed field" and $\Sigma^{\mathrm{para}}[D^{\mathrm{ph}}]$ is the electrodynamic Fan–Migdal term.

**4.3 Gauge [derived here].** In the temporal gauge the longitudinal block of the dressed field is $D^{\mathrm{e,L}}_{\alpha\beta}=-\frac{c}{\omega^2}\partial_\alpha\partial'_\beta D^{\mathrm{e}}_{\rho\rho}$ with $e^2D^{\mathrm{e}}_{\rho\rho}=W^{\mathrm{e}}$ (E.436–437, (B6.9); for $c\to\infty$, E.441 and E.445, (B9.6)–(B9.7)), i.e. $D^{\mathrm{L}}_{\alpha\beta}=-\frac{c}{e^2\omega^2}q_\alpha q_\beta\,W$ in Fourier space. Combined with continuity (E.428′) this gives $K^{\mathrm{L}}=(\Delta E/\hbar\omega')^2|\rho|^2W$, with $\Delta E=E_{n\mathbf{k}}-E_{m\mathbf{k}+\mathbf{q}}$. On shell ($\hbar\omega=E_{n\mathbf{k}}$) the excess over the Coulomb kernel, $[(\Delta E/\hbar\omega')^2-1]|\rho|^2W$, cancels the pole of $G$ and leaves $\big(\Delta E/(\hbar\omega')^2+1/\hbar\omega'\big)\,|\rho|^2W$.

- **First term.** The Bethe sum rule for one Bloch state (local $H$) is $\sum_m(E_{m\mathbf{k}+\mathbf{q}}-E_{n\mathbf{k}})|\rho_{n\mathbf{k},m\mathbf{k}+\mathbf{q}}|^2=\hbar^2q^2/2m_{\mathrm{e}}+\hbar\mathbf{q}\cdot\mathbf{v}_{n\mathbf{k}}$. The linear term vanishes at band extrema, and in general after $\mathbf{q}$ integration, because $W(\mathbf{q},\omega')=W(-\mathbf{q},\omega')$ (time reversal). What remains, $-\frac{\mathrm{i}\hbar}{2m_{\mathrm{e}}}\int q^2W/\omega'^2$, cancels $\Sigma^{\mathrm{dia,L}}$ exactly.
- **Second term.** $\sum_m|\rho|^2W/\hbar\omega'$ is independent of $n$, since $\sum_m|\rho_{n\mathbf{k},m\mathbf{k}+\mathbf{q}}|^2=1$. With a symmetric principal-value prescription it vanishes; otherwise it is a rigid shift that cannot change gaps.

On shell, therefore, para + dia is gauge invariant up to that $n$-independent constant. This is the self-energy analogue of E.310 ((B5.25)). It requires complete band sums and a local Hamiltonian, and off shell the two gauges differ.

**Proposal:** evaluate $\Sigma$ in the Coulomb gauge, building the propagators from EPP's gauge-invariant $\Pi$ (E.281, E.311; (B5.19), (B5.26)):

$$\Sigma^{\mathrm{EPP}}\simeq\Sigma^{\mathrm{L}}[W]+\Sigma^{\mathrm{T}}\big[D^{\mathrm{T}}-D^{0\mathrm{T}}\big]+\Sigma^{\mathrm{dia,T}}\big[D^{\mathrm{T}}-D^{0\mathrm{T}}\big],\qquad W=\frac{4\pi e^2}{\mathbf{q}\cdot\boldsymbol{\varepsilon}(\mathbf{q},\omega)\cdot\mathbf{q}},\quad \varepsilon_{\alpha\beta}=\delta_{\alpha\beta}+\frac{4\pi c}{\omega^2}\Pi_{\alpha\beta}.$$

$\Sigma^{\mathrm{L}}$ is a $GW$-type self-energy with the density vertices $\rho$ of §4.1; with a KS reference, $V^{\mathrm{xc}}$ is subtracted as in §4.1. For isotropic $\Pi$ the split is exact for the full $D$ at any $c$, with $W=4\pi e^2/(q^2\varepsilon^{\mathrm{L}})$ ((B6.17)); in general it is exact as $c\to\infty$ (E.441, E.445).

**4.4 The phonon pole of $W$ [derived here; checked against the Lyddane–Sachs–Teller (LST) form].**

$$W^{\mathrm{ph}}(\mathbf{q},\omega)=W-W^{\mathrm{e}}=\frac{\big[W^{\mathrm{e}}(\mathbf{q},\omega)\big]^2}{e^2}\sum_{\nu\nu'}\frac{(\mathbf{q}\cdot\mathbf{j}_\nu)(\mathbf{q}\cdot\mathbf{j}_{\nu'})^*}{2\Omega_{\mathrm{V}}\sqrt{\omega_{\mathbf{q}\nu}\omega_{\mathbf{q}\nu'}}}\,\hbar F_{\nu\nu'}(\mathbf{q},\omega).$$

Here $\mathbf{j}_\nu$ is the ionic current vertex (B5.22) at $\mathbf{G}=0$ and $\omega_{\mathbf{q}\nu}$ the bare frequency. For one isotropic mode with static screening, this is $\frac{4\pi e^2}{q^2}\frac{(\omega^{\mathrm{ion}}_{\mathrm{p}})^2/\varepsilon_\infty^2}{\omega^2-\omega_{\mathrm{LO}}^2}$: Fröhlich's $|g^{\mathrm{lr}}|^2$ times the LO propagator (§A.5, (A5.30)). Inserted into $\Sigma^{\mathrm{L}}$ at $\mathbf{G}=\mathbf{G}'=0$, it gives (G.157) with the vertex of E.448a ((B8.16)–(B8.17)). Three things are new compared with (G.157):

- $\omega_{\mathrm{LO}}$ is an output (E.454; (B9.10)), not an input;
- the screening is $\varepsilon^{\mathrm{e}}(\mathbf{q},\omega)$, so approximation (v), the static vertex (5 in primer [§A.7](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#a7-working)), is dropped and doped crystals get the $\mathrm{L}_\pm$ modes of E.453. Primer v2 now works these out ([§B.9](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b9-dynamic), (B9.16)–(B9.19)); the EPP manuscript stops at the static E.454 (see [plasmon.md](plasmon.md));
- the residue $\omega_{\mathbf{q}\nu}/\bar\omega_{\mathbf{q}\nu}$ of $\hbar F$ fixes the normalization (§X.2, (X2.12)).

**4.5 Temperature [proposal; this is the standard spectral form].** Do the Matsubara sum once, with $\mathcal{B}=-\frac\hbar\pi\mathrm{Im}\,K^{\mathrm{R}}$:

$$\Sigma^{\mathrm{R}}_{n\mathbf{k}}(\omega)=\sum_m\!\int\!\mathrm{d}^3q\!\int_0^\infty\!\mathrm{d}\omega'\,\mathcal{B}_{nm}(\mathbf{k},\mathbf{q},\omega')\Big[\frac{1-f_m+n_{\mathrm{B}}(\omega')}{\hbar(\omega+\mathrm{i}\eta)-E_m-\hbar\omega'}+\frac{f_m+n_{\mathrm{B}}(\omega')}{\hbar(\omega+\mathrm{i}\eta)-E_m+\hbar\omega'}\Big]+\Sigma^{\mathrm{x}}.$$

Here $\Sigma^{\mathrm{x}}$ is the static exchange part. For $\mathcal{B}=|g|^2\delta(\omega'-\omega_{\mathbf{q}\nu})$ this is (G.158). $D$ itself keeps EPP's rule, since its only temperature dependence is $f(T)$ in $\Pi^{\mathrm{e}}$ (E.248). The Bose factors enter only here.

**4.6 Debye–Waller analogue [assessment].** In the Coulomb gauge $\Sigma^{\mathrm{dia}}$ involves only $\langle(\mathbf{A}^{\mathrm{T}})^2\rangle$, whose ionic part lives at $q\lesssim q_{\mathrm{light}}(\omega)$, on micrometre scales. With EPP's macroscopic ($\mathbf{G}=\mathbf{G}'=0$) $D$, any $A^2$ tadpole is a rigid shift of all bands and cannot change a gap. A real DW needs the second-order coupling of (G.159), i.e. a short-range channel in the Lagrangian. Its FM and DW parts must then be added together so that (G.193)–(G.194) hold.

**4.7 Retardation [assessment].**

- **Magnitude.** $\Sigma^{\mathrm{T}}$ is of order $(v_{n\mathbf{k}}/c)^2\approx10^{-7}$–$10^{-5}$ of the Coulomb term for electronic (Breit-type) exchange, with band velocities $v_{n\mathbf{k}}=10^5$–$10^6$ m/s, and of order $(q_{\mathrm{light}}(\omega)/q)^2\approx10^{-6}$–$10^{-4}$ for phonon-frequency exchange at the polaron scale (§X.1).
- **Photon emission within a band.** On the photon-like branch, real emission is forbidden by the Cherenkov condition, since $v_{n\mathbf{k}}\ll c/\sqrt{\varepsilon}$ [derived here]. With the ionic $\Pi$, $\varepsilon(\omega)$ diverges just below $\omega_{\mathrm{TO}}$, so on the phonon-like lower polariton branch near $\omega_{\mathrm{TO}}$ emission is kinematically allowed. That is TO-phonon emission, with a transverse coupling $\propto(\omega/cq)^2$, and it is already contained in the transverse part of $\Sigma^{\mathrm{para}}[D^{\mathrm{ph}}]$.
- **Interband emission belongs to the optics (BSE) topic.** It is radiative recombination, which needs a hole, so it is a two-particle (exciton) process.
- **Vacuum subtraction.** The vacuum part $D^{0\mathrm{T}}$ is UV divergent (mass renormalization) and already contained in $m_{\mathrm{e}}$. It is subtracted in both transverse terms, $\Sigma^{\mathrm{T}}[D^{\mathrm{T}}-D^{0\mathrm{T}}]$ and $\Sigma^{\mathrm{dia,T}}[D^{\mathrm{T}}-D^{0\mathrm{T}}]$: $\langle(\mathbf{A}^{\mathrm{T}})^2\rangle_{\mathrm{vac}}$ is also UV divergent, though band independent. This is a vacuum subtraction. It is not the $\Pi(\omega)-\Pi(0)$ of E.281/E.310, which is the gauge-invariance cancellation of the diamagnetic term.

## 5. Pitfalls

**Double counting.**

- **Force constants.** Bare $\omega_{\mathbf{q}\nu}$ must come from short-range force constants $\Phi^{\mathrm{sr}}$ ([§B.1](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b1-range)). A DFPT $\omega_{\mathrm{LO}}$ in $F$, dressed again by $\Pi^{\mathrm{ph}}$, counts the LO–TO splitting twice.
- **KS reference.** Subtract $V^{\mathrm{xc}}$, since $\Sigma^{\mathrm{L}}[W^{\mathrm{e}}]$ is the electronic $GW$ self-energy. Only with a $GW$ reference ($G^{\mathrm{cn}}$ from $GW$) does one keep just $D^{\mathrm{ph}}$ and $D^{\mathrm{e,T}}-D^{0\mathrm{T}}$.
- **The polar vertex.** If DFPT $g$ is used for FM, subtract $g^{\mathrm{lr}}$ (same $\boldsymbol{\varepsilon}_\infty$, $Z^*$) before adding EPP's polar part. Subtracting $g^{\mathrm{lr}}$ also removes the sr–lr cross term $2\,\mathrm{Re}(g^{\mathrm{sr}}g^{\mathrm{lr}*})$, which $W^{\mathrm{ph}}$ does not contain, since EPP has no short-range channel. Either keep it explicitly as FM$[g^{\mathrm{sr}},g^{\mathrm{lr}}]$ cross terms evaluated with the same $\boldsymbol{\varepsilon}_\infty$ and $Z^*$, or use EPP's $W^{\mathrm{ph}}$ only as a replacement for the $|g^{\mathrm{lr}}|^2$ part and state that the cross term is kept from DFPT.
- **Static $\varepsilon_{\mathrm{s}}$.** Never add FM to a $W$ screened with $\varepsilon_{\mathrm{s}}$.
- **Hartree terms.** These are already in $V^{\mathrm{KS}}$.

**Screening.** EPP's $\varepsilon^{\mathrm{e}}$ is the head of $\varepsilon^{\mathrm{e}}_{\mathrm{RPA}}$, without local fields or xc ([§B.9](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b9-survivor)). $GW$ codes include local fields and DFPT uses $\varepsilon^{\mathrm{e}}_{\mathrm{RPA+xc}}$. A mismatch distorts $W^{\mathrm{ph}}$ at intermediate $q$.

**Vertices.** FM takes two screened vertices and the phonon self-energy one bare and one screened. Mixing them up repeats the miscount of [§B.8](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b8-counting).

**Gauge.** The temporal gauge needs the $A^2$ term, complete band sums and a symmetric $\omega'\to0$ prescription (its $1/\omega'^2$ pole). Truncated bands leave spurious terms, the analogue of the fake Drude weight of [§B.5](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b5-truncation). Nonlocal pseudopotentials and scissors break E.428′ ((X2.8)) and the Bethe sum rule.

**Sum rules.**

- the acoustic sum rule for Born charges, $\sum_\kappa\mathbf{Z}^*_\kappa=0$ ((A4.4));
- the basis phase of the ionic vertex (E p.86a′; (B5.22)), which removes spurious acoustic coupling;
- the translational sum rules (G.193–194) once DW is included.

**Convergence.**

- **3D.** The $1/q^2$ singularity needs nonadiabatic denominators, analytic small-$q$ integration and interpolated grids. The error decays only as $1/N_1$, with $N_1$ points per direction ($\propto N_{\mathbf{q}}^{-1/3}$ in the total number of points $N_{\mathbf{q}}$; primer §A.8).
- **2D.** The adiabatic integral diverges logarithmically; use nonadiabatic denominators. Independently, compute $Z^*$, $\varepsilon$ and $g$ with a 2D-truncated Coulomb interaction, so that the small-$Q_\parallel$ vertex is the true 2D one (Sohier et al.; unverified, §6). Integrate the region $Q_\parallel\lesssim\varepsilon_{\mathrm{env}}/r_0$ analytically or with mini-BZ averaging, since $W^{\mathrm{ph}}$ is concentrated there.
- **Broadening.** Replace $\hbar\eta$ by the physical linewidth $\Gamma^{\tau}_{n\mathbf{k}}$ (G.166).
- **Light line.** The pole of $D^{\mathrm{T}}$ at $q<q_{\mathrm{light}}(\omega)$ needs analytic integration.
- **Strong coupling.** For satellites use a second-order cumulant (Nery 2018); Dyson–Migdal gives better quasiparticle masses (Kandolf 2022). For $\alpha_{\mathrm{F}}\gtrsim2$–3, or self-trapping, use polaron methods (Lafuente-Bartolome 2022; unverified, §6).

## 6. State of the art (verified)

1. F. Giustino, "Electron-phonon interactions from first principles", Rev. Mod. Phys. **89**, 015003 (2017), doi:10.1103/RevModPhys.89.015003. The reference framework.
2. P. B. Allen, V. Heine, "Theory of the temperature dependence of electronic band structures", J. Phys. C **9**, 2305 (1976) (confirmed via Giustino's bibliography; DOI not confirmed). Founds Allen–Heine theory, with DW from translational invariance.
3. C. Verdi, F. Giustino, "Fröhlich electron-phonon vertex from first principles", Phys. Rev. Lett. **115**, 176401 (2015), doi:10.1103/PhysRevLett.115.176401. The polar vertex $g^{\mathrm{lr}}$, EPP's $c\to\infty$ limit.
4. S. Botti, M. A. L. Marques, "Strong renormalization of the electronic band gap due to lattice polarization in the GW formalism", Phys. Rev. Lett. **110**, 226404 (2013), doi:10.1103/PhysRevLett.110.226404. Lattice polarization in $GW$'s $W$; more than 15% of the gap in highly polar compounds.
5. W. R. L. Lambrecht, C. Bhandari, M. van Schilfgaarde, Phys. Rev. Materials **1**, 043802 (2017), doi:10.1103/PhysRevMaterials.1.043802. Revisits (4) and confines the polar term in $q$.
6. G. Antonius, S. Poncé, P. Boulanger, M. Côté, X. Gonze, "Many-body effects on the zero-point renormalization of the band structure", Phys. Rev. Lett. **112**, 215501 (2014), doi:10.1103/PhysRevLett.112.215501. $GW$-enhanced coupling, diamond ZPR (628 meV, G Fig. 13).
7. J. P. Nery, P. B. Allen, "Influence of Fröhlich polaron coupling on renormalized electron bands in polar semiconductors: Results for zinc-blende GaN", Phys. Rev. B **94**, 115135 (2016), doi:10.1103/PhysRevB.94.115135. Analytic small-$q$ Fröhlich correction: −0.05 eV of the −0.15 eV gap ZPR.
8. A. Miglio et al., npj Comput. Mater. **6**, 167 (2020), doi:10.1038/s41524-020-00434-z. Nonadiabatic ZPR survey; ZPRs above 0.3 eV, up to 0.7 eV.
9. J. P. Nery, P. B. Allen, G. Antonius, L. Reining, A. Miglio, X. Gonze, "Quasiparticles and phonon satellites in spectral functions of semiconductors and insulators: Cumulants applied to the full first-principles theory and the Fröhlich polaron", Phys. Rev. B **97**, 115145 (2018), doi:10.1103/PhysRevB.97.115145. Dyson versus cumulant satellites (MgO, LiF).
10. N. Kandolf, C. Verdi, F. Giustino, "Many-body Green's function approaches to the doped Fröhlich solid: Exact solutions and anomalous mass enhancement", Phys. Rev. B **105**, 085148 (2022), doi:10.1103/PhysRevB.105.085148. Analytic Dyson–Migdal and cumulant solutions of the doped (3D) Fröhlich model; anomalous mass enhancement at finite doping.

Unverified (could not be re-confirmed in the last check; confirm before publishing):

- J. Lafuente-Bartolome et al., Phys. Rev. Lett. **129**, 076402 (2022), doi:10.1103/PhysRevLett.129.076402. Allen–Heine theory and polaron self-trapping in one framework.
- S. Poncé, W. Li, S. Reichardt, F. Giustino, Rep. Prog. Phys. **83**, 036501 (2020), doi:10.1088/1361-6633/ab6a43. Mobilities in bulk and 2D materials.
- T. Sohier, M. Calandra, F. Mauri, Phys. Rev. B **94**, 085415 (2016), doi:10.1103/PhysRevB.94.085415. The 2D Fröhlich interaction in TMD monolayers (listed as verified in the 2D topic note; not re-checked here).
- F. Macheda, P. Barone, F. Mauri, Phys. Rev. Lett. **129**, 185902 (2022), doi:10.1103/PhysRevLett.129.185902. Doping screens the Fröhlich vertex and the LO–TO splitting.
- M. R. Filip, J. B. Haber, J. B. Neaton, Phys. Rev. Lett. **127**, 067401 (2021), doi:10.1103/PhysRevLett.127.067401. Phonon screening in the BSE, 3D halide perovskites (listed as verified in the optics and 2D topic notes; not re-checked here).

## 7. First steps

1. **Derive (1–2 weeks).** Write out Sec. 4.1–4.4 with all prefactors, including the $V^{\mathrm{xc}}$ subtraction. Prove the $D^{\mathrm{ph}}$ identity and the on-shell gauge equivalence. Check that $c\to\infty$ gives (G.157) with E.448a and the residue $\omega_{\mathbf{q}\nu}/\bar\omega_{\mathbf{q}\nu}$.
2. **Analytic benchmark (days).** For an isotropic Fröhlich crystal (GaAs, GaN, MgO), recover the band-edge shift $-\alpha_{\mathrm{F}}\hbar\omega_{\mathrm{LO}}$, the §A.8 rate formula and the Bose-factor $T$ dependence. Then confirm that the finite-$c$ corrections are of order $10^{-3}$ or less and come only from the transverse channel.
3. **2D test (2–4 weeks).** Add a head-only $W^{\mathrm{ph}}(Q_\parallel,\omega)$ to your $GW$ code, built from E.469–472 with the expansion $\varepsilon^{\mathrm{2D}}\approx\varepsilon_{\mathrm{env}}+r_0Q_\parallel$ of §B.10 ((B10.8)). Use monolayer hBN as the main test.
   - Use the code's truncated-Coulomb head and mini-BZ averaging for $Q_\parallel\lesssim\varepsilon_{\mathrm{env}}/r_0$.
   - Show the adiabatic log divergence and the nonadiabatic convergence.
   - Benchmark the vertex $|g(Q_\parallel)|$ implied by $W^{\mathrm{ph}}$ against Sohier 2016 and the sheet vertex of [§B.10](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b10-coupling), and the 2D-Fröhlich ZPR against the analytic nonadiabatic integral.
   - Vary the substrate through $\varepsilon_{\mathrm{env}}$ at fixed $r_0$, with $W^{\mathrm{e,2D}}=2\pi e^2/[Q_\parallel(\varepsilon_{\mathrm{env}}+r_0Q_\parallel)]$ ((B10.9)). This is an addition beyond EPP's E.461, which is for a freestanding sheet.
   - Treat MoS$_2$ only as a weak-polar control (the Fröhlich part of the e–ph shift only). Its e–ph renormalization is dominated by short-range couplings, so a comparison with full ZPRs would mislead.
4. **Doped crystals (weeks).** Use the dynamical $\varepsilon^{\mathrm{e}}(\mathbf{q},\omega)$ with $\mathrm{L}_\pm$ modes in n-GaAs ([§B.9](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b9-dynamic)) or gated MoS$_2$ ([§B.10](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b10-doped)), and compute lifetimes and satellites. Compare with Kandolf 2022 (a 3D model benchmark) and Macheda 2022. The coupled modes themselves are the topic of [plasmon.md](plasmon.md).
5. **Toward BSE (months).** Add the short-range channel with DW, a cumulant and coupling to BSE.

## 8. Open questions

- **Which gauge?** Off shell the gauges differ. Which one gives the measured spectral function, given that ARPES couples through $\mathbf{A}$?
- **A consistent DW.** What is the minimal extension of the Lagrangian (short range, second order in $\mathbf{u}$) that keeps gauge invariance and the translational sum rules?
- **Gated 2D crystals.** There the acoustic $\sqrt{Q_\parallel}$ plasmon overlaps the phonons ([§B.10](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b10-doped); [plasmon.md](plasmon.md)), so is the split $\Sigma^{\mathrm{e}}+\Sigma^{\mathrm{ph}}$ still meaningful?
- **Head-only $W^{\mathrm{ph}}$.** How large is its error next to a $W^{\mathrm{e}}$ with local fields?
- **dGW.** How large is it?
- **[speculation] A sheet's own light cone.** Does it ($Q_\parallel\lesssim\omega/c$ for a freestanding sheet) give any single-particle effect, or only excitonic radiative widths?

## 9. Review responses

All 23 points were checked against EPP (E.237–250, E.281, E.374, E.428, E.436–445, E.461–472), Giustino (pp. 30–32, 43–44 and the bibliography) and primer v1 §A.8, §B.10 and §X.1. Every point was accepted. This pass could not re-search citations: the shared web-search budget was used up, and doi.org, Crossref and the publishers were unreachable.

- **Fröhlich share of the ZPR (major).** Accepted. Giustino p. 44 confirms 45 of 150 meV for zb-GaN. §2 now separates the band-edge model value ($-\alpha_{\mathrm{F}}\hbar\omega_{\mathrm{LO}}$, GaN CBM ≈ −40 meV; recomputed $\alpha_{\mathrm{F}}$ = 0.44 from the §A.8 parameters) from the gap term $-(\alpha^{\mathrm{c}}_{\mathrm{F}}+\alpha^{\mathrm{v}}_{\mathrm{F}})\hbar\omega_{\mathrm{LO}}$ and from Nery & Allen's analytic correction. The CBM/VBM split (−0.02/+0.03 eV) is taken from the reviewer's reading of their abstract. It is consistent with the confirmed −0.05 eV total. Miglio 2020 is now quoted as "for many materials".
- **$K^{0\mathrm{L}}$ off shell (major).** Accepted and rederived: $-\frac1c\mathbf{j}^*D^{0\mathrm{L}}\mathbf{j}$ with E.374 and E.428 gives $(\Delta E/\hbar\omega')^2v|\rho|^2$. §4.1 now agrees with §4.3.
- **$V^{\mathrm{xc}}$ double counting (major).** Accepted: the Dyson equation is now $G^{\mathrm{KS}}(\Sigma^{\mathrm{EPP}}-V^{\mathrm{xc}})G$, with a new "KS reference" bullet in §5 that is merged with the $GW$ reference bullet.
- **sr–lr cross term (major).** Accepted; the two options are stated in §5.
- **Substrate through $\varepsilon_{\mathrm{env}}$, not $r_0$ (major).** Accepted and checked against §B.10 and E.469–470: E.461 has no environment, and E.470 is $1-v^{\mathrm{2D}}\chi^0$ ((B10.12)). §2 item 3 and §7 step 3 have been rewritten.
- **2D convergence vs truncation (major).** Accepted: these are now two separate requirements.
- **Bethe sum rule $+\mathbf{q}\cdot\mathbf{v}_{n\mathbf{k}}$, $n$-independent remainder, E.436–437 citation.** Accepted and rederived ($\mathrm{e}^{-\mathrm{i}\mathbf{q}\cdot\mathbf{r}}H\mathrm{e}^{\mathrm{i}\mathbf{q}\cdot\mathbf{r}}-H=-\mathrm{i}\hbar^2\mathbf{q}\cdot\nabla/m_{\mathrm{e}}+\hbar^2q^2/2m_{\mathrm{e}}$; completeness gives $\sum_m|\rho|^2=1$). §4.3 has been restructured accordingly.
- **Attribution of the $W^{\mathrm{ph}}$ pole, cumulant citation (G §V.B.6, p. 32), the linear-in-$\mathbf{u}$ statement (the $\rho^{\mathrm{ion}}\rho^{\mathrm{ion}}$ term of E.237 checked), the $n_{\mathrm{B}}$ rule, $N_1$ per direction (checked in §A.8: $q_{\mathrm{min}}=q_{\mathrm{D}}/N_1$).** All accepted.
- **$(v_{n\mathbf{k}}/c)^2$ vs $(q_{\mathrm{light}}(\omega)/q)^2$; light-cone fraction as an upper bound ($(2/\pi)q_{\mathrm{light}}(\omega_{\mathrm{LO}})/q_{\mathrm{pol}}=1.5\times10^{-3}$ recomputed); Cherenkov on the lower polariton branch; vacuum subtraction in $\Sigma^{\mathrm{dia,T}}$ (E.281 checked: it is the dia/para(0) cancellation).** All accepted.
- **§7 step 3 and step 4, strong-coupling advice, Kandolf description, MgO ($\alpha_{\mathrm{F}}$ = 1.725 recomputed).** All accepted. hBN is kept as the main test.
- **Citations.** Accepted. The two errata, Leibbrandt 1987, Palummo 2015 (with its "few ps at 4 K" figure) and Sio & Giustino 2023 were dropped, since none could be confirmed and none carried a load-bearing claim. Five core references are kept but marked unverified in §6 and in the text. The verified list now has 10 papers. The reviewer's "11" names only these 10. Titles for refs 2, 3, 4, 6 and 7 were checked here against Giustino's bibliography; Allen & Heine is counted as verified on that basis, with its DOI still unconfirmed. Refs 5 and 8–10 rest on the earlier search-based confirmation, which the reviewer repeated. Nery 2018 now has its full title, and Kandolf 2022 has its title and a corrected description.
