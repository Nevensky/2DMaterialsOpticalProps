# Ordered phases: superconductivity, charge density waves and condensates

Labels: **[G]** Giustino RMP 2017; **[EPP]** the manuscript (draft v2, with the three corrections that primer v2 builds into Part B and states in footnotes of §B.1, §B.2 and §B.5–§B.8; primed and lettered E tags such as E.428′ and E.454a are the corrected or added equations of [corrections/EPP_corrections.tex](../../corrections/EPP_corrections.tex)); **[lit]** published work in §6; **[derived here]** derived for this note, unchecked; **[proposal]**; **[speculation]**. Notation and units follow the [conventions](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#c-conventions) of primer v2: Gaussian units with $e>0$, so $v=4\pi e^2/q^2$ and $W$ carry $e^2$; (A$n$.$m$) and (B$n$.$m$) are equations of primer v2. Fermionic and bosonic Matsubara frequencies are $\omega^{\mathrm{F}}_p$ and $\omega^{\mathrm{B}}_p$ ([§A.7](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#a7-working)), and $\mathrm{i}\omega$ is a point on the imaginary axis, as in (A9.20).

EPP as written cannot produce any ordered phase: it is a normal-state theory of three bosonic propagators, and its static lattice self-energy only stabilises. It is, however, a natural language for pairing in polar and dilute superconductors, and it already enforces the gauge constraint against uniform photon condensation.

| Phase | Instability of | EPP has | EPP lacks |
|---|---|---|---|
| Superconductivity | Cooper channel of $G$, interaction $D$ | polar part of $D=D^{\mathrm{e}}+D^{\mathrm{ph}}$, dynamically screened | $G$ with anomalous (Nambu–Gor'kov) parts; short-range $g$ |
| CDW, ferroelectric | $F^{-1}=0$: $\bar\omega^2_{\mathbf{q}\nu}\to0$ | $F=F^0+F^0\Pi^{\mathrm{ph}}F$ (E.398) | short-range electron–lattice coupling, the $\mathbf{G}\neq0$ dipole sum, anharmonicity |
| Photon or polariton condensate | $D^{-1}(\mathbf{q},0)=0$, or driven occupation of a dispersion minimum | $\Pi(\omega)-\Pi(0)$ (E.281, E.311), which enforces the uniform no-go theorem but also sets $\chi_{\mathrm{orb}}\equiv0$ | static transverse (orbital) response, excitons (BSE), cavity $D^0$, kinetics |

**Verdict.**
- **Superconductivity.** With a Nambu–Gor'kov extension (4.1), EPP could predict $T_{\mathrm{c}}$ for polar pairing with dynamical screening and LO–plasmon modes ([plasmon.md](plasmon.md)). That problem has been solved in model form since Takada (1980). EPP would add first-principles $Z^*$, $\varepsilon^{\mathrm{e}}(\mathbf{q},\mathrm{i}\omega)$ and anisotropy, all from one propagator. Best tests: doped SrTiO$_3$ and FeSe/SrTiO$_3$. Retarded transverse photons are irrelevant in bulk [derived here], and cavity enhancement is speculative.
- **CDW.** With the $\mathbf{G}\neq0$ dipole sum (E.454a; (B9.15)), EPP gives Cochran (ferroelectric) softening. With DFPT short-range force constants as input, however, this only reproduces DFPT. $2k_{\mathrm{F}}$-type CDWs (NbSe$_2$, TaS$_2$) need short-range coupling and anharmonicity, which EPP does not improve, and TiSe$_2$ needs excitons as well.
- **Condensates.** Phonons and photons have no equilibrium BEC. Uniform photon condensation is forbidden by gauge invariance, and a finite-$\mathbf{q}$ condensate needs a static orbital response that EPP's $\Pi(\omega)-\Pi(0)$ sets to zero (4.5). Polariton condensates need BSE, a cavity and kinetics; EPP-type theory would supply the dispersion and the relaxation vertices.

## 1. What is missing in EPP

- **No fermion propagator.** EPP solves three Dyson equations:
  - $D=D^0+D^0\Pi D$, with $D^0=D^{0\mathrm{L}}+D^{0\mathrm{T}}$ (E.375, E.379–380; (B6.11), (B6.14));
  - $D^{\mathrm{e}}$ (E.388; (B5.13));
  - $F$, with $\Pi^{\mathrm{ph}}=-J^*D^{\mathrm{e}}J$ (E.398–399; (B7.5), (B7.7)).

  Electrons appear only as the RPA bubble (E.282–283; [§B.5](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b-photon-se)–[§B.7](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b-phonon-se)). There is no $G$, no $\Sigma$, no pair amplitude $\langle\psi_\uparrow\psi_\downarrow\rangle$, hence no gap equation.
- **Only the polar channel.** Ions couple only through $\frac1c\dot{\mathbf{p}}_{\ell\kappa}\cdot\mathbf{A}$ with $\mathbf{p}_{\ell\kappa}=e\mathbf{Z}^*_\kappa\cdot\mathbf{u}_{\ell\kappa}$ (E.134–135; (B2.6), (B2.10)), on a frozen crystal potential (E.136; (B2.7)) and a harmonic lattice. Deformation-potential coupling, which dominates $\lambda_{\mathrm{ep}}$ in conventional superconductors and in CDW metals, is absent ([§B.8](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b-frohlich), [§X.2](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#x2-missing)). The bare vertex $g^0$ appears in EPP Sec. 2 (E.100–109, with the double-counting caveat of pp. 20–21; primer [§B.1](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b1-bare)) but is dropped from E.134. EPP does keep harmonic short-range force constants ($\Phi^{\mathrm{sr}}$, EPP's $k^{\mathcal{S}}$ of p. 21). What it lacks is short-range electron–lattice coupling and anharmonicity.
- **The static lattice self-energy only hardens; no zero-frequency root [derived here].** At $c\to\infty$, $\Pi^{\mathrm{ph}}_{\nu\nu}=\hbar W^{\mathrm{e}}|\mathbf{q}\cdot\mathbf{j}_\nu|^2/2e^2\Omega_{\mathrm{V}}\omega_{\mathbf{q}\nu}$ (E.452; (B9.8)). In RPA, $\chi^0(\mathbf{q},0)\le0$, hence $\varepsilon^{\mathrm{e}}(\mathbf{q},0)\ge1$ and $W^{\mathrm{e}}(\mathbf{q},0)>0$, so E.454 ((B9.10)) can never drive $\omega^2$ to zero. (The sign also holds with crystal local fields, because the static $-\chi^0_{\mathbf{G}\mathbf{G}'}$ is positive semidefinite.) This sign argument is specific to RPA. Beyond it, with exchange–correlation, $\varepsilon(\mathbf{q}\neq0,0)<0$ is admissible [Dolgov, Kirzhnits, Maksimov, Rev. Mod. Phys. **53**, 81 (1981)]. Dynamically, E.453 ((B9.9)) gives an $\mathrm{L}_-$ branch below $\omega_{\mathbf{q}\nu}$ in doped crystals (worked out in [§B.9](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b9-dynamic)), but a zero-frequency root would need $W^{\mathrm{e}}(\mathbf{q},0)<0$, which RPA forbids. Cochran softening, $\omega^2_{\mathrm{TO}}=\omega_0^2-(\omega^{\mathrm{ion}}_{\mathrm{p}})^2/3$ (E.160; (B3.14) with $\zeta=1/3$), comes from the Lorentz field (E.151; (B3.13)) of the classical model. The Dyson route lacks it because E.380–399 keep only $\mathbf{G}=0$ in the ionic vertex and in $D^0$ (macroscopic field only). The $\mathbf{G}\neq0$ dipole sum (E.454a, Gonze–Lee; (B9.15)) restores the Lorentz/Ewald term. Electronic local fields ("CLFs are mostly neglected", EPP p. 95) are a separate omission: they change $\varepsilon^{\mathrm{e}}$, while $Z^*$, an input, already contains them ([§B.3](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b3-which)).
- **Temperature and ladders.** $T$ enters only through $f_{n\mathbf{k}}$ in $\Pi^{\mathrm{e}}$ ($n_{\mathrm{B}}=0$ for bosons, E.246–250). There is no anharmonicity, and no particle–hole (exciton) or particle–particle (Cooper) ladder.

## 2. Why it matters

- **Light cone.** Retardation acts only for $q\lesssim q_{\mathrm{light}}(\omega)\sim10^{-5}$–$10^{-4}$ Å$^{-1}$ ([§X.1](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#x1-cone)), while pairing samples $2k_{\mathrm{F}}\sim0.2$–$2$ Å$^{-1}$. At larger $q$ the transverse (magnetic, Breit/Darwin) interaction persists at relative order $(v_{\mathrm{F}}/c)^2$, about $10^{-5}$ of the Coulomb term for $v_{\mathrm{F}}=10^6$ m/s. The Fermi-surface fraction with $|\mathbf{k}-\mathbf{k}'|<q_{\mathrm{light}}(\omega)$ is $(q_{\mathrm{light}}(\omega)/2k_{\mathrm{F}})^2\approx2.5\times10^{-11}$–$2.5\times10^{-7}$ in 3D, and $q_{\mathrm{light}}(\omega)/(\pi k_{\mathrm{F}})\approx3\times10^{-6}$–$3\times10^{-4}$ for a 2D Fermi line [derived here]. For SrTiO$_3$ at $5.5\times10^{17}$ cm$^{-3}$ (one spin-degenerate band), $2k_{\mathrm{F}}\approx0.05$ Å$^{-1}$, $v_{\mathrm{F}}\approx1.3$–$2\times10^4$ m/s and $(v_{\mathrm{F}}/c)^2\approx2$–$4\times10^{-9}$ [derived here].
- **Dilute polar superconductors.** SrTiO$_3$ superconducts down to $5.5\times10^{17}$ cm$^{-3}$, where $E_{\mathrm{F}}\approx1.1$ meV (Lin 2013), far below its polar phonons. Migdal's condition $\hbar\omega_{\mathrm{ph}}/E_{\mathrm{F}}\ll1$ fails there. So does the separation of scales behind $\mu^*$ (G.204), because the plasma frequency falls among the phonons. The mechanism is still debated (Gastiasoro 2020): LO phonons (Gor'kov 2016), plasmons (Ruhman and Lee 2016; see [plasmon.md](plasmon.md)) and soft TO modes have all been proposed.
- **Interfacial Fröhlich coupling.** Monolayer FeSe on SrTiO$_3$ opens gaps near 70 K (bulk FeSe: about 8 K). ARPES shows replica bands about 100 meV below the main bands, attributed to small-$q$ coupling to substrate oxygen phonons (Lee 2014). Estimates of the coupling are $\lambda_{\mathrm{ep}}\approx0.5$ (attributed to Lee 2014; not confirmed in the primary text) and $0.19\pm0.02$ (Faeth 2021). An extrinsic origin of the replicas, photoelectron energy loss to Fuchs–Kliewer phonons, has been argued (Li and Sawatzky, preprint), and Faeth et al. argue against it. In model RPA the FeSe electrons screen the dipolar coupling almost completely (Zhou and Millis 2017).
- **Cavities.** Models give a pair-density-wave state with $T_{\mathrm{c}}$ in the low-kelvin range for a GaAs-like 2D gas (Schlawin 2019; about 0.5 K in the preprint). They give no enhancement of forward-scattering pairing, because mode softening offsets the coupling gain (Sentef 2018). A reported rise of Rb$_3$C$_{60}$ from 30 to 45 K (Thomas et al. 2025, figures as given in the preprint) is debated.
- **CDW.** In NbSe$_2$, $T_{\mathrm{CDW}}$ from SSCHA is about 59 K in bulk and 73 K in the monolayer, with ionic fluctuations dominating (Bianco 2020). In TiSe$_2$ an electronic mode softens to zero at finite $q$ near 190 K (Kogar 2017). The excitonic interpretation is debated (Lian, Ali and Wong, preprint).
- **Condensates.** Polariton condensation has been seen at 19 K in CdTe (Kasprzak 2006) and in monolayer MoSe$_2$ (Anton-Solanas 2021). Signals of interlayer-exciton condensation persist above 100 K in MoSe$_2$/WSe$_2$ (Wang 2019).

**For an RPA/GW/BSE practitioner.**
1. The Coulomb part of pairing, $V_{n\mathbf{k},m\mathbf{k}'}$ in G.209, is a Fermi-surface matrix element of your $W$. Use the full $W_{\mathbf{G}\mathbf{G}'}(\mathbf{q},\mathrm{i}\omega)$ of your code, not EPP's head-only $W^{\mathrm{e}}$: the matrix elements involve $q$ up to $2k_{\mathrm{F}}$ and umklapp terms. Add $W^{\mathrm{ph}}$ to the head ($\mathbf{G}=\mathbf{G}'=0$) only.
2. In dilute or 2D polar systems $W^{\mathrm{e}}$ and $W^{\mathrm{ph}}$ cannot be separated by retardation, and EPP builds both from one propagator.
3. The substrate sets, through $\varepsilon_{\mathrm{env}}$, both the exciton screening length $r_0/\varepsilon_{\mathrm{env}}$ ([§B.10](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b10-keldysh)) and the interfacial pairing.
4. Excitonic insulators and interlayer exciton condensates are BSE problems with an anomalous electron–hole self-energy.

## 3. What standard theory does

- **Pairing [G §XI].** McMillan–Allen–Dynes (G.204–207) and anisotropic Migdal–Eliashberg theory in Nambu–Gor'kov form (G.208–210), Wannier-interpolated (Margine 2013). Here $\lambda(\mathrm{i}\omega)$ is $W^{\mathrm{ph}}=W^{\mathrm{e}}\chi^{\mathrm{ion}}W^{\mathrm{e}}$ projected on pairs ([§A.9](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#a9-eliashberg); (A9.20)), and $\mu^*$ is either empirical or taken from the RPA $W$. SCDFT (G.211; Lüders 2005) builds $\mathcal{K}^{\mathrm{ee}}$ from the screened Coulomb interaction; Akashi and Arita (2013) made that kernel dynamical, so plasmons pair too.
- **Dilute polar pairing [lit].**
  - Takada (1980) solved the gap equation for n-SrTiO$_3$ with a dynamically screened polar coupling and the plasmon–LO modes. He used the dielectric formulation of pairing of Kirzhnits, Maksimov and Khomskii.
  - Gor'kov (2016) argued for nonadiabatic LO-phonon pairing, and Ruhman and Lee (2016) for plasmon pairing.
  - Gastiasoro, Chubukov and Fernandes (2019) solved phonon pairing at $E_{\mathrm{F}}\lesssim\hbar\omega_{\mathrm{ph}}$ with Coulomb repulsion, both dressed by the electron polarizability.
  - An interaction that is repulsive at every Matsubara frequency can still pair, with a gap $\Delta(\mathrm{i}\omega)$ that changes sign (Pimenov and Chubukov 2022).
- **Lattice instabilities [G §V.A, §VII].**
  - The DFPT force constants $\Pi^{\mathrm{ph,A}}$ (G.128) already contain Kohn anomalies and the imaginary frequencies of CDWs, and $\Pi^{\mathrm{ph,NA}}$ (G.141, G.145) adds the dynamics ([§A.6](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#a-phonon-se)).
  - Nesting predicts CDWs poorly: the $\mathbf{q}$ dependence of the coupling decides (Johannes and Mazin 2008).
  - The distorted phase and $T_{\mathrm{CDW}}$ need anharmonicity (SSCHA; Bianco 2020).
  - Ferroelectric soft modes balance dipolar against short-range forces (Cochran 1960). They are computed with DFPT plus the Gonze–Lee dipole sum (cf. E.454a; (B9.15)).
- **Condensates.**
  - The $A^2$ term forbids the superradiant transition (Rzążewski 1975; Nataf and Ciuti 2010).
  - Gauge invariance forbids uniform photon condensation (Andolina 2019). A finite-$q$ condensate requires $\chi_{\mathrm{orb}}(\mathbf{q})>1/4\pi$ in 3D (Andolina 2020).
  - Photon BEC needed a thermalizing dye (Klaers 2010), and polariton condensates are driven and dissipative.

## 4. How it could enter the EPP formalism

**4.1 Pairing kernel from $D$ [proposal].** With Nambu spinors $\Psi_{\mathbf{k}}=(c_{\mathbf{k}\uparrow},c^\dagger_{-\mathbf{k}\downarrow})^{\mathsf{T}}$ and Pauli matrices $\hat\tau_i$ in Nambu space, EPP's vertices become:
- the paramagnetic current $-\frac1c\mathbf{j}\cdot\mathbf{A}\,\hat\tau_0$;
- the diamagnetic term $\frac{e^2}{2m_{\mathrm{e}}c^2}n^{\mathrm{e}}A^2\,\hat\tau_3$;
- the charge $-e\rho\,\hat\tau_3$, via continuity (E.428′; (X2.8)). E.428′ holds for the normal-state $H^{\mathrm{e}}_0$ (p. 94), so use it only in the linearized gap equation ($T\ge T_{\mathrm{c}}$). Below $T_{\mathrm{c}}$, work in Coulomb gauge from the start, with the scalar potential coupling through $\hat\tau_3$.

In Coulomb-gauge form, which is exact for $c\to\infty$ (E.441–445), the self-energy is (schematic; prefactors to be fixed against G.208–210)

$$\hat\Sigma(\mathbf{k},\mathrm{i}\omega^{\mathrm{F}}_p)=-\frac{k_{\mathrm{B}}T}{N_{\mathbf{k}}}\sum_{\mathbf{k}'p'}\Big[\hat\tau_3\hat G(\mathbf{k}',\mathrm{i}\omega^{\mathrm{F}}_{p'})\hat\tau_3\,|\rho_{\mathbf{k}\mathbf{k}'}|^2W(\mathbf{q},\mathrm{i}\omega^{\mathrm{B}}_{p-p'})+\hat\tau_0\hat G(\mathbf{k}',\mathrm{i}\omega^{\mathrm{F}}_{p'})\hat\tau_0\,K^{\mathrm{T}}_{\mathbf{k}\mathbf{k}'}(\mathrm{i}\omega^{\mathrm{B}}_{p-p'})\Big],$$

$$W=\frac{4\pi e^2}{\mathbf{q}\cdot\boldsymbol{\varepsilon}\cdot\mathbf{q}},\qquad\boldsymbol{\varepsilon}=\boldsymbol{\varepsilon}^{\mathrm{e}}+\boldsymbol{\varepsilon}^{\mathrm{ion}}-1,\qquad W=W^{\mathrm{e}}+W^{\mathrm{ph}}$$

(head-only form; the Coulomb part needs $W_{\mathbf{G}\mathbf{G}'}$).

Here $\mathbf{q}=\mathbf{k}'-\mathbf{k}$, $\omega^{\mathrm{B}}_{p-p'}=\omega^{\mathrm{F}}_p-\omega^{\mathrm{F}}_{p'}$, and $\rho_{\mathbf{k}\mathbf{k}'}$ is the number-density matrix element (the charge vertex is $-e\rho$). The additivity of $\boldsymbol{\varepsilon}$ is EPP p. 33 (cf. E.432, E.445). The transverse kernel $K^{\mathrm{T}}=\frac1c\mathbf{j}^*\cdot D^{\mathrm{T}}\cdot\mathbf{j}$ is of relative order $(v_{\mathrm{F}}/c)^2$. Subtract $D^{0\mathrm{T}}$ only in the normal part, where it is a mass renormalization; in the anomalous part the vacuum transverse exchange is the physical Breit/Darwin interaction. For time-reversed pairs that interaction is repulsive (Ampère: antiparallel currents repel).

For one isotropic mode, $W^{\mathrm{ph}}(\mathbf{q},\mathrm{i}\omega)=-\frac{4\pi e^2}{q^2}\frac{(\omega^{\mathrm{ion}}_{\mathrm{p}})^2/\varepsilon_\infty^2}{\omega^2+\bar\omega_{\mathrm{LO}}^2}$, which is attractive [derived here; cf. self-energy topic §4.4]. With $\varepsilon_\infty\to\varepsilon^{\mathrm{e}}(\mathbf{q},\mathrm{i}\omega)$ the formula holds pointwise in $\mathrm{i}\omega$, with $\bar\omega^2_{\mathrm{LO}}(\mathbf{q},\mathrm{i}\omega)=\omega^2_{\mathrm{TO}}+(\omega^{\mathrm{ion}}_{\mathrm{p}})^2/\varepsilon^{\mathrm{e}}(\mathbf{q},\mathrm{i}\omega)$. On the real axis the poles of the full $W$ are the $\mathrm{L}_\pm$ roots of E.453, which primer v2 now derives for Drude carriers ((B9.18); [plasmon.md](plasmon.md)). For pairing, use $W(\mathbf{q},\mathrm{i}\omega)=4\pi e^2/\{q^2[\varepsilon^{\mathrm{e}}(\mathbf{q},\mathrm{i}\omega)+(\omega^{\mathrm{ion}}_{\mathrm{p}})^2/(\omega^2+\omega^2_{\mathrm{TO}})]\}$ without splitting. Linearized, G.209 holds with

$$\lambda_{n\mathbf{k},m\mathbf{k}'}(\mathrm{i}\omega)-N_{\mathrm{F}}V_{n\mathbf{k},m\mathbf{k}'}\;\to\;-N_{\mathrm{F}}\langle n\mathbf{k},m\mathbf{k}'|W^{\mathrm{e}}+W^{\mathrm{ph}}|n\mathbf{k},m\mathbf{k}'\rangle(\mathrm{i}\omega).$$

This is the dielectric formulation of pairing (Kirzhnits, Maksimov and Khomskii; journal reference unverified). Takada (1980) applied it to polar n-SrTiO$_3$ with the plasmon and the LO phonon.
- **New relative to G.209 (but known in model form since Takada 1980):** a dynamically screened Fröhlich vertex (E.448a; (B8.16)–(B8.17)), $\omega_{\mathrm{LO}}$ as an output, and coupled LO–plasmon modes. EPP's own contribution would be first-principles $Z^*$, $\varepsilon^{\mathrm{e}}(\mathbf{q},\mathrm{i}\omega)$ and anisotropy.
- **Still missing:** the short-range $g$, and vertex corrections for $\hbar\omega_{\mathrm{LO}}/E_{\mathrm{F}}\gtrsim1$.

**4.2 2D and cavities [proposal].**
- **Substrates.** For a sheet on a polar substrate, replace the single-sheet $W^{\mathrm{e}}(\mathbf{Q}_\parallel,\omega;z,z')$ (E.461; (B10.11)) by a layered Green's function with $\varepsilon_{\mathrm{sub}}(\omega)=\varepsilon_\infty+\sum$ (ionic poles from the substrate's $\Pi^{\mathrm{ion}}$). This is a new ingredient: EPP's E.455–472 cover only a freestanding sheet with $\delta(z)$ ions, while the substrate is a half-space. Primer v2 adds only a dielectric environment, $\varepsilon_{\mathrm{env}}=\frac12(\varepsilon_{\mathrm{sub}}+\varepsilon_{\mathrm{cap}})$, taken at the frequency of the process ((B10.8)). The FeSe/SrTiO$_3$ coupling is then a Fuchs–Kliewer surface-phonon pole of $W(\mathbf{Q}_\parallel,\mathrm{i}\omega;0,0)$, small-$Q_\parallel$ by construction.
- **Cavities.** Replace $D^0$ by the cavity Green's dyadic $D^0_{\mathrm{cav}}$; E.379 and E.398 keep their form.
  - Sentef's mechanism is $\Pi^{\mathrm{ph}}$ and $g^{\mathrm{dyn}}$ (E.399, E.416; (B7.5), (B8.6)) evaluated with $D^{\mathrm{e}}_{\mathrm{cav}}$.
  - Schlawin's is $K^{\mathrm{T}}[D_{\mathrm{cav}}]$ at finite pair momentum, which gives a pair-density wave.
  - Curtis's driven enhancement needs nonequilibrium kinetics, outside EPP.

**4.3 Meissner kernel [derived here; standard].** With Nambu propagators in $\Pi^{\mathrm{e}}$, the identity $\Pi^{\mathrm{dia}}=-\Pi^{\mathrm{para}}(\omega=0)$ behind E.281 ((B5.15)) fails: $\Pi^{\mathrm{T}}(\mathbf{q}\to0,0)\propto n_{\mathrm{s}}/m^*\neq0$ (the London kernel). Superconductivity then appears in EPP's photon sector as a photon mass, and E.279 must replace E.281. This applies to the transverse sector. In the longitudinal sector the bare Nambu bubble violates $\mathbf{q}\cdot\Pi(\mathbf{q},0)\cdot\mathbf{q}=0$, which gives a spurious longitudinal photon mass. The vertex must satisfy Nambu's Ward identity [Y. Nambu, Phys. Rev. **117**, 648 (1960), doi:10.1103/PhysRev.117.648]. That brings in the phase mode, which $D^{0\mathrm{L}}$ pushes up to the plasma frequency $\bar\omega_{\mathrm{p}}$.

**4.4 Soft modes [proposal].** Reinstate EPP's own $g^0$ (E.107; (B1.4)) in E.134, minus its long-range $Z_\kappa/|\mathbf{q}+\mathbf{G}|^2$ part, which duplicates the dipole channel of E.134. This gives a bare short-range vertex $g^{\mathrm{sr},0}$, with screened partner $g^{\mathrm{sr}}$, and the lattice Dyson equation becomes

$$F=F^0+F^0\big(\Pi^{\mathrm{ph,lr}}+\Pi^{\mathrm{ph,sr}}+\Pi^{\mathrm{ph,lr\text{-}sr}}+\Pi^{\mathrm{ph,sr\text{-}lr}}\big)F,\qquad\Pi^{\mathrm{ph,sr}}_{\nu\nu'}=\sum_{\mathbf{k}nm}\big[g^{\mathrm{sr},0}_{mn\nu}\big]^*\,\chi^0_{n\mathbf{k},m\mathbf{k}+\mathbf{q}}(\omega)\,g^{\mathrm{sr}}_{mn\nu'}(\omega).$$

Here $\Pi^{\mathrm{ph,lr}}$ is EPP's polar self-energy (E.399; (B7.5)). As in E.412 and G.141 ((A6.17)), one vertex is bare and one screened; the cross terms $\Pi^{\mathrm{ph,lr\text{-}sr}}$ and $\Pi^{\mathrm{ph,sr\text{-}lr}}$ mix the polar and nonpolar channels. Only two partitions are consistent, and mixing them counts the valence electrons twice or three times:
- **(A) Screened/adiabatic (EPP's own choice, p. 21).** $F^0$ comes from DFPT or SSCHA short-range force constants, with the Gonze–Lee dipole–dipole part $\Phi^{\mathrm{dd}}$ removed. $J$ is built with $Z^*$ and $\varepsilon_\infty$ (E.454a; (B9.15)). The short-range channel enters only through its nonadiabatic remainder $\Pi^{\mathrm{ph,sr}}(\omega)-\Pi^{\mathrm{ph,sr}}(0)$ (the G.145 structure), with $g^{\mathrm{sr}}=g^{\mathrm{DFPT}}-g^{\mathrm{lr}}$. The CDW instability is then an input through $F^0$. EPP adds the dynamics and the polar or substrate channel.
- **(B) Bare (the structure of G.119, whose nonadiabatic part is G.141).** $F^0$ comes from the bare ion–ion plus frozen-density terms. $J$ is built with the bare valence charges $Z_\kappa$, not $Z^*$, and $g^{\mathrm{sr},0}$ is E.107 with its $Z_\kappa/|\mathbf{q}+\mathbf{G}|^2$ tail removed. $Z^*$, $\varepsilon_\infty$ and the Kohn anomalies must then come out of $\chi^0$ with local fields. SSCHA cannot be combined with (B), because it works on the Born–Oppenheimer surface, which already contains the electronic screening. Giustino notes that bare force constants give dispersions very different from the real ones (p. 015003-21), so (B) is a long project.

The instability criterion, with $\Pi^{\mathrm{ph}}$ now the sum of the four terms, is

$$\bar\omega^2_{\mathbf{q}\nu}(T)=\omega^2_{\mathbf{q}\nu}(T)+\frac{2\omega_{\mathbf{q}\nu}}{\hbar}\,\mathrm{Re}\,\Pi^{\mathrm{ph}}_{\nu\nu}(\mathbf{q},\bar\omega_{\mathbf{q}\nu};T)=0,\qquad\mathbf{q}=\mathbf{q}_{\mathrm{CDW}}.$$

Here $T$ enters through $f_{n\mathbf{k}}$ in $\chi^0$ and, in partition (A) only, through an anharmonic $\omega^2_{\mathbf{q}\nu}(T)$ from SSCHA. Keeping $\mathbf{G}\neq0$ in $D^{\mathrm{e}}$ (E.454a; (B9.15), with the on-site Ewald term removed) restores the Lorentz field, and with it Cochran-type softening.

**4.5 Condensates [assessment].**
- **Phonons.** A "phonon condensate" $\langle Q_{\mathbf{q}\nu}\rangle\neq0$ at $\mathbf{q}=\mathbf{q}_{\mathrm{CDW}}$ is just the frozen distortion of 4.4.
- **Photons.** Photon condensation needs a zero transverse eigenvalue of $D^{-1}(\mathbf{q},0)=[D^0]^{-1}-\Pi(\mathbf{q},0)$. E.281 and E.311 set $\Pi(\mathbf{q},0)=0$ at every $\mathbf{q}$ ((B5.19), (B5.26)). This enforces the uniform no-go theorem, but it also deletes the static orbital response at finite $\mathbf{q}$, as primer v2 notes in [§B.5](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b5-truncation). For the point-dipole ions of E.311 that is exact (E.309–310; (B5.25)); for electrons it is not. A finite-$\mathbf{q}$ condensate requires $\chi_{\mathrm{orb}}(\mathbf{q})>1/4\pi$ in 3D, Gaussian units (Andolina et al. 2020). To test this in EPP, keep E.279 (explicit dia + para, complete band sums) for the transverse part, and use E.281 only for the longitudinal part.
- **Polaritons.** Put $\Pi^{\mathrm{e,BSE}}$ (long-range exchange removed, as in [optics.md](optics.md) §4.1) and $D^0_{\mathrm{cav}}$ into E.379. The poles of $D$ then give the polariton dispersion, the Hopfield weights and the radiative widths. Relaxation (the bottleneck) needs exciton–phonon vertices, and condensation itself needs a kinetic or Gross–Pitaevskii layer.

## 5. Pitfalls

- **Double counting.**
  - Never combine $Z^*$ or $\varepsilon_\infty$ with an $F^0$ that has no electron contribution, or DFPT/SSCHA $\omega^2$ with a full static $\Pi^{\mathrm{ph,sr}}$. Use partition (A) or (B) of 4.4.
  - Subtract $g^{\mathrm{lr}}$ from the DFPT $g$.
  - Add no empirical $\mu^*$ on top of $W^{\mathrm{e}}$.
  - Coulomb is $D^{0\mathrm{L}}$ (E.364; (B6.10)); never add $v$ separately.
  - Normal self-energy: either start from bare (DFT) bands and keep the full $W^{\mathrm{e}}+W^{\mathrm{ph}}$ normal part, or start from $GW$ bands and keep only $W^{\mathrm{ph}}$ and the anomalous part of $W^{\mathrm{e}}$. Drop the $W^{\mathrm{ph}}$ normal part too if the input bands already contain the Fröhlich mass renormalization.
- **Gauge.** $\Pi(\omega)-\Pi(0)$ is exact for the longitudinal part at all $q$, and for the transverse part only as $q\to0$. At finite transverse $q$, or in the paired state, use dia + para with band sums complete enough to satisfy the f-sum rule. In the paired state the longitudinal vertex must also satisfy Nambu's Ward identity (4.3). Truncated band sums break the f-sum rule and fake a Drude weight or a $q=0$ photon instability.
- **Sum rules.** Impose the acoustic sum rules for $Z^*$ and for the force constants.
- **Migdal.** Include vertex corrections when $\hbar\omega_{\mathrm{ph}}/E_{\mathrm{F}}\gtrsim1$.
- **Convergence.**
  - Fermi-surface integrals need up to $10^6$ $\mathbf{k}$ points.
  - Small-$q$ forward scattering needs analytic integration.
  - Matsubara cutoffs renormalize $\mu_{\mathrm{C}}$ to $\mu^*$ ((A9.19)). At $T_{\mathrm{c}}\sim0.1$ K with $\omega^{\mathrm{B}}_p$ up to several $\omega_{\mathrm{LO}}$, use an intermediate-representation (sparse) Matsubara basis.
  - $\mathbf{q}_{\mathrm{CDW}}$ must lie on the grid, and smearing is not the physical $T$.
- **Regularization.** The $1/\omega^2$ of $D^{0\mathrm{L}}$ and the noncommuting limits $\mathbf{q},\omega\to0$ need care ([§B.9](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b-electrostatic)). Subtract $D^{0\mathrm{T}}$ only in the normal part (4.1).

## 6. State of the art (verified)

1. F. Giustino, Rev. Mod. Phys. **89**, 015003 (2017), doi:10.1103/RevModPhys.89.015003. Eliashberg, SCDFT, nonadiabatic phonons.
2. E. R. Margine, F. Giustino, Phys. Rev. B **87**, 024505 (2013), doi:10.1103/PhysRevB.87.024505. Anisotropic Eliashberg with Wannier interpolation.
3. M. Lüders *et al.*, Phys. Rev. B **72**, 024545 (2005), doi:10.1103/PhysRevB.72.024545. SCDFT.
4. R. Akashi, R. Arita, Phys. Rev. Lett. **111**, 057006 (2013), doi:10.1103/PhysRevLett.111.057006. Dynamical Coulomb kernel; plasmon-assisted pairing.
5. X. Lin *et al.*, Phys. Rev. X **3**, 021002 (2013), doi:10.1103/PhysRevX.3.021002. SrTiO$_3$ superconducting at $E_{\mathrm{F}}\approx1.1$ meV.
6. J. J. Lee *et al.*, Nature **515**, 245 (2014), doi:10.1038/nature13894. FeSe/SrTiO$_3$ replica bands.
7. M. A. Sentef, M. Ruggenthaler, A. Rubio, Sci. Adv. **4**, eaau6969 (2018), doi:10.1126/sciadv.aau6969. Cavity phonon polaritons and pairing.
8. F. Schlawin, A. Cavalleri, D. Jaksch, Phys. Rev. Lett. **122**, 133602 (2019), doi:10.1103/PhysRevLett.122.133602 (DOI by APS pattern). Cavity-photon pairing; a pair-density wave with $T_{\mathrm{c}}$ in the low-kelvin range (about 0.5 K in arXiv:1804.07142v1).
9. J. B. Curtis *et al.*, Phys. Rev. Lett. **122**, 167002 (2019), doi:10.1103/PhysRevLett.122.167002. Driven cavity Eliashberg effect.
10. R. Bianco *et al.*, Phys. Rev. Lett. **125**, 106101 (2020), doi:10.1103/PhysRevLett.125.106101. Anharmonic CDW in NbSe$_2$.
11. A. Kogar *et al.*, Science **358**, 1314 (2017), doi:10.1126/science.aam6432 (DOI from secondary sources only). Exciton condensation in TiSe$_2$; interpretation debated.
12. G. M. Andolina *et al.*, Phys. Rev. B **100**, 121109(R) (2019), doi:10.1103/PhysRevB.100.121109. No-go theorem.
13. G. M. Andolina *et al.*, Phys. Rev. B **102**, 125137 (2020), doi:10.1103/PhysRevB.102.125137 (DOI from the publisher's file name). Finite-$q$ condensation; criterion $\chi_{\mathrm{orb}}(\mathbf{q})>1/4\pi$ in 3D.
14. J. Kasprzak *et al.*, Nature **443**, 409 (2006), doi:10.1038/nature05131 (DOI from the RePEc record). Polariton BEC.
15. C. Anton-Solanas *et al.*, Nat. Mater. **20**, 1233 (2021), doi:10.1038/s41563-021-01000-8. Polariton condensate in MoSe$_2$.
16. Y. Nambu, Phys. Rev. **117**, 648 (1960), doi:10.1103/PhysRev.117.648. Ward identity and a gauge-invariant Meissner kernel in the paired state.
17. Y. Takada, J. Phys. Soc. Jpn. **49**, 1267 (1980) (DOI not confirmed). Polar n-SrTiO$_3$ with plasmon and LO phonon in the dielectric formulation.
18. M. N. Gastiasoro, A. V. Chubukov, R. M. Fernandes, Phys. Rev. B **99**, 094524 (2019), doi:10.1103/PhysRevB.99.094524. Phonon pairing at $E_{\mathrm{F}}\lesssim\hbar\omega_{\mathrm{ph}}$ with Coulomb repulsion.
19. J. Ruhman, P. A. Lee, Phys. Rev. B **94**, 224515 (2016), doi:10.1103/PhysRevB.94.224515. Plasmon pairing in dilute SrTiO$_3$.
20. L. P. Gor'kov, Proc. Natl. Acad. Sci. USA **113**, 4646 (2016), doi:10.1073/pnas.1604145113. Nonadiabatic LO-phonon pairing in dilute SrTiO$_3$.
21. D. Pimenov, A. V. Chubukov, npj Quantum Mater. **7**, 45 (2022), "Quantum phase transition in a clean superconductor with repulsive dynamical interaction" (DOI not confirmed). Pairing from a dynamically repulsive interaction.
22. Y. Zhou, A. J. Millis, Phys. Rev. B **96**, 054516 (2017), doi:10.1103/PhysRevB.96.054516; arXiv:1703.04021. RPA screening of the SrTiO$_3$ dipolar coupling by the FeSe electrons.
23. G. Marchese *et al.*, Nat. Phys. **20**, 88 (2024), doi:10.1038/s41567-023-02203-3. Dynamical Born charges in metals.
24. C. Faeth *et al.*, Phys. Rev. Lett. **127**, 016803 (2021), doi:10.1103/PhysRevLett.127.016803. Intrinsic replica bands, $\lambda_{\mathrm{ep}}=0.19\pm0.02$.

Also cited:
- O. V. Dolgov, D. A. Kirzhnits, E. G. Maksimov, Rev. Mod. Phys. **53**, 81 (1981), doi:10.1103/RevModPhys.53.81. Admissible sign of $\varepsilon(\mathbf{q},0)$.
- Cochran, Adv. Phys. **9**, 387 (1960), doi:10.1080/00018736000101229.
- Nataf, Ciuti, Nat. Commun. **1**, 72 (2010), doi:10.1038/ncomms1069.
- Klaers *et al.*, Nature **468**, 545 (2010), doi:10.1038/nature09567.
- Z. Wang *et al.*, Nature **574**, 76 (2019), doi:10.1038/s41586-019-1591-7.
- Rzążewski, Wódkiewicz, Żakowicz, Phys. Rev. Lett. **35**, 432 (1975); DOI not confirmed.
- Gastiasoro, Ruhman, Fernandes, Ann. Phys. 168107 (2020), arXiv:1912.01509; volume and DOI unverified.
- Johannes, Mazin, arXiv:0708.1744; journal reference (Phys. Rev. B **77**, 165135, 2008) seen only in a citing paper, unverified.
- A. Thomas *et al.*, J. Chem. Phys. **162** (13) (2025), doi:10.1063/5.0231202 (article number unverified); arXiv:1911.01459. Rb$_3$C$_{60}$ 30→45 K as reported in the preprint.
- F. Li, G. A. Sawatzky, arXiv:1710.10795 (journal reference unverified). Extrinsic (energy-loss) origin of the FeSe/SrTiO$_3$ replicas.
- C. Lian, Z. A. Ali, B. M. Wong, arXiv:1901.11223 (journal reference unverified). Linear-response TDDFT reading of the TiSe$_2$ soft mode.
- D. A. Kirzhnits, E. G. Maksimov, D. I. Khomskii, Lebedev Physics Institute preprint No. 108 (1970), as cited in a reference list; journal reference unverified. Dielectric formulation of pairing.

## 7. First steps

1. **Derive (1–2 weeks).** Write EPP's vertices in Nambu form. Show that 4.1 reduces to G.208–210 with $|g|^2\to|g^{\mathrm{lr}}|^2$ (E.448a; (B8.17)) and $W^{\mathrm{ph}}=W^{\mathrm{e}}\chi^{\mathrm{ion}}W^{\mathrm{e}}$. Show that integrating out $|\omega^{\mathrm{B}}_p|>\omega_{\mathrm{cut}}$ reproduces Morel–Anderson ((A9.19)) only for $\bar\omega_{\mathrm{p}}\gg\omega_{\mathrm{LO}}$. Prove 4.3, including the Ward identity in the longitudinal sector, and check the $(v_{\mathrm{F}}/c)^2$ and light-cone estimates of §2.
2. **Benchmark (2–3 weeks).** Reproduce Takada (1980) and Gastiasoro et al. (2019) from EPP's single propagator. Then add what they lack: first-principles $Z^*$, $\varepsilon^{\mathrm{e}}(\mathbf{q},\mathrm{i}\omega)$ from your RPA code, and anisotropy. Solve the linearized gap equation with the full Lindhard $W(\mathbf{q},\mathrm{i}\omega)$ of 4.1, and find where the split $\lambda_{\mathrm{ep}}/\mu^*$ treatment fails.
   - Use an intermediate-representation (sparse) Matsubara basis. $T_{\mathrm{c}}\approx0.1$–$0.4$ K against $\hbar\omega^{\mathrm{B}}_p$ up to 0.1–0.5 eV means about $5\times10^2$–$10^4$ positive Matsubara frequencies [derived here], and the Coulomb part of $W$ does not decay in $\omega^{\mathrm{B}}_p$.
   - Treat a generic polar semiconductor with $E_{\mathrm{F}}/\hbar\omega_{\mathrm{LO}}$ from 0.01 to 1. A one-LO-mode model is not SrTiO$_3$-like. For SrTiO$_3$, include the soft TO$_1$ (about 1–2 meV at low $T$, $\varepsilon_{\mathrm{s}}\approx2\times10^4$) and all three LO modes (about 21, 59 and 98 meV) [lit; e.g. Gastiasoro 2020].
3. **Interface (1–2 months).** Start from Zhou and Millis (2017), who found near-complete RPA screening, and test it with first-principles geometry and local fields. For an FeSe/SrTiO$_3$-like sheet, build $W(\mathbf{Q}_\parallel,\mathrm{i}\omega;0,0)$ from §B.10 plus a layered substrate Green's function (4.2). Compare $\lambda_{\mathrm{ep}}$ and the replica weights with Faeth's $\lambda_{\mathrm{ep}}=0.19\pm0.02$, and vary the substrate.
4. **Soft modes (consistency check, about 1 week).** In partition (A), recover Cochran softening for SrTiO$_3$/BaTiO$_3$ by adding the $\mathbf{G}\neq0$ dipole sum (E.454a; (B9.15)) to $\Pi^{\mathrm{ph}}$. With DFPT short-range force constants this reproduces DFPT by construction. Optional: show with partition (A) that the polar channel is negligible for metallic NbSe$_2$ (cf. Bianco et al. 2020). As a metal it has only dynamical Born charges (Marchese et al. 2024), and its CDW is short-range and adiabatic.
5. **Polaritons (months, exploratory).** Put MoSe$_2$ in a planar cavity, with $\Pi^{\mathrm{e,BSE}}$ and $D^0_{\mathrm{cav}}$. Verify that $\Pi(\omega)-\Pi(0)$ gives no $q=0$ instability while the paramagnetic part alone does.

## 8. Open questions

- With $\hbar\omega_{\mathrm{LO}}\gg E_{\mathrm{F}}$, which vertex corrections does the Ward identity fix in EPP's current-vertex language?
- Can a dynamical-$W$ replica calculation decide whether the FeSe/SrTiO$_3$ coupling is intrinsic, and how much does the sheet's own plasmon screen it? Zhou and Millis answered the screening part within model RPA.
- Can a substrate shift $T_{\mathrm{CDW}}$ through $\Pi^{\mathrm{ph}}$, given that $\mathbf{q}_{\mathrm{CDW}}$ is large?
- **[speculation]** Can a cavity change $T_{\mathrm{c}}$ in the thermodynamic limit other than through near-field environment modes, which amount to a substrate effect?
- **[speculation]** Does any crystal satisfy $\chi_{\mathrm{orb}}(\mathbf{q})>1/4\pi$?
- Can one EPP+BSE theory treat TiSe$_2$'s combined excitonic and lattice instability?

## Review responses

All 17 points were accepted, four of them with refinements. EPP pp. 18–21 (E.100–109, $k^{\mathcal{S}}$), Giustino p. 015003-21 and the study notes E6 (item 14), E7 (E.309–311) and E9 (E.428, E.452–454), together with §B.3 and the corrections section of primer v1, were re-read. Every citation was re-searched.

1. **Lattice double counting (major): accepted, refined.** §4.4 now gives partitions (A) and (B), §5 has the new rule, and §7 step 4 uses (A). EPP pp. 18–21 confirm $g^0$ (E.107), the Hamiltonian E.109 with its double-counting caveat, and DFT short-range force constants. Refinement: (B) is labelled as the structure of G.119 (bare $Z_\kappa$, frozen density, full $W^{\mathrm{e}}$), because G.141 is only its nonadiabatic part. Giustino's remark that bare force constants are not useful (p. 015003-21) is now cited.
2. **$\Pi(\mathbf{q},0)=0$ at all $\mathbf{q}$ (major): accepted, refined.** The study note E6 already flagged that this deletes $\chi_{\mathrm{orb}}$. Refinement: for the point-dipole ions of E.311 the cancellation is exact at all $\mathbf{q}$ (E.309–310). Only the electronic E.281/E.383 lose orbital magnetism. §4.5, §5, the table and the verdict were updated.
3. **Longitudinal Nambu bubble and E.428′ (major): accepted.** Nambu 1960 was verified (it gives the generalized Ward identities and a gauge-invariant Meissner kernel).
4. **Dolgov misattribution (major): accepted, refined.** DKM 1981 was verified ("On an admissible sign of the static dielectric function of matter"). For Kirzhnits–Maksimov–Khomskii only the 1970 Lebedev preprint (in a reference list) and a 1972 translation listing were found, so the journal reference is marked unverified. Takada's use of their theory was confirmed in a citing paper. Refinement: the §1 sign argument holds in RPA even with crystal local fields, because the static $-\chi^0_{\mathbf{G}\mathbf{G}'}$ is positive semidefinite. The qualifier is therefore "RPA (no xc)", not "RPA without local fields".
5. **Missing prior work (major): accepted.** All were verified except two. Zhou and Millis is now confirmed (PRB 96, 054516, with DOI). Pimenov and Chubukov's npj Quantum Mater. 7, 45 is titled "Quantum phase transition in a clean superconductor with repulsive dynamical interaction". Their related "Twists and turns…" is an Ann. Phys. review, not this paper. Li and Sawatzky remain unverified. A "Dilute polar pairing" bullet was added to §3, and steps 2–3 and §8 Q2 were rewritten.
6. **Unrealistic plan (major): accepted.** The Matsubara counts were recomputed ($5\times10^2$–$10^4$). The SrTiO$_3$ data were confirmed: TO$_1$ about 1–2 meV at low $T$, $\varepsilon_{\mathrm{s}}\approx2\times10^4$, LO about 21, 59 and 98 meV. Marchese et al. 2024 was verified. The one disagreement is the effort figure: with step 2 now 2–3 weeks, steps 1–2 take 3–5 weeks, not 2–4.
7. **Origin of the Lorentz field (minor): accepted.** This is consistent with the corrections section of primer v1 ("From E.380 on only $\mathbf{G}=\mathbf{G}'=0$ is kept"; in v2, the note on (B9.15) in §B.9) and §B.3 ($Z^*$ already contains the local field).
8. **"Only hardens" (minor): accepted.** Heading renamed, and the $\mathrm{L}_-$ sentence added.
9. **Light-cone numbers (minor): accepted.** Recomputed: $2.5\times10^{-11}$–$2.5\times10^{-7}$ (3D) and $3\times10^{-6}$–$3\times10^{-4}$ (2D). For SrTiO$_3$, $k_{\mathrm{F}}=0.025$ Å$^{-1}$; $v_{\mathrm{F}}=1.3\times10^4$ m/s from $E_{\mathrm{F}}=1.1$ meV and up to $2\times10^4$ m/s for $m^*=1.5\,m_{\mathrm{e}}$; $(v_{\mathrm{F}}/c)^2\approx2$–$4\times10^{-9}$.
10. **$\bar\omega_{\mathrm{LO}}$ (minor): accepted.** Rederived: $W-W^{\mathrm{e}}=-\frac{4\pi e^2}{q^2}\frac{(\omega^{\mathrm{ion}}_{\mathrm{p}})^2/(\varepsilon^{\mathrm{e}})^2}{\omega^2+\omega^2_{\mathrm{TO}}+(\omega^{\mathrm{ion}}_{\mathrm{p}})^2/\varepsilon^{\mathrm{e}}}$.
11. **$g^0$ is in EPP (minor): accepted.** Checked on pp. 18–21.
12. **Full $W_{\mathbf{G}\mathbf{G}'}$ (minor): accepted.**
13. **Normal self-energy (minor): accepted.**
14. **$K^{\mathrm{T}}$, Ampère, Schlawin (minor): accepted.** The published abstract (pair-density wave, "low-kelvin") and the v1 abstract (about 0.5 K) were confirmed. §5 "Regularization" was aligned.
15. **Substrate half-space (minor): accepted.**
16. **$g\to g^{\mathrm{lr}}$ (minor): accepted.**
17. **Citation status (minor): accepted.** Faeth's DOI was confirmed, and Faeth moved to the verified list. The Thomas record (JCP 162(13), 2025, doi:10.1063/5.0231202) was confirmed. Lee's $\lambda_{\mathrm{ep}}\approx0.5$ was not found. Andolina's criterion was confirmed. Kasprzak's DOI turned up in a RePEc record and is now listed. Kogar's DOI is still flagged, and Lian, Ali and Wong (arXiv:1901.11223; authors confirmed, journal not) is cited for the debate. Johannes–Mazin appears in a citing paper but remains unverified. For Gastiasoro (Ann. Phys.) the article number and year were confirmed, the volume was not.
