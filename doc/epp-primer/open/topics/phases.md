# Ordered phases: superconductivity, charge density waves and condensates

Labels: **[G]** Giustino RMP 2017; **[EPP]** the manuscript (draft v2, corrected as in [§Z.3](../../sections/Z03-corrections.html)); **[lit]** published work in §6; **[derived here]** derived for this note, unchecked; **[proposal]**; **[speculation]**. Gaussian units.

EPP as written cannot produce any ordered phase: it is a normal-state theory of three bosonic propagators, and its static lattice self-energy only stabilises. It is, however, a natural language for pairing in polar and dilute superconductors, and it already enforces the gauge constraint against uniform photon condensation.

| Phase | Instability of | EPP has | EPP lacks |
|---|---|---|---|
| Superconductivity | Cooper channel of $G$, interaction $D$ | polar part of $D=D^{\rm e}+D^{\rm ph}$, dynamically screened | $G$ with anomalous (Nambu–Gor'kov) parts; short-range $g$ |
| CDW, ferroelectric | $F^{-1}=0$: $\tilde\omega^2_\sigma(\mathbf Q)\to0$ | $F=F^0+F^0\Omega^{\rm e}F$ (E.398) | short-range electron–lattice coupling, the $\bG\neq0$ dipole sum, anharmonicity |
| Photon or polariton condensate | $D^{-1}(\bq,0)=0$, or driven occupation of a dispersion minimum | $\Pi(\omega)-\Pi(0)$ (E.281, E.311), which enforces the uniform no-go theorem but also sets $\chi_{\rm orb}\equiv0$ | static transverse (orbital) response, excitons (BSE), cavity $D^0$, kinetics |

**Verdict.**
- **Superconductivity.** With a Nambu–Gor'kov extension (4.1), EPP could predict $T_c$ for polar pairing with dynamical screening and LO–plasmon modes. That problem has been solved in model form since Takada (1980). EPP would add first-principles $Z^*$, $\epsilon^{\rm e}(\bq,i\nu)$ and anisotropy, all from one propagator. Best tests: doped SrTiO$_3$ and FeSe/SrTiO$_3$. Retarded transverse photons are irrelevant in bulk [derived here], and cavity enhancement is speculative.
- **CDW.** With the $\bG\neq0$ dipole sum (E.454a), EPP gives Cochran (ferroelectric) softening. With DFPT short-range force constants as input, however, this only reproduces DFPT. $2k_F$-type CDWs (NbSe$_2$, TaS$_2$) need short-range coupling and anharmonicity, which EPP does not improve, and TiSe$_2$ needs excitons as well.
- **Condensates.** Phonons and photons have no equilibrium BEC. Uniform photon condensation is forbidden by gauge invariance, and a finite-$\bq$ condensate needs a static orbital response that EPP's $\Pi(\omega)-\Pi(0)$ sets to zero (4.5). Polariton condensates need BSE, a cavity and kinetics; EPP-type theory would supply the dispersion and the relaxation vertices.

## 1. What is missing in EPP

- **No fermion propagator.** EPP solves three Dyson equations:
  - $D=D^0+D^0\Pi_{jj}D$, with $D^0=D^{0l}+D^{0t}$ (E.375, E.379–380);
  - $D^{\rm e}$ (E.388);
  - $F$, with $\Omega^{\rm e}=-J^*D^{\rm e}J$ (E.398–399).

  Electrons appear only as the RPA bubble (E.282–283; [§B.5](../../sections/B05-photon-se.html)–[§B.7](../../sections/B07-phonon-se.html)). There is no $G$, no $\Sigma$, no pair amplitude $\langle\psi_\uparrow\psi_\downarrow\rangle$, hence no gap equation.
- **Only the polar channel.** Ions couple only through $\frac1c\dot{\mathbf p}\cdot\bA$ with $\mathbf p=Z^*Q$ (E.134–135), on a frozen crystal potential (E.136) and a harmonic lattice. Deformation-potential coupling, which dominates $\lambda$ in conventional superconductors and in CDW metals, is absent ([§B.8](../../sections/B08-frohlich.html), [§Z.1](../../sections/Z01-outlook.html)). The bare vertex $g^0$ appears in EPP Sec. 2 (E.100–109, with the double-counting caveat of pp. 20–21) but is dropped from E.134. EPP does keep harmonic short-range force constants ($k^{\mathcal S}$, p. 21). What it lacks is short-range electron–lattice coupling and anharmonicity.
- **The static lattice self-energy only hardens; no zero-frequency root [derived here].** At $c\to\infty$, $\Omega^{\rm e}_{\sigma\sigma}=\hbar W^{\rm e}|\bq\cdot\bj_\sigma|^2/2\Omega_{\rm uc}\omega_\sigma$ (E.452). In RPA, $\chi^0(\bq,0)\le0$, hence $\epsilon^{\rm e}(\bq,0)\ge1$ and $W^{\rm e}(\bq,0)>0$, so E.454 can never drive $\omega^2$ to zero. (The sign also holds with crystal local fields, because the static $-\chi^0_{\bG\bG'}$ is positive semidefinite.) This sign argument is specific to RPA. Beyond it, with exchange–correlation, $\epsilon(\bq\neq0,0)<0$ is admissible [Dolgov, Kirzhnits, Maksimov, Rev. Mod. Phys. **53**, 81 (1981)]. Dynamically, E.453 gives an $L_-$ branch below $\omega_\sigma$ in doped crystals, but a zero-frequency root would need $W^{\rm e}(\bq,0)<0$, which RPA forbids. Cochran softening, $\tilde\omega^2_{\rm TO}=\omega^2_O-\omega_p^2/3$ (E.160), comes from the Lorentz field (E.151) of the classical model. The Dyson route lacks it because E.380–399 keep only $\bG=0$ in the ionic vertex and in $D^0$ (macroscopic field only). The $\bG\neq0$ dipole sum (E.454a, Gonze–Lee) restores the Lorentz/Ewald term. Electronic local fields ("CLFs are mostly neglected", EPP p. 95) are a separate omission: they change $\epsilon^{\rm e}$, while $Z^*$, an input, already contains them ([§B.3](../../sections/B03-dielectric.html)).
- **Temperature and ladders.** $T$ enters only through $f_{n\bk}$ in $\Pi^{\rm e}$ ($n_B=0$ for bosons, E.246–250). There is no anharmonicity, and no particle–hole (exciton) or particle–particle (Cooper) ladder.

## 2. Why it matters

- **Light cone.** Retardation acts only for $q\lesssim q_c\sim10^{-5}$–$10^{-4}$ Å$^{-1}$ ([§X.1](../../sections/X01-light-cone.html)), while pairing samples $2k_F\sim0.2$–$2$ Å$^{-1}$. At larger $q$ the transverse (magnetic, Breit/Darwin) interaction persists at relative order $(v/c)^2$, about $10^{-5}$ of the Coulomb term for $v_F=10^6$ m/s. The Fermi-surface fraction with $|\bk-\bk'|<q_c$ is $(q_c/2k_F)^2\approx2.5\times10^{-11}$–$2.5\times10^{-7}$ in 3D, and $q_c/(\pi k_F)\approx3\times10^{-6}$–$3\times10^{-4}$ for a 2D Fermi line [derived here]. For SrTiO$_3$ at $5.5\times10^{17}$ cm$^{-3}$ (one spin-degenerate band), $2k_F\approx0.05$ Å$^{-1}$, $v_F\approx1.3$–$2\times10^4$ m/s and $(v_F/c)^2\approx2$–$4\times10^{-9}$ [derived here].
- **Dilute polar superconductors.** SrTiO$_3$ superconducts down to $5.5\times10^{17}$ cm$^{-3}$, where $E_F\approx1.1$ meV (Lin 2013), far below its polar phonons. Migdal's condition $\omega_{\rm ph}/E_F\ll1$ fails there. So does the separation of scales behind $\mu^*$ (G.204), because the plasma frequency falls among the phonons. The mechanism is still debated (Gastiasoro 2020): LO phonons (Gor'kov 2016), plasmons (Ruhman and Lee 2016) and soft TO modes have all been proposed.
- **Interfacial Fröhlich coupling.** Monolayer FeSe on SrTiO$_3$ opens gaps near 70 K (bulk FeSe: about 8 K). ARPES shows replica bands about 100 meV below the main bands, attributed to small-$q$ coupling to substrate oxygen phonons (Lee 2014). Estimates of the coupling are $\lambda\approx0.5$ (attributed to Lee 2014; not confirmed in the primary text) and $0.19\pm0.02$ (Faeth 2021). An extrinsic origin of the replicas, photoelectron energy loss to Fuchs–Kliewer phonons, has been argued (Li and Sawatzky, preprint), and Faeth et al. argue against it. In model RPA the FeSe electrons screen the dipolar coupling almost completely (Zhou and Millis 2017).
- **Cavities.** Models give a pair-density-wave state with $T_c$ in the low-kelvin range for a GaAs-like 2D gas (Schlawin 2019; about 0.5 K in the preprint). They give no enhancement of forward-scattering pairing, because mode softening offsets the coupling gain (Sentef 2018). A reported rise of Rb$_3$C$_{60}$ from 30 to 45 K (Thomas et al. 2025, figures as given in the preprint) is debated.
- **CDW.** In NbSe$_2$, $T_{\rm CDW}$ from SSCHA is about 59 K in bulk and 73 K in the monolayer, with ionic fluctuations dominating (Bianco 2020). In TiSe$_2$ an electronic mode softens to zero at finite $q$ near 190 K (Kogar 2017). The excitonic interpretation is debated (Lian, Ali and Wong, preprint).
- **Condensates.** Polariton condensation has been seen at 19 K in CdTe (Kasprzak 2006) and in monolayer MoSe$_2$ (Anton-Solanas 2021). Signals of interlayer-exciton condensation persist above 100 K in MoSe$_2$/WSe$_2$ (Wang 2019).

**For an RPA/GW/BSE practitioner.**
1. The Coulomb part of pairing, $V_{n\bk,n'\bk'}$ in G.209, is a Fermi-surface matrix element of your $W$. Use the full $W_{\bG\bG'}(\bq,i\nu)$ of your code, not EPP's head-only $W^{\rm e}$: the matrix elements involve $q$ up to $2k_F$ and umklapp terms. Add $W^{\rm ph}$ to the head ($\bG=\bG'=0$) only.
2. In dilute or 2D polar systems $W^{\rm e}$ and $W^{\rm ph}$ cannot be separated by retardation, and EPP builds both from one propagator.
3. The substrate sets both the exciton screening length $r_0$ ([§B.10](../../sections/B10-2d.html)) and the interfacial pairing.
4. Excitonic insulators and interlayer exciton condensates are BSE problems with an anomalous electron–hole self-energy.

## 3. What standard theory does

- **Pairing [G §XI].** McMillan–Allen–Dynes (G.204–207) and anisotropic Migdal–Eliashberg theory in Nambu–Gor'kov form (G.208–210), Wannier-interpolated (Margine 2013). Here $\lambda(i\omega)$ is $W_{\rm ph}=W_{\rm e}\mathcal DW_{\rm e}$ projected on pairs ([§A.9](../../sections/A09-applications.html)), and $\mu^*$ is either empirical or taken from the RPA $W$. SCDFT (G.211; Lüders 2005) builds $\mathcal K^{\rm ee}$ from the screened Coulomb interaction; Akashi and Arita (2013) made that kernel dynamical, so plasmons pair too.
- **Dilute polar pairing [lit].**
  - Takada (1980) solved the gap equation for n-SrTiO$_3$ with a dynamically screened polar coupling and the plasmon–LO modes. He used the dielectric formulation of pairing of Kirzhnits, Maksimov and Khomskii.
  - Gor'kov (2016) argued for nonadiabatic LO-phonon pairing, and Ruhman and Lee (2016) for plasmon pairing.
  - Gastiasoro, Chubukov and Fernandes (2019) solved phonon pairing at $E_F\lesssim\omega_{\rm ph}$ with Coulomb repulsion, both dressed by the electron polarizability.
  - An interaction that is repulsive at every Matsubara frequency can still pair, with a gap $\Delta(i\omega)$ that changes sign (Pimenov and Chubukov 2022).
- **Lattice instabilities [G §V.A, §VII].**
  - The DFPT force constants $\Pi^{\rm A}$ (G.128) already contain Kohn anomalies and the imaginary frequencies of CDWs, and $\Pi^{\rm NA}$ (G.141, G.145) adds the dynamics ([§A.6](../../sections/A06-phonon-se.html)).
  - Nesting predicts CDWs poorly: the $\bq$ dependence of the coupling decides (Johannes and Mazin 2008).
  - The distorted phase and $T_{\rm CDW}$ need anharmonicity (SSCHA; Bianco 2020).
  - Ferroelectric soft modes balance dipolar against short-range forces (Cochran 1960). They are computed with DFPT plus the Gonze–Lee dipole sum (cf. E.454a).
- **Condensates.**
  - The $A^2$ term forbids the superradiant transition (Rzążewski 1975; Nataf and Ciuti 2010).
  - Gauge invariance forbids uniform photon condensation (Andolina 2019). A finite-$q$ condensate requires $\chi_{\rm orb}(\bq)>1/4\pi$ in 3D (Andolina 2020).
  - Photon BEC needed a thermalizing dye (Klaers 2010), and polariton condensates are driven and dissipative.

## 4. How it could enter the EPP formalism

**4.1 Pairing kernel from $D$ [proposal].** With Nambu spinors $\Psi_\bk=(c_{\bk\uparrow},c^\dagger_{-\bk\downarrow})^T$, EPP's vertices become:
- the paramagnetic current $-\frac1c\bj\cdot\bA\,\tau_0$;
- the diamagnetic term $\frac{e^2}{2mc^2}\rho A^2\,\tau_3$;
- the charge $e\rho\,\tau_3$, via continuity (E.428′). E.428′ holds for the normal-state $H_0^{\rm e}$ (p. 94), so use it only in the linearized gap equation ($T\ge T_c$). Below $T_c$, work in Coulomb gauge from the start, with the scalar potential coupling through $\tau_3$.

In Coulomb-gauge form, which is exact for $c\to\infty$ (E.441–445), the self-energy is (schematic; prefactors to be fixed against G.208–210)

$$\hat\Sigma(\bk,i\omega_n)=-\frac{k_BT}{N}\sum_{\bk'm}\Big[\tau_3\hat G(\bk',i\omega_m)\tau_3\,|e\rho_{\bk\bk'}|^2W(\bq,i\nu)+\tau_0\hat G(\bk',i\omega_m)\tau_0\,K^{\rm T}_{\bk\bk'}(i\nu)\Big],$$

$$W=\frac{4\pi}{\bq\cdot\hat\epsilon\cdot\bq},\qquad\hat\epsilon=\hat\epsilon^{\rm e}+\hat\epsilon^{\rm ion}-1,\qquad W=W^{\rm e}+W^{\rm ph}$$

(head-only form; the Coulomb part needs $W_{\bG\bG'}$).

Here $\bq=\bk'-\bk$ and $i\nu=i\omega_n-i\omega_m$. The additivity of $\hat\epsilon$ is EPP p. 33 (cf. E.432, E.445). The transverse kernel $K^{\rm T}=\frac1c\bj^*\cdot D^t\cdot\bj$ is of relative order $(v/c)^2$. Subtract $D^{0t}$ only in the normal part, where it is a mass renormalization; in the anomalous part the vacuum transverse exchange is the physical Breit/Darwin interaction. For time-reversed pairs that interaction is repulsive (Ampère: antiparallel currents repel).

For one isotropic mode, $W^{\rm ph}(\bq,i\nu)=-\frac{4\pi}{q^2}\frac{\omega_p^2/\epsilon^2}{\nu^2+\tilde\omega_{\rm LO}^2}$, which is attractive [derived here; cf. self-energy topic §4.4]. With $\epsilon\to\epsilon^{\rm e}(\bq,i\nu)$ the formula holds pointwise in $i\nu$, with $\tilde\omega^2_{\rm LO}(\bq,i\nu)=\omega^2_{\rm TO}+\omega_p^2/\epsilon^{\rm e}(\bq,i\nu)$. On the real axis the poles of the full $W$ are the $L_\pm$ roots of E.453. For pairing, use $W(\bq,i\nu)=4\pi/\{q^2[\epsilon^{\rm e}(\bq,i\nu)+\omega_p^2/(\nu^2+\omega^2_{\rm TO})]\}$ without splitting. Linearized, G.209 holds with

$$\lambda_{n\bk,n'\bk'}(i\nu)-N_FV_{n\bk,n'\bk'}\;\to\;-N_F\langle n\bk,n'\bk'|W^{\rm e}+W^{\rm ph}|n\bk,n'\bk'\rangle(i\nu).$$

This is the dielectric formulation of pairing (Kirzhnits, Maksimov and Khomskii; journal reference unverified). Takada (1980) applied it to polar n-SrTiO$_3$ with the plasmon and the LO phonon.
- **New relative to G.209 (but known in model form since Takada 1980):** a dynamically screened Fröhlich vertex (E.448a), $\omega_{\rm LO}$ as an output, and coupled LO–plasmon modes. EPP's own contribution would be first-principles $Z^*$, $\epsilon^{\rm e}(\bq,i\nu)$ and anisotropy.
- **Still missing:** the short-range $g$, and vertex corrections for $\omega_{\rm LO}/E_F\gtrsim1$.

**4.2 2D and cavities [proposal].**
- **Substrates.** For a sheet on a polar substrate, replace the single-sheet $W^{\rm e}(\mathbf Q,\omega,z,z')$ (E.461) by a layered Green's function with $\epsilon_{\rm sub}(\omega)=\epsilon_\infty+\sum$ (ionic poles from the substrate's $\Pi^{\rm ion}$). This is a new ingredient: EPP's E.455–472 cover only a freestanding sheet with $\delta(z)$ ions, while the substrate is a half-space. The FeSe/SrTiO$_3$ coupling is then a Fuchs–Kliewer surface-phonon pole of $W(\mathbf Q,i\nu;0,0)$, small-$Q$ by construction.
- **Cavities.** Replace $D^0$ by the cavity Green's dyadic $D^0_{\rm cav}$; E.379 and E.398 keep their form.
  - Sentef's mechanism is $\Omega^{\rm e}$ and $g^{\rm dyn}$ (E.399, E.416) evaluated with $D^{\rm e}_{\rm cav}$.
  - Schlawin's is $K^{\rm T}[D_{\rm cav}]$ at finite pair momentum, which gives a pair-density wave.
  - Curtis's driven enhancement needs nonequilibrium kinetics, outside EPP.

**4.3 Meissner kernel [derived here; standard].** With Nambu propagators in $\Pi^{\rm e}$, the identity $\Pi^{\rm dia}=-\Pi^{\rm para}(\omega=0)$ behind E.281 fails: $\Pi^t_{jj}(\bq\to0,0)\propto n_s/m\neq0$ (the London kernel). Superconductivity then appears in EPP's photon sector as a photon mass, and E.279 must replace E.281. This applies to the transverse sector. In the longitudinal sector the bare Nambu bubble violates $\bq\cdot\Pi(\bq,0)\cdot\bq=0$, which gives a spurious longitudinal photon mass. The vertex must satisfy Nambu's Ward identity [Y. Nambu, Phys. Rev. **117**, 648 (1960), doi:10.1103/PhysRev.117.648]. That brings in the phase mode, which $D^{0l}$ pushes up to $\omega_{\rm pl}$.

**4.4 Soft modes [proposal].** Reinstate EPP's own $g^0$ (E.107) in E.134, minus its long-range $Z_{\rm ion}/|\bq+\bG|^2$ part, which duplicates the dipole channel of E.134. This gives a bare short-range vertex $g^{\rm S,b}$, with screened partner $g^{\rm S}$, and the lattice Dyson equation becomes

$$F=F^0+F^0\big(\Omega^{\rm e}+\Omega^{\rm S}+\Omega^{\rm LS}+\Omega^{\rm SL}\big)F,\qquad\Omega^{\rm S}_{\sigma\sigma'}=\sum_{\bk nm}g^{\rm S,b}_{nm\sigma}\,\chi^0_{n\bk,m\bk+\bq}(\omega)\,\big[g^{\rm S}_{nm\sigma'}(\omega)\big]^*.$$

As in E.412 and G.141, one vertex is bare and one screened; the cross terms $\Omega^{\rm LS}$ and $\Omega^{\rm SL}$ mix the polar and nonpolar channels. Only two partitions are consistent, and mixing them counts the valence electrons twice or three times:
- **(A) Screened/adiabatic (EPP's own choice, p. 21).** $F^0$ comes from DFPT or SSCHA short-range force constants, with the Gonze–Lee dipole–dipole part removed. $J$ is built with $Z^*$ and $\epsilon_\infty$ (E.454a). The short-range channel enters only through its nonadiabatic remainder $\Omega^{\rm S}(\omega)-\Omega^{\rm S}(0)$ (the G.145 structure), with $g^{\rm S}=g^{\rm DFPT}-g^{\rm L}$. The CDW instability is then an input through $F^0$. EPP adds the dynamics and the polar or substrate channel.
- **(B) Bare (the structure of G.119, whose nonadiabatic part is G.141).** $F^0$ comes from the bare ion–ion plus frozen-density terms. $J$ is built with the bare valence charges $Z_{\rm ion}$, not $Z^*$, and $g^{\rm S,b}$ is E.107 with its $Z_{\rm ion}/|\bq+\bG|^2$ tail removed. $Z^*$, $\epsilon_\infty$ and the Kohn anomalies must then come out of $\chi^0$ with local fields. SSCHA cannot be combined with (B), because it works on the Born–Oppenheimer surface, which already contains the electronic screening. Giustino notes that bare force constants give dispersions very different from the real ones (p. 015003-21), so (B) is a long project.

The instability criterion is

$$\tilde\omega^2_\sigma(\mathbf Q;T)=\omega^2_\sigma(\mathbf Q;T)+\frac{2\omega_\sigma}{\hbar}\,{\rm Re}\,\Omega_{\sigma\sigma}(\mathbf Q,\tilde\omega;T)=0.$$

Here $T$ enters through $f_{n\bk}$ in $\chi^0$ and, in partition (A) only, through an anharmonic $\omega^2_\sigma(T)$ from SSCHA. Keeping $\bG\neq0$ in $D^{\rm e}$ (E.454a, with the on-site Ewald term removed) restores the Lorentz field, and with it Cochran-type softening.

**4.5 Condensates [assessment].**
- **Phonons.** A "phonon condensate" $\langle Q_\sigma(\mathbf Q)\rangle\neq0$ is just the frozen distortion of 4.4.
- **Photons.** Photon condensation needs a zero transverse eigenvalue of $D^{-1}(\bq,0)=[D^0]^{-1}-\Pi_{jj}(\bq,0)$. E.281 and E.311 set $\Pi_{jj}(\bq,0)=0$ at every $\bq$. This enforces the uniform no-go theorem, but it also deletes the static orbital response at finite $\bq$. For the point-dipole ions of E.311 that is exact (E.309–310); for electrons it is not. A finite-$\bq$ condensate requires $\chi_{\rm orb}(\bq)>1/4\pi$ in 3D, Gaussian units (Andolina et al. 2020). To test this in EPP, keep E.279 (explicit dia + para, complete band sums) for the transverse part, and use E.281 only for the longitudinal part.
- **Polaritons.** Put $\Pi^{\rm BSE}_{jj}$ (long-range exchange removed, [§Z.1](../../sections/Z01-outlook.html)) and $D^0_{\rm cav}$ into E.379. The poles of $D$ then give the polariton dispersion, the Hopfield weights and the radiative widths. Relaxation (the bottleneck) needs exciton–phonon vertices, and condensation itself needs a kinetic or Gross–Pitaevskii layer.

## 5. Pitfalls

- **Double counting.**
  - Never combine $Z^*$ or $\epsilon_\infty$ with an $F^0$ that has no electron contribution, or DFPT/SSCHA $\omega^2$ with a full static $\Omega^{\rm S}$. Use partition (A) or (B) of 4.4.
  - Subtract $g^{\rm L}$ from the DFPT $g$.
  - Add no empirical $\mu^*$ on top of $W^{\rm e}$.
  - Coulomb is $D^{0l}$ (E.364); never add $v$ separately.
  - Normal self-energy: either start from bare (DFT) bands and keep the full $W^{\rm e}+W^{\rm ph}$ normal part, or start from $GW$ bands and keep only $W^{\rm ph}$ and the anomalous part of $W^{\rm e}$. Drop the $W^{\rm ph}$ normal part too if the input bands already contain the Fröhlich mass renormalization.
- **Gauge.** $\Pi(\omega)-\Pi(0)$ is exact for the longitudinal part at all $q$, and for the transverse part only as $q\to0$. At finite transverse $q$, or in the paired state, use dia + para with band sums complete enough to satisfy the f-sum rule. In the paired state the longitudinal vertex must also satisfy Nambu's Ward identity (4.3). Truncated band sums break the f-sum rule and fake a Drude weight or a $q=0$ photon instability.
- **Sum rules.** Impose the acoustic sum rules for $Z^*$ and for the force constants.
- **Migdal.** Include vertex corrections when $\omega_{\rm ph}/E_F\gtrsim1$.
- **Convergence.**
  - Fermi-surface integrals need up to $10^6$ $\bk$ points.
  - Small-$q$ forward scattering needs analytic integration.
  - Matsubara cutoffs renormalize $\mu$ to $\mu^*$. At $T_c\sim0.1$ K with $\nu$ up to several $\omega_{\rm LO}$, use an intermediate-representation (sparse) Matsubara basis.
  - $\mathbf Q_{\rm CDW}$ must lie on the grid, and smearing is not the physical $T$.
- **Regularization.** The $1/\omega^2$ of $D^{0l}$ and the noncommuting limits $\bq,\omega\to0$ need care ([§B.9](../../sections/B09-electrostatic.html)). Subtract $D^{0t}$ only in the normal part (4.1).

## 6. State of the art (verified)

1. F. Giustino, Rev. Mod. Phys. **89**, 015003 (2017), doi:10.1103/RevModPhys.89.015003. Eliashberg, SCDFT, nonadiabatic phonons.
2. E. R. Margine, F. Giustino, Phys. Rev. B **87**, 024505 (2013), doi:10.1103/PhysRevB.87.024505. Anisotropic Eliashberg with Wannier interpolation.
3. M. Lüders *et al.*, Phys. Rev. B **72**, 024545 (2005), doi:10.1103/PhysRevB.72.024545. SCDFT.
4. R. Akashi, R. Arita, Phys. Rev. Lett. **111**, 057006 (2013), doi:10.1103/PhysRevLett.111.057006. Dynamical Coulomb kernel; plasmon-assisted pairing.
5. X. Lin *et al.*, Phys. Rev. X **3**, 021002 (2013), doi:10.1103/PhysRevX.3.021002. SrTiO$_3$ superconducting at $E_F\approx1.1$ meV.
6. J. J. Lee *et al.*, Nature **515**, 245 (2014), doi:10.1038/nature13894. FeSe/SrTiO$_3$ replica bands.
7. M. A. Sentef, M. Ruggenthaler, A. Rubio, Sci. Adv. **4**, eaau6969 (2018), doi:10.1126/sciadv.aau6969. Cavity phonon polaritons and pairing.
8. F. Schlawin, A. Cavalleri, D. Jaksch, Phys. Rev. Lett. **122**, 133602 (2019), doi:10.1103/PhysRevLett.122.133602 (DOI by APS pattern). Cavity-photon pairing; a pair-density wave with $T_c$ in the low-kelvin range (about 0.5 K in arXiv:1804.07142v1).
9. J. B. Curtis *et al.*, Phys. Rev. Lett. **122**, 167002 (2019), doi:10.1103/PhysRevLett.122.167002. Driven cavity Eliashberg effect.
10. R. Bianco *et al.*, Phys. Rev. Lett. **125**, 106101 (2020), doi:10.1103/PhysRevLett.125.106101. Anharmonic CDW in NbSe$_2$.
11. A. Kogar *et al.*, Science **358**, 1314 (2017), doi:10.1126/science.aam6432 (DOI from secondary sources only). Exciton condensation in TiSe$_2$; interpretation debated.
12. G. M. Andolina *et al.*, Phys. Rev. B **100**, 121109(R) (2019), doi:10.1103/PhysRevB.100.121109. No-go theorem.
13. G. M. Andolina *et al.*, Phys. Rev. B **102**, 125137 (2020), doi:10.1103/PhysRevB.102.125137 (DOI from the publisher's file name). Finite-$q$ condensation; criterion $\chi_{\rm orb}(\bq)>1/4\pi$ in 3D.
14. J. Kasprzak *et al.*, Nature **443**, 409 (2006), doi:10.1038/nature05131 (DOI from the RePEc record). Polariton BEC.
15. C. Anton-Solanas *et al.*, Nat. Mater. **20**, 1233 (2021), doi:10.1038/s41563-021-01000-8. Polariton condensate in MoSe$_2$.
16. Y. Nambu, Phys. Rev. **117**, 648 (1960), doi:10.1103/PhysRev.117.648. Ward identity and a gauge-invariant Meissner kernel in the paired state.
17. Y. Takada, J. Phys. Soc. Jpn. **49**, 1267 (1980) (DOI not confirmed). Polar n-SrTiO$_3$ with plasmon and LO phonon in the dielectric formulation.
18. M. N. Gastiasoro, A. V. Chubukov, R. M. Fernandes, Phys. Rev. B **99**, 094524 (2019), doi:10.1103/PhysRevB.99.094524. Phonon pairing at $E_F\lesssim\omega_{\rm ph}$ with Coulomb repulsion.
19. J. Ruhman, P. A. Lee, Phys. Rev. B **94**, 224515 (2016), doi:10.1103/PhysRevB.94.224515. Plasmon pairing in dilute SrTiO$_3$.
20. L. P. Gor'kov, Proc. Natl. Acad. Sci. USA **113**, 4646 (2016), doi:10.1073/pnas.1604145113. Nonadiabatic LO-phonon pairing in dilute SrTiO$_3$.
21. D. Pimenov, A. V. Chubukov, npj Quantum Mater. **7**, 45 (2022), "Quantum phase transition in a clean superconductor with repulsive dynamical interaction" (DOI not confirmed). Pairing from a dynamically repulsive interaction.
22. Y. Zhou, A. J. Millis, Phys. Rev. B **96**, 054516 (2017), doi:10.1103/PhysRevB.96.054516; arXiv:1703.04021. RPA screening of the SrTiO$_3$ dipolar coupling by the FeSe electrons.
23. G. Marchese *et al.*, Nat. Phys. **20**, 88 (2024), doi:10.1038/s41567-023-02203-3. Dynamical Born charges in metals.
24. C. Faeth *et al.*, Phys. Rev. Lett. **127**, 016803 (2021), doi:10.1103/PhysRevLett.127.016803. Intrinsic replica bands, $\lambda=0.19\pm0.02$.

Also cited:
- O. V. Dolgov, D. A. Kirzhnits, E. G. Maksimov, Rev. Mod. Phys. **53**, 81 (1981), doi:10.1103/RevModPhys.53.81. Admissible sign of $\epsilon(\bq,0)$.
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

1. **Derive (1–2 weeks).** Write EPP's vertices in Nambu form. Show that 4.1 reduces to G.208–210 with $|g|^2\to|g^{\rm L}|^2$ (E.448a) and $W_{\rm ph}=W_{\rm e}\mathcal DW_{\rm e}$. Show that integrating out $|\nu|>\omega_c$ reproduces Morel–Anderson only for $\omega_{\rm pl}\gg\omega_{\rm LO}$. Prove 4.3, including the Ward identity in the longitudinal sector, and check the $(v/c)^2$ and light-cone estimates of §2.
2. **Benchmark (2–3 weeks).** Reproduce Takada (1980) and Gastiasoro et al. (2019) from EPP's single propagator. Then add what they lack: first-principles $Z^*$, $\epsilon^{\rm e}(\bq,i\nu)$ from your RPA code, and anisotropy. Solve the linearized gap equation with the full Lindhard $W(\bq,i\nu)$ of 4.1, and find where the split $\lambda/\mu^*$ treatment fails.
   - Use an intermediate-representation (sparse) Matsubara basis. $T_c\approx0.1$–$0.4$ K against $\nu$ up to 0.1–0.5 eV means about $5\times10^2$–$10^4$ positive Matsubara frequencies [derived here], and the Coulomb part of $W$ does not decay in $\nu$.
   - Treat a generic polar semiconductor with $E_F/\omega_{\rm LO}$ from 0.01 to 1. A one-LO-mode model is not SrTiO$_3$-like. For SrTiO$_3$, include the soft TO$_1$ (about 1–2 meV at low $T$, $\epsilon_0\approx2\times10^4$) and all three LO modes (about 21, 59 and 98 meV) [lit; e.g. Gastiasoro 2020].
3. **Interface (1–2 months).** Start from Zhou and Millis (2017), who found near-complete RPA screening, and test it with first-principles geometry and local fields. For an FeSe/SrTiO$_3$-like sheet, build $W(\mathbf Q,i\nu;0,0)$ from §B.10 plus a layered substrate Green's function (4.2). Compare $\lambda$ and the replica weights with Faeth's $\lambda=0.19\pm0.02$, and vary the substrate.
4. **Soft modes (consistency check, about 1 week).** In partition (A), recover Cochran softening for SrTiO$_3$/BaTiO$_3$ by adding the $\bG\neq0$ dipole sum (E.454a) to $\Omega^{\rm e}$. With DFPT short-range force constants this reproduces DFPT by construction. Optional: show with partition (A) that the polar channel is negligible for metallic NbSe$_2$ (cf. Bianco et al. 2020). As a metal it has only dynamical Born charges (Marchese et al. 2024), and its CDW is short-range and adiabatic.
5. **Polaritons (months, exploratory).** Put MoSe$_2$ in a planar cavity, with $\Pi^{\rm BSE}$ and $D^0_{\rm cav}$. Verify that $\Pi(\omega)-\Pi(0)$ gives no $q=0$ instability while the paramagnetic part alone does.

## 8. Open questions

- With $\omega_{\rm LO}\gg E_F$, which vertex corrections does the Ward identity fix in EPP's current-vertex language?
- Can a dynamical-$W$ replica calculation decide whether the FeSe/SrTiO$_3$ coupling is intrinsic, and how much does the sheet's own plasmon screen it? Zhou and Millis answered the screening part within model RPA.
- Can a substrate shift $T_{\rm CDW}$ through $\Omega^{\rm e}$, given that $\mathbf Q_{\rm CDW}$ is large?
- **[speculation]** Can a cavity change $T_c$ in the thermodynamic limit other than through near-field environment modes, which amount to a substrate effect?
- **[speculation]** Does any crystal satisfy $\chi_{\rm orb}(\bq)>1/4\pi$?
- Can one EPP+BSE theory treat TiSe$_2$'s combined excitonic and lattice instability?

## Review responses

All 17 points were accepted, four of them with refinements. EPP pp. 18–21 (E.100–109, $k^{\mathcal S}$), Giustino p. 015003-21 and the study notes E6 (item 14), E7 (E.309–311) and E9 (E.428, E.452–454), together with §B.3 and §Z.3, were re-read. Every citation was re-searched.

1. **Lattice double counting (major): accepted, refined.** §4.4 now gives partitions (A) and (B), §5 has the new rule, and §7 step 4 uses (A). EPP pp. 18–21 confirm $g^0$ (E.107), the Hamiltonian E.109 with its double-counting caveat, and DFT short-range force constants. Refinement: (B) is labelled as the structure of G.119 (bare $Z_\kappa$, frozen density, full $W_{\rm e}$), because G.141 is only its nonadiabatic part. Giustino's remark that bare force constants are not useful (p. 015003-21) is now cited.
2. **$\Pi(\bq,0)=0$ at all $\bq$ (major): accepted, refined.** The study note E6 already flagged that this deletes $\chi_{\rm orb}$. Refinement: for the point-dipole ions of E.311 the cancellation is exact at all $\bq$ (E.309–310). Only the electronic E.281/E.383 lose orbital magnetism. §4.5, §5, the table and the verdict were updated.
3. **Longitudinal Nambu bubble and E.428′ (major): accepted.** Nambu 1960 was verified (it gives the generalized Ward identities and a gauge-invariant Meissner kernel).
4. **Dolgov misattribution (major): accepted, refined.** DKM 1981 was verified ("On an admissible sign of the static dielectric function of matter"). For Kirzhnits–Maksimov–Khomskii only the 1970 Lebedev preprint (in a reference list) and a 1972 translation listing were found, so the journal reference is marked unverified. Takada's use of their theory was confirmed in a citing paper. Refinement: the §1 sign argument holds in RPA even with crystal local fields, because the static $-\chi^0_{\bG\bG'}$ is positive semidefinite. The qualifier is therefore "RPA (no xc)", not "RPA without local fields".
5. **Missing prior work (major): accepted.** All were verified except two. Zhou and Millis is now confirmed (PRB 96, 054516, with DOI). Pimenov and Chubukov's npj Quantum Mater. 7, 45 is titled "Quantum phase transition in a clean superconductor with repulsive dynamical interaction". Their related "Twists and turns…" is an Ann. Phys. review, not this paper. Li and Sawatzky remain unverified. A "Dilute polar pairing" bullet was added to §3, and steps 2–3 and §8 Q2 were rewritten.
6. **Unrealistic plan (major): accepted.** The Matsubara counts were recomputed ($5\times10^2$–$10^4$). The SrTiO$_3$ data were confirmed: TO$_1$ about 1–2 meV at low $T$, $\epsilon_0\approx2\times10^4$, LO about 21, 59 and 98 meV. Marchese et al. 2024 was verified. The one disagreement is the effort figure: with step 2 now 2–3 weeks, steps 1–2 take 3–5 weeks, not 2–4.
7. **Origin of the Lorentz field (minor): accepted.** This is consistent with §Z.3 ("From E.380 on only $\bG=\bG'=0$ is kept") and §B.3 ($Z^*$ already contains the local field).
8. **"Only hardens" (minor): accepted.** Heading renamed, and the $L_-$ sentence added.
9. **Light-cone numbers (minor): accepted.** Recomputed: $2.5\times10^{-11}$–$2.5\times10^{-7}$ (3D) and $3\times10^{-6}$–$3\times10^{-4}$ (2D). For SrTiO$_3$, $k_F=0.025$ Å$^{-1}$; $v_F=1.3\times10^4$ m/s from $E_F=1.1$ meV and up to $2\times10^4$ m/s for $m^*=1.5$; $(v_F/c)^2\approx2$–$4\times10^{-9}$.
10. **$\tilde\omega_{\rm LO}$ (minor): accepted.** Rederived: $W-W^{\rm e}=-\frac{4\pi}{q^2}\frac{\omega_p^2/\epsilon^{\rm e2}}{\nu^2+\omega^2_{\rm TO}+\omega_p^2/\epsilon^{\rm e}}$.
11. **$g^0$ is in EPP (minor): accepted.** Checked on pp. 18–21.
12. **Full $W_{\bG\bG'}$ (minor): accepted.**
13. **Normal self-energy (minor): accepted.**
14. **$K^{\rm T}$, Ampère, Schlawin (minor): accepted.** The published abstract (pair-density wave, "low-kelvin") and the v1 abstract (about 0.5 K) were confirmed. §5 "Regularization" was aligned.
15. **Substrate half-space (minor): accepted.**
16. **$g\to g^{\rm L}$ (minor): accepted.**
17. **Citation status (minor): accepted.** Faeth's DOI was confirmed, and Faeth moved to the verified list. The Thomas record (JCP 162(13), 2025, doi:10.1063/5.0231202) was confirmed. Lee's $\lambda\approx0.5$ was not found. Andolina's criterion was confirmed. Kasprzak's DOI turned up in a RePEc record and is now listed. Kogar's DOI is still flagged, and Lian, Ali and Wong (arXiv:1901.11223; authors confirmed, journal not) is cited for the debate. Johannes–Mazin appears in a citing paper but remains unverified. For Gastiasoro (Ann. Phys.) the article number and year were confirmed, the volume was not.
