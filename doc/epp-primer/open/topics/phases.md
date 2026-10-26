# Ordered phases: superconductivity, charge density waves and condensates

Labels: **[G]** Giustino RMP 2017; **[EPP]** the manuscript (draft v2, corrected as in [§Z.3](../../sections/Z03-corrections.html)); **[lit]** published work in §6; **[derived here]** derived for this note, unchecked; **[proposal]**; **[speculation]**. Gaussian units.

EPP as written cannot produce any ordered phase: it is a normal-state theory of three bosonic propagators, and its lattice self-energy only stabilises. It is, however, a natural language for pairing in polar and dilute superconductors, and it already enforces the gauge constraint against photon condensation.

| Phase | Instability of | EPP has | EPP lacks |
|---|---|---|---|
| Superconductivity | Cooper channel of $G$, interaction $D$ | polar part of $D=D^{\rm e}+D^{\rm ph}$, dynamically screened | $G$ with anomalous (Nambu–Gor'kov) parts; short-range $g$ |
| CDW, ferroelectric | $F^{-1}=0$: $\tilde\omega^2_\sigma(\mathbf Q)\to0$ | $F=F^0+F^0\Omega^{\rm e}F$ (E.398) | short-range channel, local fields, anharmonicity |
| Photon or polariton condensate | $D^{-1}(\bq,0)=0$, or driven occupation of a dispersion minimum | $\Pi(\omega)-\Pi(0)$ (E.281, E.311) | excitons (BSE), cavity $D^0$, kinetics |

**Verdict.**
- **Superconductivity.** With a Nambu–Gor'kov extension (4.1), EPP could predict $T_c$ for polar pairing with dynamical screening and LO–plasmon modes. Best tests: doped SrTiO$_3$ and FeSe/SrTiO$_3$. Retarded transverse photons are irrelevant in bulk [derived here], and cavity enhancement is speculative.
- **CDW.** With local fields, EPP gives Cochran (ferroelectric) softening. $2k_F$-type CDWs (NbSe$_2$, TaS$_2$) need short-range coupling and anharmonicity, which EPP does not improve, and TiSe$_2$ needs excitons as well.
- **Condensates.** Phonons and photons have no equilibrium BEC, and uniform photon condensation is forbidden by gauge invariance. Polariton condensates need BSE, a cavity and kinetics; EPP-type theory would supply the dispersion and the relaxation vertices.

## 1. What is missing in EPP

- **No fermion propagator.** EPP solves three Dyson equations:
  - $D=D^0+D^0\Pi_{jj}D$, with $D^0=D^{0l}+D^{0t}$ (E.375, E.379–380);
  - $D^{\rm e}$ (E.388);
  - $F$, with $\Omega^{\rm e}=-J^*D^{\rm e}J$ (E.398–399).

  Electrons appear only as the RPA bubble (E.282–283; [§B.5](../../sections/B05-photon-se.html)–[§B.7](../../sections/B07-phonon-se.html)). There is no $G$, no $\Sigma$, no pair amplitude $\langle\psi_\uparrow\psi_\downarrow\rangle$, hence no gap equation.
- **Only the polar channel.** Ions couple only through $\frac1c\dot{\mathbf p}\cdot\bA$ with $\mathbf p=Z^*Q$ (E.134–135), on a frozen crystal potential (E.136) and a harmonic lattice. Deformation-potential coupling, which dominates $\lambda$ in conventional superconductors and in CDW metals, is absent ([§B.8](../../sections/B08-frohlich.html), [§Z.1](../../sections/Z01-outlook.html)).
- **The lattice self-energy only hardens [derived here].** At $c\to\infty$, $\Omega^{\rm e}_{\sigma\sigma}=\hbar W^{\rm e}|\bq\cdot\bj_\sigma|^2/2\Omega_{\rm uc}\omega_\sigma$ (E.452). Without local fields, RPA gives $\chi^0(\bq,0)\le0$, hence $\epsilon^{\rm e}(\bq,0)\ge1$ and $W^{\rm e}(\bq,0)>0$, so E.454 can never drive $\omega^2$ to zero. Cochran softening, $\tilde\omega^2_{\rm TO}=\omega^2_O-\omega_p^2/3$ (E.160), comes from the Lorentz field (E.151) of the classical model. The Dyson route has no such field, because local fields are dropped after E.380 ("CLFs are mostly neglected", EPP p. 95).
- **Temperature and ladders.** $T$ enters only through $f_{n\bk}$ in $\Pi^{\rm e}$ ($n_B=0$ for bosons, E.246–250); there is no anharmonicity, and no particle–hole (exciton) or particle–particle (Cooper) ladder.

## 2. Why it matters

- **Light cone.** The transverse field acts only for $q\lesssim q_c\sim10^{-5}$–$10^{-4}$ Å$^{-1}$ ([§X.1](../../sections/X01-light-cone.html)), while pairing samples $2k_F\sim0.2$–$2$ Å$^{-1}$. Transverse exchange is of order $(v_F/c)^2\approx10^{-5}$ of the Coulomb term for $v_F=10^6$ m/s, and the Fermi-surface fraction with $|\bk-\bk'|<q_c$ is of order $(q_c/2k_F)^2\sim10^{-9}$–$10^{-7}$ [derived here].
- **Dilute polar superconductors.** SrTiO$_3$ superconducts down to $5.5\times10^{17}$ cm$^{-3}$, where $E_F\approx1.1$ meV (Lin 2013), far below its polar phonons. Migdal's condition $\omega_{\rm ph}/E_F\ll1$ fails there. So does the separation of scales behind $\mu^*$ (G.204), because the plasma frequency falls among the phonons. The mechanism (LO phonons, soft TO modes, plasmons) is still debated (Gastiasoro 2020).
- **Interfacial Fröhlich coupling.** Monolayer FeSe on SrTiO$_3$ opens gaps near 70 K (bulk FeSe: about 8 K). ARPES shows replica bands about 100 meV below the main bands, attributed to small-$q$ coupling to substrate oxygen phonons (Lee 2014). The coupling is estimated at $\lambda\approx0.5$ (Lee 2014) or $0.19\pm0.02$ (Faeth 2021), and an extrinsic origin of the replicas has also been argued.
- **Cavities.** Models give $T_c\sim0.5$ K for a GaAs-like 2D gas (Schlawin 2019), and no enhancement of forward-scattering pairing, since mode softening offsets the coupling gain (Sentef 2018). A reported rise of Rb$_3$C$_{60}$ from 30 to 45 K (Thomas, preprint) is debated.
- **CDW.** In NbSe$_2$, $T_{\rm CDW}$ from SSCHA is about 59 K in bulk and 73 K in the monolayer, with ionic fluctuations dominating (Bianco 2020). In TiSe$_2$ an electronic mode softens to zero at finite $q$ near 190 K (Kogar 2017).
- **Condensates.** Polariton condensation has been seen at 19 K in CdTe (Kasprzak 2006) and in monolayer MoSe$_2$ (Anton-Solanas 2021). Signals of interlayer-exciton condensation persist above 100 K in MoSe$_2$/WSe$_2$ (Wang 2019).

**For an RPA/GW/BSE practitioner.**
1. The Coulomb part of pairing, $V_{n\bk,n'\bk'}$ in G.209, is a Fermi-surface matrix element of your $W$.
2. In dilute or 2D polar systems $W^{\rm e}$ and $W^{\rm ph}$ cannot be separated by retardation, and EPP builds both from one propagator.
3. The substrate sets both the exciton screening length $r_0$ ([§B.10](../../sections/B10-2d.html)) and the interfacial pairing.
4. Excitonic insulators and interlayer exciton condensates are BSE problems with an anomalous electron–hole self-energy.

## 3. What standard theory does

- **Pairing [G §XI].** McMillan–Allen–Dynes (G.204–207) and anisotropic Migdal–Eliashberg theory in Nambu–Gor'kov form (G.208–210), Wannier-interpolated (Margine 2013). Here $\lambda(i\omega)$ is $W_{\rm ph}=W_{\rm e}\mathcal DW_{\rm e}$ projected on pairs ([§A.9](../../sections/A09-applications.html)), and $\mu^*$ is either empirical or taken from the RPA $W$. SCDFT (G.211; Lüders 2005) builds $\mathcal K^{\rm ee}$ from the screened Coulomb interaction; Akashi and Arita (2013) made that kernel dynamical, so plasmons pair too.
- **Lattice instabilities [G §V.A, §VII].**
  - The DFPT force constants $\Pi^{\rm A}$ (G.128) already contain Kohn anomalies and the imaginary frequencies of CDWs, and $\Pi^{\rm NA}$ (G.141, G.145) adds the dynamics ([§A.6](../../sections/A06-phonon-se.html)).
  - Nesting predicts CDWs poorly: the $\bq$ dependence of the coupling decides (Johannes and Mazin 2008).
  - The distorted phase and $T_{\rm CDW}$ need anharmonicity (SSCHA; Bianco 2020).
  - Ferroelectric soft modes balance dipolar against short-range forces (Cochran 1960). They are computed with DFPT plus the Gonze–Lee dipole sum (cf. E.454a).
- **Condensates.**
  - The $A^2$ term forbids the superradiant transition (Rzążewski 1975; Nataf and Ciuti 2010).
  - Gauge invariance forbids uniform photon condensation (Andolina 2019); a finite-$q$ condensate requires $\chi_{\rm orb}(\bq)>1/4\pi$ (Andolina 2020).
  - Photon BEC needed a thermalizing dye (Klaers 2010), and polariton condensates are driven and dissipative.

## 4. How it could enter the EPP formalism

**4.1 Pairing kernel from $D$ [proposal].** With Nambu spinors $\Psi_\bk=(c_{\bk\uparrow},c^\dagger_{-\bk\downarrow})^T$, EPP's vertices become:
- the paramagnetic current $-\frac1c\bj\cdot\bA\,\tau_0$;
- the diamagnetic term $\frac{e^2}{2mc^2}\rho A^2\,\tau_3$;
- the charge $e\rho\,\tau_3$, via continuity (E.428′).

In Coulomb-gauge form, which is exact for $c\to\infty$ (E.441–445), the self-energy is (schematic; prefactors to be fixed against G.208–210)

$$\hat\Sigma(\bk,i\omega_n)=-\frac{k_BT}{N}\sum_{\bk'm}\Big[\tau_3\hat G(\bk',i\omega_m)\tau_3\,|e\rho_{\bk\bk'}|^2W(\bq,i\nu)+\tau_0\hat G(\bk',i\omega_m)\tau_0\,K^{\rm T}_{\bk\bk'}(i\nu)\Big],$$

$$W=\frac{4\pi}{\bq\cdot\hat\epsilon\cdot\bq},\qquad\hat\epsilon=\hat\epsilon^{\rm e}+\hat\epsilon^{\rm ion}-1,\qquad W=W^{\rm e}+W^{\rm ph}.$$

Here $\bq=\bk'-\bk$ and $i\nu=i\omega_n-i\omega_m$. The additivity of $\hat\epsilon$ is EPP p. 33 (cf. E.432, E.445). $K^{\rm T}=\frac1c\bj^*\cdot(D^t-D^{0t})\cdot\bj$ is of relative order $(v/c)^2$.

For one isotropic mode, $W^{\rm ph}(\bq,i\nu)=-\frac{4\pi}{q^2}\frac{\omega_p^2/\epsilon^2}{\nu^2+\tilde\omega_{\rm LO}^2}$, which is attractive [derived here; cf. self-energy topic §4.4]. Here $\epsilon\to\epsilon^{\rm e}(\bq,i\nu)$, and $\tilde\omega_{\rm LO}$ are the roots of E.453, the $L_\pm$ modes when the crystal is doped. Linearized, G.209 holds with

$$\lambda_{n\bk,n'\bk'}(i\nu)-N_FV_{n\bk,n'\bk'}\;\to\;-N_F\langle n\bk,n'\bk'|W^{\rm e}+W^{\rm ph}|n\bk,n'\bk'\rangle(i\nu),$$

which is the dielectric formulation of pairing (Dolgov 1981).
- **New relative to G.209:** a dynamically screened Fröhlich vertex (E.448a), $\omega_{\rm LO}$ as an output, and coupled LO–plasmon modes.
- **Still missing:** the short-range $g$, and vertex corrections for $\omega_{\rm LO}/E_F\gtrsim1$.

**4.2 2D and cavities [proposal].**
- **Substrates.** For a sheet on a polar substrate, use the layered-stack propagator in place of $W^{\rm e}(\mathbf Q,\omega,z,z')$ (E.461), with the substrate's $\Pi^{\rm ion}$ in its Dyson equation (E.455–472 machinery). The FeSe/SrTiO$_3$ coupling is then a surface-phonon pole of $W(\mathbf Q,i\nu;0,0)$, small-$Q$ by construction.
- **Cavities.** Replace $D^0$ by the cavity Green's dyadic $D^0_{\rm cav}$; E.379 and E.398 keep their form.
  - Sentef's mechanism is $\Omega^{\rm e}$ and $g^{\rm dyn}$ (E.399, E.416) evaluated with $D^{\rm e}_{\rm cav}$.
  - Schlawin's is $K^{\rm T}[D_{\rm cav}]$.
  - Curtis's driven enhancement needs nonequilibrium kinetics, outside EPP.

**4.3 Meissner kernel [derived here; standard].** With Nambu propagators in $\Pi^{\rm e}$, the identity $\Pi^{\rm dia}=-\Pi^{\rm para}(\omega=0)$ behind E.281 fails: $\Pi^t_{jj}(\bq\to0,0)\propto n_s/m\neq0$ (the London kernel). Superconductivity then appears in EPP's photon sector as a photon mass, and E.279 must replace E.281.

**4.4 Soft modes [proposal].** Add to E.134 the term of $V_{\rm cr}$ linear in $Q$. This gives a bare short-range vertex $g^{\rm S,b}$, and $F^0$ must then be built from electron-free force constants:

$$F=F^0+F^0\big(\Omega^{\rm e}+\Omega^{\rm S}+\Omega^{\rm LS}+\Omega^{\rm SL}\big)F,\qquad\Omega^{\rm S}_{\sigma\sigma'}=\sum_{\bk nm}g^{\rm S,b}_{nm\sigma}\,\chi^0_{n\bk,m\bk+\bq}(\omega)\,\big[g^{\rm S}_{nm\sigma'}(\omega)\big]^*.$$

As in E.412 and G.141, one vertex is bare and one screened; the cross terms $\Omega^{\rm LS}$ and $\Omega^{\rm SL}$ mix the polar and nonpolar channels. The instability criterion is

$$\tilde\omega^2_\sigma(\mathbf Q;T)=\omega^2_\sigma(\mathbf Q;T)+\frac{2\omega_\sigma}{\hbar}\,{\rm Re}\,\Omega_{\sigma\sigma}(\mathbf Q,\tilde\omega;T)=0,$$

with $T$ entering through $f_{n\bk}$ in $\chi^0$ and through an anharmonic $\omega^2_\sigma(T)$ from SSCHA. Keeping $\bG\neq0$ in $D^{\rm e}$ (E.454a, with the on-site Ewald term removed) restores the Lorentz field, and with it Cochran-type softening.

**4.5 Condensates [assessment].**
- **Phonons.** A "phonon condensate" $\langle Q_\sigma(\mathbf Q)\rangle\neq0$ is just the frozen distortion of 4.4.
- **Photons.** Photon condensation needs a zero transverse eigenvalue of $D^{-1}(\bq,0)=[D^0]^{-1}-\Pi_{jj}(\bq,0)$. E.281 and E.311 give $\Pi_{jj}(\bq\to0,0)=0$, so the uniform no-go theorem holds identically. At finite $\bq$ the condition is $\chi_{\rm orb}(\bq)>1/4\pi$ [derived here; prefactors to check].
- **Polaritons.** Put $\Pi^{\rm BSE}_{jj}$ (long-range exchange removed, [§Z.1](../../sections/Z01-outlook.html)) and $D^0_{\rm cav}$ into E.379. The poles of $D$ then give the polariton dispersion, the Hopfield weights and the radiative widths. Relaxation (the bottleneck) needs exciton–phonon vertices, and condensation itself needs a kinetic or Gross–Pitaevskii layer.

## 5. Pitfalls

- **Double counting.**
  - DFPT $\omega_\sigma$ already contain Kohn anomalies and CDW softening. Either drop $\Omega^{\rm S}$ or use bare force constants.
  - Subtract $g^{\rm L}$ from the DFPT $g$.
  - Add no empirical $\mu^*$ on top of $W^{\rm e}$.
  - Coulomb is $D^{0l}$ (E.364); never add $v$ separately.
- **Gauge.** Use $\Pi(\omega)-\Pi(0)$ in the normal state, and explicit dia + para in the paired state. Truncated band sums break the f-sum rule and fake a Drude weight or a $q=0$ photon instability.
- **Sum rules.** Impose the acoustic sum rules for $Z^*$ and for the force constants.
- **Migdal.** Include vertex corrections when $\omega_{\rm ph}/E_F\gtrsim1$.
- **Convergence.**
  - Fermi-surface integrals need up to $10^6$ $\bk$ points.
  - Small-$q$ forward scattering needs analytic integration.
  - Matsubara cutoffs renormalize $\mu$ to $\mu^*$.
  - $\mathbf Q_{\rm CDW}$ must lie on the grid, and smearing is not the physical $T$.
- **Regularization.** The $1/\omega^2$ of $D^{0l}$ and the noncommuting limits $\bq,\omega\to0$ need care ([§B.9](../../sections/B09-electrostatic.html)). Subtract $D^{0t}$.

## 6. State of the art (verified)

1. F. Giustino, Rev. Mod. Phys. **89**, 015003 (2017), doi:10.1103/RevModPhys.89.015003. Eliashberg, SCDFT, nonadiabatic phonons.
2. E. R. Margine, F. Giustino, Phys. Rev. B **87**, 024505 (2013), doi:10.1103/PhysRevB.87.024505. Anisotropic Eliashberg with Wannier interpolation.
3. M. Lüders *et al.*, Phys. Rev. B **72**, 024545 (2005), doi:10.1103/PhysRevB.72.024545. SCDFT.
4. R. Akashi, R. Arita, Phys. Rev. Lett. **111**, 057006 (2013), doi:10.1103/PhysRevLett.111.057006. Dynamical Coulomb kernel; plasmon-assisted pairing.
5. X. Lin *et al.*, Phys. Rev. X **3**, 021002 (2013), doi:10.1103/PhysRevX.3.021002. SrTiO$_3$ superconducting at $E_F\approx1.1$ meV.
6. J. J. Lee *et al.*, Nature **515**, 245 (2014), doi:10.1038/nature13894. FeSe/SrTiO$_3$ replica bands.
7. M. A. Sentef, M. Ruggenthaler, A. Rubio, Sci. Adv. **4**, eaau6969 (2018), doi:10.1126/sciadv.aau6969. Cavity phonon polaritons and pairing.
8. F. Schlawin, A. Cavalleri, D. Jaksch, Phys. Rev. Lett. **122**, 133602 (2019), doi:10.1103/PhysRevLett.122.133602 (DOI by APS pattern). Cavity-photon pairing.
9. J. B. Curtis *et al.*, Phys. Rev. Lett. **122**, 167002 (2019), doi:10.1103/PhysRevLett.122.167002. Driven cavity Eliashberg effect.
10. R. Bianco *et al.*, Phys. Rev. Lett. **125**, 106101 (2020), doi:10.1103/PhysRevLett.125.106101. Anharmonic CDW in NbSe$_2$.
11. A. Kogar *et al.*, Science **358**, 1314 (2017), doi:10.1126/science.aam6432. Exciton condensation in TiSe$_2$.
12. G. M. Andolina *et al.*, Phys. Rev. B **100**, 121109(R) (2019), doi:10.1103/PhysRevB.100.121109. No-go theorem.
13. G. M. Andolina *et al.*, Phys. Rev. B **102**, 125137 (2020), doi:10.1103/PhysRevB.102.125137 (DOI from the publisher's file name). Finite-$q$ condensation.
14. J. Kasprzak *et al.*, Nature **443**, 409 (2006), doi:10.1038/nature05131. Polariton BEC.
15. C. Anton-Solanas *et al.*, Nat. Mater. **20**, 1233 (2021), doi:10.1038/s41563-021-01000-8. Polariton condensate in MoSe$_2$.

Also cited:
- Dolgov, Kirzhnits, Maksimov, Rev. Mod. Phys. **53**, 81 (1981), doi:10.1103/RevModPhys.53.81.
- Cochran, Adv. Phys. **9**, 387 (1960), doi:10.1080/00018736000101229.
- Nataf, Ciuti, Nat. Commun. **1**, 72 (2010), doi:10.1038/ncomms1069.
- Klaers *et al.*, Nature **468**, 545 (2010), doi:10.1038/nature09567.
- Z. Wang *et al.*, Nature **574**, 76 (2019), doi:10.1038/s41586-019-1591-7.
- Rzążewski *et al.*, Phys. Rev. Lett. **35**, 432 (1975); DOI not confirmed.
- Faeth *et al.*, Phys. Rev. Lett. **127**, 016803 (2021); DOI not confirmed.
- Gastiasoro, Ruhman, Fernandes, Ann. Phys. 168107 (2020), arXiv:1912.01509; volume and DOI unverified.
- Johannes, Mazin, arXiv:0708.1744; journal reference (Phys. Rev. B **77**, 165135, 2008) unverified.
- Thomas *et al.*, arXiv:1911.01459; reported as J. Chem. Phys. **162** (2025), unverified.

## 7. First steps

1. **Derive (1–2 weeks).** Write EPP's vertices in Nambu form. Show that 4.1 reduces to G.208–210, with $W_{\rm ph}=W_{\rm e}\mathcal DW_{\rm e}$ and the Morel–Anderson $\mu^*$. Prove 4.3 and check the $(v/c)^2$ and $(q_c/2k_F)^2$ estimates.
2. **Model (1–2 weeks).** Take one parabolic band and one LO mode (SrTiO$_3$-like, $E_F=1$–$50$ meV). Solve the linearized gap equation with the full Lindhard $W(\bq,i\nu)$, and find where the split $\lambda/\mu^*$ treatment fails.
3. **Interface (1–2 months).** For an FeSe/SrTiO$_3$-like sheet, build $W(\mathbf Q,i\nu;0,0)$ from §B.10 plus a substrate surface-phonon pole. Compare $\lambda$ and the replica weights with the measured range 0.19–0.5, and vary the substrate.
4. **Soft modes (weeks, then months).** First recover Cochran softening by adding $\bG\neq0$ to $\Omega^{\rm e}$ for SrTiO$_3$/BaTiO$_3$. Then treat monolayer NbSe$_2$: add $\Omega^{\rm S}$ with bare force constants, reproduce DFPT at equal smearing, and add SSCHA.
5. **Polaritons (months, exploratory).** Put MoSe$_2$ in a planar cavity, with $\Pi^{\rm BSE}$ and $D^0_{\rm cav}$. Verify that $\Pi(\omega)-\Pi(0)$ gives no $q=0$ instability while the paramagnetic part alone does.

## 8. Open questions

- With $\omega_{\rm LO}\gg E_F$, which vertex corrections does the Ward identity fix in EPP's current-vertex language?
- Can a dynamical-$W$ replica calculation decide whether the FeSe/SrTiO$_3$ coupling is intrinsic, and how much does the sheet's own plasmon screen it?
- Can a substrate shift $T_{\rm CDW}$ through $\Omega^{\rm e}$, given that $\mathbf Q_{\rm CDW}$ is large?
- **[speculation]** Can a cavity change $T_c$ in the thermodynamic limit other than through near-field environment modes, which amount to a substrate effect?
- **[speculation]** Does any crystal satisfy $\chi_{\rm orb}(\bq)>1/4\pi$?
- Can one EPP+BSE theory treat TiSe$_2$'s combined excitonic and lattice instability?
