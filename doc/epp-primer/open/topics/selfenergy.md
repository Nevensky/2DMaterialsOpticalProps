# Electron self-energies with photon and phonon exchange

Labels used below: **[G]**, Giustino RMP 2017; **[EPP]**, the manuscript (draft v2, with the corrections of §Z.3); **[lit]**, published work in §6; **[derived here]**, a short derivation made for this note that still needs an independent check; **[proposal]**, a suggestion for the manuscript that appears in neither source; **[speculation]**. Units are Gaussian as in EPP. In the self-energy formulas $\hbar=1$, and overall volume normalizations still have to be fixed against (G.157).

## 1. What is missing in EPP

- **No electron Dyson equation.** EPP solves three bosonic Dyson equations: the total field $D=D^0+D^0\Pi_{jj}D$ (E.379–381, §6), the electron-dressed field $D^{\rm e}$ (E.388) and the phonons $F=F^0+F^0\Omega^{\rm e}F$ with $\Omega^{\rm e}=-J^*D^{\rm e}J$ (E.398–399, §7). Electrons appear only inside bubbles: $\Pi^{\rm e}_{jj}$ (E.383) and $\chi^0$ (p. 92). There is no $\Sigma$ and no dressed $G$ (primer §X.2, §Z.1).
- **No Fan–Migdal (FM), no Debye–Waller (DW).** Nothing corresponds to (G.151)–(G.161). The same goes for everything built on them: the zero-point renormalization (ZPR) and Varshni shift (G.166–168), lifetimes (G.170), spectral functions and the cumulant (G.172–174; §A.7, §A.8).
- **The derived vertex is used only once.** The electrodynamic Fröhlich vertices $g^{\rm dyn,0}$ and $g^{\rm dyn}$ (E.415–416; their $c\to\infty$ limits are E p.95′, E p.98′ and E.448a) appear only inside the phonon self-energy (E.413), never on an electron line (§B.8).
- **No second order in the displacements.** Ions couple only through Born dipoles linear in $Q$ (E.237, E.241). The crystal potential is frozen, so there is no $\partial^2V/\partial\tau\partial\tau'$ and hence no DW. The $A^2$ terms of E.237 are second order in the field, not in the displacements. There is no short-range channel either.
- **Temperature.** The "$n_B=0$" rule (EPP §4.5, E.246–250) is exact for harmonic boson self-energies. It fails for $\Sigma$, whose numerators carry $n_{\bq\nu}$ (G.158; §B.5).

## 2. Why it matters

**Magnitudes** [lit]:

- **Diamond.** The direct-gap ZPR is 628 meV with $GW$-corrected couplings, more than 40% above DFT (Antonius 2014; G Fig. 13).
- **Survey of many materials.** ZPRs often exceed 0.3 eV, and reach 0.7 eV, in light-element materials. In polar crystals, agreement with experiment requires the nonadiabatic Fröhlich part (Miglio 2020).
- **Weak-coupling Fröhlich shift, $-\alpha\hbar\omega_{\rm LO}$.** This is about −2.5 meV in GaAs ($\alpha\approx0.07$) and −40 meV in GaN ($\alpha\approx0.44$; primer §A.8). Nery & Allen find −0.05 eV out of a −0.15 eV gap ZPR in zinc-blende GaN. MgO and LiF lie beyond weak coupling (Nery 2018).
- **Lattice screening in $W$.** It can lower $GW$ gaps by more than 15% in highly polar compounds (Botti & Marques 2013). Lambrecht et al. (2017) object that the polar term was not confined in $q$, a sign of how strongly the result depends on dynamics.
- **Light cone.** For GaAs, $q_c=\sqrt{\epsilon^\infty}\,\omega_{\rm LO}/c\approx6\times10^{-5}$ Å$^{-1}$, against the polaron wave vector $q_p=(2m^*\omega_{\rm LO}/\hbar)^{1/2}\approx0.025$ Å$^{-1}$. Only a fraction $\approx(2/\pi)q_c/q_p\sim10^{-3}$ of the Fröhlich shift comes from inside the cone (§A.8, §X.1).

**For an RPA/$GW$/BSE practitioner:**

1. **A $GW$ gap is a clamped-nuclei gap.** Without the ZPR, agreement with experiment at the 0.1–0.5 eV level can be fortuitous, especially in light-element and polar crystals such as hBN.
2. **Temperature-dependent BSE needs input.** It requires $\Sigma^{\rm ep}$ energies and widths (BSE+Allen–Heine, G p. 42). Phonon screening of excitons, its two-particle counterpart, lowers binding energies (Filip 2021).
3. **In 2D the polar channel behaves differently.** The 2D Fröhlich vertex tends to a constant as $Q\to0$ (Sohier 2016; cf. E.472). The adiabatic integrand $\propto|g|^2/Q^2\times Q\,dQ$ then diverges logarithmically, while the nonadiabatic integral converges [derived here]. The same $r_0$ sets the Fröhlich strength and the Rytova–Keldysh exciton screening (§B.10), so encapsulation changes both at once.
4. **The add-on is cheap.** Your code already has $W^{\rm e}$. The lattice adds one pole per polar mode to its head, and EPP derives that pole (Sec. 4.4).

## 3. What standard theory does

- **Exact regrouping** [G]. $G=G^{\rm cn}+G^{\rm cn}\Sigma^{\rm ep}G$ with $\Sigma^{\rm ep}=\Sigma^{\rm FM}+\Sigma^{\rm DW}+\Sigma^{\rm dGW}$ (G.151–155). FM is $i\hbar G\Gamma W_{\rm ph}$ with $W_{\rm ph}=W_{\rm e}\mathcal DW_{\rm e}$ (G.105, G.110), screened at both ends (G.156). DW is the static potential of the smeared nuclei (G.154, G.159). dGW has never been evaluated (G p. 27).
- **Working formulas** [G]. Approximations (i)–(vii) lead to (G.157)–(G.161) and to Allen–Heine (G.164–167). The rigid-ion DW (G.195–198) follows from the translational sum rules (G.193–194), and FM and DW partly cancel.
- **Polar materials** [lit]. The vertex is split as $g=g^{\rm S}+g^{\rm L}$ (Verdi & Giustino 2015) and Wannier-interpolated. The adiabatic formula diverges (G p. 43). The cures are nonadiabatic denominators, analytic small-$q$ integration (Nery & Allen 2016) and spectral functions.
- **Combining with $GW$.** Route (a) takes $G^{\rm cn}$ from $GW$ and adds FM+DW, optionally with $GW$-corrected $g$ (Antonius 2014). Route (b) puts lattice screening into $W$ (Botti & Marques; Lambrecht), which is G.104 used literally.
- **Beyond Dyson.** Options are the cumulant (Nery 2018; Kandolf 2022 for doped crystals), polaron self-trapping (Lafuente-Bartolome 2022) and 2D polarons (Sio & Giustino 2023). Mobilities are reviewed by Poncé (2020); free-carrier screening of the vertex is treated by Macheda (2022).

## 4. How it could enter the EPP formalism

**4.1 Electron Dyson equation [proposal].** The interaction $V'$ (E.237) gives two electron self-energies: a paramagnetic one from two $-\frac1c\bj\cdot\bA$ vertices and a diamagnetic tadpole from $\frac{e^2}{2mc^2}\rho_{\rm e}A^2$:

$$G=G^{\rm KS}+G^{\rm KS}\Sigma^{\rm EPP}G,\qquad \Sigma^{\rm EPP}=\Sigma^{\rm para}+\Sigma^{\rm dia},$$

$$\Sigma^{\rm para}_{nn'\bk}(\omega)=i\sum_m\!\int\!\frac{d\bq}{(2\pi)^3}\!\int\!\frac{d\omega'}{2\pi}e^{i\omega'\eta}\,G_{m\bk+\bq}(\omega-\omega')\Big[-\frac1c\,j^{\mu}_{n\bk,m\bk+\bq}D_{\mu\nu}(\bq,\omega')\big(j^{\nu}_{n'\bk,m\bk+\bq}\big)^*\Big],$$

$$\Sigma^{\rm dia}(\br)=\frac{e^2}{2mc^2}\langle A^2(\br)\rangle,\qquad \langle A_\mu A_\mu\rangle=-ic\!\int\!\frac{d\omega}{2\pi}e^{i\omega\eta}D_{\mu\mu}(\br,\br;\omega).$$

The kernel $K=-\frac1c\,\bj^*\!\cdot D\cdot\bj$ is fixed by the normalization: with $D^{0l}$ (E.374) and continuity it becomes the bare Coulomb $v$ between transition charges.

**4.2 Vector form of $W=W_{\rm e}+W_{\rm ph}$ [derived here from E.379, E.387–388, E.393, E.398].** Resumming the ions gives exactly, within EPP's RPA,

$$D=D^{\rm e}+D^{\rm ph},\qquad D^{\rm ph}=-D^{\rm e}JFJ^{\dagger}D^{\rm e}.$$

The phonon part of $K$ is therefore $\frac1c(\bj^*D^{\rm e}J)\,F\,(J^\dagger D^{\rm e}\bj)$, with screened vertices at both ends, as in (G.156). The phonon self-energy (E.412) instead pairs one bare vertex with one screened vertex. $\Sigma^{\rm para}[D^{\rm e}]$ is "$GW$ with the total, electron-dressed field" and $\Sigma^{\rm para}[D^{\rm ph}]$ is the electrodynamic Fan–Migdal term.

**4.3 Gauge [derived here].** In the temporal gauge, $D^{ll}=-\frac{c}{\omega^2}\bq\bq\,W$ (E.441) combined with continuity (E.428′) gives $K^{ll}=(\Delta E/\omega')^2e^2|\rho|^2W$, with $\Delta E=\varepsilon_{n\bk}-\varepsilon_{m\bk+\bq}$. On shell ($\omega=\varepsilon_{n\bk}$) the excess over the Coulomb kernel cancels the pole of $G$. The Bethe sum rule $\sum_m(\varepsilon_{m\bk+\bq}-\varepsilon_{n\bk})|\rho|^2=q^2/2m$ then turns it into $-\frac{ie^2}{2m}\int q^2W/\omega'^2$, which cancels $\Sigma^{\rm dia,ll}$ exactly; the remainder is odd in $\omega'$ and vanishes.

On shell, therefore, para + dia is gauge invariant. This is the self-energy analogue of E.310. It requires complete band sums and a local Hamiltonian, and off shell the two gauges differ.

**Proposal:** evaluate $\Sigma$ in the Coulomb gauge, building the propagators from EPP's gauge-invariant $\Pi_{jj}$ (E.281, E.311):

$$\Sigma^{\rm EPP}\simeq\Sigma^{\rm L}[W]+\Sigma^{\rm T}\big[D^{t}-D^{0t}\big]+\Sigma^{\rm dia,T},\qquad W=\frac{4\pi}{\bq\cdot\hat\epsilon(\bq,\omega)\cdot\bq},\quad \hat\epsilon=1+\frac{4\pi c}{\omega^2}\Pi_{jj}.$$

$\Sigma^{\rm L}$ is a $GW$-type self-energy with charge vertices $e\rho$. This form is exact for isotropic crystals at any $c$, and in general as $c\to\infty$ (E.445).

**4.4 The phonon pole of $W$ [derived here; checked against the Lyddane–Sachs–Teller (LST) form].**

$$W^{\rm ph}(\bq,\omega)=W-W^{\rm e}=\big[W^{\rm e}(\bq,\omega)\big]^2\sum_{\sigma\sigma'}\frac{(\bq\cdot\bj_\sigma)(\bq\cdot\bj_{\sigma'})^*}{2\Omega_{\rm uc}\sqrt{\omega_\sigma\omega_{\sigma'}}}\,F_{\sigma\sigma'}(\bq,\omega).$$

For one isotropic mode with static screening, this is $\frac{4\pi}{q^2}\frac{\Omega_p^2/(\epsilon^\infty)^2}{\omega^2-\omega_{\rm LO}^2}$: Fröhlich's $|g|^2$ times the LO propagator (§A.5). Inserted into $\Sigma^{\rm L}$ at $\bG=\bG'=0$, it gives (G.157) with the vertex of E.448a. Three things are new compared with (G.157):

- $\omega_{\rm LO}$ is an output (E.454), not an input;
- the screening is $\hat\epsilon^{\rm e}(\bq,\omega)$, so approximation (v) is dropped and doped crystals get the $L_\pm$ modes of E.453;
- the residue $\omega_\sigma/\tilde\omega$ of $F$ fixes the normalization (§X.2).

**4.5 Temperature [proposal; this is the standard spectral form].** Do the Matsubara sum once, with $B=-\frac1\pi{\rm Im}\,K^{R}$:

$$\Sigma^R_{n\bk}(\omega)=\sum_m\!\int\!d\bq\!\int_0^\infty\!d\omega'\,B_{nm}(\bk,\bq,\omega')\Big[\frac{1-f_m+n_B(\omega')}{\omega-\varepsilon_m-\omega'+i\eta}+\frac{f_m+n_B(\omega')}{\omega-\varepsilon_m+\omega'+i\eta}\Big]+\Sigma^{\rm x}.$$

Here $\Sigma^{\rm x}$ is the static exchange part. For $B=|g|^2\delta(\omega'-\omega_{\bq\nu})$ this is (G.158). $D$ itself keeps EPP's rule, since its only temperature dependence is $f(T)$ in $\Pi^{\rm e}$ (E.248). The Bose factors enter only here.

**4.6 Debye–Waller analogue [assessment].** In the Coulomb gauge $\Sigma^{\rm dia}$ involves only $\langle A_t^2\rangle$, whose ionic part lives at $q\lesssim q_c$, on micrometre scales. With EPP's macroscopic ($\bG=\bG'=0$) $D$, any $A^2$ tadpole is a rigid shift of all bands and cannot change a gap. A real DW needs the second-order coupling of (G.159), i.e. a short-range channel in the Lagrangian. Its FM and DW parts must then be added together so that (G.193)–(G.194) hold.

**4.7 Retardation [assessment].**

- **Magnitude.** $\Sigma^{\rm T}$ is of order $(v/c)^2\sim10^{-6}$–$10^{-4}$ of the Coulomb term (§X.1).
- **No intraband photon emission.** Real emission within a band is forbidden by the Cherenkov condition $v<c/\sqrt{\epsilon}$ [derived here].
- **Interband emission belongs to the BSE topic.** It is radiative recombination, which needs a hole, so it is a two-particle (exciton) process. TMD exciton radiative lifetimes are a few ps at 4 K (Palummo 2015).
- **Vacuum subtraction.** The vacuum part $D^{0t}$ is UV divergent (mass renormalization) and is already contained in $m_e$. Subtracting it, $D^t-D^{0t}$, mirrors $\Pi(\omega)-\Pi(0)$.

## 5. Pitfalls

**Double counting.**

- **Force constants.** Bare $\omega_\sigma$ must come from short-range force constants. A DFPT $\omega_{\rm LO}$ in $F$, dressed again by $\Omega^{\rm e}$, counts the LO–TO splitting twice.
- **The polar vertex.** If DFPT $g$ is used for FM, subtract $g^{\rm L}$ (same $\epsilon^\infty$, $Z^*$) before adding EPP's polar part.
- **A $GW$ reference.** If $G^{\rm cn}$ is $GW$, keep only $D^{\rm ph}$ and $D^{{\rm e},t}-D^{0t}$.
- **Static $\epsilon^0$.** Never add FM to a $W$ screened with $\epsilon^0$.
- **Hartree terms.** These are already in $V^{\rm KS}$.

**Screening.** EPP's $\hat\epsilon^{\rm e}$ is the RPA head, without local fields or xc. $GW$ codes include local fields and DFPT uses RPA+xc. A mismatch distorts $W^{\rm ph}$ at intermediate $q$.

**Vertices.** FM takes two screened vertices and the phonon self-energy one bare and one screened. Mixing them up repeats the miscount of §B.8.

**Gauge.** The temporal gauge needs the $A^2$ term, complete band sums and a symmetric $\omega'\to0$ prescription (its $1/\omega'^2$ pole; noncovariant gauges are reviewed by Leibbrandt 1987). Truncated bands leave spurious terms, the analogue of the fake Drude weight of §B.5. Nonlocal pseudopotentials and scissors break E.428′ and the Bethe sum rule.

**Sum rules.**

- the acoustic sum rule for Born charges, $\sum_aZ^*_a=0$;
- the basis phase of E p.86a′, which removes spurious acoustic coupling;
- the translational sum rules (G.193–194) once DW is included.

**Convergence.**

- **3D.** The $1/q^2$ singularity needs nonadiabatic denominators, analytic small-$q$ integration and interpolated grids. The error decays only as $1/N_q$.
- **2D.** The adiabatic integral diverges logarithmically. Use a truncated Coulomb interaction (Sohier 2016).
- **Broadening.** Replace $\eta$ by the physical $\Gamma$ (G.166).
- **Light line.** The pole of $D^t$ at $q<q_c$ needs analytic integration.
- **Strong coupling.** Use a cumulant beyond $\alpha\sim1$.

## 6. State of the art (verified)

1. F. Giustino, Rev. Mod. Phys. **89**, 015003 (2017), doi:10.1103/RevModPhys.89.015003. The reference framework.
2. P. B. Allen, V. Heine, J. Phys. C **9**, 2305 (1976). Founds Allen–Heine theory, with DW from translational invariance (DOI not confirmed).
3. C. Verdi, F. Giustino, Phys. Rev. Lett. **115**, 176401 (2015), doi:10.1103/PhysRevLett.115.176401. The polar vertex $g^{\rm L}$, EPP's $c\to\infty$ limit.
4. S. Botti, M. A. L. Marques, Phys. Rev. Lett. **110**, 226404 (2013), doi:10.1103/PhysRevLett.110.226404. Lattice polarization in $GW$'s $W$.
5. W. R. L. Lambrecht, C. Bhandari, M. van Schilfgaarde, Phys. Rev. Materials **1**, 043802 (2017), doi:10.1103/PhysRevMaterials.1.043802. Revisits (4) and confines the polar term in $q$.
6. G. Antonius, S. Poncé, P. Boulanger, M. Côté, X. Gonze, Phys. Rev. Lett. **112**, 215501 (2014), doi:10.1103/PhysRevLett.112.215501. $GW$-enhanced coupling, diamond ZPR.
7. J. P. Nery, P. B. Allen, Phys. Rev. B **94**, 115135 (2016), doi:10.1103/PhysRevB.94.115135. Analytic Fröhlich correction (GaN).
8. A. Miglio et al., npj Comput. Mater. **6**, 167 (2020), doi:10.1038/s41524-020-00434-z. Nonadiabatic ZPR survey.
9. J. P. Nery et al., Phys. Rev. B **97**, 115145 (2018), doi:10.1103/PhysRevB.97.115145. Dyson versus cumulant satellites (MgO, LiF).
10. N. Kandolf, C. Verdi, F. Giustino, Phys. Rev. B **105**, 085148 (2022), doi:10.1103/PhysRevB.105.085148. Exact doped Fröhlich solid.
11. J. Lafuente-Bartolome et al., Phys. Rev. Lett. **129**, 076402 (2022), doi:10.1103/PhysRevLett.129.076402 (erratum: **129**, 249901). Allen–Heine meets polaron self-trapping.
12. S. Poncé, W. Li, S. Reichardt, F. Giustino, Rep. Prog. Phys. **83**, 036501 (2020), doi:10.1088/1361-6633/ab6a43. Mobilities in bulk and 2D materials.
13. T. Sohier, M. Calandra, F. Mauri, Phys. Rev. B **94**, 085415 (2016), doi:10.1103/PhysRevB.94.085415 (erratum: PRB **96**, 159904). The 2D Fröhlich interaction in TMDs.
14. F. Macheda, P. Barone, F. Mauri, Phys. Rev. Lett. **129**, 185902 (2022), doi:10.1103/PhysRevLett.129.185902. Doping screens the Fröhlich vertex and the LO–TO splitting.
15. M. R. Filip, J. B. Haber, J. B. Neaton, Phys. Rev. Lett. **127**, 067401 (2021), doi:10.1103/PhysRevLett.127.067401. Phonon screening in the BSE.

Also cited: G. Leibbrandt, Rev. Mod. Phys. **59**, 1067 (1987), doi:10.1103/RevModPhys.59.1067; M. Palummo, M. Bernardi, J. C. Grossman, Nano Lett. **15**, 2794 (2015), doi:10.1021/nl503799t; W. H. Sio, F. Giustino, Nat. Phys. (2023), doi:10.1038/s41567-023-01953-4 (volume and pages unverified).

## 7. First steps

1. **Derive (1–2 weeks).** Write out Sec. 4.1–4.4 with all prefactors. Prove the $D^{\rm ph}$ identity and the on-shell gauge equivalence. Check that $c\to\infty$ gives (G.157) with E.448a and the residue $\omega_\sigma/\tilde\omega$.
2. **Analytic benchmark (days).** For an isotropic Fröhlich crystal (GaAs, GaN, MgO), recover the ZPR $-\alpha\hbar\omega_{\rm LO}$, the §A.8 rate formula and the Bose-factor $T$ dependence. Then confirm that the finite-$c$ corrections are of order $10^{-3}$ or less.
3. **2D test (2–4 weeks).** Add a head-only $W^{\rm ph}(Q,\omega)$, with $r_0$ (E.472), to your $GW$ code for monolayer hBN and MoS$_2$. Show the adiabatic log divergence and the nonadiabatic convergence, and compare with Sohier and Sio & Giustino. Vary the substrate through $r_0$.
4. **Doped crystals (weeks).** Use the dynamical $\hat\epsilon^{\rm e}(\bq,\omega)$ with $L_\pm$ modes in n-GaAs or gated MoS$_2$, and compute lifetimes and satellites against Kandolf and Macheda.
5. **Toward BSE (months).** Add the short-range channel with DW, a cumulant and coupling to BSE.

## 8. Open questions

- **Which gauge?** Off shell the gauges differ. Which one gives the measured spectral function, given that ARPES couples through $\bA$?
- **A consistent DW.** What is the minimal extension of the Lagrangian (short range, second order in $Q$) that keeps gauge invariance and the translational sum rules?
- **Gated 2D crystals.** There the acoustic $\sqrt{Q}$ plasmon overlaps the phonons, so is the split $\Sigma^{\rm e}+\Sigma^{\rm ph}$ still meaningful?
- **Head-only $W^{\rm ph}$.** How large is its error next to a $W^{\rm e}$ with local fields?
- **dGW.** How large is it?
- **[speculation] A sheet's own light cone.** Does it ($q_c=\omega/c$) give any single-particle effect, or only excitonic radiative widths?
