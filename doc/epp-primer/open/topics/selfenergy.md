# Electron self-energies with photon and phonon exchange

Labels used below: **[G]**, Giustino RMP 2017; **[EPP]**, the manuscript (draft v2, with the corrections of §Z.3); **[lit]**, published work in §6; **[derived here]**, a short derivation made for this note that still needs an independent check; **[proposal]**, a suggestion for the manuscript that appears in neither source; **[speculation]**. Units are Gaussian as in EPP. In the self-energy formulas $\hbar=1$, and overall volume normalizations still have to be fixed against (G.157). References marked "unverified" could not be confirmed in the last check (§6).

## 1. What is missing in EPP

- **No electron Dyson equation.** EPP solves three bosonic Dyson equations: the total field $D=D^0+D^0\Pi_{jj}D$ (E.379–381, §6), the electron-dressed field $D^{\rm e}$ (E.388) and the phonons $F=F^0+F^0\Omega^{\rm e}F$ with $\Omega^{\rm e}=-J^*D^{\rm e}J$ (E.398–399, §7). Electrons appear only inside bubbles: $\Pi^{\rm e}_{jj}$ (E.383) and $\chi^0$ (p. 92). There is no $\Sigma$ and no dressed $G$ (primer §X.2, §Z.1).
- **No Fan–Migdal (FM), no Debye–Waller (DW).** Nothing corresponds to (G.151)–(G.161). The same goes for everything built on them: the zero-point renormalization (ZPR) and Varshni shift (G.166–168), lifetimes (G.170), spectral functions (G.171–174) and the cumulant (G §V.B.6, p. 015003-32; primer §A.7, §A.8).
- **The derived vertex is used only once.** The electrodynamic Fröhlich vertices $g^{\rm dyn,0}$ and $g^{\rm dyn}$ (E.415–416; their $c\to\infty$ limits are E p.95′, E p.98′ and E.448a) appear only inside the phonon self-energy (E.413), never on an electron line (§B.8).
- **No electron–two-phonon vertex.** The electron–ion coupling is linear in the displacements: ions enter through Born dipoles, $\mathbf P$ linear in $Q$ (E.237, E.241). The only $Q^2$ term, $\rho_{\rm ion}\rho_{\rm ion}$ in the Coulomb part of E.237 (with $\rho_{\rm ion}=-\nabla\cdot\mathbf P$), is ion–ion and dresses the phonons, not the electrons. The crystal potential seen by the electrons is frozen, so there is no $\partial^2V/\partial\tau\partial\tau'$ and hence no DW. The $A^2$ terms of E.237 are second order in the field, not in the displacements. There is no short-range channel either.
- **Temperature.** The "$n_B=0$" rule (EPP §4.5, E.246–250) is exact for the retarded RPA boson propagators EPP computes, where $T$ enters only through $f$ in the bubbles (E.248). It fails for $\Sigma$, whose numerators carry $n_{\bq\nu}$ (G.158; §B.5).

## 2. Why it matters

**Magnitudes** [lit]:

- **Diamond.** The direct-gap ZPR is 628 meV with $GW$-corrected couplings, more than 40% above DFT (Antonius 2014; G Fig. 13).
- **Survey of many materials.** ZPRs often exceed 0.3 eV, and reach 0.7 eV, in light-element materials. In infrared-active crystals a large share, for many materials the majority, comes from the long-range Fröhlich term, which must be treated nonadiabatically (Miglio 2020).
- **Weak-coupling Fröhlich shift.** For one band edge the weak-coupling estimate is $-\alpha\hbar\omega_{\rm LO}$: about −2.5 meV for the GaAs CBM ($\alpha\approx0.07$) and −40 meV for the GaN CBM ($\alpha\approx0.44$), both from model parameters (primer §A.8). The gap term is $-(\alpha_{\rm e}+\alpha_{\rm h})\hbar\omega_{\rm LO}$, and $\alpha_{\rm h}$ is usually the larger, since hole masses are heavier. For zinc-blende GaN, Nery & Allen find a Fröhlich gap correction of −0.05 eV out of a total −0.15 eV (G p. 44: 45 of 150 meV), split as −0.02 eV (CBM) and +0.03 eV (VBM). These are analytic small-$q$ corrections added to an adiabatic grid calculation, so they are not the same quantity as the model value $-\alpha\hbar\omega_{\rm LO}$.
- **Strong coupling.** MgO ($\alpha\approx1.7$ from the model parameters of primer §A.8) lies outside weak coupling. Nery et al. (2018) use MgO and LiF to compare cumulant and Dyson satellites.
- **Lattice screening in $W$.** It can lower $GW$ gaps by more than 15% in highly polar compounds (Botti & Marques 2013). Lambrecht et al. (2017) object that the polar term was not confined in $q$, a sign of how strongly the result depends on dynamics.
- **Light cone.** For GaAs, $q_c=\sqrt{\epsilon^\infty}\,\omega_{\rm LO}/c\approx6\times10^{-5}$ Å$^{-1}$, against the polaron wave vector $q_p=(2m^*\omega_{\rm LO}/\hbar)^{1/2}\approx0.025$ Å$^{-1}$. The cone holds a fraction $\approx(2/\pi)q_c/q_p\approx1.5\times10^{-3}$ of the Fröhlich-shift integrand, so finite $c$ can affect at most about $10^{-3}$ of the shift. Even inside the cone it changes only the transverse (TO-polariton) channel, because the LO kernel is $c$-independent (primer §A.8, §B.8, §X.1).

**For an RPA/$GW$/BSE practitioner:**

1. **A $GW$ gap is a clamped-nuclei gap.** Without the ZPR, agreement with experiment at the 0.1–0.5 eV level can be fortuitous, especially in light-element and polar crystals such as hBN.
2. **Temperature-dependent BSE needs input.** It requires $\Sigma^{\rm ep}$ energies and widths (BSE+Allen–Heine, G p. 42). Phonon screening of excitons, its two-particle counterpart, lowers binding energies (Filip 2021; unverified, §6).
3. **In 2D the polar channel behaves differently.** The 2D Fröhlich vertex tends to a constant as $Q\to0$ (Sohier 2016; unverified, §6; cf. E.469–472). The adiabatic integrand $\propto|g|^2/Q^2\times Q\,dQ$ then diverges logarithmically, while the nonadiabatic integral converges [derived here]. The sheet's static screening, $\epsilon_{\rm 2D}(Q)\approx\epsilon_{\rm env}+r_0Q$ with $\epsilon_{\rm env}=(\epsilon_{\rm a}+\epsilon_{\rm b})/2$ (§B.10), sets both the Fröhlich vertex and the Rytova–Keldysh exciton potential. Encapsulation changes $\epsilon_{\rm env}$, which enters both at once; $r_0=2\pi\alpha_{\rm 2D}$ is a property of the sheet.
4. **The add-on is cheap.** Your code already has $W^{\rm e}$. The lattice adds one pole per polar mode to its head, and that pole follows from EPP's equations (E.388, E.393, E.398; derived in §4.4 of this note).

## 3. What standard theory does

- **Exact regrouping** [G]. $G=G^{\rm cn}+G^{\rm cn}\Sigma^{\rm ep}G$ with $\Sigma^{\rm ep}=\Sigma^{\rm FM}+\Sigma^{\rm DW}+\Sigma^{\rm dGW}$ (G.151–155). FM is $i\hbar G\Gamma W_{\rm ph}$ with $W_{\rm ph}=W_{\rm e}\mathcal DW_{\rm e}$ (G.105, G.110), screened at both ends (G.156). DW is the static potential of the smeared nuclei (G.154, G.159). dGW has never been evaluated (G p. 27).
- **Working formulas** [G]. Approximations (i)–(vii) lead to (G.157)–(G.161) and to Allen–Heine (G.164–167). The rigid-ion DW (G.195–198) follows from the translational sum rules (G.193–194), and FM and DW partly cancel.
- **Polar materials** [lit]. The vertex is split as $g=g^{\rm S}+g^{\rm L}$ (Verdi & Giustino 2015) and Wannier-interpolated. The adiabatic formula diverges (G p. 43). The cures are nonadiabatic denominators, analytic small-$q$ integration (Nery & Allen 2016) and spectral functions.
- **Combining with $GW$.** Route (a) takes $G^{\rm cn}$ from $GW$ and adds FM+DW, optionally with $GW$-corrected $g$ (Antonius 2014). Route (b) puts lattice screening into $W$ (Botti & Marques; Lambrecht), which is G.104 used literally.
- **Beyond Dyson.** Options are the cumulant (Nery 2018; Kandolf 2022 for the doped Fröhlich model) and polaron self-trapping (Lafuente-Bartolome 2022; unverified, §6). For the 2D long-range vertex and 2D polarons, see the 2D topic note. Mobilities are reviewed by Poncé (2020), and free-carrier screening of the vertex is treated by Macheda (2022) (both unverified, §6).

## 4. How it could enter the EPP formalism

**4.1 Electron Dyson equation [proposal].** The interaction $V'$ (E.237) gives two electron self-energies: a paramagnetic one from two $-\frac1c\bj\cdot\bA$ vertices and a diamagnetic tadpole from $\frac{e^2}{2mc^2}\rho_{\rm e}A^2$:

$$G=G^{\rm KS}+G^{\rm KS}\big(\Sigma^{\rm EPP}-v^{\rm xc}\big)G,\qquad \Sigma^{\rm EPP}=\Sigma^{\rm para}+\Sigma^{\rm dia},$$

$$\Sigma^{\rm para}_{nn'\bk}(\omega)=i\sum_m\!\int\!\frac{d\bq}{(2\pi)^3}\!\int\!\frac{d\omega'}{2\pi}e^{i\omega'\eta}\,G_{m\bk+\bq}(\omega-\omega')\Big[-\frac1c\,j^{\mu}_{n\bk,m\bk+\bq}D_{\mu\nu}(\bq,\omega')\big(j^{\nu}_{n'\bk,m\bk+\bq}\big)^*\Big],$$

$$\Sigma^{\rm dia}(\br)=\frac{e^2}{2mc^2}\langle A^2(\br)\rangle,\qquad \langle A_\mu A_\mu\rangle=-ic\!\int\!\frac{d\omega}{2\pi}e^{i\omega\eta}D_{\mu\mu}(\br,\br;\omega).$$

$\Sigma^{\rm EPP}$ contains the full electronic exchange (through $D^{0l}$) and the RPA correlation (through the longitudinal $D^{\rm e}$, i.e. $\Sigma^{\rm L}[W^{\rm e}]$ of §4.3), so the $v^{\rm xc}$ already in $G^{\rm KS}$ must be subtracted. The kernel is $K=-\frac1c\,\bj^*\!\cdot D\cdot\bj$. With $D^{0l}$ (E.374) and continuity (E.428′) it becomes $(\Delta E/\hbar\omega')^2\,v\,e^2|\rho|^2$, which equals the Coulomb kernel $v\,e^2|\rho|^2$ only at $\omega'=\Delta E/\hbar$. Off that shell the temporal-gauge para term differs from Coulomb exchange, and only para + dia reproduces it on shell (§4.3).

**4.2 Vector form of $W=W_{\rm e}+W_{\rm ph}$ [derived here from E.379, E.387–388, E.393, E.398].** Resumming the ions gives exactly, within EPP's RPA,

$$D=D^{\rm e}+D^{\rm ph},\qquad D^{\rm ph}=-D^{\rm e}JFJ^{\dagger}D^{\rm e}.$$

The phonon part of $K$ is therefore $\frac1c(\bj^*D^{\rm e}J)\,F\,(J^\dagger D^{\rm e}\bj)$, with screened vertices at both ends, as in (G.156). The phonon self-energy (E.412) instead pairs one bare vertex with one screened vertex. $\Sigma^{\rm para}[D^{\rm e}]$ is "$GW$ with the total, electron-dressed field" and $\Sigma^{\rm para}[D^{\rm ph}]$ is the electrodynamic Fan–Migdal term.

**4.3 Gauge [derived here].** In the temporal gauge the longitudinal block of the dressed field is $D^{{\rm e},ll}=-\frac{c}{\omega^2}\partial\partial'D^{\rm e}_{11}$ with $D^{\rm e}_{11}\equiv W^{\rm e}$ (E.436–437; for $c\to\infty$, E.441 and E.445), i.e. $D^{ll}=-\frac{c}{\omega^2}\bq\bq\,W$ in Fourier space. Combined with continuity (E.428′) this gives $K^{ll}=(\Delta E/\omega')^2e^2|\rho|^2W$, with $\Delta E=\varepsilon_{n\bk}-\varepsilon_{m\bk+\bq}$. On shell ($\omega=\varepsilon_{n\bk}$) the excess over the Coulomb kernel, $[(\Delta E/\omega')^2-1]e^2|\rho|^2W$, cancels the pole of $G$ and leaves $(\Delta E/\omega'^2+1/\omega')\,e^2|\rho|^2W$.

- **First term.** The Bethe sum rule for one Bloch state (local $H$) is $\sum_m(\varepsilon_{m\bk+\bq}-\varepsilon_{n\bk})|\rho_{n\bk,m\bk+\bq}|^2=q^2/2m+\bq\cdot\mathbf v_{n\bk}$. The linear term vanishes at band extrema, and in general after $\bq$ integration, because $W(\bq,\omega')=W(-\bq,\omega')$ (time reversal). What remains, $-\frac{ie^2}{2m}\int q^2W/\omega'^2$, cancels $\Sigma^{\rm dia,ll}$ exactly.
- **Second term.** $\sum_m|\rho|^2W/\omega'$ is independent of $n$, since $\sum_m|\rho_{n\bk,m\bk+\bq}|^2=1$. With a symmetric principal-value prescription it vanishes; otherwise it is a rigid shift that cannot change gaps.

On shell, therefore, para + dia is gauge invariant up to that $n$-independent constant. This is the self-energy analogue of E.310. It requires complete band sums and a local Hamiltonian, and off shell the two gauges differ.

**Proposal:** evaluate $\Sigma$ in the Coulomb gauge, building the propagators from EPP's gauge-invariant $\Pi_{jj}$ (E.281, E.311):

$$\Sigma^{\rm EPP}\simeq\Sigma^{\rm L}[W]+\Sigma^{\rm T}\big[D^{t}-D^{0t}\big]+\Sigma^{\rm dia,T}\big[D^{t}-D^{0t}\big],\qquad W=\frac{4\pi}{\bq\cdot\hat\epsilon(\bq,\omega)\cdot\bq},\quad \hat\epsilon=1+\frac{4\pi c}{\omega^2}\Pi_{jj}.$$

$\Sigma^{\rm L}$ is a $GW$-type self-energy with charge vertices $e\rho$; with a KS reference, $v^{\rm xc}$ is subtracted as in §4.1. For isotropic $\Pi$ the split is exact for the full $D$ at any $c$, with $W=4\pi/(q^2\epsilon_L)$; in general it is exact as $c\to\infty$ (E.441, E.445).

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

- **Magnitude.** $\Sigma^{\rm T}$ is of order $(v/c)^2\approx10^{-7}$–$10^{-5}$ of the Coulomb term for electronic (Breit-type) exchange, with band velocities $10^5$–$10^6$ m/s, and of order $(q_c/q)^2\approx10^{-6}$–$10^{-4}$ for phonon-frequency exchange at the polaron scale (§X.1).
- **Photon emission within a band.** On the photon-like branch, real emission is forbidden by the Cherenkov condition, since $v\ll c/\sqrt{\epsilon}$ [derived here]. With the ionic $\Pi$, $\epsilon(\omega)$ diverges just below $\omega_{\rm TO}$, so on the phonon-like lower polariton branch near $\omega_{\rm TO}$ emission is kinematically allowed. That is TO-phonon emission, with a transverse coupling $\propto(\omega/cq)^2$, and it is already contained in the transverse part of $\Sigma^{\rm para}[D^{\rm ph}]$.
- **Interband emission belongs to the optics (BSE) topic.** It is radiative recombination, which needs a hole, so it is a two-particle (exciton) process.
- **Vacuum subtraction.** The vacuum part $D^{0t}$ is UV divergent (mass renormalization) and already contained in $m_e$. It is subtracted in both transverse terms, $\Sigma^{\rm T}[D^t-D^{0t}]$ and $\Sigma^{\rm dia,T}[D^t-D^{0t}]$: $\langle A_t^2\rangle_{\rm vac}$ is also UV divergent, though band independent. This is a vacuum subtraction. It is not the $\Pi(\omega)-\Pi(0)$ of E.281/E.310, which is the gauge-invariance cancellation of the diamagnetic term.

## 5. Pitfalls

**Double counting.**

- **Force constants.** Bare $\omega_\sigma$ must come from short-range force constants. A DFPT $\omega_{\rm LO}$ in $F$, dressed again by $\Omega^{\rm e}$, counts the LO–TO splitting twice.
- **KS reference.** Subtract $v^{\rm xc}$, since $\Sigma^{\rm L}[W^{\rm e}]$ is the electronic $GW$ self-energy. Only with a $GW$ reference ($G^{\rm cn}$ from $GW$) does one keep just $D^{\rm ph}$ and $D^{{\rm e},t}-D^{0t}$.
- **The polar vertex.** If DFPT $g$ is used for FM, subtract $g^{\rm L}$ (same $\hat\epsilon^\infty$, $Z^*$) before adding EPP's polar part. Subtracting $g^{\rm L}$ also removes the S–L cross term $2{\rm Re}(g^{\rm S}g^{{\rm L}*})$, which $W^{\rm ph}$ does not contain, since EPP has no short-range channel. Either keep it explicitly as FM$[g^{\rm S},g^{\rm L}]$ cross terms evaluated with the same $\hat\epsilon$ and $Z^*$, or use EPP's $W^{\rm ph}$ only as a replacement for the $|g^{\rm L}|^2$ part and state that the cross term is kept from DFPT.
- **Static $\epsilon^0$.** Never add FM to a $W$ screened with $\epsilon^0$.
- **Hartree terms.** These are already in $V^{\rm KS}$.

**Screening.** EPP's $\hat\epsilon^{\rm e}$ is the RPA head, without local fields or xc. $GW$ codes include local fields and DFPT uses RPA+xc. A mismatch distorts $W^{\rm ph}$ at intermediate $q$.

**Vertices.** FM takes two screened vertices and the phonon self-energy one bare and one screened. Mixing them up repeats the miscount of §B.8.

**Gauge.** The temporal gauge needs the $A^2$ term, complete band sums and a symmetric $\omega'\to0$ prescription (its $1/\omega'^2$ pole). Truncated bands leave spurious terms, the analogue of the fake Drude weight of §B.5. Nonlocal pseudopotentials and scissors break E.428′ and the Bethe sum rule.

**Sum rules.**

- the acoustic sum rule for Born charges, $\sum_aZ^*_a=0$;
- the basis phase of E p.86a′, which removes spurious acoustic coupling;
- the translational sum rules (G.193–194) once DW is included.

**Convergence.**

- **3D.** The $1/q^2$ singularity needs nonadiabatic denominators, analytic small-$q$ integration and interpolated grids. The error decays only as $1/N_q$, with $N_q$ points per direction ($\propto N_{\rm tot}^{-1/3}$; primer §A.8).
- **2D.** The adiabatic integral diverges logarithmically; use nonadiabatic denominators. Independently, compute $Z^*$, $\epsilon$ and $g$ with a 2D-truncated Coulomb interaction, so that the small-$Q$ vertex is the true 2D one (Sohier et al.; unverified, §6). Integrate the region $Q\lesssim\epsilon_{\rm env}/r_0$ analytically or with mini-BZ averaging, since $W^{\rm ph}$ is concentrated there.
- **Broadening.** Replace $\eta$ by the physical $\Gamma$ (G.166).
- **Light line.** The pole of $D^t$ at $q<q_c$ needs analytic integration.
- **Strong coupling.** For satellites use a second-order cumulant (Nery 2018); Dyson–Migdal gives better quasiparticle masses (Kandolf 2022). For $\alpha\gtrsim2$–3, or self-trapping, use polaron methods (Lafuente-Bartolome 2022; unverified, §6).

## 6. State of the art (verified)

1. F. Giustino, "Electron-phonon interactions from first principles", Rev. Mod. Phys. **89**, 015003 (2017), doi:10.1103/RevModPhys.89.015003. The reference framework.
2. P. B. Allen, V. Heine, "Theory of the temperature dependence of electronic band structures", J. Phys. C **9**, 2305 (1976) (confirmed via Giustino's bibliography; DOI not confirmed). Founds Allen–Heine theory, with DW from translational invariance.
3. C. Verdi, F. Giustino, "Fröhlich electron-phonon vertex from first principles", Phys. Rev. Lett. **115**, 176401 (2015), doi:10.1103/PhysRevLett.115.176401. The polar vertex $g^{\rm L}$, EPP's $c\to\infty$ limit.
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

1. **Derive (1–2 weeks).** Write out Sec. 4.1–4.4 with all prefactors, including the $v^{\rm xc}$ subtraction. Prove the $D^{\rm ph}$ identity and the on-shell gauge equivalence. Check that $c\to\infty$ gives (G.157) with E.448a and the residue $\omega_\sigma/\tilde\omega$.
2. **Analytic benchmark (days).** For an isotropic Fröhlich crystal (GaAs, GaN, MgO), recover the band-edge shift $-\alpha\hbar\omega_{\rm LO}$, the §A.8 rate formula and the Bose-factor $T$ dependence. Then confirm that the finite-$c$ corrections are of order $10^{-3}$ or less and come only from the transverse channel.
3. **2D test (2–4 weeks).** Add a head-only $W^{\rm ph}(Q,\omega)$ to your $GW$ code, built from E.469–472 with the expansion $\epsilon_{\rm 2D}\approx\epsilon_{\rm env}+r_0Q$ of §B.10. Use monolayer hBN as the main test.
   - Use the code's truncated-Coulomb head and mini-BZ averaging for $Q\lesssim\epsilon_{\rm env}/r_0$.
   - Show the adiabatic log divergence and the nonadiabatic convergence.
   - Benchmark the vertex $|g(Q)|$ implied by $W^{\rm ph}$ against Sohier 2016, and the 2D-Fröhlich ZPR against the analytic nonadiabatic integral.
   - Vary the substrate through $\epsilon_{\rm env}$ at fixed $r_0$, with $W^{\rm e,2D}=2\pi/[Q(\epsilon_{\rm env}+r_0Q)]$. This is an addition beyond EPP's E.461, which is for a freestanding sheet.
   - Treat MoS$_2$ only as a weak-polar control (the Fröhlich part of the e–ph shift only). Its e–ph renormalization is dominated by short-range couplings, so a comparison with full ZPRs would mislead.
4. **Doped crystals (weeks).** Use the dynamical $\hat\epsilon^{\rm e}(\bq,\omega)$ with $L_\pm$ modes in n-GaAs or gated MoS$_2$, and compute lifetimes and satellites. Compare with Kandolf 2022 (a 3D model benchmark) and Macheda 2022.
5. **Toward BSE (months).** Add the short-range channel with DW, a cumulant and coupling to BSE.

## 8. Open questions

- **Which gauge?** Off shell the gauges differ. Which one gives the measured spectral function, given that ARPES couples through $\bA$?
- **A consistent DW.** What is the minimal extension of the Lagrangian (short range, second order in $Q$) that keeps gauge invariance and the translational sum rules?
- **Gated 2D crystals.** There the acoustic $\sqrt{Q}$ plasmon overlaps the phonons, so is the split $\Sigma^{\rm e}+\Sigma^{\rm ph}$ still meaningful?
- **Head-only $W^{\rm ph}$.** How large is its error next to a $W^{\rm e}$ with local fields?
- **dGW.** How large is it?
- **[speculation] A sheet's own light cone.** Does it ($q_c=\omega/c$) give any single-particle effect, or only excitonic radiative widths?

## 9. Review responses

All 23 points were checked against EPP (E.237–250, E.281, E.374, E.428, E.436–445, E.461–472), Giustino (pp. 30–32, 43–44 and the bibliography) and primer §A.8, §B.10 and §X.1. Every point was accepted. This pass could not re-search citations: the shared web-search budget was used up, and doi.org, Crossref and the publishers were unreachable.

- **Fröhlich share of the ZPR (major).** Accepted. Giustino p. 44 confirms 45 of 150 meV for zb-GaN. §2 now separates the band-edge model value ($-\alpha\hbar\omega_{\rm LO}$, GaN CBM ≈ −40 meV; recomputed α = 0.44 from the §A.8 parameters) from the gap term $-(\alpha_{\rm e}+\alpha_{\rm h})\hbar\omega_{\rm LO}$ and from Nery & Allen's analytic correction. The CBM/VBM split (−0.02/+0.03 eV) is taken from the reviewer's reading of their abstract. It is consistent with the confirmed −0.05 eV total. Miglio 2020 is now quoted as "for many materials".
- **$K^{0l}$ off shell (major).** Accepted and rederived: $-\frac1c\bj^*D^{0l}\bj$ with E.374 and E.428 gives $(\Delta E/\hbar\omega')^2v|\rho|^2$. §4.1 now agrees with §4.3.
- **$v^{\rm xc}$ double counting (major).** Accepted: the Dyson equation is now $G^{\rm KS}(\Sigma^{\rm EPP}-v^{\rm xc})G$, with a new "KS reference" bullet in §5 that is merged with the $GW$ reference bullet.
- **S–L cross term (major).** Accepted; the two options are stated in §5.
- **Substrate through $\epsilon_{\rm env}$, not $r_0$ (major).** Accepted and checked against §B.10 and E.469–470: E.461 has no environment, and E.470 is $1-v_Q\Pi^{\rm e,0}_{\rho\rho}$. §2 item 3 and §7 step 3 have been rewritten.
- **2D convergence vs truncation (major).** Accepted: these are now two separate requirements.
- **Bethe sum rule $+\bq\cdot\mathbf v_{n\bk}$, $n$-independent remainder, E.436–437 citation.** Accepted and rederived ($e^{-i\bq\br}He^{i\bq\br}-H=\bq\cdot\mathbf p/m+q^2/2m$; completeness gives $\sum_m|\rho|^2=1$). §4.3 has been restructured accordingly.
- **Attribution of the $W^{\rm ph}$ pole, cumulant citation (G §V.B.6, p. 32), the linear-in-$Q$ statement (the $\rho_{\rm ion}\rho_{\rm ion}$ term of E.237 checked), the $n_B$ rule, $N_q$ per direction (checked in §A.8: $q_{\rm min}=q_D/N_q$).** All accepted.
- **$(v/c)^2$ vs $(q_c/q)^2$; light-cone fraction as an upper bound ($(2/\pi)q_c/q_p=1.5\times10^{-3}$ recomputed); Cherenkov on the lower polariton branch; vacuum subtraction in $\Sigma^{\rm dia,T}$ (E.281 checked: it is the dia/para(0) cancellation).** All accepted.
- **§7 step 3 and step 4, strong-coupling advice, Kandolf description, MgO (α = 1.725 recomputed).** All accepted. hBN is kept as the main test.
- **Citations.** Accepted. The two errata, Leibbrandt 1987, Palummo 2015 (with its "few ps at 4 K" figure) and Sio & Giustino 2023 were dropped, since none could be confirmed and none carried a load-bearing claim. Five core references are kept but marked unverified in §6 and in the text. The verified list now has 10 papers. The reviewer's "11" names only these 10. Titles for refs 2, 3, 4, 6 and 7 were checked here against Giustino's bibliography; Allen & Heine is counted as verified on that basis, with its DOI still unconfirmed. Refs 5 and 8–10 rest on the earlier search-based confirmation, which the reviewer repeated. Nery 2018 now has its full title, and Kandolf 2022 has its title and a corrected description.
