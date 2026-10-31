# SYNTHESIS: Giustino RMP 89, 015003 (2017) → EPP ("Theory of the Electron-Phonon-Photon coupling")

Sources and tags used below:
- Note files: [G1]–[G6] = Giustino notes (journal page 015003-N = PDF page N); [E1]–[E9] = EPP notes (PDF page = printed page).
- "G Eq. n" = equation n of Giustino; "E Eq. n" = equation n of the EPP manuscript; "E p. n" = unnumbered equation on EPP page n.
- [derived] = derived or checked by the note authors / by me, not printed in either paper. [lit] = literature value not taken from either paper (flag on the page). [ctx] = context equation (e.g. Gonze–Lee) that the paper cites but does not print.
- PDF pages re-checked against the notes for this synthesis: EPP pp. 20–22, 24–26, 33–35, 83–87, 89–91, 95–102; Giustino pp. 3–5, 17–19, 22–25, 35. All matched the notes; the corrections listed in Sec. 6 are confirmed.
- Audience: grad students and researchers in computational condensed matter. The reader computes RPA/GW/BSE optics of 2D materials, so W, χ⁰, ε and Rytova–Keldysh screening are familiar; phonon field theory and photon propagators may not be.

---

## 1. Part A narrative: Giustino (2017), step by step

**A1. The object (G Eq. 1) [G1].** The whole review is about the canonical Hamiltonian
$$\hat H=\sum_{n\mathbf k}\varepsilon_{n\mathbf k}\hat c^\dagger_{n\mathbf k}\hat c_{n\mathbf k}+\sum_{\mathbf q\nu}\hbar\omega_{\mathbf q\nu}\big(\hat a^\dagger_{\mathbf q\nu}\hat a_{\mathbf q\nu}+\tfrac12\big)+N_p^{-1/2}\sum_{\mathbf k\mathbf q mn\nu}g_{mn\nu}(\mathbf k,\mathbf q)\,\hat c^\dagger_{m\mathbf k+\mathbf q}\hat c_{n\mathbf k}\big(\hat a_{\mathbf q\nu}+\hat a^\dagger_{-\mathbf q\nu}\big)+[\text{DW term}]$$
It has no explicit e–e interaction and no electromagnetic field. Coulomb physics is hidden inside $\varepsilon_{n\mathbf k}$, $\omega_{\mathbf q\nu}$ and the *screened* $g$. The review's two goals are (i) a recipe for the parameters (DFT/DFPT) and (ii) a field-theoretic justification of when Eq. 1 is right.

**A2. Lattice dynamics (G Eqs. 12–22, App. B) [G1, G2].** Force constants $C_{\kappa\alpha p,\kappa'\alpha'p'}$ give the dynamical matrix $D^{\rm dm}(\mathbf q)=(M_\kappa M_{\kappa'})^{-1/2}\sum_pC_{\kappa\alpha0,\kappa'\alpha'p}e^{i\mathbf q\cdot\mathbf R_p}$ (G Eq. 14). Its eigenpairs $\omega^2_{\mathbf q\nu}$, $e_{\kappa\alpha,\nu}(\mathbf q)$ define normal modes. Quantization gives $\Delta\tau_{\kappa\alpha p}=(M_0/N_pM_\kappa)^{1/2}\sum_{\mathbf q\nu}e^{i\mathbf q\cdot\mathbf R_p}e_{\kappa\alpha,\nu}(\mathbf q)\,l_{\mathbf q\nu}(\hat a_{\mathbf q\nu}+\hat a^\dagger_{-\mathbf q\nu})$ (G Eq. 20), with zero-point length $l_{\mathbf q\nu}=[\hbar/(2M_0\omega_{\mathbf q\nu})]^{1/2}$ (G Eq. 21) and $\hat H_p$ (G Eq. 22).

**A3. A century of "g" (G Eqs. 2–7, 11) [G1].** The matrix element is "displacement × gradient of an effective potential", $g=-i(\hbar/2N_pM\omega)^{1/2}\,\mathbf q\cdot\mathbf e\,V_0$ (G Eq. 2). History is the story of $V_0$:
- Bloch: a constant.
- Nordheim: the bare ion $V_\kappa(q)\propto q^{-2}$, too strong.
- Mott–Jones: truncated at the WS cell.
- Bardeen: $V_\kappa/\epsilon(q)$ with Lindhard screening (G Eqs. 3–4).
- Bardeen–Pines: dynamic $\epsilon(q,\omega_{\mathbf q\nu})$ (G Eq. 5).
- Deformation potential $E_1=\Omega\,\partial\varepsilon/\partial\Omega$ (G Eq. 6).
- Fröhlich's $1/q$ polar vertex (G Eq. 7).
- The DFT matrix element (G Eq. 11).

**A4. DFT and DFPT (G Eqs. 23–61) [G2].** $g_{mn\nu}(\mathbf k,\mathbf q)=\langle u_{m\mathbf k+\mathbf q}|\Delta_{\mathbf q\nu}v^{\rm KS}|u_{n\mathbf k}\rangle_{\rm uc}$ (G Eq. 38), with $\Delta_{\mathbf q\nu}v^{\rm KS}$ from G Eqs. 33–35. The screened perturbation solves the Sternheimer problem (G Eq. 47), which is equivalent to $\partial v^{\rm KS}=(\epsilon^{{\rm H}xc})^{-1}\partial v^{\rm en}$ (G Eq. 53), $\epsilon^{{\rm H}xc}=1-(v^{\rm C}+f^{xc})\chi^0$ (G Eq. 54, "test electron"), versus $\epsilon^{\rm H}=1-v^{\rm C}\chi^0$ (G Eq. 55, "test charge"/RPA). In jellium this reproduces Bardeen (G Eqs. 59–61). DFPT screening is static, i.e. adiabatic.

**A5. Many-body foundations: Hedin–Baym (G Sec. IV, Eqs. 62–122, Table I) [G3].**
- Interactions are instantaneous Coulomb, $v(12)=v(\mathbf r_1,\mathbf r_2)\delta(t_1-t_2)$. Schwinger functional derivatives give $W=\epsilon^{-1}v=v+vPW$ (G Eqs. 86, 90), with $P=\delta\langle\hat n\rangle/\delta V_{\rm tot}$ containing electrons *and* nuclei.
- Clamped nuclei give Hedin's GW equations, with $W_{\rm e}=v+vP_{\rm e}W_{\rm e}$ (G Eq. 94).
- Mobile nuclei give $W=W_{\rm e}+W_{\rm ph}$ with $W_{\rm ph}=W_{\rm e}DW_{\rm e}$ (G Eqs. 104–105). The vertex is the dynamically screened gradient of the nuclear potential, $\epsilon_{\rm e}^{-1}\nabla V_\kappa$ (G Eq. 110).
- Phonons obey $[\mathbf M\omega^2-\Pi(\omega)]\mathbf D=1$ (G Eq. 118), with $\Pi=Z_\kappa Z_{\kappa''}\partial\partial'[W_{\rm e}(\omega)-W_{\rm e}(0)]$-type structure (G Eq. 121) and the acoustic sum rule (G Eq. 120).
- Setting $D=0$ recovers Hedin. "Retardation" here means only the frequency dependence of $\epsilon_{\rm e}$.

**A6. Phonons from first principles (G Sec. V.A, VII; Eqs. 123–146) [G4, G5].**
- $\mathbf D=[\mathbf M\omega^2-\Pi]^{-1}$ (G Eq. 123) with $\Pi=\Pi^{\rm A}+\Pi^{\rm NA}(\omega)$ (G Eq. 125). $\Pi^{\rm A}$ is the DFPT force-constant matrix (G Eq. 128). In the adiabatic eigenmode basis $D^{\rm A}_{\mathbf q\nu}=2\omega_{\mathbf q\nu}/(\omega^2-\omega_{\mathbf q\nu}^2)$ (G Eqs. 132, 134).
- The nonadiabatic Dyson equation (G Eqs. 135–136) gives $\tilde\Omega^2=\omega^2+2\omega\Pi^{\rm NA}(\tilde\Omega)$ (G Eq. 137) and the linewidth (G Eqs. 138, 146).
- The ab initio self-energy (G Eq. 141) has **one bare vertex $g^{\rm b}$ and one screened vertex $g$**, with $\chi^0(\omega)-\chi^0(0)$ in between; approximations (i)–(v) take it to G Eq. 145.
- Sec. VII: the nonadiabatic Kohn anomaly of doped graphene.

**A7. Electrons (G Sec. V.B; Eqs. 147–174) [G4].**
- Exact rewrite: $\Sigma^{\rm ep}=\Sigma^{\rm FM}+\Sigma^{\rm DW}+\Sigma^{\rm dGW}$ (G Eqs. 152–155). Fan–Migdal is $GW_{\rm ph}$ with **two screened vertices** (G Eq. 156). It becomes G Eqs. 157–158 (finite $T$, retarded).
- The Debye–Waller term is G Eqs. 159–161.
- Consequences: AH band renormalization (G Eq. 166), Varshni (G Eq. 168), lifetimes (G Eq. 170), spectral functions with kinks and satellites (G Eqs. 171–174, Fig. 4), and the cumulant expansion.

**A8. Interpolation and the polar problem (G Sec. VI; Eqs. 175–185) [G5].** Wannier interpolation of $g$ fails for polar materials, because $g\sim1/q$. The fix is to subtract and add the analytic long-range part (Verdi–Giustino; G Eq. 180):
$$g^{\mathcal L}_{mn\nu}(\mathbf k,\mathbf q)=i\frac{4\pi}{\Omega}\frac{e^2}{4\pi\varepsilon_0}\sum_\kappa\Big(\frac{\hbar}{2M_\kappa\omega_{\mathbf q\nu}}\Big)^{1/2}\sum_{\mathbf G\neq-\mathbf q}\frac{(\mathbf q+\mathbf G)\cdot\mathbf Z^*_\kappa\cdot\mathbf e_{\kappa\nu}(\mathbf q)}{(\mathbf q+\mathbf G)\cdot\boldsymbol\epsilon^\infty\cdot(\mathbf q+\mathbf G)}\langle\psi_{m\mathbf k+\mathbf q}|e^{i(\mathbf q+\mathbf G)\cdot(\mathbf r-\boldsymbol\tau_\kappa)}|\psi_{n\mathbf k}\rangle_{\rm sc}$$
($N_p$ dropped; see Sec. 6.) Its sibling is the Gonze–Lee nonanalytic dynamical-matrix term that produces LO–TO splitting [ctx]. Fermi-surface harmonics (G Eq. 181 ff.) and the double-delta linewidth (G Eq. 185) follow.

**A9. Photoemission and metals (G Sec. VIII; Eqs. 186–192) [G5, G6].** These cover the spectral representation of $\Sigma^{\rm FM}$, $\alpha^2F$, $\lambda$, velocity renormalization $1/(1+\lambda)$, and the heat-capacity enhancement.

**A10. Applications (G Secs. IX–XII; Eqs. 193–212) [G6].**
- Allen–Heine in practice, and the divergence of adiabatic AH in polar materials ("Fröhlich catastrophe").
- Frozen-phonon and Williams–Lax approaches (G Eq. 201).
- Phonon-assisted absorption (G Eq. 200): the photon is a **classical** perturbing field and only electrons couple to it.
- Mobility and the Boltzmann transport equation (G Eqs. 202–203).
- Superconductivity: Migdal–Eliashberg, McMillan–Allen–Dynes, SCDFT (G Eqs. 204–207).
- $g$ beyond LDA.

**Takeaway for Part A.** All coupling is mediated by the instantaneous, longitudinal Coulomb interaction. Retardation enters only via $\epsilon_{\rm e}(\omega)$. Photons are external and classical. Polar coupling and LO–TO splitting are added analytically through $Z^*$ and $\epsilon^\infty$.

---

## 2. Part B narrative: EPP, step by step

**B1. Lattice (E Sec. 1, Eqs. 1–91) [E1, E2].** This is the same harmonic lattice as Giustino. The elastic matrix gives the dynamical matrix (E Eq. 15), its eigenvectors, and ladder operators. The displacement operator $Q_{ia\mu}=N^{-1/2}\sum\sqrt{\hbar/2m_a\omega_\alpha}\,e_{a\mu\alpha}(a+a^\dagger_{-\mathbf q})e^{i\mathbf q\mathbf R_i^0}$ (E Eqs. 83/242) equals G Eq. 20 without $M_0$.

**B2. Rigid-ion e–ph coupling (E Sec. 2, Eqs. 92–109) [E3].** Expanding rigid local pseudopotentials gives the **bare** vertex $g^0$ (E Eq. 107, boxed) and the "bare" model Hamiltonian (E Eq. 109). Page 21 states the double-counting issue: DFT force constants already contain the screened short-range interactions, so only short-range (analytic) force constants may be combined with explicitly mediated long-range coupling.

**B3. Polar coupling electrostatically (E Sec. 3, Eqs. 110–131) [E3].**
- Born charges: E Eqs. 110–111. In EPP $Z^*$ has the dimension of charge.
- Dipole potential in an anisotropic dielectric: $\phi(\mathbf q)=-\frac{i}{\epsilon_0}\frac{\mathbf p\cdot\mathbf q}{\mathbf q\cdot\hat\epsilon^\infty\cdot\mathbf q}$ (E Eq. 120).
- Mapping onto Sec. 2 gives $g^{\mathcal L}$ (E Eq. 131), stated to be "exactly equal" to Verdi–Giustino Eq. 4 (= G Eq. 180).
- $\hat\epsilon^\infty$ is inserted ad hoc. It is regenerated later from $\Pi^{\rm e}$.

**B4. One Lagrangian for everything (E Sec. 4, Eqs. 132–143) [E3, E4].**
- Temporal gauge $\Phi=0$, $\mathbf E=-\frac1c\partial_t\mathbf A$ (E Eqs. 132–133).
- The Lagrangian (E Eq. 134) has six parts: electrons in the frozen $V_{\rm cr}$; harmonic ions; the free field $\frac1{8\pi}\int(\mathbf E^2-\mathbf B^2)$; $\frac ec\mathbf v_e\cdot\mathbf A$; and $\frac1c\dot{\mathbf p}_{ia}\cdot\mathbf A(\mathbf R_i^0)$ with $\mathbf p=\hat Z^*\mathbf Q$ (E Eq. 135).
- Euler–Lagrange gives the Lorentz force for electrons and $m_a\ddot Q+kQ=Z^{*T}\mathbf E$ for ions (E Eq. 143).
- There is no explicit e–ph term. Electrons and phonons talk only through $\mathbf A$.

**B5. Classical polar crystal (E Eqs. 144–165) [E4].**
- Diatomic, long-wavelength: the bare optical mode is $\omega_O^2=k/M$, and the ionic plasma frequency is $\omega_p^2=4\pi Z^{*2}/(\Omega_{\rm uc}M)$ (E Eq. 152).
- With the Lorentz local field: $\epsilon(\omega)=\epsilon_\infty-\omega_p^2/(\omega^2-\omega_O^2+\omega_p^2/3)$ (E Eq. 158), factorized as $\epsilon_\infty(\omega^2-\tilde\omega_{LO}^2)/(\omega^2-\tilde\omega_{TO}^2)$ (E Eqs. 159–160).
- LO modes are zeros of $\epsilon$ (longitudinal); TO modes are poles (transverse). LST follows (E Eq. 165).

**B6. Fields, Hamiltonian, photons (E Eqs. 166–237) [E5].**
- Field Euler–Lagrange gives Maxwell with source $\mathbf j_e+\dot{\mathbf P}$.
- Canonical momenta (E Eqs. 182–184) give the minimal-coupling Hamiltonian and its four couplings: electron para/dia and ion para/dia (E Eq. 186).
- Quantization in a triad $\mathbf e_{1,2,3}(\mathbf q)$ gives transverse photons $\hbar cq$ (E Eq. 214). The longitudinal sector $H^0_{\rm Lph}=\frac1{8\pi}\int\mathbf E_l^2$ (E Eq. 227) becomes, through Gauss's law, the full instantaneous Coulomb interaction of electrons and ionic bound charges (E Eq. 235).
- Final split: $H_0$ = electrons + bare phonons + transverse photons (E Eq. 236). $V'$ = $-\frac1c\mathbf j_e\cdot\mathbf A+\frac{e^2}{2m c^2}\rho_e\mathbf A^2-\frac1c Z^*\pi\cdot\mathbf A/m_a+\frac1{2c^2}Z^*Z^*AA/m_a+$Coulomb (E Eq. 237).

**B7. Free propagators (E Secs. 4.5–4.7, Eqs. 246–269) [E6].**
- The "$n_B=0$" rule (Sec. 4.5).
- Phonon: $F_\sigma=\frac1\hbar\frac{2\omega_\sigma}{\omega^2-\omega_\sigma^2+i\eta}$ (E Eq. 256).
- Transverse photon: $D^{0t}_{\mu\nu}=-4\pi c\,[\delta_{\mu\nu}-q_\mu q_\nu/q^2]/(\omega^2-c^2q^2+i\eta)$ (E Eq. 261).
- Longitudinal photon in the interaction picture: $D^{0l}=0$ (E Eq. 263). Its physics is "hidden" in E Eq. 235.

**B8. RPA photon self-energy (E Sec. 5, Eqs. 270–311) [E6, E7].**
- Transverse Dyson equation (E Eq. 277); $\Pi=\Pi^{\rm e}+\Pi^{\rm ion}$ (E Eq. 278), each dia + para (E Eqs. 279–280).
- Gauge invariance: $\Pi=\Pi^{\rm para}(\omega)-\Pi^{\rm para}(0)$ (E Eq. 281).
- Electronic para term with current vertices (E Eqs. 282–284) = Adler–Wiser with $\rho\to\mathbf j$.
- Ionic current vertex $j_{\alpha\sigma}=\sum_{a\mu}Z^*_{a\alpha\mu}e_{a\mu\sigma}/\sqrt{m_a}$ (E Eq. 295); ionic para term (E Eqs. 297/299); dia term (E Eq. 304).
- Completeness of the eigenvectors gives dia $=-$para$(0)$ (E Eqs. 309–310) and $\Pi^{\rm ion}=-\frac{\omega^2}{c\Omega_{\rm uc}}\sum_\sigma\frac{j_{\alpha\sigma}j^*_{\beta\sigma}}{\omega^2-\omega_\sigma^2+i\eta}$ (E Eq. 311).
- 2D forms: E Eqs. 300, 305.

**B9. Charge conservation and the total propagator (E Sec. 6, Eqs. 312–385) [E7, E8].**
- A $4\times4$ propagator (scalar $v$ plus $3\times3$ transverse; E Eqs. 317–319).
- Connection formulas from continuity (E Eqs. 341, 344, 347–350), e.g. $\Pi_{\rho\rho}=-\frac{c}{\omega^2}\partial_\mu\partial'_\nu\Pi_{jj\mu\nu}$. These turn every density vertex into a current vertex and $v$ into $D^{0l}=-\frac{4\pi c}{q^2}\frac{q_\mu q_\nu}{\omega^2}$ (E Eqs. 373–374).
- Total free propagator: $D^0=D^{0l}+D^{0t}$ = vacuum Maxwell propagator (E Eq. 375).
- Dyson: $D=D^0+D^0\Pi_{jj}D$ (E Eqs. 379–385).

**B10. Phonons dressed by photons (E Sec. 7, Eqs. 386–454) [E9].**
- Resum electrons into $D^{\rm e}=D^0+D^0\Pi^{\rm e,0}D^{\rm e}$ (E Eq. 388).
- Dynamical ionic vertex $J_{\mu\sigma}=\sqrt{\hbar/(2c\Omega_{\rm uc}\omega_\sigma)}\,\omega\,j_{\mu\sigma}$ (E Eq. 392).
- Phonon Dyson $F=F^0+F^0\Omega^{\rm e}F$ (E Eq. 398), $\Omega^{\rm e}=-J^*D^{\rm e}J$ (E Eq. 399), diagonal form (E Eqs. 404–405).
- Sec. 7.1: at $q=0$ without electrons, $\tilde\omega^2=\omega_\sigma^2+\frac{4\pi}{\Omega_{\rm uc}}\sum_\alpha|j_{\alpha\sigma}|^2$ (E Eq. 410).
- Sec. 7.2: frequency-dependent bare and screened e–ph vertices (E Eqs. 415–416), $\Delta\Omega^{\rm e}=\sum g^{\rm dyn}\chi^0g^{\rm dyn,0}$ (E p. 92). Their $c\to\infty$ limits are Verdi–Giustino-like vertices (E pp. 95, 98), with $\epsilon^\infty\to\hat\epsilon^{\rm e}(\mathbf q,\omega)$ (E Eqs. 431–432, 446).
- Sec. 7.3: $\tilde\omega^2=\omega_\sigma^2+W^{\rm e}(\mathbf q,0)|\mathbf q\cdot\mathbf j_\sigma|^2/\Omega_{\rm uc}$ (E Eq. 454) = first-order Gonze–Lee.

**B11. 2D crystals (E Sec. 8, Eqs. 455–472) [E9].** A sheet of point dipoles at $z=0$ couples to $D^{\rm e}(\mathbf Q,\omega;0,0)$. In the electrostatic limit $\tilde\omega^2=\omega_\sigma^2+\frac{1}{S_{\rm uc}}\frac{v_Q}{\epsilon^{\rm e}_{2D}(\mathbf Q,0)}|\mathbf Q\cdot\mathbf j_\sigma|^2$ (E Eq. 472). This gives no LO–TO splitting at $\Gamma$, consistent with Sohier et al., Nano Lett. 17, 3758 (2017).

**Takeaway for Part B.** One Lagrangian, one propagator ($D$), one self-energy ($\Pi_{jj}$). Coulomb is the longitudinal photon. LO–TO splitting, Fröhlich coupling and $\epsilon^\infty$ are outputs, not inputs. Retardation now also means finite $c$.

---

## 3. Bridge: electron–phonon → electron–phonon–photon

### 3.1 Concept map

| Concept | Giustino (sec, eq) | EPP (sec, eq) | What changes / what is new |
|---|---|---|---|
| Degrees of freedom / Hamiltonian | Sec. I, Eq. 1; many-body $H$, Sec. IV.A, Eqs. 62–69 [G1, G3] | Sec. 2, Eq. 109 (bare); Sec. 4, Eq. 134 (Lagrangian), Eqs. 236–237 [E3–E5] | Adds transverse photons. No explicit $g$ in $V'$: e–ph coupling is field-mediated (Coulomb $\rho_e\rho_{\rm ion}$ cross term plus photon exchange $\mathbf j_e\cdot\mathbf A\leftrightarrow Z^*\pi\cdot\mathbf A$). No second-order (DW) vertex. |
| Lattice dynamics | Sec. III.A, Eqs. 12–22, App. B [G1, G2] | Sec. 1, Eqs. 1–91 [E1, E2] | Same harmonic theory. EPP's "bare" phonons must come from **short-range (analytic)** force constants; the long-range dipole–dipole part is generated by $D$ (Sec. 7). |
| Bare vertex | Eqs. 57–58, 142 ($g^{\rm b}=\Delta_{\mathbf q\nu}V^{\rm en}$) [G2, G4] | Sec. 2, Eq. 107 ($g^0$, rigid ion) [E3] | Same object. EPP puts $N^{-1/2}$ inside $g$. |
| Screened vertex | Eq. 38 (DFPT, static RPA+xc); Eqs. 53–56; Eq. 143 (dynamic $\epsilon_{\rm e}^{-1}g^{\rm b}$) [G2, G4] | Sec. 7, $D^{\rm e}$ (Eq. 388), $\hat\epsilon^{\rm e}$ (Eqs. 431–432) [E9] | Screening becomes a dynamic RPA **tensor** acting on photon lines. No xc, no local fields. |
| Coulomb interaction | instantaneous $v(\mathbf r,\mathbf r')\delta(t-t')$, Eqs. 64, 74 [G3] | longitudinal photon: Eq. 235; $D^{0l}=-\frac{4\pi c}{q^2}\frac{q_\mu q_\nu}{\omega^2}$, Eq. 374 [E5, E8] | Derived, not postulated: Gauss's law turns $\frac1{8\pi}\int E_l^2$ into Coulomb. |
| Transverse photons | absent; classical field only in Eq. 200 [G6] | Eq. 214 ($\hbar cq$), Eq. 261 ($D^{0t}$) [E5, E6] | New, quantized, finite $c$. |
| How matter couples | scalar potential, density vertices $\rho$ (Eqs. 74–80) [G3] | $-\frac1c\mathbf j\cdot\mathbf A$ plus $A^2$ terms, Eqs. 186, 237, 273 [E5, E6] | Current vertices. Diamagnetic terms are needed for gauge invariance. |
| Polarization / boson self-energy | $P=P_{\rm e}+P_{\rm n}$ (Eq. 89); $P_{\rm e}=GG\Gamma$ (Eq. 92); Lindhard (Eqs. 4, 59) [G1–G3] | $\Pi_{jj}=\Pi^{\rm e}+\Pi^{\rm ion}$ (Eqs. 278, 382); $\Pi^{\rm e}$ (Eqs. 282/383); $\Pi^{\rm ion}$ (Eqs. 297, 311, 385) [E6–E8] | Current–current correlators, gauge-invariant as $\Pi(\omega)-\Pi(0)$. Ions enter as a **polarization current**, not as a moving density. |
| Density ↔ current | — | connection formulas, Eqs. 341, 344, 347–350; vertex form, Eq. 428 [E7–E9] | Ward identities. $\Pi_{\rho\rho}=-\frac{c}{\omega^2}\mathbf q\cdot\Pi_{jj}\cdot\mathbf q$ maps Giustino's density language onto EPP's longitudinal sector. |
| Screened interaction | $W=v+vPW$ (Eq. 90); $W=W_{\rm e}+W_{\rm ph}$ (Eqs. 104–105, 110) [G3] | $D=D^0+D^0\Pi_{jj}D$ (Eqs. 379–380); $D^{\rm e}$ (Eq. 388); $D^{\rm e}_{11}\equiv W^{\rm e}$ (Eq. 437) [E8, E9] | The scalar $W$ becomes a $3\times3$ tensor $D$. $W$ is its $ll$ sector: $D^{ll}=-\frac{c}{\omega^2}q_\mu q_\nu W$ (Eq. 441). $W_{\rm e}DW_{\rm e}$ ↔ $\Pi^{\rm ion}$ dressed by $D^{\rm e}$ (Eq. 387). |
| Phonon propagator | $D^{\rm A}_{\mathbf q\nu}=2\omega/(\omega^2-\omega^2_{\mathbf q\nu})$ (Eqs. 132, 134) [G4] | $F_\sigma=\frac1\hbar\frac{2\omega_\sigma}{\omega^2-\omega_\sigma^2+i\eta}$ (Eqs. 256, 390) [E6, E9] | Identical up to $1/\hbar$; EPP's $Q=a+a^\dagger$ is dimensionless. |
| Phonon Dyson equation | $\mathbf D=\mathbf D^{\rm A}+\mathbf D^{\rm A}\Pi^{\rm NA}\mathbf D$ (Eq. 135); $\tilde\Omega^2=\omega^2+2\omega\Pi^{\rm NA}$ (Eq. 137) [G4] | $F=F^0+F^0\Omega^{\rm e}F$ (Eq. 398); poles $\omega^2=\omega_\sigma^2+2\omega_\sigma\Omega^{\rm e}/\hbar$ (Eq. 405) [E9] | Same form, with $\Omega^{\rm e}=\hbar\Pi$. |
| Phonon self-energy | $\Pi=ZZ'\partial\partial'[W_{\rm e}(\omega)-W_{\rm e}(0)]$ (Eqs. 121, 140); $g^{\rm b}P_{\rm e}g$ (Eq. 141); Eq. 145 [G3, G4] | $\Omega^{\rm e}=-J^*D^{\rm e}J$ (Eq. 399); $\Omega^0$ (Eq. 401); $\Delta\Omega^{\rm e}=-J^*D^{\rm e}\Pi^{\rm e,0}D^0J$ (Eq. 412) = $\sum g^{\rm dyn}\chi^0g^{\rm dyn,0}$ (p. 92) [E9] | Photon exchange (L + T) replaces Coulomb exchange. Born charges $Z^*$ replace bare nuclear $Z_\kappa$. Long-range channel only. Vertices are $\omega$-dependent. Both theories keep **one bare and one screened** vertex. |
| Adiabatic split / double counting | $\Pi=\Pi^{\rm A}+\Pi^{\rm NA}$ (Eq. 125); $\chi^0(\omega)-\chi^0(0)$ in Eq. 145 [G4] | short-range DFT force constants only (p. 21); static $W^{\rm e}(\mathbf q,0)$ in Eq. 454 [E3, E9] | Giustino subtracts the static bubble. EPP subtracts the long-range part from the force constants and regenerates it. |
| Fröhlich / long-range vertex | Eq. 7; Eq. 180 (VG) [G1, G5] | Eq. 131 (static); Eq. 415 → p. 95 (bare, $c\to\infty$); Eq. 416 → p. 98 (screened, $c\to\infty$) [E3, E9] | $\epsilon^\infty\to\hat\epsilon^{\rm e}(\mathbf q,\omega)$, nonadiabatic. TO phonons couple near the light cone. The LO kernel is $c$-independent. |
| Born charges, $\epsilon^\infty$ | inputs of Eq. 180 (DFPT) [G5] | Eqs. 110–111, 135; $\Pi^{\rm ion,dia}$ (Eq. 304); $\epsilon^\infty$ from $\Pi^{\rm e}$ (Eqs. 383, 431–432) [E3, E7, E9] | EPP's $Z^*$ has charge dimension. $\epsilon^\infty$ is generated (no longer put in by hand). |
| LO–TO splitting | implicit: $\Pi^{\rm A}$ plus Gonze–Lee term (cited in Sec. VI.A.3) [G5] | classical (Eqs. 159–165); retarded $q=0$ (Eq. 410); electrostatic (Eq. 454) [E4, E9] | Splitting becomes a self-energy from photon exchange. Order of limits matters ($q\to0$ vs $c\to\infty$). Phonon polaritons emerge. |
| Dielectric function | $\epsilon_{\rm e}=1-vP_{\rm e}$ (Eq. 95); $\epsilon^{-1}$ (Eq. 80); LST via Botti–Marques (p. 43) [G3, G6] | $\epsilon(\omega)$ (Eqs. 158–159); $\epsilon=1+4\pi c\Pi/\omega^2$; $\hat\epsilon^{\rm e}$ (Eq. 432) [E4, E9] | Electronic and ionic responses add inside one $\Pi$. |
| Electron self-energy | $\Sigma^{\rm FM}$, two screened vertices (Eqs. 153, 156–158); DW (Eqs. 154, 159–161) [G4] | not derived in EPP | Missing piece: electron self-energy from photon + phonon exchange. |
| Optics | Eq. 200 (HBB, classical field, velocity vertices) [G6] | quantized $\mathbf A$; $\Pi^{\rm ion}$ gives IR response (Eq. 311 → Gonze–Lee IR $\epsilon$ [derived]) [E7] | Light is part of the system. Phonon-assisted absorption is not yet worked out in EPP. |
| 2D materials | not treated | Sec. 8, Eqs. 455–472 [E9] | Sheet propagator and Rytova–Keldysh $W$. No LO–TO splitting at $\Gamma$. |
| Electrostatic limit | Giustino *is* this limit | Eqs. 433–454 (bulk), 460–472 (2D) [E9] | $c\to\infty$: only the $ll$ sector survives, and Giustino-type formulas are recovered. |
| Temperature | Matsubara; $n_B$, $f$ in Eqs. 158, 161, 170 [G4] | "$n_B=0$" rule (Sec. 4.5) [E6] | Valid for harmonic photon/phonon self-energies, not for electron self-energies. |
| Second-order terms | DW: Eqs. 39–40, 154, 159–161 [G2, G4] | $A^2$ diamagnetic terms: Eqs. 186, 271, 304 [E5, E6] | Analogy only. Both are "contact" terms needed for sum rules. |
| Sum rules | acoustic sum rule (Eq. 120) [G3] | f-sum / oscillator-strength sum rule (Eqs. 309–311) [E7] | Translational invariance vs gauge invariance, both enforced by subtracting static parts. |

### 3.2 Explanation: from electron–phonon to electron–phonon–photon

**What Giustino assumes.** In Giustino's field theory every interaction is the instantaneous longitudinal Coulomb interaction $v(12)=v(\mathbf r_1-\mathbf r_2)\delta(t_1-t_2)$. The only "retardation" is the frequency dependence of electronic screening, $W_{\rm e}=\epsilon_{\rm e}^{-1}(\omega)v$. Moving nuclei add a second channel, $W=W_{\rm e}+W_{\rm ph}$ with $W_{\rm ph}=W_{\rm e}DW_{\rm e}$ (G Eqs. 104–105). Its vertex is the screened gradient of the nuclear potential, $\epsilon_{\rm e}^{-1}\nabla V_\kappa$ (G Eq. 110). Phonons are dressed by the same object, $\Pi=ZZ'\partial\partial'[W_{\rm e}(\omega)-W_{\rm e}(0)]$ (G Eqs. 121, 140). In practice this is G Eq. 145, with one bare and one screened vertex. In polar crystals the dipolar part of the vertex is Fröhlich/Verdi–Giustino, $\propto\mathbf q\cdot\mathbf Z^*\cdot\mathbf e/(\mathbf q\cdot\boldsymbol\epsilon^\infty\cdot\mathbf q)$ (G Eqs. 7, 180). The same kernel gives the nonanalytic LO–TO term. Light appears only as a classical perturbation (G Eq. 200).

**What EPP replaces.** EPP replaces "Coulomb plus phonons" by "the full electromagnetic field plus bare short-range phonons". Ions couple through their Born dipoles, $\frac1c\dot{\mathbf p}\cdot\mathbf A$ (E Eq. 134), which gives the ionic current vertex $j_{\alpha\sigma}=\sum Z^*e/\sqrt m$ (E Eq. 295). In the temporal gauge the longitudinal field reproduces the instantaneous Coulomb interaction exactly (E Eq. 235). The free propagator is therefore $D^0=D^{0l}+D^{0t}$ (E Eq. 375):
- $D^{0l}=-\frac{4\pi c}{q^2}\frac{q_\mu q_\nu}{\omega^2}$ is Coulomb in vector form;
- $D^{0t}$ is the retarded photon.

All matter enters through one current–current self-energy, $\Pi_{jj}=\Pi^{\rm e}+\Pi^{\rm ion}$, made gauge-invariant as $\Pi^{\rm para}(\omega)-\Pi^{\rm para}(0)$ (E Eqs. 281, 311). Charge conservation (E Eqs. 341–350) maps density vertices onto current vertices. Giustino's objects are therefore the longitudinal projections of EPP's: $W^{\rm e}=D^{\rm e}_{11}$ and $\Pi_{\rho\rho}=-\frac c{\omega^2}\mathbf q\cdot\Pi_{jj}\cdot\mathbf q$.

**Phonons dressed by photons.** A phonon is dressed by emitting and reabsorbing an electron-screened photon: $\Omega^{\rm e}=-J^*D^{\rm e}J$ and $F=F^0+F^0\Omega^{\rm e}F$ (E Eqs. 398–399). This is G Eq. 136 with $\Omega^{\rm e}=\hbar\Pi$. There are two limits, and the order of limits matters:
- Without electrons, at $q=0$: $\omega^2=\omega_\sigma^2+4\pi|j|^2/\Omega_{\rm uc}$ (E Eq. 410). This is the retarded, polariton limit, with an unscreened ionic plasma frequency.
- Taking $c\to\infty$ first: $D^{\rm e}\to-\frac{4\pi c}{\omega^2}\frac{\mathbf q\mathbf q}{\mathbf q\cdot\hat\epsilon^{\rm e}\cdot\mathbf q}$ (E Eq. 446). The shift becomes $W^{\rm e}(\mathbf q,0)|\mathbf q\cdot\mathbf j|^2/\Omega_{\rm uc}$, the Gonze–Lee term (E Eq. 454).

**Frequency-dependent vertices.** The e–ph vertices become $g^{\rm dyn,0}=iJ^{{\rm e}*}D^0J$ and $g^{\rm dyn}=iJ^*D^{\rm e}J^{\rm e}$ (E Eqs. 415–416). For $c\to\infty$ they reduce to Verdi–Giustino with $\epsilon^\infty\to\hat\epsilon^{\rm e}(\mathbf q,\omega)$ (E pp. 95, 98). The self-energy $\Delta\Omega^{\rm e}=\sum g^{\rm dyn}\chi^0g^{\rm dyn,0}$ keeps Giustino's bare–screened structure.

**What is genuinely new:**
- TO phonons couple through transverse photons for $q\lesssim\sqrt{\epsilon_\infty}\,\omega/c\sim10^{-4}$ Å$^{-1}$.
- Polaritons, LO-phonon–plasmon modes and $\epsilon^\infty$ all come out of one Dyson equation.
- In 2D, $W^{\rm 2D}=v_Q/\epsilon_{2D}$, so there is no splitting at $\Gamma$ (E Eq. 472).

**Still missing in EPP:** short-range deformation-potential coupling, electron self-energies, and xc and local fields.

---

## 4. Notation dictionary

| Quantity | Giustino | EPP | Notes |
|---|---|---|---|
| Cell index | $p$ | $i$ | |
| Atom in cell | $\kappa$ | $a$ (and $b$) | |
| Cartesian index | $\alpha,\alpha'$ | $\mu,\nu$ (Secs. 1–4); $\alpha,\beta$ also Cartesian in Secs. 5–7 | EPP reuses $\alpha,\beta$ for branches in Eq. 287. |
| Branch | $\nu$ | $\alpha$ (Secs. 1–6), $\sigma$ (from p. 61 and Sec. 7) | $\sigma$ also labels photon polarization (Eq. 219) and $l/t$ (Eq. 324). |
| Atomic mass | $M_\kappa$; reference $M_0$ | $m_a$; no reference mass ($m_0$ in Eq. 271 cancels) | |
| Number of cells | $N_p$ | $N$ | |
| Cell volume | $\Omega$ | $\Omega_{\rm uc}$ (also $\Omega_0$) | EPP's $\Omega$ is the crystal volume $=N\Omega_{\rm uc}$. EPP also uses $\Omega^2$ for the dynamical matrix and $\Omega^{\rm e}$ for the phonon self-energy. |
| Displacement | $\Delta\tau_{\kappa\alpha p}$ | $Q_{ia\mu}$ | |
| Force constants | $C_{\kappa\alpha p,\kappa'\alpha'p'}$ | $k_{ia\mu,jb\nu}$ | |
| Dynamical matrix | $D^{\rm dm}(\mathbf q)$, G Eq. 14 | $\tilde\Omega^2(\mathbf q)$, E Eq. 15 | E Eq. 15 = $D^{\rm dm}(-\mathbf q)$ (phase convention; Sec. 6). |
| Eigenvector | $e_{\kappa\alpha,\nu}(\mathbf q)$ | $e_{a\mu\alpha}(\mathbf q)$ | Same normalization. |
| Band energy / Bloch state | $\varepsilon_{n\mathbf k}$, $\psi_{n\mathbf k}=e^{i\mathbf k\mathbf r}u_{n\mathbf k}$ (normalized on supercell, $u$ on cell) | $E_{n\mathbf k}$, $\phi_{n\mathbf k}=e^{i\mathbf k\mathbf r}u/\sqrt\Omega$ | |
| e–ph term | $N_p^{-1/2}g_{mn\nu}(\mathbf k,\mathbf q)c^\dagger_{m\mathbf k+\mathbf q}c_{n\mathbf k}(a_{\mathbf q\nu}+a^\dagger_{-\mathbf q\nu})$ | $g_{nm\alpha}(\mathbf k,\mathbf q)c^\dagger_{n\mathbf k}c_{m\mathbf k-\mathbf q}(\ldots)$, $N^{-1/2}$ inside $g$ | $g^{\rm EPP}_{nm}(\mathbf k,\mathbf q)=N^{-1/2}g^{\rm G}_{nm}(\mathbf k-\mathbf q,\mathbf q)$ (labels shifted). |
| Born charge | $\mathbf Z^*_\kappa$ dimensionless (units of $e$) | $Z^*_{a\mu\nu}$ with dimension of charge | $Z^*_{\rm EPP}=|e|Z^*_{\rm G}$. |
| Coulomb kernel | $v^{\rm C}(\mathbf q)=e^2/(\varepsilon_0\Omega q^2)$ (SI, per cell) | $v(q)=4\pi/q^2$ (Gaussian, charges in the vertices) | $e^2/4\pi\varepsilon_0=14.3996$ eV Å. |
| Units | SI throughout | Sec. 3 SI; Sec. 4 onward Gaussian; Eq. 131 hybrid | Gaussian→SI: $e^2\to e^2/4\pi\varepsilon_0$, $4\pi\to1/\varepsilon_0$ in Poisson; $\omega_p^2=4\pi Z^{*2}/(\Omega M)\to Z^{*2}e^2/(\varepsilon_0\Omega M)$. |
| Electron charge | $-e$ explicit | $-e$ in Sec. 3, signed $e$ in Sec. 4+ | |
| Spatial Fourier | crystal: $f(\mathbf G)=\Omega^{-1}\int_{\rm uc}e^{-i\mathbf G\mathbf r}f$ | $f(\mathbf q)=\int d\mathbf r\,e^{-i\mathbf q\mathbf r}f(\mathbf r)$; field $\mathbf A=\Omega^{-1/2}\sum_{\mathbf q}\mathbf A_{\mathbf q}e^{i\mathbf q\mathbf r}$; crystal kernels via E Eq. 294 | Lattice phase $e^{+i\mathbf q\cdot\mathbf R}$ in both. |
| Time Fourier | $f(\omega)=\int dt\,e^{i\omega t}f(t)$ | same (E Eq. 255) | Time dependence $e^{-i\omega t}$ in both. |
| Phonon propagator | $D=-\frac i\hbar\langle T\Delta\tau\Delta\tau\rangle$ (G Eq. 109); mode form $D^{\rm A}_{\mathbf q\nu}$ in 1/frequency | $F=-\frac i\hbar\langle TQQ\rangle$, $Q$ dimensionless (E Eq. 251), 1/energy | $D^{\rm A}_{\mathbf q\nu}=\hbar F_\sigma$. |
| Phonon self-energy | $\Pi$ (frequency units in mode basis; force-constant units in Cartesian basis) | $\Omega^{\rm e}_{\sigma\sigma'}$ (energy) | $\Pi^{\rm G}=\Omega^{\rm e}/\hbar$. |
| "$D$" | nuclear density correlator (G Eq. 103) and phonon GF (G Eq. 109) | photon propagator, $D=+\frac{i}{\hbar c}\langle TAA\rangle$ (E Eqs. 257, 275) | Name clash: rename on the page (see Sec. 6). |
| "$\Pi$" | phonon self-energy | photon self-energy $\Pi_{jj}=+\frac{i}{\hbar c}\langle Tjj\rangle+$dia; $\Pi_{\rho\rho}=-\frac i\hbar\langle T\rho\rho\rangle$; $\Pi_{j\rho}$ with $-\frac i\hbar$; $\Pi_{\rho j}$ with $+\frac{i}{\hbar c}$ (E Eqs. 320, 345–346) | Name clash. |
| Screened interaction | $W$, $W_{\rm e}$, $W_{\rm ph}$ | $D_{11}$; $W^{\rm e}\equiv D^{\rm e}_{11}$ (E Eq. 437) | |
| Electronic polarizability | $P_{\rm e}$ (G Eq. 92), $\chi^0$ (G Eq. 50) | $\Pi^{\rm e}_{\rho\rho}$; $\chi^0_{n\mathbf k,m\mathbf k+\mathbf q}$ (E p. 92, not plain Lindhard) | |
| Field relations | — | $\mathbf E=\frac{i\omega}{c}\mathbf A$ ($\Phi=0$); $\nabla\cdot\mathbf j=i\omega\rho$ (E Eqs. 334–335) | $\epsilon=1+4\pi c\,\Pi/\omega^2$; $\sigma=-\frac{ic}{\omega}\Pi$ [derived from E Eqs. 279, 334]. |
| Mode oscillator strength | $\mathbf Z^*$, $\mathbf e$ in Eq. 180 | $j_{\alpha\sigma}j^*_{\beta\sigma}$ (E Eq. 295) | Equals the Gonze–Lee $S_{\sigma,\alpha\beta}$ (up to $e^2$). |

---

## 5. Web page outline (18 sections)

Global page conventions:
- State units per panel: meV, eV, Å, amu, cm$^{-1}$.
- Use the Giustino convention without $N_p$ (G Eq. 38) for $g$ in Part A, and Gaussian units with charge-dimensioned $Z^*$ in Part B. Convert to SI numbers via $e^2/4\pi\varepsilon_0=14.3996$ eV Å.
- Constants for JS: $\hbar=6.582120\times10^{-16}$ eV s; $\hbar c=1973.27$ eV Å; $\hbar^2/2m_e=3.80998$ eV Å$^2$; $k_B=0.0861733$ meV/K; 1 eV = 8065.54 cm$^{-1}$; $\hbar^2/({\rm amu}\cdot{\rm Å}^2)=4.1803$ meV.
- Rename on the page: Giustino's phonon GF $\mathcal D$, phonon self-energy $\Pi^{\rm ph}$; EPP's photon propagator $D$, photon self-energy $\Pi_{jj}$, phonon self-energy $\Omega^{\rm e}$.
- Every viz uses retarded denominators $\omega^2-\omega_0^2+i\gamma\omega$.

### Part A: Electron–phonon physics (Giustino 2017)

#### §1. Lattice dynamics and phonons
- **Learning goals:**
  - Go from force constants to the dynamical matrix to normal modes.
  - Distinguish acoustic and optical branches.
  - Quantize the modes and understand zero-point amplitude and thermal displacement.
  - Recognize that Giustino and EPP use the same theory.
- **Key equations:**
  - $D^{\rm dm}_{\kappa\alpha,\kappa'\alpha'}(\mathbf q)=(M_\kappa M_{\kappa'})^{-1/2}\sum_pC_{\kappa\alpha0,\kappa'\alpha'p}e^{i\mathbf q\cdot\mathbf R_p}$ (G Eq. 14; = E Eq. 15 up to $\mathbf q\to-\mathbf q$).
  - $\sum_{\kappa'\alpha'}D^{\rm dm}_{\kappa\alpha,\kappa'\alpha'}e_{\kappa'\alpha',\nu}=\omega^2_{\mathbf q\nu}e_{\kappa\alpha,\nu}$ (G Eq. 16).
  - $\Delta\tau_{\kappa\alpha p}=(M_0/N_pM_\kappa)^{1/2}\sum_{\mathbf q\nu}e^{i\mathbf q\cdot\mathbf R_p}e_{\kappa\alpha,\nu}(\mathbf q)\,l_{\mathbf q\nu}(\hat a_{\mathbf q\nu}+\hat a^\dagger_{-\mathbf q\nu})$ (G Eq. 20).
  - $l_{\mathbf q\nu}=[\hbar/(2M_0\omega_{\mathbf q\nu})]^{1/2}$ (G Eq. 21).
  - $\hat H_p=\sum\hbar\omega_{\mathbf q\nu}(\hat a^\dagger\hat a+\tfrac12)$ (G Eq. 22; = E Eq. 72/78).
  - $\langle u^2\rangle=\frac{\hbar}{2M\omega}\coth\frac{\hbar\omega}{2k_BT}$ (E Eq. 83 [derived]).
- **Interactive viz "Diatomic chain lab":**
  - Plot (a) $\hbar\omega(q)$ for $qa\in[-\pi,\pi]$, two branches. Plot (b) an animated chain, 12 cells, $u_\kappa(R,t)=A(M_0/M_\kappa)^{1/2}{\rm Re}[e_\kappa e^{i(qR-\omega t)}]$. Plot (c) rms displacement vs $T$.
  - Controls: $m_1,m_2$ (1–200 amu); springs $K_1,K_2$ (0.1–30 eV/Å$^2$); clickable $q$ cursor; branch selector; $T$ (0–1000 K).
  - JS: $A=(K_1+K_2)/m_1$, $C=(K_1+K_2)/m_2$, $D_{12}=-(K_1+K_2e^{-iqa})/\sqrt{m_1m_2}$, $|B|^2=(K_1^2+K_2^2+2K_1K_2\cos qa)/(m_1m_2)$, $\omega^2_\pm=\frac{A+C}2\pm\sqrt{(\frac{A-C}2)^2+|B|^2}$.
  - Unit conversion: $\hbar\omega\,[{\rm meV}]=64.654\sqrt{\omega^2[{\rm eV\,Å^{-2}amu^{-1}}]}$. Eigenvector $(B^*\!,\ \lambda-A)$ normalized, component 1 real positive (Giustino phase). $l^2[{\rm Å^2}]=2.0901/(m[{\rm amu}]\,\hbar\omega[{\rm meV}])$; rms uses $\coth(\hbar\omega/2k_BT)$.
  - Presets:
    - GaAs-like: 69.72/74.92 amu, $K_1=K_2=4.79$ eV/Å$^2$, giving optical $\Gamma$ = 33.3 meV.
    - hBN-like: 10.81/14.01 amu, $K=21$ eV/Å$^2$, giving 169.6 meV.
    - EPP worked example (E Eq. 30): $m=1$, $K_1=0.5$, $K_2=2.5$, $qa=\pi$, giving $[[3,2],[2,3]]$, $\omega^2=1,5$, $\mathbf e=(1,\mp1)/\sqrt2$.
  - What to notice: the gap at the zone edge when $m_1\neq m_2$; in-phase vs out-of-phase motion at $\Gamma$; heavier atom still in the optical mode at $q=\pi/a$; a zero-point plateau below $T\sim\hbar\omega/3k_B$; light atoms vibrate more.
- **Static figure:** chain schematic with masses, springs $K_1,K_2$, cell length $a$; a 3D BZ sketch with $\mathbf q$.
- **Sources:** [G1] III.A, [G2] V1, [E1] Eq. 30, [E2] V3–V4.

#### §2. The electron–phonon Hamiltonian and a century of "g"
- **Learning goals:**
  - Read Eq. 1 term by term, including what is *not* in it (e–e, light).
  - $g$ = zero-point amplitude × gradient of an effective potential.
  - Screening turns Nordheim's $1/q^2$ catastrophe into Bardeen's $-\tfrac23\varepsilon_F$.
  - Know the deformation-potential and Fröhlich special cases.
- **Key equations:**
  - G Eq. 1 (as in Part A, A1).
  - $g_{mn\nu}(\mathbf k,\mathbf q)=-i\big(\frac{\hbar}{2M_\kappa\omega_{\mathbf q\nu}}\big)^{1/2}\mathbf q\cdot\mathbf e_{\kappa\nu}(\mathbf q)V_0$ (G Eq. 2; $N_p$ dropped).
  - $V_0\to V_\kappa(\mathbf q)/\epsilon(q)$, $\epsilon(q)=1+(k_{TF}/q)^2F(q/2k_F)$, $F(x)=\frac12+\frac{1-x^2}{4x}\ln\left|\frac{1+x}{1-x}\right|$ (G Eqs. 3–4).
  - $V_0\to V_\kappa/\epsilon(q,\omega_{\mathbf q\nu})$ (G Eq. 5).
  - $V_0\to E_1=\Omega\,\partial\varepsilon/\partial\Omega$ (G Eq. 6).
  - $V_0\to-\big[\frac{e^2M\omega^2}{\varepsilon_0\Omega}(\frac1{\epsilon^\infty}-\frac1{\epsilon^0})\big]^{1/2}\frac1{q^2}$ (G Eq. 7).
  - $g_{mn\nu}=\langle u_{m\mathbf k+\mathbf q}|\Delta_{\mathbf q\nu}v^{\rm KS}|u_{n\mathbf k}\rangle_{\rm uc}$ (G Eq. 38).
- **Interactive viz "From Bloch to Bardeen":**
  - Plot (a) $V_0(q)$ in eV vs $q/k_F\in[0.02,3]$ for five models. Plot (b) $|g(q)|$ in meV on a log axis.
  - Models:
    - Bloch: $-1.7007$ eV.
    - Nordheim: $V_N=-4\pi Z\cdot14.3996/(\Omega q^2)$.
    - Mott–Jones: $V_N[1-\cos(qR_{\rm WS})]$.
    - Thomas–Fermi: $V_N/(1+k_{TF}^2/q^2)$.
    - Bardeen–Lindhard: $V_N/\epsilon(q)$.
  - Coupling: $|g|=\sqrt{2.0901/(M\hbar\omega[{\rm meV}])}\;q\,|V|$ with $\hbar\omega[{\rm meV}]=6.582\times10^{-3}v_s[{\rm m/s}]\,q[{\rm Å^{-1}}]$.
  - Jellium: $k_Fa_0=1.91916/r_s$, $\varepsilon_F=(k_Fa_0)^2\cdot13.6057$ eV, $(k_{TF}a_0)^2=4k_Fa_0/\pi$, $\Omega=Z\frac{4\pi}3(r_sa_0)^3$, $R_{\rm WS}=Z^{1/3}r_sa_0$, $a_0=0.529177$ Å.
  - Controls: $r_s$ (1.5–6), $Z$ (1–4), $M$ (amu), $v_s$ (1000–8000 m/s), model checkboxes.
  - Presets:
    - Al: $r_s=2.07$, $Z=3$, 26.98 amu, 6400 m/s.
    - Na: 3.93, 1, 22.99, 3200.
    - K: 4.86, 1, 39.10, 2000 [lit].
  - Marker $-\frac23\varepsilon_F$: Al $-7.80$, Na $-2.16$, K $-1.41$ eV (Cu $-4.69$).
  - What to notice: Nordheim $|g|\sim q^{-3/2}$ vs screened $|g|\sim q^{1/2}$; every screened model tends to $-\tfrac23\varepsilon_F$; a kink in the Lindhard curve at $2k_F$ (seed of Kohn anomalies); Mott–Jones is "screening by truncation".
- **Static figure:** a timeline (Bloch 1929 → Nordheim 1931 → Mott–Jones 1936 → Bardeen 1937 → Fröhlich 1950 / Bardeen–Shockley 1950 → Bardeen–Pines 1955 → DFPT 1985–87 → Wannier 2007 → Verdi–Giustino 2015). Feynman vertex: straight electron lines $n\mathbf k\to m\mathbf k+\mathbf q$, wavy phonon $\mathbf q\nu$, dot = $g$.
- **Sources:** [G1] Eqs. 1–7, V1–V2; [G2] Eqs. 38, 57–61.

#### §3. DFPT: the vertex is screened by static RPA(+xc)
- **Learning goals:**
  - Linear response as a self-consistent screening loop.
  - $\partial v^{\rm KS}=(\epsilon^{{\rm H}xc})^{-1}\partial v^{\rm en}$.
  - Test-charge vs test-electron dielectric functions.
  - DFPT is adiabatic (static $\chi^0$). Make the link to the reader's RPA $\chi^0$.
- **Key equations:**
  - G Eqs. 33–35 ($\Delta_{\mathbf q\nu}v^{\rm KS}$).
  - $(\hat H^{\rm KS}_{\mathbf k+\mathbf q}-\varepsilon_{v\mathbf k})\partial u_{v\mathbf k,\mathbf q}=-\partial_{\kappa\alpha,\mathbf q}v^{\rm KS}u_{v\mathbf k}$ (G Eq. 47, with projector form).
  - $\partial n=\chi^0_{\mathbf q}\partial v^{\rm KS}$ (G Eq. 50).
  - $\partial v^{\rm KS}=(\epsilon^{{\rm H}xc})^{-1}\partial v^{\rm en}$ (G Eq. 53).
  - $\epsilon^{{\rm H}xc}=1-(v^{\rm C}+f^{xc})\chi^0$ (G Eq. 54); $\epsilon^{\rm H}=1-v^{\rm C}\chi^0$ (G Eq. 55).
  - Jellium: $\chi^0_{\mathbf G\mathbf G'}(\mathbf q)=-\frac{m_ek_F}{\pi^2\hbar^2\Omega}F(|\mathbf q+\mathbf G|/2k_F)\delta_{\mathbf G\mathbf G'}$ (G Eq. 59).
- **Interactive viz "The SCF loop as a geometric series":**
  - Plot (a) $|\delta v_n/\delta v_{\rm en}|$ vs iteration $n$ (0–40, log axis) with the exact fixed point $1/\epsilon$. Plot (b) the contraction factor $r(q)$ vs $q/k_F$.
  - JS: $v\chi^0=-(k_{TF}/q)^2F(q/2k_F)$; iteration $\delta v_{n+1}=(1-\beta)\delta v_n+\beta[\delta v_{\rm en}+(v+f)\chi^0\delta v_n]$.
  - Fixed point $\delta v_{\rm en}/\epsilon$ with $\epsilon=1-(v+f)\chi^0$; $r=1-\beta\epsilon$.
  - Kerker option: $\beta(q)=\beta_0q^2/(q^2+q_0^2)$, $q_0=k_{TF}$.
  - Toy xc toggle: $(v+f_x)\chi^0=v\chi^0(1-q^2/4k_F^2)$.
  - Controls: $q/k_F$ (0.05–3), $r_s$ (1–6), $\beta$ (0.05–1.5), Kerker on/off, xc on/off.
  - What to notice: plain mixing with $\beta=1$ diverges whenever $\epsilon>2$, i.e. at small $q$ in metals ("charge sloshing"); Kerker makes $\beta\epsilon\approx\beta_0$, uniform in $q$; xc slightly weakens screening (test electron vs test charge).
- **Static figure:** flowchart "$\partial v^{\rm en}$ → Sternheimer → $\partial n$ → $\partial v^{\rm H}+\partial v^{xc}$ → $\partial v^{\rm KS}$ → loop". Diagram: dressed vertex = bare dot + bubble chain (Dyson for $g$).
- **Sources:** [G2] B.3–B.5, V3, V5.

#### §4. Polar materials: Fröhlich, Verdi–Giustino and the LO–TO splitting
- **Learning goals:**
  - Why $g\sim1/q$ in polar crystals.
  - What sets the strength: $1/\epsilon^\infty-1/\epsilon^0$, $Z^*$.
  - Why Wannier interpolation needs the analytic long-range part.
  - The same kernel gives the nonanalytic LO–TO term.
  - Where electrostatics must fail (light cone, preview of EPP).
- **Key equations:**
  - G Eq. 7; G Eq. 180 (as in Part A, A8).
  - Reduction to Fröhlich: $|g^{\mathcal L}|=\frac1q\big[\frac{e^2}{4\pi\varepsilon_0}\frac{4\pi}{\Omega}\frac{\hbar\omega_{\rm LO}}2(\frac1{\epsilon^\infty}-\frac1{\epsilon^0})\big]^{1/2}$ ([G5] [derived]).
  - Gonze–Lee: $\tilde C^{\rm NA}_{\kappa\alpha,\kappa'\beta}=\frac{4\pi}{\Omega}\frac{e^2}{4\pi\varepsilon_0}\frac{(\mathbf q\cdot\mathbf Z^*_\kappa)_\alpha(\mathbf q\cdot\mathbf Z^*_{\kappa'})_\beta}{\mathbf q\cdot\boldsymbol\epsilon^\infty\cdot\mathbf q}$ [ctx].
  - LST: $\omega_{\rm LO}^2/\omega_{\rm TO}^2=\epsilon^0/\epsilon^\infty$.
  - EPP's static version: E Eq. 131.
- **Interactive viz "Fröhlich vertex explorer":**
  - Plot (a) $|g(q)|$ in meV vs $q$ (Å$^{-1}$, log–log, $10^{-6}$ to 1) for the LO Fröhlich and acoustic deformation-potential vertices. Plot (b) a bar of the LO–TO splitting.
  - JS:
    - $\sqrt{N}|g|=\frac1q\sqrt{180.95\,\hbar\omega_{\rm LO}[{\rm eV}]/(2\Omega)\,(1/\epsilon^\infty-1/\epsilon^0)}$ eV.
    - $|g_{\rm ac}|=E_1q\sqrt{2.0901/(M_{\rm cell}\cdot6.582\times10^{-3}v_sq)}$ (meV inside the root).
    - $(\hbar\omega_{\rm LO})^2=(\hbar\omega_{\rm TO})^2+0.7564\,Z^{*2}/(\epsilon^\infty\Omega\mu)$ eV$^2$ ($\Omega$ Å$^3$, $\mu$ amu).
    - $\alpha=(1/\epsilon^\infty-1/\epsilon^0)\sqrt{13.606\,m^*/\hbar\omega_{\rm LO}[{\rm eV}]}$.
    - Shaded "light cone" $q<q_c=\sqrt{\epsilon^\infty}\,\omega_{\rm LO}/c=\sqrt{\epsilon^\infty}\,\hbar\omega_{\rm LO}/1973.27$ Å$^{-1}$.
  - Controls: material preset; sliders $\epsilon^\infty$ (1–15), $\epsilon^0$ ($\ge\epsilon^\infty$, to 20), $\hbar\omega_{\rm LO}$ (10–200 meV), $\Omega$ (10–100 Å$^3$), $Z^*$ (0–3), $E_1$ (0–20 eV).
  - Presets ($\epsilon^\infty$, $\epsilon^0$, $\hbar\omega_{\rm LO}$/$\hbar\omega_{\rm TO}$ meV, $\Omega$ Å$^3$, $m^*$, $Z^*$, $\mu$):
    - GaAs: 10.9, 12.9, 36.2/33.3, 45.2, 0.067, 2.07, 36.1 → $\sqrt N|g|q=0.032$ eV Å, $\alpha=0.071$.
    - GaN: 5.35, 9.5, 92, 45.7, 0.2 → $\alpha\approx0.44$.
    - MgO: 2.96, 9.8, 89, 18.7, 0.35, $Z^*=1.98$, $\mu=9.65$ → 0.32 eV Å.
    - LiF: 1.93, 9.0, 82, 16.4 → 0.43 eV Å.
    - NaCl: 2.34, 5.9, 32.7, 44.9, $Z^*=1.12$, $\mu=13.9$ → 0.13 eV Å.
  - What to notice: LO coupling grows as $1/q$ and dominates the acoustic coupling ($\propto q^{1/2}$) at small $q$; $\alpha$ spans weak (GaAs) to intermediate (oxides, fluorides); the LO–TO splitting uses the same $Z^*$ and $\epsilon^\infty$; the light cone ($\sim10^{-5}$–$10^{-4}$ Å$^{-1}$) is where instantaneous electrostatics stops being exact (teaser for §10).
  - Optional tab "1D Wannier interpolation lab": $g^{\mathcal S}(q)=A\sin^2(qa/2)+B\sin^2(qa)$ and $g^{\mathcal L}(q)=C\sum_ne^{-\alpha(q+nb)^2}/(q+nb)$. Fourier-interpolate from $N_c=4$–16 coarse points with and without subtracting $g^{\mathcal L}$, and show ringing vs convergence.
- **Static figure:** a lattice of point dipoles $e\mathbf Z^*\cdot\Delta\boldsymbol\tau$ modulated by $e^{i\mathbf q\cdot\mathbf R}$; the potential falls as $1/r^2$. Workflow "subtract $g^{\mathcal L}$ → Wannier-interpolate $g^{\mathcal S}$ → add $g^{\mathcal L}$".
- **Sources:** [G1] V3, [G5] VI.A.3, Viz 1–2; [E3] V3.

#### §5. Hedin–Baym: the lattice inside W ($W=W_{\rm e}+W_{\rm ph}$)
- **Learning goals:**
  - The reader's GW $W$ is $W_{\rm e}$, at clamped nuclei.
  - Moving nuclei add $W_{\rm ph}=W_{\rm e}\mathcal DW_{\rm e}$.
  - The e–ph vertex is $\epsilon_{\rm e}^{-1}\nabla V_\kappa$.
  - Phonons themselves are nuclei interacting via $W_{\rm e}$.
  - Table I closes the set; $\mathcal D=0$ recovers Hedin.
  - Every $v$ is instantaneous.
- **Key equations:**
  - $W=v+vPW$, $P=P_{\rm e}+P_{\rm n}$ (G Eqs. 89–90).
  - $W_{\rm e}=v+vP_{\rm e}W_{\rm e}$ (G Eq. 94).
  - $W=W_{\rm e}+W_{\rm ph}$, $W_{\rm ph}=W_{\rm e}\mathcal DW_{\rm e}$ (G Eqs. 104–105).
  - $W_{\rm ph}(12)=\sum\int d(34)\,\epsilon_{\rm e}^{-1}(13)\nabla_{3\alpha}V_\kappa\,\mathcal D_{\kappa\alpha p,\kappa'\alpha'p'}(t_3t_4)\,\epsilon_{\rm e}^{-1}(24)\nabla_{4\alpha'}V_{\kappa'}$ (G Eq. 110).
  - $\Pi^{\rm ph}=\sum Z_\kappa Z_{\kappa''}\partial_\alpha\partial'_{\alpha'}[\delta_{\kappa'p',\kappa''p''}W_{\rm e}(\omega)-\delta_{\kappa p,\kappa'p'}W_{\rm e}(0)]$ (G Eq. 121).
  - Long-wavelength polar model [derived, E4 Concept 4.10]: $W_{\rm ph}=v[1/\epsilon(\omega)-1/\epsilon^\infty]=|g_{\rm Fr}|^2\mathcal D_{\rm LO}$.
- **Interactive viz A "Lattice screening in a polar insulator":**
  - Plot Re and Im of $W/v$, $W_{\rm e}/v$, $W_{\rm ph}/v$ vs $\hbar\omega$ (0–$2\hbar\omega_{\rm LO}$).
  - JS: $\epsilon=\epsilon^\infty(E_{\rm LO}^2-E^2-i\gamma E)/(E_{\rm TO}^2-E^2-i\gamma E)$; $w=1/\epsilon$; $w_{\rm e}=1/\epsilon^\infty$; $w_{\rm ph}=\frac1{\epsilon^\infty}\frac{E_{\rm TO}^2-E_{\rm LO}^2}{E_{\rm LO}^2-E^2-i\gamma E}$.
  - Overlay the Fröhlich form $w_{\rm ph}^{\rm F}=\frac{g^2}{v}\frac{2E_{\rm LO}}{E^2-E_{\rm LO}^2}$ with $\frac{g^2}{v}=\frac{E_{\rm LO}}2(\frac1{\epsilon^\infty}-\frac1{\epsilon^0})$. It is identical for $\gamma=0$ via LST.
  - Toggle "$\mathcal D=0$" (Hedin limit).
  - Controls: $\gamma$ (0.1–5 meV) and the presets below.
  - Presets ($E_{\rm TO}$/$E_{\rm LO}$ meV, $\epsilon^\infty$): GaAs 33.3/36.2/10.9; MgO 49.1/89.0/2.95; hBN 169.9/199.6/4.87.
  - What to notice: $W(0)=v/\epsilon^0<W_{\rm e}=v/\epsilon^\infty$ (the lattice screens statically); $W$ changes sign between TO and LO (Reststrahlen); $W_{\rm ph}$ has a pole at LO, not TO; the Fröhlich constant is just the residue.
- **Interactive viz B "Jellium with mobile ions (Bardeen–Pines)":**
  - Plot ${\rm Re}\,W/v$ vs $\hbar\omega$.
  - JS: $W/v=\frac1{\epsilon_{\rm e}}\frac{\omega^2}{\omega^2-\omega_q^2+i\gamma\omega}$, $\omega_q^2=\Omega_p^2/\epsilon_{\rm e}(q)$, $\epsilon_{\rm e}=1+(k_{TF}/q)^2F$, $(\hbar\Omega_p)^2=0.7564\,Z^2/(\Omega M)$ eV$^2$ (Na: 29.6 meV).
  - What to notice: attraction ($W<0$) below the screened ion-plasma (phonon) frequency, the root of BCS pairing.
- **Interactive viz C "Clickable Hedin–Baym map":** an 11-node graph of Table I ($G$, $\Sigma$, $\Gamma$, $P_{\rm e}$, $\epsilon_{\rm e}$, $W_{\rm e}$, $W_{\rm ph}$, $\mathcal D$, $\Pi^{\rm ph}$, $\langle n_{\rm e}\rangle$, $\langle n_{\rm n}\rangle$, $V_{\rm tot}$). Clicking a node shows its equation. Toggle "$\mathcal D=0$" greys out the phonon half. Toggle "where is $v(12)=v\,\delta(t_1-t_2)$?" highlights every instantaneous Coulomb line; these are what EPP replaces.
- **Static figure:** diagram $W_{\rm ph}$ = (double-dashed $W_{\rm e}$)–(wavy $\mathcal D$)–(double-dashed $W_{\rm e}$), equivalently gray disk $g$ – wavy – gray disk $g$.
- **Sources:** [G3] Eqs. 88–122, Table I, Viz 1, 2, 4; [E4] 4.10.

#### §6. Phonon self-energy: adiabatic vs nonadiabatic
- **Learning goals:**
  - Split $\Pi^{\rm ph}=\Pi^{\rm A}+\Pi^{\rm NA}$; $\Pi^{\rm A}$ is the DFPT force-constant matrix.
  - The Dyson equation gives shifts, widths and new poles.
  - Eq. 145 has one bare and one screened vertex, plus the $\chi^0(\omega)-\chi^0(0)$ double-counting subtraction.
  - Landau damping and the nonadiabatic Kohn anomaly; when the adiabatic approximation fails ($\varepsilon_F\sim\hbar\omega$).
- **Key equations:**
  - $\mathcal D(\omega)=[\mathbf M\omega^2-\Pi^{\rm ph}(\omega)]^{-1}$ (G Eq. 123).
  - $\Pi^{\rm ph}=\Pi^{\rm A}+\Pi^{\rm NA}$ (G Eq. 125).
  - $\mathcal D^{\rm A}_{\mathbf q\nu}=\frac{1}{\omega-\omega_{\mathbf q\nu}+i\eta}-\frac{1}{\omega+\omega_{\mathbf q\nu}-i\eta}$ (G Eq. 134).
  - $\mathcal D^{-1}_{\mathbf q\nu\nu'}=\frac1{2\omega_{\mathbf q\nu}}[\delta_{\nu\nu'}(\omega^2-\omega^2_{\mathbf q\nu})-2\omega_{\mathbf q\nu}\Pi^{\rm NA}_{\mathbf q\nu\nu'}]$ (G Eq. 136).
  - $\tilde\Omega^2=\omega^2+2\omega\Pi^{\rm NA}(\tilde\Omega)$ (G Eq. 137).
  - $\hbar\Pi^{\rm NA}=g^{\rm b}P_{\rm e}(\omega)g^{\rm cc}(\omega)-g^{\rm b}P_{\rm e}(0)g^*(0)$ (G Eq. 141).
  - $\hbar\Pi^{\rm NA}_{\mathbf q\nu\nu'}(\omega)=2\sum_{mn}\int\frac{d\mathbf k}{\Omega_{\rm BZ}}g^{\rm b}_{mn\nu}g^*_{mn\nu'}\Big[\frac{f_{m\mathbf k+\mathbf q}-f_{n\mathbf k}}{\varepsilon_{m\mathbf k+\mathbf q}-\varepsilon_{n\mathbf k}-\hbar(\omega+i\eta)}-\frac{f_{m\mathbf k+\mathbf q}-f_{n\mathbf k}}{\varepsilon_{m\mathbf k+\mathbf q}-\varepsilon_{n\mathbf k}}\Big]$ (G Eq. 145).
  - Linewidth: G Eq. 146.
- **Interactive viz "Nonadiabatic phonon in a doped 2D semiconductor":**
  - Plots (a) heat map of the phonon spectral function $B(Q,\omega)=-\frac1\pi{\rm Im}\frac{2\omega_0}{\omega^2-\omega_0^2-2\omega_0\Pi^{\rm NA}}$ on $Q=q/k_F\in[0,3]$, $\hbar\omega\in[0,2\hbar\omega_0]$, with particle–hole edges $\varepsilon_F(Q^2\pm2Q)$; (b) $\Omega(Q)$ and $\gamma(Q)$; (c) $\tau^{\rm ph}=1/2\gamma$ in ps.
  - JS: $\tilde\chi(Q,W)=-[1+\frac{s(a)+s(b)}{2Q^2}]$, $a=W-Q^2$, $b=-W-Q^2$, $s(z)=z\sqrt{1-4Q^2/z^2}$ (principal root, complex $W=(\hbar\omega+i\eta)/\varepsilon_F$).
  - $\hbar\Pi^{\rm NA}=G[\tilde\chi(Q,W)-\tilde\chi(Q,0)]$; $\Omega=\sqrt{\omega_0^2+2\omega_0{\rm Re}\Pi}$; $\gamma=-(\omega_0/\Omega){\rm Im}\Pi$.
  - Scales: $k_F=\sqrt{2\pi n}$, $\varepsilon_F=3.80998\,k_F^2/(m^*/m_e)$ eV.
  - Controls: $n$ ($10^{12}$–$10^{14}$ cm$^{-2}$, log); $m^*$ (0.1–1); $\hbar\omega_0$ (10–200 meV, default 50); $G$ (0–15 meV, default 5); $\eta$ (0.5–5 meV); first-order vs self-consistent.
  - Default $n=10^{13}$ cm$^{-2}$, $m^*=0.5$ gives $k_F=0.0793$ Å$^{-1}$, $\varepsilon_F=47.9$ meV.
  - What to notice: above the continuum the phonon hardens by about $G$ (nonadiabatic stiffening); inside it, Landau damping $\hbar\gamma\approx G\hbar\omega_0/(\varepsilon_FQ\sqrt{4-Q^2})$; for $Q>2$, $\Pi^{\rm NA}\to0$; larger $\varepsilon_F$ shrinks the nonadiabatic window.
  - Second tab "Jellium phonons": $\omega^2(q)=\Omega_p^2/\epsilon_{\rm e}(q)$ with Lindhard $\epsilon_{\rm e}$; Bohm–Staver $v_s=\Omega_p/k_{TF}$; plot $d\omega/dq$ to reveal the Kohn anomaly at $2k_F$.
  - Optional tab: graphene G band vs doping using literature parameters [lit, flagged].
- **Static figures:** G Fig. 1 (Dyson for $\mathcal D$, self-energy bubble with one dot and one disk, dressed vertex). A "double counting" cartoon: DFT phonons already contain $\chi^0(0)$.
- **Sources:** [G3] 3.8, [G4] V.A, V1, V5; [G5] VII, Viz 3.

#### §7. Electron self-energy: Fan–Migdal, Debye–Waller, kinks and satellites
- **Learning goals:**
  - $\Sigma^{\rm ep}=\Sigma^{\rm FM}+\Sigma^{\rm DW}+\Sigma^{\rm dGW}$, an exact rewrite.
  - FM = "$GW_{\rm ph}$", with **two screened** vertices (vs one in the phonon self-energy).
  - Retarded finite-$T$ form; quasiparticle renormalization $Z\approx1/(1+\lambda)$; ARPES kinks.
  - Dyson vs cumulant for satellites (same story as GW vs GW+C).
- **Key equations:**
  - G Eqs. 152–155.
  - $\Sigma^{\rm FM}=i\hbar G\Gamma W_{\rm ph}$ (G Eq. 153).
  - $\Sigma^{\rm FM}_{nn'\mathbf k}(\omega)=\frac1\hbar\sum_{m\nu}\int\frac{d\mathbf q}{\Omega_{\rm BZ}}g^*_{mn\nu}g_{mn'\nu}\Big[\frac{1-f_{m\mathbf k+\mathbf q}+n_{\mathbf q\nu}}{\omega-\varepsilon_{m\mathbf k+\mathbf q}/\hbar-\omega_{\mathbf q\nu}+i\eta}+\frac{f_{m\mathbf k+\mathbf q}+n_{\mathbf q\nu}}{\omega-\varepsilon_{m\mathbf k+\mathbf q}/\hbar+\omega_{\mathbf q\nu}+i\eta}\Big]$ (G Eqs. 157–158).
  - $\Sigma^{\rm DW}_{nn\mathbf k}=\sum_\nu\int\frac{d\mathbf q}{\Omega_{\rm BZ}}g^{\rm DW}_{nn\nu\nu}(\mathbf k,\mathbf q,-\mathbf q)(2n_{\mathbf q\nu}+1)$ (G Eqs. 160–161).
  - $A=-\frac1\pi{\rm Im}\,G$; $Z=[1-\partial_E{\rm Re}\Sigma]^{-1}$ (G Eqs. 172–174).
- **Interactive viz "Kinks and satellites":**
  - Plots (a) heat map $A(k,E)$, $k\in[0,1.5k_F]$, $E\in[-\varepsilon_F-3\hbar\omega_0,\,2\hbar\omega_0]$, clipped at 3 eV$^{-1}$; (b) side panel Re/Im $\Sigma(E)$ and the line $E-\varepsilon_k$.
  - JS: $\Sigma(E)=\frac{\lambda\hbar\omega_0}2\Big\{\ln\frac{z-\hbar\omega_0}{z-\hbar\omega_0-\varepsilon_c}+\ln\frac{z+\hbar\omega_0+\varepsilon_F}{z+\hbar\omega_0}\Big\}-{\rm Re}\Sigma(0)$, $z=E+i\eta$, logs via $\ln|N|-\ln|D|+i[\arg N-\arg D]$.
  - $A=-\frac1\pi{\rm Im}[E+i\eta-\varepsilon_k-\Sigma]^{-1}$, $\varepsilon_k=-\varepsilon_F+3.80998\,k^2/(m^*/m_e)$ eV.
  - Overlays: Brillouin–Wigner roots of $E-\varepsilon_k-{\rm Re}\Sigma(E)$ (bisection, reject log singularities); Rayleigh–Schrödinger $\varepsilon_k+{\rm Re}\Sigma(\varepsilon_k)$; $Z$ at $k_F$.
  - Controls: $\varepsilon_F/\hbar\omega_0$ (0.5–15; presets 10 and 2), $\lambda$ (0–2, default 0.8), $\eta$ (5–40 meV, default 20); fixed $\hbar\omega_0=100$ meV, $\varepsilon_c=5$ eV, $m^*=0.1$.
  - What to notice: the kink at $-\hbar\omega_0$; velocity reduced by $1/(1+\lambda)$; broadening only below $-\hbar\omega_0$; for $\varepsilon_F=2\hbar\omega_0$ a satellite below the band bottom that BW misses.
  - Tab "Dyson vs cumulant (single level)":
    - Dyson: $\Sigma=|g|^2/(E-\varepsilon_0+\hbar\omega_0+i\eta)$, poles $x_\pm=\frac{\hbar\omega_0}2(-1\pm\sqrt{1+4S})$, weights $Z_\pm=[1+S(\hbar\omega_0)^2/(x_\pm+\hbar\omega_0)^2]^{-1}$.
    - Cumulant: $A=\sum_ne^{-S}\frac{S^n}{n!}L_\eta(E-\varepsilon_0-S\hbar\omega_0+n\hbar\omega_0)$.
    - Controls: $S$ (0–4), $\hbar\omega_0$ (20–200 meV), $\eta$.
    - What to notice: Dyson gives one misplaced satellite; the cumulant gives a Poisson ladder spaced by $\hbar\omega_0$.
- **Static figures:** G Fig. 2: FM (rainbow phonon on electron line, gray disks at both ends) and DW (tadpole loop at one time). Labels: "two screened vertices" vs §6 "one bare + one screened".
- **Sources:** [G4] V.B, V3, V4; [G5] VIII.

#### §8. Temperature-dependent band gaps, zero-point renormalization, lifetimes and the Fröhlich catastrophe
- **Learning goals:**
  - AH band renormalization (FM + DW); Varshni curve; extracting the ZPR from the high-$T$ asymptote.
  - Lifetimes from Im $\Sigma$.
  - Why the adiabatic AH diverges in polar materials, and why this is an electrostatic artifact of dropping $\hbar\omega$ (not of dropping $c$).
  - How small the light cone is.
- **Key equations:**
  - G Eq. 166 (AH, finite $T$).
  - $E_g(T)=E_g^{\rm cn}-|\Delta E_g^{\rm ZP}|[1+2n(\hbar\omega_0/k_BT)]$ (G Eq. 168).
  - $\frac1{\tau_{n\mathbf k}}=\frac{2\pi}\hbar\sum_{m\nu}\int\frac{d\mathbf q}{\Omega_{\rm BZ}}|g_{mn\nu}|^2[(1-f+n)\delta(\varepsilon_{n\mathbf k}-\hbar\omega-\varepsilon_{m\mathbf k+\mathbf q})+(f+n)\delta(\varepsilon_{n\mathbf k}+\hbar\omega-\varepsilon_{m\mathbf k+\mathbf q})]$ (G Eq. 170).
  - AH in practice and the polar divergence: G Eqs. 195–199.
- **Interactive viz A "Varshni":**
  - Plots: $E_g(T)$ with asymptote and ZPR arrow; $\Gamma(T)$ and $\tau(T)$.
  - JS: $F=\coth(\hbar\omega_0/2k_BT)$; $E_g=E_g^{\rm cn}-{\rm ZPR}\cdot F$; asymptote $E_g^{\rm cn}-{\rm ZPR}\cdot2k_BT/\hbar\omega_0$; $\Gamma=\Gamma^{\rm ZP}F$; $\tau=658.212/(2\Gamma[{\rm meV}])$ fs.
  - Controls: $E_g^{\rm cn}$ (0.5–8 eV), ZPR (0–700 meV), $\hbar\omega_0$ (10–200 meV), $\Gamma^{\rm ZP}$ (0–100 meV), $T_{\max}$ (300–1000 K), normal/inverse-Varshni toggle.
  - Presets:
    - Diamond: $E_g^{\rm cn}=7.73$ eV, ZPR 0.628 eV (GW) or 0.40 (LDA), $\hbar\omega_0=125$ meV.
    - Si: 1.23 eV, 0.06 eV, 30 meV.
    - Fig.-3 toy: 1 eV, 100 meV, $\hbar\omega_0=50$ meV (matches the printed curves; caption says 100, show the label), $\Gamma^{\rm ZP}=50$ meV → 6.6 fs.
  - What to notice: zero-point plateau, linear high-$T$ regime, intercept $=E_g^{\rm cn}$, width and shift share the Bose factor.
- **Interactive viz B "Fröhlich catastrophe":**
  - Plot $\Delta E$ (meV) vs $N_q$ (q-grid size, 4–1000, log axis).
  - JS: $\Delta E_{\rm NA}=-\alpha\hbar\omega_{\rm LO}\frac2\pi[\arctan(q_D/q_p)-\arctan(q_{\min}/q_p)]$; $\Delta E_{\rm A}=-\alpha\hbar\omega_{\rm LO}\frac2\pi q_p(1/q_{\min}-1/q_D)$; $q_p=\sqrt{(m^*/m_e)\hbar\omega_{\rm LO}[{\rm eV}]/3.80998}$ Å$^{-1}$; $q_D=(6\pi^2/\Omega)^{1/3}$; $q_{\min}=q_D/N_q$.
  - Also mark the light-cone momentum $q_c$ and the fraction $(q_c/q_D)$.
  - GaAs: $q_p=0.0252$, $q_D=1.094$ Å$^{-1}$, $q_c\approx6\times10^{-5}$ Å$^{-1}$.
  - What to notice: nonadiabatic converges to $-\alpha\hbar\omega_{\rm LO}$; adiabatic diverges linearly in $N_q$; the light cone is $\sim10^{-4}$ of the BZ, so retardation is not what cures the divergence; dynamics ($\hbar\omega$) is.
- **Interactive viz C "LO-phonon scattering rate":**
  - JS: $1/\tau=2\alpha\omega_{\rm LO}\sqrt{\hbar\omega/\varepsilon}\,[n\,{\rm asinh}\sqrt{\varepsilon/\hbar\omega}+(n+1){\rm asinh}\sqrt{\varepsilon/\hbar\omega-1}\,\Theta(\varepsilon-\hbar\omega)]$, vs $\varepsilon$ (0–5$\hbar\omega_{\rm LO}$) and $T$.
  - GaAs: $2\alpha\omega_{\rm LO}\approx7.8\times10^{12}$ s$^{-1}$ ($\alpha=0.071$).
  - What to notice: the emission threshold at $\hbar\omega_{\rm LO}$.
- **Static figure:** G Fig. 3 redrawn (gap and lifetime vs $T$).
- **Sources:** [G4] V.B.3–V.B.4, V2; [G6] IX.A, Viz 1, 2, 4.

#### §9. Where Giustino meets light (classically), transport and superconductivity
- **Learning goals:**
  - Indirect absorption needs a phonon: two second-order paths (G Eq. 200), emission vs absorption onsets.
  - The photon is a classical external field coupled only to electrons; the lattice IR response is absent.
  - Temperature-dependent spectra (Williams–Lax).
  - Brief: mobility (Boltzmann transport), $\alpha^2F\to\lambda\to T_c$.
- **Key equations:**
  - G Eq. 200 (HBB) with its two energy denominators and weights $n+\frac12\pm\frac12$.
  - G Eq. 201 (Williams–Lax $\epsilon_2(\omega;T)$).
  - $\lambda=2\int d\omega\,\alpha^2F(\omega)/\omega$ and $\omega_{\log}$ (G Eqs. 186–192).
  - Allen–Dynes: $k_BT_c=\frac{\hbar\omega_{\log}}{1.2}\exp\big[-\frac{1.04(1+\lambda)}{\lambda-\mu^*(1+0.62\lambda)}\big]$ (G Sec. XI).
- **Interactive viz A "Indirect absorption edge":**
  - Plot $\alpha(\hbar\omega)$ (cm$^{-1}$, log) for $\hbar\omega\in[E_g-0.15,\,E_g+0.3]$ eV.
  - JS: $\alpha=\frac{C}{\hbar\omega}\sum_\nu w_\nu[(n_\nu+1)(\hbar\omega-E_g-\hbar\omega_\nu)^2\Theta+n_\nu(\hbar\omega-E_g+\hbar\omega_\nu)^2\Theta]$.
  - Toggle "AH gap renormalization" with $E_g(T)$ from §8.
  - Si preset: $E_g^{\rm cn}=1.23$ eV, ZPR 0.06 eV, $\hbar\omega_{\rm eff}=30$ meV; TO 57.8 meV ($w=1$), TA 18.4 meV ($w=0.3$) [lit]; $C\approx1.2\times10^4$ cm$^{-1}$eV$^{-1}$ (illustrative).
  - Controls: $T$ (4–600 K).
  - What to notice: absorption tail below $E_g$ appears with $T$ (phonon absorption), two-slope edge, red shift.
- **Interactive viz B "$\alpha^2F$ to $T_c$":**
  - Sum of Gaussians, each specified by center (meV) and partial $\lambda$; numerical $\lambda$ and $\omega_{\log}$; Allen–Dynes $T_c$ for $\mu^*\in[0.08,0.15]$.
  - Presets: Al (20, 0.25; 37, 0.18); Pb (4.4, 0.9; 8.5, 0.65); MgB$_2$ isotropic (65, 0.6; 35, 0.1; 90, 0.03); H$_3$S (30, 0.5; 120, 1.0; 170, 0.6) [lit].
- **Static figure:** two second-order paths of G Eq. 200 as Feynman diagrams: photon (external classical, dashed arrow) – electron – phonon (wavy), in both orders. Caption: "light is external here; EPP makes it internal".
- **Sources:** [G6] IX.B, X, XI, Viz 3, 5; [G5] VIII.

### Bridge

#### §10. From instantaneous Coulomb to the full electromagnetic field
- **Learning goals:**
  - Every line of Part A is a longitudinal, instantaneous Coulomb line.
  - EPP replaces it by the photon propagator $D^0=D^{0l}+D^{0t}$.
  - "Retardation" has two meanings: $\epsilon_{\rm e}(\omega)$ (Giustino) vs finite $c$ (EPP).
  - Momentum scales show where each matters.
  - Read the concept table (Sec. 3 of this synthesis).
- **Key equations:**
  - $W=W_{\rm e}+W_{\rm e}\mathcal DW_{\rm e}$ (G Eqs. 104–105) vs $D=D^0+D^0\Pi_{jj}D$ (E Eq. 379).
  - $D^0_{\mu\nu}(\mathbf q,\omega)=\frac{4\pi c}{\omega^2}\frac{-q_\mu q_\nu+\frac{\omega^2}{c^2}\delta_{\mu\nu}}{q^2-\frac{\omega^2}{c^2}-i\eta}$ (E Eq. 375) $=-\frac{4\pi c}{\omega^2}\hat q_\mu\hat q_\nu+\frac{4\pi}{c}\frac{\delta_{\mu\nu}-\hat q_\mu\hat q_\nu}{q^2-\omega^2/c^2}$ [derived].
  - $\Pi_{\rho\rho}=-\frac c{\omega^2}\mathbf q\cdot\Pi_{jj}\cdot\mathbf q$ (E Eq. 350 in $\mathbf q$-space).
- **Interactive viz A "How big is the light cone?":**
  - Plot $|D^0_T/D^0_L|=1/|1-x^2|$ vs $q$ (log axis $10^{-7}$–1 Å$^{-1}$), $x=cq/\omega$, with $D^0_L\propto-1$ and $D^0_T\propto-1/(1-x^2+i\eta)$.
  - Controls: $\hbar\omega$ (1 meV–5 eV, log); "speed of light multiplier" $\lambda_c$ ($c\to\lambda_cc$, $10^{-2}$–$10^{6}$).
  - Markers: 50 meV phonon, hBN TO 170 meV, 2 eV exciton, BZ edge $\sim1$ Å$^{-1}$.
  - Scale ladder (static bar): $q_{\rm light}\sim10^{-5}$–$10^{-3}$ ≪ $q_p\sim0.02$–0.07 ≈ $\omega/v_F\sim0.03$ ≪ $k_{TF}\sim1$–2 ≈ $q_{\rm BZ}\sim1$ Å$^{-1}$.
  - What to notice: outside the cone the transverse part is suppressed by $(\omega/cq)^2$; as $\lambda_c\to\infty$ only $D^0_L$ survives (the Giustino world); inside the cone L and T are comparable and polaritons live there.
- **Interactive viz B "Diagram dictionary" (side by side, hover to highlight the counterpart):**
  - Giustino glyphs: straight $G$; wavy phonon $\mathcal D$; dashed $v$; double-dashed $W$; dot $g^{\rm b}$; gray disk $g$.
  - EPP glyphs: dashed $v$; red double-dashed $D^{0l}$; blue wavy $D^{0t}$; magenta $D^0$; blue double $D^{\rm e}$; zigzag phonon $F$; circle = charge vertex $\rho$; square = current vertex $j$; teal square = ionic $J$; blue square = electronic $J^{\rm e}$.
- **Static figure:** the concept table (Sec. 3.1) rendered as an interactive table with filters ("same / generalized / new").
- **Sources:** [G3] 3.9, [E8] 3.5, S1–S2, V1, V5; [E9] Concepts.

### Part B: The electron–phonon–photon (EPP) formalism

#### §11. One Lagrangian: electrons, ions and the electromagnetic field
- **Learning goals:**
  - Temporal gauge; why ions couple via $\dot{\mathbf p}\cdot\mathbf A$ (Born dipoles, not positions).
  - Lorentz force and ionic EOM from one Lagrangian.
  - Minimal coupling with para and dia terms for electrons *and* ions.
  - Field quantization: two transverse photons; the longitudinal field = Coulomb via Gauss's law.
  - No explicit $g$ in $V'$.
- **Key equations:**
  - $\mathbf E=-\frac1c\partial_t\mathbf A$, $\Phi=0$ (E Eqs. 132–133).
  - $L=\sum_e[\frac{m_e\mathbf v_e^2}2-V_{\rm cr}]+\sum\frac{m_a\dot Q^2}2-\frac12\sum kQQ+\int\frac{d\mathbf r}{8\pi}(\mathbf E^2-\mathbf B^2)+\frac ec\sum_e\mathbf v_e\cdot\mathbf A(\mathbf r_e)+\frac1c\sum_{ia}\dot{\mathbf p}_{ia}\cdot\mathbf A(\mathbf R^0_i)$ (E Eq. 134), $p_{ia\mu}=\sum_\nu Z^*_{a\mu\nu}Q_{ia\nu}$ (E Eq. 135).
  - $m_a\ddot Q_{ia\mu}+\sum k\,Q=\sum_\nu Z^*_{a\nu\mu}E_\nu(\mathbf R^0_i)$ (E Eq. 143).
  - $\pi_{ia\mu}=m_a\dot Q+\frac1c\sum_\alpha Z^*_{a\alpha\mu}A_\alpha$ (E Eq. 183).
  - $H^0_{\rm Tph}=\sum_{\mathbf q,\alpha=1,2}\hbar cq(a^\dagger a+\frac12)$ (E Eq. 214).
  - $\mathbf A(\mathbf r)=\sum_{\mathbf q\alpha}\sqrt{\frac{2\pi c\hbar}{q\Omega}}\mathbf e_\alpha(\mathbf q)(a_{\mathbf q\alpha}+a^\dagger_{-\mathbf q\alpha})e^{i\mathbf q\mathbf r}$ (E Eq. 216).
  - $H^0_{\rm Lph}=\frac1{8\pi}\int\mathbf E_l^2=\frac12\int\!\!\int\frac{[\rho_e+\rho_{\rm ion}][\rho_e'+\rho_{\rm ion}']}{|\mathbf r-\mathbf r'|}$, $\rho_{\rm ion}=-\nabla\cdot\mathbf P$ (E Eqs. 227, 234–235).
  - $H_0$ (E Eq. 236) and $V'$ (E Eq. 237) printed in full.
- **Interactive viz A "Helmholtz splitter":**
  - Choose $\mathbf q$ (two angles) and an arbitrary $\tilde{\mathbf A}$ (three sliders). Show the triad $\mathbf e_1=\mathrm{normalize}(\mathbf n-(\mathbf n\cdot\hat{\mathbf q})\hat{\mathbf q})$, $\mathbf e_2=-\hat{\mathbf q}\times\mathbf e_1$, $\mathbf e_3=\hat{\mathbf q}$ and the projections.
  - Magnetic energy $\propto q^2(\tilde A_1^2+\tilde A_2^2)$ (the longitudinal part has no $\mathbf B$).
  - What to notice: only $\mathbf A_t$ costs magnetic energy, so only it oscillates.
- **Interactive viz B "Oscillator vs drift":**
  - Phase space $(X,P)$ for the transverse mode $H=\frac{\hbar cq}4(X^2+P^2)$ (circles) and the longitudinal mode $H=\frac{\hbar cq}4P^2$ (straight drift, no restoring force).
  - Controls: $q$, initial point.
  - What to notice: longitudinal "photons" are not particles; their energy is fixed by the charges (Gauss) and reappears as Coulomb.
- **Interactive viz C "Coupling map":** nodes electrons, ions, $\mathbf A_t$, $\mathbf A_l$ (Coulomb). Edges are labelled with E Eq. 237 terms. Toggle "electrostatic limit ($c\to\infty$)" removes $\mathbf A_t$ edges and leaves Giustino's Coulomb-only graph.
- **Static figure:** E Fig. 6 (polarization triad); a table of the four minimal-coupling terms (electron $\mathbf p\cdot\mathbf A$, $A^2$; ion $Z^*\pi\cdot\mathbf A$, $Z^*Z^*AA$).
- **Sources:** [E3] Sec. 4 start, [E4] Eqs. 134–143, [E5] Eqs. 166–237, V1–V2.

#### §12. The polar crystal as a dielectric (classical EPP)
- **Learning goals:**
  - The diatomic crystal in a field gives a Lorentz oscillator with ionic plasma frequency $\omega_p^2=4\pi Z^{*2}/(\Omega_{\rm uc}M)$.
  - LO = zeros, TO = poles of $\epsilon$; LST.
  - "LO–TO splitting" vs "LO shift".
  - Only $\omega_{LO}^2-\omega_{TO}^2=\omega_p^2/\epsilon_\infty$ is model-independent; the local-field model affects individual frequencies.
  - Link to $W=W_{\rm e}+W_{\rm ph}$.
- **Key equations:**
  - $P=-\frac{Z^{*2}}{\Omega_{\rm uc}M}\frac{E}{\omega^2-\omega_O^2}$ (E Eq. 150).
  - $\omega_p^2=\frac{4\pi Z^{*2}}{\Omega_{\rm uc}M}$ (E Eq. 152).
  - $\epsilon(\omega)=\epsilon_\infty-\frac{\omega_p^2}{\omega^2-\omega_O^2+\frac13\omega_p^2}$ (E Eq. 158) $=\epsilon_\infty\frac{\omega^2-\tilde\omega_{LO}^2}{\omega^2-\tilde\omega_{TO}^2}$ (E Eq. 159).
  - $\tilde\omega^2_{LO}=\omega_O^2+(\frac1{\epsilon_\infty}-\frac13)\omega_p^2$, $\tilde\omega^2_{TO}=\omega_O^2-\frac13\omega_p^2$ (E Eq. 160).
  - $\frac{\epsilon(0)}{\epsilon_\infty}=\frac{\tilde\omega_{LO}^2}{\tilde\omega_{TO}^2}$ (E Eq. 165).
  - $W/v=1/\epsilon(\omega)$, $W_{\rm ph}/v=1/\epsilon-1/\epsilon_\infty$ [derived, E4 4.10].
- **Interactive viz "Reststrahlen lab":**
  - Plots (a) Re/Im $\epsilon(\omega)$; (b) reflectance $R=|(\sqrt\epsilon-1)/(\sqrt\epsilon+1)|^2$; (c) loss $-{\rm Im}\,1/\epsilon$; (d) $W/v$ split into $1/\epsilon_\infty$ and the lattice part.
  - Frequency axis: wavenumber (cm$^{-1}$) and meV.
  - JS: $\epsilon=\epsilon_\infty(\omega_{LO}^2-\omega^2-i\gamma\omega)/(\omega_{TO}^2-\omega^2-i\gamma\omega)$.
  - Controls: $\omega_{TO}$, $\omega_{LO}$, $\epsilon_\infty$, $\gamma$ (1–20 cm$^{-1}$). Second mode "build from microscopic": $Z^*$, $\Omega$, $\mu$, $\omega_O$, local-field on/off.
  - Presets (cm$^{-1}$: TO/LO/$\epsilon_\infty$/$\gamma$): GaAs 267/292/10.9/2.5; NaCl 164/264/2.34/8; SiC 796/972/6.52/5; hBN 1370/1610/4.87/5.
  - LST-implied $\Omega_p=\sqrt{\epsilon_\infty(\omega_{LO}^2-\omega_{TO}^2)}$: GaAs ≈390, NaCl 317, SiC 1424, hBN 1866 cm$^{-1}$. Microscopic: $\Omega_p^2[{\rm cm^{-2}}]=4.9206\times10^7Z^{*2}/(\Omega[{\rm Å^3}]\mu[{\rm u}])$; GaAs with $Z^*=2.07$ gives 359 cm$^{-1}$.
  - Level diagram: bare $\omega_O$, $\tilde\omega_{TO}$, $\tilde\omega_{LO}$ vs $\omega_p$. Warn when $\tilde\omega_{LO}<\omega_O$ ($\epsilon_\infty>3$, Born-charge local field) and when $\omega_p^2>3\omega_O^2$ (Cochran soft-mode instability).
  - What to notice: $R\approx1$ between TO and LO; loss peak at LO; $\epsilon(0)/\epsilon_\infty$ = LST; the difference of squares is independent of the local-field model.
- **Static figure:** E Figs. 3–5 (chain with first neighbours, two masses + one spring, charge-neutral dipole); arrows for LO ($\mathbf E=-4\pi\mathbf P$) and TO ($\mathbf E=0$).
- **Sources:** [E4] Eqs. 144–165, Concepts 4.3–4.10.

#### §13. RPA photon self-energy: electrons and ions on the same footing
- **Learning goals:**
  - Matter enters only via $\Pi_{jj}=\Pi^{\rm e}+\Pi^{\rm ion}$.
  - Dia + para, and why the dia term is needed: gauge invariance, $\Pi(\omega\to0)\to0$ for insulators.
  - Electronic $\Pi$ = Adler–Wiser with current vertices.
  - Ionic $\Pi$ = phonon propagator between Born-charge current vertices, exact in the harmonic approximation.
  - The oscillator-strength sum rule from completeness.
- **Key equations:**
  - $D^t=D^{0t}+D^{0t}\Pi^{\rm RPA}D^t$ (E Eq. 277).
  - $\Pi^{\rm e}=-\frac{e^2}{mc}\rho^{\rm e}\delta\delta_{\alpha\beta}+\frac{i}{\hbar c}\int dt\,e^{i\omega t}\langle Tj^{\rm e}_\alpha j^{\rm e}_\beta\rangle$ (E Eq. 279).
  - $\Pi^{\rm RPA,e}=\Pi^{\rm para,e}(\omega)-\Pi^{\rm para,e}(0)$ (E Eq. 281).
  - $\Pi^{\rm para,e}_{\alpha\beta\mathbf G\mathbf G'}=-\frac2{c\Omega}\sum_{\mathbf k nm}\frac{f_{n\mathbf k}-f_{m\mathbf k+\mathbf q}}{\hbar\omega+E_{n\mathbf k}-E_{m\mathbf k+\mathbf q}+i\eta}j^\alpha(\mathbf G)[j^\beta(\mathbf G')]^*$ (E Eq. 282), $j^\alpha=\frac{e\hbar}{2im}\{\phi^*\partial_\alpha\phi-(\partial_\alpha\phi^*)\phi\}$ (E Eq. 284).
  - $j_{\alpha\sigma}(\mathbf q)=\sum_{a\mu}\frac{Z^*_{a\alpha\mu}}{\sqrt{m_a}}e_{a\mu\sigma}(\mathbf q)$ (E Eq. 295).
  - $\Pi^{\rm ion,para}=-\frac1{c\Omega_{\rm uc}}\sum_\sigma\frac{\omega_\sigma^2}{\omega^2-\omega_\sigma^2+i\eta}j_{\alpha\sigma}j^*_{\beta\sigma}$ (E Eq. 297).
  - $\Pi^{\rm ion,dia}=-\frac1{\Omega_{\rm uc}c}\sum_{a\mu}\frac{Z^*_{a\alpha\mu}Z^*_{a\beta\mu}}{m_a}$ (E Eq. 304).
  - $\sum_\sigma j_{\alpha\sigma}j^*_{\beta\sigma}=\sum_{a\mu}Z^*_{a\alpha\mu}Z^*_{a\beta\mu}/m_a$ (E Eq. 309).
  - $\Pi^{\rm ion,RPA}=-\frac{\omega^2}{c\Omega_{\rm uc}}\sum_\sigma\frac{j_{\alpha\sigma}j^*_{\beta\sigma}}{\omega^2-\omega_\sigma^2+i\eta}$ (E Eq. 311).
  - $\epsilon=\epsilon_\infty+\frac{4\pi c}{\omega^2}\Pi^{\rm ion}=\epsilon_\infty+\frac{4\pi}{\Omega_{\rm uc}}\sum_\sigma\frac{j_{\alpha\sigma}j^*_{\beta\sigma}}{\omega_\sigma^2-\omega^2}$ (Gonze–Lee IR form [derived]).
- **Interactive viz A "Para + dia = gauge invariant":**
  - Notation $P\equiv-4\pi c\,\Pi$: $P^{\rm dia}=\Omega_p^2$; $P^{\rm para}=\Omega_p^2\omega_{TO}^2/(\omega^2-\omega_{TO}^2+i\gamma\omega)$; $P^{\rm RPA}=\Omega_p^2\omega^2/(\omega^2-\omega_{TO}^2+i\gamma\omega)$; $\epsilon=\epsilon_\infty+\Omega_p^2/(\omega_{TO}^2-\omega^2-i\gamma\omega)$.
  - Toggle "omit dia": $\epsilon_{\rm para}=\epsilon+\Omega_p^2/\omega^2$, a wrong-sign Drude term ($\to+\infty$ as $\omega\to0$).
  - Plot all three $P$ and $\epsilon$ vs $\omega$.
- **Interactive viz B "Electrons: truncated band sums":**
  - Model: $N_{\rm tr}$ transitions $\Delta_i=\Delta_1+(i-1)\delta$, strengths $f_i\propto e^{-(i-1)/3}$ normalized over 20 transitions; $\hbar\omega_p=16.6$ eV (Si valence). $\Delta_1=3.4$ eV, $\delta=1.2$ eV gives $\epsilon(0)\approx12.9$ [derived]; tune $\Delta_1$ for 11.7.
  - Gauge-invariant form: $\epsilon_{\rm GI}=1+\omega_p^2\sum_{i\le N_t}f_i/(\Delta_i^2-\omega^2)$.
  - "dia + truncated para": $\epsilon_{\rm tr}=\epsilon_{\rm GI}-(1-S)\omega_p^2/\omega^2$ with $S=\sum_{i\le N_t}f_i$, a fake Drude metal ($\to-\infty$).
  - Slider $N_t$ = 1–20.
  - Show the f-sum $(2/\pi)\int\omega\,{\rm Im}\,\epsilon\,d\omega$ vs $\omega_p^2$.
- **Interactive viz C "Oscillator-strength sum rule":**
  - 3-atom chain (masses and $Z^*$ sliders, $\sum Z^*=0$ enforced). Bar chart of $|j_\sigma|^2$ per mode at chosen $q$, plus the constant total $\sum_aZ_a^{*2}/m_a$.
  - What to notice: force constants redistribute strength among modes but never change the total.
- **Static figure:** E Fig. 7-type bubbles: electronic bubble (two square current vertices, particle–hole loop) + contact dia; ionic: square – zigzag phonon $F$ – square + contact. Caption: "ionic RPA is exact (harmonic)".
- **Sources:** [E6] 3.5–3.7, V2–V3; [E7] 3.2–3.4, V1–V2; [E8] S5.

#### §14. Charge conservation, the total propagator and polaritons
- **Learning goals:**
  - The $4\times4$ (scalar + vector) structure.
  - Connection formulas (Ward identities) turn density into current vertices and $v$ into $D^{0l}$.
  - $D^0=D^{0l}+D^{0t}$ is the vacuum Maxwell propagator.
  - One Dyson equation for everything; in an isotropic medium L and T decouple.
  - Plasmons and LO from $\epsilon_L=0$; polaritons from $c^2q^2=\omega^2\epsilon_T$.
- **Key equations:**
  - $\Pi_{\rho\rho}=-\frac i\omega\partial_\mu\Pi_{j\rho\mu}$ (E Eq. 341).
  - $\Pi_{\rho j\nu}=-\frac i\omega\partial_\mu\Pi_{jj\mu\nu}$ (E Eq. 344).
  - $\Pi_{\rho\rho}=-\frac c{\omega^2}\partial_\mu\partial'_\nu\Pi_{jj\mu\nu}$ (E Eq. 350).
  - $D^{0l}_{\mu\nu}=-\frac{4\pi c}{q^2}\frac{q_\mu q_\nu}{\omega^2}$ (E Eq. 374).
  - E Eq. 375 (as in §10).
  - $D=D^0+D^0\Pi_{jj}D$, $\Pi_{jj}=\Pi^{\rm e}+\Pi^{\rm ion}$ (E Eqs. 379–382).
  - Isotropic [derived, E8 S1]: $D_L=-\frac{4\pi c}{\omega^2\epsilon_L}$, $D_T=-\frac{4\pi c}{\omega^2\epsilon_T-c^2q^2}$, $\epsilon_{L,T}=1+\frac{4\pi c\Pi_{L,T}}{\omega^2}$.
- **Interactive viz "Photon spectral map":**
  - Heat map of $-{\rm Im}\,D_T$ and $-{\rm Im}\,D_L$ (switch) on $q\in[0,3\,q_{\rm res}]$ (nm$^{-1}$; light line $q=E[{\rm eV}]/197.327$) × $\hbar\omega$.
  - Models for $4\pi c\Pi$:
    - vacuum (0);
    - Drude $-\omega_p^2\omega/(\omega+i\gamma)$;
    - insulator two-level $\omega^2\omega_p^2f/(\Delta^2-\omega^2)$;
    - polar crystal $\epsilon_\infty-1$ background + $\Omega_p^2\omega^2/(\omega_{TO}^2-\omega^2-i\gamma\omega)$ (i.e. $\epsilon_T=\epsilon_\infty+\Omega_p^2/(\omega_{TO}^2-\omega^2-i\gamma\omega)$).
  - Analytic branches (polar): $\omega_\pm^2=\frac{B\pm\sqrt{B^2-4\epsilon_\infty c^2q^2\omega_{TO}^2}}{2\epsilon_\infty}$, $B=\epsilon_\infty\omega_{TO}^2+\Omega_p^2+c^2q^2$.
  - Controls: model, parameters (presets from §12; Drude $\hbar\omega_p$ 1–10 eV), $\gamma$.
  - What to notice: anticrossing of light with TO, a Reststrahlen gap between $\omega_{TO}$ and $\omega_{LO}$, the upper branch at $q=0$ starts at $\omega_{LO}$, $D_L$ peaks at $\omega_{LO}$ for every $q$ (no dispersion, no retardation); Drude gives the plasma edge.
  - Tab "Diagram transformer" (E Figs. 9→10→11): click to (a) replace each $\rho$ vertex by a $j$ vertex and $v$ by $D^{0l}$ (E Eq. 377); (b) merge $D^{0l}+D^{0t}\to D^0$; (c) merge $\Pi^{\rm e}+\Pi^{\rm ion}\to\Pi_{jj}$.
  - Tab "$\Pi_{\rho\rho}$ from $\Pi_{jj}$": compute $\Pi^{\rm ion}_{\rho\rho}=\frac1{\Omega_{\rm uc}}\sum_\sigma|\mathbf q\cdot\mathbf j_\sigma|^2/(\omega^2-\omega_\sigma^2)$ with and without the dia term; only the gauge-invariant version gives a negative static density response.
- **Static figure:** the $4\times4$ block $D^{1,i}=\begin{bmatrix}v\Pi_{\rho\rho}v&v\Pi_{\rho j}D^{0t}\\D^{0t}\Pi_{j\rho}v&D^{0t}\Pi_{jj}D^{0t}\end{bmatrix}$ (E Eq. 319) with arrows to the connection formulas.
- **Sources:** [E7] 3.5–3.6, V5; [E8] 3.2–3.8, V1–V3, V5.

#### §15. Electrodynamic phonon self-energy and the LO–TO splitting
- **Learning goals:**
  - Phonons dressed by an electron-screened photon: $\Omega^{\rm e}=-J^*D^{\rm e}J$.
  - The dynamical ionic vertex $J\propto\omega j$.
  - The same Dyson structure as G Eq. 137.
  - Without electrons at $q=0$: the shift $4\pi|j|^2/\Omega_{\rm uc}$ (E Eq. 410) is the retarded upper polariton with an unscreened ionic plasma frequency.
  - The electrostatic LO–TO splitting needs $c\to\infty$ first (§17).
  - Order of limits.
- **Key equations:**
  - $D^{\rm e}=D^0+D^0\Pi^{\rm e,0}D^{\rm e}$ (E Eq. 388).
  - $J_{\mu\sigma}(\mathbf q,\omega)=\sqrt{\frac{\hbar}{2c\Omega_{\rm uc}\omega_\sigma(\mathbf q)}}\,\omega\,j_{\mu\sigma}(\mathbf q)$ (E Eq. 392).
  - $F=F^0+F^0\Omega^{\rm e}F$ (E Eq. 398); $\Omega^{\rm e}_{\sigma\sigma'}=-\sum_{\alpha\beta}J^*_{\alpha\sigma}D^{\rm e}_{\alpha\beta}J_{\beta\sigma'}$ (E Eq. 399).
  - $F_\sigma=\frac1\hbar\frac{2\omega_\sigma}{\omega^2-\omega_\sigma^2-2\omega_\sigma\Omega^{\rm e}_{\sigma\sigma}/\hbar}$ (E Eq. 405).
  - $\tilde\omega^2_\sigma(0)=\omega^2_\sigma(0)+\frac{4\pi}{\Omega_{\rm uc}}\sum_\alpha|j_{\alpha\sigma}(0)|^2$ (E Eq. 410).
  - General $\mathbf q$ without electrons [derived, E9]: $\omega^2=\omega_\sigma^2+\omega_{p\parallel}^2+\omega_{p\perp}^2\frac{\omega^2}{\omega^2-c^2q^2}$.
- **Interactive viz "Order of limits":**
  - Plot dispersion $\omega(q)$, $q\in[10^{-7},10^{-2}]$ Å$^{-1}$ (log), with LO, TO and polariton branches.
  - JS: $\omega_L^2=\omega_T^2+\omega_p^2/\epsilon_\infty$ with $\omega_p^2=\epsilon_\infty(\omega_L^2-\omega_T^2)$ unscreened; $\omega_\pm^2=\frac12[(\omega_L^2+c^2q^2/\epsilon_\infty)\pm\sqrt{(\omega_L^2+c^2q^2/\epsilon_\infty)^2-4\omega_T^2c^2q^2/\epsilon_\infty}]$.
  - Toggle "electrons off ($\epsilon_\infty=1$)" reproduces E Eq. 410 at $q=0$, e.g. GaAs ≈56 meV ($Z^*=2.07$) vs measured LO 36.2 meV.
  - Speed-of-light multiplier $\lambda_c$ ($1$–$10^6$) shows the polariton region shrinking to $q=0$ and the plain LO/TO pair emerging.
  - Presets ($\omega_T$/$\omega_L$ meV, $\epsilon_\infty$): GaAs 33.2/36.2/10.9; NaCl 20.3/32.7/2.34; hBN 170/200/4.9.
  - What to notice: at $q\to0$ with finite $c$, the "LO" is the upper polariton and TO is photon-like; for $cq\gg\omega$ you get LO at $\omega_L$ and TO at $\omega_T$; E Eq. 410 is not the electrostatic splitting (no $1/\epsilon_\infty$, isotropic sum).
  - Tab "Why resum?": $F_N=F^0\sum_{n=0}^{N}(F^0\Omega)^n$ vs exact $F^0/(1-F^0\Omega)$ for one mode with $\Omega=$ const (E Eq. 409). Controls $N$ (0–20). Notice that no finite order moves the pole.
- **Static figure:** E Fig. 13 (phonon Dyson chain: zigzag $F^0$ – teal $J$ – blue $D^{\rm e}$ – teal $J^*$ – zigzag); E Fig. 12-type expansion of $\Pi^{\rm ion}$ (E Eqs. 386–387).
- **Sources:** [E9] 3.2–3.3, V1, V5.

#### §16. The electrodynamic Fröhlich vertex
- **Learning goals:**
  - The e–ph vertex is generated by photon exchange between an electronic current vertex and an ionic current vertex.
  - It depends on frequency, unlike the adiabatic $g(\mathbf k,\mathbf q)$.
  - Bare ($D^0$) and screened ($D^{\rm e}$) versions; the self-energy keeps one of each.
  - For $c\to\infty$ it reduces to Verdi–Giustino with $\epsilon^\infty\to\hat\epsilon^{\rm e}(\mathbf q,\omega)$.
  - TO phonons couple only near the light cone.
  - The charge–current vertex identity (continuity).
- **Key equations:**
  - $J^\mu_{n\mathbf k,m\mathbf k+\mathbf q}(\omega)=\sqrt{\frac1{c\Omega}}\frac{\hbar\omega}{E_{n\mathbf k}-E_{m\mathbf k+\mathbf q}}j^\mu_{n\mathbf k,m\mathbf k+\mathbf q}$ (E Eq. 414).
  - $g^{\rm dyn,0}_{nm\sigma}=i\sum_{\alpha\beta}[J^\alpha_{n\mathbf k,m\mathbf k+\mathbf q}]^*D^0_{\alpha\beta}J_{\beta\sigma}$ (E Eq. 415).
  - $g^{\rm dyn}_{nm\sigma}=i\sum_{\alpha\beta}J^*_{\alpha\sigma}D^{\rm e}_{\alpha\beta}J^\beta_{n\mathbf k,m\mathbf k+\mathbf q}$ (E Eq. 416).
  - $\Delta\Omega^{\rm e}_{\sigma\sigma'}=\sum_{\mathbf k nm}g^{\rm dyn}_{nm\sigma}\chi^0_{n\mathbf k,m\mathbf k+\mathbf q}g^{\rm dyn,0}_{nm\sigma'}$ (E p. 92).
  - $\rho_{n\mathbf k,m\mathbf k+\mathbf q}=-\sum_\mu\hbar q_\mu\frac{j^\mu}{E_{n\mathbf k}-E_{m\mathbf k+\mathbf q}}$ (E Eq. 428; read $e\rho$ on the left).
  - $\epsilon^{\rm e}_{\mu\nu}=\delta_{\mu\nu}-\sum_\alpha D^0_{\mu\alpha}\Pi^{\rm e,0}_{jj\alpha\nu}$ (E Eq. 432).
  - $\lim_{c\to\infty}D^{\rm e}_{\mu\nu}=-\frac{4\pi c}{\omega^2}\frac{q_\mu q_\nu}{\mathbf q\cdot\hat\epsilon^{\rm e}\cdot\mathbf q}$ (E Eq. 446).
  - $\lim_{c\to\infty}g^{\rm dyn,0}=\frac{4\pi i}{\Omega_{\rm uc}}\sum_{a\mu\nu}\sqrt{\frac{\hbar}{2Nm_a\omega_\sigma}}\frac{q_\mu Z^*_{a\mu\nu}e_{a\nu\sigma}}{q^2}[\rho_{n\mathbf k,m\mathbf k+\mathbf q}]^*$ (E p. 95).
  - $\lim_{c\to\infty}g^{\rm dyn}=\frac{4\pi i}{\Omega_{\rm uc}}\sum_{a\mu\nu}\sqrt{\frac{\hbar}{2Nm_a\omega_\sigma}}\frac{q_\mu Z^*_{a\mu\nu}e^*_{a\nu\sigma}}{\mathbf q\cdot\hat\epsilon^{\rm e}(\mathbf q,\omega)\cdot\mathbf q}\rho_{n\mathbf k,m\mathbf k+\mathbf q}$ (E p. 98).
  - Relation [derived, E9]: $g^{\rm dyn}=-[g^{\rm dyn,0}]^*\,q^2/(\mathbf q\cdot\hat\epsilon^{\rm e}\cdot\mathbf q)$.
- **Interactive viz "LO vs TO coupling across the light line":**
  - Plot normalized kernels vs $q/q_0$ (log $10^{-3}$–$10^{3}$), $q_0=\sqrt{\epsilon_\infty}\,\omega/c$:
    - $K_{LO}=-\frac{\omega^2}{4\pi c}D^{\rm e}_L=1/\epsilon_\infty$, independent of $c$;
    - $|K_{TO}|=\frac{1/\epsilon_\infty}{|(q/q_0)^2-1-i\gamma|}$ [derived, E8 S1, E9 V2].
  - Second panel: bare vs screened static vertex magnitudes, $|g^{\rm dyn,0}|=\epsilon_\infty|g^{\rm dyn}|$, with VG numbers from §4 (GaAs ratio ≈10.9).
  - Controls: $\hbar\omega$ (phonon, 10–200 meV), $\epsilon_\infty$ (1–15), $\gamma$ (0.001–0.1), $\lambda_c$.
  - What to notice: LO coupling is the same at every $q$ (Coulomb); TO coupling is resonant at the light line and decays as $(q_0/q)^2$; as $c\to\infty$, $q_0\to0$ and TO decouples. The screened static vertex is VG; the bare one is $\epsilon_\infty$ times larger, which is why using two bare vertices would double count.
- **Static figure:** E Fig. 14: phonon self-energy bubble = teal $J$ – blue $D^{\rm e}$ – blue $J^{\rm e}$ – $\chi^0$ loop – blue $J^{\rm e}$ – magenta $D^0$ – teal $J$, with brackets labelling $g^{\rm dyn}$ (left, screened) and $g^{\rm dyn,0}$ (right, bare). Side by side with G Fig. 1(c) (dot $g^{\rm b}$ and disk $g$).
- **Sources:** [E9] 3.4–3.6, V2; [E3] V3; [G4] V5.

#### §17. The electrostatic limit: recovering Giustino, plus LO-phonon–plasmon coupling
- **Learning goals:**
  - As $c\to\infty$ only the $ll$ sector survives, so $D^{\rm e}\to-\frac{c}{\omega^2}qqW^{\rm e}$ with $W^{\rm e}=4\pi/(\mathbf q\cdot\hat\epsilon^{\rm e}\cdot\mathbf q)$ = RPA, no local fields.
  - The phonon self-energy becomes $W^{\rm e}|\mathbf q\cdot\mathbf j|^2$; its static version is Gonze–Lee.
  - Keeping the dynamic version gives coupled LO-phonon–plasmon modes in doped semiconductors.
  - The adiabatic replacement $W^{\rm e}(\omega)\to W^{\rm e}(0)$ is valid only when plasmons are far from phonons.
- **Key equations:**
  - $\lim_{c\to\infty}D^{\rm e}_{\mu\nu}=D^{{\rm e},ll}_{\mu\nu}$ (E Eq. 434).
  - $D^{\rm e}_{11}\equiv W^{\rm e}$ (E Eq. 437).
  - $\lim_{c\to\infty}D^{\rm e}_{\mu\nu}(\mathbf q,\omega)=-\frac{c}{\omega^2}q_\mu q_\nu W^{\rm e}(\mathbf q,\omega)$ (E Eq. 441).
  - $W^{\rm e}=\frac{4\pi}{\mathbf q\cdot\hat\epsilon^{\rm e}\cdot\mathbf q}$ (E Eq. 445).
  - $\lim\Omega^{\rm e}_{\sigma\sigma}=\frac{\hbar W^{\rm e}(\mathbf q,\omega)}{2\Omega_{\rm uc}\omega_\sigma}|\sum_\alpha q_\alpha j_{\alpha\sigma}|^2$ (E Eq. 452).
  - $\omega^2=\omega_\sigma^2+\frac{W^{\rm e}(\mathbf q,\omega)}{\Omega_{\rm uc}}|\mathbf q\cdot\mathbf j_\sigma|^2$ (E Eq. 453).
  - $\tilde\omega_\sigma^2=\omega_\sigma^2+\frac{W^{\rm e}(\mathbf q,0)}{\Omega_{\rm uc}}|\mathbf q\cdot\mathbf j_\sigma|^2$ (E Eq. 454).
  - With $W^{\rm e}=4\pi/(\mathbf q\cdot\epsilon_\infty\cdot\mathbf q)$: $\tilde\omega^2=\omega^2_\sigma+\sum\frac{e^*_{a\mu\sigma}}{\sqrt{m_a}}\big[\frac{4\pi}{\Omega_{\rm uc}}\frac{(\mathbf q\cdot Z^*_a)_\mu(\mathbf q\cdot Z^*_b)_\nu}{\mathbf q\cdot\epsilon_\infty\cdot\mathbf q}\big]\frac{e_{b\nu\sigma}}{\sqrt{m_b}}$ = first-order Gonze–Lee [derived, E9].
- **Interactive viz "LO-phonon–plasmon coupled modes":**
  - Plot (a) mode energies vs carrier density $n$ ($10^{16}$–$10^{19}$ cm$^{-3}$, log) for E Eq. 453 (dynamic) vs E Eq. 454 with static Thomas–Fermi screening; (b) loss function $-{\rm Im}\,1/\epsilon_{\rm tot}$ heat map.
  - Dynamic: pole of $\epsilon_{\rm e}(\omega)+\Omega_p^2/(\omega_T^2-\omega^2)=0$, $\epsilon_{\rm e}=\epsilon_\infty-\omega_{pl}^2/(\omega(\omega+i\gamma))$. Closed form $\omega_\pm^2=\frac12[(\omega_L^2+\tilde\omega_{pl}^2)\pm\sqrt{(\omega_L^2+\tilde\omega_{pl}^2)^2-4\tilde\omega_{pl}^2\omega_T^2}]$, $\tilde\omega_{pl}^2=\omega_{pl}^2/\epsilon_\infty$, $\omega_{pl}^2=ne^2/(\varepsilon_0m^*)$ (SI).
  - Static TF: $\tilde\omega^2=\omega_T^2+\Omega_p^2q^2/(\epsilon_\infty q^2+k_{TF}^2)$, $k_{TF}^2=3ne^2/(2\varepsilon_0E_F)$ (SI), with $q$ slider (0.001–0.1 Å$^{-1}$).
  - Controls: $n$, $m^*$ (0.02–1), $\gamma$, $q$.
  - GaAs preset: 33.3/36.2 meV, $\epsilon_\infty=10.9$, $m^*=0.067$; $n=10^{18}$ cm$^{-3}$ gives $\hbar\tilde\omega_{pl}=43.4$ meV; anticrossing near $5\times10^{17}$ cm$^{-3}$ (Mooradian–Wright [lit]).
  - What to notice: the static approximation only screens the LO back toward TO; the dynamic one gives two branches ($L_-\to\omega_T$, $L_+\to$ plasmon-like). The "3D screening is practically static" argument (E p. 99) fails near the anticrossing.
- **Static figure:** limit chain "$D^{\rm e}\to D^{{\rm e},ll}\to W^{\rm e}\to$ Gonze–Lee"; E Fig. 10 (four sectors $ll,lt,tl,tt$, with $lt,tl,tt$ greyed for $c\to\infty$).
- **Sources:** [E9] 3.6–3.7, V3; [G3] (Bardeen–Pines link), [G4] (coupled plasmon–phonon remark, p. 23).

#### §18. 2D crystals: no LO–TO splitting at Γ
- **Learning goals:**
  - A sheet of point dipoles couples to $D^{\rm e}(\mathbf Q,\omega;z=0,z'=0)$.
  - In the electrostatic limit, $W^{\rm 2D}=v_Q/\epsilon_{2D}$ makes the LO shift vanish linearly at $Q\to0$ (no splitting at $\Gamma$), with the same $r_0$ that controls Keldysh exciton screening.
  - The doped case, and retardation via the sheet propagator.
  - Connect to the reader's 2D RPA/BSE work (same $W$).
- **Key equations:**
  - $\Pi^{\rm ion}_{jj}(\mathbf Q,\omega,z,z')=\frac{\Omega_{\rm uc}}{S_{\rm uc}}\Pi^{\rm ion}_{jj}(\mathbf Q,\omega)\delta(z)\delta(z')$ (E Eq. 455).
  - $\Omega^{\rm e}_{\sigma\sigma'}(\mathbf Q,\omega)=-\sum J^*_{\alpha\sigma}D^{\rm e}_{\alpha\beta}(\mathbf Q,\omega,0,0)J_{\beta\sigma'}$, $J=\sqrt{\frac{\hbar}{2cS_{\rm uc}\omega_\sigma}}\omega j$ (E Eqs. 458–459).
  - $W^{\rm e}(\mathbf Q,\omega,z,z')=v_Qe^{-Q|z-z'|}+v_Q^2e^{-Q(|z|+|z'|)}\Pi^{\rm e}_{\rho\rho}$, $v_Q=2\pi/Q$ (E Eqs. 461–462).
  - $W^{\rm e,2D}=\frac{v_Q}{\epsilon^{\rm e}_{2D}}$, $\epsilon^{\rm e}_{2D}=1-v_Q\Pi^{\rm e,0}_{\rho\rho}$ (E Eqs. 469–470).
  - $\tilde\omega_\sigma^2(\mathbf Q)=\omega_\sigma^2+\frac1{S_{\rm uc}}\frac{v_Q}{\epsilon^{\rm e}_{2D}(\mathbf Q,0)}|\sum_{\alpha=x,y}j_{\alpha\sigma}Q_\alpha|^2$ (E Eq. 472).
  - Model: $\epsilon_{2D}=\epsilon_{\rm env}+r_0Q$ (Rytova–Keldysh, [lit] parameterization).
  - Sheet propagator [derived, E8 S4]: $D^0_{\mu\nu}(\mathbf Q,\omega;0,0)=\frac{2\pi c}{\omega^2\beta}(\frac{\omega^2}{c^2}\delta_{\mu\nu}-Q_\mu Q_\nu)$, $\beta=\sqrt{Q^2-\omega^2/c^2}$.
- **Interactive viz "2D LO dispersion":**
  - Plot $\hbar\omega_{LO}(Q)$ for $Q\in[0,0.5]$ Å$^{-1}$ with TO flat and the 3D-bulk LO for comparison.
  - JS: $(\hbar\omega)^2=(\hbar\omega_T)^2+C\,Q/(\epsilon_{\rm env}+r_0Q)$, $C=2\pi\cdot14.40\cdot4.18\times10^{-3}Z^{*2}/(S_{\rm uc}\mu)=0.3782\,Z^{*2}/(S_{\rm uc}[{\rm Å^2}]\mu[{\rm u}])$ eV$^2$Å.
  - hBN preset: $S=5.43$ Å$^2$, $Z^*=2.7$, $\mu=6.10$ → $C=0.0832$; $\hbar\omega_T=170$ meV, $r_0=10$ Å, $\epsilon_{\rm env}=1$; $Q=0.1$ Å$^{-1}$ → 181.8 meV; slope ≈0.25 eV Å at $Q\to0$; saturation $\sqrt{\omega_T^2+C/r_0}\approx193$ meV.
  - MoS$_2$-like preset (weakly polar): $Z^*\approx1$, $\hbar\omega_T\approx48$ meV, $r_0\approx40$ Å [lit].
  - Controls: $\epsilon_{\rm env}$ (1–10; substrate/encapsulation), $r_0$ (5–80 Å), $Z^*$, $\hbar\omega_T$.
  - Doped tab: solve $(\epsilon_{\rm env}+r_0Q)\omega^4-\omega^2[(\epsilon_{\rm env}+r_0Q)\omega_T^2+\omega_{pl}^2+C'Q]+\omega_{pl}^2\omega_T^2=0$ with $C'=C/\hbar^2$ and $(\hbar\omega_{pl})^2=689.4\,n[{\rm Å^{-2}}]\,Q[{\rm Å^{-1}}]/(m^*/m_e)$ eV$^2$.
  - Retardation tab [derived, flagged]: replace $v_Q\to\frac{2\pi}{Q}\sqrt{1-\omega^2/c^2Q^2}$ near the light line; radiative width of order 0.04 meV for hBN [derived].
  - What to notice: no LO–TO splitting at $\Gamma$, linear onset with slope $\propto Z^{*2}/\epsilon_{\rm env}$; saturation set by $r_0$ (the same $r_0$ as in Keldysh exciton binding); encapsulation reduces the slope; doping produces 2D plasmon–phonon anticrossing.
- **Static figure:** sheet of dipoles at $z=0$ with field lines leaking into the environment; side-by-side 3D (constant LO shift) vs 2D (linear in $Q$). Cite Sohier et al. 2017 (agreement stated as "probably identical").
- **Sources:** [E9] 3.8, V4; [E8] S4, V4; [E7] V4.

---

## 6. Subtleties and risks (do not propagate)

1. **$N_p$ placement in Giustino.** G Eqs. 1, 37, 38 put $N_p^{-1/2}$ outside an $N_p$-independent $g$. G Eqs. 2, 57, 58, 61 and 180 put $N_p$ inside the square root, which is inconsistent. On the page, drop $N_p$ (Eq. 38 convention) in Part A. EPP puts $N^{-1/2}$ inside $g$ and nothing in $H$.
2. **Gaussian vs SI.**
   - EPP Sec. 3 is SI, Sec. 4 onward is Gaussian (unannounced), and Eq. 131 is a hybrid.
   - Giustino is SI.
   - Pick one system per panel. Convert with $e^2\to e^2/4\pi\varepsilon_0$ and $4\pi\to1/\varepsilon_0$.
   - Numeric constants for the page: $e^2/4\pi\varepsilon_0=14.3996$ eV Å; $e^2/\varepsilon_0=180.95$ eV Å.
3. **Born-charge dimension.** EPP's $Z^*$ includes $e$; Giustino's is dimensionless. This plus the missing $e$ in EPP's $\rho_e$ explains the author's Croatian margin note (p. 96: "compared to Giustino you are missing $e^2$"): one $e$ is in $Z^*$, one is missing from E Eqs. 420/428.
4. **$\rho_e$ number vs charge density.** E Eqs. 234, 235, 237 (last term), 312, 420, 428 use $\Psi^\dagger\Psi$ where a charge density $e\Psi^\dagger\Psi$ is required. E Eq. 428 should read $e\rho=-\hbar\mathbf q\cdot\mathbf j/\Delta E$.
5. **Two different $\mathbf j_e$.** E Eq. 166 is the kinetic current $e\sum\mathbf v\delta$. E Eq. 238 is the paramagnetic current. Do not mix them.
6. **Missing basis phase.** EPP places dipoles at the cell origin. E Eqs. 121–131 should carry $e^{-i(\mathbf q+\mathbf G)\cdot\mathbf r^0_a}$. The same factor is missing in Verdi–Giustino Eq. 4 (Giustino footnote 7; G Eq. 180 has it). EPP cites "F. Giustino, PRL 115, 176401"; the correct citation is Verdi & Giustino, PRL 115, 176401 (2015).
7. **E Eq. 410 is not the LO–TO splitting.**
   - E Eq. 410 is the $q=0$, finite-$c$, electron-free result: the upper polariton, with an isotropic $\sum_\alpha|j_\alpha|^2$ and an **unscreened** $\omega_p^2$ (GaAs ≈56 meV vs LO 36 meV). It also shifts TO-like modes.
   - The electrostatic LO–TO splitting is E Eq. 454, obtained with $c\to\infty$ first. It has $1/\epsilon_\infty$ and $|\hat{\mathbf q}\cdot\mathbf j|^2$.
   - The text's comparison of Eq. 410 to the Quantum Espresso nonanalytic term needs this caveat.
8. **Wording slips in EPP.**
   - "Retarded limit $c\to\infty$" (p. 94) means *non-retarded/instantaneous*.
   - "Photon self-energy/photon energy" before E Eqs. 406–407 means *phonon*.
   - E Eq. 391 should read $\Pi^{\rm ion,0}$.
   - p. 64 cites Eq. 354 where Eq. 291 is meant.
9. **Local field and Born charges.** EPP Sec. 4.2.1 uses the Lorentz local field with Born (not Szigeti) charges and claims validity for all symmetries (true only for cubic sites). This gives $\tilde\omega_{LO}<\omega_O$ when $\epsilon_\infty>3$. Only $\tilde\omega_{LO}^2-\tilde\omega_{TO}^2=\omega_p^2/\epsilon_\infty$ is robust. The Dyson-equation route (Secs. 5–7) has no local field, so there TO = bare $\omega_\sigma$. State which model each panel uses.
10. **E Eq. 15 phase.** EPP's dynamical matrix equals $D^{\rm dm}(-\mathbf q)$, and is called "symmetric" although it is only Hermitian. Use Giustino's $D_{12}=-(K_1+K_2e^{-iqa})/\sqrt{m_1m_2}$ in viz. E Eq. 4 lacks the on-site term (acoustic sum rule).
11. **E Eq. 260 factor.** The printed transverse propagator has a duplicated $2\omega_{\mathbf q\alpha}$. Use E Eq. 261, which is correct.
12. **Longitudinal propagator, two statements.** $D^{0l}=0$ (E Eq. 263) holds in the interaction picture, $\omega\neq0$. $D^{0l}=-\frac{4\pi c}{q^2}\frac{qq}{\omega^2}$ (E Eq. 374) is the physical Coulomb propagator, recovered via E Eq. 235 and the connection formulas. Both are valid; explain why, or readers will see a contradiction.
13. **E Eq. 276** sums from $n=1$; it should start at $n=0$ to include $D^{0t}$.
14. **E Eq. 301** ($\nabla\cdot\mathbf j^{\rm ion}=0$, "no free charges") is wrong as stated: ionic bound charge $-\nabla\cdot\mathbf P$ exists. The result dia $=-$para$(0)$ is still correct, via completeness (E Eqs. 307–309).
15. **E Eq. 440** power counting: $D^{lt}\sim c^{-1}$, not $c^0$. The conclusion (only $ll$ survives) holds.
16. **Local fields.** EPP Eqs. 380 onward drop $\mathbf G,\mathbf G'$ (no crystal local fields). Electronic $\Pi$ is RPA without xc. For the 2D reader: no LFE in $z$ means a strictly 2D model.
17. **2D "$\Omega\to S$" is insufficient.** A sheet needs the $q_z$-integrated propagator $D^0(\mathbf Q,\omega;0,0)$ (E8 S4). E Eqs. 463–472 drop out-of-plane (ZO, $Z^*_{zz}$) dipoles, although E Eq. 300 contains them. Label the 2D panel "in-plane modes, electrostatic".
18. **EPP $\chi^0$ (p. 92) is not the Lindhard function.** It differs by a term that vanishes only for a complete band set. The vertex $J^{\rm e}\propto1/\Delta E$ is singular for degenerate/intraband transitions. Do not implement it naively.
19. **"$n_B=0$" rule (E Sec. 4.5).** It is valid for free bosons and for harmonic photon/phonon self-energies, where the response is $T$-independent. It is not valid for electron self-energies (G Eq. 158 needs $n_B$).
20. **What EPP omits.** The Lagrangian (E Eq. 134) has no short-range deformation-potential coupling, since $V_{\rm cr}$ is frozen at equilibrium. Only the polar channel is extended. Bare phonons must be built from short-range (analytic) force constants to avoid double counting (p. 21). EPP does not derive electron self-energies.
21. **Which quantity is screened.**
    - Giustino's phonon self-energy (G Eqs. 141, 145) uses one bare and one screened vertex. His FM electron self-energy (G Eq. 156) uses two screened vertices. EPP's $\Delta\Omega^{\rm e}$ uses one bare ($D^0$) and one screened ($D^{\rm e}$) photon line.
    - G Eq. 145 with $g^{\rm b}g^*(0)$ vs G Eq. 146 with $|g|^2$ differ by $\epsilon(\mathbf q,0)$. Only a consistent RPA treatment settles it.
    - Fröhlich uses $1/\epsilon^\infty-1/\epsilon^0$ (ionic part). $W_{\rm e}$ uses $\epsilon^\infty$, $W$ uses $\epsilon(\omega)$.
    - EPP puts $\epsilon^\infty$ in ad hoc in Sec. 3 and regenerates it in Sec. 7.
22. **Giustino corrections.**
    - G Eq. 127: the sign should be $+Z_\kappa$ with $\epsilon_{\rm e}^{-1}(\mathbf r',\mathbf r;0)$.
    - G Eq. 100: the index order should be $\epsilon_{\rm e}^{-1}(31)$.
    - G Eq. 108 silently truncates at one phonon.
    - Fig. 3: the caption says $\hbar\omega_0=100$ meV, but the curves match 50 meV. The lifetime axis is inconsistent with $\tau=\hbar/2\Gamma$ (6.6 fs for 50 meV).
    - G Eq. B6 should read $x/(\sqrt2l)$.
    - The semiclassical current $J$ is missing $1/\Omega$.
    - BTE G Eqs. 202–203 have normalization issues.
    - G Eqs. 210/200 have index-order issues.
    - The Saitta formula is a Fermi-surface average.
    - The $\gamma$ in G Eq. 185 is a HWHM.
    - $N_F$ is per spin.
23. **Name clashes.** $D$ (nuclear correlator / phonon GF vs photon propagator), $\Pi$ (phonon vs photon self-energy), $\Omega$ (cell volume vs crystal volume vs dynamical matrix vs phonon self-energy), $\sigma$ (branch, photon polarization, $l/t$ label, conductivity). Rename on the page.
24. **Time-ordered vs retarded.** EPP and Giustino write time-ordered $+i\eta$ forms. Viz must use retarded $\omega^2-\omega_0^2+i\gamma\omega$ (E p. 88 digression; G Eq. 158).
25. **Corrected sign in the "omit diamagnetic" demo.** Omitting the dia term gives $\epsilon_{\rm para}=\epsilon_{\rm correct}+\Omega_p^2/\omega^2$, which tends to $+\infty$ (wrong-sign Drude). The E7 note's V1 statement ($-\Omega_p^2/\omega^2$) is wrong; E8 S5 is right. The $-(1-S)\omega_p^2/\omega^2$ "fake metal" belongs to the different case "dia + truncated para".
26. **Sohier et al. agreement.** EPP says "very similar (probably identical)" for E Eqs. 454/471–472. Present it as agreement in form, not as a proven identity.
27. **Fig. 14 colours** in EPP do not match the text. Heisenberg evolution on p. 93 is missing $\hbar$.
28. **Two meanings of "retardation".** Giustino means $\omega$-dependent electronic screening (Bardeen–Pines). EPP means finite $c$ (and also $\omega$-dependent $\hat\epsilon^{\rm e}$). Always say which.
29. **Bardeen–Pines Im $\epsilon$ and Eq. 145.** The dynamic vertex is complex. Approximation (iv) in Giustino, $g(\omega)\to g(0)$, is "ubiquitous but never mentioned". EPP's vertices keep $\omega$-dependence, and this is a selling point.
30. **Literature numbers.** Diamond/Si ZPR, H$_3$S/MgB$_2$ $\alpha^2F$, graphene G band, MoS$_2$ $r_0$, Mooradian–Wright and Si phonon energies are [lit]. Flag them as illustrative.
