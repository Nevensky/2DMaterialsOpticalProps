# 2D crystals beyond in-plane modes: out-of-plane dipoles, environment and radiative widths

Labels used below: **[est]** established in the literature; **[derived]** algebra done for this note from EPP's own equations, not printed in EPP; **[proposal]** a suggested extension of EPP; **[spec]** speculation. Numbers marked *illustrative* use the hBN preset of the primer ($S_{\rm uc}=5.43$ Å², $\mu=6.10$ u, $Z^*_\parallel=2.7$, $\hbar\omega_{\rm T}=170$ meV, $r_0=10$ Å).

## 1. What is missing in EPP

EPP §8 (E.455–E.472, pp. 100–102; primer [§B.10](../../sections/B10-2d.html)) treats one **planar** sheet. The ionic bubble is pinned to $z=0$ (E.455), the photon line obeys a $z$-resolved Dyson equation (p. 100, unnumbered), and then three restrictions are made at once:

1. **Out-of-plane dipoles are dropped.** E.300 and E.305 (pp. 64, 67) contain all $\alpha,\beta\in\{x,y,z\}$, hence $Z^*_{zz}$. Before E.463 the "divergent local term $\sim\delta(z-z')$" is discarded and only $\alpha,\beta=x,y$ are kept, so ZO-type modes (hBN ZO, MoS₂ A₂″) get no self-energy.
2. **Atomic heights are lost.** The printed E.300/E.305 put every dipole at $z=0$; the corrected E.300′/E.305′ ([§Z.3](../../sections/Z03-corrections.html)) carry $\boldsymbol\rho^0_a$ and $z_a$, but §8 is evaluated only for $z_a=0$. Buckled sheets (MoS₂, InSe), multilayers and heterostructures are absent.
3. **Strictly 2D, electrostatic, freestanding.** The bubble is squeezed into $\delta(z)$ (E p.101′, RPA, no local fields, no $f_{xc}$), $c\to\infty$ (E.460), and there is no substrate, so $\epsilon_{\rm env}$ and finite thickness are missing. The $\epsilon_{\rm env}+r_0Q$ form in the primer is a literature parameterization, not EPP.

Also missing: the finite-$c$ sheet propagator (the primer derives only its in-plane block), radiative widths, 2D phonon polaritons, and any 2D electron–phonon vertex. E.415–E.416 are never carried to 2D, so there is no 2D Fröhlich coupling, mobility or polaron ([§Z.1](../../sections/Z01-outlook.html)).

## 2. Why it matters

- **LO dispersion and environment [est].** LO and TO are degenerate at Γ. The LO slope scales with $Z^{*2}$, is reduced by the environment, and grows linearly with layer number (Sohier *et al.* 2017). For illustrative hBN the freestanding slope is 0.245 eV Å, and 0.054 eV Å for $\epsilon_{\rm env}=4.5$ (primer). EELS on monolayer hBN on Cu finds a slope below the freestanding prediction, attributed to substrate screening (Li *et al.* 2024).
- **Out-of-plane dipoles in transport [est].** In fluorinated graphene, omitting the dipolar coupling overestimates the 300 K hole mobility by 99%, and a 3D treatment underestimates it by 33% (Deng *et al.* 2021, preprint figures). Quadrupoles change 2D mobilities by more than 75% (Poncé *et al.* 2023). 2D Fröhlich coupling in TMDs is much stronger than earlier ab initio estimates (Sohier *et al.* 2016).
- **ZO dispersion [derived].** The dropped $zz$ term is not only a $Q$-independent contact shift. Its nonlocal part *softens* ZO modes linearly in $Q$ (§4). With an illustrative $Z^*_{zz}=0.5$ and $\hbar\omega_{\rm ZO}=0.1$ eV, the slope is $-0.014$ eV Å, about $-1.4$ meV at $Q=0.1$ Å⁻¹.
- **Radiative widths [derived].** In-plane polar modes of freestanding hBN have $\hbar\gamma_0\approx42$ μeV at $Q=0$ (primer). ZO modes are dark at normal incidence and radiate only in p-polarization, with a width $\propto x^2$ (§4). The same angular laws govern bright and $z$-polarized 2D excitons, whose radiative coupling is far larger (Wang *et al.* 2018).
- **For an RPA/GW/BSE practitioner.** The phonon kernel *is* your $W$: $r_0$, $\epsilon_{\rm env}$ and Coulomb truncation enter the LO dispersion exactly as they enter Rytova–Keldysh exciton binding. Without truncation a supercell restores a spurious LO–TO splitting (Sohier *et al.* 2017), the artifact that also spoils untruncated 2D BSE. A wrong $r_0$ convention (§5) moves the LO plateau $\sqrt{\omega_{\rm T}^2+C'/r_0}$.

## 3. What standard theory does

Giustino's review is three-dimensional: the Verdi–Giustino vertex (G.180) with $\bq\cdot\epsilon^\infty\cdot\bq$, the Gonze–Lee term (cited, not printed), $\Pi^{\rm A}$ (G.125) and $W_{\rm e}=v+vP_{\rm e}W_{\rm e}$ (G.94) ([§A.4](../../sections/A04-polar.html), [§A.6](../../sections/A06-phonon-se.html)). Nothing in it is specific to 2D.

Established 2D first-principles practice [est]:

- **DFPT with a Coulomb cutoff along $z$** removes periodic images and handles gating (Sohier, Calandra, Mauri 2016, 2017).
- **Long-range force constants.** The 2D dipole term (Sohier *et al.* 2017) was made exact by Royo and Stengel (2021): 2×2 (in-plane/out-of-plane) screening, $z$-dependent macroscopic potentials, out-of-plane dynamical charges, quadrupoles.
- **Long-range e–ph vertices.** 2D dipoles with out-of-plane charges (Deng *et al.* 2021); a unified 3D↔2D Fröhlich vertex without truncation (Sio and Giustino 2022); dipoles, quadrupoles and a Berry-connection term for mobilities of MoS₂, InSe, hBN (Poncé *et al.* 2023).
- **Environment.** Layer-resolved screening (QEH, Andersen *et al.* 2015), remote screening of Fröhlich scattering (Sohier, Gibertini, Verstraete 2021), $r_0$ from the sheet polarizability (Cudazzo *et al.* 2011; Berkelbach *et al.* 2013).
- **Electrodynamics and excitons.** The 2D phonon polariton is the LO phonon (Rivera, Christensen, Narang 2019); phonon screening of the BSE kernel (Filip, Haber, Neaton 2021).

## 4. How it could enter the EPP formalism

**(a) Full sheet propagator [derived from E.375].** Write $\tilde\partial=(iQ_x,iQ_y,\partial_z)$, $\tilde\partial'=(-iQ_x,-iQ_y,\partial_{z'})$, $k^2=\epsilon\omega^2/c^2$ and $\beta=\sqrt{Q^2-k^2-i\eta}$ with ${\rm Re}\,\beta>0$. Then, in a homogeneous medium $\epsilon$,
$$D^0_{\mu\nu}(\bQ,\omega;z,z')=\frac{4\pi c}{\epsilon\,\omega^2}\Big[k^2\delta_{\mu\nu}-\tilde\partial_\mu\tilde\partial'_\nu\Big]\frac{e^{-\beta|z-z'|}}{2\beta}.$$
In-plane it gives $\frac{2\pi c}{\epsilon\omega^2\beta}(k^2\delta_{\mu\nu}-Q_\mu Q_\nu)$ (the primer's result); $xz\propto-iQ_x\,{\rm sgn}(z-z')\,e^{-\beta|z-z'|}$, nonzero only between different heights; and $zz=\frac{4\pi c}{\epsilon\omega^2}\big[\frac{Q^2}{2\beta}e^{-\beta|z-z'|}-\delta(z-z')\big]$, whose $\delta$ is EPP's "divergent local term". For $c\to\infty$ it reduces to E.460 with $W\to v_Qe^{-Q|z-z'|}$. A layered environment replaces $e^{-\beta|z-z'|}/2\beta$ by the Fresnel (Sommerfeld) Green's function with $r_s$, $r_p$ **[proposal]**.

**(b) Ionic sheet self-energy [proposal].** Keep E.300′/E.305′ with all three Cartesian components and every height $z_a$, in the subtracted form (the analogue of E.311):
$$\Pi^{\rm ion}_{jj,\alpha\beta}(\bQ,\omega;z,z')=-\frac{\omega^2}{cS_{\rm uc}}\sum_{\sigma=1}^{3s}\sum_{ab}\frac{j^a_{\alpha\sigma}(\bQ)\,j^{b*}_{\beta\sigma}(\bQ)}{\omega^2-\omega_\sigma^2(\bQ)}\,\delta(z-z_a)\,\delta(z'-z_b).$$

**(c) Phonon self-energy [proposal, generalizes E.458].** Use atom-resolved vertices $J^a_{\alpha\sigma}=\sqrt{\hbar/(2cS_{\rm uc}\omega_\sigma)}\,\omega\,j^a_{\alpha\sigma}$:
$$\Omega^{\rm e}_{\sigma\sigma'}(\bQ,\omega)=-\sum_{ab}\sum_{\alpha\beta\in\{x,y,z\}}J^{a*}_{\alpha\sigma}\,D^{\rm e}_{\alpha\beta}(\bQ,\omega;z_a,z_b)\,J^b_{\beta\sigma'},$$
with $D^{\rm e}$ from the $z$-resolved Dyson equation of p. 100 built on (a). The phonon Dyson equation E.467 is unchanged.

**(d) Electrostatic kernel [derived].** As $c\to\infty$ the photon line becomes $-\frac{c}{\omega^2}K$, with $K_{\alpha\beta}(z,z')=\tilde\partial_\alpha\tilde\partial'_\beta W^{\rm e}(\bQ,\omega;z,z')$. For the bare kernel,
$$\partial_z\partial_{z'}\,v_Qe^{-Q|z-z'|}=4\pi\,\delta(z-z')-Q^2v_Qe^{-Q|z-z'|}.$$
Outside the contact, $K_{xx}+K_{yy}+K_{zz}=0$, which is Laplace's equation. Whatever the in-plane longitudinal mode gains, the $z$ mode loses. For a planar, $\sigma_h$-symmetric sheet, in-plane and ZO modes decouple. The in-plane modes recover E.472. The ZO modes obey, to leading order in $Q$,
$$\omega^2_{\rm ZO}(\bQ)=\bar\omega^2_{\rm ZO}-\frac{2\pi Q}{S_{\rm uc}}\,|j_{z\sigma}|^2 ,$$
where $\bar\omega_{\rm ZO}$ contains the contact (depolarization) term. The in-plane sheet bubble does not screen this term, because a $z$-dipole sheet has $\phi(z=0)=0$. Electronic out-of-plane screening needs a $z$-dipole channel in $\Pi^{\rm e}$: $\Pi^{\rm e}(z,z')\approx\Pi_{\rho\rho}\delta(z)\delta(z')+\alpha_{zz}\,\delta'(z)\delta'(z')$ **[proposal]**. This is EPP's version of Royo and Stengel's 2×2 screening.

**(e) Radiative widths [derived].** Insert (a) with $\beta=-i\kappa$ inside the cone and $x=cQ/(\sqrt{\epsilon}\,\omega)$. The in-plane LO-like and TO-like modes have widths $\gamma_0\sqrt{1-x^2}$ and $\gamma_0/\sqrt{1-x^2}$ (primer). The ZO mode has
$$\gamma_{\rm ZO}=\gamma_{0,z}\,\frac{x^2}{\sqrt{1-x^2}},\qquad\hbar\gamma_{0,z}=\frac{C_z}{\hbar c\sqrt{\epsilon}},\qquad C_z=\frac{2\pi e^2\hbar^2Z^{*2}_{zz}}{S_{\rm uc}\,\mu}\ \ (\text{diatomic}),$$
with $Z^*_{zz}$ in units of $|e|$, as for the primer's $C$. Illustratively, $Z^*_{zz}=0.5$ gives $\hbar\gamma_{0,z}\approx1.4$ μeV. The p-polarized poles of the full sheet $D$ are the 2D phonon polaritons; outside the cone they are the confined LO branch itself (Rivera *et al.*).

**(f) 2D Fröhlich vertex [proposal].** Carry E.416 to the sheet with $D^{\rm e}(\bQ,\omega;z,z_a)$. In the static limit, up to the phase convention of E p.98′,
$$g^{\rm 2D}_{mn\sigma}(\bk,\bQ)\propto\frac{ie}{S_{\rm uc}}\sum_a\Big\langle m\bk+\bQ\Big|\big[\bQ\cdot\bj^a_{\sigma\parallel}-i\,j^a_{z\sigma}\partial_{z_a}\big]W^{\rm e}(\bQ,0;z,z_a)\Big|n\bk\Big\rangle .$$
With $W=v_Qe^{-Q|z-z_a|}$ both terms stay finite as $Q\to0$; Sohier's 2D Fröhlich vertex is the in-plane one. The $z$ term is odd in $z-z_a$, so it vanishes within a $\sigma_h$-even band. It survives for Janus or buckled layers and for remote ZO scattering of carriers in an adjacent layer.

**(g) Multilayers [proposal].** Use $W^{\rm e}(\bQ,\omega;z_l,z_{l'})$ from a layer-resolved RPA, as in QEH. Then $\Omega^{\rm e}$ becomes a matrix in layer modes, and the in-phase combination should reproduce the slope $\propto N$.

## 5. Pitfalls

- **Double counting with DFPT force constants.** 2D-cutoff DFPT already contains the $Q$-linear dipole terms and the analytic ZO contact shift. Subtract the 2D nonanalytic part (Sohier 2017; Royo–Stengel 2021) before adding $\Omega^{\rm e}$. Assign the contact term once, either to the short-range constants or to $\Omega^{\rm e}$ with a regularized form factor.
- **Which $Z^*_{zz}$.** The open-circuit (fixed $D_z$) and short-circuit (fixed $E_z$) charges differ by $\epsilon_{zz}$. If $\epsilon_{zz}$ is also in the kernel, using the screened charge double-counts the screening.
- **Double counting with $D$.** The Coulomb interaction is $D^{0l}$; never add $v$ on top. In a BSE built on this $W$, drop the $\bG=0$ exchange, which is that photon line ([§B.5](../../sections/B05-photon-se.html)).
- **Gauge invariance.** E.310 (dia $=-$para(0)) holds only with all $3s$ modes and all three components. Keeping only in-plane modes in para but the full $Z^*Z^*/m$ in dia breaks it, so use the subtracted form. For electrons, apply E.281 per $(G_z,G_z')$.
- **Sum rules.** $\sum_aZ^*_{a,zz}=0$ must hold. A violation produces a spurious monopole $\propto v_Q\sim1/Q$.
- **Conventions.** In Gaussian units $r_0=2\pi\alpha_{\rm 2D}$, and from a supercell of height $L$, $r_0=L(\epsilon^{\rm sc}_\parallel-1)/2$ (Berkelbach *et al.* 2013). Some authors write $\epsilon_{\rm env}(1+r_0'Q)$. In SI, $r_0=\alpha_{\rm 2D}/2\varepsilon_0$. Use $\epsilon_{\rm env}$ at the phonon frequency, not its static value, for polar substrates.
- **Regularization and convergence.** The $\delta(0)$ and the ${\rm sgn}$ kink at $z_a=z_b$ need atomic form factors or many $G_z$. The $1/\sqrt{1-x^2}$ divergence at the light line is integrable but needs a finite $\eta$. The $\bQ\to0$ sampling must resolve $Q^*=\epsilon_{\rm env}/r_0$.

## 6. State of the art (verified)

1. P. Cudazzo, I. V. Tokatly, A. Rubio, PRB **84**, 085406 (2011), [10.1103/PhysRevB.84.085406](https://doi.org/10.1103/PhysRevB.84.085406). Exact 2D screened potential; $r_0$ from the sheet polarizability.
2. T. C. Berkelbach, M. S. Hybertsen, D. R. Reichman, PRB **88**, 045318 (2013), [10.1103/PhysRevB.88.045318](https://doi.org/10.1103/PhysRevB.88.045318). Keldysh-screened exciton and trion model with ab initio $r_0$.
3. K. Andersen, S. Latini, K. S. Thygesen, Nano Lett. **15**, 4616 (2015), [10.1021/acs.nanolett.5b01251](https://doi.org/10.1021/acs.nanolett.5b01251). QEH model: layer-resolved screening of heterostructures.
4. T. Sohier, M. Calandra, F. Mauri, PRB **94**, 085415 (2016), [10.1103/PhysRevB.94.085415](https://doi.org/10.1103/PhysRevB.94.085415). 2D Fröhlich coupling in TMD monolayers with cutoff DFPT. Erratum: PRB 96, 159904.
5. T. Sohier, M. Calandra, F. Mauri, PRB **96**, 075448 (2017), [10.1103/PhysRevB.96.075448](https://doi.org/10.1103/PhysRevB.96.075448). 2D DFPT with Coulomb truncation and gating.
6. T. Sohier, M. Gibertini, M. Calandra, F. Mauri, N. Marzari, Nano Lett. **17**, 3758 (2017), [10.1021/acs.nanolett.7b01090](https://doi.org/10.1021/acs.nanolett.7b01090). Breakdown of LO–TO splitting; slope $\propto N$, set by the environment.
7. N. Rivera, T. Christensen, P. Narang, Nano Lett. **19**, 2653 (2019), [10.1021/acs.nanolett.9b00518](https://doi.org/10.1021/acs.nanolett.9b00518). The 2D phonon polariton is the LO phonon; confinement and losses in hBN.
8. T. Deng *et al.*, PRB **103**, 075410 (2021), [10.1103/PhysRevB.103.075410](https://doi.org/10.1103/PhysRevB.103.075410). Ab initio 2D dipolar e–ph coupling including out-of-plane dipoles.
9. M. Royo, M. Stengel, PRX **11**, 041027 (2021), [10.1103/PhysRevX.11.041027](https://doi.org/10.1103/PhysRevX.11.041027). Exact quasi-2D long-range screening and force constants, with out-of-plane charges and quadrupoles.
10. T. Sohier, M. Gibertini, M. J. Verstraete, PRMaterials **5**, 024004 (2021), [10.1103/PhysRevMaterials.5.024004](https://doi.org/10.1103/PhysRevMaterials.5.024004). Remote free-carrier screening suppresses Fröhlich scattering.
11. M. R. Filip, J. B. Haber, J. B. Neaton, PRL **127**, 067401 (2021), [10.1103/PhysRevLett.127.067401](https://doi.org/10.1103/PhysRevLett.127.067401). Phonon screening of the BSE kernel.
12. W. H. Sio, F. Giustino, PRB **105**, 115414 (2022), [10.1103/PhysRevB.105.115414](https://doi.org/10.1103/PhysRevB.105.115414). Unified 2D/3D Fröhlich vertex without truncation.
13. S. Poncé, M. Royo, M. Stengel, N. Marzari, M. Gibertini, PRB **107**, 155424 (2023), [10.1103/PhysRevB.107.155424](https://doi.org/10.1103/PhysRevB.107.155424). Dipoles and quadrupoles in 2D mobilities (MoS₂, InSe, hBN, …).
14. J. Li *et al.*, Nat. Commun. **15**, 1938 (2024), [10.1038/s41467-024-46229-4](https://doi.org/10.1038/s41467-024-46229-4). EELS of the nonanalytic LO in monolayer hBN; substrate-reduced slope.
15. G. Wang *et al.*, RMP **90**, 021001 (2018), [10.1103/RevModPhys.90.021001](https://doi.org/10.1103/RevModPhys.90.021001). Review of 2D excitons: screening, environment, light–matter coupling.

## 7. First steps

1. **Derive and check (days).** Write out §4(a)–(e). Verify the limit E.460 → E.472, the tracelessness of $K$, and the primer's numbers (0.245 eV Å, 42 μeV). Then add the ZO branch.
2. **hBN monolayer benchmark (1–2 weeks).** Run cutoff DFPT (Quantum ESPRESSO's 2D Coulomb cutoff) for $Z^*_\parallel$, $Z^*_{zz}$, $\epsilon$ and short-range IFCs. Compare §4(d) for LO and ZO with direct DFPT at small $Q$ and with the EELS slope on Cu (Li 2024). Compute $r_0$ with your own RPA code and compare.
3. **Use your GW code's $W$ (2–4 weeks).** Project the truncated $W_{G_zG_z'}(\bQ,\omega)$ onto atomic heights to get $K(z_a,z_b)$. Test on buckled MoS₂ (E′ and A₂″), InSe and bilayer hBN (the $N$ scaling), against the dipole and dipole-plus-quadrupole models.
4. **Environment (weeks).** Replace $D^0$ by the Fresnel Green's function with $\epsilon_{a,b}(\omega)$ for SiO₂ and hBN, which gives hybrid sheet–substrate modes and remote Fröhlich scattering.
5. **Electrodynamics (months).** Solve the full finite-$c$ sheet Dyson equation, and link it to exciton–phonon and BSE phonon screening.

## 8. Open questions

- Does the ZO contact term belong in the short-range force constants or in $\Omega^{\rm e}$, and does EPP's Lagrangian decide?
- Are point dipoles at heights $z_a$ enough at $O(Q^2)$, or are electronic quadrupoles unavoidable (Poncé 2023)? **[spec]** Probably needed for InSe.
- Polar substrate: background $\epsilon(\omega)$ or an explicit second $\Pi^{\rm ion}$ layer, and when does it matter?
- How does the light-cone physics cross over from $N$ layers to the bulk LO–TO splitting?
- Are 1–40 μeV radiative widths measurable against anharmonic and inhomogeneous broadening? **[spec]** Only in the cleanest encapsulated samples.
