# Anisotropic crystals and polaritons beyond the zone centre

Labels: **[est]** established; **[derived]** algebra done for this note from EPP's equations and checked numerically; **[proposal]** a suggested EPP extension; **[spec]** speculation. Notation and units are those of primer v2 ([§C.0](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#c-conventions)): Gaussian units (as in EPP §§4–8) with $e>0$, $\hbar$ and $c$ explicit; $Z^*$ is dimensionless and $eZ^*$ is the charge of the dipole (EPP's $Z^*$ carries the charge); for a uniaxial crystal $z$ is the optic axis and $x,y$ lie in the plane, as in [§B.9](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b-electrostatic). (B$n$.$m$) are equations of primer v2. Primer v1's E.410a–h, a proposed EPP §7.1.1, are built into v2's [§B.7](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b-phonon-se): E.410a is (B7.14), with its tensor form in the box after (B7.19); E.410b–c are (B7.15)–(B7.16); E.410d is (B7.19); E.410e–f are (B7.20)–(B7.21); E.410g is (B7.17) and E.410h is (B7.22). hBN parameters (a widely used rounded set, also the primer preset of [§A.4](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#a4-direction); its attribution to [7], also made there, is unverified, and published fits differ by 5–20 cm⁻¹): in-plane ($x$) TO/LO $1370/1610$ cm⁻¹ with $\varepsilon_{\infty,xx}=4.87$, and out-of-plane ($z$) $780/830$ cm⁻¹ with $\varepsilon_{\infty,zz}=2.95$. Sensitivity [derived]: each entry of the §2 table lies between its band's TO and LO, so it moves by an amount between the shifts of the two band edges. With the alternative fit 4.90/2.95, 1360/1614 and 760/825 cm⁻¹ (reported by the review; source not re-verified), the upper row moves by −8 to +4 cm⁻¹ and the lower row by −8 to −20 cm⁻¹.

## 1. What is missing in EPP

- **The formalism is tensorial.** E.380 is $3\times3$ at $\mathbf{G}=\mathbf{G}'=0$, E.385 (with the basis phase of (B5.22)) carries full $\mathbf{Z}^*$ tensors, E.431–E.432 define an electronic dielectric tensor ($\mathbf{S}^{\mathrm{e}}$ in [§B.8](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b-frohlich)), and E.445–E.454 contain $\mathbf{q}\cdot\boldsymbol{\varepsilon}^{\mathrm{e}}\cdot\mathbf{q}$. The direction-dependent LO–TO term therefore exists for $c\to\infty$; (B9.12)–(B9.15) give its matrix form ([§B.9](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b-electrostatic)).
- **What is evaluated is isotropic.** §4.2.1 (E.144–E.165, pp. 28–35) assumes diagonal $Z^*$ and force constants (p. 29), treats one Cartesian direction at a time and assumes $\nabla\times\mathbf{E}=0$. §7.1 evaluates only $\mathbf{q}=0$ without electrons (E.406–E.410), where $D^0$ is isotropic (E.408).
  - §4.2.1 also adds the cubic Lorentz field $\frac{4\pi}{3}\mathbf{P}$ (E.151, p. 32; (B3.13)) and states that "the same expression is valid for all other crystal symmetries". For a uniaxial, layered crystal such as hBN the Lorentz factor is a tensor, not $\frac13\mathbf{1}$.
  - In the first-principles route the local field is already inside $Z^*$ and the analytic force constants (the $\mathbf{G}\neq0$ dipole sum), so E.151 must not be added on top ([§B.3](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b-dielectric): the local-field factor is 0 "when the local field is counted inside the force constants, as in density-functional perturbation theory"; [§B.7](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b-phonon-se): "There is no Lorentz local field").
  - It also sets $\mathbf{Z}^*_1\approx-\mathbf{Z}^*_2$ (p. 32). The acoustic sum rule makes this exact for two inequivalent sublattices, as in hBN, but not for α-MoO₃ or β-Ga₂O₃.
- **What is missing at finite $c$:**
  - L–T mixing in $D$;
  - the extraordinary (Fresnel) dispersion and the hyperbolic regime $\varepsilon_{xx}\varepsilon_{zz}<0$;
  - slabs, surfaces and heterostructures. §8 has one $\delta(z)$ ionic sheet (E.455) with a formally retarded, $z$-resolved $D^{\mathrm{e}}$ (p. 100). It evaluates only the electrostatic, in-plane limit (§8.1, E.460–E.472: $\alpha=x,y$, with the contact term discarded);
  - damping beyond $\mathrm{i}\eta$;
  - low-symmetry crystals;
  - a protocol for first-principles inputs.
- **Primer status.** [§B.6](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b-propagator) decouples $D^{\mathrm{L}}$ and $D^{\mathrm{T}}$ for isotropic media only, and notes that in uniaxial hBN the mixed block is proportional to the difference of the in-plane and out-of-plane components of $\Pi$, $\hat{\mathbf{q}}\cdot\Pi\cdot\mathbf{t}_1\propto(\Pi_{xx}-\Pi_{zz})\sin\theta\cos\theta$, and that this mixing produces the hyperbolic polaritons. [§B.7](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b-phonon-se) derives the isotropic branches and treats hBN as isotropic.
  - [§B.9](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b9-sectors) ("Four sectors, one survivor") derives $D^{\mathrm{e}}=\frac{4\pi c}{\omega^2}(x^2\mathcal{P}^{\mathrm{T}}-\boldsymbol{\varepsilon})^{-1}$ ((B9.3)), with $x=cq/\omega$ and the transverse projector $\mathcal{P}^{\mathrm{T}}=\mathbf{1}-\hat{\mathbf{q}}\hat{\mathbf{q}}^{\mathsf{T}}$, for a constant uniaxial $\boldsymbol{\varepsilon}$. It plots the L, LT and T sectors at angle $\theta$ and shows that LT and T fall as $x^{-2}$.
  - [§B.1](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b-foundations) plots the anisotropic $\mathbf{G}=0$ kernel of the long-range vertex for hBN.
  - New in v2, in the electrostatic limit: the [§A.4](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#a4-direction) figure "The modes at Γ depend on the direction of approach" solves the $3\times3$ problem for bulk hBN with the parameters above, i.e. the $q\gg q_{\mathrm{light}}(\omega)$ limit behind the §2 table, and [§B.9](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b9-matrix) ("First order, or the full matrix") shows for wurtzite GaN that the first-order (diagonal) formula fails for oblique $\hat{\mathbf{q}}$.
  - What is missing is $\Pi^{\mathrm{ion}}$ inside the tensor. The polariton and hyperbolic dispersion therefore appear only through the determinant (B7.19) ([§B.7](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b-phonon-se)).

## 2. Why it matters

- **Hyperbolic bands [est].** hBN has $\varepsilon_{zz}<0<\varepsilon_{xx}$ at 780–830 cm⁻¹ (Type I) and $\varepsilon_{xx}<0<\varepsilon_{zz}$ at 1370–1610 cm⁻¹ (Type II).
  - In hBN, volume-confined modes reach 1/86 of the free-space wavelength with quality factors up to 283 [7], the film dispersion tunes with thickness down to a few layers [6], and isotopic purification triples the lifetime [11].
  - In α-MoO₃, the in-plane polaritons are up to 60 times shorter than light and have amplitude lifetimes of about 8 ps [12].
  - In monoclinic β-Ga₂O₃, shear polaritons have frequency-dependent axes [15].
- **hBN numbers [derived].**
  - $\varepsilon_{\mathrm{s},xx}=6.73$ and $\varepsilon_{\mathrm{s},zz}=3.34$ (Lyddane–Sachs–Teller).
  - $\omega_{\mathrm{LO}}^2-\omega_{\mathrm{TO}}^2=4\pi e^2Z^{*2}/(\Omega_{\mathrm{V}}M_{\mathrm{r}}\varepsilon_\infty)$ implies $Z^*_{\mathrm{B},xx}=2.80$ and $Z^*_{\mathrm{B},zz}=0.73$ ($\Omega_{\mathrm{V}}=18.1$ Å³ per BN pair, $M_{\mathrm{r}}=6.10$ amu; §A.4 fits 2.79 and 0.730).
  - The light cone is at $q_{\mathrm{light}}(\omega_{\mathrm{TO}})=1.9\times10^{-4}$ Å⁻¹ in-plane and $0.84\times10^{-4}$ Å⁻¹ out-of-plane.
  - For $q\gg q_{\mathrm{light}}(\omega)$ the directional LO modes are, with $\theta$ the angle between $\hat{\mathbf{q}}$ and the optic axis $z$:

  | $\theta$ | 0° | 15° | 30° | 45° | 60° | 75° | 90° |
  |---|---|---|---|---|---|---|---|
  | upper (cm⁻¹) | 1370 (T) | 1399 | 1463 | 1527 | 1574 | 1601 | 1610 |
  | lower (cm⁻¹) | 830 | 822 | 807 | 794 | 786 | 781 | 780 (T) |

- **Slab regime [derived].** Take a free-standing 10 nm film, with in-plane wavevector $Q_\parallel$ and slab modes $Q_{\parallel,l}$ (P.6).
  - At 1500 cm⁻¹, $Q_{\parallel,0}=4.4\times10^{-3}$ Å⁻¹ $=46\,\omega/c$ (polariton wavelength ≈144 nm), and the rays run at 38° to the plane.
  - At 1400 cm⁻¹, $Q_{\parallel,0}=5.4\times10^{-4}$ Å⁻¹ $=6.1\,\omega/c$; that this is $\approx\sqrt{|\varepsilon_{xx}|}\,\omega/c$ is a coincidence. The exact retarded TM solution ($q_z^2=\varepsilon_{xx}(\omega^2/c^2-Q_\parallel^2/\varepsilon_{zz})$, $\zeta=\sqrt{Q_\parallel^2-\omega^2/c^2}$, $q_zd=\pi l-2\arctan(|\varepsilon_{xx}|\zeta/q_z)$) shifts it by only +1.3% (+0.03% at 1500 cm⁻¹).
  - The shift grows to 3% at 1390 cm⁻¹ and 13% at 1380 cm⁻¹, where $Q_{\parallel,0}\approx2\omega/c$. For a 10 nm film, retardation matters only within ~10–20 cm⁻¹ of $\omega_{\mathrm{TO},x}$. Since $Q_{\parallel,0}\propto1/d$, the window is wider for thicker films.
  - At 1600 cm⁻¹, $Q_{\parallel,0}=0.063$ and $Q_{\parallel,2}=0.28$ Å⁻¹, so phonon dispersion (nonlocality) matters. Gubbin and De Liberato analyse this nonlocal phonon response for polar films of ~10 nm (Phys. Rev. X 2020; unverified, §6).
- **For 2D GW/BSE practitioners.** hBN is the standard encapsulant.
  - The effective environment $\sqrt{\varepsilon_{xx}\varepsilon_{zz}}$ rises from 3.79 ($\varepsilon_\infty$) to 4.74 ($\varepsilon_{\mathrm{s}}$), about 25% [derived], across the 97–200 meV phonon window. Which of the two values $W$ should see is a frequency question, and EPP's total $D$ answers it.
  - Between the two limits the screening is not an interpolation. Inside 97–103 and 170–200 meV, $\varepsilon_{xx}\varepsilon_{zz}<0$ (−12.5 at 1500 cm⁻¹), so $W(\mathbf{Q}_\parallel,\omega)$ has a hyperbolic-polariton continuum instead of a real $\varepsilon_{\mathrm{eff}}$ [derived]. This is where a frequency-dependent $W$ matters most.
  - Hyperbolic hBN modes cool hot graphene carriers in picoseconds [10].
  - The monolayer limit is the 2D LO dispersion without LO–TO splitting ([9], [§B.10](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b-2d)).

## 3. What standard theory does

- **Giustino.** (G.180) is "the generalization of Fröhlich's model to the case of anisotropic crystalline lattices and multiple phonon modes", with $\mathbf{Z}^*_\kappa$ and $\boldsymbol{\varepsilon}_\infty$ tensors and the cutoff $\mathrm{e}^{-a|\mathbf{q}+\mathbf{G}|^2}$ (p. 015003-35). For the phonon side, §VI.A.3 cites Gonze–Lee.
  - $W=W^{\mathrm{e}}+W^{\mathrm{ph}}$ (G.104–G.105) is dynamically screened. $W^{\mathrm{ph}}(12)=\int W^{\mathrm{e}}(13)\,\chi^{\mathrm{ion}}(34)\,W^{\mathrm{e}}(24)$ contains the nuclear density–density correlation function $\chi^{\mathrm{ion}}$ (G.103; [§A.5](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#a-hedin-baym)), so $W^{\mathrm{ph}}(\omega)$ carries the phonon poles, and $W^{\mathrm{e}}(\omega)$ is frequency dependent too (G.94–G.95).
  - It is built on the instantaneous Coulomb $v$ ($c\to\infty$), so it contains no transverse photons and no retardation.
  - Light enters only as a classical field with an input refractive index $n_{\mathrm{r}}(\omega)$ (G.200, p. 015003-45). Polaritons appear only as the analogy for coupled plasmon–phonon modes, "the electronic analog of photon polaritons" (p. 015003-23); [§B.9](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b9-dynamic) derives these modes, and [plasmon.md](plasmon.md) takes them further. There is no polariton theory.
- **DFPT [4, 5] [est].** DFPT supplies $\mathbf{Z}^*_\kappa$, $\boldsymbol{\varepsilon}_\infty$ and the analytic force constants. The nonanalytic term is the $\mathbf{q}\to0$ limit of the dipole–dipole force constants, $\Phi^{\mathrm{dd}}_{\kappa\alpha,\kappa'\beta}=\frac{4\pi e^2}{\Omega_{\mathrm{V}}}\frac{(\hat{\mathbf{q}}\cdot\mathbf{Z}^*_\kappa)_\alpha(\hat{\mathbf{q}}\cdot\mathbf{Z}^*_{\kappa'})_\beta}{\hat{\mathbf{q}}\cdot\boldsymbol{\varepsilon}_\infty\cdot\hat{\mathbf{q}}}$ ([§B.9](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b9-dd)), and $\boldsymbol{\varepsilon}(\omega)=\boldsymbol{\varepsilon}_\infty+\frac{4\pi e^2}{\Omega_{\mathrm{V}}}\sum_\nu\mathbf{S}_\nu/(\omega_{\mathbf{0}\nu}^2-\omega^2)$ with mode oscillator strengths $\mathbf{S}_\nu$ ((B5.27)). The quasi-LO/TO angular dispersion in uniaxial crystals is classic Raman physics (Loudon 1964, unverified).
- **Polaritons in practice [est].** A two-step route: an ab initio or fitted $\boldsymbol{\varepsilon}(\omega)$ goes into Maxwell's equations (Fresnel, transfer matrices, slab formulas [6, 7, 13]). Classical roots: [1, 2]. Royo–Stengel [14] give exact long-range 2D electrostatics, with $2\times2$ dielectric matrices and the out-of-plane dynamical-charge coupling, but no retardation.

## 4. How it could enter the EPP formalism

**P.A Closed-form propagator [proposal; derived].** At $\mathbf{G}=\mathbf{G}'=0$, E.380 inverts algebraically. With $(D^0)^{-1}=\frac{c}{4\pi}\mathsf{M}^0$ from E.375, where $\mathsf{M}^0=q^2\mathbf{1}-\mathbf{q}\mathbf{q}^{\mathsf{T}}-\frac{\omega^2}{c^2}\mathbf{1}$ is the Maxwell operator of [§X.1](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#x-light-cone),

$$D(\mathbf{q},\omega)=\frac{4\pi}{c}\,\mathsf{M}^{-1},\qquad \mathsf{M}=q^2\mathbf{1}-\mathbf{q}\mathbf{q}^{\mathsf{T}}-\frac{\omega^2}{c^2}\boldsymbol{\varepsilon}(\mathbf{q},\omega),\qquad \boldsymbol{\varepsilon}=\mathbf{1}+\frac{4\pi c}{\omega^2}\big(\Pi^{\mathrm{e}}+\Pi^{\mathrm{ion}}\big),\tag{P.1}$$

and E.385, with the vertex (B5.22), gives

$$\varepsilon^{\mathrm{ion}}_{\alpha\beta}(\mathbf{q},\omega)=\frac{4\pi}{\Omega_{\mathrm{V}}}\sum_{\nu}\frac{j_{\alpha\nu}(\mathbf{q})\,j^*_{\beta\nu}(\mathbf{q})}{\omega^2_{\mathbf{q}\nu}-\omega^2-\mathrm{i}\eta\omega}.\tag{P.2}$$

P.2 is Gonze–Lee's $\boldsymbol{\varepsilon}(\omega)$ with $e^2\mathbf{S}_\nu=\mathbf{j}_\nu\mathbf{j}_\nu^\dagger$, extended to finite $\mathbf{q}$; at $\mathbf{q}\to0$ it is the primer's (B5.27). $\boldsymbol{\varepsilon}$ is the material tensor; EPP's E.432 tensor ($\mathbf{S}^{\mathrm{e}}$ in §B.8) contains $D^0$ and agrees with it only in the longitudinal projection, as §B.8 states. (B7.14) is P.1 with $\boldsymbol{\varepsilon}=\varepsilon_\infty\mathbf{1}$. With $\Pi^{\mathrm{ion}}=0$ and a constant uniaxial $\boldsymbol{\varepsilon}^{\mathrm{e}}$, P.1 is the closed form (B9.3) of [§B.9](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b-electrostatic).

**P.B Pole condition [derived].** Sylvester's identity gives $\det[\text{matrix of (B7.19)}]=\prod_\nu(\omega^2-\omega^2_{\mathbf{q}\nu})\,\det\mathsf{M}[\boldsymbol{\varepsilon}]/\det\mathsf{M}[\boldsymbol{\varepsilon}^{\mathrm{e}}]$, so the modes solve $\det\mathsf{M}=0$, the Fresnel equation. For a uniaxial crystal the ordinary wave has $c^2q^2=\omega^2\varepsilon_{xx}$, and the extraordinary wave has

$$\frac{c^2q^2}{\omega^2}=\frac{\varepsilon_{xx}(\omega)\,\varepsilon_{zz}(\omega)}{\varepsilon_{xx}(\omega)\sin^2\theta+\varepsilon_{zz}(\omega)\cos^2\theta}.\tag{P.3}$$

- **$\mathbf{q}\to0$ at finite $c$:** $\det\boldsymbol{\varepsilon}=0$, so every polar mode sits at its LO frequency, isotropically. This is the anisotropic (B7.20) (hBN: 1610 cm⁻¹ twice and 830 cm⁻¹).
- **$q\gg q_{\mathrm{light}}(\omega)$:** $\varepsilon^{\mathrm{L}}\equiv\hat{\mathbf{q}}\cdot\boldsymbol{\varepsilon}\cdot\hat{\mathbf{q}}=0$. This is the $c\to\infty$ limit of (B7.19), i.e. the $\mathbf{G}=0$ term (B9.12) projected on the modes, (B9.14); it gives the table in §2.
  - For oblique $\hat{\mathbf{q}}$ both E₁u and A₂u have $\hat{\mathbf{q}}\cdot\mathbf{j}_\nu\neq0$ and mix. Diagonalizing the $2\times2$ matrix reproduces the table [derived].
  - The diagonal (B7.21)/E.454 with a tensor $\boldsymbol{\varepsilon}_\infty$ (first order, as in (B9.13)) misplaces the oblique modes. At $\theta=30^\circ$ it gives 1459.7/812.6 cm⁻¹ instead of 1462.8/807.1 cm⁻¹, an error of up to 3 cm⁻¹ (upper) and 5.5 cm⁻¹ (lower, 11% of the Type I band) [derived]. §B.9 now makes the same point for wurtzite GaN.

At large $q$, bulk hyperbolic polaritons are these oblique LO modes ($\nabla\cdot(\boldsymbol{\varepsilon}\cdot\nabla A_0)=0$).

**P.C L–T mixing in $W$ [derived].** In the basis $(\hat{\mathbf{q}},\mathbf{t}_1)$, with $\mathbf{t}_1\perp\hat{\mathbf{q}}$ in the plane of $\hat{\mathbf{q}}$ and the optic axis, $\varepsilon^{\mathrm{LT}}=(\varepsilon_{xx}-\varepsilon_{zz})\sin\theta\cos\theta$ and $\varepsilon^{\mathrm{T}}=\varepsilon_{xx}\cos^2\theta+\varepsilon_{zz}\sin^2\theta$. With $x=cq/\omega$ and $W\equiv-\frac{e^2\omega^2}{c\,q^2}\hat{\mathbf{q}}\cdot D\cdot\hat{\mathbf{q}}$, i.e. $W=e^2D_{\rho\rho}$ (§B.6),

$$W(\mathbf{q},\omega)=\frac{4\pi e^2}{q^2}\,\frac{\varepsilon^{\mathrm{T}}-x^2}{\varepsilon^{\mathrm{L}}(\varepsilon^{\mathrm{T}}-x^2)-(\varepsilon^{\mathrm{LT}})^2}\;\xrightarrow{x\to\infty}\;\frac{4\pi e^2}{\mathbf{q}\cdot\boldsymbol{\varepsilon}\cdot\mathbf{q}}\ (\text{E.445; (B9.7)}).\tag{P.4}$$

- At $q=0.01$ Å⁻¹ and $\theta=45^\circ$, the retarded L–T mixing moves the oblique LO poles by −0.03 cm⁻¹ (at 1526.5 cm⁻¹) and −0.004 cm⁻¹ (at 794.2 cm⁻¹).
- Away from resonances the relative correction $(\varepsilon^{\mathrm{LT}})^2/(\varepsilon^{\mathrm{L}}x^2)$ is ~10⁻⁵–10⁻⁴; with the static $\boldsymbol{\varepsilon}_{\mathrm{s}}$, valid for $\hbar\omega\ll97$ meV, it is $4\times10^{-5}$. Next to a TO pole (170 meV) it reaches $4\times10^{-2}$ undamped and $2\times10^{-2}$ with $\gamma\approx5$ cm⁻¹.
- It is negligible for transport and of order one for $q\lesssim q_{\mathrm{light}}(\omega)$.

**P.D Phonons and vertices [proposal].** E.398–E.399 with $D^{\mathrm{e}}=\frac{4\pi}{c}(\mathsf{M}^{\mathrm{e}})^{-1}$, where $\mathsf{M}^{\mathrm{e}}$ is $\mathsf{M}$ with $\boldsymbol{\varepsilon}^{\mathrm{e}}$. The off-diagonal $\Pi^{\mathrm{ph}}_{\nu\nu'}$ mixes E₁u and A₂u for oblique $\hat{\mathbf{q}}$, and its $c\to\infty$ limit is (B9.12) in the matrix condition (B9.14). The retarded vertex E.416/E.430 ($g^{\mathrm{dyn}}$ of (B8.6)) with the full $D^{\mathrm{e}}$ couples electrons in an adjacent layer to the hyperbolic modes. Its transverse part matters only for $q\lesssim$ a few $\sqrt{|\varepsilon|}\,\omega/c$; at the $q\sim k_{\mathrm{F}}$ that dominate transport and cooling, the LT and T sectors are suppressed by $(q_{\mathrm{light}}(\omega)/q)^2$ (P.C; [§B.9](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b-electrostatic)).

**P.E Slabs and heterostructures [proposal].** Generalize the $z$-resolved Dyson equation of EPP p. 100 ((B10.18)) to the whole field:

$$D(\mathbf{Q}_\parallel,\omega;z,z')=D^0+\!\int\!\mathrm{d}z_1\,\mathrm{d}z_2\,D^0(z,z_1)\big[\Pi^{\mathrm{e}}+\Pi^{\mathrm{ion}}\big](z_1,z_2)\,D(z_2,z'),\qquad \Pi^{\mathrm{ion}}=\sum_{\ell\ell'}\tfrac{\Omega_{\mathrm{V}}}{\Omega_{\mathrm{S}}}\Pi^{\mathrm{ion}}_{\ell\ell'}(\mathbf{Q}_\parallel,\omega)\,f_\ell(z)f_{\ell'}(z')\tag{P.5}$$

- $f_\ell(z)$ is a normalized layer profile ($\int f_\ell\,\mathrm{d}z=1$) of width ~3.3 Å, centred on layer $\ell$. Point sheets $\delta(z-z_\ell)$ are not allowed (next bullet).
- **Contact term.** $D^0$ is the $q_z$ transform of E.375. It is $\propto\mathrm{e}^{-\zeta|z-z'|}/2\zeta$ with $\zeta=\sqrt{Q_\parallel^2-\omega^2/c^2}$, as in (B10.27), and its $zz$ component also contains the contact term $-\frac{4\pi c}{\omega^2}\delta(z-z')$.
  - With point sheets, the $\ell=\ell'$ $zz$ terms would give $\delta(0)$. EPP discards exactly this "divergent local term" before E.463 and keeps $\alpha=x,y$.
  - Here the term must be kept. It is the depolarization field that sets the A₂u frequency and the whole Type I band ($\varepsilon_{zz}<0$).
  - Regularize it with the profiles $f_\ell$, or remove the on-layer self term Ewald-style and add the contact term back analytically.
  - Require that the $N\to\infty$ continuum limit reproduce $\varepsilon_{zz}(\omega)$ and the 780–830 cm⁻¹ band.
- The ionic part is a layer-resolved E.455; $\ell\neq\ell'$ terms come from interlayer force constants.
- The continuum quasi-static limit gives the slab modes of [6] [derived]. For a film of thickness $d$ between media $\varepsilon_1$ and $\varepsilon_3$,

$$Q_{\parallel,l}=\frac{1}{\mathcal{K}d}\Big[\pi l+\arctan\frac{\varepsilon_1}{\varepsilon_{zz}\mathcal{K}}+\arctan\frac{\varepsilon_3}{\varepsilon_{zz}\mathcal{K}}\Big],\qquad \mathcal{K}=\sqrt{-\varepsilon_{xx}/\varepsilon_{zz}},\quad Q_{\parallel,l}>0.\tag{P.6}$$

**Gauge invariance.** Every $\Pi$ must have the form $\Pi^{\mathrm{para}}(\omega)-\Pi^{\mathrm{para}}(0)$ (E.383, E.385; (B5.15)).
- This removes the spurious $1/\omega^2$ Drude term, so $\boldsymbol{\varepsilon}^{\mathrm{ion}}(\omega\to0)$ is finite.
- It makes the static longitudinal response vanish (a static longitudinal $\mathbf{A}$ is pure gauge), consistent with continuity. E.341–E.350 constrain only the longitudinal projections.
- Applied componentwise at finite $\mathbf{q}$, it also drops the $O(q^2)$ static orbital-magnetic response, which is negligible for phonon polaritons.

## 5. Pitfalls

- **Phonon double counting.**
  - Bulk: $\omega_{\mathbf{q}\nu}$ is the DFPT dynamical matrix minus only its $\mathbf{G}=0$ nonanalytic term, because E.380 re-adds only $\mathbf{G}=0$. The $\mathbf{G}\neq0$ dipolar part stays in, so the bare $\omega_{\mathbf{0}\nu}$ is the zero-field TO. Phonons interpolated with Gonze–Lee would count LO–TO twice. The subtraction must use the same $\mathbf{Z}^*$ and $\boldsymbol{\varepsilon}_\infty$ as the force constants.
  - Slab: compute the bare force constants with DFPT and a 2D Coulomb truncation. Then remove the layer-resolved, $\mathbf{Q}_\parallel$-dependent 2D long-range term (in the form of Sohier [9] or Royo–Stengel [14], with the same $Z^*_\ell$ and electronic screening) that P.5 adds back. The bulk $\mathbf{G}=0$ rule does not apply, and a periodic-image supercell gives wrong $\ell\ell'$ dipole couplings.
- **Electron double counting.** Use either the model $\varepsilon_\infty$ ((B7.14)) or the RPA $\Pi^{\mathrm{e}}$, not both.
  - The $\Pi^{\mathrm{e}}$ that enters $D$ at $\mathbf{G}=\mathbf{G}'=0$ must be proper with respect to the macroscopic Coulomb line: $\varepsilon_{\mathrm{M}}=1-v_0\bar\chi_{00}$, with $\bar\chi$ computed with the $\mathbf{G}=0$ Hartree term removed. Without local fields, $\bar\chi_{00}=\chi^0_{00}$. Inserting the full RPA head $\chi_{00}$ counts the $\mathbf{G}=0$ Coulomb twice; this is the RPA analogue of the BSE remark below.
  - A density-response code gives only $\hat{\mathbf{q}}\cdot\Pi^{\mathrm{e}}\cdot\hat{\mathbf{q}}$ (through E.350). At $q\to0$ in an insulator the tensor follows from heads along several directions. Transverse components at finite $q$ need a current–current bubble.
  - RPA without local fields or $f^{\mathrm{xc}}$ is not DFPT's $\varepsilon_\infty$, and a 10% error in $\varepsilon_{\infty,xx}$ moves $\omega_{\mathrm{LO},x}$ by 2.8 meV (22 cm⁻¹), 9% of the Type II band [derived].
  - [est] Local fields across the vdW gap matter most for $\varepsilon_{zz}$, so a $\mathbf{G}=0$ RPA may misplace the Type I band. For graphite, TDDFT reproduces the anisotropic loss spectra for $\mathbf{q}$ off the basal planes only with local fields included (Marinopoulos et al. 2002; unverified, §6).
- **Coulomb double counting.** $D$ already contains $v$ (E.364, E.374). In a BSE kernel, drop the $\mathbf{G}=0$ exchange.
- **Local-field double counting.** Do not add the Lorentz field of E.151 ((B3.13)) to DFPT inputs. It is cubic-only, and $\mathbf{Z}^*$ and the analytic force constants already contain the local field (§1).
- **Sum rules.** Use the phase-corrected vertex (B5.22) and enforce $\sum_\kappa Z^*_{\kappa,\alpha\beta}=0$. Without either, the acoustic modes acquire spurious, cell-dependent dipoles ([§B.1](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b1-acoustic), (B1.17)–(B1.19)).
- **Convergence and regularization.**
  - Keep $\hat{\mathbf{q}}$ in the head of $\Pi^{\mathrm{e}}$.
  - (B9.15) needs an Ewald cutoff and removal of the self term. P.5 needs the regularized $zz$ contact term (§4 P.E).
  - Anharmonic and isotopic damping [11] belong in $F$ (E.398), not in $\eta$ or $\Pi^{\mathrm{e}}$.
- **Approximations.** E.404 fails for mixed modes; use (B7.19) or P.B instead (§4 P.B quantifies the error). P.6 is local; for $Q_\parallel\gtrsim0.05$ Å⁻¹ keep $\omega_{\mathbf{q}\nu}$ and $\mathbf{j}_\nu(\mathbf{q})$ in P.2.

## 6. State of the art (verified)

1. K. Huang, "On the interaction between the radiation field and ionic crystals", Proc. R. Soc. Lond. A **208**, 352 (1951), doi:10.1098/rspa.1951.0166. First phonon-polariton theory.
2. R. Fuchs, K. L. Kliewer, "Optical modes of vibration in an ionic crystal slab", Phys. Rev. **140**, A2076 (1965), doi:10.1103/PhysRev.140.A2076. Slab and surface modes.
3. R. Geick, C. H. Perry, G. Rupprecht, "Normal modes in hexagonal boron nitride", Phys. Rev. **146**, 543 (1966), doi:10.1103/PhysRev.146.543. Polarized IR modes of hBN.
4. X. Gonze, C. Lee, "Dynamical matrices, Born effective charges, dielectric permittivity tensors, and interatomic force constants from density-functional perturbation theory", Phys. Rev. B **55**, 10355 (1997), doi:10.1103/PhysRevB.55.10355. Nonanalytic term, $\mathbf{S}_\nu$, $\boldsymbol{\varepsilon}(\omega)$.
5. S. Baroni, S. de Gironcoli, A. Dal Corso, P. Giannozzi, "Phonons and related crystal properties from density-functional perturbation theory", Rev. Mod. Phys. **73**, 515 (2001), doi:10.1103/RevModPhys.73.515. DFPT $\mathbf{Z}^*$ and $\boldsymbol{\varepsilon}_\infty$.
6. S. Dai *et al.*, "Tunable phonon polaritons in atomically thin van der Waals crystals of boron nitride", Science **343**, 1125 (2014), doi:10.1126/science.1246833. Thickness scaling.
7. J. D. Caldwell *et al.*, "Sub-diffractional volume-confined polaritons in the natural hyperbolic material hexagonal boron nitride", Nat. Commun. **5**, 5221 (2014), doi:10.1038/ncomms6221. Type I and II modes.
8. D. N. Basov, M. M. Fogler, F. J. García de Abajo, "Polaritons in van der Waals materials", Science **354**, aag1992 (2016), doi:10.1126/science.aag1992. Review.
9. T. Sohier, M. Gibertini, M. Calandra, F. Mauri, N. Marzari, "Breakdown of optical phonons' splitting in two-dimensional materials", Nano Lett. **17**, 3758 (2017), doi:10.1021/acs.nanolett.7b01090. 2D limit.
10. A. Principi *et al.*, "Super-Planckian electron cooling in a van der Waals stack", Phys. Rev. Lett. **118**, 126804 (2017), doi:10.1103/PhysRevLett.118.126804. Electrons coupled to hyperbolic modes, treated electrostatically.
11. A. J. Giles *et al.*, "Ultralow-loss polaritons in isotopically pure boron nitride", Nat. Mater. **17**, 134 (2018), doi:10.1038/nmat5047. Damping.
12. W. Ma *et al.*, "In-plane anisotropic and ultra-low-loss polaritons in a natural van der Waals crystal", Nature **562**, 557 (2018), doi:10.1038/s41586-018-0618-9. Biaxial α-MoO₃.
13. G. Álvarez-Pérez *et al.*, "Infrared permittivity of the biaxial van der Waals semiconductor α-MoO₃ from near- and far-field correlative studies", Adv. Mater. **32**, 1908176 (2020), doi:10.1002/adma.201908176. Benchmark tensor.
14. M. Royo, M. Stengel, "Exact long-range dielectric screening and interatomic force constants in quasi-two-dimensional crystals", Phys. Rev. X **11**, 041027 (2021), doi:10.1103/PhysRevX.11.041027. 2D electrostatic benchmark.
15. N. C. Passler *et al.*, "Hyperbolic shear polaritons in low-symmetry crystals", Nature **602**, 595 (2022), doi:10.1038/s41586-021-04328-y. Non-co-diagonal oscillators.

Unverified:
- R. Loudon, Adv. Phys. **13**, 423 (1964). The page and DOI were not confirmed.
- N. Rivera, T. Christensen, P. Narang, "Phonon polaritonics in two-dimensional materials", arXiv:1901.05308. The journal version was not confirmed.
- The next three were supplied by the review with full bibliographic data. They were not re-checked in this revision: the web-search budget was used up and the publisher and Crossref hosts were not reachable.
  - K.-J. Tielrooij *et al.*, Nat. Nanotechnol. **13**, 41 (2018), doi:10.1038/s41565-017-0008-8. Measured out-of-plane cooling of graphene carriers into hBN hyperbolic modes.
  - A. G. Marinopoulos *et al.*, "Anisotropy and interplane interactions in the dielectric response of graphite", Phys. Rev. Lett. **89**, 076402 (2002), doi:10.1103/PhysRevLett.89.076402. Local fields and out-of-plane loss.
  - C. R. Gubbin, S. De Liberato, "Optical nonlocality in polar dielectrics", Phys. Rev. X **10**, 021027 (2020), doi:10.1103/PhysRevX.10.021027. Nonlocal phonon response in nanoscale polar films.

## 7. First steps

1. **Bulk hBN model (1–3 days).** Start from the `inv3`/`sectors` code of [§B.9](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#b9-sectors) (its figure "Four sectors, one survivor") and replace $\boldsymbol{\varepsilon}$ by $\boldsymbol{\varepsilon}^{\mathrm{e}}+\boldsymbol{\varepsilon}^{\mathrm{ion}}(\omega)$ (P.2).
   - Implement P.1–P.3 and reproduce the §2 table (from the $2\times2$ mode matrix, not the diagonal (B7.21)) and the $q\to0$ limit; the [§A.4](https://claude.ai/artifact/VpqcSjXyAFNMcgBmiv3Gid#a4-direction) figure already gives the $q\gg q_{\mathrm{light}}(\omega)$ row values for comparison.
   - Check (B7.19) against $\det\mathsf{M}$, then plot the isofrequency surfaces across $q_{\mathrm{light}}(\omega)$.
2. **First-principles inputs (1–2 weeks).**
   - DFPT: $\mathbf{Z}^*$, $\boldsymbol{\varepsilon}_\infty$ and the zero-field Γ modes of bulk hBN.
   - Your RPA code: take $\varepsilon_{\mathrm{M}}(\hat{\mathbf{q}})=1-v_0\bar\chi_{00}(\hat{\mathbf{q}})$, with the $\mathbf{G}=0$ Hartree term removed, with and without $\mathbf{G}\neq0$ local fields. Do this along $x$ and $z$, and along 45° as a consistency check (uniaxial: it must equal $(\varepsilon_{xx}+\varepsilon_{zz})/2$).
   - Reconstruct the $q\to0$ tensor and insert it as $\frac{4\pi c}{\omega^2}\Pi^{\mathrm{e}}=\boldsymbol{\varepsilon}_{\mathrm{M}}-\mathbf{1}$. Compare with DFPT. For low-symmetry crystals use at least 4–6 directions, including one oblique in the $xz$ plane for $\varepsilon_{xz}$ in the monoclinic case.
   - The transverse components at finite $q$ cannot come from the charge–charge correlator $\Pi^{\mathrm{e}}_{\rho\rho}$. They need a current–current implementation with $\Pi^{\mathrm{para}}(\omega)-\Pi^{\mathrm{para}}(0)$.
   - Propagate the differences to the LO frequencies with P.3.
3. **Slab (4–8 weeks).** Solve P.5 for $N$-layer hBN with the retarded $D^0$.
   - Regularize the $zz$ contact term (profiles $f_\ell$ or Ewald self-term removal) and check the $N\to\infty$ limit against $\varepsilon_{zz}(\omega)$ and the Type I band. Take the bare force constants from 2D-truncated DFPT minus the layer-resolved long-range term (§5). Together these add about 1–2 weeks.
   - Benchmark against P.6 and against the data of [6].
   - Test $N=1$ against EPP §8, [9] and [14].
   - Map where retardation ($Q_\parallel\lesssim3$–$5\,\omega/c$; check against the exact TM root of §2) and nonlocality ($Q_\parallel\gtrsim0.05$ Å⁻¹) matter.
4. **Low symmetry (1–3 months).** α-MoO₃ [12, 13], then β-Ga₂O₃ [15].
   - Start with fitted $\gamma_{\mathbf{q}\nu}$ in $F$ (E.398) and benchmark against [13]. The dispersion, and the frequency-dependent axis rotation of shear polaritons [15], already follow from P.2 with non-co-diagonal $\mathbf{j}_\nu\mathbf{j}_\nu^\dagger$ and any damping.
   - Ab initio anharmonic and isotope linewidths are a separate, longer follow-up. It overlaps the damping-beyond-$\mathrm{i}\eta$ gap.
5. **Electrons (months).**
   - (a) Electrostatic tensor limit. Use the $z$-resolved P.5 in the $c\to\infty$ limit (interlayer $W$ including the hBN $\Pi^{\mathrm{ion}}$). Benchmark the graphene→hBN cooling power against [10] and the measurement of Tielrooij et al. (unverified, §6). This needs the electron energy-loss rate, i.e. an electron self-energy with $D$, from the self-energy topic.
   - (b) Retardation can be tested only where $Q_\parallel\lesssim$ a few $\omega/c$: far-field emission or IR absorption of the stack, and polariton launching.
   - Put the dynamical hBN $W$ into the exciton screening.

## 8. Open questions

- How large are the $\mathbf{G}\neq0$ local-field corrections to $\varepsilon_{\infty,zz}$ of hBN? Do they move the Type I band by more than the 5–20 cm⁻¹ spread of published fits? For graphite, local fields are known to be essential for $\mathbf{q}$ off the basal plane (Marinopoulos et al. 2002; unverified, §6).
- Interface Born charges differ from bulk ones. Should P.5 use layer-dependent $Z^*_\ell$ from slab DFPT?
- How do dynamical quadrupoles [14] enter at finite $c$, at order $O(q)$?
- [spec] Do exciton polaritons and the hyperbolic continuum hybridize measurably in hBN/TMD stacks?

## Review responses

Every reviewer number was re-derived with a short script (exact retarded TM slab root, $2\times2$ mode matrix, extraordinary condition at $q=0.01$ Å⁻¹). The EPP and Giustino page and equation references were checked in the PDFs, and the primer claims in §B.1, §B.3, §B.7, §B.9 and the outlook and corrections sections of primer v1.

1. **Slab retardation (major): accepted.** At 1400 cm⁻¹ the exact correction is +1.34% (+0.025% at 1500, 3.2% at 1390, 13.2% at 1380 cm⁻¹, where $Q_\parallel\approx2.1\,\omega/c$). The §2 bullet and §7 step 3 were rewritten, and "quasi-statics fails" was removed.
2. **P.5 contact term and slab double counting (major): accepted.** EPP p. 101 drops the "divergent local term $\sim\delta(z-z')$" and keeps $\alpha,\beta=x,y$ (E.463–E.464). P.5 now uses profiles $f_\ell(z)$ and keeps the $zz$ term, with an $N\to\infty$ check. §5 has a slab rule, and step 3 is now 4–8 weeks.
3. **RPA head (major): accepted.** Step 2 and §5 now use $\varepsilon_{\mathrm{M}}=1-v_0\bar\chi_{00}$, the reconstruction of the tensor from several directions, and the current–current requirement at finite $q$. Wording refined: without local fields, $\bar\chi_{00}$ reduces to $\chi^0_{00}$. What must not be used is the full reducible head $\chi_{00}$.
4. **Giustino's $W$ (major): accepted.** Checked G.94–G.95, G.103–G.105 (p. 015003-17), G.200 with $n_{\mathrm{r}}(\omega)$ (p. 015003-45) and the polariton analogy (p. 015003-23).
5. **Step 5 cooling benchmark (major): accepted.** It is split into (a) the electrostatic tensor limit, which depends on the self-energy topic, and (b) retardation tests at $Q_\parallel\lesssim$ a few $\omega/c$. P.D now states where the transverse part matters. Tielrooij et al. was added but is marked unverified (see 14).
6. **Diagonal E.410f (minor): accepted.** Reproduced 1459.7/812.6 vs 1462.8/807.1 cm⁻¹ at 30° and 1523.9/799.2 vs 1526.5/794.2 cm⁻¹ at 45°. The $2\times2$ matrix gives the exact values. (E.410f is (B7.21) in v2.)
7. **P.C numbers (minor): accepted.** Reproduced $4.2\times10^{-5}$ (static), $4.1\times10^{-2}$ (undamped, 170 meV), $1.7\times10^{-2}$ ($\gamma=5/4$ cm⁻¹) and pole shifts of −0.028 and −0.004 cm⁻¹.
8. **Primer status (minor): accepted, with one correction.** §B.9 (v1 lines 229–250, `inv3`/`sectors`) and §B.1 (the kernel $K(\hat{\mathbf{q}})$ of v1) were confirmed. The suggested formula had the wrong sign: §B.9 defines $\frac{\omega^2}{4\pi c}D^{\mathrm{e}}=(x^2\mathcal{P}^{\mathrm{T}}-\boldsymbol{\varepsilon})^{-1}$ ((B9.3) in v2), so $D^{\mathrm{e}}=+\frac{4\pi c}{\omega^2}(x^2\mathcal{P}^{\mathrm{T}}-\boldsymbol{\varepsilon})^{-1}$, which also agrees with P.1. Step 1 now reuses that code (1–3 days).
9. **Lorentz field (minor): accepted, refined.** EPP p. 32 (E.151 and the "all other crystal symmetries" sentence) was confirmed. Added that $\mathbf{Z}^*_1\approx-\mathbf{Z}^*_2$ is exact by the acoustic sum rule for hBN's two sublattices but fails for the low-symmetry targets. Added a pitfall line.
10. **§8 description (minor): accepted.** EPP p. 100 (E.455–E.459 with a $z$-resolved Dyson equation for $D^{\mathrm{e}}$) and §8.1 (E.460–E.472) were checked.
11. **Gauge invariance (minor): accepted.** E.341–E.350 are continuity relations for the longitudinal projections. The paragraph was rewritten.
12. **Sum rules (minor): accepted.** The corrections section of primer v1 (in v2, §B.1) attributes the spurious acoustic dipoles to the phase-less vertex; both conditions are now stated.
13. **Parameter attribution (minor): accepted, refined.** Attribution to [7] is marked unverified. "Rigid shift" is not quite right: each oblique mode moves by an amount between the shifts of its band edges. With the review's alternative set, the computed shifts are −8 to +4 cm⁻¹ (upper) and −8 to −20 cm⁻¹ (lower), so the header now gives this range.
14. **Local fields and nonlocality citations (minor): accepted in substance.** The local-field remark is now [est], and both papers are cited. Their bibliographic data could not be re-checked in this revision (web-search budget exhausted, publisher and Crossref hosts blocked). They are therefore listed as unverified, together with Tielrooij et al., not as verified references.
15. **Reststrahlen interpolation (minor): accepted.** $\varepsilon_{xx}\varepsilon_{zz}=-12.5$ at 1500 cm⁻¹ was confirmed, and the bullet was added.
16. **Step 4 (minor): accepted.** Fitted damping first; ab initio linewidths as a separate follow-up.
