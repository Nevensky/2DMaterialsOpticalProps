# Anisotropic crystals and polaritons beyond the zone centre

Labels: **[est]** established; **[derived]** algebra done for this note from EPP's equations and checked numerically; **[proposal]** a suggested EPP extension; **[spec]** speculation. Units are Gaussian (EPP §§4–8), and $Z^*$ carries the charge. E.410a–h are the primer's proposed EPP §7.1.1 ([§Z.3](../../sections/Z03-corrections.html)). hBN parameters (standard set, cf. [7]; primer preset; provenance not re-verified): in-plane (t) TO/LO $1370/1610$ cm⁻¹ with $\epsilon^\infty_t=4.87$, and out-of-plane (z) $780/830$ cm⁻¹ with $\epsilon^\infty_z=2.95$.

## 1. What is missing in EPP

- **The formalism is tensorial.** E.380 is $3\times3$ at $\bG=\bG'=0$, E.385 (with E p.86a′) carries full $Z^*$ tensors, E.431–E.432 define an electronic dielectric tensor, and E.445–E.454 contain $\bq\cdot\hat\epsilon^{\rm e}\cdot\bq$. The direction-dependent LO–TO term therefore exists for $c\to\infty$; E.454a gives its matrix form ([§B.9](../../sections/B09-electrostatic.html)).
- **What is evaluated is isotropic.** §4.2.1 (E.144–E.165, pp. 28–35) assumes diagonal $Z^*$ and force constants (p. 29), treats one Cartesian direction at a time and assumes $\nabla\times\bE=0$. §7.1 evaluates only $\bq=0$ without electrons (E.406–E.410), where $D^0$ is isotropic (E.408).
- **What is missing at finite $c$:**
  - L–T mixing in $D$;
  - the extraordinary (Fresnel) dispersion and the hyperbolic regime $\epsilon_t\epsilon_z<0$;
  - slabs, surfaces and heterostructures (§8 is one electrostatic $\delta(z)$ sheet with in-plane modes only);
  - damping beyond $i\eta$;
  - low-symmetry crystals;
  - a protocol for first-principles inputs.
- **Primer status.** [§B.6](../../sections/B06-propagator.html) decouples $D_L$ and $D_T$ for isotropic media only, and notes that in hBN $\hat\bq\cdot\Pi_{jj}\cdot\hat{\mathbf t}\propto(\Pi_\parallel-\Pi_\perp)\sin\theta\cos\theta$. [§B.7](../../sections/B07-phonon-se.html) derives the isotropic branches and treats hBN as isotropic. Anisotropy enters only formally, through the determinant E.410d ([§Z.1](../../sections/Z01-outlook.html)).

## 2. Why it matters

- **Hyperbolic bands [est].** hBN has $\epsilon_z<0<\epsilon_t$ at 780–830 cm⁻¹ (Type I) and $\epsilon_t<0<\epsilon_z$ at 1370–1610 cm⁻¹ (Type II).
  - In hBN, volume-confined modes reach $\lambda/86$ with $Q\le283$ [7], the film dispersion tunes with thickness down to a few layers [6], and isotopic purification triples the lifetime [11].
  - In α-MoO₃, the in-plane polaritons are up to 60 times shorter than light and have amplitude lifetimes of about 8 ps [12].
  - In monoclinic β-Ga₂O₃, shear polaritons have frequency-dependent axes [15].
- **hBN numbers [derived].**
  - $\epsilon^0_t=6.73$ and $\epsilon^0_z=3.34$ (Lyddane–Sachs–Teller).
  - $\omega_{\rm LO}^2-\omega_{\rm TO}^2=4\pi Z^{*2}/(\Omega_{\rm uc}\mu\epsilon^\infty)$ implies $Z^*_t=2.80|e|$ and $Z^*_z=0.73|e|$ ($\Omega_{\rm uc}=18.1$ Å³ per BN pair, $\mu=6.10$ u).
  - The light cone is at $q_c=1.9\times10^{-4}$ Å⁻¹ in-plane and $0.84\times10^{-4}$ Å⁻¹ out-of-plane.
  - For $q\gg q_c$ the directional LO modes are, with $\theta=\angle(\bq,c)$:

  | $\theta$ | 0° | 15° | 30° | 45° | 60° | 75° | 90° |
  |---|---|---|---|---|---|---|---|
  | upper (cm⁻¹) | 1370 (T) | 1399 | 1463 | 1527 | 1574 | 1601 | 1610 |
  | lower (cm⁻¹) | 830 | 822 | 807 | 794 | 786 | 781 | 780 (T) |

- **Slab regime [derived].** Take a free-standing 10 nm film.
  - At 1500 cm⁻¹, $q_0=4.4\times10^{-3}$ Å⁻¹ $=46\,k_0$ ($\lambda_p\approx144$ nm), and the rays run at 38° to the plane.
  - At 1400 cm⁻¹, $q_0\approx6k_0\approx k_0\sqrt{|\epsilon_t|}$. Quasi-statics fails there, and the retarded $D$ is needed.
  - At 1600 cm⁻¹, $q_0=0.063$ and $q_2=0.28$ Å⁻¹, so phonon dispersion (nonlocality) matters.
- **For 2D GW/BSE practitioners.** hBN is the standard encapsulant.
  - The effective environment $\sqrt{\epsilon_t\epsilon_z}$ rises from 3.79 ($\epsilon^\infty$) to 4.74 ($\epsilon^0$), about 25% [derived], across the 97–200 meV phonon window. Which of the two values $W$ should see is a frequency question, and EPP's total $D$ answers it.
  - Hyperbolic hBN modes cool hot graphene carriers in picoseconds [10].
  - The monolayer limit is the 2D LO dispersion without LO–TO splitting ([9], [§B.10](../../sections/B10-2d.html)).

## 3. What standard theory does

- **Giustino.** (G.180) is "the generalization of Fröhlich's model to the case of anisotropic crystalline lattices and multiple phonon modes", with $\mathbf Z_\kappa$ and $\boldsymbol\epsilon^\infty$ tensors and the cutoff $e^{-a|\bq+\bG|^2}$ (p. 015003-35). For the phonon side, §VI.A.3 cites Gonze–Lee. $W=W_{\rm e}+W_{\rm ph}$ (G.104–G.105) is instantaneous, and light enters only as a classical field with an input refractive index (G.200). There are no polaritons.
- **DFPT [4, 5] [est].** DFPT supplies $Z^*_\kappa$, $\epsilon^\infty$ and the analytic force constants. The nonanalytic term is $C^{\rm NA}_{\kappa\alpha,\kappa'\beta}=\frac{4\pi}{\Omega}\frac{(\hat\bq\cdot Z^*_\kappa)_\alpha(\hat\bq\cdot Z^*_{\kappa'})_\beta}{\hat\bq\cdot\epsilon^\infty\cdot\hat\bq}$, and $\epsilon(\omega)=\epsilon^\infty+\frac{4\pi}{\Omega}\sum_mS_m/(\omega_m^2-\omega^2)$ with mode oscillator strengths $S_m$. The quasi-LO/TO angular dispersion in uniaxial crystals is classic Raman physics (Loudon 1964, unverified).
- **Polaritons in practice [est].** A two-step route: an ab initio or fitted $\epsilon(\omega)$ goes into Maxwell's equations (Fresnel, transfer matrices, slab formulas [6, 7, 13]). Classical roots: [1, 2]. Royo–Stengel [14] give exact long-range 2D electrostatics, with $2\times2$ dielectric matrices and the out-of-plane dynamical-charge coupling, but no retardation.

## 4. How it could enter the EPP formalism

**P.A Closed-form propagator [proposal; derived].** At $\bG=\bG'=0$, E.380 inverts algebraically. With $(D^0)^{-1}=-\frac{1}{4\pi c}[\omega^2\mathbb 1-c^2(q^2\mathbb 1-\bq\bq)]$ from E.375,

$$D(\bq,\omega)=-4\pi c\,\Lambda^{-1},\qquad \Lambda=\omega^2\hat\epsilon(\bq,\omega)-c^2(q^2\mathbb 1-\bq\bq),\qquad \hat\epsilon=\mathbb 1+\frac{4\pi c}{\omega^2}\big(\Pi^{\rm e}_{jj}+\Pi^{\rm ion}_{jj}\big),\tag{P.1}$$

and E.385 with E p.86a′ gives

$$\hat\epsilon^{\rm ion}_{\mu\nu}(\bq,\omega)=\frac{4\pi}{\Omega_{\rm uc}}\sum_{\sigma}\frac{j_{\mu\sigma}(\bq)\,j^*_{\nu\sigma}(\bq)}{\omega^2_\sigma(\bq)-\omega^2-i\eta\omega}.\tag{P.2}$$

P.2 is Gonze–Lee's $\epsilon(\omega)$ with $S_\sigma=\bj_\sigma\bj_\sigma^\dagger$, extended to finite $\bq$. $\hat\epsilon$ is the material tensor; EPP's E.432 tensor contains $D^0$ and agrees with it only in the longitudinal projection. E.410a is P.1 with $\hat\epsilon=\epsilon_\infty\mathbb 1$.

**P.B Pole condition [derived].** Sylvester's identity gives $\det(\text{E.410d})=\prod_\sigma(\omega^2-\omega_\sigma^2)\,\det\Lambda[\hat\epsilon]/\det\Lambda[\hat\epsilon^{\rm e}]$, so the modes solve $\det\Lambda=0$, the Fresnel equation. For a uniaxial crystal the ordinary wave has $c^2q^2=\omega^2\epsilon_t$, and the extraordinary wave has

$$\frac{c^2q^2}{\omega^2}=\frac{\epsilon_t(\omega)\,\epsilon_z(\omega)}{\epsilon_t(\omega)\sin^2\theta+\epsilon_z(\omega)\cos^2\theta}.\tag{P.3}$$

- **$\bq\to0$ at finite $c$:** $\det\hat\epsilon=0$, so every polar mode sits at its LO frequency, isotropically. This is the anisotropic E.410e (hBN: 1610 cm⁻¹ twice and 830 cm⁻¹).
- **$q\gg q_c$:** $\epsilon_{LL}\equiv\hat\bq\cdot\hat\epsilon\cdot\hat\bq=0$, which is E.410f and E.454 with a tensor (the table in §2).

At large $q$, bulk hyperbolic polaritons are these oblique LO modes ($\nabla\cdot\epsilon\nabla\phi=0$).

**P.C L–T mixing in $W$ [derived].** In the basis $(\hat\bq,\hat{\mathbf t}_1)$, with $\hat{\mathbf t}_1$ in the $\bq$–$c$ plane, $\epsilon_{LT}=(\epsilon_t-\epsilon_z)\sin\theta\cos\theta$ and $\epsilon_{TT}=\epsilon_t\cos^2\theta+\epsilon_z\sin^2\theta$. With $n=cq/\omega$ and $W\equiv-\frac{\omega^2}{cq^2}\hat\bq\cdot D\cdot\hat\bq$ (§B.6),

$$W(\bq,\omega)=\frac{4\pi}{q^2}\,\frac{\epsilon_{TT}-n^2}{\epsilon_{LL}(\epsilon_{TT}-n^2)-\epsilon_{LT}^2}\;\xrightarrow{n\to\infty}\;\frac{4\pi}{\bq\cdot\hat\epsilon\cdot\bq}\ (\text{E.445}).\tag{P.4}$$

The correction $\approx\epsilon_{LT}^2/(\epsilon_{LL}n^2)$ is $4\times10^{-5}$ for static hBN at $\theta=45^\circ$, $q=0.01$ Å⁻¹ and $\hbar\omega=170$ meV. It is negligible for transport and of order one for $q\lesssim q_c$.

**P.D Phonons and vertices [proposal].** E.398–E.399 with $D^{\rm e}=-4\pi c\Lambda_{\rm e}^{-1}$. The off-diagonal $\Omega^{\rm e}_{\sigma\sigma'}$ mixes E₁u and A₂u for oblique $\hat\bq$, and its $c\to\infty$ limit is E.454a. The retarded vertex E.416/E.430 with the full $D^{\rm e}$, transverse part kept, couples electrons in an adjacent layer to the hyperbolic modes.

**P.E Slabs and heterostructures [proposal].** Generalize the $z$-resolved Dyson equation of EPP p. 100 to the whole field:

$$D(\bQ,\omega;z,z')=D^0+\!\int\! dz_1dz_2\,D^0(z,z_1)\big[\Pi^{\rm e}_{jj}+\Pi^{\rm ion}_{jj}\big](z_1,z_2)\,D(z_2,z'),\qquad \Pi^{\rm ion}_{jj}=\sum_{\ell\ell'}\tfrac{\Omega_{\rm uc}}{S_{\rm uc}}\Pi^{\rm ion}_{\ell\ell'}(\bQ,\omega)\,\delta(z-z_\ell)\delta(z'-z_{\ell'})\tag{P.5}$$

- The ionic part is a layer-resolved E.455; $\ell\neq\ell'$ terms come from interlayer force constants.
- $D^0$ is the $q_z$ transform of E.375, $\propto e^{-\kappa_0|z-z'|}/2\kappa_0$ with $\kappa_0=\sqrt{Q^2-\omega^2/c^2}$, plus a contact term.
- The continuum quasi-static limit gives the slab modes of [6] [derived]. For a film of thickness $d$ between media $\epsilon_1$ and $\epsilon_3$,

$$q_l=\frac{1}{\kappa d}\Big[\pi l+\arctan\frac{\epsilon_1}{\epsilon_z\kappa}+\arctan\frac{\epsilon_3}{\epsilon_z\kappa}\Big],\qquad \kappa=\sqrt{-\epsilon_t/\epsilon_z},\quad q_l>0.\tag{P.6}$$

**Gauge invariance.** Every $\Pi$ must have the form $\Pi^{\rm para}(\omega)-\Pi^{\rm para}(0)$ (E.383, E.385). Only then is $\hat\epsilon^{\rm ion}(\omega\to0)$ finite, and only then does E.350 hold for every tensor component.

## 5. Pitfalls

- **Phonon double counting.** $\omega_\sigma(\bq)$ is the DFPT dynamical matrix minus only its $\bG=0$ nonanalytic term, because E.380 re-adds only $\bG=0$. The $\bG\neq0$ dipolar part stays in, so the bare $\omega_\sigma(0)$ is the zero-field TO. Phonons interpolated with Gonze–Lee would count LO–TO twice. The subtraction must use the same $Z^*$ and $\epsilon^\infty$ as the force constants.
- **Electron double counting.** Use either the model $\epsilon_\infty$ (E.410a) or the RPA $\Pi^{\rm e}$, not both.
  - RPA without local fields or $f^{xc}$ is not DFPT's $\epsilon^\infty$, and a 10% error in $\epsilon^\infty_t$ moves $\omega_{{\rm LO},t}$ by 2.8 meV (22 cm⁻¹), 9% of the Type II band [derived].
  - [spec] Local fields across the vdW gap matter most for $\epsilon_z$, so a $\bG=0$ RPA may misplace the Type I band.
- **Coulomb double counting.** $D$ already contains $v$ (E.364, E.374). In a BSE kernel, drop the $\bG=0$ exchange.
- **Sum rules.** Enforce the tensor acoustic sum rule $\sum_aZ^*_{a,\mu\nu}=0$; otherwise the acoustic modes acquire spurious dipoles (§Z.3).
- **Convergence and regularization.**
  - Keep $\hat\bq$ in the head of $\Pi^{\rm e}$.
  - E.454a needs an Ewald cutoff and removal of the self term.
  - Anharmonic and isotopic damping [11] belong in $F$ (E.398), not in $\eta$ or $\Pi^{\rm e}$.
- **Approximations.** E.404 fails for mixed modes; use E.410d or P.B instead. P.6 is local; for $q\gtrsim0.05$ Å⁻¹ keep $\omega_\sigma(\bq)$ and $\bj_\sigma(\bq)$ in P.2.

## 6. State of the art (verified)

1. K. Huang, "On the interaction between the radiation field and ionic crystals", Proc. R. Soc. Lond. A **208**, 352 (1951), doi:10.1098/rspa.1951.0166. First phonon-polariton theory.
2. R. Fuchs, K. L. Kliewer, "Optical modes of vibration in an ionic crystal slab", Phys. Rev. **140**, A2076 (1965), doi:10.1103/PhysRev.140.A2076. Slab and surface modes.
3. R. Geick, C. H. Perry, G. Rupprecht, "Normal modes in hexagonal boron nitride", Phys. Rev. **146**, 543 (1966), doi:10.1103/PhysRev.146.543. Polarized IR modes of hBN.
4. X. Gonze, C. Lee, "Dynamical matrices, Born effective charges, dielectric permittivity tensors, and interatomic force constants from density-functional perturbation theory", Phys. Rev. B **55**, 10355 (1997), doi:10.1103/PhysRevB.55.10355. Nonanalytic term, $S_m$, $\epsilon(\omega)$.
5. S. Baroni, S. de Gironcoli, A. Dal Corso, P. Giannozzi, "Phonons and related crystal properties from density-functional perturbation theory", Rev. Mod. Phys. **73**, 515 (2001), doi:10.1103/RevModPhys.73.515. DFPT $Z^*$ and $\epsilon^\infty$.
6. S. Dai *et al.*, "Tunable phonon polaritons in atomically thin van der Waals crystals of boron nitride", Science **343**, 1125 (2014), doi:10.1126/science.1246833. Thickness scaling.
7. J. D. Caldwell *et al.*, "Sub-diffractional volume-confined polaritons in the natural hyperbolic material hexagonal boron nitride", Nat. Commun. **5**, 5221 (2014), doi:10.1038/ncomms6221. Type I and II modes.
8. D. N. Basov, M. M. Fogler, F. J. García de Abajo, "Polaritons in van der Waals materials", Science **354**, aag1992 (2016), doi:10.1126/science.aag1992. Review.
9. T. Sohier, M. Gibertini, M. Calandra, F. Mauri, N. Marzari, "Breakdown of optical phonons' splitting in two-dimensional materials", Nano Lett. **17**, 3758 (2017), doi:10.1021/acs.nanolett.7b01090. 2D limit.
10. A. Principi *et al.*, "Super-Planckian electron cooling in a van der Waals stack", Phys. Rev. Lett. **118**, 126804 (2017), doi:10.1103/PhysRevLett.118.126804. Electrons coupled to hyperbolic modes.
11. A. J. Giles *et al.*, "Ultralow-loss polaritons in isotopically pure boron nitride", Nat. Mater. **17**, 134 (2018), doi:10.1038/nmat5047. Damping.
12. W. Ma *et al.*, "In-plane anisotropic and ultra-low-loss polaritons in a natural van der Waals crystal", Nature **562**, 557 (2018), doi:10.1038/s41586-018-0618-9. Biaxial α-MoO₃.
13. G. Álvarez-Pérez *et al.*, "Infrared permittivity of the biaxial van der Waals semiconductor α-MoO₃ from near- and far-field correlative studies", Adv. Mater. **32**, 1908176 (2020), doi:10.1002/adma.201908176. Benchmark tensor.
14. M. Royo, M. Stengel, "Exact long-range dielectric screening and interatomic force constants in quasi-two-dimensional crystals", Phys. Rev. X **11**, 041027 (2021), doi:10.1103/PhysRevX.11.041027. 2D electrostatic benchmark.
15. N. C. Passler *et al.*, "Hyperbolic shear polaritons in low-symmetry crystals", Nature **602**, 595 (2022), doi:10.1038/s41586-021-04328-y. Non-co-diagonal oscillators.

Unverified:
- R. Loudon, Adv. Phys. **13**, 423 (1964). The page and DOI were not confirmed.
- N. Rivera, T. Christensen, P. Narang, "Phonon polaritonics in two-dimensional materials", arXiv:1901.05308. The journal version was not confirmed.

## 7. First steps

1. **Bulk hBN model (2–4 days).** Implement P.1–P.3 and reproduce the §2 table and the $q\to0$ limit. Check E.410d against $\det\Lambda$, then plot the isofrequency surfaces across $q_c$.
2. **First-principles inputs (1–2 weeks).**
   - DFPT: $Z^*$, $\epsilon^\infty$ and the zero-field Γ modes of bulk hBN.
   - Your RPA code: the head of $\Pi^{\rm e}$ along $a$, $c$ and 45°, without and with local fields, compared with DFPT.
   - Propagate the differences to the LO frequencies with P.3.
3. **Slab (3–6 weeks).** Solve P.5 for $N$-layer hBN with the retarded $D^0$.
   - Benchmark against P.6 and against the data of [6].
   - Test $N=1$ against EPP §8, [9] and [14].
   - Map where retardation ($q\lesssim k_0\sqrt{|\epsilon|}$) and nonlocality ($q\gtrsim0.05$ Å⁻¹) matter.
4. **Low symmetry (1–3 months).** α-MoO₃ [12, 13], then β-Ga₂O₃ [15]. Both need an anharmonic phonon self-energy.
5. **Electrons (months).** Build a retarded anisotropic Fröhlich vertex for graphene on hBN, test it against the cooling rate of [10], and put the dynamical hBN $W$ into the exciton screening.

## 8. Open questions

- Is a $\bG=0$ RPA adequate for $\epsilon_z$ of a layered crystal, or are $\bG\neq0$ local fields needed before the Type I band can be trusted?
- Interface Born charges differ from bulk ones. Should P.5 use layer-dependent $Z^*_\ell$ from slab DFPT?
- How do dynamical quadrupoles [14] enter at finite $c$, at order $O(q)$?
- [spec] Do exciton polaritons and the hyperbolic continuum hybridize measurably in hBN/TMD stacks?
