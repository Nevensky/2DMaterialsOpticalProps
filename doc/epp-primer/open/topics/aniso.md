# Anisotropic crystals and polaritons beyond the zone centre

Labels: **[est]** established; **[derived]** algebra done for this note from EPP's equations and checked numerically; **[proposal]** a suggested EPP extension; **[spec]** speculation. Units are Gaussian (EPP §§4–8), and $Z^*$ carries the charge. E.410a–h are the primer's proposed EPP §7.1.1 ([§Z.3](../../sections/Z03-corrections.html)). hBN parameters (a widely used rounded set, also the primer preset; its attribution to [7] is unverified, and published fits differ by 5–20 cm⁻¹): in-plane (t) TO/LO $1370/1610$ cm⁻¹ with $\epsilon^\infty_t=4.87$, and out-of-plane (z) $780/830$ cm⁻¹ with $\epsilon^\infty_z=2.95$. Sensitivity [derived]: each entry of the §2 table lies between its band's TO and LO, so it moves by an amount between the shifts of the two band edges. With the alternative fit 4.90/2.95, 1360/1614 and 760/825 cm⁻¹ (reported by the review; source not re-verified), the upper row moves by −8 to +4 cm⁻¹ and the lower row by −8 to −20 cm⁻¹.

## 1. What is missing in EPP

- **The formalism is tensorial.** E.380 is $3\times3$ at $\bG=\bG'=0$, E.385 (with E p.86a′) carries full $Z^*$ tensors, E.431–E.432 define an electronic dielectric tensor, and E.445–E.454 contain $\bq\cdot\hat\epsilon^{\rm e}\cdot\bq$. The direction-dependent LO–TO term therefore exists for $c\to\infty$; E.454a gives its matrix form ([§B.9](../../sections/B09-electrostatic.html)).
- **What is evaluated is isotropic.** §4.2.1 (E.144–E.165, pp. 28–35) assumes diagonal $Z^*$ and force constants (p. 29), treats one Cartesian direction at a time and assumes $\nabla\times\bE=0$. §7.1 evaluates only $\bq=0$ without electrons (E.406–E.410), where $D^0$ is isotropic (E.408).
  - §4.2.1 also adds the cubic Lorentz field $\frac{4\pi}{3}\mathbf P$ (E.151, p. 32) and states that "the same expression is valid for all other crystal symmetries". For a uniaxial, layered crystal such as hBN the Lorentz factor is a tensor, not $\frac13\mathbb 1$.
  - In the first-principles route the local field is already inside $Z^*$ and the analytic force constants (the $\bG\neq0$ dipole sum), so E.151 must not be added on top ([§B.3](../../sections/B03-dielectric.html); [§B.7](../../sections/B07-phonon-se.html): "the Dyson route has no Lorentz local field").
  - It also sets $Z^*_{\rm a}\approx-Z^*_{\rm b}$ (p. 32). The acoustic sum rule makes this exact for two inequivalent sublattices, as in hBN, but not for α-MoO₃ or β-Ga₂O₃.
- **What is missing at finite $c$:**
  - L–T mixing in $D$;
  - the extraordinary (Fresnel) dispersion and the hyperbolic regime $\epsilon_t\epsilon_z<0$;
  - slabs, surfaces and heterostructures. §8 has one $\delta(z)$ ionic sheet (E.455) with a formally retarded, $z$-resolved $D^{\rm e}$ (p. 100). It evaluates only the electrostatic, in-plane limit (§8.1, E.460–E.472: $\alpha=x,y$, with the contact term discarded);
  - damping beyond $i\eta$;
  - low-symmetry crystals;
  - a protocol for first-principles inputs.
- **Primer status.** [§B.6](../../sections/B06-propagator.html) decouples $D_L$ and $D_T$ for isotropic media only, and notes that in hBN $\hat\bq\cdot\Pi_{jj}\cdot\hat{\mathbf t}\propto(\Pi_\parallel-\Pi_\perp)\sin\theta\cos\theta$. [§B.7](../../sections/B07-phonon-se.html) derives the isotropic branches and treats hBN as isotropic.
  - [§B.9](../../sections/B09-electrostatic.html) ("Which sector survives $c\to\infty$?") derives $D^{\rm e}=\frac{4\pi c}{\omega^2}(x^2P_T-\epsilon)^{-1}$, with $x=cq/\omega$, for a constant uniaxial $\epsilon$. It plots the $ll$, $lt$ and $tt$ sectors at angle $\theta$ and shows that $lt$ and $tt$ fall as $x^{-2}$.
  - [§B.1](../../sections/B01-foundations.html) plots the anisotropic $\bG=0$ kernel $K(\hat\bq)$.
  - What is missing is $\Pi^{\rm ion}$ inside the tensor. The polariton and hyperbolic dispersion therefore appear only through the determinant E.410d ([§Z.1](../../sections/Z01-outlook.html)).

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
  - At 1400 cm⁻¹, $q_0=5.4\times10^{-4}$ Å⁻¹ $=6.1\,k_0$; that this is $\approx k_0\sqrt{|\epsilon_t|}$ is a coincidence. The exact retarded TM solution ($k_z^2=\epsilon_t(k_0^2-q^2/\epsilon_z)$, $\kappa_{\rm out}=\sqrt{q^2-k_0^2}$, $k_zd=\pi m-2\arctan(|\epsilon_t|\kappa_{\rm out}/k_z)$) shifts it by only +1.3% (+0.03% at 1500 cm⁻¹).
  - The shift grows to 3% at 1390 cm⁻¹ and 13% at 1380 cm⁻¹, where $q_0\approx2k_0$. For a 10 nm film, retardation matters only within ~10–20 cm⁻¹ of $\omega_{{\rm TO},t}$. Since $q_0\propto1/d$, the window is wider for thicker films.
  - At 1600 cm⁻¹, $q_0=0.063$ and $q_2=0.28$ Å⁻¹, so phonon dispersion (nonlocality) matters. Gubbin and De Liberato analyse this nonlocal phonon response for polar films of ~10 nm (Phys. Rev. X 2020; unverified, §6).
- **For 2D GW/BSE practitioners.** hBN is the standard encapsulant.
  - The effective environment $\sqrt{\epsilon_t\epsilon_z}$ rises from 3.79 ($\epsilon^\infty$) to 4.74 ($\epsilon^0$), about 25% [derived], across the 97–200 meV phonon window. Which of the two values $W$ should see is a frequency question, and EPP's total $D$ answers it.
  - Between the two limits the screening is not an interpolation. Inside 97–103 and 170–200 meV, $\epsilon_t\epsilon_z<0$ (−12.5 at 1500 cm⁻¹), so $W(\bQ,\omega)$ has a hyperbolic-polariton continuum instead of a real $\epsilon_{\rm eff}$ [derived]. This is where a frequency-dependent $W$ matters most.
  - Hyperbolic hBN modes cool hot graphene carriers in picoseconds [10].
  - The monolayer limit is the 2D LO dispersion without LO–TO splitting ([9], [§B.10](../../sections/B10-2d.html)).

## 3. What standard theory does

- **Giustino.** (G.180) is "the generalization of Fröhlich's model to the case of anisotropic crystalline lattices and multiple phonon modes", with $\mathbf Z_\kappa$ and $\boldsymbol\epsilon^\infty$ tensors and the cutoff $e^{-a|\bq+\bG|^2}$ (p. 015003-35). For the phonon side, §VI.A.3 cites Gonze–Lee.
  - $W=W_{\rm e}+W_{\rm ph}$ (G.104–G.105) is dynamically screened. $W_{\rm ph}(12)=\int W_{\rm e}(13)D(34)W_{\rm e}(24)$ contains the nuclear density–density correlation function $D$ (G.103), so $W_{\rm ph}(\omega)$ carries the phonon poles, and $W_{\rm e}(\omega)$ is frequency dependent too (G.94–G.95).
  - It is built on the instantaneous Coulomb $v$ ($c\to\infty$), so it contains no transverse photons and no retardation.
  - Light enters only as a classical field with an input refractive index $n_{\rm r}(\omega)$ (G.200, p. 015003-45). Polaritons appear only as the analogy for coupled plasmon–phonon modes, "the electronic analog of photon polaritons" (p. 015003-23). There is no polariton theory.
- **DFPT [4, 5] [est].** DFPT supplies $Z^*_\kappa$, $\epsilon^\infty$ and the analytic force constants. The nonanalytic term is $C^{\rm NA}_{\kappa\alpha,\kappa'\beta}=\frac{4\pi}{\Omega}\frac{(\hat\bq\cdot Z^*_\kappa)_\alpha(\hat\bq\cdot Z^*_{\kappa'})_\beta}{\hat\bq\cdot\epsilon^\infty\cdot\hat\bq}$, and $\epsilon(\omega)=\epsilon^\infty+\frac{4\pi}{\Omega}\sum_mS_m/(\omega_m^2-\omega^2)$ with mode oscillator strengths $S_m$. The quasi-LO/TO angular dispersion in uniaxial crystals is classic Raman physics (Loudon 1964, unverified).
- **Polaritons in practice [est].** A two-step route: an ab initio or fitted $\epsilon(\omega)$ goes into Maxwell's equations (Fresnel, transfer matrices, slab formulas [6, 7, 13]). Classical roots: [1, 2]. Royo–Stengel [14] give exact long-range 2D electrostatics, with $2\times2$ dielectric matrices and the out-of-plane dynamical-charge coupling, but no retardation.

## 4. How it could enter the EPP formalism

**P.A Closed-form propagator [proposal; derived].** At $\bG=\bG'=0$, E.380 inverts algebraically. With $(D^0)^{-1}=-\frac{1}{4\pi c}[\omega^2\mathbb 1-c^2(q^2\mathbb 1-\bq\bq)]$ from E.375,

$$D(\bq,\omega)=-4\pi c\,\Lambda^{-1},\qquad \Lambda=\omega^2\hat\epsilon(\bq,\omega)-c^2(q^2\mathbb 1-\bq\bq),\qquad \hat\epsilon=\mathbb 1+\frac{4\pi c}{\omega^2}\big(\Pi^{\rm e}_{jj}+\Pi^{\rm ion}_{jj}\big),\tag{P.1}$$

and E.385 with E p.86a′ gives

$$\hat\epsilon^{\rm ion}_{\mu\nu}(\bq,\omega)=\frac{4\pi}{\Omega_{\rm uc}}\sum_{\sigma}\frac{j_{\mu\sigma}(\bq)\,j^*_{\nu\sigma}(\bq)}{\omega^2_\sigma(\bq)-\omega^2-i\eta\omega}.\tag{P.2}$$

P.2 is Gonze–Lee's $\epsilon(\omega)$ with $S_\sigma=\bj_\sigma\bj_\sigma^\dagger$, extended to finite $\bq$. $\hat\epsilon$ is the material tensor; EPP's E.432 tensor contains $D^0$ and agrees with it only in the longitudinal projection. E.410a is P.1 with $\hat\epsilon=\epsilon_\infty\mathbb 1$. With $\Pi^{\rm ion}=0$ and a constant uniaxial $\hat\epsilon^{\rm e}$, P.1 is the closed form of [§B.9](../../sections/B09-electrostatic.html).

**P.B Pole condition [derived].** Sylvester's identity gives $\det(\text{E.410d})=\prod_\sigma(\omega^2-\omega_\sigma^2)\,\det\Lambda[\hat\epsilon]/\det\Lambda[\hat\epsilon^{\rm e}]$, so the modes solve $\det\Lambda=0$, the Fresnel equation. For a uniaxial crystal the ordinary wave has $c^2q^2=\omega^2\epsilon_t$, and the extraordinary wave has

$$\frac{c^2q^2}{\omega^2}=\frac{\epsilon_t(\omega)\,\epsilon_z(\omega)}{\epsilon_t(\omega)\sin^2\theta+\epsilon_z(\omega)\cos^2\theta}.\tag{P.3}$$

- **$\bq\to0$ at finite $c$:** $\det\hat\epsilon=0$, so every polar mode sits at its LO frequency, isotropically. This is the anisotropic E.410e (hBN: 1610 cm⁻¹ twice and 830 cm⁻¹).
- **$q\gg q_c$:** $\epsilon_{LL}\equiv\hat\bq\cdot\hat\epsilon\cdot\hat\bq=0$. This is the $c\to\infty$ limit of E.410d, i.e. the $\bG=0$ term of E.454a projected on the modes; it gives the table in §2.
  - For oblique $\hat\bq$ both E₁u and A₂u have $\hat\bq\cdot\bj_\sigma\neq0$ and mix. Diagonalizing the $2\times2$ matrix reproduces the table [derived].
  - The diagonal E.410f/E.454 with a tensor $\epsilon^\infty$ misplaces the oblique modes. At $\theta=30^\circ$ it gives 1459.7/812.6 cm⁻¹ instead of 1462.8/807.1 cm⁻¹, an error of up to 3 cm⁻¹ (upper) and 5.5 cm⁻¹ (lower, 11% of the Type I band) [derived].

At large $q$, bulk hyperbolic polaritons are these oblique LO modes ($\nabla\cdot\epsilon\nabla\phi=0$).

**P.C L–T mixing in $W$ [derived].** In the basis $(\hat\bq,\hat{\mathbf t}_1)$, with $\hat{\mathbf t}_1$ in the $\bq$–$c$ plane, $\epsilon_{LT}=(\epsilon_t-\epsilon_z)\sin\theta\cos\theta$ and $\epsilon_{TT}=\epsilon_t\cos^2\theta+\epsilon_z\sin^2\theta$. With $n=cq/\omega$ and $W\equiv-\frac{\omega^2}{cq^2}\hat\bq\cdot D\cdot\hat\bq$ (§B.6),

$$W(\bq,\omega)=\frac{4\pi}{q^2}\,\frac{\epsilon_{TT}-n^2}{\epsilon_{LL}(\epsilon_{TT}-n^2)-\epsilon_{LT}^2}\;\xrightarrow{n\to\infty}\;\frac{4\pi}{\bq\cdot\hat\epsilon\cdot\bq}\ (\text{E.445}).\tag{P.4}$$

- At $q=0.01$ Å⁻¹ and $\theta=45^\circ$, the retarded L–T mixing moves the oblique LO poles by −0.03 cm⁻¹ (at 1526.5 cm⁻¹) and −0.004 cm⁻¹ (at 794.2 cm⁻¹).
- Away from resonances the relative correction $\epsilon_{LT}^2/(\epsilon_{LL}n^2)$ is ~10⁻⁵–10⁻⁴; with the static $\epsilon^0$, valid for $\hbar\omega\ll97$ meV, it is $4\times10^{-5}$. Next to a TO pole (170 meV) it reaches $4\times10^{-2}$ undamped and $2\times10^{-2}$ with $\gamma\approx5$ cm⁻¹.
- It is negligible for transport and of order one for $q\lesssim q_c$.

**P.D Phonons and vertices [proposal].** E.398–E.399 with $D^{\rm e}=-4\pi c\Lambda_{\rm e}^{-1}$. The off-diagonal $\Omega^{\rm e}_{\sigma\sigma'}$ mixes E₁u and A₂u for oblique $\hat\bq$, and its $c\to\infty$ limit is E.454a. The retarded vertex E.416/E.430 with the full $D^{\rm e}$ couples electrons in an adjacent layer to the hyperbolic modes. Its transverse part matters only for $q\lesssim$ a few $k_0\sqrt{|\epsilon|}$; at the $q\sim k_F$ that dominate transport and cooling, the $lt$ and $tt$ sectors are suppressed by $(q_c/q)^2$ (P.C; [§B.9](../../sections/B09-electrostatic.html)).

**P.E Slabs and heterostructures [proposal].** Generalize the $z$-resolved Dyson equation of EPP p. 100 to the whole field:

$$D(\bQ,\omega;z,z')=D^0+\!\int\! dz_1dz_2\,D^0(z,z_1)\big[\Pi^{\rm e}_{jj}+\Pi^{\rm ion}_{jj}\big](z_1,z_2)\,D(z_2,z'),\qquad \Pi^{\rm ion}_{jj}=\sum_{\ell\ell'}\tfrac{\Omega_{\rm uc}}{S_{\rm uc}}\Pi^{\rm ion}_{\ell\ell'}(\bQ,\omega)\,f_\ell(z)f_{\ell'}(z')\tag{P.5}$$

- $f_\ell(z)$ is a normalized layer profile ($\int f_\ell\,dz=1$) of width ~3.3 Å, centred on layer $\ell$. Point sheets $\delta(z-z_\ell)$ are not allowed (next bullet).
- **Contact term.** $D^0$ is the $q_z$ transform of E.375. It is $\propto e^{-\kappa_0|z-z'|}/2\kappa_0$ with $\kappa_0=\sqrt{Q^2-\omega^2/c^2}$, and its $zz$ component also contains the contact term $-\frac{4\pi c}{\omega^2}\delta(z-z')$.
  - With point sheets, the $\ell=\ell'$ $zz$ terms would give $\delta(0)$. EPP discards exactly this "divergent local term" before E.463 and keeps $\alpha=x,y$.
  - Here the term must be kept. It is the depolarization field that sets the A₂u frequency and the whole Type I band ($\epsilon_z<0$).
  - Regularize it with the profiles $f_\ell$, or remove the on-layer self term Ewald-style and add the contact term back analytically.
  - Require that the $N\to\infty$ continuum limit reproduce $\epsilon_z(\omega)$ and the 780–830 cm⁻¹ band.
- The ionic part is a layer-resolved E.455; $\ell\neq\ell'$ terms come from interlayer force constants.
- The continuum quasi-static limit gives the slab modes of [6] [derived]. For a film of thickness $d$ between media $\epsilon_1$ and $\epsilon_3$,

$$q_l=\frac{1}{\kappa d}\Big[\pi l+\arctan\frac{\epsilon_1}{\epsilon_z\kappa}+\arctan\frac{\epsilon_3}{\epsilon_z\kappa}\Big],\qquad \kappa=\sqrt{-\epsilon_t/\epsilon_z},\quad q_l>0.\tag{P.6}$$

**Gauge invariance.** Every $\Pi$ must have the form $\Pi^{\rm para}(\omega)-\Pi^{\rm para}(0)$ (E.383, E.385).
- This removes the spurious $1/\omega^2$ Drude term, so $\hat\epsilon^{\rm ion}(\omega\to0)$ is finite.
- It makes the static longitudinal response vanish (a static longitudinal $\mathbf A$ is pure gauge), consistent with continuity. E.341–E.350 constrain only the longitudinal projections.
- Applied componentwise at finite $\bq$, it also drops the $O(q^2)$ static orbital-magnetic response, which is negligible for phonon polaritons.

## 5. Pitfalls

- **Phonon double counting.**
  - Bulk: $\omega_\sigma(\bq)$ is the DFPT dynamical matrix minus only its $\bG=0$ nonanalytic term, because E.380 re-adds only $\bG=0$. The $\bG\neq0$ dipolar part stays in, so the bare $\omega_\sigma(0)$ is the zero-field TO. Phonons interpolated with Gonze–Lee would count LO–TO twice. The subtraction must use the same $Z^*$ and $\epsilon^\infty$ as the force constants.
  - Slab: compute the bare force constants with DFPT and a 2D Coulomb truncation. Then remove the layer-resolved, $\bQ$-dependent 2D long-range term (in the form of Sohier [9] or Royo–Stengel [14], with the same $Z^*_\ell$ and electronic screening) that P.5 adds back. The bulk $\bG=0$ rule does not apply, and a periodic-image supercell gives wrong $\ell\ell'$ dipole couplings.
- **Electron double counting.** Use either the model $\epsilon_\infty$ (E.410a) or the RPA $\Pi^{\rm e}$, not both.
  - The $\Pi^{\rm e}$ that enters $D$ at $\bG=\bG'=0$ must be proper with respect to the macroscopic Coulomb line: $\epsilon_M=1-v_0\bar\chi_{00}$, with $\bar\chi$ computed with the $\bG=0$ Hartree term removed. Without local fields, $\bar\chi_{00}=\chi^0_{00}$. Inserting the full RPA head $\chi_{00}$ counts the $\bG=0$ Coulomb twice; this is the RPA analogue of the BSE remark below.
  - A density-response code gives only $\hat\bq\cdot\Pi^{\rm e}\cdot\hat\bq$ (through E.350). At $q\to0$ in an insulator the tensor follows from heads along several directions. Transverse components at finite $q$ need a current–current bubble.
  - RPA without local fields or $f^{xc}$ is not DFPT's $\epsilon^\infty$, and a 10% error in $\epsilon^\infty_t$ moves $\omega_{{\rm LO},t}$ by 2.8 meV (22 cm⁻¹), 9% of the Type II band [derived].
  - [est] Local fields across the vdW gap matter most for $\epsilon_z$, so a $\bG=0$ RPA may misplace the Type I band. For graphite, TDDFT reproduces the anisotropic loss spectra for $\bq$ off the basal planes only with local fields included (Marinopoulos et al. 2002; unverified, §6).
- **Coulomb double counting.** $D$ already contains $v$ (E.364, E.374). In a BSE kernel, drop the $\bG=0$ exchange.
- **Local-field double counting.** Do not add the Lorentz field of E.151 to DFPT inputs. It is cubic-only, and $Z^*$ and the analytic force constants already contain the local field (§1).
- **Sum rules.** Use the phase-corrected vertex E p.86a′ and enforce $\sum_{\rm a}Z^*_{{\rm a},\mu\nu}=0$. Without either, the acoustic modes acquire spurious, cell-dependent dipoles ([§Z.3](../../sections/Z03-corrections.html), E.131b).
- **Convergence and regularization.**
  - Keep $\hat\bq$ in the head of $\Pi^{\rm e}$.
  - E.454a needs an Ewald cutoff and removal of the self term. P.5 needs the regularized $zz$ contact term (§4 P.E).
  - Anharmonic and isotopic damping [11] belong in $F$ (E.398), not in $\eta$ or $\Pi^{\rm e}$.
- **Approximations.** E.404 fails for mixed modes; use E.410d or P.B instead (§4 P.B quantifies the error). P.6 is local; for $q\gtrsim0.05$ Å⁻¹ keep $\omega_\sigma(\bq)$ and $\bj_\sigma(\bq)$ in P.2.

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

1. **Bulk hBN model (1–3 days).** Start from the `inv3`/`sectors` code of [§B.9](../../sections/B09-electrostatic.html) and replace $\epsilon$ by $\hat\epsilon^{\rm e}+\hat\epsilon^{\rm ion}(\omega)$ (P.2).
   - Implement P.1–P.3 and reproduce the §2 table (from the $2\times2$ mode matrix, not the diagonal E.410f) and the $q\to0$ limit.
   - Check E.410d against $\det\Lambda$, then plot the isofrequency surfaces across $q_c$.
2. **First-principles inputs (1–2 weeks).**
   - DFPT: $Z^*$, $\epsilon^\infty$ and the zero-field Γ modes of bulk hBN.
   - Your RPA code: take $\epsilon_M(\hat\bq)=1-v_0\bar\chi_{00}(\hat\bq)$, with the $\bG=0$ Hartree term removed, with and without $\bG\neq0$ local fields. Do this along $a$ and $c$, and along 45° as a consistency check (uniaxial: it must equal $(\epsilon_t+\epsilon_z)/2$).
   - Reconstruct the $q\to0$ tensor and insert it as $\frac{4\pi c}{\omega^2}\Pi^{\rm e}=\hat\epsilon_M-1$. Compare with DFPT. For low-symmetry crystals use at least 4–6 directions, including one oblique in the $ac$ plane for $\epsilon_{xz}$ in the monoclinic case.
   - The transverse components at finite $q$ cannot come from $\chi_{\rho\rho}$. They need a current–current implementation with $\Pi^{\rm para}(\omega)-\Pi^{\rm para}(0)$.
   - Propagate the differences to the LO frequencies with P.3.
3. **Slab (4–8 weeks).** Solve P.5 for $N$-layer hBN with the retarded $D^0$.
   - Regularize the $zz$ contact term (profiles $f_\ell$ or Ewald self-term removal) and check the $N\to\infty$ limit against $\epsilon_z(\omega)$ and the Type I band. Take the bare force constants from 2D-truncated DFPT minus the layer-resolved long-range term (§5). Together these add about 1–2 weeks.
   - Benchmark against P.6 and against the data of [6].
   - Test $N=1$ against EPP §8, [9] and [14].
   - Map where retardation ($q\lesssim3$–$5\,k_0$; check against the exact TM root of §2) and nonlocality ($q\gtrsim0.05$ Å⁻¹) matter.
4. **Low symmetry (1–3 months).** α-MoO₃ [12, 13], then β-Ga₂O₃ [15].
   - Start with fitted $\gamma_\sigma$ in $F$ (E.398) and benchmark against [13]. The dispersion, and the frequency-dependent axis rotation of shear polaritons [15], already follow from P.2 with non-co-diagonal $\bj_\sigma\bj_\sigma^\dagger$ and any damping.
   - Ab initio anharmonic and isotope linewidths are a separate, longer follow-up. It overlaps the damping-beyond-$i\eta$ gap.
5. **Electrons (months).**
   - (a) Electrostatic tensor limit. Use the $z$-resolved P.5 in the $c\to\infty$ limit (interlayer $W$ including the hBN $\Pi^{\rm ion}$). Benchmark the graphene→hBN cooling power against [10] and the measurement of Tielrooij et al. (unverified, §6). This needs the electron energy-loss rate, i.e. an electron self-energy with $D$, from the self-energy topic.
   - (b) Retardation can be tested only where $q\lesssim$ a few $k_0$: far-field emission or IR absorption of the stack, and polariton launching.
   - Put the dynamical hBN $W$ into the exciton screening.

## 8. Open questions

- How large are the $\bG\neq0$ local-field corrections to $\epsilon^\infty_z$ of hBN? Do they move the Type I band by more than the 5–20 cm⁻¹ spread of published fits? For graphite, local fields are known to be essential for $\bq$ off the basal plane (Marinopoulos et al. 2002; unverified, §6).
- Interface Born charges differ from bulk ones. Should P.5 use layer-dependent $Z^*_\ell$ from slab DFPT?
- How do dynamical quadrupoles [14] enter at finite $c$, at order $O(q)$?
- [spec] Do exciton polaritons and the hyperbolic continuum hybridize measurably in hBN/TMD stacks?

## Review responses

Every reviewer number was re-derived with a short script (exact retarded TM slab root, $2\times2$ mode matrix, extraordinary condition at $q=0.01$ Å⁻¹). The EPP and Giustino page and equation references were checked in the PDFs, and the primer claims in §B.1, §B.3, §B.7, §B.9, §Z.1 and §Z.3.

1. **Slab retardation (major): accepted.** At 1400 cm⁻¹ the exact correction is +1.34% (+0.025% at 1500, 3.2% at 1390, 13.2% at 1380 cm⁻¹, where $q\approx2.1k_0$). The §2 bullet and §7 step 3 were rewritten, and "quasi-statics fails" was removed.
2. **P.5 contact term and slab double counting (major): accepted.** EPP p. 101 drops the "divergent local term $\sim\delta(z-z')$" and keeps $\mu,\nu=x,y$ (E.463–E.464). P.5 now uses profiles $f_\ell(z)$ and keeps the $zz$ term, with an $N\to\infty$ check. §5 has a slab rule, and step 3 is now 4–8 weeks.
3. **RPA head (major): accepted.** Step 2 and §5 now use $\epsilon_M=1-v_0\bar\chi_{00}$, the reconstruction of the tensor from several directions, and the current–current requirement at finite $q$. Wording refined: without local fields, $\bar\chi_{00}$ reduces to $\chi^0_{00}$. What must not be used is the full reducible head $\chi_{00}$.
4. **Giustino's $W$ (major): accepted.** Checked G.94–G.95, G.103–G.105 (p. 015003-17), G.200 with $n_{\rm r}(\omega)$ (p. 015003-45) and the polariton analogy (p. 015003-23).
5. **Step 5 cooling benchmark (major): accepted.** It is split into (a) the electrostatic tensor limit, which depends on the self-energy topic, and (b) retardation tests at $q\lesssim$ a few $k_0$. P.D now states where the transverse part matters. Tielrooij et al. was added but is marked unverified (see 14).
6. **Diagonal E.410f (minor): accepted.** Reproduced 1459.7/812.6 vs 1462.8/807.1 cm⁻¹ at 30° and 1523.9/799.2 vs 1526.5/794.2 cm⁻¹ at 45°. The $2\times2$ matrix gives the exact values.
7. **P.C numbers (minor): accepted.** Reproduced $4.2\times10^{-5}$ (static), $4.1\times10^{-2}$ (undamped, 170 meV), $1.7\times10^{-2}$ ($\gamma=5/4$ cm⁻¹) and pole shifts of −0.028 and −0.004 cm⁻¹.
8. **Primer status (minor): accepted, with one correction.** §B.9 (lines 229–250, `inv3`/`sectors`) and §B.1 ($K(\hat\bq)$) were confirmed. The suggested formula had the wrong sign: §B.9 defines $d^{\rm e}=\frac{\omega^2}{4\pi c}D^{\rm e}=(x^2P_T-\epsilon)^{-1}$, so $D^{\rm e}=+\frac{4\pi c}{\omega^2}(x^2P_T-\epsilon)^{-1}$, which also agrees with P.1. Step 1 now reuses that code (1–3 days).
9. **Lorentz field (minor): accepted, refined.** EPP p. 32 (E.151 and the "all other crystal symmetries" sentence) was confirmed. Added that $Z^*_{\rm a}\approx-Z^*_{\rm b}$ is exact by the acoustic sum rule for hBN's two sublattices but fails for the low-symmetry targets. Added a pitfall line.
10. **§8 description (minor): accepted.** EPP p. 100 (E.455–E.459 with a $z$-resolved Dyson equation for $D^{\rm e}$) and §8.1 (E.460–E.472) were checked.
11. **Gauge invariance (minor): accepted.** E.341–E.350 are continuity relations for the longitudinal projections. The paragraph was rewritten.
12. **Sum rules (minor): accepted.** §Z.3 attributes the spurious acoustic dipoles to the phase-less vertex; both conditions are now stated.
13. **Parameter attribution (minor): accepted, refined.** Attribution to [7] is marked unverified. "Rigid shift" is not quite right: each oblique mode moves by an amount between the shifts of its band edges. With the review's alternative set, the computed shifts are −8 to +4 cm⁻¹ (upper) and −8 to −20 cm⁻¹ (lower), so the header now gives this range.
14. **Local fields and nonlocality citations (minor): accepted in substance.** The local-field remark is now [est], and both papers are cited. Their bibliographic data could not be re-checked in this revision (web-search budget exhausted, publisher and Crossref hosts blocked). They are therefore listed as unverified, together with Tielrooij et al., not as verified references.
15. **Reststrahlen interpolation (minor): accepted.** $\epsilon_t\epsilon_z=-12.5$ at 1500 cm⁻¹ was confirmed, and the bullet was added.
16. **Step 4 (minor): accepted.** Fitted damping first; ab initio linewidths as a separate follow-up.
