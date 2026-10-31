# Primer v2 · Notation and conventions

This file fixes one notation for the whole of primer v2. Every section uses it in the main text.
The symbols of the two sources appear only in footnotes and in the source map at the end of each section:
F. Giustino, Rev. Mod. Phys. **89**, 015003 (2017), cited as **G**, and the EPP manuscript, cited as **E**.
Primer v1 (`../sections/`) keeps the source notations and stays as the reference.

Section 9 maps the source symbols. Sections 10 and 11 set the writing rules that make a section self-contained.

---

## 1. Typography

| Rule | Example |
|---|---|
| Running indices and variables are italic. | $\kappa$, $\alpha$, $\nu$, $n$, $\bk$, $\omega$, $t$ |
| Labels that abbreviate a word are upright, the **whole** word. | $\Pi^{\mathrm{ion}}$, $\Pi^{\mathrm{ph}}$, $\varepsilon^{\mathrm e}$, $\Omega_{\mathrm V}^{\mathrm{SC}}$, $M_{\mathrm r}$, $k_{\mathrm B}$, $\omega_{\mathrm{LO}}$, $V^{\mathrm{KS}}$ |
| Euler's number, the imaginary unit and the differential are upright. | $\mathrm e^{\mathrm i\bq\cdot\bR_\ell}$, $\mathrm d^3r$ |
| Vectors and Cartesian tensors are bold; their components are italic with indices. | $\bu_{\ell\kappa}$ and $u_{\ell\kappa\alpha}$; $\boldsymbol\varepsilon$ and $\varepsilon_{\alpha\beta}$; $\mathbf Z^*_\kappa$ and $Z^*_{\kappa,\alpha\beta}$ |
| A hat marks the unit vector of a vector; a vector defined as a unit vector (a bond direction $\bn$) needs none; otherwise hats mark operators where ambiguous. | $\hat{\bq}=\bq/q$ |
| Operators carry hats only where an operator and a number of the same name meet (A.1 quantization, expectation values). | $\hat u_{\ell\kappa\alpha}$, $\langle\hat u^2\rangle$ |
| A comma separates index groups that belong to different objects. | $e_{\kappa\alpha,\nu}$ (row $\kappa\alpha$, column $\nu$); $\Phi_{\ell\kappa\alpha,\ell'\kappa'\beta}$ |
| A tilde means mass-scaled. | $\tilde\Phi=M^{-1/2}\Phi M^{-1/2}$ |
| Complex conjugate $^*$, Hermitian conjugate $^\dagger$, transpose $^{\mathsf T}$; time ordering $\mathcal T$. | |

Reason for the upright $\mathrm e$: the exponential sits next to the eigenvector components $e_{\kappa\alpha,\nu}$
and the elementary charge $e$ in many equations, for example $\mathrm e^{\mathrm i\bq\cdot\bR_\ell}e_{\kappa\alpha,\nu}(\bq)$.

MathJax macros (v2 head): `\ee` = $\mathrm e$, `\ii` = $\mathrm i$, `\dd` = $\mathrm d$, `\btau` = $\boldsymbol\tau$,
`\ba` = $\mathbf a$, `\bb` = $\mathbf b$, `\bF` = $\mathbf F$, `\bn` = $\mathbf n$,
`\OV` = $\Omega_{\mathrm V}$, `\OS` = $\Omega_{\mathrm S}$, `\OVSC` = $\Omega_{\mathrm V}^{\mathrm{SC}}$,
`\OSSC` = $\Omega_{\mathrm S}^{\mathrm{SC}}$, `\OBZ` = $\Omega_{\mathrm{BZ}}$, `\rme` = $\mathrm e$ (label),
`\ion` = $\mathrm{ion}$, `\ph` = $\mathrm{ph}$, `\SC` = $\mathrm{SC}$, `\kB` = $k_{\mathrm B}$.
Counts such as $N_\ell$ and $N_\kappa$ are written in full. The v1 macros (`\bk`, `\bq`, `\br`, `\bR`, `\bG`, `\bA`, `\bE`,
`\bB`, `\bj`, `\bJ`, `\bP`, `\bQ`, `\bu`, `\be`, `\bp`, `\bZ`, `\Tr`, `\ket`, `\bra`, `\braket`, `\mel`, `\hc`) remain.

---

## 2. Index alphabet

| Index | Meaning | Range |
|---|---|---|
| $\ell,\ell'$ | unit cell; also the integer triple $(\ell_1,\ell_2,\ell_3)$ | $N_\ell$ cells in the Born–von Kármán supercell |
| $\kappa,\kappa'$ | atom of the basis | $1,\dots,N_\kappa$ |
| $\alpha,\beta$ | Cartesian component, free (outer) indices | $x,y,z$ |
| $\gamma,\delta$ | Cartesian component, summed (inner) indices | $x,y,z$ |
| $\nu,\nu'$ | phonon branch | $1,\dots,3N_\kappa$ |
| $n,m$ | electron band | |
| $\lambda,\lambda'$ | photon polarization | $1,2$ transverse, $3$ longitudinal |
| $j,j'$ | primitive-vector label | $1,2,3$ |

Relative to a surface or a 2D sheet, $z$ is the normal direction and $x,y$ lie in the plane:
$\varepsilon_{zz}$ is out of plane, $\varepsilon_{xx},\varepsilon_{yy}$ are in plane.
$\mu$ is reserved for the chemical potential and is never a Cartesian index. $i$ is never an index.

---

## 3. Crystal, supercell, reciprocal space

| Symbol | Meaning | Note |
|---|---|---|
| $\ba_j$, $a$ | primitive lattice vectors; lattice constant | |
| $\bR_\ell=\ell_1\ba_1+\ell_2\ba_2+\ell_3\ba_3$ | lattice vector of cell $\ell$ | $\bR_0=\mathbf 0$ |
| $\btau_\kappa$ | equilibrium position of atom $\kappa$ inside the cell | |
| $\btau_{\ell\kappa}=\bR_\ell+\btau_\kappa$ | equilibrium position of atom $\kappa$ of cell $\ell$ | |
| $\bu_{\ell\kappa}$, $u_{\ell\kappa\alpha}$ | displacement from $\btau_{\ell\kappa}$ and its components | instantaneous position $\btau_{\ell\kappa}+\bu_{\ell\kappa}$ |
| $\br$ | electron position | |
| $N_1,N_2,N_3$; $N_\ell=N_1N_2N_3$ | supercell size along $\ba_j$; number of cells | BvK: $\bu$ repeats after $N_j$ cells along $\ba_j$ |
| $N_\kappa$ | atoms per cell | $3N_\kappa$ branches, $3N_\kappa N_\ell$ coordinates |
| $N_\bq$, $N_\bk$ | number of phonon and electron wavevectors | $N_\bq=N_\ell$ |
| $\bb_j$ | reciprocal primitive vectors, $\ba_j\cdot\bb_{j'}=2\pi\delta_{jj'}$ | |
| $\bG$ | reciprocal lattice vector | |
| $\bq$, $\bk$ | phonon (or photon) and electron wavevectors | on the grid $\bq=\sum_j (s_j/N_j)\,\bb_j$ |
| $\Omega_{\mathrm V}$, $\Omega_{\mathrm S}$ | volume of the unit cell (3D); area of the unit cell (2D) | |
| $\Omega_{\mathrm V}^{\mathrm{SC}}=N_\ell\Omega_{\mathrm V}$, $\Omega_{\mathrm S}^{\mathrm{SC}}=N_\ell\Omega_{\mathrm S}$ | supercell volume and area = the whole crystal under BvK | |
| $\Omega_{\mathrm{BZ}}=(2\pi)^3/\Omega_{\mathrm V}$ | Brillouin-zone volume; $(2\pi)^2/\Omega_{\mathrm S}$ in 2D | |

Lattice sums: $\sum_{\ell}\mathrm e^{\mathrm i\bq\cdot\bR_\ell}=N_\ell\sum_{\bG}\delta_{\bq,\bG}$ and
$\sum_{\bq\in\mathrm{BZ}}\mathrm e^{\mathrm i\bq\cdot\bR_\ell}=N_\ell\,\delta_{\ell,0}$.
Grid sums become integrals as $\frac1{N_\bq}\sum_\bq\to\frac{\Omega_{\mathrm V}}{(2\pi)^3}\int_{\mathrm{BZ}}\mathrm d^3q$.

**2D.** In-plane wavevectors are capitals with $\parallel$, out-of-plane components lowercase:
$\bq+\bG\to(\bQ_\parallel+\bG_\parallel,\ q_z+G_z)$; likewise $\mathbf K_\parallel$ for electrons.
In-plane positions are $\br_\parallel$, the normal coordinate is $z$.

---

## 4. Lattice dynamics

| Symbol | Meaning | Units |
|---|---|---|
| $M_\kappa$ | mass of atom $\kappa$ | amu |
| $M_{\mathrm c}=\sum_\kappa M_\kappa$, $M_{\mathrm r}$ | cell mass; reduced mass of a two-atom cell, $1/M_{\mathrm r}=1/M_1+1/M_2$ | amu |
| $U$, $U_0$ | Born–Oppenheimer energy surface (electrons in their ground state, nuclei clamped); its minimum | eV |
| $\bF_{\ell\kappa}$, $F_{\ell\kappa\alpha}$ | force on atom $\ell\kappa$ | eV/Å |
| $\Phi_{\ell\kappa\alpha,\ell'\kappa'\beta}=\partial^2U/\partial u_{\ell\kappa\alpha}\partial u_{\ell'\kappa'\beta}$ | interatomic force constants | eV/Å² |
| $\Phi_{\kappa\alpha,\kappa'\beta}(\bq)=\sum_{\ell}\Phi_{0\kappa\alpha,\ell\kappa'\beta}\,\mathrm e^{\mathrm i\bq\cdot\bR_\ell}$ | their lattice Fourier transform | eV/Å² |
| $\tilde\Phi_{\kappa\alpha,\kappa'\beta}(\bq)=\Phi_{\kappa\alpha,\kappa'\beta}(\bq)/\sqrt{M_\kappa M_{\kappa'}}$ | dynamical matrix (mass-scaled) | eV Å⁻² amu⁻¹ (= $\omega^2$) |
| $e_{\kappa\alpha,\nu}(\bq)$; $\be_{\kappa\nu}(\bq)$ | eigenvector of $\tilde\Phi(\bq)$: component; Cartesian 3-vector of atom $\kappa$ | dimensionless |
| $\omega_{\bq\nu}$ | phonon frequency, $\omega^2_{\bq\nu}$ eigenvalue of $\tilde\Phi(\bq)$ | $\hbar\omega$ in meV |
| $\xi$, $\xi_\kappa$, $\boldsymbol\xi_\kappa$ | complex amplitude of a lattice wave (physical displacement amplitude of atom $\kappa$), $\boldsymbol\xi_\kappa\propto\be_{\kappa\nu}/\sqrt{M_\kappa}$ for a normalized eigenvector; in phasor pictures the complex displacement $u_\ell=\xi\,\mathrm e^{\mathrm i(q\ell a-\omega t)}$, physical part $\mathrm{Re}\,u_\ell$ | Å |
| $K$, $K_1$, $K_2$; $K_{\mathrm{eff}}$ | spring stiffness in model crystals; effective stiffness felt by an atom in a wave | eV/Å² |
| $\bd$ | a common displacement of every atom (rigid translation) | Å |
| $\Phi^{\mathrm{sr}}$, $\Phi^{\mathrm{dd}}$ | short-range (analytic) and long-range dipole–dipole parts of $\Phi$ | |
| $Q_{\bq\nu}$, $P_{\bq\nu}$ | normal coordinate and its conjugate momentum | $\sqrt{\text{mass}}\times$length |
| $a_{\bq\nu}$, $a^\dagger_{\bq\nu}$ | phonon annihilation and creation operators | |
| $\phi_{\bq\nu}=a_{\bq\nu}+a^\dagger_{-\bq\nu}$ | phonon field (dimensionless displacement of mode $\bq\nu$) | |
| $n_{\bq\nu}=n_{\mathrm B}(\omega_{\bq\nu})$, $n_{\mathrm B}(\omega)=1/(\mathrm e^{\hbar\omega/k_{\mathrm B}T}-1)$ | Bose–Einstein occupation | |
| $(\hbar/2M_\kappa\omega_{\bq\nu})^{1/2}$ | zero-point amplitude of atom $\kappa$ in mode $\bq\nu$, always written out | Å |
| LA, TA, LO, TO, ZA, ZO | branch labels, upright | |

Key relations, in this sign convention:

$$\sum_{\kappa'\beta}\tilde\Phi_{\kappa\alpha,\kappa'\beta}(\bq)\,e_{\kappa'\beta,\nu}(\bq)=\omega^2_{\bq\nu}\,e_{\kappa\alpha,\nu}(\bq),\qquad
\sum_{\ell'\kappa'}\Phi_{\ell\kappa\alpha,\ell'\kappa'\beta}=0\ \ \text{(acoustic sum rule)},$$

$$\hat u_{\ell\kappa\alpha}=\sum_{\bq\nu}\left(\frac{\hbar}{2N_\ell M_\kappa\,\omega_{\bq\nu}}\right)^{1/2}e_{\kappa\alpha,\nu}(\bq)\,\mathrm e^{\mathrm i\bq\cdot\bR_\ell}\,\hat\phi_{\bq\nu},\qquad
e_{\kappa\alpha,\nu}(-\bq)=e^*_{\kappa\alpha,\nu}(\bq).$$

The phase $\mathrm e^{+\mathrm i\bq\cdot\bR_\ell}$ in $\Phi(\bq)$ belongs to the cell of the second atom. This is what a displacement pattern
$\propto\mathrm e^{+\mathrm i\bq\cdot\bR_\ell}$ requires. E Eq. (15) uses the opposite assignment.

---

## 5. Electrons

| Symbol | Meaning |
|---|---|
| $E_{n\bk}$ | band energy (G writes $\varepsilon_{n\bk}$; $\varepsilon$ is kept for the dielectric function) |
| $\psi_{n\bk}(\br)=\mathrm e^{\mathrm i\bk\cdot\br}u_{n\bk}(\br)/\sqrt{N_\ell}$ | Bloch state normalized on the supercell; $u_{n\bk}$ normalized on the cell (scalar, band labels: not a displacement) |
| $c_{n\bk}$, $c^\dagger_{n\bk}$ | electron operators; spin is suppressed and a factor 2 is written where it enters |
| $f_{n\bk}$, $\mu$, $E_{\mathrm F}$ | Fermi–Dirac occupation, chemical potential, Fermi energy |
| $m_{\mathrm e}$, $-e$ | electron mass and charge, with $e>0$ |
| $V^{\mathrm{KS}}$, $V^{\mathrm{H}}$, $V^{\mathrm{xc}}$ | Kohn–Sham, Hartree and exchange–correlation potential energies |
| $g_{mn\nu}(\bk,\bq)$ | electron–phonon matrix element for $n\bk\to m\,\bk+\bq$ with absorption of $\bq\nu$ (or emission of $-\bq\nu$) |
| $g^{\mathrm{DW}}_{mn\nu\nu'}(\bk,\bq,\bq')$ | second-order (Debye–Waller) matrix element |

The electron–phonon Hamiltonian is
$\hat H^{\mathrm{e\text{-}ph}}=N_\ell^{-1/2}\sum_{\bk\bq mn\nu}g_{mn\nu}(\bk,\bq)\,c^\dagger_{m\bk+\bq}c_{n\bk}\,\phi_{\bq\nu}$
(the EPW and Giustino convention; E puts $N^{-1/2}$ inside $g$ and shifts the labels, see Section 9).

---

## 6. Fields, charges, response, propagators

| Symbol | Meaning |
|---|---|
| $\bE$, $\bB$, $\bA$, $A_0$ | electric field, magnetic field, vector potential, scalar potential; temporal gauge is $A_0=0$ |
| $n$, $\rho$, $\bj$ | number density, charge density ($\rho=-e\,n_{\mathrm e}+\rho^{\mathrm{ion}}$), current density |
| $\bP$ | polarization density |
| $Z_\kappa$; $\mathbf Z^*_\kappa$, $Z^*_{\kappa,\alpha\beta}$ | ionic valence; Born effective charge tensor, **dimensionless** (the charge is $eZ^*$) |
| $\varepsilon$, $\boldsymbol\varepsilon$, $\varepsilon_{\alpha\beta}$ | dielectric function: scalar, tensor, components |
| $\varepsilon^{\mathrm e}$ | electrons only, ions clamped; no superscript = total (electrons + ions) |
| $\varepsilon_\infty$, $\varepsilon_{\mathrm s}$ | electronic (clamped-ion) dielectric constant, i.e. $\varepsilon^{\mathrm e}$ well above the phonon and below the gap frequencies; total static dielectric constant |
| $\varepsilon_{\mathrm M}$, $\varepsilon_{\bG\bG'}$ | macroscopic and microscopic dielectric functions |
| $\sigma_{\alpha\beta}$ | conductivity |
| $c$, $\alpha_{\mathrm{fs}}$ | speed of light, fine-structure constant |
| $G$ | electron propagator |
| $F$, $F_{\nu\nu'}(\bq,\omega)$ | phonon propagator, built from $\phi_{\bq\nu}$ (a force always carries atom labels: $F_{\ell\kappa\alpha}$) |
| $D$, $D_{\alpha\beta}$, $D_{\lambda\lambda'}$ | photon propagator; $D^0=D^{0\mathrm L}+D^{0\mathrm T}$ (free, longitudinal + transverse); $D^{\mathrm e}$ dressed by electrons only |
| $v$, $v(\bq)=4\pi e^2/q^2$ | bare Coulomb interaction (Gaussian, per unit volume) |
| $W$, $W^{\mathrm e}$, $W^{\mathrm{ph}}$ | screened interaction: total, electronic, phonon part |
| $\Sigma$ | electron self-energy |
| $\Pi_{\alpha\beta}$; $\Pi^{\mathrm e}$, $\Pi^{\mathrm{ion}}$; $\Pi^{\mathrm{para}}$, $\Pi^{\mathrm{dia}}$ | photon self-energy (Cartesian indices): electronic and ionic parts; paramagnetic and diamagnetic parts |
| $\Pi^{\mathrm{ph}}_{\nu\nu'}(\bq,\omega)$ | phonon self-energy (branch indices) |
| $\be_\lambda(\bq)$ | photon polarization basis, $\be_3=\hat\bq$; $D_{\alpha\beta}=\sum_{\lambda\lambda'}e_{\lambda\alpha}D_{\lambda\lambda'}e_{\lambda'\beta}$ |
| $\Xi$ | reserved for the BSE kernel |
| $\alpha^2F(\omega)$, $\lambda_{\mathrm{ep}}$ | Eliashberg function (one traditional name, not a product with $\alpha$ or $F$); electron–phonon coupling strength |

Current vertices of electrons and ions are written $\mathbf J^{\mathrm e}$ and $\mathbf J^{\mathrm{ion}}$; the exact form is fixed in Part B together with $\Pi^{\mathrm e}$.

---

## 7. Units

- **Gaussian units with $e$, $\hbar$, $m_{\mathrm e}$ and $c$ written explicitly in every equation.** Atomic units follow by setting
  $e=\hbar=m_{\mathrm e}=1$ and $c=1/\alpha_{\mathrm{fs}}\approx137.036$.
- **SI footnotes:** every key equation that contains $e$, $c$ or a Coulomb kernel carries a footnote with its SI form.
  Rules: $e^2\to e^2/4\pi\varepsilon_0$; $4\pi$ in Poisson's and Maxwell's equations $\to1/\varepsilon_0$; $\bA_{\mathrm{Gauss}}=c\,\bA_{\mathrm{SI}}$;
  minimal coupling $\mathbf p+\frac ec\bA\to\mathbf p+e\bA$ for an electron of charge $-e$.
  Here $\varepsilon_0$ is the vacuum permittivity, used only in SI footnotes.
- **Numbers** in eV, meV, Å, amu, K. Conversion box:

| Quantity | Value |
|---|---|
| $e^2$ (Gaussian) $=e^2/4\pi\varepsilon_0$ (SI) | 14.3996 eV Å |
| $\hbar c$ | 1973.27 eV Å |
| $\hbar^2/2m_{\mathrm e}$ | 3.80998 eV Å² |
| $\hbar^2/(\mathrm{amu}\cdot\text{Å}^2)$ | 4.1802 meV, so $\hbar\omega=64.654\ \text{meV}\times\sqrt{K[\text{eV/Å}^2]/M[\text{amu}]}$ |
| $k_{\mathrm B}$ | 0.0861733 meV/K |
| 1 meV | 8.06554 cm⁻¹ = 0.241799 THz |
| $\hbar/(2M\omega)$ | 2.0901 Å² / ($M$[amu] · $\hbar\omega$[meV]) |

---

## 8. Fourier and sign conventions

- **Time:** $f(\omega)=\int\mathrm dt\,\mathrm e^{\mathrm i\omega t}f(t)$, $f(t)=\int\frac{\mathrm d\omega}{2\pi}\mathrm e^{-\mathrm i\omega t}f(\omega)$.
  Time dependence $\mathrm e^{-\mathrm i\omega t}$; retarded functions use $\omega+\mathrm i\eta$; damped oscillators $\omega^2-\omega_0^2+\mathrm i\gamma\omega$.
- **Lattice:** displacement patterns $\propto\mathrm e^{+\mathrm i\bq\cdot\bR_\ell}$; force constants transformed with the phase on the second atom's cell (Section 4).
- **Normal coordinates:** symmetric $N_\ell^{-1/2}$ in both directions,
  $u_{\ell\kappa\alpha}=(N_\ell M_\kappa)^{-1/2}\sum_{\bq\nu}\mathrm e^{\mathrm i\bq\cdot\bR_\ell}e_{\kappa\alpha,\nu}(\bq)\,Q_{\bq\nu}$.
- **Continuum fields in a crystal:** kernels per unit volume, so sums over states carry $1/\Omega_{\mathrm V}^{\mathrm{SC}}$
  (for example $\chi^0(\bq,\omega)=\frac{2}{\Omega_{\mathrm V}^{\mathrm{SC}}}\sum\dots$). The normalization of the field expansion of $\bA$ is fixed in Part B.

---

## 9. Source map (for footnotes)

| Quantity | v2 | Giustino (G) | EPP (E) |
|---|---|---|---|
| cell index; lattice vector | $\ell$; $\bR_\ell$ | $p$; $\bR_p$ | $i$; $\bR^0_i$ |
| number of cells | $N_\ell=N_1N_2N_3$ | $N_p=N_1N_2N_3$ | $N$ |
| atom index; atoms per cell | $\kappa$; $N_\kappa$ | $\kappa$; $M$ (a count!) | $a,b$; $s$ |
| basis position; equilibrium site | $\btau_\kappa$; $\btau_{\ell\kappa}$ | $\btau_\kappa$; $\btau^0_{\kappa p}$ ($\btau_{\kappa p}=\bR_p+\btau_\kappa$ is the moving position) | $\mathbf r^0_a$; $\bR^0_{ia}$ |
| Cartesian index | $\alpha,\beta$ (outer), $\gamma,\delta$ (inner) | $\alpha,\alpha'$ | $\mu,\nu$ (§§1–4); $\alpha,\beta,\mu,\nu$ from p. 61 (§5.2.1) on |
| displacement | $u_{\ell\kappa\alpha}$ | $\Delta\tau_{\kappa\alpha p}$ | $Q_{ia\mu}$ |
| mass | $M_\kappa$ | $M_\kappa$; reference mass $M_0$ | $m_a$ |
| force constants | $\Phi_{\ell\kappa\alpha,\ell'\kappa'\beta}$ | $C_{\kappa\alpha p,\kappa'\alpha'p'}$ | $k_{ia\mu,jb\nu}$ (elastic matrix) |
| Fourier-transformed force constants | $\Phi_{\kappa\alpha,\kappa'\beta}(\bq)$ | (not named) | $\Omega^2_{a\mu b\nu}(\bk)$, opposite phase |
| dynamical matrix | $\tilde\Phi_{\kappa\alpha,\kappa'\beta}(\bq)$ | $D^{\mathrm{dm}}_{\kappa\alpha,\kappa'\alpha'}(\bq)$ | $\tilde\Omega^2_{a\mu b\nu}(\bk)$, opposite phase ($=\tilde\Phi(\bq)^*$) |
| branch | $\nu$ | $\nu$ | $\alpha$ (§§1–4; $\sigma$ already in §§4.5–4.6), $\sigma$ from p. 61 (§5.2.1) on |
| eigenvector | $e_{\kappa\alpha,\nu}(\bq)$ | $e_{\kappa\alpha,\nu}(\bq)$ | $e_{a\mu\alpha}(\bk)$, complex conjugate of ours (opposite phase) |
| phonon wavevector | $\bq$ | $\bq$ | $\bk$ in §1, $\bq$ from §2 |
| normal coordinate; momentum | $Q_{\bq\nu}$; $P_{\bq\nu}=\dot Q^*_{\bq\nu}$ | $z_{\bq\nu}=Q_{\bq\nu}/\sqrt{M_0}$ | $q_{\bk\alpha}$; $p_{\bk\alpha}=\dot q_{\bk\alpha}$ (our $P_{-\bq\nu}$) |
| zero-point amplitude | $(\hbar/2M_\kappa\omega_{\bq\nu})^{1/2}$ | $(M_0/M_\kappa)^{1/2}l_{\bq\nu}$ | written out |
| phonon field | $\phi_{\bq\nu}$ | $\hat a_{\bq\nu}+\hat a^\dagger_{-\bq\nu}$ | $a_{\bq\sigma}+a^+_{-\bq\sigma}$ |
| cell volume; crystal volume | $\Omega_{\mathrm V}$; $\Omega_{\mathrm V}^{\mathrm{SC}}$ | $\Omega$; $N_p\Omega$ | $\Omega_{\mathrm{uc}}$ (also $\Omega_0$); $\Omega$ |
| band energy | $E_{n\bk}$ | $\varepsilon_{n\bk}$ | $E_{n\bk}$ |
| e–ph matrix element | $g_{mn\nu}(\bk,\bq)$, $N_\ell^{-1/2}$ outside | same | $g_{nm\alpha}(\bk,\bq)=N^{-1/2}g^{\mathrm G}_{nm}(\bk-\bq,\bq)$ |
| Born charge | $\mathbf Z^*_\kappa$ dimensionless | $\mathbf Z^*_\kappa$ dimensionless | $Z^*_{a\mu\nu}$ with dimension of charge, $=eZ^*$ |
| dielectric constants | $\varepsilon_\infty$, $\varepsilon_{\mathrm s}$, $\varepsilon^{\mathrm e}$ | $\epsilon^\infty$, $\epsilon^0$ | $\epsilon_\infty$, $\epsilon(0)$; $\epsilon^{\mathrm e}$ |
| phonon propagator | $F$ | $D$ (G Eq. 109) | $F$ |
| phonon self-energy | $\Pi^{\mathrm{ph}}$ | $\Pi$ | $\Omega^{\mathrm e}$ |
| photon propagator; free parts | $D$; $D^{0\mathrm L}$, $D^{0\mathrm T}$ | none | $D$; $D^{0l}$, $D^{0t}$ |
| photon self-energy | $\Pi$, $\Pi^{\mathrm e}$, $\Pi^{\mathrm{ion}}$ | none ($P$ is the density response) | $\Pi_{jj}$, $\Pi^{\mathrm e}$, $\Pi^{\mathrm{ion}}$ |
| screened interaction | $W$, $W^{\mathrm e}$, $W^{\mathrm{ph}}$ | $W$, $W_{\mathrm e}$, $W_{\mathrm{ph}}$ | $D_{11}$, $W^{\mathrm e}=D^{\mathrm e}_{11}$ |
| units | Gaussian, constants explicit | SI | SI in §3, Gaussian from §4 |

---

## 10. Writing rules for a self-contained section

1. **Picture, then model, then general case, then sources.** Each section opens with the physical picture.
   Next comes the simplest concrete model, with a geometric figure, and then the general formalism.
   The source mapping goes in footnotes.
2. **Define every symbol at first use**, in words, with its range and units.
   A display equation that brings in new symbols is followed by a "where" key (`<p class="where">`).
3. **Derivations motivate every step.** Use `<ol class="steps">`. Each step says what is done and why (`<span class="why">`).
   A step is never "it follows that".
4. **Geometric figures** go wherever an index, a vector or a matrix element has a spatial meaning.
   A figure comes before the formalism that uses it.
5. **Source symbols only in footnotes** (`<sup class="fn">`). The same goes for source equation numbers.
   The main text never needs the papers. It may say "the review" or "the EPP manuscript" when it credits an idea.
6. **Background boxes** (`<details class="background">`) teach material the reader may already know.
   They are closed by default and the main text never depends on their content.
   Lattice dynamics, Born charges and photon quantization are taught from scratch in the main text.
   Second quantization, Green's functions, Dyson equations and RPA/GW are assumed, with background boxes.
7. **Glossary.** Every symbol a section defines gets an entry in the section's `<dl class="glossary">`. The entry gives
   its meaning, units, where it is first defined, and its source symbols. The floating "Symbols" panel collects all entries.
8. **Equation numbers** are the primer's own, (A1.1), (A1.2), …; source equation numbers live in footnotes.
9. **SI footnotes** go on key equations with $e$, $c$ or Coulomb kernels (Section 7).

## 11. Linting (machine-checkable)

`v2/build/lint.py` flags:
- source-only symbols in the main text (outside footnotes and source maps): `\Delta\tau`, `C_{\kappa`, `D^{\rm dm}`, `N_p`, `M_0`,
  `l_{\bq\nu}`, `\varepsilon_{n`, `\Omega^2`, `\Omega_{\rm uc}`, `Q_{ia`, `m_a`, `\mu` as a Cartesian index;
- italic labels: `^{\rm `, `_{\rm ` (use `\mathrm`), `e^{i`, `e^{-i`, `e^{\ii` (use `\ee^{\ii`), `^{ion}`, `^{ph}`, `^{e}`;
- a symbol inventory (first use in reading order) for the naive-reader check, plus glossary entries never used and symbols with no glossary entry.
