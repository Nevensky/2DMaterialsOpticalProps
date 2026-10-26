# Short-range coupling and higher multipoles

Labels: **[est]** established in the literature; **[derived]** algebra done for this note from EPP's own equations, not printed in EPP; **[proposal]** a suggested extension of EPP; **[spec]** speculation. Numbers marked *illustrative* come from the primer's toy models or literature presets.

## 1. What is missing in EPP

**No deformation-potential vertex.** EPP derives the rigid-ion vertex $g^0$ in §2 (E.100–E.108), then drops it in §4. In the Lagrangian (E.134) the electrons feel the frozen potential $V_{\rm cr}(\br)=\sum_{ia}v_a(\br-\bR^0_{ia})$ (E.136). The ions reach the field only through the Born-dipole current $\frac1c\dot{\mathbf p}_{ia}\cdot\bA(\bR^0_{ia})$, $\mathbf p=\hat Z^*\mathbf Q$ (E.135). No term contains both $\br_e$ and $\mathbf Q_{ia}$ ([§B.2](../../sections/B02-lagrangian.html)), so "plausibly closed the scheme" (p. 26) holds only for the polar channel. As a result:

- In Si, diamond or graphene ($Z^*=0$), EPP's electron–phonon coupling vanishes.
- Non-polar optical modes (MoS₂ A$_1'$) and the acoustic deformation potential are absent.
- Intervalley coupling reduces to the screened dipole tail.

The table in [§B.1](../../sections/B01-foundations.html) lists $g^{\mathcal S}$ as "DFPT $g$ minus $g^{\mathcal L}$; not used."

**Only the dipole moment.** Point dipoles at the sites (E.241′, 314′) give the ionic current vertex $j_{\alpha\sigma}(\bq+\bG)=\sum_{a\mu}Z^*_{a\alpha\mu}e_{a\mu\sigma}(\bq)e^{-i(\bq+\bG)\cdot\br^0_a}/\sqrt{m_a}$ (E.295′). This is the zeroth moment of the current a moving atom carries.

- **Optical modes.** The vertex has the full $O(q^{-1})$ Fröhlich term.
- **Acoustic modes.** The $O(q^0)$ term is incomplete. With the basis phase (E.131b, [§Z.3](../../sections/Z03-corrections.html)) only the internal-strain dipole survives. The clamped-ion response, carried by dynamical quadrupoles, is missing.
- **GaAs check.** The internal-strain term alone gives ≈1.2 C/m² along [111] (primer toy), against a measured $(2/\sqrt3)|e_{14}|\approx0.19$ C/m² (*illustrative*).

**Force constants.** EPP p. 21 assumes force constants that vanish beyond one cell ($k^{\mathcal S}$). DFPT force constants carry dipolar and multipolar tails, and EPP does not say what to subtract.

**Local fields.** From E.380 on, only $\bG=\bG'=0$ is kept, and $\Pi^{\rm e}$ is RPA without $f^{xc}$. A $\bG$-resolved ionic $\Pi$ (E.454a of the correction) needs Ewald removal of the on-site self term. The draft has neither.

## 2. Why it matters

- **Transport.** In MoS₂ the room-temperature mobility (~410 cm²/Vs) is limited by optical deformation potentials. They are 4.1×10⁸ eV/cm for the homopolar intravalley mode and 2.6×10⁸ eV/cm for the intervalley LO mode (Kaasbjerg *et al.* 2012). EPP has neither. In GaAs, intervalley scattering dominates hot-carrier relaxation (Sjakste *et al.* 2015).
- **Quadrupoles [est].**
  - Si, GaAs, GaP: neglecting quadrupoles gives large mobility errors (Brunin *et al.* 2020), up to ~30% per the authors' institution.
  - GaN: with quadrupoles, acoustic scattering of electrons falls from about half the total rate to a negligible share, and for holes from ~75% to below 50% (Jhalani *et al.* 2020).
  - MoS₂: neglecting quadrupoles gives 23% and 76% errors in the electron and hole Hall mobilities (Poncé *et al.* 2023).
- **2D power counting [derived].** In 2D ($v_Q=2\pi/Q$) the dipole vertex is finite at $Q\to0$ ([§B.10](../../sections/B10-2d.html)). The piezoelectric/quadrupole acoustic coupling is $O(Q)$, like the acoustic deformation potential. Optical deformation potentials are $O(Q^0)$, like 2D Fröhlich. So in 2D, unlike 3D, short-range coupling is leading order, and EPP §8 alone cannot give 2D carrier or exciton scattering.
- **For an RPA/GW/BSE practitioner.** Exciton linewidths in WS₂ and MoSe₂ come from intravalley phonon scattering and from intervalley scattering into dark states, at $q\approx K$ (Selig *et al.* 2016). That is short-range territory, as is phonon-assisted absorption. DFPT screens with local fields (G.53) and EPP's $\Pi^{\rm e}$ has none. Reconciling the two is the $\epsilon_M$ / $\bar v$ split of your BSE code.

## 3. What standard theory does

- **DFPT vertex.** $g_{mn\nu}(\bk,\bq)=\langle u_{m\bk+\bq}|\Delta_{\bq\nu}v^{\rm KS}|u_{n\bk}\rangle$ (G.38), with $\partial v^{\rm KS}=(\epsilon^{{\rm H}xc})^{-1}\partial v^{\rm en}$ (G.53). It contains all orders, statically screened ([§A.3](../../sections/A03-dfpt.html)).
- **Decay in real space.** The potential of a displaced atom decays as $R^{-3}$ or faster in nonpolar insulators, so $g$ is finite at $q\to0$. In polar insulators it decays as $R^{-2}$ and $g\sim1/q$ (G §VI.A.2–3).
- **Split for interpolation [est].** Subtract (G.180) on the coarse grid, Wannier-interpolate $g^{\mathcal S}$ (G.177–G.179), and add $g^{\mathcal L}$ back, Gaussian-filtered (footnote 7; Verdi and Giustino 2015; Sjakste *et al.* 2015). Force constants are treated the same way with the Ewald dipole–dipole term (Gonze and Lee 1997).
- **Quadrupoles [est].** Long-wave DFPT gives the dynamical quadrupoles $Q^{\beta\gamma}_{\kappa\alpha}$ (Royo and Stengel 2019). With $\bk=\bq+\bG$ and $\rho(\bk)=-i\bk\cdot\mathbf p-\tfrac12k_\beta k_\gamma Q^{\beta\gamma}$, the long-range vertex is, schematically,
  $$g^{\mathcal L}_{mn\nu}\propto\sum_{\kappa}\sum_{\bG\neq-\bq}\frac{\big[ik_\beta Z^*_{\kappa,\beta\alpha}+\tfrac12k_\beta k_\gamma Q^{\beta\gamma}_{\kappa\alpha}\big]e_{\kappa\alpha\nu}(\bq)}{\bk\cdot\epsilon^\infty\cdot\bk}\,\langle\psi_{m\bk_e+\bq}|e^{i\bk\cdot(\br-\boldsymbol\tau_\kappa)}[1+\dots]|\psi_{n\bk_e}\rangle .$$
  The quadrupole term is $O(q^0)$ and direction dependent. It spoils interpolation even in Si (Park *et al.* 2020; Brunin *et al.* 2020). The "$\dots$" stands for two corrections. One is a local-field term in Brunin *et al.*, the Hxc response to a macroscopic field (its form was not re-checked here). The other, at quadrupolar order, is a Berry-connection term needed for Wannier-gauge covariance (Poncé *et al.* 2023).
- **Multipolar force constants.** The same multipole orders in the force constants fix the sound velocity of piezoelectrics (Royo, Hahn, Stengel 2020).
- **2D.** Truncated-Coulomb DFPT (Sohier *et al.* 2017) and exact quasi-2D electrostatics (Royo and Stengel 2021).

## 4. How it could enter the EPP formalism

**Route A** unfreezes $V_{\rm cr}$ and replaces the Born dipoles by bare pseudo-ion charges. $Z^*$, $Q$, $\epsilon^\infty$ and $g^{\mathcal S}$ must then emerge from a $\bG\bG'$-resolved $\Pi^{\rm e}$ with $f^{xc}$. This is Hedin–Baym (G.121) in vector form: exact, but expensive. **Route B**, P0–P5 below, keeps EPP's structure and adds DFPT inputs under subtraction rules. All of it is **[proposal]** unless marked [derived].

**P0 · Partition.** Microscopic dressing ($\bG\neq0$ Coulomb and $f^{xc}$) goes into the inputs. EPP's lines carry only the macroscopic longitudinal photon and the transverse photon, so the electronic self-energy must be macroscopically proper:
$$\bar\chi=\chi^0+\chi^0(\bar v+f^{xc})\bar\chi,\quad \bar v_{\bG}=v(\bq+\bG)(1-\delta_{\bG0}),\qquad \Pi^{\rm e}_{\rho\rho}\to e^2\bar\chi_{00}(\bq,\omega).\tag{P0}$$
Since $1-v_{\bq}e^2\bar\chi_{00}=\epsilon_M$ [est], EPP then regenerates (via E.350) the $\epsilon^\infty$ that DFPT uses, and $Z^*$ and $Q$ are already zero-macroscopic-field quantities.

**P1 · Multipolar ionic current.** Replace term (6) of E.134 by $\frac1c\int\mathbf j^{\rm ion}\cdot\bA$, with $\mathbf j^{\rm ion}(\br)=\sum_{ia}\dot Q_{ia\mu}\boldsymbol{\mathcal J}_{a\mu}(\br-\bR^0_{ia})$. Here $\boldsymbol{\mathcal J}$ is the current-density response to the atomic velocity (Stengel 2013). Its moments are $\int\mathcal J_{a\mu,\alpha}=Z^*_{a\alpha\mu}$ and $\int r_\beta\mathcal J_{a\mu,\alpha}=\frac12\mathcal Q_{a,\alpha\beta\mu}$. Then E.295′ becomes
$$j_{\alpha\sigma}(\bk)=\sum_{a\mu}\Big[Z^*_{a\alpha\mu}-\tfrac i2k_\beta\,\mathcal Q_{a,\alpha\beta\mu}+O(k^2)\Big]\frac{e_{a\mu\sigma}(\bq)\,e^{-i\bk\cdot\br^0_a}}{\sqrt{m_a}},\qquad\bk=\bq+\bG.\tag{P1}$$
- **[derived]** $\bk\cdot\mathbf j$ reproduces the dipole plus quadrupole charge of §3, so the symmetric part of $\mathcal Q$ is fixed by $Q$.
- The antisymmetric part couples to $\bB$ **[spec: negligible]**.
- **[derived]** P1 is linear in $\dot Q$, so E.183′, E.304′ and the gauge sum rule E.310 carry over with $\hat Z^*_a\to\tilde Z_a(\bk)=\hat Z^*_a-\frac i2\bk\cdot\mathcal Q_a$.

**P2 · Bare phonons.** Take DFPT minus exactly what EPP regenerates statically:
$$\tilde\Omega^2_0=\tilde C^{\rm DFPT}-\tilde C^{\mathcal L},\qquad \tilde C^{\mathcal L}_{a\mu,b\nu}(\bq)=\frac{4\pi}{\Omega_{\rm uc}\sqrt{m_am_b}}\sum_{\bG}{}^{\rm Ewald}\frac{[\bk\cdot\tilde Z_a]^*_\mu[\bk\cdot\tilde Z_b]_\nu}{\bk\cdot\epsilon^\infty\cdot\bk}\,e^{i\bk\cdot(\br^0_a-\br^0_b)}.\tag{P2}$$
This is E.454a with $Z^*\to\tilde Z$: dipole–dipole, dipole–quadrupole and quadrupole–quadrupole terms, with the self term removed. **[derived]** EPP's $\bG=0$-only term is not periodic in $\bq$, so it jumps at the zone boundary. Keep the Gaussian-filtered $\bG$ sum (as in G.180), in $\tilde C^{\mathcal L}$ and in EPP alike. What the filter removes is analytic and stays in $\tilde\Omega^2_0$.

**P3 · Vertex.** Fix $g^{\mathcal S}$ by matching to DFPT in the static, instantaneous limit:
$$g^{\rm tot}(\omega)=g^{\mathcal S}+\delta g^{\mathcal S}(\omega)+g^{\mathcal L,{\rm dyn}}(\omega),\qquad g^{\mathcal S}\equiv g^{\rm DFPT}-g^{\mathcal L,{\rm dyn}}\big|_{c\to\infty,\,\omega=0}.\tag{P3}$$
$g^{\mathcal L,{\rm dyn}}$ is E.416 with P0 and P1, in Giustino's labels with the E.448a sign convention. Its static limit is (G.180) plus quadrupoles, so $g^{\rm tot}(c\to\infty,0)=g^{\rm DFPT}$ by construction.

**P4 · Dynamic screening of $g^{\mathcal S}$.** Only the frequency-dependent remainder is added, as in EPP's own $\Pi(\omega)-\Pi(0)$:
$$\delta g^{\mathcal S}_{mn\sigma}(\omega)=\rho^*_{n\bk_e,m\bk_e+\bq}\big[W^{\rm e}(\bq,\omega)\Lambda_\sigma(\bq,\omega)-W^{\rm e}(\bq,0)\Lambda_\sigma(\bq,0)\big].\tag{P4}$$
$\Lambda_\sigma$ is the bubble with one charge vertex and one proper vertex $\bar g^{\mathcal S}$ (screened by $\bar v+f^{xc}$ only). P4 matters in doped sheets ([§B.9](../../sections/B09-electrostatic.html)), not in gapped insulators.

**P5 · Dyson equations.**
$$F=F^0+F^0\Omega^{\rm tot}F,\qquad \Omega^{\rm tot}=\Omega^{\rm e}(\omega)+\big[\Sigma^{\mathcal S}(\omega)-\Sigma^{\mathcal S}(0)\big].\tag{P5}$$
- $\Omega^{\rm e}$ (E.399a with P1) enters in full, because P2 removed its static part.
- $\Sigma^{\mathcal S}$, every bubble with a $g^{\mathcal S}$ vertex including the $g^{\mathcal S}$–$J$ cross terms, enters only nonadiabatically (EPP p. 21; G.125).
- Inside $\Pi^{\rm ion}$ (E.352′), dress the phonons only by $\Sigma^{\mathcal S}(\omega)-\Sigma^{\mathcal S}(0)$, so that $D$ (E.386–E.387) adds the dipolar dressing once.
- **[derived]** That dressing vanishes at $\omega=0$, so $\Pi^{\rm para}(0)$ and E.310 survive.

## 5. Pitfalls

- **Inconsistent $\epsilon$.** If the subtraction uses DFPT's $\epsilon^\infty$ but EPP regenerates the RPA head without local fields, about 10% larger in sp semiconductors (*illustrative*, §B.9), $g^{\mathcal S}$ keeps a tail $\propto(1/\epsilon_{\rm EPP}-1/\epsilon_{\rm DFPT})/q$. That tail is double-counted and breaks interpolation.
- **Double screening.** DFPT's $g$ is already screened at $\bG=0$, so routing it through $D^{\rm e}$ again screens it twice. A screened $g^{\mathcal S}$ on both sides of a bubble overcounts by $\epsilon$: this is the G.145 versus G.146 ambiguity.
- **Routing and periodicity.** Every nonanalytic $\bG=0$ charge (dipole, quadrupole, piezoelectric) goes through $A_l$, never into $g^{\mathcal S}$ or $\tilde\Omega^2_0$. A $\bG=0$-only $g^{\mathcal L,{\rm dyn}}$ makes $g^{\mathcal S}$ jump at the zone boundary, so use the filtered $\bG$ sum.
- **Force constants.** Subtract the same multipole orders you regenerate. Keep the acoustic-sum-rule on-site block (missing from E.4) and check that results do not depend on the Ewald parameter.
- **Basis phase.** Without $e^{-i\bk\cdot\br^0_a}$ and Giustino's sign of E.15, a spurious acoustic coupling appears (GaAs 0.6–1.8 C/m², *illustrative*), comparable to the quadrupole term itself.
- **Nonlocal pseudopotentials.** Current vertices need the $i[V_{\rm NL},\br]$ term, as in velocity-gauge optics (Dreyer *et al.* 2018).
- **2D.** Use one sheet propagator in DFPT, in the subtraction and in EPP §8. Include $Z^*_{zz}$ and out-of-plane quadrupoles.

## 6. State of the art (verified)

1. S. de Gironcoli, S. Baroni, R. Resta, PRL **62**, 2853 (1989), [10.1103/PhysRevLett.62.2853](https://doi.org/10.1103/PhysRevLett.62.2853). Linear-response piezoelectric constants of III–V semiconductors.
2. X. Gonze, C. Lee, PRB **55**, 10355 (1997), [10.1103/PhysRevB.55.10355](https://doi.org/10.1103/PhysRevB.55.10355). $Z^*$, $\epsilon^\infty$ and the Ewald dipole–dipole split of force constants.
3. M. Stengel, PRB **88**, 174106 (2013), [10.1103/PhysRevB.88.174106](https://doi.org/10.1103/PhysRevB.88.174106). Current response to long acoustic phonons. Nonlocal-pseudopotential implementation: C. E. Dreyer, M. Stengel, D. Vanderbilt, PRB **98**, 075153 (2018), [10.1103/PhysRevB.98.075153](https://doi.org/10.1103/PhysRevB.98.075153).
4. C. Verdi, F. Giustino, PRL **115**, 176401 (2015), [10.1103/PhysRevLett.115.176401](https://doi.org/10.1103/PhysRevLett.115.176401). Long-range Fröhlich vertex and the $g^{\mathcal S}+g^{\mathcal L}$ interpolation.
5. J. Sjakste, N. Vast, M. Calandra, F. Mauri, PRB **92**, 054307 (2015), [10.1103/PhysRevB.92.054307](https://doi.org/10.1103/PhysRevB.92.054307). Polar Wannier interpolation in GaAs; Fröhlich versus intervalley scattering.
6. M. Royo, M. Stengel, PRX **9**, 021050 (2019), [10.1103/PhysRevX.9.021050](https://doi.org/10.1103/PhysRevX.9.021050). Long-wave DFPT for dynamical quadrupoles.
7. G. Brunin *et al.*, PRL **125**, 136601 (2020), [10.1103/PhysRevLett.125.136601](https://doi.org/10.1103/PhysRevLett.125.136601), and PRB **102**, 094308 (2020), [10.1103/PhysRevB.102.094308](https://doi.org/10.1103/PhysRevB.102.094308). Quadrupolar vertex; mobilities of Si, GaAs, GaP.
8. V. A. Jhalani, J.-J. Zhou, J. Park, C. E. Dreyer, M. Bernardi, PRL **125**, 136602 (2020), [10.1103/PhysRevLett.125.136602](https://doi.org/10.1103/PhysRevLett.125.136602). Piezoelectric coupling from quadrupoles in GaN.
9. J. Park, J.-J. Zhou, V. A. Jhalani, C. E. Dreyer, M. Bernardi, PRB **102**, 125203 (2020), [10.1103/PhysRevB.102.125203](https://doi.org/10.1103/PhysRevB.102.125203). Quadrupole vertex tested against DFPT in Si and PbTiO₃.
10. M. Royo, K. R. Hahn, M. Stengel, PRL **125**, 217602 (2020), [10.1103/PhysRevLett.125.217602](https://doi.org/10.1103/PhysRevLett.125.217602). Multipolar force constants and piezoelectric sound velocities.
11. M. Royo, M. Stengel, PRX **11**, 041027 (2021), [10.1103/PhysRevX.11.041027](https://doi.org/10.1103/PhysRevX.11.041027). Exact quasi-2D long-range electrostatics.
12. S. Poncé, M. Royo, M. Gibertini, N. Marzari, M. Stengel, PRL **130**, 166301 (2023), [10.1103/PhysRevLett.130.166301](https://doi.org/10.1103/PhysRevLett.130.166301), and PRB **107**, 155424 (2023), [10.1103/PhysRevB.107.155424](https://doi.org/10.1103/PhysRevB.107.155424). Gauge-covariant quadrupoles; 2D mobilities.
13. T. Sohier, M. Calandra, F. Mauri, PRB **96**, 075448 (2017), [10.1103/PhysRevB.96.075448](https://doi.org/10.1103/PhysRevB.96.075448). 2D DFPT with Coulomb truncation and gating.
14. K. Kaasbjerg, K. S. Thygesen, K. W. Jacobsen, PRB **85**, 115317 (2012), [10.1103/PhysRevB.85.115317](https://doi.org/10.1103/PhysRevB.85.115317). MoS₂ mobility limited by optical deformation potentials.
15. M. Selig *et al.*, Nat. Commun. **7**, 13279 (2016), [10.1038/ncomms13279](https://doi.org/10.1038/ncomms13279). TMD exciton linewidths from intravalley and intervalley phonon scattering.

## 7. First steps

1. **Derive and check (days).** Verify E.310 with P1 in the primer's three-atom chain ([§B.5](../../sections/B05-photon-se.html)). Show that P5 at $c\to\infty$, $\omega=0$ gives the multipolar force constants (Royo, Hahn, Stengel 2020).
2. **Bulk benchmarks (1–2 weeks).** Si and GaAs, with $Z^*$, $Q$, $\epsilon^\infty$ from ABINIT long-wave DFPT. Build $g^{\mathcal S}$ by P3, check its real-space decay, and compare LA[111] with $e_{14}$.
3. **$\epsilon$ consistency (1–2 weeks).** Compute the head with your RPA code, with and without local fields and $f^{xc}$, and measure the residual $1/q$ tail in $g^{\mathcal S}$.
4. **2D (1–2 months).** Monolayer MoS₂ (A$_1'$, LO, LA at $K$) with truncated DFPT and the §8 propagator, plus P4 for gated samples. Compare with Poncé *et al.* and Kaasbjerg *et al.*
5. **Hand-off (months).** Pass $g^{\rm tot}$ to the electron self-energy and exciton–phonon topics.

## 8. Open questions

- Is Route A tractable: does a $\bG\bG'$ RPA+$f^{xc}$ self-energy reproduce $Z^*$, $Q$, $g^{\mathcal S}$?
- Should $Q$ and $g^{\mathcal S}$ be screened dynamically by carriers, as the dipole term is?
- Which vertex, bare or screened, belongs on each side of $\Sigma^{\mathcal S}$ when the screening is dynamic?
- Acoustic modes never enter the light cone ($v_s\ll c/\sqrt\epsilon$) **[derived]**. Is finite $c$ ever relevant to short-range physics?
- **[spec]** Does the magnetic part of $\mathcal Q$ matter for chiral phonons in 2D?
