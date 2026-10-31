# Short-range coupling and higher multipoles

Labels: **[est]** established in the literature; **[derived]** algebra done for this note from EPP's own equations, not printed in EPP; **[proposal]** a suggested extension of EPP; **[spec]** speculation. Numbers marked *illustrative* come from the primer's toy models or literature presets. References marked † in §6 could not be re-verified in this revision.

## 1. What is missing in EPP

**No deformation-potential vertex.** EPP derives the rigid-ion vertex $g^0$ in §2 (E.100–E.108), then drops it in §4. In the Lagrangian (E.134) the electrons feel the frozen potential $V_{\rm cr}(\br)=\sum_{ia}v_a(\br-\bR^0_{ia})$ (E.136). The ions reach the field only through the Born-dipole current $\frac1c\dot{\mathbf p}_{ia}\cdot\bA(\bR^0_{ia})$, $\mathbf p=\hat Z^*\mathbf Q$ (E.135). This is E.134 as corrected in place in [§Z.3](../../sections/Z03-corrections.html); the draft prints $\bA(\bR^0_i)$ (p. 26). No term contains both $\br_e$ and $\mathbf Q_{ia}$ ([§B.2](../../sections/B02-lagrangian.html)), so "plausibly closed the scheme" (p. 26) holds only for the polar channel. As a result:

- In Si, diamond or graphene ($Z^*=0$), EPP's electron–phonon coupling vanishes.
- Non-polar optical modes (MoS₂ A$_1'$) and the acoustic deformation potential are absent.
- Intervalley coupling reduces to the screened dipole tail.

[§B.1](../../sections/B01-foundations.html) lists $g^{\mathcal S}$ as "not used."

**Only the dipole moment.** Point dipoles at the sites (E.241′, 314′) give the ionic current vertex $j_{\alpha\sigma}(\bq+\bG)=\sum_{a\mu}Z^*_{a\alpha\mu}e_{a\mu\sigma}(\bq)e^{-i(\bq+\bG)\cdot\br^0_a}/\sqrt{m_a}$ (E.295′). This is the zeroth moment of the current a moving atom carries.

- **Optical modes.** The vertex has the full $O(q^{-1})$ Fröhlich term.
- **Acoustic modes.** The $O(q^0)$ term is incomplete. With the basis phase (E.131b, [§Z.3](../../sections/Z03-corrections.html)) only the internal-strain dipole survives. The clamped-ion response, carried by dynamical quadrupoles, is missing.
- **GaAs check.** The internal-strain term alone gives ≈1.2 C/m² along [111] in the primer's central-force toy, against a measured $(2/\sqrt3)|e_{14}|\approx0.19$ C/m² (*illustrative*). The toy has Kleinman parameter $\zeta=1$. Its next-nearest-neighbour springs join atoms of the same sublattice and do not couple strain to the internal displacement, so $\zeta$ is set by the central nearest-neighbour springs alone **[derived]**. In general $e_{14}^{\rm int}=Z^*|e|\zeta/a^2=1.04\,\zeta$ C/m² for GaAs ($Z^*=2.07$, $a=5.65$ Å), or $1.20\,\zeta$ C/m² along [111] **[derived]**. With a realistic $\zeta\approx0.5$–$0.75$ (literature range, unverified) the internal-strain term is about 0.6–0.9 C/m² along [111]. The clamped-ion (quadrupolar) term still has to cancel most of it. The split into the two parts should be quoted from de Gironcoli, Baroni and Resta (1989)† once that paper is checked.

**Force constants.** EPP p. 21 assumes force constants that vanish beyond one cell ($k^{\mathcal S}$). DFPT force constants carry dipolar and multipolar tails, and EPP does not say what to subtract.

**Local fields.** From E.380 on, only $\bG=\bG'=0$ is kept, and $\Pi^{\rm e}$ is RPA without $f^{xc}$. A $\bG$-resolved ionic $\Pi$ (E.454a of the correction) needs Ewald removal of the on-site self term. The draft has neither. Nor does it dress the electron's coupling to the macroscopic field with local fields (§3, P0).

## 2. Why it matters

- **Transport.** The room-temperature mobility of n-type monolayer MoS₂ (~410 cm²/Vs) is dominated by optical-phonon scattering through intra- and intervalley deformation potentials and the Fröhlich interaction (Kaasbjerg *et al.* 2012; figure and summary as recalled from the abstract, not re-checked here). EPP has only the Fröhlich part. Above the L-valley threshold, hot electrons in GaAs also relax by intervalley ($\Gamma\to$L, X) scattering, which needs short-range coupling. Ab initio hot-electron thermalization in GaAs is due to Sjakste, Vast and Tyuterev (2007) and Bernardi *et al.* (2015), who resolved carrier lifetimes valley by valley (G §X.A, p. 015003-46). Intervalley scattering times in Ge were computed the same way (Tyuterev *et al.* 2011).
- **Quadrupoles [est].**
  - Si, GaAs, GaP: neglecting quadrupoles changes the computed mobilities significantly (Brunin *et al.* 2020†; magnitude to be quoted from the paper).
  - GaN: quadrupoles reshape the piezoelectric (acoustic) scattering of electrons and holes (Jhalani *et al.* 2020†). The earlier draft of this note said the acoustic share falls from about half of the total rate to a negligible one for electrons, and from ~75% to below 50% for holes. Those numbers are not re-checked.
  - MoS₂: neglecting quadrupoles gives 23% and 76% errors in the electron and hole Hall mobilities (Poncé *et al.* 2023†; not re-checked in this revision).
- **2D power counting [derived].** In a sheet with electronic screening length $r_0$ the macroscopic interaction is $W(Q)=2\pi e^2/[Q(1+r_0Q)]$ ([§B.10](../../sections/B10-2d.html)). Per unit displacement $u$, the macroscopic charges are:
  - optical dipole, $\rho_{\rm dip}\sim Z^*Qu$;
  - quadrupole, $\rho_{\rm quad}\sim \mathcal Q\,Q^2u$;
  - acoustic piezoelectric (internal strain plus clamped ion), $\rho_{\rm piezo}\sim Q\cdot(e\,Qu)$, i.e. $Q$ times a strain $\sim Qu$.

  The vertex is $g\sim W\rho$. For the optical dipole that gives $g\sim 2\pi e^2Z^*u/(1+r_0Q)$. For $Q\ll1/r_0$ this is $O(Q^0)$, like the optical deformation potential. The acoustic piezoelectric and quadrupole terms are $O(Q)$, like the acoustic deformation potential. So in 2D, unlike 3D, short-range coupling is leading order. For $Q\gg1/r_0$ the dipole term is further suppressed by $1/(r_0Q)$. In MoS₂, $r_0\approx40$ Å (the §B.10 preset, *illustrative*; not checked against the literature), so $1/r_0\approx0.025$ Å⁻¹. Thermal electrons at 300 K with $m^*\approx0.5\,m_e$ have $k_{\rm th}=\sqrt{2m^*k_BT}/\hbar\approx0.058$ Å⁻¹, and emission of a 48 meV optical phonon (§B.10 preset, *illustrative*) needs $q\approx0.08$ Å⁻¹. Over the relevant 0.06–0.12 Å⁻¹, $r_0Q\approx2.5$–$5$ and $1/(1+r_0Q)\approx0.17$–$0.29$. The finite 2D Fröhlich vertex is that of Sohier, Calandra and Mauri (2016)†. Short-range coupling therefore dominates even more, and EPP §8 alone cannot give 2D carrier or exciton scattering.
- **For an RPA/GW/BSE practitioner.** Exciton linewidths in monolayer TMDs come from intravalley phonon scattering. In W-based TMDs they also come from scattering into dark intervalley excitons, K–K′ and K–Λ, so the phonon momenta are near $q\approx K$ and $q\approx|K-\Lambda|$, both deep in the zone (Selig *et al.* 2016†; which materials were measured and which only modelled is to be checked against the paper). That is short-range territory, as is phonon-assisted absorption (G §IX.B, p. 015003-44). DFPT screens with local fields (G.53) and EPP's $\Pi^{\rm e}$ has none. Reconciling the two is the $\epsilon_M$ / $\bar v$ split of your BSE code.

## 3. What standard theory does

- **DFPT vertex.** $g_{mn\nu}(\bk,\bq)=\langle u_{m\bk+\bq}|\Delta_{\bq\nu}v^{\rm KS}|u_{n\bk}\rangle$ (G.38), with $\partial v^{\rm KS}=(\epsilon^{{\rm H}xc})^{-1}\partial v^{\rm en}$ (G.53). It contains all orders, statically screened ([§A.3](../../sections/A03-dfpt.html); Baroni *et al.* 2001).
- **Decay in real space.** The potential of a displaced atom decays as $R^{-3}$ or faster in nonpolar insulators, a consequence of the analytic properties of the dielectric matrix (Pick, Cohen and Martin 1970). So $g$ is finite at $q\to0$ (Vogl 1976; G p. 015003-34). In polar insulators it decays as $R^{-2}$ and $g\sim1/q$ (G §VI.A.2–3). Vogl (1976) is also the original microscopic multipole treatment of the long-range vertex.
- **Split for interpolation [est].** Subtract (G.180) on the coarse grid, Wannier-interpolate $g^{\mathcal S}$ (G.177–G.179), and add $g^{\mathcal L}$ back, Gaussian-filtered (footnote 7; Verdi and Giustino 2015; Sjakste *et al.* 2015). Force constants are treated the same way with the Ewald dipole–dipole term (Gonze and Lee 1997).
- **Quadrupoles [est].** Long-wave DFPT gives the dynamical quadrupoles $Q^{\beta\gamma}_{\kappa\alpha}$ (Royo and Stengel 2019†). With $\tilde\bq=\bq+\bG$ and $\rho(\tilde\bq)=-i\tilde\bq\cdot\mathbf p-\tfrac12\tilde q_\beta \tilde q_\gamma Q^{\beta\gamma}$, the long-range vertex is, schematically,
  $$g^{\mathcal L}_{mn\nu}\propto\sum_{\kappa}\sum_{\bG\neq-\bq}\frac{\big[i\tilde q_\beta Z^*_{\kappa,\beta\alpha}+\tfrac12\tilde q_\beta \tilde q_\gamma Q^{\beta\gamma}_{\kappa\alpha}\big]e_{\kappa\alpha\nu}(\bq)}{\tilde\bq\cdot\epsilon^\infty\cdot\tilde\bq}\,\langle\psi_{m\bk+\bq}|e^{i\tilde\bq\cdot(\br-\boldsymbol\tau_\kappa)}[1+\dots]|\psi_{n\bk}\rangle .$$
  The quadrupole term is $O(q^0)$ and direction dependent. It spoils interpolation even in Si (Park *et al.* 2020†; Brunin *et al.* 2020†). The "$\dots$" stands for two further $O(q^0)$ corrections.
  - **Local fields of the macroscopic field.** The macroscopic field $E_\gamma e^{i\bq\cdot\br}$ of the dipole induces a periodic Hxc potential $V^{{\rm Hxc},E_\gamma}(\br)$, the Hxc potential per unit macroscopic field (Brunin *et al.*'s $V^{{\rm Hxc},\mathcal E}$ term†). From the dressed vertex P0b below **[derived]**, $1+\sum_{\bG\neq0}e^{i\bG\cdot\br}\big((\bar v+f^{xc})\bar\chi\big)_{\bG0}(\bq,0)=1-iq_\gamma V^{{\rm Hxc},E_\gamma}(\br)+O(q^2)$. Times the dipole numerator, this adds $q_\beta q_\gamma Z^*_{\kappa,\beta\alpha}e_{\kappa\alpha\nu}\langle\psi_{m\bk+\bq}|e^{i\bq\cdot(\br-\boldsymbol\tau_\kappa)}V^{{\rm Hxc},E_\gamma}|\psi_{n\bk}\rangle/(\bq\cdot\epsilon^\infty\cdot\bq)$ to the $\bG=0$ term: $O(q^0)$ and direction dependent, the same order as the quadrupole. Map the factor and sign to Brunin *et al.*'s convention before comparing.
  - **Berry connection.** In the Wannier gauge the overlap $\langle u_{m\bk+\bq}|u_{n\bk}\rangle$ must be expanded to $O(q)$. The $O(q)$ term is the Berry connection, needed for gauge covariance at quadrupolar order (Poncé *et al.* 2023†).
- **Implementations.** The static subtraction with quadrupoles is available in ABINIT (Brunin *et al.*†), Perturbo (Park *et al.*†, Jhalani *et al.*†) and EPW (Poncé *et al.*†). This is as reported in those papers and was not re-checked here.
- **Multipolar force constants.** The same multipole orders in the force constants fix the sound velocity of piezoelectrics (Royo, Hahn, Stengel 2020†).
- **Martin relation.** The clamped-ion piezoelectric tensor is the sublattice sum of the first moments of the induced charge, that is, of the dynamical quadrupoles (Martin 1972†).
- **2D.** Truncated-Coulomb DFPT (Sohier *et al.* 2017†) and exact quasi-2D electrostatics (Royo and Stengel 2021†).

## 4. How it could enter the EPP formalism

**Route A** unfreezes $V_{\rm cr}$ and replaces the Born dipoles by bare pseudo-ion charges. $Z^*$, $Q$, $\epsilon^\infty$ and $g^{\mathcal S}$ must then emerge from a $\bG\bG'$-resolved $\Pi^{\rm e}$ with $f^{xc}$. This is Hedin–Baym (G.121) in vector form: exact, but expensive. To avoid double counting, split $v_a=v_a^{\rm SR}-Z_ae^2/r$. Only $v_a^{\rm SR}$ moves with the ion in $V_{\rm cr}$. The point charges $Z_a$ couple to $\bA$, and the bare ion–ion Coulomb is removed from $k$, because the longitudinal photon regenerates both. **Route B**, P0–P5 below, keeps EPP's structure and adds DFPT inputs under subtraction rules. All of it is **[proposal]** unless marked [derived].

**P0 · Partition.** Microscopic dressing ($\bG\neq0$ Coulomb and $f^{xc}$) goes into the inputs. EPP's lines carry only the macroscopic ($\bG=0$) longitudinal and transverse photon, so everything attached to them must be macroscopically proper.

- **P0a, polarization.** Replace $\Pi^{\rm e}_{jj}$ by $\bar\Pi^{\rm e}_{jj}$. This is the $\bG=\bG'=0$ element of the current–current response that is irreducible with respect to $\bG=0$ photon lines, longitudinal and transverse, with $\bar v+f^{xc}$ in the longitudinal kernel:
  $$\bar\chi=\chi^0+\chi^0(\bar v+f^{xc})\bar\chi,\qquad \bar v_{\bG}=v(\bq+\bG)(1-\delta_{\bG0}).\tag{P0a}$$
  E.350 then gives $\Pi^{\rm e}_{\rho\rho}=e^2\bar\chi_{00}(\bq,\omega)$ automatically, so continuity and the L/T consistency of the photon self-energy are kept. Replacing only the $\rho\rho$ slot would break both. The transverse local-field kernel needs TDCDFT. Alternatively, it follows from the $q\to0$ equality of the longitudinal and transverse heads in cubic crystals.
- **P0b, vertex.** The electron's coupling to the macroscopic photon is dressed by the same kernel:
  $$\langle m\bk+\bq|e^{i\bq\cdot\br}|n\bk\rangle\;\to\;\langle m\bk+\bq|e^{i\bq\cdot\br}\Big[1+\sum_{\bG\neq0}e^{i\bG\cdot\br}\big((\bar v+f^{xc})\bar\chi\big)_{\bG0}(\bq,\omega)\Big]|n\bk\rangle .\tag{P0b}$$
  This is the EPP counterpart of DFPT's $V^{{\rm Hxc},\mathcal E}$ term (§3). Use the same dressed vertex in $g^{\mathcal L,{\rm dyn}}$ and in the P3 subtraction.
- **Consequences.** Since $1-v_{\bq}e^2\bar\chi_{00}=\epsilon_M$ [est], EPP regenerates (via E.350) the $\epsilon^\infty$ that DFPT uses. $Z^*$ and $Q$ are already zero-macroscopic-field quantities.

**P1 · Multipolar ionic current.** Replace term (6) of E.134 by $\frac1c\int\mathbf j^{\rm ion}\cdot\bA$, with $\mathbf j^{\rm ion}(\br)=\sum_{ia}\dot Q_{ia\mu}\boldsymbol{\mathcal J}_{a\mu}(\br-\bR^0_{ia})$. Here $\boldsymbol{\mathcal J}$ is the current-density response to the atomic velocity (Stengel 2013†). Its moments are $\int\mathcal J_{a\mu,\alpha}=Z^*_{a\alpha\mu}$ and $\int r_\beta\mathcal J_{a\mu,\alpha}=\frac12\mathcal Q_{a,\alpha\beta\mu}$. Then E.295′ becomes
$$j_{\alpha\sigma}(\tilde\bq)=\sum_{a\mu}\Big[Z^*_{a\alpha\mu}-\tfrac i2\tilde q_\beta\,\mathcal Q_{a,\alpha\beta\mu}+O(\tilde q^2)\Big]\frac{e_{a\mu\sigma}(\bq)\,e^{-i\tilde\bq\cdot\br^0_a}}{\sqrt{m_a}},\qquad\tilde\bq=\bq+\bG.\tag{P1}$$
- **[derived]** $\tilde\bq\cdot\mathbf j$ reproduces the dipole plus quadrupole charge of §3, so the symmetric part of $\mathcal Q$ is fixed by $Q$.
- The antisymmetric part couples to $\mathbf B$ **[spec: negligible]**.
- **[derived]** P1 is linear in $\dot Q$, so E.183′, E.304′ and the gauge sum rule E.310 carry over with $\hat Z^*_a\to\tilde Z_a(\tilde\bq)=\hat Z^*_a-\frac i2\tilde\bq\cdot\mathcal Q_a$. This holds by construction (completeness of the eigenvectors), so it is not a test (§7).

**P2 · Bare phonons.** Take DFPT minus exactly what EPP regenerates statically:
$$\tilde\Omega^2_0=\tilde C^{\rm DFPT}-\tilde C^{\mathcal L},\qquad \tilde C^{\mathcal L}_{a\mu,b\nu}(\bq)=\frac{4\pi}{\Omega_{\rm uc}\sqrt{m_am_b}}\sum_{\bG}{}^{\rm Ewald}\frac{[\tilde\bq\cdot\tilde Z_a]^*_\mu[\tilde\bq\cdot\tilde Z_b]_\nu}{\tilde\bq\cdot\epsilon^\infty\cdot\tilde\bq}\,e^{i\tilde\bq\cdot(\br^0_a-\br^0_b)}.\tag{P2}$$
This is E.454a with $Z^*\to\tilde Z$: dipole–dipole, dipole–quadrupole and quadrupole–quadrupole terms, with the self term removed.

- **Order [derived].** Dipole–dipole is $O(1)$, dipole–quadrupole $O(q)$ and quadrupole–quadrupole $O(q^2)$. P2 is complete through $O(q)$ in all blocks. It is complete through $O(q^2)$ only in the acoustic–acoustic block after projection with the acoustic sum rule, which is enough for sound velocities. There, $\sum Z^*=0$ makes the acoustic charge $O(q^2)$, and the dipole–octupole and $\epsilon(q)$ corrections enter only at $O(q^3)$ and $O(q^4)$. The optical and mixed blocks at $O(q^2)$ also need octupoles and the $q^2$ dispersion of $\epsilon^\infty$.
- **[derived]** EPP's $\bG=0$-only term is not periodic in $\bq$, so it jumps at the zone boundary. Keep the Gaussian-filtered $\bG$ sum (as in G.180), in $\tilde C^{\mathcal L}$ and in EPP alike. What the filter removes is analytic and stays in $\tilde\Omega^2_0$.

**P3 · Vertex.** Fix $g^{\mathcal S}$ by matching to DFPT in the static, instantaneous limit:
$$g^{\rm tot}(\omega)=g^{\mathcal S}+\delta g^{\mathcal S}(\omega)+g^{\mathcal L,{\rm dyn}}(\omega),\qquad g^{\mathcal S}\equiv g^{\rm DFPT}-g^{\mathcal L,{\rm dyn}}\big|_{c\to\infty,\,\omega=0}.\tag{P3}$$
$g^{\mathcal L,{\rm dyn}}$ is E.416 with P0a, P0b and P1, in Giustino's labels with the E.448a sign convention. Its static limit is (G.180) plus quadrupoles plus the local-field term of §3, which is the long-range part DFPT codes subtract. Then $g^{\rm tot}(c\to\infty,0)=g^{\rm DFPT}$ by construction. In the Wannier gauge, also expand the overlap to $O(q)$ (Berry connection, §3). If P0b is not used, $g^{\mathcal S}$ keeps the $O(q^0)$ nonanalytic remainder $q_\beta q_\gamma Z^*_{\kappa,\beta\alpha}\langle V^{{\rm Hxc},E_\gamma}\rangle/(\bq\cdot\epsilon^\infty\cdot\bq)$, and its Wannier interpolation fails.

**P4 · Dynamic screening of $g^{\mathcal S}$.** Let $\Lambda_\sigma(\bq,\omega)$ be the $\bG=0$ charge induced by the proper vertex $\bar g^{\mathcal S}$: the bubble with one macroscopic charge vertex and one $\bar g^{\mathcal S}$ vertex, with $\bar g^{\mathcal S}$ screened by $\bar v+f^{xc}$ only. The $Z^*$ and $Q$ in $J$ are total static (ionic plus electronic) charges, so $J$ already contains $\Lambda(0)$. Only the change $\Lambda(\omega)-\Lambda(0)$ may be added:
$$\delta g^{\mathcal S}_{mn\sigma}(\omega)=\tilde\rho^*_{n\bk,m\bk+\bq}\,W^{\rm e}(\bq,\omega)\big[\Lambda_\sigma(\bq,\omega)-\Lambda_\sigma(\bq,0)\big],\tag{P4}$$
where $\tilde\rho$ is the dressed vertex P0b. Equivalently, carry the dynamic part as a frequency-dependent ionic current $J\to J+\Delta J(\omega)$, where $\Delta J$ carries, by continuity, the charge $\Lambda(\omega)-\Lambda(0)$, and $\Delta J(0)=0$.

- **Check [derived].** The exact vertex is $\bar g^{\mathcal S}(\omega)+\tilde\rho^*W^{\rm e}(\omega)[\rho_{\rm ion}+\Lambda(\omega)]$. P3 gives $g^{\mathcal S}=\bar g^{\mathcal S}(0)$ and $g^{\mathcal L,{\rm dyn}}=\tilde\rho^*W^{\rm e}(\omega)[\rho_{\rm ion}+\Lambda(0)]$. Adding P4 reproduces the exact vertex, up to $\bar g^{\mathcal S}(\omega)-\bar g^{\mathcal S}(0)$, the microscopic dynamic screening, which is small in gapped systems.
- P4 matters in doped sheets ([§B.9](../../sections/B09-electrostatic.html)), not in gapped insulators.

**P5 · Dyson equations.**
$$F=F^0+F^0\Omega^{\rm tot}F,\qquad \Omega^{\rm tot}=-(J+\Delta J)^*D^{\rm e}(J+\Delta J)+\big[\bar\Sigma^{\mathcal S}(\omega)-\bar\Sigma^{\mathcal S}(0)\big].\tag{P5}$$
- The first term is E.399a with P0 and P1 plus $\Delta J$. It enters in full, because P2 removed its static part and $\Delta J(0)=0$.
- $\bar\Sigma^{\mathcal S}$ contains only bubbles whose $g^{\mathcal S}$ ends are joined by $\bar v+f^{xc}$, with no $\bG=0$ photon line. Every $\bG=0$ charge is in $J+\Delta J$ and every $\bG=0$ line in $D^{\rm e}$. Its static part is already in the DFPT force constants, so it enters only nonadiabatically (EPP p. 21; G.125).
- Inside $\Pi^{\rm ion}$ (E.352′), use $J+\Delta J$ and dress the phonons only by $\bar\Sigma^{\mathcal S}(\omega)-\bar\Sigma^{\mathcal S}(0)$, so that $D$ (E.386–E.387) adds the macroscopic dressing once.
- **[derived]** Both $\Delta J$ and $\bar\Sigma^{\mathcal S}(\omega)-\bar\Sigma^{\mathcal S}(0)$ vanish at $\omega=0$, so $\Pi^{\rm para}(0)$ and E.310 survive. Gauge invariance at finite $\omega$ with an $\omega$-dependent $\Delta J$ is not checked (§8).

## 5. Pitfalls

- **Inconsistent $\epsilon$.** Suppose $g^{\mathcal L}$ is subtracted with $\epsilon_{\rm DFPT}$ but regenerated with $\epsilon_{\rm EPP}$. Then $g^{\rm tot}$ has a spurious Fröhlich tail $\propto(1/\epsilon_{\rm EPP}-1/\epsilon_{\rm DFPT})/q$. If it is subtracted with $\epsilon_{\rm EPP}$, the same tail sits in $g^{\mathcal S}$ and spoils interpolation. Local fields lower the RPA head by about 10% in sp semiconductors (*illustrative*, [§B.9](../../sections/B09-electrostatic.html)), and ALDA $f^{xc}$ raises it again. Measure the net difference (§7 step 3) instead of assuming 10%. P0a removes the mismatch.
- **Double screening.** DFPT's $g$ is already screened at $\bG=0$, so routing it through $D^{\rm e}$ again screens it twice. A screened $g^{\mathcal S}$ on both sides of a bubble screens it once too often, a factor $1/\epsilon$ too small. This is the G.145 ($g^bg^*$) versus G.146 ($|g|^2$) issue (G p. 015003-25).
- **Routing and periodicity.** Every nonanalytic $\bG=0$ charge (dipole, quadrupole, piezoelectric, local-field) goes through $A_l$, never into $g^{\mathcal S}$ or $\tilde\Omega^2_0$. A $\bG=0$-only $g^{\mathcal L,{\rm dyn}}$ makes $g^{\mathcal S}$ jump at the zone boundary, so use the filtered $\bG$ sum.
- **Force constants.** Subtract the same multipole orders you regenerate. Keep the acoustic-sum-rule on-site block (missing from E.4) and check that results do not depend on the Ewald parameter.
- **Basis phase.** Without $e^{-i\tilde\bq\cdot\br^0_a}$ and Giustino's sign of E.15, a spurious acoustic coupling appears (GaAs 0.6–1.8 C/m² along [111], *illustrative*). It is comparable to or larger than the clamped-ion (quadrupole) term itself, which is about 0.4–0.7 C/m² along [111] if $\zeta=0.5$–$0.75$ (§1; $\zeta$ unverified).
- **Convention mapping.** Map $\mathcal Q$ to ABINIT's $Q$ explicitly. Check the Fourier sign (EPP uses $e^{-i\tilde\bq\cdot\br^0_a}$; this note writes $\rho=-i\tilde\bq\cdot\mathbf p-\frac12\tilde\bq\tilde\bq Q$), the factor $\frac12$, the conversion from e·bohr to EPP units, and $Z^*=|e|Z^*_{\rm ABINIT}$ ([§Z.3](../../sections/Z03-corrections.html)). A mismatch flips the sign of the dipole–quadrupole interference: the leading correction for optical modes and the clamped-ion piezoelectric term for acoustic modes. Check the sign of the dQ term against Brunin *et al.*'s published Si or GaAs $g$ before any physics run.
- **Nonlocal pseudopotentials.** Current vertices need the $i[V_{\rm NL},\br]$ term, as in velocity-gauge optics (Dreyer *et al.* 2018†).
- **2D.** Use one sheet propagator in DFPT, in the subtraction and in EPP §8. Include $Z^*_{zz}$ and out-of-plane quadrupoles.

## 6. State of the art

**Status.** Entries without a mark were confirmed against Giustino's reference list (RMP 89, 015003, 2017) in this revision. Entries marked † are bibliographic data that could not be re-verified in this revision: web search was exhausted, and Crossref, doi.org and the publisher hosts were unreachable. Check them before citing. DOIs of APS papers are built from journal, volume and first page.

Confirmed:

1. S. Baroni, S. de Gironcoli, A. Dal Corso, P. Giannozzi, RMP **73**, 515–562 (2001), [10.1103/RevModPhys.73.515](https://doi.org/10.1103/RevModPhys.73.515). DFPT review.
2. R. M. Pick, M. H. Cohen, R. M. Martin, PRB **1**, 910–920 (1970), [10.1103/PhysRevB.1.910](https://doi.org/10.1103/PhysRevB.1.910). Analytic dielectric matrix; quadrupolar decay in nonpolar insulators.
3. P. Vogl, PRB **13**, 694–704 (1976), [10.1103/PhysRevB.13.694](https://doi.org/10.1103/PhysRevB.13.694). Microscopic multipole theory of the electron–phonon vertex in insulators and semiconductors.
4. X. Gonze, C. Lee, PRB **55**, 10355–10368 (1997), [10.1103/PhysRevB.55.10355](https://doi.org/10.1103/PhysRevB.55.10355). Ewald dipole–dipole split of force constants.
5. C. Verdi, F. Giustino, PRL **115**, 176401 (2015), [10.1103/PhysRevLett.115.176401](https://doi.org/10.1103/PhysRevLett.115.176401). The $g^{\mathcal S}+g^{\mathcal L}$ interpolation.
6. J. Sjakste, N. Vast, M. Calandra, F. Mauri, PRB **92**, 054307 (2015), [10.1103/PhysRevB.92.054307](https://doi.org/10.1103/PhysRevB.92.054307). Wannier interpolation of the polar-optical coupling in GaAs.
7. J. Sjakste, N. Vast, V. Tyuterev, PRL **99**, 236405 (2007), [10.1103/PhysRevLett.99.236405](https://doi.org/10.1103/PhysRevLett.99.236405). Electron–phonon scattering times in GaAs and GaP; hot-electron thermalization.
8. M. Bernardi, D. Vigil-Fowler, C. S. Ong, J. B. Neaton, S. G. Louie, PNAS **112**, 5291–5296 (2015). Ab initio study of hot electrons in GaAs; valley-resolved lifetimes.
9. V. G. Tyuterev, S. V. Obukhov, N. Vast, J. Sjakste, PRB **84**, 035201 (2011), [10.1103/PhysRevB.84.035201](https://doi.org/10.1103/PhysRevB.84.035201). Electron–phonon scattering times in Ge (intervalley).
10. K. Kaasbjerg, K. S. Thygesen, K. W. Jacobsen, PRB **85**, 115317 (2012), [10.1103/PhysRevB.85.115317](https://doi.org/10.1103/PhysRevB.85.115317). Phonon-limited mobility of n-type single-layer MoS₂.

Not re-verified (†):

11. † S. de Gironcoli, S. Baroni, R. Resta, PRL **62**, 2853 (1989), [10.1103/PhysRevLett.62.2853](https://doi.org/10.1103/PhysRevLett.62.2853). Linear-response piezoelectric constants of III–V semiconductors.
12. † R. M. Martin, PRB **5**, 1607 (1972), [10.1103/PhysRevB.5.1607](https://doi.org/10.1103/PhysRevB.5.1607). Piezoelectricity; the clamped-ion term as a sum of first moments.
13. † M. Stengel, PRB **88**, 174106 (2013), [10.1103/PhysRevB.88.174106](https://doi.org/10.1103/PhysRevB.88.174106). Flexoelectricity from DFPT; current response to atomic motion. Nonlocal-pseudopotential implementation: C. E. Dreyer, M. Stengel, D. Vanderbilt, PRB **98**, 075153 (2018), [10.1103/PhysRevB.98.075153](https://doi.org/10.1103/PhysRevB.98.075153).
14. † M. Royo, M. Stengel, PRX **9**, 021050 (2019), [10.1103/PhysRevX.9.021050](https://doi.org/10.1103/PhysRevX.9.021050). Long-wave DFPT for dynamical quadrupoles.
15. † G. Brunin *et al.*, PRL **125**, 136601 (2020), [10.1103/PhysRevLett.125.136601](https://doi.org/10.1103/PhysRevLett.125.136601), and PRB **102**, 094308 (2020), [10.1103/PhysRevB.102.094308](https://doi.org/10.1103/PhysRevB.102.094308). Quadrupolar vertex with the $V^{{\rm Hxc},\mathcal E}$ term; mobilities of Si, GaAs, GaP.
16. † V. A. Jhalani, J.-J. Zhou, J. Park, C. E. Dreyer, M. Bernardi, PRL **125**, 136602 (2020), [10.1103/PhysRevLett.125.136602](https://doi.org/10.1103/PhysRevLett.125.136602). Piezoelectric coupling from quadrupoles in GaN.
17. † J. Park, J.-J. Zhou, V. A. Jhalani, C. E. Dreyer, M. Bernardi, PRB **102**, 125203 (2020), [10.1103/PhysRevB.102.125203](https://doi.org/10.1103/PhysRevB.102.125203). Quadrupole vertex in Si and PbTiO₃.
18. † M. Royo, K. R. Hahn, M. Stengel, PRL **125**, 217602 (2020), [10.1103/PhysRevLett.125.217602](https://doi.org/10.1103/PhysRevLett.125.217602). Multipolar force constants and piezoelectric sound velocities.
19. † M. Royo, M. Stengel, PRX **11**, 041027 (2021), [10.1103/PhysRevX.11.041027](https://doi.org/10.1103/PhysRevX.11.041027). Exact quasi-2D long-range electrostatics.
20. † S. Poncé, M. Royo, M. Gibertini, N. Marzari, M. Stengel, PRL **130**, 166301 (2023), [10.1103/PhysRevLett.130.166301](https://doi.org/10.1103/PhysRevLett.130.166301). Gauge-covariant quadrupoles; Hall mobilities of 2D materials.
21. † S. Poncé, M. Royo, M. Stengel, N. Marzari, M. Gibertini, PRB **107**, 155424 (2023), [10.1103/PhysRevB.107.155424](https://doi.org/10.1103/PhysRevB.107.155424). Long-range contributions to couplings and mobilities in 2D and bulk materials. (Author orders of 20 and 21 as recalled. They agree with [twod.md](twod.md) entry 14; confirm on the journal pages.)
22. † T. Sohier, M. Calandra, F. Mauri, PRB **94**, 085415 (2016), [10.1103/PhysRevB.94.085415](https://doi.org/10.1103/PhysRevB.94.085415). Finite 2D Fröhlich vertex in TMD monolayers.
23. † T. Sohier, M. Calandra, F. Mauri, PRB **96**, 075448 (2017), [10.1103/PhysRevB.96.075448](https://doi.org/10.1103/PhysRevB.96.075448). 2D DFPT with Coulomb truncation and gating.
24. † M. Selig *et al.*, Nat. Commun. **7**, 13279 (2016), [10.1038/ncomms13279](https://doi.org/10.1038/ncomms13279). Phonon-limited exciton linewidths in TMDs.

## 7. First steps

1. **Derive and check (days).** Use checks that can fail. Do them in 3D zincblende directly. The primer's three-atom chain ([§B.5](../../sections/B05-photon-se.html)) is 1D and longitudinal. It cannot hold the off-diagonal $Q^{yz}_x$-type components, the only ones the $T_d$ site symmetry of Si and zincblende allows, nor the antisymmetric part of $\mathcal Q$. E.310 with P1 holds by construction and is not a test.
   - (i) Continuity: $\tilde\bq\cdot\mathbf j(\tilde\bq)/\omega$ from P1 reproduces $-i\tilde\bq\cdot Z^*\cdot u-\frac12\tilde\bq\tilde\bq Q\cdot u$ in EPP's Fourier convention.
   - (ii) With the basis phase, the $O(q)$ acoustic vertex does not depend on the choice of cell.
   - (iii) Martin relation: the acoustic $O(q)$ vertex equals internal strain plus the sublattice-summed $Q$. Compare it with ABINIT's clamped-ion plus relaxed-ion $e_{14}$ for GaAs, not only with the measured value.
   - (iv) For P2, compare $\tilde C^{\mathcal L}$ with an independent long-range force-constant implementation for Si and GaAs. An example is ABINIT anaddb's dipole–quadrupole and quadrupole–quadrupole options; check the variable names in the current docs.
2. **Bulk benchmarks (1–2 weeks).**
   - (a) Reproduce Brunin *et al.*'s published Si and GaAs $g$ and mobilities with ABINIT. This is a regression test that fixes the conventions (units and signs of $Q$, the $V^{{\rm Hxc},\mathcal E}$ term).
   - (b) Then swap in EPP's $g^{\mathcal L,{\rm dyn}}$ with P0 and P1, and measure the difference. Build $g^{\mathcal S}$ by P3 and check its real-space decay. An $O(q^0)$ direction-dependent remainder shows up as a slow $R^{-3}$ tail, instead of the faster decay a complete subtraction leaves.
   - (c) Compare LA[111] with the DFPT $e_{14}$ (clamped-ion plus relaxed) as well as with the measured one.
   - Before choosing pseudopotentials, check the long-wave DFPT restrictions on functional and pseudopotential. They have historically been LDA, norm-conserving, and no nonlinear core correction; check the current version.
3. **$\epsilon$ consistency (1–2 weeks).** Compute the head with your RPA code, with and without local fields and $f^{xc}$, and measure the net residual $1/q$ tail in $g^{\mathcal S}$.
4. **2D (1–2 months).** First reproduce Poncé *et al.*'s published MoS₂ couplings and mobilities. Then treat monolayer MoS₂ (A$_1'$, LO, LA at $K$ and at $|K-\Lambda|$) with truncated DFPT and the §8 propagator, plus P4 for gated samples. Compare with Kaasbjerg *et al.*
5. **Hand-off (months).** Pass $g^{\rm tot}$ to the electron self-energy and exciton–phonon topics.

## 8. Open questions

- Is Route A tractable: does a $\bG\bG'$ RPA+$f^{xc}$ self-energy reproduce $Z^*$, $Q$, $g^{\mathcal S}$?
- Should $Q$ and $g^{\mathcal S}$ be screened dynamically by carriers, as the dipole term is? P4 says how, but $\bar g^{\mathcal S}(\omega)-\bar g^{\mathcal S}(0)$ is dropped.
- Does an $\omega$-dependent $\Delta J$ (P4–P5) keep gauge invariance at finite $\omega$, beyond E.310 at $\omega=0$?
- How good is the transverse local-field kernel from the equal-heads argument in non-cubic and 2D crystals, or is TDCDFT needed?
- Which vertex, bare or screened, belongs on each side of $\bar\Sigma^{\mathcal S}$ when the screening is dynamic?
- Acoustic modes never enter the light cone ($v_s\ll c/\sqrt\epsilon$) **[derived]**. Is finite $c$ ever relevant to short-range physics?
- **[spec]** Does the magnetic part of $\mathcal Q$ matter for chiral phonons in 2D?

## 9. Review responses

All 20 points were checked. EPP was checked for E.134 on p. 26 and for E.350. Giustino was checked for the bibliography, p. 015003-34 (Vogl), p. 015003-25 (G.145 and G.146), p. 015003-46 (§X.A, hot electrons) and p. 015003-44 (§IX.B). The numbers were recomputed: $Z^*|e|/a^2$, $k_{\rm th}$, $r_0Q$. Every point was accepted. Web search was exhausted and Crossref was unreachable, so no new citation could be searched.

1. **P4 double counting (major): accepted.** The algebra was reproduced: the old P4 left a surplus $[W^{\rm e}(\omega)-W^{\rm e}(0)]\Lambda(0)$. P4 now adds only $W^{\rm e}(\omega)[\Lambda(\omega)-\Lambda(0)]$, with the requested sentence and the $\Delta J$ form. P5 is now $-(J+\Delta J)^*D^{\rm e}(J+\Delta J)+[\bar\Sigma^{\mathcal S}(\omega)-\bar\Sigma^{\mathcal S}(0)]$, with $\bar\Sigma^{\mathcal S}$ restricted to $\bar v+f^{xc}$ links. An exactness check was added. The dropped $\bar g^{\mathcal S}(\omega)-\bar g^{\mathcal S}(0)$ and gauge invariance at finite $\omega$ were moved to §8.
2. **Local-field vertex (major): accepted.** Added P0b. §3 now gives the explicit $O(q^0)$ term $q_\beta q_\gamma Z^*V^{{\rm Hxc},E_\gamma}/(\bq\epsilon\bq)$ derived from it, with a mapping caveat for Brunin's factor. P3 now states the remainder if P0b is omitted. The Berry-connection expansion was added to P3.
3. **GaAs attribution (major): accepted.** Giustino p. 015003-46 credits Sjakste–Vast–Tyuterev 2007 and Bernardi *et al.* 2015. Neither text, as Giustino cites it, says "dominates", so the word was dropped. Sjakste 2015 now appears only in §3. Tyuterev 2011 (Ge intervalley, confirmed) was added.
4. **"~30%" (major): accepted.** Replaced in the note and the overview by "significantly; magnitude to be quoted from the paper". The Jhalani and Poncé numbers are now flagged as not re-checked as well.
5. **Verified count (minor): accepted.** §6 is now split into 10 entries confirmed against Giustino's bibliography and 14 marked †. The overview count is 10. The † entries were kept, not dropped. They carry the established results that §§2–4 rely on, the reviewer's fix asks for marking, and dropping them would leave those claims unsourced.
6. **Poncé author orders (minor): accepted.** The entry was split into 20 and 21 with the two orders. They match twod.md, but are still marked † until checked on the journal pages.
7. **Kaasbjerg (minor): accepted.** The wording now includes the Fröhlich interaction. The two deformation-potential values were removed because they could not be checked against the paper's table.
8. **Selig (minor): accepted.** Now reads "$q\approx K$ and $q\approx|K-\Lambda|$". The material names were removed pending a check.
9. **P2 order (minor): accepted.** Verified by power counting and added as the "Order" bullet. The self-check was replaced by an anaddb comparison (step 1 iv).
10. **E.310 check (minor): accepted.** Replaced by the continuity, cell-independence and Martin-relation checks, done in 3D zincblende. The chain's limits are stated.
11. **Toy $\zeta$ (minor): accepted.** Confirmed: next-nearest-neighbour springs in the correction's model are intra-sublattice, so $\zeta=1$. $1.04\times2/\sqrt3=1.20$ C/m². $\zeta$-scaled values were added. The literature $\zeta$ range is marked unverified.
12. **Inconsistent $\epsilon$ (minor): accepted.** Rewritten with both subtraction choices and the opposing local-field and ALDA shifts.
13. **P0 tensor (minor): accepted.** P0a is now stated for $\bar\Pi^{\rm e}_{jj}$, so E.350 yields $e^2\bar\chi_{00}$. The transverse kernel caveat was added.
14. **Route A (minor): accepted.** Added the $v_a^{\rm SR}$ split and the removal of the bare ion–ion term from $k$.
15. **2D power counting (minor): accepted.** Added the derivation, the $r_0Q$ regime and the MoS₂ numbers ($r_0$ unverified). Sohier 2016 was added (†).
16. **Convention mapping (minor): accepted.** New pitfall.
17. **Double screening (minor): accepted.** Reworded. Giustino p. 015003-25 was checked: G.145 has $g^bg^*$, G.146 $|g|^2$.
18. **E.134 (minor): accepted.** The printed $\bA(\bR^0_i)$ was confirmed on EPP p. 26, and the qualifier was added.
19. **Steps 2 and 4 (minor): accepted.** Step 2 is recast as regression test then swap, with long-wave restrictions and DFPT $e_{14}$. Step 4 starts by reproducing Poncé *et al.*
20. **Missing references (minor): accepted.** Vogl 1976 was added (confirmed, G p. 015003-34), along with Pick–Cohen–Martin 1970 (confirmed). Martin 1972 and Sohier 2016 were added as †.
