# Issue report: the missing basis phase (Born dipoles at $\mathbf R_i^0+\mathbf r_a^0$, not at the cell origin)

The paste-ready correction is in `corrections/issue-phase.tex`. The numerical checks are in `test/corr-phase/check_phase.py`; its output is in `test/corr-phase/results.txt` (30 checks, all PASS). The reviewers' script `test/corr-phase-verify/verify_phase.py` was re-run unchanged (17 checks, all PASS).

Revision 2 answers the two skeptic reviews. The responses, one per objection, are in §9.

Sources were verified against the page images:
- EPP pp. 1–38, 39–72 and 73–102 (all of Secs. 1–8). For revision 2, pp. 26–28, 35–39, 48–49, 55–56, 60–70, 75–78, 84–102 were re-read.
- Giustino RMP 89, 015003 (2017), pp. 015003-7, -8, -34 and -35 (Eqs. 12–22, 177–180, footnote 7).

Notes E1, E3, E4, E5 and E7 were cross-read. No equation below is taken from the notes without re-derivation.

---

## 1. Summary

1. **The error.** EPP §3 (Eqs. 121–124) places the dipole $\mathbf p_{ia}=\hat Z^*_a\mathbf Q_{ia}$ of atom $a$ at the cell origin $\mathbf R^0_i$, not at the atomic site $\mathbf R^0_{ia}=\mathbf R^0_i+\mathbf r^0_a$ (E.94). The polar vertex (E.131) therefore lacks the factor $e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}$ inside the atom sum. The same placement recurs in:
   - the Lagrangian (E.134) and everything derived from it: E.141–143, E.183, the Hamiltonians on pp. 38–39, E.186, E.237, E.270;
   - the polarization and ionic densities: E.167, E.241, E.271–272, E.274, E.314–315;
   - the ionic photon self-energy: E.286–311 and the unnumbered displays between them (pp. 61, 62, 65, 66);
   - the ionic $\rho\rho$ correlator: E.352–354 and the displays on pp. 76–77;
   - through the vertex $j_{\alpha\sigma}(\mathbf q)$ of E.295/E.385, every equation of §§6.1–8 (E.389–472). The vertex is written out explicitly, so the phase must be inserted by hand, in the two displays after E.417, in E.418, E.430, E.447, E.448 and in the displays on pp. 95, 96 and 98.
2. **Why it is wrong.** EPP's normal modes use the *cell* phase, $Q_{ia\nu}\propto e_{a\nu\alpha}(\mathbf q)e^{i\mathbf q\mathbf R^0_i}$ (E.7, E.83, E.242), as Giustino does (G Eqs. 14, 20). With this convention, the Fourier component at $\mathbf q+\mathbf G$ of a dipole sitting at $\mathbf R^0_i+\mathbf r^0_a$ carries $e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}$. §2 keeps the factor (E.101); §3 drops it.
3. **The fix.** The minimal change is $V^{\mathcal L}_a(\mathbf q)\to e^{-i\mathbf q\mathbf r^0_a}V^{\mathcal L}_a(\mathbf q)$ (124′), exactly parallel to E.101. Then E.125–E.130 keep their form, and E.131′ equals Giustino's Eq. 180 term by term (up to the $e$ inside $Z^*$, the labels, and $\mathbf G\neq-\mathbf q$). From §4 on, every $\mathbf A(\mathbf R^0_i)$ and $\delta(\mathbf r-\mathbf R^0_i)$, in numbered and unnumbered displays alike, becomes $\mathbf A(\mathbf R^0_{ia})$ and $\delta(\mathbf r-\mathbf R^0_{ia})$. The ionic vertex becomes $\mathbf G$-resolved, $j_{\alpha\sigma}(\mathbf q+\mathbf G)$ (295′), and the ionic self-energy becomes a genuine $\mathbf G\mathbf G'$ matrix (297′, 304′). The sum rules E.309–E.311 survive as matrix identities.
4. **A precondition.** The cancellation that removes the spurious acoustic term (131b), and the cell-invariance law (131a), need the eigenvectors that solve the $e^{i\mathbf q\mathbf R^0_i}$ equations of motion, i.e. the sign of E.15/E.35 corrected (Z.2 errata item 3). With E.15 as printed, the phase would *double* the acoustic error. The two fixes must be adopted together; the `.tex` now says so in the main text.
5. **What does not change.** Everything at $\mathbf q\to0$ and $\mathbf G=0$ for optical modes:
   - the Fröhlich limit of E.131;
   - the Huang equations of §4.2.1 (E.144–165);
   - E.406–E.410 at $\mathbf q=0$;
   - the LO–TO splitting (the $q\to0$ limit of E.454);
   - the polariton region (relative change $5\times10^{-9}$ in GaAs).
6. **What changes.**
   - Finite $q$ in multi-atom cells. In GaAs, in the 3D force-constant model, $|g_{\rm LO}|^2$ changes by 0.5–0.8 % at 0.1 Å⁻¹ and by 12–22 % at 0.5 Å⁻¹, up or down depending on direction and cell. The cruder two-atom estimate (q = 0 eigenvectors) gives 0.2–1.5 % and 4–33 % and is quoted only as an illustration.
   - All $\mathbf G\neq0$ (local-field) components.
   - Acoustic modes. The printed vertex equals the corrected one plus a spurious, cell-dependent term with $|g|^2\propto1/q$, present even in centrosymmetric crystals. In GaAs the spurious term alone is equivalent to a piezoelectric constant of 0.6–1.8 C/m² along [111] (depending on the cell) and 1.04 C/m² along [100]. That is comparable to or larger than the internal-strain coupling of the corrected vertex (1.2 C/m² along [111] in the toy model, 0 along [100]). The printed total can be larger (1.8) or smaller (0.6 C/m²) than the corrected one. The measured $(2/\sqrt3)|e_{14}|\approx0.19$ C/m² needs, in addition, the clamped-ion term carried by the dynamical quadrupoles.
7. **Correction to the task statement.** The task wrote that "$\sum_aZ^*_ae^{-i\mathbf q\mathbf r_a}\approx-i\sum_aZ^*_a(\mathbf q\mathbf r_a)$ gives an $O(q)$ dipolar coupling" for acoustic modes. In EPP's (and Giustino's) cell-phase convention the acoustic eigenvector itself carries $e^{+i\mathbf q\mathbf r^0_a}$.
   - The *printed* (phase-less) vertex therefore contains the spurious term $+\frac{i}{\sqrt M}\sum_a(\mathbf q\mathbf r^0_a)\hat Z^*_a\mathbf u$.
   - The *corrected* vertex cancels that term exactly, and what remains at $O(q)$ is physical (internal strain).
   - The expression $-i\sum_aZ^*_a(\mathbf q\mathbf r_a)$ arises only if one multiplies the explicit phase by an atom-phase (real) eigenvector, i.e. if conventions are mixed. It also arises without the phase if the eigenvectors of the printed E.15 are used.
   - The magnitude estimate is right; the attribution is reversed. Verified numerically (§5.4, checks B and C3).

---

## 2. Conventions used (verified on the pages)

| Item | EPP | Where |
|---|---|---|
| Atom position | $\mathbf R^0_{ia}=\mathbf R^0_i+\mathbf r^0_a$ | E.1, E.2, E.94 (pp. 1, 17) |
| Lattice Fourier transform | $Q_{ia\mu}=N^{-1/2}\sum_{\mathbf k}e^{i\mathbf k\mathbf R^0_i}Q_{\mathbf ka\mu}$ (**cell phase**) | E.7 (p. 2), E.62, E.83 (p. 15), E.242 (p. 49) |
| Continuous Fourier transform | $f(\mathbf r)=\int\frac{d\mathbf q}{(2\pi)^3}e^{i\mathbf q\mathbf r}f(\mathbf q)$ | E.101 (p. 18), E.117 (p. 23) |
| Pseudopotential basis phase | $V_a(\mathbf q')=e^{-i\mathbf q'\mathbf r^0_a}v_a(\mathbf q')$ | E.101 (p. 18) |
| Units | SI in §3 ("here we use SI units", p. 23); Gaussian from §4 ($\mathbf E=-\frac1c\dot{\mathbf A}$, $\int(E^2-B^2)/8\pi$) | pp. 23, 25–26 |
| Born charges | $Z^*_{a\mu\nu}=\delta p_\mu/\delta R_{a\nu}$, dimension of charge: $Z^*_{\rm EPP}=\lvert e\rvert Z^*_{\rm G}$ | E.110 (p. 22) |
| Electron charge in §3 | $-e$, $e>0$ | before E.123 (p. 24) |
| Static dielectric constant | $\epsilon_\mu(0)$ | E.165 (p. 35) |
| Volumes | $\Omega$ = crystal, $\Omega_{\rm uc}=\Omega/N$; Giustino's $\Omega$ = EPP $\Omega_{\rm uc}$ | p. 20; G p. 7 |

**Eigenvector remark (Z.2 errata item 3), now part of the correction.**
- E.15/E.35 define $\tilde\Omega^2(\mathbf k)$ with $e^{+i\mathbf k\boldsymbol\delta}$, $\boldsymbol\delta=\mathbf R^0_i-\mathbf R^0_j$. That is Giustino's $D^{\rm dm}(-\mathbf k)=D^{\rm dm}(\mathbf k)^*$.
- The step after E.43 (p. 9) uses $\tilde\Omega^2_{b\nu a\mu}(\mathbf k)=\tilde\Omega^2_{a\mu b\nu}(\mathbf k)$, which is true only for real $\tilde\Omega^2$.
- E.83 is consistent only if $e_{a\mu\alpha}(\mathbf k)$ are the eigenvectors of $(m_am_b)^{-1/2}\sum_{\boldsymbol\delta}e^{-i\mathbf k\boldsymbol\delta}k(\boldsymbol\delta)$, i.e. those that solve the equations of motion for the pattern $e^{i\mathbf k\mathbf R^0_i}$. Checked numerically (section B of the script): residual $2\times10^{-16}$ for those eigenvectors, against 0.61 for eigenvectors of the printed E.35.
- The *sign* of the basis phase is fixed by E.83 together with E.117, not by E.15. Its *cancellation* against the eigenvector, (131a) and (131b), needs the correct eigenvectors. With the printed E.15 the eigenvectors are conjugated and carry $e^{-i\mathbf q\mathbf r^0_a}$ for an acoustic wave, so the phase doubles the spurious term (script check B(iii): $2.0000\times$ the spurious value). Revision 1 put this in a footnote and claimed that the correction "uses only the cell phase of E.83", which was overstated. The `.tex` now states the dependence in the main text.

---

## 3. Derivation

### 3.1 Why the basis phase must appear explicitly

Mode $(\mathbf q,\alpha)$ displaces atom $a$ of cell $i$ by $\propto e_{a\nu\alpha}(\mathbf q)e^{i\mathbf q\mathbf R^0_i}$ (E.83). The dipole $p_{ia\mu}=Z^*_{a\mu\nu}Q_{ia\nu}$ sits at $\mathbf R^0_i+\mathbf r^0_a$, so the polarization density is $\mathbf P(\mathbf r)=\sum_{ia}\mathbf p_{ia}\delta(\mathbf r-\mathbf R^0_{ia})$. With the transform of E.117 its Fourier component is

$$\mathbf P(\mathbf k)=\sum_{ia}\mathbf p_{ia}e^{-i\mathbf k(\mathbf R^0_i+\mathbf r^0_a)}\propto\sum_a\hat Z^*_a\mathbf e_{a\alpha}(\mathbf q)e^{-i\mathbf k\mathbf r^0_a}\sum_ie^{i(\mathbf q-\mathbf k)\mathbf R^0_i}=N\sum_{\mathbf G}\delta_{\mathbf k,\mathbf q+\mathbf G}\sum_a\hat Z^*_a\mathbf e_{a\alpha}(\mathbf q)e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}.$$

Hence $e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}$ multiplies $e_{a\nu\alpha}(\mathbf q)$ in every vertex (checked, script section A).

**Cell-choice invariance (131a).** Suppose atom $a$ is assigned to a neighbouring cell, $\mathbf r^0_a\to\mathbf r^0_a+\mathbf R_L$. The physical displacement of each atom is unchanged, so the cell-phase eigenvector component becomes $e_{a\nu\alpha}\to e^{i\mathbf q\mathbf R_L}e_{a\nu\alpha}$. Only the product $e_{a\nu\alpha}(\mathbf q)e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}$ is invariant, so the printed vertex depends on an arbitrary bookkeeping choice. Script check C1, GaAs model with As at $(1,1,1)a/4$ or at $(1,-1,-1)a/4$:
- corrected vertices: identical to $2\times10^{-13}$;
- printed vertices: differ by up to $0.15\,e/\sqrt{\rm amu}$ (the LO vertex itself is $0.34$).

The 1D analogue is check B(ii): reassigning atom 2 of a chain to the previous cell leaves $\arg(e_2/e_1)-\mathbf q\mathbf r_2$ unchanged to $10^{-10}$.

**Atom-phase convention.** With $Q_{ia\nu}\propto\tilde e_{a\nu\alpha}e^{i\mathbf q\mathbf R^0_{ia}}$ (so $\tilde e_a=e^{-i\mathbf q\mathbf r^0_a}e_a$), the factor reduces to $e^{-i\mathbf G\mathbf r^0_a}$. E.131 as printed is correct only in this convention and only at $\mathbf G=0$.

### 3.2 Corrected §3

Start from E.120, the potential of a dipole at the origin, which is unchanged: $\phi(\mathbf q)=-\frac{i}{\epsilon_0}\frac{\mathbf p\mathbf q}{\mathbf q\hat\epsilon^\infty\mathbf q}$. Translating the dipole to $\mathbf R^0_{ia}$ multiplies $\phi(\mathbf q)$ by $e^{-i\mathbf q\mathbf R^0_{ia}}$. This gives (121′) and (122′), with $e^{i\mathbf q(\mathbf r-\mathbf R^0_i-\mathbf r^0_a)}$ in the integrand; I also fixed the left label of (121), $\phi_{ia\mu}\to\phi_{ia\nu}$, as in Z.2. Multiplying by the electron charge $-e$ and summing over $i$ gives (123′):

$$\Delta V^{\mathcal L}_{\rm cr}(\mathbf r)=i\sum_{ia\mu\nu}Z^*_{a\mu\nu}Q_{ia\nu}\int\frac{d\mathbf q}{(2\pi)^3}\frac{e\,q_\mu}{\epsilon_0\mathbf q\hat\epsilon^\infty\mathbf q}e^{i\mathbf q(\mathbf r-\mathbf R^0_i-\mathbf r^0_a)}.$$

Define $V^{\mathcal L}_a(\mathbf q)=e^{-i\mathbf q\mathbf r^0_a}e^2/(\epsilon_0\mathbf q\hat\epsilon^\infty\mathbf q)$ (124′). With E.125, $Q^{\mathcal L}=-\hat Z^*Q/e$, we get

$$-i\sum Q^{\mathcal L}_{ia\mu}\int\frac{d\mathbf q}{(2\pi)^3}q_\mu V^{\mathcal L}_a(\mathbf q)e^{i\mathbf q(\mathbf r-\mathbf R^0_i)}=\frac{ie}{\epsilon_0}\sum Z^*Q\int\frac{d\mathbf q}{(2\pi)^3}\frac{q_\mu}{\mathbf q\hat\epsilon^\infty\mathbf q}e^{i\mathbf q(\mathbf r-\mathbf R^0_i-\mathbf r^0_a)},$$

which reproduces (123′). So E.128 holds as printed with (124′), and its identity with E.102 now includes the basis phase of E.101.

The steps E.102→E.107 are then reused verbatim. Inserting E.83, the lattice sum $\sum_ie^{i(\mathbf q-\mathbf q')\mathbf R^0_i}=N\sum_{\mathbf G}\delta_{\mathbf q',\mathbf q+\mathbf G}$ sets $\mathbf q'=\mathbf q+\mathbf G$, so $V^{\mathcal L}_a(\mathbf q+\mathbf G)$ carries $e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}$. E.130 keeps its form. The prefactor $-\frac{i}{\Omega_{\rm uc}}\cdot(-\frac1e)\cdot\frac{e^2}{\epsilon_0}=\frac{4\pi i}{\Omega_{\rm uc}}\frac{e}{4\pi\epsilon_0}$ then gives

$$g^{\mathcal L}_{nm\alpha}(\mathbf k,\mathbf q)=\frac{4\pi i}{\Omega_{\rm uc}}\frac{e}{4\pi\epsilon_0}\sum_{\mathbf G}\sum_{a\mu\nu}\sqrt{\frac{\hbar}{2Nm_a\omega_\alpha(\mathbf q)}}\frac{(q_\mu+G_\mu)Z^*_{a\mu\nu}e_{a\nu\alpha}(\mathbf q)e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}}{(\mathbf q+\mathbf G)\hat\epsilon^\infty(\mathbf q+\mathbf G)}\rho_{n\mathbf k,m\mathbf k-\mathbf q}(-\mathbf G).\qquad(131')$$

### 3.3 Agreement with Giustino Eq. (180)

G Eq. 180 reads

$$g^{\mathcal L}_{mn\nu}(\mathbf k,\mathbf q)=i\frac{4\pi}{\Omega}\frac{e^2}{4\pi\varepsilon_0}\sum_\kappa\left(\frac{\hbar}{2N_pM_\kappa\omega_{\mathbf q\nu}}\right)^{1/2}\sum_{\mathbf G\neq-\mathbf q}\frac{(\mathbf q+\mathbf G)\cdot\mathbf Z^*_\kappa\cdot\mathbf e_{\kappa\nu}(\mathbf q)}{(\mathbf q+\mathbf G)\cdot\boldsymbol\epsilon^\infty\cdot(\mathbf q+\mathbf G)}\langle\psi_{m\mathbf k+\mathbf q}|e^{i(\mathbf q+\mathbf G)\cdot(\mathbf r-\boldsymbol\tau_\kappa)}|\psi_{n\mathbf k}\rangle.$$

Map onto EPP with:
- $\Omega\to\Omega_{\rm uc}$, $N_p\to N$, $M_\kappa\to m_a$, $\boldsymbol\tau_\kappa\to\mathbf r^0_a$;
- $e\,Z^*_{\rm G}\to Z^*_{\rm EPP}$, so $\frac{e^2}{4\pi\varepsilon_0}Z^*_{\rm G}=\frac{e}{4\pi\epsilon_0}Z^*_{\rm EPP}$;
- final/initial labels $(m\mathbf k+\mathbf q,n\mathbf k)\to(n\mathbf k,m\mathbf k-\mathbf q)$, so that $\langle\psi_{m\mathbf k+\mathbf q}|e^{i(\mathbf q+\mathbf G)\mathbf r}|\psi_{n\mathbf k}\rangle\to\rho_{n\mathbf k,m\mathbf k-\mathbf q}(-\mathbf G)$ (E.106).

The phase $e^{-i(\mathbf q+\mathbf G)\cdot\boldsymbol\tau_\kappa}$ then matches. Giustino's phonon convention is the cell phase (G Eq. 20, $\Delta\tau_{\kappa\alpha p}\propto e^{i\mathbf q\cdot\mathbf R_p}e_{\kappa\alpha,\nu}(\mathbf q)$; G Eq. 14 with $e^{i\mathbf q\cdot\mathbf R_p}$), as in EPP, so no further conversion is needed. The only remaining difference is the exclusion $\mathbf G\neq-\mathbf q$.

Footnote 7 (G p. 35) says only that VG Eq. 4 lacks the factor and that it "needs to be retained in order to correctly describe the acoustic modes near $\mathbf q=0$", and that results do not depend on the cutoff $a$. It does not give the cell-phase argument; the `.tex` now attaches the citation to the statement that VG Eq. 4 lacked the factor.

### 3.4 Acoustic modes, and the correct reading of footnote 7

In the cell-phase convention a long acoustic wave displaces every atom as $\mathbf u\,e^{i\mathbf q(\mathbf R^0_i+\mathbf r^0_a)}$, plus an $O(q)$ internal relaxation. So

$$e_{a\nu\alpha}(\mathbf q)=e^{i\mathbf q\mathbf r^0_a}\left[\sqrt{m_a/M}\,u_\nu+i\,w_{a\nu}(\mathbf q)+O(q^2)\right],\qquad M=\sum_am_a,$$

with $\mathbf w_a=O(q)$ the internal strain, which does not depend on the cell choice. Script checks B(i)–(ii):
- symmetric chain ($k_1=k_2$, atom 2 at $d=a/2$, so $\mathbf w_a=0$): $\arg(e_2/e_1)=qd$ to $10^{-9}$;
- general chain ($k_1\neq k_2$, $d=1.3$): $\arg(e_2/e_1)=qak_2/(k_1+k_2)$, the uniform-stress value, so $\arg-qd=-1.57\times10^{-4}$ at $q=10^{-3}$ is the internal-strain part; it is the same for both cell assignments.

Then (131b) follows:
- printed minus corrected: $\sum_a\hat Z^*_a\mathbf e_a/\sqrt{m_a}-\sum_a\hat Z^*_a\mathbf e_ae^{-i\mathbf q\mathbf r^0_a}/\sqrt{m_a}=\frac{i}{\sqrt M}\sum_a(\mathbf q\mathbf r^0_a)\hat Z^*_a\mathbf u+O(q^2)$;
- corrected: $\sum_a\hat Z^*_a\mathbf e_ae^{-i\mathbf q\mathbf r^0_a}/\sqrt{m_a}=\frac{1}{\sqrt M}(\sum_a\hat Z^*_a)\mathbf u+i\sum_a\hat Z^*_a\mathbf w_a/\sqrt{m_a}+O(q^2)=i\sum_a\hat Z^*_a\mathbf w_a/\sqrt{m_a}+O(q^2)$.

(Revision 1 wrote the printed form as "$=\frac{i}{\sqrt M}\sum_a(\mathbf q\mathbf r^0_a)\hat Z^*_a\mathbf u+O(q)$", which says nothing because both terms are $O(q)$. The difference form is what the script checks: the remainder scales as $q^{2.00}$ for [111] in both cells and for [123], check C3.)

**Properties of the spurious term.** It is origin-independent, because $\sum_a\hat Z^*_a=0$. It does depend on the assignment of atoms to cells. With $\omega_{\rm ac}\propto q$ it gives $g^{\mathcal L}\propto\frac1q\cdot q\cdot q^{-1/2}$, i.e. $|g|^2\propto1/q$: the scaling of a piezoelectric coupling, with an apparent piezoelectric constant $e_{\rm fake}=|\sum_aZ^*_a\hat{\mathbf q}\mathbf r^0_a|/\Omega_{\rm uc}$. It is nonzero even in centrosymmetric crystals, where the physical term $\sum_a\hat Z^*_a\mathbf w_a$ vanishes (inversion pairs have equal $\hat Z^*$ and opposite $\mathbf w$). For rocksalt with Cl at $(\pm\frac12,0,0)a$, the sign of $e_{\rm fake}$ flips with the cell choice.

**Script checks C3.**
- LA[100] in zincblende is piezo-inactive. The corrected vertex scales as $q^{2.00}$; the printed one has slope $Z^*(a/4)/\sqrt M=0.2432\,e$Å$/\sqrt{\rm amu}$ for both cells.
- LA[111]: the corrected slope is $0.2809$ in the toy model for both cells (physical, from internal strain). The printed slope is $0.1404$ (cell A) or $0.4213$ (cell B). In cell A the printed total is *smaller* than the corrected one: $0.1404=|0.2809-0.4213|$, while in cell B $0.4213=0.2809+0.1404$.
- At $q=10^{-4}$ Å⁻¹, printed minus corrected equals $i\sum_aZ^*_a(\mathbf q\mathbf r_a)/\sqrt M$ to $10^{-3}$ relative; at $q=10^{-3}$ and $2\times10^{-3}$ the remainder of (131b) scales as $q^2$.
- Slopes converted to piezoelectric constants, $s\sqrt M/\Omega_{\rm uc}$: corrected 1.20 C/m² ([111]) and 0 ([100]); printed totals 0.60 ([111] A), 1.80 ([111] B), 1.04 C/m² ([100]).

**The physical $O(q)$ term.** It comes from internal strain plus the clamped-ion (electronic) piezoelectric response. The dipole vertex, corrected or not, contains only the internal-strain part. The clamped-ion part requires the dynamical quadrupoles, which enter the long-range potential at the same order. In III–V semiconductors the two parts are of opposite sign and largely cancel, which is why the measured $|e_{14}|$ of GaAs is small (literature, from memory: S. de Gironcoli, S. Baroni, R. Resta, PRL 62, 2853 (1989)). Quadrupole literature (flagged, citations from memory; authors should verify):
- G. Brunin *et al.*, Phys. Rev. Lett. 125, 136601 (2020);
- V. A. Jhalani *et al.*, Phys. Rev. Lett. 125, 136602 (2020);
- J. Park *et al.*, Phys. Rev. B 102, 125203 (2020);
- background: M. Royo and M. Stengel, Phys. Rev. X 9, 021050 (2019).

So the corrected E.131′ has the right structure at $O(q^0)$ and the right *dipole* part at $O(q)$, but it is not complete at $O(q)$. This is what Giustino's footnote 7 means by "needs to be retained in order to correctly describe the acoustic modes near $\mathbf q=0$".

### 3.5 §4 (Gaussian)

The coupling is $\frac1c\sum_{ia}\dot{\mathbf p}_{ia}\mathbf A(\mathbf R^0_{ia})$ (134′). The Euler–Lagrange equation E.139 then gives $\partial L/\partial\dot Q_{ia\mu}=m_a\dot Q_{ia\mu}+\frac1c\sum_\nu Z^*_{a\nu\mu}A_\nu(\mathbf R^0_{ia})$. Its time derivative involves $\partial_tA$ at the fixed equilibrium site (dipole approximation within each atom), so (141′)–(143′) have $\dot A_\nu$ and $E_\nu$ at $\mathbf R^0_{ia}$. In §4.2.1 the field is taken uniform over the cell ($|\mathbf q|\ll1/a$), so E.144–E.165 are unchanged.

The `.tex` now states one rule for all of §§4–6, numbered and unnumbered displays alike: $A(\mathbf R^0_i)\to A(\mathbf R^0_{ia})$, $\delta(\mathbf r-\mathbf R^0_i)\to\delta(\mathbf r-\mathbf R^0_{ia})$, with the list of displays it touches (pp. 26, 28, 35, 38, 39, 48, 49, 56, 60–62, 64–66, 68, 76). It displays (134′), (143′), (183′), (186′); (141′), (142′), (237′) and (270′) follow from the rule and are not displayed. Further:
- E.167 becomes (167′) = (241′) = (314′).
- E.175–E.180 are unchanged in form; $\mathbf P$ is now (167′).
- E.235 is unchanged in form, with $\rho_{\rm ion}=-\nabla\mathbf P$ built from (241′).

### 3.6 §5: ionic photon self-energy, $\mathbf G$-resolved

E.272′ is $j^{\rm ion}_\alpha=\sum_{ia\mu}\frac{1}{m_a}Z^*_{a\alpha\mu}\pi_{a\mu}(\mathbf r)\delta(\mathbf r-\mathbf R^0_{ia})$, written, as in the manuscript, with the continuous field $\pi_{a\mu}(\mathbf r)$. (Revision 1 wrote $\pi_{ia\mu}$; equivalent, but not the manuscript's form.) The continuous field of E.274 must read $e^{i\mathbf q(\mathbf r-\mathbf r^0_a)}$ (274′), so that $\pi_{a\mu}(\mathbf R^0_{ia})$ equals E.243. The contraction of E.287–E.289 is unchanged; the factor $e^{i\mathbf q(\mathbf r-\mathbf r')}$ becomes $e^{i\mathbf q(\mathbf R^0_i-\mathbf R^0_j)}$ on the support of the delta functions. Taking $n_{\mathbf q}=0$ gives (290′), now displayed in the `.tex`.

The lattice sums are
$$\sum_ie^{i\mathbf q\mathbf R^0_i}\delta(\mathbf r-\mathbf R^0_{ia})=\Omega_{\rm uc}^{-1}\sum_{\mathbf G}e^{i(\mathbf q+\mathbf G)(\mathbf r-\mathbf r^0_a)},\qquad\sum_je^{-i\mathbf q\mathbf R^0_j}\delta(\mathbf r'-\mathbf R^0_{jb})=\Omega_{\rm uc}^{-1}\sum_{\mathbf G'}e^{-i(\mathbf q+\mathbf G')(\mathbf r'-\mathbf r^0_b)}.$$
(Revision 1 quoted the sum without $e^{i\mathbf q\mathbf R^0_i}$; the result was right, the step incomplete.) With $(1/N)\sum_{\mathbf q}\to\Omega_{\rm uc}\int_{BZ}d\mathbf q/(2\pi)^3$ the prefactor is $\frac{1}{c\Omega_{\rm uc}}$ as in E.293, giving (293′), now displayed. Comparing with E.294 yields (295′)–(297′):

$$j_{\alpha\sigma}(\mathbf q+\mathbf G)=\sum_{a\mu}\frac{Z^*_{a\alpha\mu}}{\sqrt{m_a}}e_{a\mu\sigma}(\mathbf q)e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a},\qquad\Pi^{\rm ion,para}_{\alpha\beta,\mathbf G\mathbf G'}=-\frac{1}{c\Omega_{\rm uc}}\sum_\sigma\frac{\omega^2_\sigma}{\omega^2-\omega^2_\sigma+i\eta}j_{\alpha\sigma}(\mathbf q+\mathbf G)j^*_{\beta\sigma}(\mathbf q+\mathbf G').$$

In the cell-phase convention $e_a(\mathbf q+\mathbf G)=e_a(\mathbf q)$, because $\tilde\Omega^2$ is $\mathbf G$-periodic, so $j_{\alpha\sigma}(\mathbf q+\mathbf G)$ is a well-defined function of the full wave vector $\mathbf q+\mathbf G$. It is not periodic: the basis phase supplies the structure factor of each sublattice.

**Diamagnetic part.** Applying the same lattice sum to E.302′ gives $\frac{1}{\Omega_{\rm uc}}\sum_{\mathbf G}e^{i\mathbf G(\mathbf r-\mathbf r^0_a)}\int\frac{d\mathbf q'}{(2\pi)^3}e^{i\mathbf q'(\mathbf r-\mathbf r')}$. Decompose $\mathbf q'=\mathbf q+\mathbf G'$ and relabel $\mathbf G+\mathbf G'\to\mathbf G$ (the implicit relabelling on p. 66). This gives (304′), $\Pi^{\rm ion,dia}_{\alpha\beta,\mathbf G\mathbf G'}=-\frac{1}{\Omega_{\rm uc}c}\sum_{a\mu}\frac{Z^*_{a\alpha\mu}Z^*_{a\beta\mu}}{m_a}e^{-i(\mathbf G-\mathbf G')\mathbf r^0_a}$.

**Completeness, (306)–(311).** E.306 and E.307 become (297′) at $\omega=0$ and (307′). E.308 writes the vertex product out with $e_{a\mu\sigma}e^*_{b\nu\sigma}$ and no phases, so it does *not* keep its printed form: (308′) carries $e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}e^{+i(\mathbf q+\mathbf G')\mathbf r^0_b}$. (Revision 1 listed E.308 as "same printed form"; corrected.) In (309′), $\sum_\sigma e_{a\mu\sigma}e^*_{b\nu\sigma}=\delta_{ab}\delta_{\mu\nu}$ (E.28) keeps only $a=b$, so $e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}e^{+i(\mathbf q+\mathbf G')\mathbf r^0_a}=e^{-i(\mathbf G-\mathbf G')\mathbf r^0_a}$. Hence E.310 and E.311 hold as $\mathbf G\mathbf G'$ matrix identities once the subscripts are added (script check C4: $|\Pi^{\rm dia}+\Pi^{\rm para}(0)|<10^{-16}$ for 16 $(\mathbf G,\mathbf G')$ pairs). The diagonal $\mathbf G=\mathbf G'$ of (304′) is E.304 as printed.

**2D.** (300′) and (305′) are now displayed, with in-plane $\mathbf G$, in-plane positions $\boldsymbol\rho^0_a$ and heights $z_a$: the para part carries $\sum_{ab}j^a_{\alpha\sigma}(\mathbf Q+\mathbf G)j^{b*}_{\beta\sigma}(\mathbf Q+\mathbf G')\delta(z-z_a)\delta(z'-z_b)$, the dia part $e^{-i(\mathbf G-\mathbf G')\boldsymbol\rho^0_a}\delta(z-z_a)\delta(z'-z_a)$ (from $\delta(z-z_a)\delta(z-z')$). At $\omega=0$, completeness again makes dia $=-$para. For a planar sheet, $z_a=0$ and the atom sums collapse to the manuscript's $\delta(z)\delta(z')$ form.

### 3.7 §6

E.352 becomes the matrix (352′). The derivation on pp. 76–77, with $\delta(\mathbf r-\mathbf R^0_{ia})\delta(\mathbf r'-\mathbf R^0_{jb})$ in the two displays after E.353 and in the substitution E.354, yields (354a): $\Pi^{\rm ion}_{\rho\rho,\mathbf G\mathbf G'}=-\frac{c}{\omega^2}(q+G)_\mu(q+G')_\nu\Pi^{\rm ion}_{jj\mu\nu,\mathbf G\mathbf G'}$, using $\partial_\mu\to i(q+G)_\mu$ and $\partial'_\nu\to-i(q+G')_\nu$. So the connection formula E.351 still holds, and E.341–E.350 are general and unchanged.

### 3.8 §§6.1–8

From E.380 on, the manuscript drops $\mathbf G\mathbf G'$. That step was justified for the ions by "the ionic self-energy is a constant in $\mathbf G\mathbf G'$ space" (p. 63), which is no longer true. For the electrons it was justified by "CLFs are mostly neglected here" (p. 95). Within the $\mathbf G=\mathbf G'=0$ (macroscopic) approximation:
- the equations that contain $j_{\mu\sigma}$ only as a symbol keep their printed form with the vertex (385′): E.389, E.391–393, E.399–417, E.429, E.449–454 and, in 2D, E.455–472;
- the equations that write the vertex out need the phase inserted by hand, $e^{-i\mathbf q\mathbf r^0_a}$ next to $e_{a\mu\sigma}$ and $e^{+i\mathbf q\mathbf r^0_a}$ next to $e^*_{a\mu\sigma}$. These are the two displays after E.417 ((417a), (417b) in the `.tex`), E.418 ((418′)), the p. 95 display ((428b)), E.430, E.447, E.448 ((430′), (447′), (448′)) and the boxed p. 98 display ((448b)). Revision 1 missed E.418 and the two p. 93 displays, and wrongly listed E.430, E.447 and E.448 as "same printed form". As a result its (428b), derived "after using (428) in (418)", did not follow from a corrected (418). The chain is now consistent: (417a) → (417b) [with E.375] → (418′) [$c\to\infty$] → (428b) [with E.428]; and (430′) → (447′) [with E.446] → (448′) → (448b) [with E.428]. Each step was re-derived against the printed displays on pp. 93, 95, 96 and 98.

With local fields kept, the phonon self-energy is the matrix (399′).

**Signs and the factor $e$ in the electrostatic limits.**
- (448b) $=-[(428b)]^*\,q^2/(\mathbf q\hat\epsilon^{\rm e}\mathbf q)$ exactly, because $(4\pi i)^*=-4\pi i$; at $\hat\epsilon^{\rm e}=1$ this is *minus* the complex conjugate. Revision 1 called it "the complex conjugate"; corrected in the `.tex` and the p. 99 replacement text. Checked with random complex data and a complex symmetric $\hat\epsilon^{\rm e}$ (script section G), and by the reviewers' `verify_phase.py` section E.
- With the printed prefactor, (428b) $=g_{\rm G}/|e|$, where $g_{\rm G}$ is the $\mathbf G=0$, $\epsilon^\infty=1$ term of G.180 in Gaussian units, phase included. With the factor $e$ of `issue-efactor.tex` restored ($e$ the signed electron charge), $e\cdot$(428b) $=(e/|e|)g_{\rm G}=-g_{\rm G}$, which is that correction's (428a) with the phase. Likewise $e\cdot$(448b) $=-(e/|e|)[g_{\rm G}]^*q^2/(\mathbf q\hat\epsilon^{\rm e}\mathbf q)$, its (448a) with the phase (script section G, all three relations to $10^{-12}$).
- The `.tex` therefore says: "up to the factor $e$ and an overall sign (the electron-charge convention, cf. (428a), (448a) of the e-factor correction), (428b) is the $\mathbf G=0$, $\hat\epsilon^\infty=1$ term of G.180, and (448b) is its complex conjugate with $\hat\epsilon^\infty\to\hat\epsilon^{\rm e}$". This matches the wording of `issue-efactor.tex`, whose combined relation is $\lim g^{\rm dyn}=-[\lim g^{\rm dyn,0}]^*q^2/(\mathbf q\hat\epsilon^{\rm e}\mathbf q)$.

**Merging with `corrections/issue-efactor.tex`.** Both corrections replace the same displays on pp. 95, 96 and 98, and the same sentence on p. 99. The `.tex` now gives the merged displays (p95″), (p96″), (p98″): factor $e$ in the prefactor and phase in the atom summand. Three places in `issue-efactor.tex` must change if both are adopted (identified by content, since that file is being edited in parallel):
1. the paragraph introducing (p96′): "with $\boldsymbol\tau_\kappa=0$ as for our dipoles at $\mathbf R^0_i$ (314)" → "with $\boldsymbol\tau_\kappa=\mathbf r^0_a$, the positions of our dipoles (314′)";
2. the sentence after (448a): "($\mathbf G=0$, $\boldsymbol\tau_\kappa=0$)" → "($\mathbf G=0$)";
3. check (iv) of its Justification: "$\mathbf G=0$, $\boldsymbol\tau_\kappa=0$" → "$\mathbf G=0$, phase $e^{-i\mathbf q\boldsymbol\tau_\kappa}$ included".

Its 304 meV comparison is unaffected, because both sides of (428a)/(448a) carry the same phase. Tags do not clash: this correction uses (417a), (417b), (418′), (428b), (430′), (447′), (448′), (448b), (p95″), (p96″), (p98″); that one uses (p94′), (422′), (428′), (p95′), (p96′), (428a), (p98′), (448a), (p101′), (234′), (235′).

**E.454 and Gonze–Lee.** E.454 keeps its form with (385′). Its $\mathbf G$-resolved matrix form (454a) is the dipole–dipole dynamical matrix of Gonze & Lee (1997), with the structure factor $e^{i(\mathbf q+\mathbf G)(\mathbf r^0_a-\mathbf r^0_b)}$. Its projection on mode $\sigma$ at $\mathbf G=0$ is $e^*_\sigma\Delta\tilde\Omega^2e_\sigma=\frac{W}{\Omega_{\rm uc}}|\mathbf q\,j_\sigma|^2$, i.e. E.454. The sign of the structure factor was verified independently by an Ewald test (§5.3).

---

## 4. Equation-by-equation list

**Corrected and displayed (primed in the `.tex`):**
- §3: 121′, 122′, 123′, 124′, 131′.
- §4: 134′, 143′, 167′, 183′, 186′.
- §5: 241′, 271′, 272′, 274′, 290′, 293′, 295′, 296′, 297′, 299′, 300′, 302′, 304′, 305′, 307′, 308′, 309′.
- §6: 314′, 315′, 352′.
- §§6.1–8: 385′ (the vertex displays after E.385), 399′, 418′, 430′, 447′, 448′.

**Corrected by the stated rule, not displayed:**
- 141, 142, 237, 270 ($A$, $\dot A$ at $\mathbf R^0_{ia}$);
- 286, 289, 291, 303, 354 and the 2D substitutions on pp. 64, 66 ($\delta$ at $\mathbf R^0_{ia}$, Fourier substitution with $\mathbf r^0_a$);
- 306 ((297′) at $\omega=0$); 310, 311 (add $\mathbf G\mathbf G'$);
- unnumbered: the Hamiltonian after E.184 (p. 38), the two Hamiltonians before E.186 (p. 39), the display after E.288 (p. 61), the displays after E.303 (p. 65), the displays after E.353 (p. 76) and their successors on p. 77 (result: (354a)).

**New or replacing unnumbered displays (lettered or double-primed):**
- 131a: cell-choice invariance.
- 131b: acoustic expansion (difference form).
- 354a: $\Pi^{\rm ion}_{\rho\rho,\mathbf G\mathbf G'}$.
- 417a, 417b: the two displays after E.417 (p. 93).
- 428b: replaces the unnumbered display on p. 95.
- 448b: replaces the boxed unnumbered display on p. 98.
- p95″, p96″, p98″: the displays on pp. 95, 96, 98 with both this and the e-factor correction.
- 454a: $\mathbf G$-resolved E.454 = Gonze–Lee.

**Same printed form, new content through the corrected vertex or density:**
- 125–130 (with 124′);
- 169, 175–180 and 235 (with 167′/241′);
- 273 (with 271′, 272′);
- 313, 351;
- 386–417 and 429 (symbol $j$ only), 449–454;
- 455–472 (2D, in-plane phase; symbol $j$ only).

**Unchanged:**
- 1–120, including E.120, the potential of a single dipole;
- 144–165 (uniform field over a cell);
- 166, 168, 170–174, 181–182, 184–185, 187–234, 236 (photon field), 238–240, 242–269, 275–285 (electronic and photon parts), 316–350;
- 406–410 at $\mathbf q=0$ (E.409–410 explicitly);
- the $q\to0$ limits of 131′ (Fröhlich) and 454 (LO–TO) for optical modes.

---

## 5. Checks

### 5.1 Dimensions (Gaussian)
- The correction factor $e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}$ is dimensionless, so no printed dimension changes.
- E.131′ in Gaussian units: $[eZ^*]=$ statC$^2=$ erg·cm. Then $\Omega_{\rm uc}^{-1}$ (cm$^{-3}$) × $\sqrt{\hbar/m\omega}$ (cm) × $q/q^2$ (cm) gives erg.
- $[j]=$ statC/$\sqrt{\rm g}$, as printed.
- Script check F: the SI evaluation ($Z^*_{\rm G}$ dimensionless, $e^2/4\pi\epsilon_0$) and the Gaussian evaluation ($e^2=14.3996$ eV Å, $Z^*_{\rm EPP}=|e|Z^*_{\rm G}$) of $\sqrt N|g^{\mathcal L}|$ at $q=0.1$ Å⁻¹ agree: 303.991 vs 303.990 meV.

### 5.2 Limits
- **$q\to0$, $\mathbf G=0$.** Only relative phases $e^{-i\mathbf q(\mathbf r^0_a-\mathbf r^0_b)}$ are observable; they go to 1. The optical-triplet strength $\sum_{\rm opt}|\hat{\mathbf q}j|^2=Z^{*2}/\mu=0.11865$, identical for both forms (check C2). E.131′ equals Fröhlich (G Eq. 7) exactly when $Z^*$ is the LO–TO-implied 2.177 (check E). With the input $Z^*=2.07$ the two differ by 5 %: 0.0304 vs 0.0320 eV Å⁻¹, an inconsistency of the inputs, not of the formula. LO–TO: $\hbar\omega_p=44.56$ meV, $\hbar\omega_{\rm LO}=35.93$ meV (measured 36.2).
- **$c\to\infty$.** The chains (417a)→(428b) and (430′)→(448b) are re-derived from the printed displays; signs and the factor $e$ as in §3.8 (script section G).
- **Polariton (Huang).** $q_c=\omega_{\rm TO}\sqrt{\epsilon(0)}/c=6.06\times10^{-5}$ Å⁻¹ and $q_cd=1.5\times10^{-4}$, so $|j_{\rm LO}|^2$ changes by $5.5\times10^{-9}$. E.157–E.165 and E.410 are unaffected. (Revision 1 wrote $\epsilon_0$, which in the same fragment denotes the SI vacuum permittivity; now $\epsilon(0)$ as in E.165.)

### 5.3 Gonze–Lee: Ewald test of the $\mathbf G$-resolved structure (script section D)
- **Setup.** Energy per cell of a modulated dipole lattice $\mathbf p_{ia}=\mathbf p_ae^{i\mathbf q\mathbf R^0_i}$ at $\mathbf R^0_{ia}$ (GaAs fcc geometry, random complex $\mathbf p_a$). The reciprocal sum uses the structure factor of (295′); the real-space part uses the erfc-screened dipole tensor; the Gaussian self term is subtracted.
- **Correct sign $e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}$.** $E=-0.070411\,e^2$/Å, independent of the Ewald parameter ($\Lambda=0.45,0.6,0.8$ Å⁻¹) to $3\times10^{-14}$.
- **Wrong sign.** The result varies by 34 %.
- **Dipoles at the cell origin.** The on-site ($a=b$) term needs the self-term removal in both forms. With coincident dipoles the intracell $a\neq b$ pair diverges as well, so the $\mathbf G$ sum of the dipole–dipole energy or dynamical matrix built from the printed E.297/E.454 has no limit. The vertex $g^{\mathcal L}$ itself does not diverge: coincident dipoles just add into one dipole per cell. (Revision 1's Justification bullet implied a divergence of the vertex; restricted now.)

### 5.4 Script summary
`python3 test/corr-phase/check_phase.py` runs 30 checks, all PASS, in sections:
- A: Fourier sign;
- B: E.15 remark; acoustic eigenvector phase for (i) a symmetric chain with atom 2 explicitly at $d=a/2$, (ii) a general chain with $d=1.3$ (internal strain, cell invariance), (iii) printed E.15 plus phase doubles the spurious term;
- C: zincblende model, covering cell invariance, $q\to0$, acoustic scaling, (131b) remainder $O(q^2)$, piezoelectric equivalents, and E.310′;
- D: Ewald;
- E: unchanged limits;
- F: units;
- G: (448b) $=-[(428b)]^*q^2/\mathbf q\hat\epsilon^{\rm e}\mathbf q$; relations to G.180 with the factor $e$.

`python3 test/corr-phase-verify/verify_phase.py` (the reviewers' script, unchanged): 17 checks, all PASS.

---

## 6. Numbers (GaAs-like; $d=|\mathbf r^0_{\rm As}-\mathbf r^0_{\rm Ga}|=2.448$ Å)

Inputs: $\hbar\omega_{\rm TO}=33.3$ meV, $\hbar\omega_{\rm LO}=36.2$ meV, $\epsilon^\infty=10.9$, $Z^*=2.07\,e$, $\Omega_{\rm uc}=45.2$ Å³, masses 69.72/74.92 amu ($\mu=36.1$), $a=5.653$ Å, $v_{\rm LA}=5.40$ km/s ([111]) and 4.73 km/s ([100]).

**(a) 3D toy model (main estimate).**
- Model: zincblende with NN central springs $\alpha=7.185$ eV/Å² (fit to $\omega_{\rm TO}$) and NNN springs $\beta=0.451$ eV/Å² (fit to $v_{\rm LA}[100]$).
- Eigenvectors at finite $q$; $\mathbf G=0$; vertices in $e/\sqrt{\rm amu}$, with $Z^*/\sqrt\mu=0.3445$.

| direction, cell | q | LO corr. | LO printed | ratio² | LA corr. | LA printed |
|---|---|---|---|---|---|---|
| [111], (1,1,1)a/4 | 0.1 | 0.3433 | 0.3442 | 1.005 | 0.0282 | 0.0139 |
| [111], (1,1,1)a/4 | 0.5 | 0.3105 | 0.3398 | 1.198 | 0.1492 | 0.0562 |
| [111], (1,−1,−1)a/4 | 0.1 | 0.3433 | 0.3419 | 0.992 | 0.0282 | 0.0421 |
| [111], (1,−1,−1)a/4 | 0.5 | 0.3105 | 0.2739 | 0.778 | 0.1492 | 0.2089 |
| [100], (1,1,1)a/4 | 0.1 | 0.3445 | 0.3436 | 0.995 | 0.0001 | 0.0243 |
| [100], (1,1,1)a/4 | 0.5 | 0.3444 | 0.3232 | 0.881 | 0.0032 | 0.1191 |

LO error of the printed form: 0.5–0.8 % at 0.1 Å⁻¹, 11.9–22.2 % at 0.5 Å⁻¹, of either sign (+20 % for [111] cell A, −22 % for cell B, −12 % for [100]).

**(b) Two-atom estimate with leading-order eigenvectors (illustration only).** $\theta=\mathbf q\,(\mathbf r^0_{\rm As}-\mathbf r^0_{\rm Ga})$. It uses the $\mathbf q=0$ eigenvectors with the cell phase and so misses the $O(q)$ internal-strain phase of the optical eigenvector; for [111] cell A at 0.5 Å⁻¹ it even has the wrong sign (0.670 against 1.198 in the 3D model).

Formulas:
- LO amplitude ratio (printed/correct): $|m_{\rm As}+m_{\rm Ga}e^{i\theta}|/M$.
- Spurious LA vertex: $(2Z^*/\sqrt M)|\sin\frac\theta2|$; the correct value is 0 at this order.

| direction, cell | q (Å⁻¹) | θ | LO $\lvert g\rvert^2$ ratio | $\lvert j_{\rm LA}\rvert/\lvert j_{\rm LO}\rvert$ | $\lvert g_{\rm LA}\rvert/\lvert g_{\rm LO}\rvert$ |
|---|---|---|---|---|---|
| [111], As at (1,1,1)a/4 | 0.1 | 0.245 | 0.985 | 0.122 | 0.39 |
| [111], As at (1,1,1)a/4 | 0.5 | 1.224 | 0.670 | 0.574 | 0.82 |
| [111], As at (1,−1,−1)a/4 | 0.1 | −0.082 | 0.998 | 0.041 | 0.13 |
| [111], As at (1,−1,−1)a/4 | 0.5 | −0.408 | 0.959 | 0.202 | 0.29 |
| [100], either cell | 0.1 | 0.141 | 0.995 | 0.071 | 0.24 |
| [100], either cell | 0.5 | 0.707 | 0.880 | 0.346 | 0.53 |

How to read the tables:
- The corrected columns do not depend on the cell; the printed ones do.
- Along [111] the $O(q)$ internal-strain phase of the optical eigenvector can push the printed LO value either up or down. That is why the 3D model gives $|g|^2$ ratios of 1.005 and 0.992 at 0.1 Å⁻¹, where the leading-order estimate gave 0.985.
- The `.tex` quotes the 3D model (0.5–0.8 %, 12–22 %) and gives the two-atom values (0.2–1.5 %, 4–33 %) only as an illustration. Revision 1 merged both into "4–33 %".

**(c) Acoustic coupling expressed as piezoelectric constants.** $e_{\rm fake}=|\sum_aZ^*_a\hat{\mathbf q}\mathbf r^0_a|/\Omega_{\rm uc}$ for the spurious term alone; a slope $s$ maps to $s\sqrt M/\Omega_{\rm uc}$.

| | [111] cell A | [111] cell B | [100] |
|---|---|---|---|
| spurious term alone | 1.80 | 0.60 | 1.04 |
| corrected vertex (internal strain, toy model) | 1.20 | 1.20 | 0 |
| printed total (toy model) | 0.60 | 1.80 | 1.04 |
| measured, $(2/\sqrt3)\lvert e_{14}\rvert$ (literature $\lvert e_{14}\rvert=0.16$) | 0.185 | 0.185 | 0 |

(all in C/m²). The spurious term is comparable to or larger than the internal-strain coupling of the corrected vertex. The printed total can be larger or smaller than the corrected one, and it is nonzero along [100], where the true coupling vanishes. Neither dipole vertex reaches the measured value: that needs the clamped-ion term (dynamical quadrupoles), which largely cancels the internal-strain term in III–V semiconductors. Revision 1 compared "0.6–1.8 C/m² against 0.19 C/m²", which set the spurious term against the full measured $e_{14}$ and suggested that the corrected vertex reproduces 0.19 C/m². It does not.

---

## 7. What the primer should change

| Section | Change |
|---|---|
| **B01-foundations** (§B.1) | Replace the "Note on the source: basis phase" callout (line 149–150) with a *correction* callout. Display (E.124′) and (E.131′); tag them as corrected and cite `corrections/issue-phase.tex`, since the guide forbids inventing source numbers, so mark them "corrected". Add the one-line cell-reassignment argument (131a), and say that it and (131b) need the E.15 sign fix (Z.2 item 3). Add an interactive "cell-choice test": a GaAs-like chain (or zincblende along [111]/[100]) with toggles "printed E.131 / corrected E.131′" and "As at +d / at d−a". Plot $\lvert g_{\rm LO}(q)\rvert^2$ and $\lvert g_{\rm LA}(q)\rvert^2$. The corrected curves must not move when the cell toggles. The printed LA curve diverges as $q^{-1/2}$, while the corrected one vanishes along [100]. |
| **A04-polar** (§A.4) | Rewrite the acoustic bullet under (G.180) (line 136) precisely. In Giustino's (Maradudin–Vosko, cell-phase) convention the acoustic eigenvector carries $e^{i\mathbf q\boldsymbol\tau_\kappa}$; the explicit $e^{-i(\mathbf q+\mathbf G)\boldsymbol\tau_\kappa}$ cancels it. Without it, the vertex acquires a spurious, cell-dependent term, $\lvert g\rvert^2\propto1/q$, present even in centrosymmetric crystals. The physical $O(q)$ coupling is the internal-strain dipole plus the clamped-ion term carried by dynamical quadrupoles (Brunin/Jhalani/Park 2020, literature). Update the source callout (line 201–202) to point to the correction; footnote 7 should be cited for "VG Eq. 4 lacked the factor", not for the cell-phase reason. If the 1D lab is meant to "reproduce this logic", include the basis phase there too. |
| **B02-lagrangian** (§B.2) | In the "Dipole approximation" bullet (line 97), state the corrected coupling (E.134′) with $\mathbf A(\mathbf R^0_{ia})$ instead of "This drops the basis phase". Keep "harmless at photon wave vectors", quantified as $q_cd\approx10^{-4}$ and a relative change of $5\times10^{-9}$ for GaAs. |
| **B03-dielectric** (§4.2.1) | Note that (E.143′) puts $\mathbf E$ at $\mathbf R^0_{ia}$ and that the long-wavelength results E.146–E.165 are unchanged. |
| **B04-fields** | Replace $\bR^0_i$ by $\bR^0_{ia}$ in line 18 ((E.167), $\delta(\br-\bR^0_{ia})$), line 34 ((E.183′)), line 41 (Hamiltonian, p. 39), line 43 ((E.186′)) and lines 149–150 ((E.237) and (E.241′)). $\rho_{\rm ion}=-\nabla\cdot\bP$ acquires the structure factors at $\mathbf G\neq0$. In the double-counting callout (line 161), keep the double-counting content and update only the closing clause "where the basis phase that (E.167) and (E.241) also drop is discussed" to "where the basis phase, dropped in the printed (E.167) and (E.241) and restored in (E.167′), (E.241′), is discussed". |
| **B05-photon-se** | The largest change. Display (E.295′), (E.297′) and (E.299′) as $\mathbf G\mathbf G'$ matrices with $j(\mathbf q+\mathbf G)$. In the derivation, the lattice sum $\sum_ie^{i\mathbf q\mathbf R^0_i}\delta(\mathbf r-\mathbf R^0_{ia})=\Omega_{\rm uc}^{-1}\sum_{\mathbf G}e^{i(\mathbf q+\mathbf G)(\mathbf r-\mathbf r^0_a)}$ now gives structure factors. Replace "The $\mathbf G$-independence is an artefact … only $\mathbf G=\mathbf G'=0$ should be used" with the replacement text of the `.tex` (sublattice dipole gratings, reliable only at small $\lvert\mathbf q+\mathbf G\rvert$). Display (E.304′). In the completeness subsection, show (E.308′) with both structure factors and then that the phases cancel for $a=b$, so E.309–E.311 hold as matrices. In **viz-b-photon-se-sumrule**, give the three ions explicit positions and use (295′). Add a toggle "printed (E.295) / corrected (E.295′)": the stacked sum (E.309) is the same for both, but the per-branch weights change, and in the printed form the acoustic branch acquires spurious weight $\propto q^2$. |
| **B06-propagator** | P in line 18 becomes (E.314′). Add (E.354a), $\Pi^{\rm ion}_{\rho\rho,\mathbf G\mathbf G'}$. |
| **B07-phonon-se** | In the comparison table, change "point dipoles at cell origins" to "point dipoles at atomic sites $\mathbf R^0_{ia}$ (E.314′)". State that the vertex is (E.385′) and that E.409–E.410 ($\mathbf q=0$) are unchanged. |
| **B08-frohlich** | Replace "The basis phase of G Eq. 180 is missing" (line 87) with the corrected displays. Show (E.418′) as the step before the p. 95 display. With the e-factor correction already adopted in the primer (E p.95′, E p.98′, line 75), use the merged forms (p95″), (p98″) of `issue-phase.tex`. They equal, up to the overall sign $e/\lvert e\rvert$ of (428a), the $\mathbf G=0$ term of G.180 including $e^{-i\mathbf q\mathbf r^0_a}$, and its complex conjugate with $e^{+i\mathbf q\mathbf r^0_a}$ (line 81 keeps "complex conjugated, for real $\hat\epsilon^{\rm e}$", which is right with the factor $e$; without it the relation to (428b) is *minus* the conjugate). Any vertex plot at finite $q$ should use the phase. |
| **B09-electrostatic** | With (E.385′), E.454 equals the Gonze–Lee $\mathbf G=0$ term at all $q$, not only "at small $q$" (line 111). Display (E.454a) as the $\mathbf G$-resolved form (the Gonze–Lee dipole–dipole matrix after Ewald regularization) and mention the Ewald $\Lambda$-independence check. Give the finite-$q$ size in GaAs from the 3D model (0.5–0.8 % at 0.1 Å⁻¹, 12–22 % at 0.5 Å⁻¹). Change "Born-charge point dipoles at the cell origins" in the bridge callout to "at the atomic sites". |
| **B10-2d** | Add the in-plane basis phase $e^{-i(\mathbf Q+\mathbf G)\boldsymbol\rho^0_a}$ in (300′)/(305′); for h-BN, B and N sit at inequivalent in-plane sites. For buckled or multilayer sheets (MoS₂ with S at $\pm z$), add the $\delta(z-z_a)$ form displayed in (300′)/(305′). |
| **X02-dictionary** | Add rows: basis vector EPP $\mathbf r^0_a$ ↔ G $\boldsymbol\tau_\kappa$; phase convention cell phase in both (E.7/E.83 ↔ G.14/G.20; mind Z.2 item 3); ionic vertex (E.295′) ↔ the Verdi–Giustino numerator with structure factor. |
| **Z02-errata** | Item 14 (line 31): list the corrected equations (from E.121′ to E.454a, as in §4 above, including E.308′, E.418′, E.430′, E.447′, E.448′) and state the effect precisely: spurious acoustic $\lvert g\rvert^2\propto1/q$, equivalent to 0.6–1.8 C/m² in GaAs, comparable to the internal-strain coupling and not to be compared with the measured $e_{14}$; finite-$q$ LO error 0.5–0.8 % at 0.1 Å⁻¹ and 12–22 % at 0.5 Å⁻¹; $\mathbf G\neq0$. Link `corrections/issue-phase.tex`. Cross-reference item 3 (E.15 sign), which must be fixed together with item 14, and item 29 (point-dipole self-energy), which (E.454a) needs for its regularization. Update the intro sentence (line 9) accordingly. |

---

## 8. Open questions and caveats

1. **Task framing on the acoustic sign.** See §1, item 7, and §3.4. The spurious $O(q)$ term belongs to the printed vertex; the corrected vertex removes it. The task's expression $-i\sum_aZ^*_a(\mathbf q\mathbf r_a)$ mixes the cell-phase factor with an atom-phase eigenvector. The magnitude estimate is right.
2. **E.15 sign (Z.2 item 3).** The correction assumes that $e(\mathbf q)$ in E.83 solves the $e^{i\mathbf q\mathbf R^0_i}$ equations of motion (Giustino's G.14). If the authors kept E.15 and E.26 literally, the eigenvectors would come out conjugated, and the phase would double the acoustic error (check B(iii)). The `.tex` now states in the main text that item 3 must be fixed together with this one.
3. **Dynamical quadrupoles and the measured $e_{14}$.** The literature citations (Brunin, Jhalani, Park 2020; Royo–Stengel 2019; de Gironcoli–Baroni–Resta 1989 for the cancellation between clamped-ion and internal-strain terms in III–V semiconductors) are from memory and flagged; the authors should verify the bibliographic details. No quadrupole formula is given in the `.tex`. The toy-model internal-strain value (1.2 C/m²) depends on the force-constant model (central NN + NNN springs) and is illustrative only.
4. **Toy-model numbers.** The 3D model is a toy: it fixes the internal strain, hence the physical LA[111] slope and the [111] LO numbers. It is now the main estimate; the two-atom estimate is labelled as an illustration. The GaAs $e_{14}$ is a literature value.
5. **Regularizing the $\mathbf G$ sum.** In (454a) and in any $\mathbf G\neq0$ use of (297′) and (304′), the on-site ($a=b$) point-dipole terms diverge. They need a consistent prescription: Ewald self-term removal, the cutoff $e^{-a|\mathbf q+\mathbf G|^2}$, and avoiding double counting with the short-range force constants $k^{\mathcal S}$ and with the self-energies kept in E.235 (Z.2 item 29). This is beyond the present correction.
6. **The approximation in §§6.1–8 should be stated.** Dropping $\mathbf G\mathbf G'$ there is now a genuine approximation for the ions as well as for the electrons: it neglects ionic local fields. The manuscript should say so where it drops $\mathbf G\mathbf G'$ (p. 64 and p. 95).
7. **Overlap with the e-factor correction.** The merged displays (p95″), (p96″), (p98″) and the three edits to `issue-efactor.tex` are in §3.8. That file was being edited while this revision was made; the edits are identified by content, not line numbers. I did not edit it.
8. **2D buckled sheets.** (300′)/(305′) now display the $\delta(z-z_a)$ form. For MoS₂-type layers the out-of-plane positions also enter $D^{\rm e}(\mathbf Q,\omega;z_a,z_b)$ inside the atom sums of E.458, which the manuscript's (456) evaluates only at $z=z'=0$.

---

## 9. Review responses

Two skeptics raised 20 objections; several overlap. Each was checked against the PDF pages named and against the derivation. 19 are accepted and fixed; 1 is accepted in part.

| # | Objection (short) | Verdict | What changed |
|---|---|---|---|
| 1 | (major) E.418 and the two unnumbered displays after E.417 write $Z^*_{a\beta\mu}e_{a\mu\sigma}(\mathbf q)$ and need $e^{-i\mathbf q\mathbf r^0_a}$; "p. 93, before (417)" is the wrong place | **Accepted.** Verified on p. 93: E.417 is generic in $j_{\beta\sigma}(\mathbf q)$; the explicit vertex appears in the two displays after it and in E.418. | Added (417a), (417b), (418′). The vertex list now runs to (417) only as a symbol, and the replacement location reads "p. 93, after (417), where (295) is inserted". The chain (417a)→(417b)→(418′)→(428b) was re-derived. |
| 2 | E.308 writes the vertex product out without phases, so it is not "same printed form" | **Accepted.** Verified on p. 67. | Added (308′) with $e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}e^{+i(\mathbf q+\mathbf G')\mathbf r^0_b}$, plus (307′); removed 306–308 from "same printed form" (306 and 307 now listed under the rule). |
| 3 | (131b) left equation is vacuous: both sides $O(q)$ | **Accepted.** | (131b) now gives the difference printed − corrected $=\frac{i}{\sqrt M}\sum_a(\mathbf q\mathbf r^0_a)\hat Z^*_a\mathbf u+O(q^2)$. New check C3: remainder exponent 2.00 in three cases. |
| 4 | Headline "0.6–1.8 vs 0.19 C/m²" compares the spurious term with the full measured $e_{14}$; the printed total can be smaller than the corrected one; the corrected vertex is ≈1.2 C/m² | **Accepted.** Reproduced: corrected 1.20 C/m², printed totals 0.60/1.80 ([111]) and 1.04 ([100]). | Paragraph after (131b) rewritten as the reviewer proposed. Comparison table in §6(c). |
| 5 | (428b)/(448b) identities hold only up to sign and factor $e$; inconsistent with `issue-efactor.tex` (428a) | **Accepted.** | `.tex` now states (448b) $=-[(428b)]^*q^2/\mathbf q\hat\epsilon^{\rm e}\mathbf q$ exactly, and "up to the factor $e$ and an overall sign (electron-charge convention, cf. (428a), (448a))" for the relation to G.180. Script section G checks all three relations. |
| 6 | $\epsilon_0$ used for the static dielectric constant | **Accepted.** | $\epsilon(0)$, as in E.165, with a gloss. |
| 7 | Summary lists equations as corrected that are given only as rules; unnumbered displays not mentioned; the lattice sum lacks $e^{i\mathbf q\mathbf R_i}$ | **Accepted.** Verified pp. 26, 28, 38–39, 48, 56, 60–66, 76–77. | Explicit rule for all displays of §§4–6, numbered or not, with the list of pages. Displayed (186′), (290′), (293′), (300′), (305′). The lattice sum now reads $\sum_ie^{i\mathbf q\mathbf R^0_i}\delta(\mathbf r-\mathbf R^0_{ia})=\Omega_{\rm uc}^{-1}\sum_{\mathbf G}e^{i(\mathbf q+\mathbf G)(\mathbf r-\mathbf r^0_a)}$. §4 separates "displayed" from "by the rule". |
| 8 | "4–33 % at 0.5 Å⁻¹" mixes two estimates that disagree | **Accepted.** | The 3D model is the main estimate (0.5–0.8 %, 12–22 %, either sign); the two-atom values are labelled as an illustration. The script prints both ranges. |
| 9 | Footnote "uses only the cell phase of E.83" is overstated: (131a)/(131b) need the E.15 fix | **Accepted.** New check B(iii): with the printed E.15 the phase gives 2.0000× the spurious term. | Main-text paragraph: the sign follows from (83)+(117); the cancellation needs the corrected E.15/E.35; the fixes must be adopted together. |
| 10 | Script section B: $d=1.3$ unused; the test only shows $\arg(e_2/e_1)=qa/2$ | **Accepted.** | B(i) now places atom 2 explicitly at $d=a/2$. New B(ii) uses $d=1.3$ with $k_1\neq k_2$: $\arg(e_2/e_1)=qak_2/(k_1+k_2)$, $\arg-qd$ is the internal strain, and it is cell-independent. |
| 11 | Same as 1, plus E.430, E.447, E.448 carry the expanded vertex | **Accepted.** Verified pp. 96 and 98. | Added (430′), (447′), (448′) with $e^*e^{+i\mathbf q\mathbf r^0_a}$; general insertion rule stated; list of "symbol-only" equations corrected. |
| 12 | "(448b) is its complex conjugate" has the wrong sign | **Accepted** (same as 5). | "minus the complex conjugate" in the `.tex` sentence and in the p. 99 replacement text. |
| 13 | $A(\mathbf R^0_i)$ on pp. 38–39 and $\delta(\mathbf r-\mathbf R^0_i)$ in E.291–292, E.303, E.354, p. 76 not covered; E.308 | **Accepted** (overlaps 2, 7). | Covered by the rule and its explicit list; (308′) added. |
| 14 | Piezo comparison overstated; corrected vertex gives only the internal-strain part | **Accepted** (same as 4). The reviewer's "realistic ζ ≈ 0.5 gives ≈ 0.6 C/m²" is from memory and was not adopted; the `.tex` quotes only the toy-model value and the literature statement that the clamped-ion term largely cancels the internal-strain term. | As 4. |
| 15 | $\epsilon_0$ clash; "$e^{-i\mathbf q\mathbf r_a}=1+O(q\lvert\mathbf r_a-\mathbf r_b\rvert)$" true only for the relative phase | **Accepted.** | $\epsilon(0)$; "only the relative phases $e^{-i\mathbf q(\mathbf r^0_a-\mathbf r^0_b)}$ are observable; a common factor drops out of $\lvert g\rvert^2$ and $jj^*$". |
| 16 | (131b): make the internal strain explicit | **Accepted** (same as 3). | $e_a=e^{i\mathbf q\mathbf r^0_a}[\sqrt{m_a/M}\mathbf u+i\mathbf w_a+O(q^2)]$ displayed; the right equation of (131b) shows $i\sum_a\hat Z^*_a\mathbf w_a/\sqrt{m_a}$. |
| 17 | Justification bullet 3 overstated: the vertex does not diverge | **Accepted.** | Bullet restricted to the ionic self-energy and the dipole–dipole matrix (454a): on-site terms need regularization in both forms, and with coincident dipoles the $a\neq b$ terms diverge too; the vertex stays finite. Script note in section D updated. |
| 18 | Footnote 7 cited for the cell-phase reason, which it does not give | **Accepted.** Verified on G p. 35. | Citation moved to "Verdi and Giustino's Eq. (4) lacked it (Giustino's footnote 7, which notes that it is needed for the acoustic modes near $\mathbf q=0$)". |
| 19 | Coordination with `issue-efactor.tex`: three places with $\boldsymbol\tau_\kappa=0$; merged displays for pp. 95, 96, 98 | **Accepted.** In the current file the three places are the (p96′) paragraph, the sentence after (448a) and check (iv) (line numbers shift; that file is being edited in parallel). | Merged (p95″), (p96″), (p98″) added; the three edits named in the `.tex` and in §3.8. That file was not edited here. |
| 20 | Primer B04 pointer wrong: line 161 is the double-counting callout; the edits are at 34, 41, 43, 149–150 | **Accepted in part.** The missing edits (and line 18, E.167) are real and now listed. But line 161 does contain a clause on the dropped phase ("where the basis phase that (E.167) and (E.241) also drop is discussed"), so it needs an update, not deletion. The "drops the basis phase" sentence in B02 line 97 is covered by the B02 item. | B04 item rewritten: lines 18, 34, 41, 43, 149–150 to $\bR^0_{ia}$; line 161 keeps the double-counting content and only its closing clause is updated. |

---

## Note (25 Oct 2026): merged deliverable

The paste-ready correction for the authors is now `corrections/EPP_corrections.tex`. It merges the three
corrections (basis phase, E.410, factor e) into one standalone document and supersedes `issue-*.tex`, which
carry a header saying so. This report stays as the supporting derivation. Changes relative to this report:

- **One tag convention.** Every corrected or new display is `equation*`/`align*` with an explicit `\tag`, so
  nothing steps the equation counter. (n′) replaces numbered equation n, (pN′) replaces an unnumbered display on
  page N, (n a) is new. Pasting rule (stated once, merged §1): `\stepcounter{equation}` after each replaced
  numbered equation.
- **Merged checks.** `test/corr-merge/check_merge.py` (24 checks, all pass) tests the merged forms. The LaTeX was
  validated without TeX by `test/corr-merge/validate_tex.py` (braces, environments, math delimiters, `\left/\right`,
  tags, labels). Every math region was rendered with MathJax: `extract_math.py` builds `test/corr-merge-frag/`,
  then `build/assemble.py` and `build/check.js --eager`, and `inspect_math.js` reads the MathML, because check.js
  misses undefined macros and unclosed braces. Result: 889 regions, 63 tags, 0 problems.
- **Specific to this correction.**
  - Tags renamed: (417a), (417b) → (p93a′), (p93b′); (399′) → (399a), a new equation after E.399 on p. 89,
    because E.399 itself keeps its form at **G** = **G′** = 0.
  - (428b), (448b), (p95″), (p96″), (p98″) are retired. The final pp. 95/96/98 displays carry both the phase and
    the factor e and are tagged (p95′), (p96′), (p98′) (merged §5.1). This resolves the overlap noted in §3.8 and §8
    item 7. The three τ_κ edits to `issue-efactor.tex` listed in §3.8 are made in the merged file.
  - Combined tags split for pasting: (167′); (241′, 314′); (271′, 315′); (272′). E.315 is the ionic density tensor,
    the twin of E.271; E.272 is the ionic current.
  - New merged display (237′): the rule A(**R**⁰_ia) plus eρ_e in the Coulomb term.
  - New check M5: (410d) for cq ≫ √ε ω, written in the atom basis, equals the **G** = 0 term of (454a), with
    deviation 14.7 (ω/cq)².
  - Order after E.454 on p. 99: the E.410 paragraph first, then the paragraph with (454a) (merged §5.5).
  - Added a sentence after E.380 (p. 85) saying that **G** = **G′** = 0 now also neglects ionic local fields (§8 item 6).
