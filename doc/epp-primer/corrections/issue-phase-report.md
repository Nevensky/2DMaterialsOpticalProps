# Issue report: the missing basis phase (Born dipoles at $\mathbf R_i^0+\mathbf r_a^0$, not at the cell origin)

The paste-ready correction is in `corrections/issue-phase.tex`. The numerical checks are in `test/corr-phase/check_phase.py`; its output is in `test/corr-phase/results.txt` (22 checks, all PASS).

Sources were verified against the page images:
- EPP pp. 1–38, 39–72 and 73–102 (all of Secs. 1–8).
- Giustino RMP 89, 015003 (2017), pp. 015003-7, -8, -34 and -35 (Eqs. 12–22, 177–180, footnote 7).

Notes E1, E3, E4, E5 and E7 were cross-read. No equation below is taken from the notes without re-derivation.

---

## 1. Summary

1. **The error.** EPP §3 (Eqs. 121–124) places the dipole $\mathbf p_{ia}=\hat Z^*_a\mathbf Q_{ia}$ of atom $a$ at the cell origin $\mathbf R^0_i$, not at the atomic site $\mathbf R^0_{ia}=\mathbf R^0_i+\mathbf r^0_a$ (E.94). The polar vertex (E.131) therefore lacks the factor $e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}$ inside the atom sum. The same placement recurs in:
   - the Lagrangian (E.134) and everything derived from it: E.141–143, E.183, E.186, E.237, E.270;
   - the polarization and ionic densities: E.167, E.241, E.271–272, E.274, E.314–315;
   - the ionic photon self-energy: E.286–311;
   - the ionic $\rho\rho$ correlator: E.352–354;
   - through the vertex $j_{\alpha\sigma}(\mathbf q)$ of E.295/E.385, every equation of §§6.1–8 (E.389–472).
2. **Why it is wrong.** EPP's normal modes use the *cell* phase, $Q_{ia\nu}\propto e_{a\nu\alpha}(\mathbf q)e^{i\mathbf q\mathbf R^0_i}$ (E.7, E.83, E.242), as Giustino does (G Eqs. 14, 20). With this convention, the Fourier component at $\mathbf q+\mathbf G$ of a dipole sitting at $\mathbf R^0_i+\mathbf r^0_a$ carries $e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}$. §2 keeps the factor (E.101); §3 drops it.
3. **The fix.** The minimal change is $V^{\mathcal L}_a(\mathbf q)\to e^{-i\mathbf q\mathbf r^0_a}V^{\mathcal L}_a(\mathbf q)$ (124′), exactly parallel to E.101. Then E.125–E.130 keep their form, and E.131′ equals Giustino's Eq. 180 term by term (up to the $e$ inside $Z^*$, the labels, and $\mathbf G\neq-\mathbf q$). From §4 on, every $\mathbf A(\mathbf R^0_i)$ and $\delta(\mathbf r-\mathbf R^0_i)$ becomes $\mathbf A(\mathbf R^0_{ia})$ and $\delta(\mathbf r-\mathbf R^0_{ia})$. The ionic vertex becomes $\mathbf G$-resolved, $j_{\alpha\sigma}(\mathbf q+\mathbf G)$ (295′), and the ionic self-energy becomes a genuine $\mathbf G\mathbf G'$ matrix (297′, 304′). The sum rules E.309–E.311 survive as matrix identities.
4. **What does not change.** Everything at $\mathbf q\to0$ and $\mathbf G=0$:
   - the Fröhlich limit of E.131;
   - the Huang equations of §4.2.1 (E.144–165);
   - E.406–E.410 at $\mathbf q=0$;
   - the LO–TO splitting (the $q\to0$ limit of E.454);
   - the polariton region (relative change $5\times10^{-9}$ in GaAs).
5. **What changes.**
   - Finite $q$ in multi-atom cells: in GaAs, $|g_{\rm LO}|^2$ changes by up to 1.5 % at 0.1 Å⁻¹ and by 4–33 % at 0.5 Å⁻¹.
   - All $\mathbf G\neq0$ (local-field) components.
   - Acoustic modes: the printed vertex predicts a spurious, cell-dependent "piezoelectric" coupling with $|g|^2\propto1/q$, even in centrosymmetric crystals. For GaAs it amounts to an apparent piezoelectric constant of 0.6–1.8 C/m², against the physical 0.19 C/m² for LA[111] and exactly 0 for LA[100].
6. **Correction to the task statement.** The task wrote that "$\sum_aZ^*_ae^{-i\mathbf q\mathbf r_a}\approx-i\sum_aZ^*_a(\mathbf q\mathbf r_a)$ gives an $O(q)$ dipolar coupling" for acoustic modes. In EPP's (and Giustino's) cell-phase convention the acoustic eigenvector itself carries $e^{+i\mathbf q\mathbf r^0_a}$.
   - The *printed* (phase-less) vertex therefore contains the spurious term $+\frac{i}{\sqrt M}\sum_a(\mathbf q\mathbf r^0_a)\hat Z^*_a\mathbf u$.
   - The *corrected* vertex cancels that term exactly, and what remains at $O(q)$ is physical (internal strain).
   - The expression $-i\sum_aZ^*_a(\mathbf q\mathbf r_a)$ only arises if one multiplies the explicit phase by an atom-phase (real) eigenvector, i.e. if conventions are mixed.
   - The magnitude estimate is right; the attribution is reversed. Verified numerically (§6, checks C3).

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
| Volumes | $\Omega$ = crystal, $\Omega_{\rm uc}=\Omega/N$; Giustino's $\Omega$ = EPP $\Omega_{\rm uc}$ | p. 20; G p. 7 |

**Eigenvector remark (already Z.2 errata item 3).**
- E.15/E.35 define $\tilde\Omega^2(\mathbf k)$ with $e^{+i\mathbf k\boldsymbol\delta}$, $\boldsymbol\delta=\mathbf R^0_i-\mathbf R^0_j$. That is Giustino's $D^{\rm dm}(-\mathbf k)=D^{\rm dm}(\mathbf k)^*$.
- The step after E.43 (p. 9) uses $\tilde\Omega^2_{b\nu a\mu}(\mathbf k)=\tilde\Omega^2_{a\mu b\nu}(\mathbf k)$, which is true only for real $\tilde\Omega^2$.
- E.83 is consistent only if $e_{a\mu\alpha}(\mathbf k)$ are the eigenvectors of $(m_am_b)^{-1/2}\sum_{\boldsymbol\delta}e^{-i\mathbf k\boldsymbol\delta}k(\boldsymbol\delta)$, i.e. those that solve the equations of motion for the pattern $e^{i\mathbf k\mathbf R^0_i}$. I checked this numerically (section B of the script): residual $2\times10^{-16}$ for those eigenvectors, against 0.61 for eigenvectors of the printed E.35.
- The basis-phase correction needs only the cell phase of E.83. Its sign is fixed by E.83 together with E.117, not by E.15. I assume the item-3 fix throughout and say so in a footnote of the `.tex`.

---

## 3. Derivation

### 3.1 Why the basis phase must appear explicitly

Mode $(\mathbf q,\alpha)$ displaces atom $a$ of cell $i$ by $\propto e_{a\nu\alpha}(\mathbf q)e^{i\mathbf q\mathbf R^0_i}$ (E.83). The dipole $p_{ia\mu}=Z^*_{a\mu\nu}Q_{ia\nu}$ sits at $\mathbf R^0_i+\mathbf r^0_a$, so the polarization density is $\mathbf P(\mathbf r)=\sum_{ia}\mathbf p_{ia}\delta(\mathbf r-\mathbf R^0_{ia})$. With the transform of E.117 its Fourier component is

$$\mathbf P(\mathbf k)=\sum_{ia}\mathbf p_{ia}e^{-i\mathbf k(\mathbf R^0_i+\mathbf r^0_a)}\propto\sum_a\hat Z^*_a\mathbf e_{a\alpha}(\mathbf q)e^{-i\mathbf k\mathbf r^0_a}\sum_ie^{i(\mathbf q-\mathbf k)\mathbf R^0_i}=N\sum_{\mathbf G}\delta_{\mathbf k,\mathbf q+\mathbf G}\sum_a\hat Z^*_a\mathbf e_{a\alpha}(\mathbf q)e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}.$$

Hence $e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}$ multiplies $e_{a\nu\alpha}(\mathbf q)$ in every vertex (checked, script section A).

**Cell-choice invariance (131a).** Suppose atom $a$ is assigned to a neighbouring cell, $\mathbf r^0_a\to\mathbf r^0_a+\mathbf R_L$. The physical displacement of each atom is unchanged, so the cell-phase eigenvector component becomes $e_{a\nu\alpha}\to e^{i\mathbf q\mathbf R_L}e_{a\nu\alpha}$. Only the product $e_{a\nu\alpha}(\mathbf q)e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}$ is invariant, so the printed vertex depends on an arbitrary bookkeeping choice. Script check C1, GaAs model with As at $(1,1,1)a/4$ or at $(1,-1,-1)a/4$:
- corrected vertices: identical to $2\times10^{-13}$;
- printed vertices: differ by up to $0.15\,e/\sqrt{\rm amu}$ (the LO vertex itself is $0.34$).

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

### 3.4 Acoustic modes, and the correct reading of footnote 7

In the cell-phase convention a long acoustic wave displaces every atom as $\mathbf u\,e^{i\mathbf q(\mathbf R^0_i+\mathbf r^0_a)}$, plus $O(q)$ internal relaxation. So

$$e_{a\nu\alpha}(\mathbf q)=e^{i\mathbf q\mathbf r^0_a}\left[\sqrt{m_a/M}\,u_\nu+O(q)\right],\qquad M=\sum_am_a.$$

The script checks this: for a symmetric chain, $\arg(e_2/e_1)=qa/2$ to $10^{-9}$. Then (131b) follows:
- printed: $\sum_a\hat Z^*_a\mathbf e_a/\sqrt{m_a}=\frac{i}{\sqrt M}\sum_a(\mathbf q\mathbf r^0_a)\hat Z^*_a\mathbf u+O(q)_{\rm phys}$;
- corrected: $\sum_a\hat Z^*_a\mathbf e_ae^{-i\mathbf q\mathbf r^0_a}/\sqrt{m_a}=\frac{1}{\sqrt M}(\sum_a\hat Z^*_a)\mathbf u+O(q)_{\rm phys}=O(q)_{\rm phys}$.

**Properties of the spurious term.** It is origin-independent, because $\sum_a\hat Z^*_a=0$. It does depend on the assignment of atoms to cells. With $\omega_{\rm ac}\propto q$ it gives $g^{\mathcal L}\propto\frac1q\cdot q\cdot q^{-1/2}$, i.e. $|g|^2\propto1/q$: the scaling of a piezoelectric coupling, with an apparent piezoelectric constant $e_{\rm fake}=|\sum_aZ^*_a\hat{\mathbf q}\mathbf r^0_a|/\Omega_{\rm uc}$. It is nonzero even in centrosymmetric crystals. For rocksalt with Cl at $(\pm\frac12,0,0)a$, the sign of $e_{\rm fake}$ flips with the cell choice.

**Script checks C3.**
- LA[100] in zincblende is piezo-inactive. The corrected vertex scales as $q^{2.00}$; the printed one has slope $Z^*(a/4)/\sqrt M=0.2432\,e$Å$/\sqrt{\rm amu}$ for both cells.
- LA[111]: the corrected slope is $0.2809$ in the toy model for both cells (physical, from internal strain). The printed slope is $0.1404$ or $0.4213$, depending on the cell.
- At $q=10^{-4}$ Å⁻¹, printed minus corrected equals $i\sum_aZ^*_a(\mathbf q\mathbf r_a)/\sqrt M$ to $10^{-3}$ relative.

**The physical $O(q)$ term.** It comes from internal strain plus the clamped-ion (electronic) piezoelectric response. The latter requires the dynamical quadrupoles, which enter the long-range potential at the same order as the internal-strain dipole term. Literature (flagged, citations from memory; authors should verify):
- G. Brunin *et al.*, Phys. Rev. Lett. 125, 136601 (2020);
- V. A. Jhalani *et al.*, Phys. Rev. Lett. 125, 136602 (2020);
- J. Park *et al.*, Phys. Rev. B 102, 125203 (2020);
- background: M. Royo and M. Stengel, Phys. Rev. X 9, 021050 (2019).

So the corrected E.131′ has the right structure at $O(q^0)$ and the right *dipole* part at $O(q)$, but it is not complete at $O(q)$. This is what Giustino's footnote 7 means by "needs to be retained in order to correctly describe the acoustic modes near $\mathbf q=0$".

### 3.5 §4 (Gaussian)

The coupling is $\frac1c\sum_{ia}\dot{\mathbf p}_{ia}\mathbf A(\mathbf R^0_{ia})$ (134′). The Euler–Lagrange equation E.139 then gives $\partial L/\partial\dot Q_{ia\mu}=m_a\dot Q_{ia\mu}+\frac1c\sum_\nu Z^*_{a\nu\mu}A_\nu(\mathbf R^0_{ia})$. Its time derivative involves $\partial_tA$ at the fixed equilibrium site (dipole approximation within each atom), so (141′)–(143′) have $A_\nu$ and $E_\nu$ at $\mathbf R^0_{ia}$. In §4.2.1 the field is taken uniform over the cell ($|\mathbf q|\ll1/a$), so E.144–E.165 are unchanged.

Further changes:
- E.167 becomes (167′) = (241′) = (314′).
- E.183, E.186, E.237 and E.270 get $A(\mathbf R^0_{ia})$.
- E.175–E.180 are unchanged in form; $\mathbf P$ is now (167′).
- E.235 is unchanged in form, with $\rho_{\rm ion}=-\nabla\mathbf P$ built from (241′).

### 3.6 §5: ionic photon self-energy, $\mathbf G$-resolved

E.272′ is $j^{\rm ion}_\alpha=\sum_{ia\mu}\frac{1}{m_a}Z^*_{a\alpha\mu}\pi_{ia\mu}\delta(\mathbf r-\mathbf R^0_{ia})$, with $\pi_{ia\mu}$ from E.243 ($\propto e^{i\mathbf q\mathbf R^0_i}$). The "continuous field" of E.274 must read $e^{i\mathbf q(\mathbf r-\mathbf r^0_a)}$ (274′), so that it equals E.243 at $\mathbf r=\mathbf R^0_{ia}$. The contraction of E.287–E.289 is unchanged. Taking $n_{\mathbf q}=0$ gives

$$\Pi^{\rm ion,para}_{\alpha\beta}(\mathbf r,\mathbf r',\omega)=-\frac{1}{\hbar cN}\sum_{\sigma\mathbf q}\sum_{ia\mu,jb\nu}\frac{Z^*_{a\alpha\mu}Z^*_{b\beta\nu}}{\sqrt{m_am_b}}\frac{\hbar\omega^2_\sigma}{\omega^2-\omega^2_\sigma+i\eta}e_{a\mu\sigma}e^*_{b\nu\sigma}e^{i\mathbf q(\mathbf R^0_i-\mathbf R^0_j)}\delta(\mathbf r-\mathbf R^0_{ia})\delta(\mathbf r'-\mathbf R^0_{jb})\quad(290').$$

With $\delta(\mathbf r-\mathbf R^0_{ia})=\int\frac{d\mathbf q_1}{(2\pi)^3}e^{i\mathbf q_1(\mathbf r-\mathbf R^0_i-\mathbf r^0_a)}$, the sum over $i$ gives $\frac N\Omega\sum_{\mathbf G}e^{i(\mathbf q+\mathbf G)(\mathbf r-\mathbf r^0_a)}$, and similarly for $j$ with $\mathbf G'$. The prefactor $\frac{\Omega}{N}\frac{N^2}{\Omega^2}\frac{1}{c}=\frac{1}{c\Omega_{\rm uc}}$ is as in E.293. Comparing with E.294 yields (295′)–(297′):

$$j_{\alpha\sigma}(\mathbf q+\mathbf G)=\sum_{a\mu}\frac{Z^*_{a\alpha\mu}}{\sqrt{m_a}}e_{a\mu\sigma}(\mathbf q)e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a},\qquad\Pi^{\rm ion,para}_{\alpha\beta,\mathbf G\mathbf G'}=-\frac{1}{c\Omega_{\rm uc}}\sum_\sigma\frac{\omega^2_\sigma}{\omega^2-\omega^2_\sigma+i\eta}j_{\alpha\sigma}(\mathbf q+\mathbf G)j^*_{\beta\sigma}(\mathbf q+\mathbf G').$$

In the cell-phase convention $e_a(\mathbf q+\mathbf G)=e_a(\mathbf q)$, because $\tilde\Omega^2$ is $\mathbf G$-periodic, so $j_{\alpha\sigma}(\mathbf q+\mathbf G)$ is a well-defined function of the full wave vector $\mathbf q+\mathbf G$. It is not periodic: the basis phase supplies the structure factor of each sublattice.

**Diamagnetic part.** Applying the same lattice sum to E.302′ gives $\frac{1}{\Omega_{\rm uc}}\sum_{\mathbf G}e^{i\mathbf G(\mathbf r-\mathbf r^0_a)}\int\frac{d\mathbf q'}{(2\pi)^3}e^{i\mathbf q'(\mathbf r-\mathbf r')}$. Decompose $\mathbf q'=\mathbf q+\mathbf G'$ and relabel $\mathbf G+\mathbf G'\to\mathbf G$ (the implicit relabelling on p. 66). This gives (304′), $\Pi^{\rm ion,dia}_{\alpha\beta,\mathbf G\mathbf G'}=-\frac{1}{\Omega_{\rm uc}c}\sum_{a\mu}\frac{Z^*_{a\alpha\mu}Z^*_{a\beta\mu}}{m_a}e^{-i(\mathbf G-\mathbf G')\mathbf r^0_a}$.

**Completeness.** In (309′), $\sum_\sigma e_{a\mu\sigma}e^*_{b\nu\sigma}=\delta_{ab}\delta_{\mu\nu}$ (E.28) keeps only $a=b$, so $e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}e^{+i(\mathbf q+\mathbf G')\mathbf r^0_a}=e^{-i(\mathbf G-\mathbf G')\mathbf r^0_a}$. Hence E.310 and E.311 hold as $\mathbf G\mathbf G'$ matrix identities (script check C4: $|\Pi^{\rm dia}+\Pi^{\rm para}(0)|<10^{-16}$ for 16 $(\mathbf G,\mathbf G')$ pairs). The diagonal $\mathbf G=\mathbf G'$ of (304′) is E.304 as printed.

**2D.** For a planar sheet, (300′)/(305′) are the same with in-plane $\mathbf G$ and in-plane positions $\boldsymbol\rho^0_a$. For atoms at heights $z_a$, the factor $\delta(z)\delta(z')$ becomes $\sum_{ab}j^a_{\alpha\sigma}j^{b*}_{\beta\sigma}\delta(z-z_a)\delta(z'-z_b)$, with the atom-resolved vertex $j^a$ of (295′).

### 3.7 §6

E.352 becomes the matrix (352′). The derivation on pp. 76–77, with $\delta(\mathbf r-\mathbf R^0_{ia})\delta(\mathbf r'-\mathbf R^0_{jb})$, yields (354a): $\Pi^{\rm ion}_{\rho\rho,\mathbf G\mathbf G'}=-\frac{c}{\omega^2}(q+G)_\mu(q+G')_\nu\Pi^{\rm ion}_{jj\mu\nu,\mathbf G\mathbf G'}$, using $\partial_\mu\to i(q+G)_\mu$ and $\partial'_\nu\to-i(q+G')_\nu$. So the connection formula E.351 still holds, and E.341–E.350 are general and unchanged.

### 3.8 §§6.1–8

From E.380 on, the manuscript drops $\mathbf G\mathbf G'$. That step was justified for the ions by "the ionic self-energy is a constant in $\mathbf G\mathbf G'$ space" (p. 63), which is no longer true. For the electrons it was justified by "CLFs are mostly neglected here" (p. 95). Within the $\mathbf G=\mathbf G'=0$ (macroscopic) approximation, every equation from E.385 to E.472 keeps its printed form, provided the vertex is (385′), $j_{\mu\sigma}(\mathbf q)=\sum Z^*e\,e^{-i\mathbf q\mathbf r^0_a}/\sqrt m$. With local fields kept, the phonon self-energy is the matrix (399′).

Electrostatic limits:
- **p. 95.** $D^0\to-\frac{4\pi c}{\omega^2}\frac{q_\alpha q_\beta}{q^2}$ (E.442), combined with E.428, gives (428b). This equals the $\mathbf G=0$, $\hat\epsilon^\infty=1$ term of G Eq. 180, including the phase.
- **p. 98.** Equation (448b) carries $e^*_{a\nu\sigma}e^{+i\mathbf q\mathbf r^0_a}$, the complex conjugate, consistent with the author's remark that it "represents the conjugated complex part".
- **p. 96.** The display that EPP calls Verdi–Giustino Eq. 4 should cite G Eq. 180 and include $e^{-i\mathbf q\mathbf r^0_a}$.
- **Interplay with `corrections/issue-efactor.tex`.** That correction fixes the missing factor $e$ in the same displays and uses the tags (p95′), (p96′), (428a), (p98′) and (448a). The two corrections are independent: the factor $e$ changes $4\pi i/\Omega_{\rm uc}\to4\pi ie/\Omega_{\rm uc}$, and the phase multiplies the atom summand. To avoid clashing tags, this correction uses (428b) and (448b). If both are adopted, (p96′) there must drop "with $\boldsymbol\tau_\kappa=0$ as for our dipoles at $\mathbf R^0_i$ (314)" and keep the phase $e^{-i\mathbf q\mathbf r^0_a}$. Its relations (428a) and (448a) between our vertex and Verdi–Giustino's then hold with the same phase on both sides.

**E.454 and Gonze–Lee.** E.454 keeps its form with (385′). Its $\mathbf G$-resolved matrix form (454a) is the dipole–dipole dynamical matrix of Gonze & Lee (1997), with the structure factor $e^{i(\mathbf q+\mathbf G)(\mathbf r^0_a-\mathbf r^0_b)}$. Its projection on mode $\sigma$ at $\mathbf G=0$ is $e^*_\sigma\Delta\tilde\Omega^2e_\sigma=\frac{W}{\Omega_{\rm uc}}|\mathbf q\,j_\sigma|^2$, i.e. E.454. The sign of the structure factor was verified independently by an Ewald test (§6.2).

---

## 4. Equation-by-equation list

**Corrected (primed in the `.tex`):**
- §3: 121′, 122′, 123′, 124′, 131′.
- §4: 134′, 141′–143′ (shown: 143′), 167′, 183′, 186′, 237′ (A at $\mathbf R^0_{ia}$).
- §5: 241′, 270′, 271′, 272′, 274′, 286′–290′ (shown: rule), 293′, 295′, 296′, 297′, 299′, 300′, 302′, 304′, 305′, 309′.
- §6: 314′, 315′, 352′, 385′, 399′.

**New (lettered):**
- 131a: cell-choice invariance.
- 131b: acoustic expansion.
- 354a: $\Pi_{\rho\rho,\mathbf G\mathbf G'}$.
- 428b: replaces the unnumbered display on p. 95.
- 448b: replaces the boxed unnumbered display on p. 98.
- 454a: $\mathbf G$-resolved E.454 = Gonze–Lee.

**Same printed form, new content through the corrected vertex or density:**
- 125–130 (with 124′);
- 169, 175–180 and 235 (with 167′/241′);
- 306–308, 310, 311 (as $\mathbf G\mathbf G'$ matrices);
- 313, 351;
- 386–417, 429–430, 447–453, 454;
- 455–472 (2D, in-plane phase).

**Unchanged:**
- 1–120, including E.120, the potential of a single dipole;
- 144–165 (uniform field over a cell);
- 166, 168, 170–174, 181–182, 184–185, 187–234, 236 (photon field), 238–240, 242–269, 275–285 (electronic and photon parts), 316–350;
- 406–410 at $\mathbf q=0$ (E.409–410 explicitly);
- the $q\to0$ limits of 131′ (Fröhlich) and 454 (LO–TO).

---

## 5. Checks

### 5.1 Dimensions (Gaussian)
- The correction factor $e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}$ is dimensionless, so no printed dimension changes.
- E.131′ in Gaussian units: $[eZ^*]=$ statC$^2=$ erg·cm. Then $\Omega_{\rm uc}^{-1}$ (cm$^{-3}$) × $\sqrt{\hbar/m\omega}$ (cm) × $q/q^2$ (cm) gives erg.
- $[j]=$ statC/$\sqrt{\rm g}$, as printed.
- Script check F: the SI evaluation ($Z^*_{\rm G}$ dimensionless, $e^2/4\pi\epsilon_0$) and the Gaussian evaluation ($e^2=14.3996$ eV Å, $Z^*_{\rm EPP}=|e|Z^*_{\rm G}$) of $\sqrt N|g^{\mathcal L}|$ at $q=0.1$ Å⁻¹ agree: 303.991 vs 303.990 meV.

### 5.2 Limits
- **$q\to0$, $\mathbf G=0$.** The phase goes to 1. The optical-triplet strength $\sum_{\rm opt}|\hat{\mathbf q}j|^2=Z^{*2}/\mu=0.11865$, identical for both forms (check C2). E.131′ equals Fröhlich (G Eq. 7) exactly when $Z^*$ is the LO–TO-implied 2.177 (check E). With the input $Z^*=2.07$ the two differ by 5 %: 0.0304 vs 0.0320 eV Å⁻¹, an inconsistency of the inputs, not of the formula. LO–TO: $\hbar\omega_p=44.56$ meV, $\hbar\omega_{\rm LO}=35.93$ meV (measured 36.2).
- **$c\to\infty$.** (428b) is G Eq. 180 at $\mathbf G=0$, $\epsilon^\infty=1$, including the phase. (448b) is its conjugate with $\hat\epsilon^{\rm e}$.
- **Polariton (Huang).** $q_c=\omega_{\rm TO}\sqrt{\epsilon_0}/c=6.06\times10^{-5}$ Å⁻¹ and $q_cd=1.5\times10^{-4}$, so $|j_{\rm LO}|^2$ changes by $5.5\times10^{-9}$. E.157–E.165 and E.410 are unaffected.

### 5.3 Gonze–Lee: Ewald test of the $\mathbf G$-resolved structure (script section D)
- **Setup.** Energy per cell of a modulated dipole lattice $\mathbf p_{ia}=\mathbf p_ae^{i\mathbf q\mathbf R^0_i}$ at $\mathbf R^0_{ia}$ (GaAs fcc geometry, random complex $\mathbf p_a$). The reciprocal sum uses the structure factor of (295′); the real-space part uses the erfc-screened dipole tensor; the Gaussian self term is subtracted.
- **Correct sign $e^{-i(\mathbf q+\mathbf G)\mathbf r^0_a}$.** $E=-0.070411\,e^2$/Å, independent of the Ewald parameter ($\Lambda=0.45,0.6,0.8$ Å⁻¹) to $3\times10^{-14}$.
- **Wrong sign.** The result varies by 34 %.
- **Dipoles at the cell origin.** The intracell $a\neq b$ pair coincides and its interaction diverges, so the $\mathbf G$ sum of printed E.297 has no limit.

### 5.4 Script summary
`python3 test/corr-phase/check_phase.py` runs 22 checks, all PASS, in sections:
- A: Fourier sign;
- B: E.15 remark;
- C: zincblende model, covering cell invariance, $q\to0$, acoustic scaling and E.310′;
- D: Ewald;
- E: unchanged limits;
- F: units.

---

## 6. Numbers (GaAs-like; $d=|\mathbf r^0_{\rm As}-\mathbf r^0_{\rm Ga}|=2.448$ Å)

Inputs: $\hbar\omega_{\rm TO}=33.3$ meV, $\hbar\omega_{\rm LO}=36.2$ meV, $\epsilon^\infty=10.9$, $Z^*=2.07\,e$, $\Omega_{\rm uc}=45.2$ Å³, masses 69.72/74.92 amu ($\mu=36.1$), $a=5.653$ Å, $v_{\rm LA}=5.40$ km/s ([111]) and 4.73 km/s ([100]).

**(a) Two-atom estimate with leading-order eigenvectors.** $\theta=\mathbf q\,(\mathbf r^0_{\rm As}-\mathbf r^0_{\rm Ga})$.

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

**(b) 3D toy model.**
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

How to read the tables:
- The corrected columns do not depend on the cell; the printed ones do.
- Along [111] the $O(q)$ internal-strain phase of the optical eigenvector can push the printed LO value either up or down. That is why the model gives $|g|^2$ ratios of 1.005 and 0.992 at 0.1 Å⁻¹, where the leading-order estimate gave 0.985.
- Across (a) and (b), the LO error is 0.2–1.5 % at 0.1 Å⁻¹ and 4–33 % at 0.5 Å⁻¹.

**(c) Fake piezoelectric constant of the printed vertex**, $e_{\rm fake}=|\sum_aZ^*_a\hat{\mathbf q}\mathbf r^0_a|/\Omega_{\rm uc}$:
- 1.80 C/m² for [111] with As at (1,1,1)a/4;
- 0.60 C/m² for [111] with As at (1,−1,−1)a/4;
- 1.04 C/m² for [100].

The physical values for GaAs (literature $|e_{14}|=0.16$ C/m², illustrative) are $(2/\sqrt3)|e_{14}|=0.185$ C/m² for LA[111] and 0 for LA[100].

---

## 7. What the primer should change

| Section | Change |
|---|---|
| **B01-foundations** (§B.1) | Replace the "Note on the source: basis phase" callout with a *correction* callout. Display (E.124′) and (E.131′); tag them as corrected and cite `corrections/issue-phase.tex`, since the guide forbids inventing source numbers, so mark them "corrected". Add the one-line cell-reassignment argument (131a). Add an interactive "cell-choice test": a GaAs-like chain (or zincblende along [111]/[100]) with toggles "printed E.131 / corrected E.131′" and "As at +d / at d−a". Plot $\lvert g_{\rm LO}(q)\rvert^2$ and $\lvert g_{\rm LA}(q)\rvert^2$. The corrected curves must not move when the cell toggles. The printed LA curve diverges as $q^{-1/2}$, while the corrected one vanishes along [100]. |
| **A04-polar** (§A.4) | Rewrite the acoustic bullet under (G.180) precisely. In Giustino's (Maradudin–Vosko, cell-phase) convention the acoustic eigenvector carries $e^{i\mathbf q\boldsymbol\tau_\kappa}$; the explicit $e^{-i(\mathbf q+\mathbf G)\boldsymbol\tau_\kappa}$ cancels it. Without it, $\lvert g\rvert^2\propto1/q$ is spurious, cell-dependent, and present even in centrosymmetric crystals. The physical $O(q)$ coupling also needs dynamical quadrupoles (Brunin/Jhalani/Park 2020, literature). Update the source callout to point to the correction. If the 1D lab is meant to "reproduce this logic", include the basis phase there too. |
| **B02-lagrangian** (§B.2) | In the "Dipole approximation" bullet, state the corrected coupling (E.134′) with $\mathbf A(\mathbf R^0_{ia})$. Keep "harmless at photon wave vectors", quantified as $q_cd\approx10^{-4}$ and a relative change of $5\times10^{-9}$ for GaAs. |
| **B03-dielectric** (§4.2.1) | Note that (E.143′) puts $\mathbf E$ at $\mathbf R^0_{ia}$ and that the long-wavelength results E.146–E.165 are unchanged. |
| **B04-fields** | (E.167′)/(E.241′): dipoles at atomic sites. $\rho_{\rm ion}=-\nabla\mathbf P$ acquires the structure factors at $\mathbf G\neq0$. Remove "drop the basis phase" from line 161. |
| **B05-photon-se** | The largest change. Display (E.295′), (E.297′) and (E.299′) as $\mathbf G\mathbf G'$ matrices with $j(\mathbf q+\mathbf G)$. In the derivation, the lattice sum $\sum_i\delta(\mathbf r-\mathbf R^0_{ia})=\Omega_{\rm uc}^{-1}\sum_{\mathbf G}e^{i\mathbf G(\mathbf r-\mathbf r^0_a)}$ now gives structure factors. Replace "The $\mathbf G$-independence is an artefact … only $\mathbf G=\mathbf G'=0$ should be used" with the replacement text of the `.tex` (sublattice dipole gratings, reliable only at small $\lvert\mathbf q+\mathbf G\rvert$). Display (E.304′). In the completeness subsection, show that the phases cancel for $a=b$, so E.309–E.311 hold as matrices. In **viz-b-photon-se-sumrule**, give the three ions explicit positions and use (295′). Add a toggle "printed (E.295) / corrected (E.295′)": the stacked sum (E.309) is the same for both, but the per-branch weights change, and in the printed form the acoustic branch acquires spurious weight $\propto q^2$. |
| **B06-propagator** | P in line 18 becomes (E.314′). Add (E.354a), $\Pi^{\rm ion}_{\rho\rho,\mathbf G\mathbf G'}$. |
| **B07-phonon-se** | In the comparison table, change "point dipoles at cell origins" to "point dipoles at atomic sites $\mathbf R^0_{ia}$ (E.314′)". State that the vertex is (E.385′) and that E.409–E.410 ($\mathbf q=0$) are unchanged. |
| **B08-frohlich** | Replace "The basis phase of G Eq. 180 is missing" with the corrected limits (428b)/(448b), which equal the $\mathbf G=0$ term of G.180 including $e^{\mp i\mathbf q\mathbf r^0_a}$. Any vertex plot at finite $q$ should use the phase. |
| **B09-electrostatic** | With (E.385′), E.454 equals the Gonze–Lee $\mathbf G=0$ term at all $q$, not only "at small $q$". Display (E.454a) as the $\mathbf G$-resolved form (the Gonze–Lee dipole–dipole matrix after Ewald regularization) and mention the Ewald $\Lambda$-independence check. Give the finite-$q$ size in GaAs (up to 1.5 % at 0.1 Å⁻¹, 4–33 % at 0.5 Å⁻¹). Change "Born-charge point dipoles at the cell origins" in the bridge callout to "at the atomic sites". |
| **B10-2d** | Add the in-plane basis phase $e^{-i\mathbf Q\boldsymbol\rho^0_a}$ in (300′)/(305′); for h-BN, B and N sit at inequivalent in-plane sites. For buckled or multilayer sheets (MoS₂ with S at $\pm z$), add the $\delta(z-z_a)$ replacement. |
| **X02-dictionary** | Add rows: basis vector EPP $\mathbf r^0_a$ ↔ G $\boldsymbol\tau_\kappa$; phase convention cell phase in both (E.7/E.83 ↔ G.14/G.20; mind Z.2 item 3); ionic vertex (E.295′) ↔ the Verdi–Giustino numerator with structure factor. |
| **Z02-errata** | Item 14: list the corrected equations (from E.121′ to E.454a, as in §4 above) and state the effect precisely (spurious acoustic $\lvert g\rvert^2\propto1/q$, fake piezoelectric constant 0.6–1.8 C/m² in GaAs; finite-$q$ LO error; $\mathbf G\neq0$). Link `corrections/issue-phase.tex`. Update the intro sentence (line 9) accordingly. Cross-reference item 29 (point-dipole self-energy), which (E.454a) needs for its regularization. |

---

## 8. Open questions and caveats

1. **Task framing on the acoustic sign.** See §1, item 6, and §3.4. The spurious $O(q)$ term belongs to the printed vertex; the corrected vertex removes it. The task's expression $-i\sum_aZ^*_a(\mathbf q\mathbf r_a)$ mixes the cell-phase factor with an atom-phase eigenvector. The magnitude estimate is right.
2. **E.15 sign (Z.2 item 3).** The correction assumes that $e(\mathbf q)$ in E.83 solves the $e^{i\mathbf q\mathbf R^0_i}$ equations of motion (Giustino's G.14). If the authors kept E.15 and E.26 literally, the eigenvectors would come out conjugated, and the phase would pair with the wrong vectors. Item 3 should be fixed together with this one.
3. **Dynamical quadrupoles.** The literature citations (Brunin, Jhalani, Park 2020; Royo–Stengel 2019) are from memory and flagged; the authors should verify the bibliographic details. No quadrupole formula is given in the `.tex`.
4. **Toy-model numbers.** The 3D model (central NN plus NNN springs) is a toy: it fixes the internal strain, hence the physical LA[111] slope and the [111] LO numbers. The quoted ranges combine it with the model-free two-atom estimate. The GaAs $e_{14}$ is a literature value.
5. **Regularizing the $\mathbf G$ sum.** In (454a) and in any $\mathbf G\neq0$ use of (297′) and (304′), the on-site ($a=b$) point-dipole terms diverge. They need a consistent prescription: Ewald self-term removal, the cutoff $e^{-a|\mathbf q+\mathbf G|^2}$, and avoiding double counting with the short-range force constants $k^{\mathcal S}$ and with the self-energies kept in E.235 (Z.2 item 29). This is beyond the present correction.
6. **The approximation in §§6.1–8 should be stated.** Dropping $\mathbf G\mathbf G'$ there is now a genuine approximation for the ions as well as for the electrons: it neglects ionic local fields. The manuscript should say so where it drops $\mathbf G\mathbf G'$ (p. 64 and p. 95).
7. **Overlap with the e-factor correction.** `issue-efactor.tex` sets $\boldsymbol\tau_\kappa=0$ in its (p96′), "as for our dipoles at $\mathbf R^0_i$ (314)". This must be reconciled: with (314′) the phase stays. Tags (428b) and (448b) were chosen to avoid its (428a) and (448a).
8. **2D buckled sheets.** The $\delta(z-z_a)$ generalization is only sketched. For MoS₂-type layers the out-of-plane positions enter $D^{\rm e}(\mathbf Q,\omega;z_a,z_b)$ inside the atom sums of E.458.
