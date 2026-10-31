# Optics beyond RPA: excitons, local fields, xc and phonon-assisted absorption

Labels: **[est]** established literature; **[derived]** worked out here from EPP's equations, not printed in EPP; **[proposal]** suggested extension; **[spec]** speculation. Units are Gaussian, as in EPP. "G" is Giustino's review, "E" the EPP draft.

## 1. What is missing in EPP

- **Independent-particle bubble.** E.278–E.279 (pp. 57–58) define the photon self-energy as the lowest-order bubble plus contact (diamagnetic) term per species (primer [§B.5](../../sections/B05-photon-se.html)). E.282 and E.383 (pp. 59, 85) evaluate it with single-particle energies $E_{n\bk}$ and current vertices (E.284, E.384). No electron–hole ladder means no excitons; no self-energy insertions means no $GW$ gap and no phonon renormalization.
- **No local fields, no xc.** From E.380 on (p. 85) EPP drops the $\bG$ indices, so only $\bG=\bG'=0$ is kept for electrons. For ions, printed E.293–E.295 has no basis phase at all (after E.293 the ionic self-energy is called "a constant in $GG'$ space"); the primer's corrected vertex E.295′ carries the structure factors ([§Z.3](../../sections/Z03-corrections.html)). Hence $W^{\rm e}(\bq,\omega)=4\pi/(\bq\cdot\hat\epsilon^{\rm e}\cdot\bq)$ (E.445), with $\hat\epsilon^{\rm e}$ (E.432) the head of the RPA dielectric matrix without local fields. The $\epsilon_\infty$ implicit in $W^{\rm e}(\bq,0)$ of E.454, and the $r_0$ implicit in $\epsilon^{\rm e}_{\rm 2D}(Q,0)$ of E.472, are RPA numbers without local fields. No $f^{xc}$ appears anywhere. DFPT screening (G.54) is static but includes both ([§A.3](../../sections/A03-dfpt.html), [§B.6](../../sections/B06-propagator.html)).
- **No non-RPA electronic insertions.** Mixed $\langle j^{\rm e}j^{\rm ion}\rangle$ bubbles vanish at lowest order (E.316). Phonon-assisted absorption (G.200) and exciton–phonon coupling are instead the first non-RPA electronic insertions: self-energy and vertex diagrams of $\Pi^{\rm e}$ containing an internal $D^{\rm ph}$ line. EPP has no direct electron–ion vertex, so the phonon reaches the electrons only through this photon line. Because these diagrams are photon-irreducible, they live inside $\Pi^{\rm e}$ and are not double counted by the Dyson resummation of E.379. They are not derived ([§A.9](../../sections/A09-applications.html), [§Z.1](../../sections/Z01-outlook.html)).
- **No short-range electron–phonon coupling.** The Lagrangian freezes the crystal potential (§Z.1), so $D^{\rm ph}$ carries only the dipolar (Fröhlich) part of $g$. Any phonon-assisted process built from EPP alone is its polar part only ([shortrange.md](shortrange.md)).
- **$T=0$ boson lines.** The $n_B=0$ rule (E.249, §B.5) is safe for E.379. It fails once a phonon line sits inside $\Pi^{\rm e}$, which needs the weights $n_{\bq\nu}+\frac12\pm\frac12$ of G.200.
- **2D.** The sheet bubble is $l\,\Pi^{\rm e}_{\bG=\bG'=0}$ (p. 60). EPP writes a $z$-resolved $D^{\rm e}(Q,\omega;z,z')$ "including electron $z$-LFE" (p. 100), but keeps only $\Pi^{{\rm e},0}_{\rho\rho}(Q,\omega)\delta(z)\delta(z')$ (before E.461). It therefore has Rytova–Keldysh screening (E.469–E.470, without $\epsilon_{\rm env}$) but no electron–hole ladder, so neither intralayer nor interlayer excitons.

## 2. Why it matters

- **Excitons dominate 2D optics [est].** Monolayer WS₂ has an exciton binding energy of 0.32 eV and a non-hydrogenic Rydberg series [5]. An RPA $\Pi^{\rm e}$ misses the peaks that carry the oscillator strength.
- **Linewidths are meV [est].** hBN-encapsulated MoS₂ shows exciton linewidths of about 2 meV FWHM at 4 K, approaching the homogeneous limit [6]. What remains is radiative (transverse photon) plus exciton–phonon broadening. EPP covers the radiative channel and the polar (Fröhlich) part of the phonon channel. The acoustic and non-polar optical contributions need the short-range coupling that EPP lacks ([shortrange.md](shortrange.md)). In a MoS₂-like sheet the Fröhlich channel is weak ($r_0^{\rm ion}\approx0.5$ Å with the primer preset, §4.3) [derived].
- **Screening inputs propagate [est].** Giustino (§XII; primer §A.9) estimates that a 20% error in $\epsilon_\infty$ means about 40% in band renormalization and 10% in phonon frequencies. EPP's $\epsilon_\infty$ and $r_0$ currently lack local fields and xc.
- **Temperature [est].** Zero-point motion shifts and broadens excitons in Si and hBN [7]. The $GW$ zero-point renormalization of diamond's direct gap is −628 meV (Antonius *et al.* [16], via G Fig. 13), which moves every absorption edge.
- **Phonon screening of $W$ [est].** Adding phonon screening to the BSE kernel lowers exciton binding energies in halide perovskites [12]. In EPP this screening is already inside $D$ (§4.3).
- **Light cone [derived, §X.1].** A 2 eV exciton has $q_c\approx10^{-3}$ Å⁻¹, while an exciton of radius ~1 nm (illustrative) samples $Q\sim0.1$ Å⁻¹. The kernel stays electrostatic; radiative decay and polaritons live in the transverse sector of $D$.
- **For your code.** Your BSE already builds the photon-irreducible response. EPP would combine it with phonons and light: IR and UV $\epsilon(\omega)$ from one object, polaritons and radiative widths from $D_T$, and BSE-quality $\epsilon_\infty$, $r_0$ for E.454 and E.472.

## 3. What standard theory does

- **Giustino.** $P=GG\Gamma$ (G.89–G.92) contains excitons through the vertex, but the review uses RPA $P_{\rm e}$; excitons appear only as BSE plus Allen–Heine (Marini [7]; G p. 42). Absorption is Hall–Bardeen–Blatt (G.200): classical field, independent particles, $n_r$ put in by hand; first principles in [8]. Temperature enters through Williams–Lax (G.201), stochastically (Zacharias, Patrick, Giustino 2015, G p. 45) or with one deterministic configuration [9].
- **BSE [est].** $L=L^0+L^0(2v-W)L$. The absorption uses $\bar L$, with the long-range $v_{\bG=0}$ removed from the exchange, and $\epsilon_M=1-v_0\bar\chi_{00}$. The full $L$ gives the loss function [1].
- **TDDFT [est].** $\chi=\chi^0+\chi^0(v+f^{xc})\chi$, with local fields from $\bG\neq0$. A long-range kernel $-\alpha/q^2$ mimics continuum excitons [3].
- **Exciton–phonon [est].** Exciton self-energy [13]; exciton–phonon matrix elements and relaxation times (hBN) [11]; "elemental" (irreducible) versus reducible excitons [14]; exciton lifetimes and optical line-width profile from exciton–phonon interactions, monolayer MoS₂ [15]; phonon-assisted absorption and emission of bulk hBN by finite differences [10].
- **Polaritons [est].** Hopfield [2] mixes the exciton field with photons. In a BSE code this is the transverse pole of §4.1.
- **2D [est].** For a freestanding sheet $W^{\rm 2D}=2\pi/[Q(1+r_0Q)]$ with $r_0=2\pi\alpha_{\rm 2D}$ [4]. In a dielectric environment the standard Rytova–Keldysh generalization is $2\pi/[Q(\epsilon_{\rm env}+r_0Q)]$, used with ab initio $r_0$ for TMD excitons in [17]. EPP's E.469–E.470 give the sheet part; $\epsilon_{\rm env}$ is added in [§B.10](../../sections/B10-2d.html).

## 4. How it could enter the EPP formalism

**4.1 Photon-irreducible BSE as $\Pi^{\rm e}$ [proposal].** With $t=(n\bk,m\bk+\bq)$ and $L^0_t=(f_{n\bk}-f_{m\bk+\bq})/(\hbar\omega+E_{n\bk}-E_{m\bk+\bq}+i\eta)$ as in E.282:

$$\Pi^{\rm e,BSE}_{\mu\nu}(\bq,\omega)=-\frac{2}{c\Omega}\sum_{tt'}j^{\mu}_{t}\big[\bar L_{tt'}(\bq,\omega)-\bar L_{tt'}(\bq,0)\big]\,j^{\nu*}_{t'},\qquad \bar L=L^0+L^0\big(2\bar K^{x}-K^{d}\big)\bar L .$$

*Kernels.* Write all kernels in EPP's scalar + vector (4×4) form of E.317–E.319, $D^0={\rm diag}(v,\hat D^{0,t})$, which is Coulomb-gauge-like. The exchange is $\bar K^{x}_{tt'}=\rho_t^*\,\bar v\,\rho_{t'}$ ($\bG\neq0$ only), plus a transverse (magnetic) part $\bj_T^*\cdot D^{0t}\cdot\bj_T$. That part is of relative size $(\hbar|\bq+\bG|/m_ec)^2\approx10^{-5}$–$10^{-4}$ (§B.5: Biot–Savart is $(v/c)^2$ weaker than Coulomb), and we drop it [derived]. The temporal-gauge form $-\frac1c\,\bj^*\!\cdot D^{0l}\cdot\bj$ gives $(\Delta_t\Delta_{t'}/\hbar^2\omega^2)\,\bar v_{tt'}$ by continuity and equals $\rho^*\bar v\rho$ only on shell ([selfenergy.md](selfenergy.md) §4.3), so it is not used inside the kernel. The current form is kept only for the two external vertices of $\Pi^{\rm e}$. The direct kernel is $K^d=\rho^*D_{11}\rho$, with $D_{11}=W$ the scalar (Coulomb) sector of the total $D$ (§4.3). With $\bar K^x=K^d=0$ this reduces to E.281. The subtraction $\bar L(\omega)-\bar L(0)$ stands in for the diamagnetic term. This is exact only if the static BSE current response obeys the same Ward identity, which holds for the conserving pair of §7 step 1 but not in general for $G_0W_0$ + static-$W$ BSE.

*Length form.* At $\bq\to0$ we define $\Pi^{\rm e,BSE}$ by its length-gauge form (what BSE codes compute):

$$\Pi^{\rm e,BSE}_{\mu\nu}=\frac{\omega^2}{4\pi c}\big[\epsilon^{\bar L}_{\mu\nu}(\omega)-\delta_{\mu\nu}\big],\quad \epsilon^{\bar L}_{\mu\nu}=\delta_{\mu\nu}+\frac{8\pi e^2}{\Omega}\sum_\lambda\Big[\frac{R^\mu_\lambda R^{\nu*}_\lambda}{E_\lambda-\hbar\omega-i\eta}+\frac{R^{\mu*}_\lambda R^{\nu}_\lambda}{E_\lambda+\hbar\omega+i\eta}\Big],$$

with $R^\mu_\lambda=\sum_{vc\bk}A^\lambda_{vc\bk}\langle c\bk|\hat v_\mu|v\bk\rangle/(i\omega^{\rm KS}_{cv\bk})$ (TDA; standard [1]; for $\omega^{\rm KS}$ see §5). The velocity-gauge expression above agrees with it only if the current vertex is dressed consistently with the kernel and $\Sigma$, $\hat{\bf v}=\hat{\bf p}/m+\frac{i}{\hbar}[V_{\rm NL}+\Sigma_{\rm nl},\hat\br]$ (Ward identity). With bare vertices the two differ. Since $\epsilon^{\bar L}(0)$ is finite, $\Pi\propto\omega^2$, so gauge invariance ($\Pi(0)=0$) holds by construction for the length form. The independent check, the f-sum rule ($\Pi\to\Pi^{\rm dia}$ as $\omega\to\infty$), is then not guaranteed, and TDA violates it. Inserted into E.380 for a cubic crystal (§B.6):

- $\epsilon^{\bar L}_L=0$ gives the longitudinal excitons, which are the poles of the full-$L$ BSE [est].
- $c^2q^2=\omega^2\epsilon^{\bar L}_T$ gives exciton polaritons [2].
- For one exciton, $\omega_L^2-\omega_T^2=F/\epsilon_b$, the algebra of $\omega_{\rm LO}^2-\omega_{\rm TO}^2$. The exciton L–T splitting comes from the removed $\bG=0$ exchange, which $D^{0l}$ restores here [est; reproduced here in EPP form; 1, 2].

**4.2 Local fields and xc [proposal].** Fold the microscopic lines into a macroscopically irreducible self-energy, then keep E.380 at $\bG=\bG'=0$:

$$\bar\Pi^{\rm e}=\Pi^{\rm e,0}+\Pi^{\rm e,0}\big(\tilde D^0+F^{xc}\big)\bar\Pi^{\rm e},\quad \tilde D^0_{\bG\bG'}=(1-\delta_{\bG0})D^0(\bq+\bG,\omega)\delta_{\bG\bG'},\quad F^{xc}_{\mu\nu}=-\frac{c}{\omega^2}\partial_\mu\partial'_\nu f^{xc}.$$

$F^{xc}$ is the E.364 rule applied to $f^{xc}$. At this response-function level the current form is safe: both vertices carry the external $\omega$, and with the diamagnetic term in $\Pi^{\rm e,0}$ the longitudinal projection reproduces the density Dyson equation. The on-shell problem of §4.1 arises only between transition vertices [derived]. The long-range kernel $-\alpha/q^2$ [3] becomes $F^{xc}_{00}=+\frac{c\alpha}{\omega^2}\hat q\hat q$, a $q$-independent term that belongs in the kernel of the $\bar\Pi$ Dyson equation (the $\bG=\bG'=0$ element of $F^{xc}$), not in $D$ [derived].

**4.3 Direct kernel from $D$: phonon screening [proposal].** Take $K^d=\rho^*D_{11}\rho$, with $D_{11}=W$ the scalar (Coulomb) sector of the total $D=D^{\rm e}+D^{\rm ph}$, $D^{\rm ph}=-D^{\rm e}JFJ^\dagger D^{\rm e}$ (selfenergy.md §4.2). Its scalar sector is $W=W^{\rm e}+W^{\rm e}\Pi^{\rm ion}_{\rho\rho}W$, i.e. $W=W^{\rm e}+W^{\rm e}\mathcal D\,W^{\rm e}$, with $\mathcal D$ the reducible ionic density correlator whose poles are the screened LO frequencies of E.454 (G.104–105; §X.1). $\Pi^{\rm ion}$ itself is irreducible and built from bare short-range phonons, with poles at the unscreened $\omega_\sigma$. So EPP's kernel contains Fröhlich phonon screening [12], dynamically. In 2D the static sheet gives [derived from E.470–E.472]

$$\epsilon_{\rm 2D}(Q,0)=\epsilon_{\rm env}+\big(r_0^{\rm e}+r_0^{\rm ion}\big)Q,\qquad r_0^{\rm ion}=\frac{2\pi e^2Z^{*2}}{S_{\rm uc}\,\mu\,\omega_{\rm T}^2}=\frac{C}{(\hbar\omega_{\rm T})^2}.$$

With the primer's illustrative presets: 2.9 Å against $r_0^{\rm e}=10$ Å for hBN (+29%), 0.5 Å against 40 Å for a MoS₂-like sheet (~1%). Unlike bulk, 2D ionic screening vanishes as $Q\to0$: it changes $r_0$, not the $1/(\epsilon_{\rm env}r)$ tail. For $E_B\gg\hbar\omega_{\rm LO}$ the ions cannot follow, so the static value is an upper bound and the dynamic kernel is needed. Monolayer hBN is this case: $E_B$ is of order 2 eV (illustrative) against $\hbar\omega_{\rm LO}\lesssim0.2$ eV, so the static +29% badly overestimates the effect on the bound exciton.

**4.4 Phonon-assisted optics as the first non-RPA electronic insertion [proposal].** In EPP the phonon line is $D^{\rm ph}$. At first order in $D^{\rm ph}$, $\Pi^{\rm e}$ gains two self-energy diagrams and one vertex correction, $\Pi^{(1)}=\Pi^{\rm SE}[D^{\rm ph}]+\Pi^{\rm VC}[D^{\rm ph}]$, with the vertex $g$ of §B.8 (E.415–E.416). Golden-rule cuts through $F$ should reproduce the long-range (dipole) part of G.200, with $g\to$ the electrostatic limit of E.416 (E.448a) (to derive). The full G.200 needs the short-range vertex of [shortrange.md](shortrange.md). Ladders on both sides give the excitonic form [est, schematic; cf. 10]:

$${\rm Im}\,\epsilon^{(1)}_{\mu\mu}(\omega)\propto\frac{1}{\omega^2}\sum_{\lambda\bq\sigma s}\Big|\sum_{\lambda'}\frac{J^\mu_{\lambda'}\,\mathcal G_{\lambda'\lambda\sigma}(\bq)}{E_{\lambda'}-\hbar\omega}\Big|^2\big(n_{\bq\sigma}+\tfrac12+\tfrac s2\big)\,\delta\big(E_\lambda(\bq)+s\hbar\omega_{\bq\sigma}-\hbar\omega\big).$$

Here $J^\mu_{\lambda'}$ is the current of the bright exciton and $\mathcal G$ the exciton–phonon vertex; from EPP alone $\mathcal G$ contains only the Fröhlich part. The light is now internal: $\alpha=2\,{\rm Im}\,q$ with $c^2q^2=\omega^2\epsilon_T$, so $n_r$ becomes an output (§A.9). Luminescence needs the lesser $\Pi^<$ for a thermalized exciton population (detailed balance, as in [10]), which means a Keldysh extension of EPP. The same diagrams resummed in the exciton basis give $\Sigma^X_\lambda(\omega)$ for linewidths [13, 15].

**4.5 2D radiative widths and polaritons [proposal].** Replace $(\Omega_{\rm uc}/S_{\rm uc})\Pi^{\rm ion}$ in the finite-$c$ sheet algebra of §B.10 by the BSE sheet response. With $\beta=\sqrt{Q^2-\epsilon_{\rm env}\omega^2/c^2}$ one obtains exciton radiative widths with the $\sqrt{1-x^2}$ (p) and $1/\sqrt{1-x^2}$ (s) laws. Out-of-plane and interlayer excitons need the $zz$ block and the $z$-resolved $D(Q,\omega;z,z')$ (p. 100; [twod.md](twod.md)).

## 5. Pitfalls

- **$\bG=0$ exchange twice.** Keeping $v_0$ in the kernel while E.380 adds $D^{0l}$ doubles the L–T splitting. The transverse $\bG=0$ exchange belongs to $D$ only.
- **Temporal-gauge kernels.** Inside the BSE kernel use $\rho^*v\rho$ and $\rho^*W\rho$ (E.317–E.319), not $-\frac1c\bj^*\!\cdot D\cdot\bj$. The latter carries $(\Delta/\hbar\omega)^2$ factors and diverges for a static $W$ (§4.1).
- **Screening twice.** EPP needs bare short-range phonons (§Z.1). Phonons that already carry the Gonze–Lee LO–TO term, combined with $\Pi^{\rm ion}$, count it twice. $Z^*$ and $\epsilon_\infty$ from DFPT include local fields and xc, while EPP's $\Pi^{\rm e}$ does not: keep them consistent, or E.454 shifts. Restoring the ionic $\bG\neq0$ terms must match the Gonze–Lee Ewald split ([shortrange.md](shortrange.md)). Williams–Lax supercells already contain the static polar fields that $D^{\rm ph}$ would add again.
- **Polar-only phonon channel.** Without short-range $g$, EPP gives zero phonon-assisted absorption in a nonpolar crystal such as Si ($Z^*=0$), and misses zone-boundary processes such as those near T/K in bulk hBN.
- **Which excitons couple.** Compute $\mathcal G$ from $\bar L$ eigenvectors and let $D$ add the macroscopic field. In EPP terms, Paleari–Marini's elemental versus reducible distinction [14] is $\bar L$ versus $L$ [spec].
- **Vertex screening.** DFPT $g$ is screened by $\epsilon^{{\rm H}xc}$, EPP's by $D^{\rm e}$: use one. Self-energies and phonon kernel must be a $\Phi$-derivable pair.
- **Gauge and sum rules.** E.350 (length = velocity gauge) holds only if the current vertex is dressed consistently with the kernel and $\Sigma$ (Ward identity). Include $[V_{\rm NL},\br]$ and, with $GW$, the nonlocal-$\Sigma$ correction to velocities [est] (Del Sole–Girlanda 1993; Levine–Allan 1989; both unverified, §6). For $G_0W_0$ + static-$W$ BSE no exact Ward identity exists: the kernel is not $\delta\Sigma_{GW}/\delta G$ (the $iG\,\delta W/\delta G$ term is dropped), and the static approximation breaks the frequency structure. The length–velocity mismatch must then be measured, not assumed away.
- **KS versus QP energies in $R$.** $R^\mu_\lambda$ must use KS transition energies $\omega^{\rm KS}_{cv}$ when $\langle\hat v\rangle$ is the KS velocity (i.e. KS dipoles $\langle c|\br|v\rangle$), even when $E_\lambda$ uses QP energies. Using $\omega^{\rm QP}$ rescales oscillator strengths by $(\omega^{\rm KS}/\omega^{\rm QP})^2$ [est].
- **f-sum rule and TDA.** Check the f-sum rule with the full (coupling) BSE or at the IP/RPA level; TDA violates it. Beyond TDA, keep the coupling block for $D_L$ near plasmons [est: TDA loss spectra are unreliable there, e.g. the Si plasmon] (Olevano–Reining 2001; Sander–Maggio–Kresse 2015; both unverified, §6). Keep SE and VC together, or $\Pi^{(1)}(0)\neq0$ gives a fake Drude term (§B.5).
- **Numerics.** Compute $\Pi^{\rm e}$ from $\epsilon^{\bar L}-1$ times $\omega^2$, since the $1/\omega^2$ double pole of E.374 cancels only analytically. 2D BSE needs truncated Coulomb and dense $\bk$ meshes. $\eta$ must be below meV [6], and Fröhlich $q\to0$ singularities need principal-value care (G p. 43).
- **$n_B=0$** fails for internal phonon lines; use Matsubara or real-time occupations.

## 6. State of the art (verified)

Bibliographic data (authors, journal, volume, article number, year) of all entries below were confirmed: refs 1–15 by an earlier search pass, refs 1, 7, 8, 9 and 16 again against Giustino's reference list, ref 17 in the 2D topic note's verified list. DOIs marked † follow APS's volume.article pattern and were not displayed in search results.

1. G. Onida, L. Reining, A. Rubio, Rev. Mod. Phys. **74**, 601 (2002), doi:10.1103/RevModPhys.74.601. BSE versus TDDFT, $\bar v$, $\epsilon_M$.
2. J. J. Hopfield, Phys. Rev. **112**, 1555 (1958), doi:10.1103/PhysRev.112.1555†. Exciton polaritons.
3. L. Reining, V. Olevano, A. Rubio, G. Onida, Phys. Rev. Lett. **88**, 066404 (2002), doi:10.1103/PhysRevLett.88.066404†. xc kernel from BSE, $-\alpha/q^2$.
4. P. Cudazzo, I. V. Tokatly, A. Rubio, Phys. Rev. B **84**, 085406 (2011), doi:10.1103/PhysRevB.84.085406. 2D screening of a freestanding sheet, $r_0$ from layer polarizability.
5. A. Chernikov *et al.*, Phys. Rev. Lett. **113**, 076802 (2014), doi:10.1103/PhysRevLett.113.076802. WS₂: $E_B=0.32$ eV, non-hydrogenic series.
6. F. Cadiz *et al.*, Phys. Rev. X **7**, 021026 (2017), doi:10.1103/PhysRevX.7.021026. ~2 meV linewidths in encapsulated MoS₂.
7. A. Marini, Phys. Rev. Lett. **101**, 106405 (2008), doi:10.1103/PhysRevLett.101.106405. BSE + Allen–Heine, finite-$T$ excitons.
8. J. Noffsinger *et al.*, Phys. Rev. Lett. **108**, 167402 (2012), doi:10.1103/PhysRevLett.108.167402. Ab initio indirect absorption in Si.
9. M. Zacharias, F. Giustino, Phys. Rev. B **94**, 075125 (2016), doi:10.1103/PhysRevB.94.075125†. One-shot Williams–Lax.
10. F. Paleari, H. P. C. Miranda, A. Molina-Sánchez, L. Wirtz, Phys. Rev. Lett. **122**, 187401 (2019), doi:10.1103/PhysRevLett.122.187401†. Phonon-assisted absorption and emission with excitons (bulk hBN).
11. H.-Y. Chen, D. Sangalli, M. Bernardi, Phys. Rev. Lett. **125**, 107401 (2020), doi:10.1103/PhysRevLett.125.107401. Exciton–phonon matrix elements and relaxation times.
12. M. R. Filip, J. B. Haber, J. B. Neaton, Phys. Rev. Lett. **127**, 067401 (2021), doi:10.1103/PhysRevLett.127.067401. Phonon screening in the BSE kernel (3D halide perovskites).
13. G. Antonius, S. G. Louie, Phys. Rev. B **105**, 085111 (2022), doi:10.1103/PhysRevB.105.085111†. Exciton–phonon self-energy.
14. F. Paleari, A. Marini, Phys. Rev. B **106**, 125403 (2022), doi:10.1103/PhysRevB.106.125403. Elemental versus reducible excitons.
15. Y.-h. Chan *et al.*, Nano Lett. **23**, 3971 (2023) (DOI not confirmed). Exciton lifetimes and optical line-width profile from exciton–phonon interactions, monolayer MoS₂.
16. G. Antonius, S. Poncé, P. Boulanger, M. Côté, X. Gonze, Phys. Rev. Lett. **112**, 215501 (2014), doi:10.1103/PhysRevLett.112.215501†. Many-body ($GW$) effects on the zero-point renormalization; diamond direct gap −628 meV (G Fig. 13).
17. T. C. Berkelbach, M. S. Hybertsen, D. R. Reichman, Phys. Rev. B **88**, 045318 (2013), doi:10.1103/PhysRevB.88.045318. Keldysh-screened excitons and trions in TMD monolayers with ab initio $r_0$.

**Unverified (from memory; not counted; confirm before citing).** No search budget was available in the last pass, and publisher and Crossref hosts were blocked.

- V. Olevano, L. Reining, Phys. Rev. Lett. **86**, 5962 (2001). Excitonic effects on the Si plasmon; coupling block beyond TDA.
- T. Sander, E. Maggio, G. Kresse, Phys. Rev. B **92**, 045209 (2015). BSE beyond TDA for extended systems.
- R. Del Sole, R. Girlanda, Phys. Rev. B **48**, 11789 (1993). Optical properties with quasiparticle (nonlocal-$\Sigma$) velocities.
- Z. H. Levine, D. C. Allan, Phys. Rev. Lett. **63**, 1719 (1989). Linear optical response of Si and Ge with self-energy (scissor) corrections.

## 7. First steps

1. **Derive (1–2 weeks).** Write 4.1–4.2 with all prefactors, with kernels in the 4×4 scalar + vector form. Prove E.350 for the conserving pair $\Sigma_{\rm SEX}[W_{\rm static}$ fixed$]$ with $K=2\bar v-W_{\rm static}$ and the vertex $\hat{\bf v}=\hat{\bf p}/m+\frac{i}{\hbar}[V_{\rm NL}+\Sigma_{\rm SEX},\hat\br]$. For $G_0W_0$ + BSE, quantify the length–velocity mismatch numerically in your code. Show that the poles of $D_L$ are the full-$L$ excitons and that $D_T$ gives Hopfield polaritons.
2. **2D optics in your code (1–2 weeks).** This assumes that 2D-truncated Coulomb, SOC and $\bk$ meshes of at least 30×30 are already available. For monolayer hBN and encapsulated MoS₂ (with SOC to resolve A and B; encapsulation via $\epsilon_{\rm env}$ in the sheet propagator, not explicit hBN layers), build the sheet $\Pi$ from $\epsilon^{\bar L}$. (a) Compare $r_0$ from BSE, RPA and RPA without local fields (what E.472 inherits). (b) Compute radiative widths and polariton dispersion from the sheet propagator. (c) Check length versus velocity gauge. Check the f-sum rule with the full (coupling) BSE or at the IP/RPA level, since TDA violates it.
3. **Phonon screening (weeks).** Use $D_{11}=W^{\rm e}+W^{\rm e}\Pi^{\rm ion}_{\rho\rho}D_{11}$ (bare short-range phonons, $Z^*$) as the kernel, first static and then with a first-order dynamic correction [12]. Test bulk MgO or GaN against $\epsilon_0-\epsilon_\infty$ (Lyddane–Sachs–Teller) and the model of [12]; in 2D compare with $r_0^{\rm ion}$. For monolayer hBN go directly to the dynamical kernel, since $E_B\gg\hbar\omega_{\rm LO}$.
4. **Phonon-assisted optics (months).** Derive $\Pi^{(1)}$ at finite $T$ and recover the Fröhlich part of G.200, then the excitonic form. Benchmark the polar part on a polar crystal (e.g. GaN or MgO), then bulk hBN [10] once short-range $g$ is added ([shortrange.md](shortrange.md)). A Si comparison [8] requires the short-range vertex first.
5. **Local fields and xc (weeks).** Restore $\bG,\bG'$ and $F^{xc}$; quantify $\epsilon_\infty$ and E.454.

## 8. Open questions

- In a supercell sheet, should $\bar v$ remove only $\bG=0$ or the whole $\bG_\parallel=0$ column, given that EPP's sheet photon line integrates over $q_z$? [spec]
- With a dynamical $D$ in the kernel, the BSE becomes frequency dependent. Is there a conserving, gauge-invariant $\Phi$-derivable truncation beyond the static SEX pair of §7 step 1?
- Does EPP's structure prove Paleari–Marini's claim [14]?
- In 2D, how do radiative and phonon widths compare as $\epsilon_{\rm env}$, temperature and polar-substrate (remote-phonon) screening vary?
- Does a finite-$q$ $\bar L$ inside E.380 settle the additional-boundary-condition problem of bulk exciton polaritons?

## 9. Review responses

All 16 points were checked against the EPP text (E.293–E.295, E.316–E.319, E.380, E.432–E.445, E.454, p. 100 and E.461), Giustino (G.103–G.105, Fig. 13 and the reference list), selfenergy.md §4.3 and primer §B.5 and §Z.1. All were accepted.

1. **Temporal-gauge BSE kernels (major).** Accepted. §4.1 now writes $\bar K^x=\rho^*\bar v\rho$ and $K^d=\rho^*D_{11}\rho$ in the E.317–E.319 form. It states that $-\frac1c\bj^*D^{0l}\bj$ equals Coulomb only on shell and keeps the current form only for the external vertices. A pitfall was added to §5. §4.2 gains a note on why the current form is safe there (response level, external $\omega$ at both vertices).
2. **Phonon-screened $W$ (major).** Accepted. §4.3 now gives the Dyson form $W=W^{\rm e}+W^{\rm e}\Pi^{\rm ion}_{\rho\rho}W$, its equivalent $W^{\rm e}+W^{\rm e}\mathcal DW^{\rm e}$ (G.104–105), and the irreducible/reducible pole structure. Step 3 uses the same Dyson form.
3. **Length vs velocity Π (major).** Accepted. §4.1 now defines $\Pi^{\rm e,BSE}$ by its length form, says when the velocity form agrees (dressed vertex, Ward identity) and qualifies the $\bar L(0)$ subtraction. The f-sum rule is flagged as not guaranteed in TDA (§4.1, §5, step 2). A §5 item on $\omega^{\rm KS}$ in $R$ was added.
4. **Short-range coupling for phonon-assisted optics (major).** Accepted. §1 has a new bullet. §2 (linewidths), §4.4, a new §5 item and step 4 now restrict EPP to the polar part of G.200. Si needs the short-range vertex first; the benchmark order is polar crystal, then bulk hBN.
5. **No Ward identity for $G_0W_0$ + BSE (major).** Accepted. Step 1 now proves E.350 only for the conserving static-SEX pair and measures the mismatch for $G_0W_0$ + BSE. §5 explains why. The overview first step mirrors this.
6. **arXiv:2503.13613 and the reference count (major).** Accepted. The preprint was dropped: no search budget remained, arxiv.org and publisher hosts were blocked, and its authors and content could not be checked. The verified count in the overview is now 10 fully confirmed entries (refs 1, 4–8, 11, 12, 14, 17), with 6 DOIs inferred (†: 2, 3, 9, 10, 13, 16) and the DOI of [15] unconfirmed.
7. **Transverse-exchange size (minor).** Accepted and rechecked: $\hbar^2|\bG|^2/(m_e\Delta)\approx15$ for $|\bG|=2$ Å⁻¹ and 2 eV, and $(\hbar|\bG|/m_ec)^2\approx6\times10^{-5}$. §4.1 now quotes $10^{-5}$–$10^{-4}$ and cites §B.5.
8. **Attributions to EPP (minor).** Accepted. §1 bullets 1–2 were rewritten as proposed. The phrase "constant in $GG'$ space" after E.293 was confirmed in the EPP text.
9. **2D bullet (minor).** Accepted. The bullet now separates the screening EPP has (E.469–E.470) from the e–h ladder it lacks. The quotation "including electron $z$-LFE" (p. 100) was confirmed.
10. **"Mixed insertions" (minor).** Accepted. These are now called non-RPA electronic insertions with an internal $D^{\rm ph}$ line, and a sentence on irreducibility and double counting was added. §4.4 was retitled.
11. **Diamond ZPR (minor).** Accepted. The text now says direct gap and −628 meV, and Antonius *et al.* 2014 was added as [16]; its data were confirmed in Giustino's reference list (the DOI follows the APS pattern, †).
12. **Content claims for [11], [12], [15] (minor).** Accepted. [11] is described only as matrix elements and relaxation times (hBN). [15] has one description, close to its title. For [12], the quoted "significantly" and the "CsPbX₃" list were removed in favour of "halide perovskites", which is in the title. A new search could not be run (budget exhausted).
13. **[spec] on TDA and the unsourced Σ correction (minor).** Accepted in substance: the label is now [est]. The four suggested papers are named in the text but listed under "Unverified" in §6, not as verified references, because they could not be searched in this pass.
14. **Cudazzo vs environment form (minor).** Accepted. [4] is now attributed the freestanding form only. For the $\epsilon_{\rm env}$ generalization I did not attribute the form itself to Berkelbach *et al.*, since I could not confirm that their paper includes $\epsilon_{\rm env}$; [17] (verified in the 2D note) is cited for the ab initio $r_0$ use.
15. **L–T splitting label; $F^{xc}$ placement (minor).** Accepted. The label is now [est; reproduced here in EPP form; 1, 2]. Hopfield [2] is added because the L–T splitting is part of the polariton model. In 4.2, $F^{xc}_{00}$ is placed in the kernel of the $\bar\Pi$ Dyson equation.
16. **Steps 2–3 (minor).** Accepted. Step 2 now includes the SOC, $\epsilon_{\rm env}$ and infrastructure assumptions. Step 3 now uses LST for bulk and $r_0^{\rm ion}$ only in 2D, and goes straight to the dynamical kernel for monolayer hBN. §4.3 explains why.
