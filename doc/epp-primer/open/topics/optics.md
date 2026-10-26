# Optics beyond RPA: excitons, local fields, xc and phonon-assisted absorption

Labels: **[est]** established literature; **[derived]** worked out here from EPP's equations, not printed in EPP; **[proposal]** suggested extension; **[spec]** speculation. Units are Gaussian, as in EPP. "G" is Giustino's review, "E" the EPP draft.

## 1. What is missing in EPP

- **Independent-particle bubble.** E.278–E.279 (pp. 57–58) define the photon self-energy as "one bubble plus one contact term" per species. E.282 and E.383 (pp. 59, 85) evaluate it with single-particle energies $E_{n\bk}$ and current vertices (E.284, E.384). No electron–hole ladder means no excitons; no self-energy insertions means no $GW$ gap and no phonon renormalization (primer [§B.5](../../sections/B05-photon-se.html)).
- **No local fields, no xc.** From E.380 on (p. 85) only $\bG=\bG'=0$ is kept, for electrons and, through the structure factors of E.295′, for ions. No $f^{xc}$ appears anywhere. Hence $W^{\rm e}=D^{\rm e}_{11}$ (E.437) is the head of the RPA dielectric matrix, and the $\epsilon_\infty$ in E.454 and the $r_0$ in E.472 are RPA numbers without local fields. DFPT screening (G.54) is static but includes both ([§A.3](../../sections/A03-dfpt.html), [§B.6](../../sections/B06-propagator.html)).
- **No mixed electron–ion insertions.** Mixed bubbles vanish at lowest order (E.316), and EPP stops there. Phonon-assisted absorption (G.200) and exciton–phonon coupling are the first mixed insertions: an ionic line inside the electronic bubble. They are not derived ([§A.9](../../sections/A09-applications.html), [§Z.1](../../sections/Z01-outlook.html)).
- **$T=0$ boson lines.** The $n_B=0$ rule (E.249, §B.5) is safe for E.379. It fails once a phonon line sits inside $\Pi^{\rm e}$, which needs the weights $n_{\bq\nu}+\frac12\pm\frac12$ of G.200.
- **2D.** The sheet bubble is $l\,\Pi^{\rm e}_{\bG=\bG'=0}$ (p. 60), strictly 2D. Layer excitons (Keldysh screening, interlayer pairs) are absent.

## 2. Why it matters

- **Excitons dominate 2D optics [est].** Monolayer WS₂ has an exciton binding energy of 0.32 eV and a non-hydrogenic Rydberg series [5]. An RPA $\Pi^{\rm e}$ misses the peaks that carry the oscillator strength.
- **Linewidths are meV [est].** hBN-encapsulated MoS₂ shows exciton linewidths of about 2 meV FWHM at 4 K, approaching the homogeneous limit [6]. What remains is radiative (transverse photon) plus exciton–phonon broadening, which are exactly the two channels EPP is built to unify.
- **Screening inputs propagate [est].** Giustino (§XII; primer §A.9) estimates that a 20% error in $\epsilon_\infty$ means about 40% in band renormalization and 10% in phonon frequencies. EPP's $\epsilon_\infty$ and $r_0$ currently lack local fields and xc.
- **Temperature [est].** Zero-point motion shifts and broadens excitons in Si and hBN [7]. The diamond zero-point gap renormalization is 628 meV with $GW$ (G Fig. 13), which moves every absorption edge.
- **Phonon screening of $W$ [est].** Adding phonon screening to the BSE kernel "significantly" lowers exciton binding energies in CsPbX₃ [12]. In EPP this screening is already inside $D$ (§4.3).
- **Light cone [derived, §X.1].** A 2 eV exciton has $q_c\approx10^{-3}$ Å⁻¹, while an exciton of radius ~1 nm (illustrative) samples $Q\sim0.1$ Å⁻¹. The kernel stays electrostatic; radiative decay and polaritons live in the transverse sector of $D$.
- **For your code.** Your BSE already builds the photon-irreducible response. EPP would combine it with phonons and light: IR and UV $\epsilon(\omega)$ from one object, polaritons and radiative widths from $D_T$, and BSE-quality $\epsilon_\infty$, $r_0$ for E.454 and E.472.

## 3. What standard theory does

- **Giustino.** $P=GG\Gamma$ (G.89–G.92) contains excitons through the vertex, but the review uses RPA $P_{\rm e}$; excitons appear only as BSE plus Allen–Heine (Marini [7]; G p. 42). Absorption is Hall–Bardeen–Blatt (G.200): classical field, independent particles, $n_r$ put in by hand; first principles in [8]. Temperature enters through Williams–Lax (G.201), stochastically (Zacharias, Patrick, Giustino 2015, G p. 45) or with one deterministic configuration [9].
- **BSE [est].** $L=L^0+L^0(2v-W)L$. The absorption uses $\bar L$, with the long-range $v_{\bG=0}$ removed from the exchange, and $\epsilon_M=1-v_0\bar\chi_{00}$. The full $L$ gives the loss function [1].
- **TDDFT [est].** $\chi=\chi^0+\chi^0(v+f^{xc})\chi$, with local fields from $\bG\neq0$. A long-range kernel $-\alpha/q^2$ mimics continuum excitons [3].
- **Exciton–phonon [est].** Exciton self-energy [13]; relaxation times and phonon-assisted PL, LO-dominated in hBN [11]; "elemental" (irreducible) versus reducible excitons [14]; absorption at second order with photon and phonon on an equal footing (Toyozawa–Hopfield), off-diagonal terms included [15]; phonon-assisted absorption and emission of bulk hBN by finite differences [10].
- **Polaritons [est].** Hopfield [2] mixes the exciton field with photons. A first-principles BSE with retarded electron–hole exchange appeared recently (Mauri et al., arXiv:2503.13613, preprint).
- **2D [est].** $W^{\rm 2D}=2\pi/[Q(\epsilon_{\rm env}+r_0Q)]$ with $r_0=2\pi\alpha_{\rm 2D}$ [4]. EPP's E.469–E.470 plus $\epsilon_{\rm env}$ ([§B.10](../../sections/B10-2d.html)).

## 4. How it could enter the EPP formalism

**4.1 Photon-irreducible BSE as $\Pi^{\rm e}$ [proposal].** With $t=(n\bk,m\bk+\bq)$ and $L^0_t=(f_{n\bk}-f_{m\bk+\bq})/(\hbar\omega+E_{n\bk}-E_{m\bk+\bq}+i\eta)$ as in E.282:

$$\Pi^{\rm e,BSE}_{\mu\nu}(\bq,\omega)=-\frac{2}{c\Omega}\sum_{tt'}j^{\mu}_{t}\big[\bar L_{tt'}(\bq,\omega)-\bar L_{tt'}(\bq,0)\big]\,j^{\nu*}_{t'},\qquad \bar L=L^0+L^0\big(2\bar K^{x}-K^{d}\big)\bar L .$$

The kernels follow from EPP's current-form interaction $K=-\frac1c\,\bj^*\!\cdot D\cdot\bj$ ([selfenergy.md](selfenergy.md) §4). The exchange is $\bar K^{x}=-\frac1c\,\bj^*\!\cdot[D^0-D^0_{\bG=0}]\cdot\bj\simeq\bar v$, because at $\bG\neq0$ the transverse part is $\sim(\omega/c|\bq+\bG|)^2\approx10^{-6}$ of the longitudinal one at 2 eV [derived]. The direct term $K^d$ is the screened $W$ (§4.3). With $\bar K^x=K^d=0$ this reduces to E.281. At $\bq\to0$, equivalently,

$$\Pi^{\rm e,BSE}_{\mu\nu}=\frac{\omega^2}{4\pi c}\big[\epsilon^{\bar L}_{\mu\nu}(\omega)-\delta_{\mu\nu}\big],\quad \epsilon^{\bar L}_{\mu\nu}=\delta_{\mu\nu}+\frac{8\pi e^2}{\Omega}\sum_\lambda\Big[\frac{R^\mu_\lambda R^{\nu*}_\lambda}{E_\lambda-\hbar\omega-i\eta}+\frac{R^{\mu*}_\lambda R^{\nu}_\lambda}{E_\lambda+\hbar\omega+i\eta}\Big],$$

with $R^\mu_\lambda=\sum_{vc\bk}A^\lambda_{vc\bk}\langle c\bk|\hat v_\mu|v\bk\rangle/(i\omega_{cv\bk})$ (TDA; standard [1]). Since $\epsilon^{\bar L}(0)$ is finite, $\Pi\propto\omega^2$ and gauge invariance holds by construction. Inserted into E.380 for a cubic crystal (§B.6):

- $\epsilon^{\bar L}_L=0$ gives the longitudinal excitons, which are the poles of the full-$L$ BSE [est].
- $c^2q^2=\omega^2\epsilon^{\bar L}_T$ gives exciton polaritons [2].
- For one exciton, $\omega_L^2-\omega_T^2=F/\epsilon_b$, the algebra of $\omega_{\rm LO}^2-\omega_{\rm TO}^2$. The L–T splitting is the removed $\bG=0$ exchange, now restored by $D^{0l}$ [derived].

**4.2 Local fields and xc [proposal].** Fold the microscopic lines into a macroscopically irreducible self-energy, then keep E.380 at $\bG=\bG'=0$:

$$\bar\Pi^{\rm e}=\Pi^{\rm e,0}+\Pi^{\rm e,0}\big(\tilde D^0+F^{xc}\big)\bar\Pi^{\rm e},\quad \tilde D^0_{\bG\bG'}=(1-\delta_{\bG0})D^0(\bq+\bG,\omega)\delta_{\bG\bG'},\quad F^{xc}_{\mu\nu}=-\frac{c}{\omega^2}\partial_\mu\partial'_\nu f^{xc}.$$

$F^{xc}$ is the E.364 rule applied to $f^{xc}$. The long-range kernel $-\alpha/q^2$ [3] becomes $+\frac{c\alpha}{\omega^2}\hat q\hat q$, a $q$-independent term that belongs in $\bar\Pi$, not in $D$ [derived].

**4.3 Direct kernel from $D$: phonon screening [proposal].** Take $K^d=-\frac1c\bj^*\!\cdot D\cdot\bj$ with the total $D=D^{\rm e}+D^{\rm ph}$, $D^{\rm ph}=-D^{\rm e}JFJ^\dagger D^{\rm e}$ (selfenergy.md §4.2). Its $ll$ sector is $W=W^{\rm e}+W^{\rm e}\Pi^{\rm ion}_{\rho\rho}W^{\rm e}$ (§X.1), so EPP's kernel contains Fröhlich phonon screening [12], dynamically. In 2D the static sheet gives [derived from E.470–E.472]

$$\epsilon_{\rm 2D}(Q,0)=\epsilon_{\rm env}+\big(r_0^{\rm e}+r_0^{\rm ion}\big)Q,\qquad r_0^{\rm ion}=\frac{2\pi e^2Z^{*2}}{S_{\rm uc}\,\mu\,\omega_{\rm T}^2}=\frac{C}{(\hbar\omega_{\rm T})^2}.$$

With the primer's illustrative presets: 2.9 Å against $r_0^{\rm e}=10$ Å for hBN (+29%), 0.5 Å against 40 Å for a MoS₂-like sheet (~1%). Unlike bulk, 2D ionic screening vanishes as $Q\to0$: it changes $r_0$, not the $1/(\epsilon_{\rm env}r)$ tail. For $E_B\gg\hbar\omega_{\rm LO}$ the ions cannot follow; the static value is an upper bound and the dynamic kernel is needed.

**4.4 Phonon-assisted optics as the first mixed insertion [proposal].** In EPP the phonon line is $D^{\rm ph}$. At first order in $D^{\rm ph}$, $\Pi^{\rm e}$ gains two self-energy diagrams and one vertex correction, $\Pi^{(1)}=\Pi^{\rm SE}[D^{\rm ph}]+\Pi^{\rm VC}[D^{\rm ph}]$, with the vertex $g$ of §B.8 (E.415–E.416). Golden-rule cuts through $F$ should reproduce G.200 (to derive). Ladders on both sides give the excitonic form [est, schematic; 10, 11, 13, 15]:

$${\rm Im}\,\epsilon^{(1)}_{\mu\mu}(\omega)\propto\frac{1}{\omega^2}\sum_{\lambda\bq\sigma s}\Big|\sum_{\lambda'}\frac{J^\mu_{\lambda'}\,\mathcal G_{\lambda'\lambda\sigma}(\bq)}{E_{\lambda'}-\hbar\omega}\Big|^2\big(n_{\bq\sigma}+\tfrac12+\tfrac s2\big)\,\delta\big(E_\lambda(\bq)+s\hbar\omega_{\bq\sigma}-\hbar\omega\big).$$

Here $J^\mu_{\lambda'}$ is the current of the bright exciton and $\mathcal G$ the exciton–phonon vertex. The light is now internal: $\alpha=2\,{\rm Im}\,q$ with $c^2q^2=\omega^2\epsilon_T$, so $n_r$ becomes an output (§A.9). Luminescence needs the lesser $\Pi^<$ for a thermalized exciton population (detailed balance, as in [10]), which means a Keldysh extension of EPP. The same diagrams resummed in the exciton basis give $\Sigma^X_\lambda(\omega)$ for linewidths [13, 15].

**4.5 2D radiative widths and polaritons [proposal].** Replace $(\Omega_{\rm uc}/S_{\rm uc})\Pi^{\rm ion}$ in the finite-$c$ sheet algebra of §B.10 by the BSE sheet response. With $\beta=\sqrt{Q^2-\epsilon_{\rm env}\omega^2/c^2}$ one obtains exciton radiative widths with the $\sqrt{1-x^2}$ (p) and $1/\sqrt{1-x^2}$ (s) laws. Out-of-plane and interlayer excitons need the $zz$ block and the $z$-resolved $D(Q,\omega;z,z')$ (p. 100; [twod.md](twod.md)).

## 5. Pitfalls

- **$\bG=0$ exchange twice.** Keeping $v_0$ in the kernel while E.380 adds $D^{0l}$ doubles the L–T splitting. The transverse $\bG=0$ exchange belongs to $D$ only.
- **Screening twice.** EPP needs bare short-range phonons (§Z.1). Phonons that already carry the Gonze–Lee LO–TO term, combined with $\Pi^{\rm ion}$, count it twice. $Z^*$ and $\epsilon_\infty$ from DFPT include local fields and xc, while EPP's $\Pi^{\rm e}$ does not: keep them consistent, or E.454 shifts. Restoring the ionic $\bG\neq0$ terms must match the Gonze–Lee Ewald split ([shortrange.md](shortrange.md)). Williams–Lax supercells already contain the static polar fields that $D^{\rm ph}$ would add again.
- **Which excitons couple.** Compute $\mathcal G$ from $\bar L$ eigenvectors and let $D$ add the macroscopic field. In EPP terms, Paleari–Marini's elemental versus reducible distinction [14] is $\bar L$ versus $L$ [spec].
- **Vertex screening.** DFPT $g$ is screened by $\epsilon^{{\rm H}xc}$, EPP's by $D^{\rm e}$: use one. Self-energies and phonon kernel must be a $\Phi$-derivable pair.
- **Gauge and sum rules.** E.350 (length = velocity gauge) holds only if the current vertex is dressed consistently with the kernel (Ward identity). Include $[V_{\rm NL},\br]$ and, with $GW$, the nonlocal-$\Sigma$ correction to velocities (standard; reference not verified here). Beyond TDA, keep the coupling block for $D_L$ near plasmons [spec: TDA loss spectra are unreliable there]. Keep SE and VC together, or $\Pi^{(1)}(0)\neq0$ gives a fake Drude term (§B.5).
- **Numerics.** Compute $\Pi^{\rm e}$ from $\epsilon^{\bar L}-1$ times $\omega^2$, since the $1/\omega^2$ double pole of E.374 cancels only analytically. 2D BSE needs truncated Coulomb and dense $\bk$ meshes. $\eta$ must be below meV [6], and Fröhlich $q\to0$ singularities need principal-value care (G p. 43).
- **$n_B=0$** fails for internal phonon lines; use Matsubara or real-time occupations.

## 6. State of the art (verified)

DOIs marked † follow APS's volume.article pattern and were not displayed in search results.

1. G. Onida, L. Reining, A. Rubio, Rev. Mod. Phys. **74**, 601 (2002), doi:10.1103/RevModPhys.74.601. BSE versus TDDFT, $\bar v$, $\epsilon_M$.
2. J. J. Hopfield, Phys. Rev. **112**, 1555 (1958), doi:10.1103/PhysRev.112.1555†. Exciton polaritons.
3. L. Reining, V. Olevano, A. Rubio, G. Onida, Phys. Rev. Lett. **88**, 066404 (2002), doi:10.1103/PhysRevLett.88.066404†. xc kernel from BSE, $-\alpha/q^2$.
4. P. Cudazzo, I. V. Tokatly, A. Rubio, Phys. Rev. B **84**, 085406 (2011), doi:10.1103/PhysRevB.84.085406. 2D screening, $r_0$ from layer polarizability.
5. A. Chernikov *et al.*, Phys. Rev. Lett. **113**, 076802 (2014), doi:10.1103/PhysRevLett.113.076802. WS₂: $E_B=0.32$ eV, non-hydrogenic series.
6. F. Cadiz *et al.*, Phys. Rev. X **7**, 021026 (2017), doi:10.1103/PhysRevX.7.021026. ~2 meV linewidths in encapsulated MoS₂.
7. A. Marini, Phys. Rev. Lett. **101**, 106405 (2008), doi:10.1103/PhysRevLett.101.106405. BSE + Allen–Heine, finite-$T$ excitons.
8. J. Noffsinger *et al.*, Phys. Rev. Lett. **108**, 167402 (2012), doi:10.1103/PhysRevLett.108.167402. Ab initio indirect absorption in Si.
9. M. Zacharias, F. Giustino, Phys. Rev. B **94**, 075125 (2016), doi:10.1103/PhysRevB.94.075125†. One-shot Williams–Lax.
10. F. Paleari, H. P. C. Miranda, A. Molina-Sánchez, L. Wirtz, Phys. Rev. Lett. **122**, 187401 (2019), doi:10.1103/PhysRevLett.122.187401†. Phonon-assisted absorption and emission with excitons (hBN).
11. H.-Y. Chen, D. Sangalli, M. Bernardi, Phys. Rev. Lett. **125**, 107401 (2020), doi:10.1103/PhysRevLett.125.107401. Exciton–phonon matrix elements and relaxation times.
12. M. R. Filip, J. B. Haber, J. B. Neaton, Phys. Rev. Lett. **127**, 067401 (2021), doi:10.1103/PhysRevLett.127.067401. Phonon screening in the BSE kernel.
13. G. Antonius, S. G. Louie, Phys. Rev. B **105**, 085111 (2022), doi:10.1103/PhysRevB.105.085111†. Exciton–phonon self-energy.
14. F. Paleari, A. Marini, Phys. Rev. B **106**, 125403 (2022), doi:10.1103/PhysRevB.106.125403. Elemental versus reducible excitons.
15. Y.-h. Chan *et al.*, Nano Lett. **23**, 3971 (2023) (DOI not confirmed). Second-order linewidths, MoS₂.

## 7. First steps

1. **Derive (1–2 weeks).** Write 4.1–4.2 with all prefactors. Prove E.350 for $\Pi^{\rm e,BSE}$ with a Ward-consistent vertex. Show that the poles of $D_L$ are the full-$L$ excitons and that $D_T$ gives Hopfield polaritons.
2. **2D optics in your code (1–2 weeks).** For monolayer hBN and encapsulated MoS₂, build the sheet $\Pi$ from $\epsilon^{\bar L}$. (a) Compare $r_0$ from BSE, RPA and RPA without local fields (what E.472 inherits). (b) Compute radiative widths and polariton dispersion from the sheet propagator. Check gauges and the f-sum rule.
3. **Phonon screening (weeks).** Use $D_{11}$ with $\Pi^{\rm ion}$ (bare phonons, $Z^*$) as the kernel, first static and then with a first-order dynamic correction [12]. Test bulk MgO or GaN, then monolayer hBN and MoS₂, against the $r_0^{\rm ion}$ estimates above.
4. **Phonon-assisted optics (months).** Derive $\Pi^{(1)}$ at finite $T$, recover G.200, then the excitonic form. Benchmark on bulk hBN [10, 11]; zone-boundary phonons need short-range coupling (shortrange.md).
5. **Local fields and xc (weeks).** Restore $\bG,\bG'$ and $F^{xc}$; quantify $\epsilon_\infty$ and E.454.

## 8. Open questions

- In a supercell sheet, should $\bar v$ remove only $\bG=0$ or the whole $\bG_\parallel=0$ column, given that EPP's sheet photon line integrates over $q_z$? [spec]
- With a dynamical $D$ in the kernel, the BSE becomes frequency dependent. Is there a conserving, gauge-invariant $\Phi$-derivable truncation?
- Does EPP's structure prove Paleari–Marini's claim [14]?
- In 2D, how do radiative and phonon widths compare as $\epsilon_{\rm env}$, temperature and polar-substrate (remote-phonon) screening vary?
- Does a finite-$q$ $\bar L$ inside E.380 settle the additional-boundary-condition problem of bulk exciton polaritons?
