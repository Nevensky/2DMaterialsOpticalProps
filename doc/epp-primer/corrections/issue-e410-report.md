# Issue E.410: the q = 0 upper phonon polariton, not the LO–TO splitting

Report behind `corrections/issue-e410.tex`. It covers the full derivation, the checks, the numbers, what the primer should change, and the responses to the two skeptical reviews (§8).

- Sources: the EPP draft v2, PDF pp. 32–35, 78–102, read as page images. Giustino, RMP 89, 015003 (2017), pp. 015003-21 to -24 and -34 to -35.
- Units: EPP from §4 on is Gaussian. Born charges $Z^*$ carry the dimension of charge ($Z^*_{\rm EPP}=|e|Z^*_{\rm G}$). Giustino is SI.
- Numerical check: `test/corr-e410/check_e410.py`. Its output is `test/corr-e410/results.txt`, with 45 checks, all PASS (32 before the review; §8 of the script was rewritten and §9–§10 added).

---

## 0. Verdict

- **The equations are right.** The algebra of E.406–E.410 is correct, and I checked every step against the PDF.
- **The reading is wrong.** The section title (§7.1, "Lowest self-energy correction and LO-TO splitting") and its first sentence ("completely analogous to the non-analytic dynamical matrix correction done by Quantum Espresso", p. 90) misread E.410.
  - E.410 is taken at $\mathbf q=0$ with finite $c$. There the vacuum propagator is isotropic (E.408), so every polar mode, transverse or longitudinal, shifts by the same amount. There is no LO–TO splitting at that point.
  - Without electrons, that shift is unscreened.
  - E.410 is therefore the common $\mathbf q\to0$ frequency of the longitudinal mode and of the upper transverse phonon polariton (Huang) of an electron-free crystal.
- **The standard non-analytic term is a different limit.** QE's term (Giannozzi et al. 1991; Baroni et al. 2001; in the Gonze–Lee 1997 formulation) is $\propto|\hat{\mathbf q}\cdot\mathbf j_\sigma|^2/(\hat{\mathbf q}\cdot\hat\epsilon^\infty\cdot\hat{\mathbf q})$. It carries $1/\epsilon_\infty$ and leaves the TO modes alone. Its EPP counterpart is E.454, the $c\to\infty$-first limit, with electrons resummed into $D^{\rm e}$.
- **Two separate approximations, kept apart.**
  - *Order of limits.* The longitudinal term of the general pole condition contains no $c$. For the L mode, E.410 is the unscreened LO frequency at every $q$ and coincides with E.454 for $\epsilon_\infty=1$. The order of limits changes only the transverse modes: at $q\to0$ with finite $c$ they are degenerate with the L mode; for $cq\gg\omega$ they are unshifted.
  - *Screening.* The 56 meV vs 36.2 meV gap for GaAs comes from the electron-free approximation, which §7.1 states openly ("Suppose that we neglect the interaction with electrons"). It does not come from the order of limits. With electrons, the same $q\to0$, finite-$c$ limit gives 35.9 meV for all three polar modes.
  - The primer's erratum 47 already says this in one line: "The two agree only for a mode along $\hat{\mathbf q}$ with $\epsilon^\infty=1$."
- **One formula covers both limits.** The general-$\mathbf q$ result proposed in the task is correct as stated:
  $$\omega^2=\omega_\sigma^2+\frac{4\pi}{\Omega_{\rm uc}}\Big[\frac{|\hat{\mathbf q}\cdot\mathbf j_\sigma|^2}{\epsilon_\infty}+\frac{\omega^2|(1-\hat{\mathbf q}\hat{\mathbf q})\cdot\mathbf j_\sigma|^2}{\epsilon_\infty\omega^2-c^2q^2}\Big].$$
  I re-derived it from E.375, E.388, E.392, E.399 and E.405. It contains both limits.
- **Even with electrons, the $q\to0$ limit at finite $c$ is still isotropic.** All polar modes sit at $\omega_{\rm L}$ (35.9 meV for GaAs). Screening fixes the size of the shift. The order of limits decides whether the transverse modes share it.

---

## 1. What the manuscript says (verified on the page images)

| Location | Content |
|---|---|
| p. 33 (top, after E.152) | Ionic plasma frequency $\omega_{p\mu}^2=4\pi Z^{*2}_\mu/(\Omega_{\rm uc}M)$, with $M$ the reduced mass (p. 31). p. 32 only uses $\omega_{p\mu}$, in E.152. |
| pp. 33–34, E.158–E.160 | $\epsilon_\mu=\epsilon_\infty-\omega_{p\mu}^2/(\omega^2-\omega_{\mu{\rm O}}^2+\frac13\omega_{p\mu}^2)$; $\tilde\omega^2_{\rm TO}=\omega^2_{\rm O}-\frac13\omega_p^2$, $\tilde\omega^2_{\rm LO}=\omega^2_{\rm O}+(\frac1{\epsilon_\infty}-\frac13)\omega_p^2$ (bare $\omega_{\rm O}$ with a Lorentz term). |
| p. 82, E.372 | $D^{0,t}_{\mu\nu}(\mathbf q,\omega)=-\frac{4\pi c}{q^2}\frac{q^2\delta_{\mu\nu}-q_\mu q_\nu}{\omega^2-\omega_{\mathbf q}^2+i\eta}$ |
| p. 83, E.374–375 | $D^{0,l}=-\frac{4\pi c}{q^2}\frac{q_\mu q_\nu}{\omega^2}$; $D^0_{\mu\nu}=\frac{4\pi c}{\omega^2}\frac{-q_\mu q_\nu+\frac{\omega^2}{c^2}\delta_{\mu\nu}}{q^2-\frac{\omega^2}{c^2}-i\eta}$ (boxed). E.381 on p. 85 is identical. |
| p. 85, E.383 | $\Pi^{\rm e}_{jj\mu\nu}=\frac{2}{c\Omega}\sum_{\mathbf k nm}\frac{\hbar\omega}{E_{n\mathbf k}-E_{m\mathbf k+\mathbf q}}\frac{f_{n\mathbf k}-f_{m\mathbf k+\mathbf q}}{\hbar\omega+E_{n\mathbf k}-E_{m\mathbf k+\mathbf q}+i\eta}j^\mu j^{\nu*}$ |
| p. 86, E.385 | $\Pi^{\rm ion}_{jj\mu\nu}=-\frac{\omega^2}{c\Omega_{\rm uc}}\sum_\sigma\frac{j_{\mu\sigma}j^*_{\nu\sigma}}{\omega^2-\omega_\sigma^2+i\eta}$, with $j_{\mu\sigma}=\sum_{a\alpha}Z^*_{a\mu\alpha}e_{a\alpha\sigma}/\sqrt{m_a}$ (no basis phase; see the companion correction, Eq. (385′)) |
| p. 87, E.388 | $D^{\rm e}=D^0+D^0\Pi^{\rm e,0}_{jj}D^{\rm e}$ (boxed) |
| p. 87, E.390, E.392 | $F^0_\sigma=\frac{2\omega_\sigma}{\hbar}\frac{1}{\omega^2-\omega_\sigma^2+i\eta}$; $J_{\mu\sigma}=\sqrt{\frac{\hbar}{2c\Omega_{\rm uc}\omega_\sigma}}\,\omega\,j_{\mu\sigma}$ (boxed) |
| p. 89, E.398–403 | $F=F^0+F^0\Omega^{\rm e}F$; $\Omega^{\rm e}_{\sigma\sigma'}=-\sum_{\alpha\beta}J^*_{\alpha\sigma}D^{\rm e}_{\alpha\beta}J_{\beta\sigma'}$; split $\Omega^{\rm e}=\Omega^0+\Delta\Omega^{\rm e}$ |
| p. 90, E.404–405 | Diagonal approximation; $F_\sigma=\frac1\hbar\frac{2\omega_\sigma}{\omega^2-\omega_\sigma^2-2\omega_\sigma\Omega^{\rm e}_{\sigma\sigma}/\hbar}$ |
| p. 90, §7.1 | Title "Lowest self-energy correction and LO-TO splitting"; "This procedure is completely analogous to the non-analytic dynamical matrix correction done by Quantum Espresso." |
| p. 91, before E.406 and E.407 | "the photon self-energy is" (E.406 writes $\Omega^0_\sigma$) and "the new corrected photon energy satisfy". Both should say "phonon" (primer erratum 48). |
| p. 91, E.406–410 | $D^{\rm e}\to D^0$; E.407 is the self-consistent pole condition; E.408 $D^0(\mathbf q=0,\omega)=-\frac{4\pi c}{\omega^2}\delta$; E.409 $\Omega^0_{\sigma\sigma}(0)=\frac{2\pi\hbar}{\Omega_{\rm uc}\omega_\sigma}\sum_\alpha|j_{\alpha\sigma}|^2$; E.410 $\tilde\omega^2_\sigma(0)=\omega^2_\sigma(0)+\frac{4\pi}{\Omega_{\rm uc}}\sum_\alpha|j_{\alpha\sigma}(0)|^2$ |
| p. 96, E.431–432 | $\sum_\alpha\epsilon^{\rm e}_{\mu\alpha}D^{\rm e}_{\alpha\nu}=D^0_{\mu\nu}$, $\epsilon^{\rm e}_{\mu\nu}=\delta_{\mu\nu}-\sum_\alpha D^0_{\mu\alpha}\Pi^{\rm e,0}_{jj\alpha\nu}$ |
| p. 98, E.442–446 | $\lim_{c\to\infty}D^0=-\frac{4\pi c}{\omega^2}\frac{q_\mu q_\nu}{q^2}$; $W^{\rm e}=4\pi/(\mathbf q\cdot\hat\epsilon^{\rm e}\cdot\mathbf q)$; $\lim_{c\to\infty}D^{\rm e}=-\frac{4\pi c}{\omega^2}\frac{q_\mu q_\nu}{\mathbf q\cdot\hat\epsilon^{\rm e}\cdot\mathbf q}$ |
| p. 99, E.449–454 | $\lim_{c\to\infty}\Omega^{\rm e}_{\sigma\sigma}=\frac{\hbar W^{\rm e}}{2\Omega_{\rm uc}\omega_\sigma}|\mathbf q\cdot\mathbf j_\sigma|^2$; E.454 $\tilde\omega^2_\sigma=\omega^2_\sigma+\frac{W^{\rm e}(\mathbf q,0)}{\Omega_{\rm uc}}|\sum_\alpha q_\alpha j_{\alpha\sigma}|^2$. The sentence continues after the display: "which is very similar (probably identical) to Eq.2 published in T. Sohier, …, Nano Lett. 17, 3758 (2017)." |
| p. 100, E.455–E.459 and the $z$-dependent Dyson equation | The general 2D transcription at finite $c$: E.458 uses the full $D^{\rm e}(\mathbf Q,\omega,z{=}0,z'{=}0)$. |
| pp. 101–102, §8.1, E.460–E.472 | Electrostatic ($c\to\infty$) 2D results. |

**Where the QE comparison appears.** I searched the whole text dump (pdftotext) for "LO-TO", "Quantum Espresso", "non-analytic" and "splitting". The QE comparison appears once, on p. 90. "LO-TO" appears elsewhere only in the classical §4.2.1 (pp. 33–34, the "LO shift" discussion of E.159–E.165). That passage is electrostatic ($\nabla\times\mathbf E=0$, used on p. 34) and does not refer to E.410. The draft has no abstract, table of contents or conclusions that repeat the claim. **No later equation cites E.410.** References to (405)–(410) occur only on pp. 90–91, plus (405) on pp. 99 and 102.

---

## 2. Re-derivation from the manuscript's definitions

### 2.1 E.406–E.410 as printed (correct)
- **Vertices.** With E.392, $J^*_{\alpha\sigma}J_{\beta\sigma}=\frac{\hbar\omega^2}{2c\Omega_{\rm uc}\omega_\sigma}j^*_{\alpha\sigma}j_{\beta\sigma}$.
- **E.408.** At $\mathbf q=0$, E.375 gives $\frac{4\pi c}{\omega^2}\cdot\frac{(\omega^2/c^2)\delta}{-\omega^2/c^2}=-\frac{4\pi c}{\omega^2}\delta$. ✓
- **E.409.** $\Omega^0_{\sigma\sigma}=-\frac{\hbar\omega^2}{2c\Omega_{\rm uc}\omega_\sigma}\big(-\frac{4\pi c}{\omega^2}\big)\sum_\alpha|j_{\alpha\sigma}|^2=\frac{2\pi\hbar}{\Omega_{\rm uc}\omega_\sigma}\sum_\alpha|j_{\alpha\sigma}|^2$. It does not depend on $\omega$, so the self-consistency of E.407 is trivial. ✓
- **E.410.** $2\omega_\sigma\Omega^0_{\sigma\sigma}/\hbar=\frac{4\pi}{\Omega_{\rm uc}}\sum_\alpha|j_{\alpha\sigma}|^2$. ✓

So E.410 needs no algebraic correction. The tex only relabels it (410′), to make the order of limits explicit, and keeps the equation counter in step (§7).

### 2.2 Longitudinal–transverse decomposition of E.375
Let $\hat{\mathbf q}=\mathbf q/q$, $P^L=\hat{\mathbf q}\hat{\mathbf q}$ and $P^T=1-P^L$. Write $\delta=P^L+P^T$ in the numerator of E.375:
$$-q^2P^L+\tfrac{\omega^2}{c^2}(P^L+P^T)=-(q^2-\tfrac{\omega^2}{c^2})P^L+\tfrac{\omega^2}{c^2}P^T ,$$
$$D^0=-\frac{4\pi c}{\omega^2}P^L-\frac{4\pi c}{\omega^2-c^2q^2+i0}P^T .$$
- The longitudinal projection, $-4\pi c/\omega^2$, contains no $q$ and therefore no retardation. It equals E.374 projected.
- The transverse projection, $-4\pi c/(\omega^2-c^2q^2)$, tends to $-4\pi c/\omega^2$ for $cq\ll\omega$ (this is E.408) and to $4\pi/(cq^2)$, i.e. $O(1/c)$, for $cq\gg\omega$ (this gives E.442).

Inverse: $(D^0)^{-1}=-\frac{1}{4\pi c}\big[\omega^2P^L+(\omega^2-c^2q^2)P^T\big]=-\frac{1}{4\pi c}\big[\omega^2\,1-c^2(q^2 1-\mathbf q\mathbf q)\big]$. This is the temporal-gauge Maxwell operator.

### 2.3 Electrons as a constant $\epsilon_\infty$
For an insulator at $\hbar\omega\ll E_{\rm gap}$ and $q$ much smaller than the BZ, E.383 is $\propto\omega^2$ with a positive coefficient.
- Take an occupied–empty pair with gap $\Delta$. The two orderings give $\hbar\omega\big[\frac{1}{\Delta(\Delta-\hbar\omega)}-\frac{1}{\Delta(\Delta+\hbar\omega)}\big]=\frac{2\hbar^2\omega^2}{\Delta(\Delta^2-\hbar^2\omega^2)}>0$.
- This defines $\epsilon_\infty$ through
$$\Pi^{\rm e,0}_{jj\mu\nu}\simeq\frac{(\epsilon_\infty-1)\omega^2}{4\pi c}\delta_{\mu\nu}\qquad\big(\text{tensor: }\tfrac{\omega^2}{4\pi c}(\epsilon^\infty_{\mu\nu}-\delta_{\mu\nu})\big).$$
- Consistency with the manuscript: E.432 gives $\hat{\mathbf q}\cdot\epsilon^{\rm e}\cdot\hat{\mathbf q}=1-(-\frac{4\pi c}{\omega^2})\frac{(\epsilon_\infty-1)\omega^2}{4\pi c}=\epsilon_\infty$. E.445 then gives $W^{\rm e}(\mathbf q)=4\pi/(q^2\epsilon_\infty)$.
- Caution: E.432's $\epsilon^{\rm e}$ is not $\epsilon_\infty\delta$ in its transverse block, because it is built with the full $D^0$. Only its longitudinal head enters E.446.
- The ionic analogue, E.385 in the same projection, gives $\epsilon(\omega)=\epsilon_\infty+\omega_p^2/(\omega_\sigma^2-\omega^2)$ with $\omega_p^2=4\pi|\hat{\mathbf q}\cdot\mathbf j|^2/\Omega_{\rm uc}$. That is E.158 without the Lorentz term, as it should be: the Dyson route keeps only $\mathbf G=0$, and the $\mathbf G\neq0$ local field sits in the short-range force constants that define $\omega_\sigma$. So $\omega_\sigma(0)=\omega_{\rm TO}$ in this model (§4).

### 2.4 Solving E.388: (410a)
$D^{\rm e}=[(D^0)^{-1}-\Pi^{\rm e,0}]^{-1}$. With the model above,
$$(D^{\rm e})^{-1}=-\frac{1}{4\pi c}\big[\omega^2\hat\epsilon^\infty-c^2(q^21-\mathbf q\mathbf q)\big]\equiv-\frac{\Lambda}{4\pi c},\qquad D^{\rm e}=-4\pi c\,\Lambda^{-1}.$$
$\Lambda$ is the Maxwell matrix. It was called $M$ in the first draft; it is renamed because $M$ is the manuscript's reduced mass (pp. 31–33). For isotropic $\epsilon_\infty$ the $L$ and $T$ blocks decouple ($P^LP^T=0$):
$$D^{\rm e}_{\mu\nu}=-4\pi c\Big[\frac{\hat q_\mu\hat q_\nu}{\epsilon_\infty\omega^2}+\frac{\delta_{\mu\nu}-\hat q_\mu\hat q_\nu}{\epsilon_\infty\omega^2-c^2q^2}\Big].\tag{410a}$$
This is exactly the proposal in the task. Setting $\epsilon_\infty=1$ gives E.375 back.

### 2.5 Phonon self-energy: (410b)
E.399 on the diagonal with E.392 gives
$$\Omega^{\rm e}_{\sigma\sigma}=-\frac{\hbar\omega^2}{2c\Omega_{\rm uc}\omega_\sigma}\,\mathbf j_\sigma^*\cdot D^{\rm e}\cdot\mathbf j_\sigma=\frac{2\pi\hbar\omega^2}{\Omega_{\rm uc}\omega_\sigma}\,\mathbf j_\sigma^*\cdot \Lambda^{-1}\cdot\mathbf j_\sigma$$
$$\overset{\rm iso}{=}\frac{2\pi\hbar}{\Omega_{\rm uc}\omega_\sigma}\Big[\frac{|\hat{\mathbf q}\cdot\mathbf j_\sigma|^2}{\epsilon_\infty}+\frac{\omega^2|P^T\mathbf j_\sigma|^2}{\epsilon_\infty\omega^2-c^2q^2}\Big].\tag{410b}$$
- $\mathbf j_\sigma^*\cdot P^L\cdot\mathbf j_\sigma=|\hat{\mathbf q}\cdot\mathbf j_\sigma|^2$ and $\mathbf j_\sigma^*\cdot P^T\cdot\mathbf j_\sigma=\sum_\alpha|j_{\alpha\sigma}|^2-|\hat{\mathbf q}\cdot\mathbf j_\sigma|^2$. These hold for complex $\mathbf j_\sigma$ too, since $P^{L,T}$ are real symmetric.
- **Split check (E.400–403).** With E.403, $D^0+D^0\Pi^{\rm e}D^0=D^{\rm e}$, so $\Omega^0+\Delta\Omega^{\rm e}=-J^*D^{\rm e}J$ (E.399). (410b) at $\epsilon_\infty=1$ is $\Omega^0$ (E.406). The difference is $\Delta\Omega^{\rm e}$. At $\mathbf q=0$ it equals $\frac{2\pi\hbar}{\Omega_{\rm uc}\omega_\sigma}(\epsilon_\infty^{-1}-1)\sum_\alpha|j_{\alpha\sigma}|^2<0$: screening always softens the bare dipolar stiffening.

### 2.6 Pole condition (410c) and matrix form (410d)
- **Diagonal.** E.405 has poles at $\omega^2=\omega_\sigma^2+2\omega_\sigma\Omega^{\rm e}_{\sigma\sigma}(\mathbf q,\omega)/\hbar$. With (410b) this gives
$$\omega^2=\omega^2_\sigma(\mathbf q)+\frac{4\pi}{\Omega_{\rm uc}}\Big[\frac{|\hat{\mathbf q}\cdot\mathbf j_\sigma|^2}{\epsilon_\infty}+\frac{\omega^2|(1-\hat{\mathbf q}\hat{\mathbf q})\cdot\mathbf j_\sigma|^2}{\epsilon_\infty\omega^2-c^2q^2}\Big].\tag{410c}$$
  In general this is $\omega^2=\omega_\sigma^2+\frac{4\pi\omega^2}{\Omega_{\rm uc}}\mathbf j^*_\sigma\cdot \Lambda^{-1}\cdot\mathbf j_\sigma$.
- **Matrix.** From E.398, $F^{-1}=(F^0)^{-1}-\Omega^{\rm e}$ with $(F^0)^{-1}_{\sigma\sigma}=\hbar(\omega^2-\omega_\sigma^2)/(2\omega_\sigma)$. Multiply rows and columns by $\sqrt{2\omega_\sigma/\hbar}$. Since $\Omega^{\rm e}_{\sigma\sigma'}=\frac{2\pi\hbar\omega^2}{\Omega_{\rm uc}\sqrt{\omega_\sigma\omega_{\sigma'}}}\mathbf j_\sigma^*\cdot \Lambda^{-1}\cdot\mathbf j_{\sigma'}$, this gives
$$\det\Big[(\omega^2-\omega^2_\sigma)\delta_{\sigma\sigma'}-\frac{4\pi\omega^2}{\Omega_{\rm uc}}\mathbf j^*_\sigma\cdot \Lambda^{-1}\cdot\mathbf j_{\sigma'}\Big]=0.\tag{410d}$$
  - For $c\to\infty$, $\frac{4\pi\omega^2}{\Omega_{\rm uc}}\Lambda^{-1}\to\frac{4\pi}{\Omega_{\rm uc}}\frac{\hat{\mathbf q}\hat{\mathbf q}}{\hat{\mathbf q}\cdot\hat\epsilon^\infty\cdot\hat{\mathbf q}}$. The bracket then becomes the mode-basis form of "analytic dynamical matrix + non-analytic term", which QE diagonalizes.
  - (410d) reduces to (410c) when the off-diagonal elements vanish. That happens for a single polar mode, or for degenerate modes chosen with mutually orthogonal $\mathbf j_\sigma$ along or across $\hat{\mathbf q}$ (isotropic $\epsilon_\infty$).
  - Numerical check §6: a degenerate triplet in a random basis gives the correct roots $\{\omega_-,\omega_-,\omega_{\rm L},\omega_+,\omega_+\}$ from (410d). The diagonal approximation in that basis gives wrong upper roots (100.2–100.9 meV instead of 100.92 at $q=3q_c$).

### 2.7 The two limits
The longitudinal term of (410c) contains no $c$. The transverse term depends on $c$ and $q$ only through $cq$, so the limits do not commute.
- **$\mathbf q\to0$ at fixed $c$** ($cq\ll\sqrt{\epsilon_\infty}\omega$): $\Lambda\to\omega^2\hat\epsilon^\infty$.
$$\tilde\omega^2_\sigma=\omega_\sigma^2+\frac{4\pi}{\Omega_{\rm uc}}\mathbf j^*_\sigma\cdot(\hat\epsilon^\infty)^{-1}\cdot\mathbf j_\sigma\overset{\rm iso}{=}\omega_\sigma^2+\frac{4\pi}{\Omega_{\rm uc}\epsilon_\infty}\sum_\alpha|j_{\alpha\sigma}|^2 .\tag{410e}$$
  This is isotropic and does not depend on the direction of approach. For $\epsilon_\infty=1$ it is E.410.
- **$cq\gg\sqrt{\epsilon_\infty}\omega$ (or $c\to\infty$ first).** Use the Schur complement of the $T$ block, which grows as $c^2q^2$: $\Lambda^{-1}\to\frac{\hat{\mathbf q}\hat{\mathbf q}}{\omega^2\,\hat{\mathbf q}\cdot\hat\epsilon^\infty\cdot\hat{\mathbf q}}+O\big((cq)^{-2}\big)$.
$$\tilde\omega^2_\sigma(\mathbf q)=\omega_\sigma^2+\frac{4\pi}{\Omega_{\rm uc}}\frac{|\hat{\mathbf q}\cdot\mathbf j_\sigma|^2}{\hat{\mathbf q}\cdot\hat\epsilon^\infty\cdot\hat{\mathbf q}}=\omega^2_\sigma+\frac{W^{\rm e}(\mathbf q,0)}{\Omega_{\rm uc}}\Big|\sum_\alpha q_\alpha j_{\alpha\sigma}\Big|^2 .\tag{410f}$$
  - This is E.454 with $W^{\rm e}(\mathbf q,0)=4\pi/(\mathbf q\cdot\hat\epsilon^\infty\cdot\mathbf q)$, the same as E.446 with $\hat\epsilon^{\rm e}\to\hat\epsilon^\infty$.
  - The TO correction is $-\frac{4\pi}{\Omega_{\rm uc}}|P^T\mathbf j|^2\omega^2/(c^2q^2-\epsilon_\infty\omega^2)$, suppressed as $(q_c/q)^2$.
- **What differs between the two limits.** Only the transverse modes. For a mode with $\mathbf j_\sigma\parallel\hat{\mathbf q}$, (410e) and (410f) are the same number, and with $\epsilon_\infty=1$ both equal E.410.
- **Non-analyticity.** The direction dependence of (410f) at $\Gamma$ (the "non-analyticity") does not survive at finite $c$. Below $q_c$ it is replaced by the isotropic (410e).
- **Range of (410f).** For $q\gg q_c$ retardation is negligible, and E.454 applies within the manuscript's $\mathbf G=0$ treatment: on the diagonal for a single polar mode, and through the $c\to\infty$ limit of (410d) for mixed polar modes. Local fields ($\mathbf G\neq0$) give the full dipole–dipole sum, Eq. (454a) of the companion basis-phase correction.
- **Basis phase.** At finite $q$, $j_{\alpha\sigma}(\mathbf q)$ must carry $e^{-i\mathbf q\cdot\mathbf r^0_a}$ (companion correction, Eq. (385′)). Without it, (410c) gives an acoustic mode a spurious, cell-dependent L-term of order $q^2$, the same order as $\omega_{\rm ac}^2$ (§3.4 item 9).

### 2.8 Cubic diatomic crystal: Huang's polaritons, (410g)
- **Mode currents.** Take $Z^*_a=-Z^*_b=Z^*$ (isotropic), with eigenvectors $e_a=\sqrt{m_b/(m_a+m_b)}\,\hat{\mathbf e}$ and $e_b=-\sqrt{m_a/(m_a+m_b)}\,\hat{\mathbf e}$. Then $\mathbf j=Z^*\hat{\mathbf e}/\sqrt M$, so $\sum_\alpha|j_\alpha|^2=Z^{*2}/M$ and $\frac{4\pi}{\Omega_{\rm uc}}|\mathbf j|^2=\omega_p^2$ (defined at the top of p. 33). Choose $\hat{\mathbf e}\parallel\hat{\mathbf q}$ (L) or $\perp\hat{\mathbf q}$ (T1, T2).
- **L mode:** $\omega_{\rm L}^2=\omega_\sigma^2+\omega_p^2/\epsilon_\infty$ for every $q$ and $c$.
- **T modes:** $(\omega^2-\omega_\sigma^2)(\epsilon_\infty\omega^2-c^2q^2)=\omega_p^2\omega^2$. This is equivalent to $c^2q^2=\omega^2\epsilon(\omega)$ with $\epsilon(\omega)=\epsilon_\infty+\frac{\omega_p^2}{\omega_\sigma^2-\omega^2}=\epsilon_\infty\frac{\omega_{\rm L}^2-\omega^2}{\omega_\sigma^2-\omega^2}$, which is Huang's polariton dispersion. The roots are
$$\omega_\pm^2=\frac{B\pm\sqrt{B^2-4\epsilon_\infty c^2q^2\omega_\sigma^2}}{2\epsilon_\infty},\quad B=\epsilon_\infty\omega_{\rm L}^2+c^2q^2 .\tag{410g}$$
- **Limits of the branches:**
  - $q\to0$: $\omega_+\to\omega_{\rm L}$, and $\omega_-\approx cq/\sqrt{\epsilon(0)}\to0$, which is photon-like.
  - $cq\gg\sqrt{\epsilon_\infty}\omega_\sigma$: $\omega_-\to\omega_\sigma=\omega_{\rm TO}$, and $\omega_+\approx cq/\sqrt{\epsilon_\infty}$, which is photon-like.
- **LST:** $\epsilon(0)/\epsilon_\infty=\omega_{\rm L}^2/\omega_\sigma^2$.
- **Phonon weight.** The residue of $F_\sigma$ at each pole, in the variable $\omega^2$ and relative to the bare residue, is $Z_\pm=[1+\omega_p^2c^2q^2/(\epsilon_\infty\omega_\pm^2-c^2q^2)^2]^{-1}$, with $Z_++Z_-=1$.

### 2.9 Crossover
$$q_c=\sqrt{\epsilon_\infty}\,\omega/c,\qquad\omega\simeq\omega_\sigma .\tag{410h}$$
This is the light line in the medium at the phonon frequency. Without electrons it is $\omega_\sigma/c$, smaller by $\sqrt{\epsilon_\infty}$, and that is the crossover that applies to the electron-free (410′).

---

## 3. Checks

### 3.1 Dimensions (Gaussian)
| Quantity | Dimensions | Result |
|---|---|---|
| $e^2$ | energy·length | |
| $Z^*$ | $({\rm E\,L})^{1/2}$ | |
| $j$ | $Z^*\,{\rm M}^{-1/2}$ | |
| $\frac{4\pi}{\Omega_{\rm uc}}\lvert j\rvert^2$ | $\frac{{\rm E\,L}}{{\rm M\,L^3}}={\rm T^{-2}}$ | ✓ |
| $D^0$ | $c/\omega^2={\rm L\,T}$ | |
| $J$ | ${\rm E\,T^{1/2}L^{-3/2}M^{-1/2}}$ | |
| $J^*DJ$ | ${\rm E^2T^2L^{-2}M^{-1}}={\rm E}$ | ✓ (energy, as $\Omega^{\rm e}$ must be) |
| $\omega^2/(\epsilon_\infty\omega^2-c^2q^2)$ | dimensionless | ✓ |
| $q_c$ | ${\rm L}^{-1}$ | ✓ |

`check_e410.py` §3 does this bookkeeping in M, L, T exponents.

### 3.2 Limits
| Limit | Result | Script |
|---|---|---|
| $\epsilon_\infty=1$ | (410a) → E.375 | §1 |
| $\epsilon_\infty=1$, $q\to0$ | E.408, E.409, E.410 | §1, §2, §5 |
| $c\to\infty$ | (410a) → E.446 with deviation $\propto(\omega/cq)^2$; (410c) → E.452/E.454 | §1, §7 |
| Tensor $\hat\epsilon^\infty$, $c\to\infty$ | $4\pi\lvert\hat{\mathbf q}\cdot\mathbf j\rvert^2/(\Omega\,\hat{\mathbf q}\cdot\hat\epsilon^\infty\cdot\hat{\mathbf q})$, the standard non-analytic term | §7 |
| L mode | shift $=\omega_p^2/\epsilon_\infty$ at every $q$, $\omega$, $c$ | §5, §8 |

### 3.3 Agreement with Giustino, Verdi–Giustino, the non-analytic term and Huang
- **Giustino Eq. 137 (p. -23).** $\tilde\Omega^2=\omega^2+2\omega\Pi^{\rm NA}(\tilde\Omega)$ has the same structure as E.405/E.407, with $\Omega^{\rm e}=\hbar\Pi$.
- **Giustino Eqs. 125–128 (pp. -21, -22).** The LO–TO splitting sits in the adiabatic $\Pi^{\rm A}$ (DFPT force constants). For $c\to\infty$ and $\omega\to0$, (410f) is the dipolar part of $\Pi^{\rm A}$.
- **Giustino p. -35, §VI.A.3.** The polar Wannier strategy "is analogous to that in use for calculating the splitting between LO and TO phonons in polar materials (Gonze and Lee, 1997)". That analogy is built on the electrostatic kernel $(\mathbf q\cdot\mathbf Z\cdot\mathbf e)/(\mathbf q\cdot\epsilon^\infty\cdot\mathbf q)$, the kernel of (410f)/E.454, which has no transverse shift.
- **The non-analytic term in SI** (Giustino's dimensionless $Z_\kappa$; Giannozzi et al. 1991, Baroni et al. 2001, Gonze–Lee 1997): $\tilde C^{\rm NA}_{\kappa\alpha,\kappa'\beta}=\frac{4\pi}{\Omega}\frac{e^2}{4\pi\varepsilon_0}\frac{(\mathbf q\cdot\mathbf Z_\kappa)_\alpha(\mathbf q\cdot\mathbf Z_{\kappa'})_\beta}{\mathbf q\cdot\epsilon^\infty\cdot\mathbf q}$.
  - Projected on mode $\sigma$ (divide by $\sqrt{M_\kappa M_{\kappa'}}$ and contract with $e^*_\sigma,e_\sigma$), it gives $\frac{e^2}{\varepsilon_0\Omega}\frac{|\hat{\mathbf q}\cdot\sum_\kappa\mathbf Z_\kappa\mathbf e_{\kappa\sigma}/\sqrt{M_\kappa}|^2}{\hat{\mathbf q}\cdot\epsilon^\infty\cdot\hat{\mathbf q}}$.
  - This equals (410f) under $Z^*_{\rm EPP}=|e|Z_{\rm G}$ and $e^2\to e^2/4\pi\varepsilon_0$.
  - Numerically the same $\omega_p^2$ (script §4).
- **Verdi–Giustino, Giustino Eq. 180 (p. -35): a test of screening, not of the order of limits.** Take $\mathbf G=0$, overlap 1, $N_p=1$ and the LO mode of a diatomic, where $\sum_\kappa\mathbf q\cdot\mathbf Z_\kappa\cdot\mathbf e_\kappa/\sqrt{M_\kappa}=qZ/\sqrt M$:
  $$|g^{\mathcal L}|^2=\Big(\frac{e^2}{\varepsilon_0\Omega}\Big)^2\frac{\hbar}{2M\omega_{\rm LO}}\frac{Z^2}{q^2\epsilon_\infty^2},\qquad |g_{\rm F}|^2=\frac{e^2\hbar\omega_{\rm LO}}{2\varepsilon_0\Omega q^2}\Big(\frac1{\epsilon_\infty}-\frac1{\epsilon(0)}\Big).$$
  - With $\epsilon(0)$ from LST, $|g^{\mathcal L}|^2/|g_{\rm F}|^2=\omega_p^2/[\epsilon_\infty(\omega_{\rm LO}^2-\omega_{\rm TO}^2)]$. With the model's $\epsilon(0)=\epsilon_\infty+\omega_p^2/\omega_{\rm TO}^2$, the ratio is $(\epsilon_\infty\omega_{\rm TO}^2+\omega_p^2)/(\epsilon_\infty\omega_{\rm LO}^2)$.
  - Either way the ratio is 1 if and only if $\omega_{\rm LO}^2-\omega_{\rm TO}^2=\omega_p^2/\epsilon_\infty$. That is the screened L-mode shift of (410c), which has no $c$ and holds at every $q$. So it holds for (410e) and (410f) alike (script §8: ratio 1 for both, with either $\epsilon(0)$).
  - E.410's unscreened shift passes the test in its own electron-free model ($\epsilon_\infty=1$ in G.180, Fröhlich and $\epsilon(0)$; script §8: ratio 1).
  - Mixing the unscreened shift with $\epsilon_\infty=10.9$ gives a ratio that depends on the choice of $\epsilon(0)$: 0.0917 $=1/\epsilon_\infty$ with LST ($\epsilon(0)=30.4$), and 0.417 with the model's $\epsilon(0)=12.69$. The first draft's "differ by a factor $\epsilon_\infty$" was therefore an artefact of that choice, and is withdrawn.
  - **Conclusion:** G.180 encodes the screened, unretarded L-mode shift. It cannot tell (410e) from (410f), so it says nothing about the order of limits.
- **Huang (1951), Born–Huang.** The transverse roots satisfy $c^2q^2=\omega^2\epsilon(\omega)$ (script §5, relative deviation $<3\times10^{-9}$). The upper branch starts at $\omega_{\rm L}$, the lower one is photon-like, and the reststrahlen gap lies between $\omega_{\rm TO}$ and $\omega_{\rm L}$. LST holds: $\epsilon(0)=12.69$ for the $Z^*=2.07$ numbers.

### 3.4 Numerical script (`test/corr-e410/check_e410.py`, Python 3 + numpy)
1. $D^{\rm e}$ from E.388, solved numerically with the E.375 matrix, against (410a), isotropic and tensor. The worst relative deviation is $10^{-11}$, limited by roundoff at $cq/\omega=10^3$.
2. $\Omega^{\rm e}$ from E.392 and E.399 with random complex $\mathbf j$, against (410b): $4\times10^{-16}$.
3. Dimensions.
4. GaAs numbers (table below).
5. Self-consistent roots of (410c) by bisection, against the closed form (410g). Also Huang, LST, weights and the limits; branch tables for $\epsilon_\infty=10.9$ and for $\epsilon_\infty=1$.
6. (410d) for a degenerate triplet in a random basis, and the failure of the diagonal approximation in that basis.
7. The $c\to\infty$ limit with an anisotropic $\hat\epsilon^\infty$, converging to the non-analytic term as $(\omega/cq)^2$.
8. G.180 against Fröhlich (rewritten after review). Ratio 1 for the LO of (410e) and of (410f), each with LST and with model $\epsilon(0)$; ratio 1 for E.410 in its electron-free model. In the mixed case the ratio depends on the choice of $\epsilon(0)$ (0.0917 vs 0.417).
9. Basis phase (new). LA[111] of GaAs with cell-phase eigenvectors and no $e^{-i\mathbf q\cdot\mathbf r^0_a}$: the L-term of (410c) gives $\sqrt{\Delta\omega^2}/q=2508$ m/s with As at $\frac a4(111)$, and 5852 m/s with As moved by the lattice vector $\frac a2(110)$. Compare $v_{\rm LA}[111]=5393$ m/s from $C_{11},C_{12},C_{44}$. With the phase, the term vanishes (to $10^{-10}$ m/s).
10. Probed wave vectors and the bare frequency (new). Visible backscattering Raman: $q=4\pi n/\lambda=1.03\times10^{-2}$ Å$^{-1}$ = 184 $q_c$ (514.5 nm, $n=4.2$) and $7.7\times10^{-3}$ Å$^{-1}$ = 138 $q_c$ (632.8 nm, $n=3.86$). Electron-free crossover $\omega_{\rm TO}/c=0.30\,q_c$. E.410 with the Lorentz-shifted bare $\omega_{\rm O}$ of §4.2.1 would give 61.3 meV.

---

## 4. Numbers (GaAs; literature values, illustrative)

The inputs are $\hbar\omega_{\rm TO}=33.3$ meV, $\hbar\omega_{\rm LO}=36.2$ meV, $\epsilon_\infty=10.9$, $Z^*=2.07\,e$, $\Omega_{\rm uc}=45.2$ Å$^3$ and $M=36.1$ amu. The constants are $e^2=14.3996$ eV Å (Gaussian), $\hbar c=1973.27$ eV Å and $\hbar^2/({\rm amu\,Å^2})=4.18016$ meV. The bare frequency is $\omega_\sigma(0)=\omega_{\rm TO}$, because the $\mathbf G=0$ Dyson treatment has no Lorentz local-field term (§2.3).

| Quantity | Value |
|---|---|
| $\hbar\omega_p=\hbar\sqrt{4\pi e^2Z^{*2}/(\Omega_{\rm uc}M)}$ | 44.57 meV |
| E.410, $\sqrt{\omega_{\rm TO}^2+\omega_p^2}$ (no electrons, $q\to0$, finite $c$; all polar modes) | **55.63 meV** |
| E.410 with the bare $\omega_{\rm O}$ of §4.2.1, $\omega_{\rm O}^2=\omega_{\rm TO}^2+\omega_p^2/3$ (not the right input for the $\mathbf G=0$ model) | 61.30 meV |
| (410e), $q\to0$ at finite $c$ with $\epsilon_\infty$ (all three polar modes) | 35.93 meV |
| (410f), $cq\gg\sqrt{\epsilon_\infty}\omega$: LO / TO | 35.93 / 33.30 meV (measured 36.2 / 33.3) |
| Screening alone (E.410 → (410e)) | 55.63 → 35.93 meV |
| Order of limits alone ((410e) → (410f)), TO modes | 35.93 → 33.30 meV; LO unchanged |
| $Z^*$ implied by the measured splitting, $\omega_p^2=\epsilon_\infty(\omega_{\rm LO}^2-\omega_{\rm TO}^2)$ | 2.177 $e$, which would put E.410 at 57.50 meV (the primer's B.7 figure uses this) |
| $q_c=\sqrt{\epsilon_\infty}\omega_{\rm TO}/c$ | $5.57\times10^{-5}$ Å$^{-1}$ = 0.557 μm$^{-1}$ |
| $q_c$ at $\omega_{\rm LO}$ | $6.06\times10^{-5}$ Å$^{-1}$ |
| $q_c/(2\pi/a)$, $a=5.653$ Å | $5.0\times10^{-5}$ |
| Vacuum light wave vector $\omega_{\rm TO}/c$ (electron-free crossover) | $1.69\times10^{-5}$ Å$^{-1}$ |
| TO shift at $q=10q_c$ / $100q_c$ | $-27.6$ μeV / $-0.27$ μeV |
| $\epsilon(0)$ from LST (with $Z^*=2.07$) | 12.69 |
| Visible backscattering Raman $q$ (514.5 / 632.8 nm) | $1.03\times10^{-2}$ / $7.7\times10^{-3}$ Å$^{-1}$ = 184 / 138 $q_c$ |
| Spurious LA[111] term without basis phase, $\sqrt{\Delta\omega^2}/q$ | 2.5 km/s (cell-dependent; 5.9 km/s in another cell), vs $v_{\rm LA}=5.4$ km/s |

Polariton branches with $\epsilon_\infty=10.9$ (script §5):

| $q/q_c$ | $q$ (Å$^{-1}$) | $\hbar\omega_-$ (meV) | $\hbar\omega_+$ (meV) | $Z_-$ | $Z_+$ | $\hbar\omega_{\rm L}$ |
|---|---|---|---|---|---|---|
| 0 | 0 | 0 | 35.932 | 0 | 1 | 35.932 |
| 0.1 | 5.57e-6 | 3.084 | 35.954 | 0.0012 | 0.9988 | 35.932 |
| 0.3 | 1.67e-5 | 9.205 | 36.141 | 0.0123 | 0.9877 | 35.932 |
| 1 | 5.57e-5 | 27.228 | 40.727 | 0.4007 | 0.5993 | 35.932 |
| 3 | 1.67e-4 | 32.964 | 100.918 | 0.9775 | 0.0225 | 35.932 |
| 10 | 5.57e-4 | 33.272 | 333.276 | 0.9983 | 0.0017 | 35.932 |
| 100 | 5.57e-3 | 33.300 | 3330.0 | 1.0000 | 0.0000 | 35.932 |

Electron-free ($\epsilon_\infty=1$, §7.1 as printed), with $q_0=\omega_{\rm TO}/c$:

| $q/q_0$ | $\hbar\omega_-$ | $\hbar\omega_+$ | $Z_-$ | $Z_+$ | $\hbar\omega_{\rm L}$ |
|---|---|---|---|---|---|
| 0.001 | 0.020 | 55.635 | 0.0000 | 1.0000 | 55.635 |
| 1 | 17.784 | 62.352 | 0.2219 | 0.7781 | 55.635 |
| 10 | 33.003 | 335.999 | 0.9821 | 0.0179 | 55.635 |
| 100 | 33.297 | 3330.3 | 0.9998 | 0.0002 | 55.635 |

---

## 5. Downstream: what changes and what does not

**Equations.**
- No equation is algebraically wrong, and no later equation uses E.410, so nothing inherits a change.
- **Relabelled:** E.410 → (410′). The content is the same; the left side now states $\mathbf q\to0$ at finite $c$. The tex follows it with `\stepcounter{equation}` so that (411)–(472) keep their numbers (§7).
- **New:** (410a)–(410h), in a new §7.1.1 after E.410.
- **Unchanged:**
  - E.406–E.409 (apart from the subscript slip $\Omega^0_\sigma\to\Omega^0_{\sigma\sigma}$ in E.406, noted with the photon/phonon slips).
  - E.411–E.454: the Fröhlich vertices and the electrostatic self-energy. They do not depend on E.410.
  - E.455–E.459 and the $z$-dependent Dyson equation on p. 100: the general, finite-$c$ transcription to 2D, which is unaffected.
  - E.460–E.472 (§8.1): the electrostatic 2D results. They do not describe $Q\lesssim\omega/c$, but nothing in them changes.
  - E.159–E.165: the classical, electrostatic §4.2.1.

**Prose.**
1. **§7.1 title (p. 90):** "Lowest self-energy correction and the phonon polariton at $\mathbf q=0$".
2. **§7.1 opening (p. 90):** replace the QE sentence with the paragraph in the tex (item 2). It keeps the valid part of the analogy (the long-range field dressing short-range phonons), attributes the non-analytic term to Giannozzi et al. 1991, Baroni et al. 2001 and Gonze–Lee 1997, and says that the limit taken here is different. It also points the reader to §7.1.1 and §7.3.
3. **The two photon/phonon slips (p. 91, primer erratum 48):** "the photon self-energy is" (before E.406) and "the new corrected photon energy satisfy" (before E.407) become "phonon"; $\Omega^0_\sigma\to\Omega^0_{\sigma\sigma}$ in E.406.
4. **Before E.408 and before E.410 (p. 91):** say "$\mathbf q=0$ for finite $c$, $cq\ll\omega$, isotropic" and "the pole condition gives, for $\mathbf q\to0$ at finite $c$ and without electrons".
5. **After E.410:** the interpretive paragraph (now with the L-mode statement) and the new §7.1.1.
6. **§7.3 (p. 99):** new title "… and the LO–TO splitting". One paragraph, inserted after the sentence that follows E.454 (the one ending "Nano Lett. 17, 3758 (2017)."), identifies E.454 as the non-analytic term and as the $cq\gg\sqrt{\epsilon_\infty}\omega$ limit of (410c).

---

## 6. What the primer should change

The primer already gives the correct physics in B.7, B.9, X.1, X.2, Z.1 and Z.2 item 47. The changes align it with the correction and close the remaining gaps.

1. **`sections/B07-phonon-se.html`** (§B.7):
   - **Title.** Rename the `h2`, `data-title` and lede ("The electrodynamic phonon self-energy and the LO–TO splitting") to something like "The electrodynamic phonon self-energy, polaritons and the order of limits", so the heading does not repeat the manuscript's misreading. Update the section map in `build/BUILDER_GUIDE.md` (B.7 line) to match.
   - **Main display.** In "What E Eq. 410 is, and what it is not", promote the screened general-$q$ pole condition (correction 410c, with $\epsilon_\infty$) to the main display. It currently shows only the electron-free version $\omega_\sigma^2+\omega_{p\parallel}^2+\omega_{p\perp}^2\omega^2/(\omega^2-c^2q^2)$, and the $\epsilon_\infty$ form appears only in the derivation box. Do not tag it with an E number (BUILDER_GUIDE forbids invented source numbers). Write "derived; proposed as Eq. (410c) in the correction, §Z.2 item 47".
   - **Bullets.** Add two bullets to the existing three:
     - "Even with electrons, the $q\to0$, finite-$c$ limit is isotropic. All three polar modes sit at $\omega_{\rm L}$ (35.9 meV). Screening fixes the size of the shift (55.6 → 35.9 meV); the order of limits decides whether the TO modes share it."
     - "For the longitudinal mode, E.410 is the electron-free LO frequency at every $q$ and equals E.454 with $\epsilon_\infty=1$. Only the transverse modes depend on the order of limits."
   - **Derivation box.** Add the compact tensor form $D^{\rm e}=-4\pi c[\omega^2\hat\epsilon^\infty-c^2(q^2-\mathbf q\mathbf q)]^{-1}$ and the matrix pole condition (410d). Do not call the Maxwell matrix $M$, which is the reduced mass. Add that the diagonal approximation is exact for the cubic diatomic only in an L/T-adapted basis; in a random basis it misplaces the upper roots by up to 0.75 meV at $q=3q_c$ (script §6).
   - **Order of limits.** Next to the existing $q_c=5.6\times10^{-5}$ Å$^{-1}$, add that the TO shift at $q=10q_c$ is only $-28$ μeV. Say that the polariton region is probed by infrared spectroscopy and near-forward Raman scattering, while visible backscattering Raman ($q\approx10^{-2}$ Å$^{-1}$, more than $10^2q_c$) sees the electrostatic TO and LO.
   - **Optional.** A small static table of $\omega_\pm$, $Z_\pm$ at $q/q_c=0.1, 0.3, 1, 3, 10$ (§4 above) as a figure caption or readout. The interactive figure already uses the right model (constant $\epsilon_\infty$ in $\Pi^{\rm e,0}$) and needs no change.
   - **New connect callout (corrected after review).** "G.180 equals the Fröhlich coupling if and only if $\omega_{\rm LO}^2-\omega_{\rm TO}^2=\omega_p^2/\epsilon_\infty$, the screened, unretarded L-mode shift, which is the same at every $q$. The unscreened shift of E.410 is consistent with it only for $\epsilon_\infty=1$. The check tests screening, not the order of limits." Do not write "holds only with the E.454 splitting" or "differ by $\epsilon_\infty$".
2. **`sections/B09-electrostatic.html`** (§B.9):
   - After the GaAs numbers (35.9 meV, already consistent), add one sentence: E.454 is the $cq\gg\sqrt{\epsilon_\infty}\omega$ limit of the general pole condition of §B.7, which the correction proposes to state in §7.3. Retardation is negligible except for $q\lesssim q_c$.
   - Mention that the correction retitles EPP §7.3 "… and the LO–TO splitting".
3. **`sections/X01-light-cone.html`** (§X.1): "$q_c\approx6\times10^{-5}$ Å$^{-1}$" should name the frequency used. Write $5.6\times10^{-5}$ at $\omega_{\rm TO}$ or $6.1\times10^{-5}$ at $\omega_{\rm LO}$, consistent with B.7 (which uses $\omega_\sigma=\omega_{\rm TO}$). Without electrons the crossover is $\omega/c$, not $q_c$.
4. **`sections/Z02-errata.html`**, item 47: expand "Suggested correction" to the four actions of the tex:
   1. retitle §7.1;
   2. replace the QE sentence;
   3. add (410a)–(410h) after E.410: the general-$q$ pole condition, its two limits and $q_c$;
   4. identify E.454 as the non-analytic term in §7.3.

   Mention that E.410 itself is algebraically correct, and that its 56 vs 36.2 meV mismatch is the electron-free approximation, not the order of limits. Keep the existing sentence "The two agree only for a mode along $\hat{\mathbf q}$ with $\epsilon^\infty=1$", which is right. Cross-reference item 48 (the same passage). The intro paragraph (line 9) is already right.
5. **`sections/X02-dictionary.html`** ("Phonon polaritons" row) and **`sections/Z01-outlook.html`** ("Polaritons beyond the zone centre"): both say the manuscript works out only the $q=0$ end point. Add "(the proposed correction, §Z.2 item 47, adds the full pole condition)". The LO–TO row of X.2 is already right.
6. **Notes.**
   - `notes/E9.md` §3.3 "General-$\mathbf q$ form": add the $\epsilon_\infty$ version and the matrix form, plus the basis-phase caveat for acoustic modes at finite $q$.
   - `notes/SYNTHESIS.md` §6 item 7 and §15 "Key equations": add the screened general-$q$ formula and a pointer to `corrections/issue-e410.tex`.

---

## 7. Caveats and open points

- **The electron model is the insulating, adiabatic, $\mathbf G=0$ limit.** It uses a constant $\epsilon_\infty$, the head of the RPA dielectric matrix without local fields or xc, consistent with E.380 and E.432–E.446.
  - Spatial dispersion of $\Pi^{\rm e}$ is negligible at $q\lesssim q_c$.
  - With free carriers, $D^{\rm e}$ needs its Drude/plasmon pole, which gives LO-phonon–plasmon modes (the dynamical E.453). This is outside this issue.
- **Basis phases.** The missing basis phase $e^{-i\mathbf q\cdot\mathbf r^0_a}$ (erratum Z.2 #14; companion correction `issue-phase.tex`, Eq. (385′)) does not matter for optical modes at $q\lesssim q_c$. The general-$q$ formulas (410c), (410d) and (410f) need it at finite $q$. Without it, acoustic modes acquire a spurious, cell-dependent dipolar term comparable to the sound velocity (§3.4 item 9). The tex says so in one clause.
- **Damping and real crystals.** Damping is not included. Real GaAs values (E.410 = 55.6 vs 57.5 meV) depend on which $Z^*$ is used: 2.07 from the literature, or 2.18 implied by the measured splitting with this model's $\epsilon_\infty$. They also depend on the bare frequency: $\omega_{\rm TO}$ in the $\mathbf G=0$ model; the bare $\omega_{\rm O}$ of §4.2.1 would give 61.3 meV.
- **2D.** The 2D analogue (§8) at $Q\lesssim\omega/c$ involves the sheet propagator with a radiative part. It is not derived here. E.455–E.459 are general, and §8.1's electrostatic results need no change.
- **Labelling and the equation counter.** amsmath's `equation` environment steps the counter even when `\tag` is given, while `equation*` never does (checked in `amsmath.dtx`, v2.17). The tex therefore writes (410′) as `equation*` + `\tag{410$'$}` followed by `\stepcounter{equation}`, and puts the new (410a)–(410h) in `equation*`. The authors may prefer to leave E.410 unprimed and change only the prose; the relabelling is cosmetic.
- **Citations to check.** K. Huang, Proc. R. Soc. Lond. A 208, 352 (1951); P. Giannozzi, S. de Gironcoli, P. Pavone and S. Baroni, Phys. Rev. B 43, 7231 (1991); S. Baroni, S. de Gironcoli, A. Dal Corso and P. Giannozzi, Rev. Mod. Phys. 73, 515 (2001); X. Gonze and C. Lee, Phys. Rev. B 55, 10355 (1997).

---

## 8. Review responses

Two skeptics sent 20 objections (8 and 12). Several are duplicates (marked). Each was checked against the PDF page images and the derivation. All are accepted in substance. One number in R2-8 is corrected, and R1-7 is accepted with a precision added.

### Reviewer 1

| # | Objection | Verdict | Action |
|---|---|---|---|
| R1-1 (major) | "Differ by a factor $\epsilon_\infty$" depends on the choice of $\epsilon(0)$; the G.180 check tests screening, not the order of limits | **Accepted.** Recomputed: with the model's $\epsilon(0)=12.69$ the mixed ratio is 0.417, not 0.0917. The L mode is unretarded, so (410e) passes the test as well as (410f). | Tex Justification (v) rewritten: "if and only if $\omega_{\rm LO}^2-\omega_{\rm TO}^2=\omega_p^2/\epsilon_\infty$ … for (410e) and (410f) alike … only for $\epsilon_\infty=1$ … tests screening, not the order of limits". Report §3.3 and the primer callout (§6) changed to match. Script §8 rewritten: 5 ratio-1 cases plus a check that the mixed ratio depends on $\epsilon(0)$. The old assertion is removed. |
| R1-2 | $\omega_p$ is defined on p. 33, not p. 32 | **Accepted** (page image: top of p. 33, after E.152). | Tex: "defined after Eq. (152), p. 33". Report §1 and §2.8 fixed. |
| R1-3 | E.455–E.459 are general, finite-$c$ | **Accepted** (p. 100: E.458 uses the full $D^{\rm e}(\mathbf Q,\omega,z,z')$; only §8.1 takes $c\to\infty$). | Tex "What changes" and report §5 rewritten. |
| R1-4 | `equation*` + `\tag` drops a counter step | **Accepted.** `equation*` never steps the counter. (amsmath's `equation` does step it, even with `\tag`.) | `\stepcounter{equation}` after (410′), with a LaTeX comment giving the alternative and saying to delete it if the manuscript numbers by hand. |
| R1-5 | Second photon→phonon slip before E.407 | **Accepted** (page image, p. 91). | Parenthetical in item 2 now lists both slips and the $\Omega^0_\sigma$ subscript (erratum 48). |
| R1-6 | Electron-free crossover is $\omega/c$; "probed by optics" too broad | **Accepted.** Script §10: backscattering Raman probes 184 $q_c$ (514.5 nm) and 138 $q_c$ (632.8 nm). | Tex: "polariton region $q\lesssim q_c$ ($q\lesssim\omega/c$ without electrons), probed by infrared spectroscopy and near-forward (polariton) Raman scattering. Visible-light backscattering Raman probes $q\approx10^{-2}$ Å$^{-1}$, more than $10^2q_c$, and sees the split TO and LO of (410f)." |
| R1-7 | "Holds over essentially the whole BZ" overstated | **Accepted, with a precision.** Within the manuscript's $\mathbf G=0$ treatment, the $\mathbf G=0$ dipole term at finite $q$ is exactly $4\pi|\hat{\mathbf q}\cdot\mathbf j|^2/(\Omega\hat{\mathbf q}\epsilon^\infty\hat{\mathbf q})$, provided $\mathbf j$ carries the basis phase. Local fields add the $\mathbf G\neq0$ terms, Eq. (454a) of the companion correction. | Tex: "Except for $q\lesssim q_c$, retardation is therefore negligible, and the electrostatic result (454) applies within the manuscript's $\mathbf G=0$ treatment (in the diagonal approximation; for mixed polar modes, the $c\to\infty$ limit of (410d))." |
| R1-8 | Report: (405) is cited on pp. 99 and 102, not 101 | **Accepted** (p. 102: "the same equation (405) applies in 2D case"). | Report §1 fixed. |

### Reviewer 2

| # | Objection | Verdict | Action |
|---|---|---|---|
| R2-1 (major) | The Problem paragraph and check (v) blame the missing screening on the order of limits; G.180 cannot distinguish the limits | **Accepted.** In (410c) the L term has no $c$. E.410 is the electron-free LO at every $q$ and equals E.454 for $\epsilon_\infty=1$. The 56 vs 36.2 meV gap is screening (55.6 → 35.9 meV), and the order of limits moves only the TO modes (35.9 → 33.3 meV). | Problem paragraph now separates (a) order of limits, with the L-mode statement, from (b) screening, which is where the 56 vs 36.2 meV gap goes. Item 4 and §7.1.1 add "for the longitudinal mode alone the two agree when $\epsilon_\infty=1$" and "the difference between 55.6 and 35.9 meV is electronic screening". (v) rewritten as for R1-1. Report §0, §2.7, §3.3, §4 and §6 changed to match. |
| R2-2 | $M$ is used for both the Maxwell matrix and the reduced mass | **Accepted.** | Maxwell matrix renamed $\Lambda_{\mu\nu}$ in the tex, at (410a) and (410d), and throughout the report. $\Lambda$ is not used in EPP (checked in the text dump). |
| R2-3 | The paste point after E.454 splits a sentence | **Accepted** (p. 99: "which is very similar (probably identical) …" continues the sentence of E.454). | Tex item 5: insert after "… Nano Lett. 17, 3758 (2017).". |
| R2-4 | = R1-3 | Accepted | As R1-3. |
| R2-5 | = R1-7 | Accepted | As R1-7. |
| R2-6 | Basis phase at finite $q$; spurious acoustic shift (2.5 km/s vs $v_{\rm LA}\approx5.2$ km/s) | **Accepted.** Reproduced: 2508 m/s for LA[111] with As at $\frac a4(111)$, and 5852 m/s in another cell choice. The reference velocity from elastic constants is 5.39 km/s, slightly above the reviewer's 5.2; the conclusion is unchanged. | Tex: one clause citing Eq. (385′) of the companion basis-phase correction. Script §9 added. Report §2.7 and §7 updated. |
| R2-7 | = R1-2 | Accepted | As R1-2. |
| R2-8 | "Probed by optics" loose | **Accepted; one number corrected.** The reviewer's "visible backscattering probes $q\sim10^{-3}$ Å$^{-1}$, about 20 $q_c$" is low by about a factor 8. Backscattering gives $q=4\pi n/\lambda\approx1.0\times10^{-2}$ Å$^{-1}$ (514.5 nm, $n\approx4.2$), about 180 $q_c$, in line with Reviewer 1's figure. The conclusion ($q\gg q_c$, electrostatic TO/LO) holds either way. | Same wording as R1-6. |
| R2-9 | = R1-5 | Accepted | As R1-5. |
| R2-10 | Say which bare frequency enters E.410 | **Accepted.** Reproduced: with §4.2.1's bare $\omega_{\rm O}$, E.410 would be 61.3 meV. | Tex: "With $\omega_\sigma(0)=\omega_{\rm TO}$, since the $\mathbf G=0$ treatment has no Lorentz local-field term, …". The 61.3 meV figure is in report §4 and script §10. |
| R2-11 | "Quantum Espresso (Gonze–Lee) term" misattributes | **Accepted.** | Tex now says "the standard non-analytic term that Quantum Espresso adds", citing Giannozzi et al. 1991, Baroni et al. 2001 and Gonze–Lee 1997 in item 2, and refers to them from item 5. "Gonze–Lee term" removed from (410f) and (iii). |
| R2-12 | = R1-8 | Accepted | As R1-8. |

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
  - (410′) no longer has `\stepcounter{equation}` in the fragment. The merged file states the pasting rule once,
    and every corrected display follows it.
  - In (410b)–(410h), j_ασ(**q**) is stated to be the phase-corrected vertex (385′). The last paragraph of §7.1.1
    refers to (385′) and (131b). Merged §5.3 writes out the two strengths with the phase and states the consequences:
    at q → 0 nothing changes (55.6, 35.9, 33.3 meV); at finite q the strengths are cell-invariant, and a rigid
    acoustic translation gets no dipolar shift (check_merge M4). It also adds that (410d) for cq ≫ √ε ω, in the atom
    basis, is the **G** = 0 term of (454a) (check_merge M5).
  - The §7.3 paragraph is the first of two insertions after E.454; the paragraph with (454a) follows (merged §5.5).
  - The subscript slip Ω⁰_σ → Ω⁰_σσ in E.406 is fixed in place (no tag), per the merged convention.
