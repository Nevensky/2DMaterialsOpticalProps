# Integration decisions for primer v2

These decisions resolve the cross-section conflicts found by the section reviews. They are binding for the integration
pass and extend NOTATION §12.

## 1. One home per symbol

A symbol is defined once, in its **home**: the earliest section in reading order (C00, A.1…A.9, X.1, X.2, B.1…B.10)
whose glossary defines it with that meaning. Each later section:
- reminds the reader of its meaning in a few words at its first use in the section;
- links it with `<a class="sym" href="#gl-…">` to the home entry;
- removes the symbol from its own glossary. Delete the entry if nothing else is in it, and keep the rest of a
  multi-symbol `dt`.

Never delete a home entry, or an entry that another section links to; `gloss.json` lists `linked_from`.

When the home entry's wording is narrower than the later uses, the home section widens it. An example is a light-cone
entry that must become $q_{\mathrm{light}}(\omega)$.

The same symbol with **different meanings** in different sections is allowed only in two cases:
- a model-local variable that is defined in a where-list and whose glossary entry says "in the model of …";
- the cases listed in §3, where one side is renamed.

## 2. Homes (glossary ids other sections link to)

| Symbol | Home id | Notes |
|---|---|---|
| $e$, $m_{\mathrm e}$ | `gl-a2-me` | |
| $Z_\kappa$ (ionic valence) | `gl-a2-Z` | |
| $V_\kappa(\bq)$ (all-space Fourier transform of one ion's potential) | `gl-a2-Vq` | |
| $E_{\mathrm{def}}$ | `gl-a2-Edef` | |
| $n_0$, $k_{\mathrm F}$, $r_{\mathrm s}$ (upright s), $a_0$, $R_{\mathrm{WS}}$ | `gl-a2-jel`, `gl-a2-RWS`, `gl-a2-a0` | write $r_{\mathrm s}$ everywhere, figures included |
| $k_{\mathrm{TF}}$, jellium $\varepsilon^{\mathrm e}(q)$, $\varepsilon^{\mathrm e}(q,\omega)$ | `gl-a2-eps`, `gl-a2-epsw` | |
| $F_{\mathrm L}(x)$, $v_{\mathrm F}$, $V^{\mathrm{scr}}_\kappa$, $n^{\mathrm{ind}}$ | `gl-a2-FL`, `gl-a2-vF`, `gl-a2-scr` | |
| $\nabla$, $\partial_\alpha$ | `gl-a2-nabla` | |
| $\eta$ (positive infinitesimal) | `gl-a3-eta` (A.3 splits it out of `gl-a3-chi0` into its own entry) | A.1's local $\eta$ (shift eigenvalue) stays as it is |
| $\OVSC$ | `gl-a3-OVSC` | |
| $c$ (speed of light) | `gl-a4-c` (A.4 creates this entry, separate from the light-cone entry) | |
| $\bE$, $A_0$ | `gl-a4-E` | B.2's `gl-b2-A` keeps $\bA$ and the temporal gauge |
| $q_{\mathrm{light}}(\omega)=\sqrt{\varepsilon_\infty}\,\omega/c$ | `gl-a4-qc` (A.4 rewrites its entry and text to define $q_{\mathrm{light}}(\omega)$, and its value at $\omega_{\mathrm{LO}}$) | replaces $q_{\mathrm c}$ (A.4, B.7, B.8), $q_{\mathrm{lc}}$ (A.8) and X.1's own entry; B.7's crossover is $q_{\mathrm{light}}(\omega_{\mathrm{TO}})$ |
| $M_{\mathrm c}$ | `gl-a4-Mc` | |
| $m^*$ | `gl-a4-mstar` | |
| $\mathbf p$, $\mathbf p_{\ell\kappa}$ (electric dipole) | `gl-a4-p` | see §3 for momenta |
| $\hat\psi(\br)$ (field operator) | `gl-a5-psi` | |
| $\mathcal T$ (time ordering) | `gl-a5-Tord` | |
| $\omega_{\mathrm p}$, $\omega^{\mathrm{ion}}_{\mathrm p}$ | `gl-a5-wp` | see §4 |
| $g^{\mathrm b}$ (bare, unscreened vertex) | A.6's entry | B.1 renames its $g^0$ to $g^{\mathrm b}$; its footnote says EPP writes $g^0$ |
| $N_{\mathrm F}$ (density of states per spin and per cell at $E_{\mathrm F}$) | A.6's entry (`gl-a6-NF`) | A.3 writes its per-volume, both-spin quantity as $\partial n_0/\partial\mu$, without the name $N_{\mathrm F}$ |
| $q_{\mathrm D}$ | A.2 or A.6, whichever is earlier with this meaning | |
| $\mathcal V_{\bq\nu}$ (screened vertex) | earliest of A.5, A.6 | |
| $E^{\mathrm{qp}}_{n\bk}$, $E_\bk$, $\omega_0$ (Einstein phonon of the models) | earliest entry, A.6 or A.7 | |
| Matsubara frequencies $\omega^{\mathrm F}_p$, $\omega^{\mathrm B}_p$ | `gl-a7-Matsubara` | |
| $\Sigma^{\mathrm{ep}}$ | A.7's entry | |
| $n_{\mathrm c}$ (density of free carriers in a doped semiconductor) | A.9's entry | B.8 renames its carrier density $n_{\mathrm e}$ to $n_{\mathrm c}$. $n_{\mathrm e}$ is the total electron number density (NOTATION §6) |
| $\sigma_{\alpha\beta}$ (conductivity) | A.9's entry | X.1 links it |
| $\be_\lambda$ (light polarization) | A.9's entry | B.4 keeps an entry for the triad $\be_1,\be_2,\be_3(\bq)=\hat\bq$ and links A.9 |
| $\bB$ | `gl-x1-B` | |
| $D^0$, $D^{0\mathrm L}$, $D^{0\mathrm T}$ | `gl-x1-D0` | B.5 and B.6 link; B.5's derivation is linked from that entry with `a.eqref` |
| $D^{\mathrm L}$, $D^{\mathrm T}$, $\varepsilon^{\mathrm L}$, $\varepsilon^{\mathrm T}$, L and T parts | `gl-x1-LT` | B.6 links |
| $\bj^{\mathrm{ext}}$, $\bj^{\mathrm{ind}}$, $\rho^{\mathrm{ind}}$, $n^{\mathrm{ext}}$, $n^{\mathrm{ind}}$ | `gl-x1-ext` | |
| $\vartheta$ (gauge function) | earliest of X.1, B.2 | |
| $g^{\mathrm{dyn}}$ | `gl-x2-gdyn` | B.8 links |
| $\Pi^{\mathrm{ph},0}$ | `gl-b7-Piph0` | B.8 links |
| $\Omega_{\mathrm S}^{\mathrm{SC}}$ | earliest entry (B.1 or B.10) | |

## 3. Renames

| Where | From | To | Why |
|---|---|---|---|
| A.5, B.4 (and any section that writes a nuclear or ionic momentum) | $\hat\bp_{\ell\kappa}$ | $\hat{\mathbf P}_{\ell\kappa}$ | $\mathbf p_{\ell\kappa}$ is the Born dipole of A.4 and B.1. A capital letter, as for $\bR$ and $M$; $\bP$ without atom labels stays the polarization density |
| A.7 | $E_{\mathrm c}$ (cutoff of intermediate states) | $E_{\mathrm{cut}}$ | $E_{\mathrm c}$ is the conduction-band minimum in A.8 and A.9 |
| X.1 | $\mathcal M_{\alpha\gamma}$ (Maxwell operator) | $\mathsf M_{\alpha\gamma}$ | A.9's $\mathcal M^s_{mn\nu}$ is the two-path amplitude |
| B.8 | $\varepsilon^{\mathrm e}_{\alpha\beta}=\delta_{\alpha\beta}-\sum_\gamma D^0_{\alpha\gamma}\Pi^{\mathrm e}_{\gamma\beta}$ | a distinct symbol for this screening matrix, plus a sentence saying that its longitudinal element equals $\varepsilon^{\mathrm{e,L}}$ | the dielectric tensor is $\varepsilon^{\mathrm e}_{\alpha\beta}=\delta_{\alpha\beta}+(4\pi c/\omega^2)\Pi^{\mathrm e}_{\alpha\beta}$ (B.5, B.9); the transverse elements differ at finite $c$ |
| A.4, B.7, B.8, A.8, X.1 | $q_{\mathrm c}$, $q_{\mathrm{lc}}$, X.1's own entry | $q_{\mathrm{light}}(\omega)$ with its argument | §2 |
| A.3, figures of A.3 and A.5 | $r_s$ | $r_{\mathrm s}$ | §2 |
| A.8 | $\langle u_{\kappa\alpha}u_{\kappa\beta}\rangle$ | $\langle\hat u_{\kappa\alpha}\hat u_{\kappa\beta}\rangle$ | as A.3 and (A1.28) |

## 4. Numbers and conventions

- **Plasma frequency.** $\omega_{\mathrm p}^2=4\pi n_{\mathrm c}e^2/m^*$ (unscreened), so that
  $\varepsilon(\omega)=\varepsilon_\infty-\omega_{\mathrm p}^2/\omega^2$ for free carriers. B.9 rewrites its screened
  definition in this convention. Where the screened frequency is needed, write $\omega_{\mathrm p}/\sqrt{\varepsilon_\infty}$.
- **GaAs.** Use one parameter set everywhere:
  - $\hbar\omega_{\mathrm{TO}}=33.3$ meV, $\hbar\omega_{\mathrm{LO}}=36.2$ meV (measured);
  - $\varepsilon_\infty=10.9$, $\varepsilon_{\mathrm s}=12.9$;
  - $Z^*=2.07$, $M_{\mathrm r}=36.11$ amu, $\OV=45.17$ Å³;
  - hence $\hbar\omega^{\mathrm{ion}}_{\mathrm p}=44.6$ meV.

  A.5 changes its 46.9 meV to 44.6 meV and says that $Z^*\approx2.18$ would reproduce the measured splitting exactly,
  as B.3 does.
- **Sodium.** Use the free-electron values of Ashcroft and Mermin: $n=2.65\times10^{22}$ cm⁻³, $r_{\mathrm s}=3.93$,
  $\hbar\omega_{\mathrm p}=6.05$ eV. B.6 changes 5.9 eV to 6.05 eV, and may add that the measured plasmon is about 5.7 eV.
- **Linewidths.** A.6's $\gamma_{\bq\nu}$ entry states its width convention, half or full width, and its relation to the
  damping rate $\gamma$ of the oscillator models of B.2 and B.3.
- **Photon-propagator sectors** (NOTATION §12.4 is updated to match B.6 and B.9):
  - $D_{\rho\rho}$ is normalized so that $D^0_{\rho\rho}=v/e^2$ in vacuum, and $W^{\mathrm e}=e^2D^{\mathrm e}_{\rho\rho}$;
  - the current–current sector is $D^{\mathrm T}_{\alpha\beta}$;
  - $D_{\alpha\beta}$ is the full 3×3 propagator.

  B.5, B.7, B.8, B.9, X.1 and X.2 must agree with this. Wherever a section uses another normalization, convert it.
- **Operator hats.** Part A writes hats on operators where an operator and a number of the same name meet, as in A.1.
  Part B works with operators throughout and drops the hats. Each Part B section that does so says it once, and C00
  states the rule.
- **Comma rule.** A comma may also separate a band index from a wavevector expression, as in $E_{1,k+q}$ and
  $\psi_{m,k+q}$.

## 5. Layout

- Inline math must not be wider than the text column at 390 px. The checker lists offenders as CLIPPED INLINE MATH.
  Shorten the formula or move it to a display equation.
- Display equations: at 390 px, the body plus tag of each display should fit its box. The A.1 pilot already scrolls
  some equations sideways. Split every display whose overflow exceeds 40 px (`eqwidth.json` lists them). Smaller
  overflows may stay.
