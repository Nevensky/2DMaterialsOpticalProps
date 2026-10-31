# The missing factor e in the electrostatic Fröhlich vertices (EPP draft v2, pp. 95 and 98)

Full derivation, checks and numbers behind `corrections/issue-efactor.tex`.

- Paste-ready correction: `doc/epp-primer/corrections/issue-efactor.tex`
- Numerical checks: `doc/epp-primer/test/corr-efactor/check_efactor.py` (run with `python3 -I check_efactor.py`; all assertions pass, re-run after the review round with the new checks §2e longitudinal part, §2f sign convention and the 2k_F range of §4) and its output `results.json`. The `test/` folder is git-ignored.
- Conventions: EPP §4 on is Gaussian, Born charges carry the charge (`Z*_EPP = |e| Z*_G`), and `e` is the signed electron charge (`e = -|e|`; Lagrangian E.134, momenta E.182). Giustino's RMP is SI with dimensionless `Z*`. Numbers use `e^2 = 14.3996 eV Å` (Gaussian, the SI `e^2/4πε0`), `ħ^2/(amu Å^2) = 4.1803 meV`, `ħ^2/2m_e = 3.80998 eV Å^2` and `ħc = 1973.27 eV Å`.
- Every equation quoted below was checked against the page images of the PDFs (EPP pp. 19–26, 35–49, 55–58, 60–102; the PDF has 102 pages; Giustino pp. 34–35). During the review round, EPP pp. 23–25, 55–58, 68–69, 83, 91–99 and 101–102 were read again.

---

## 1. What the draft says

| Where | What is printed | Kind of density |
|---|---|---|
| E.106 (p. 19) | `ρ⁰_{nk,mk+q}(G) = ⟨φ⁰_nk\| e^{-i(q+G)r} \|φ⁰_mk+q⟩` | number (dimensionless matrix element) |
| E.166, p. 38, E.179–E.180 | `j_e = e Σ v_e δ(r - r_e)`; continuity `∇j_e = -dρ_e/dt`; Gauss `∇D = 4πρ_e` | here `ρ_e` is a **charge** density (classical, never defined explicitly) |
| E.234 (p. 47) | `∇E = 4π(ρ_e + ρ_ion)`, "where `ρ_e = Ψ⁺Ψ` is electronic charge density operator" | `Ψ⁺Ψ` is a number density, used as a charge |
| E.235, E.237 last term, E.312 last term | `½∫∫ [ρ_e + ρ_ion][ρ_e' + ρ_ion'] / \|r - r'\|` | needs charge; `ρ_ion = -∇P` is a charge density (Z* has charge) |
| E.237, E.270, E.273, E.312 (A² terms); E.279, E.321 (Π^{e,dia}) | `e²/(2m_e c²) ∫ ρ_e A²`; `Π^{e,dia} = -(e²/mc) ρ^e δ` | number (correct as printed) |
| E.239, E.420 | `ρ_e = Ψ⁺Ψ = Σ φ*_nk φ_mk' c⁺_nk c_mk'` | number |
| E.284, E.384, E.419 | `j^μ_{nk,mk+q} = (eħ/2im_e) ∫_Ω dr e^{-iqr} {φ*_nk ∂_μ φ_mk+q - (∂_μ φ*_nk) φ_mk+q}` | **charge** current (contains e) |
| p. 94 (unnumbered), E.422 | `ρ̇_e = -∇j_e`; `[ρ_e, H⁰_e]_- = -iħ∇j_e` | number on the left, charge on the right |
| E.426 | `∫_Ω dr e^{-iqr} φ*_nk φ_mk' = δ_{k',k+q} ρ_{nk,mk+q}` | number |
| E.428 (p. 95, boxed) | `ρ_{nk,mk+q} = -Σ_μ ħ q_μ j^μ_{nk,mk+q} / (E_nk - E_mk+q)` | number = charge: **inconsistent** |
| p. 95 (unnumbered) | `lim g^{dyn,0}_{nmσ} = (4πi/Ω_uc) Σ_a Σ_μν √(ħ/2Nm_aω_σ) q_μ Z*_{aμν} e_{aνσ}/q² [ρ_{nk,mk+q}]^*` | missing e |
| p. 96 (display + blue note) | VG Eq. 4 as `(4πe²i/Ω_uc) Σ √(…) q_μ Z*_{aμν} e_{aνσ}/q² ⟨φ_mk+q\|e^{iqr}\|φ_nk⟩`; note: *"U odnosu na Giustino-a fali ti taj e². Vjerojatno je u interakciji D ili u Z*. Također njegova je oznaka tog matričnog elemeta g⁰_{mnσ}, a naša g⁰_{nmσ}, međutim to je nebitno."* ("Compared with Giustino you are missing that e². It is probably in the interaction D or in Z*. Also, his notation for this matrix element is g⁰_mnσ and ours g⁰_nmσ, but that is irrelevant.") | with the manuscript's Z* (charge) the display has e³ |
| p. 98 (boxed) | `lim g^{dyn}_{nmσ} = (4πi/Ω_uc) Σ √(…) q_μ Z*_{aμν} e*_{aνσ}/(q·ε^e(q,ω)·q) ρ_{nk,mk+q}` | missing e |
| p. 99 | "almost identical to … Verdi and Giustino … the conjugated complex part of the electron-phonon vertex" | |
| p. 101 (unnumbered) | `Π^{e,0}_ρρ(Q,ω) = (2/S) Σ_K Σ_nm (f_nK - f_mK+Q)/(ħω + E_nK - E_mK+Q + iη) \|ρ_{nK,mK+Q}\|²`, then `ε^e_2D = 1 - v_Q Π^{e,0}_ρρ` (E.470), `v_Q = 2π/Q` | missing e² |
| E.131 (p. 25, §3, SI) | `g^L_{nmα}(k,q) = (4πi/Ω_uc)(e/4πε0) Σ_G Σ_aμν √(ħ/2Nm_aω_α) (q_μ+G_μ) Z*_{aμν} e_{aνα}/((q+G)ε^∞(q+G)) ρ_{nk,mk-q}(-G)` | one explicit e, correct (the text says the second e "is hidden in Z*") |

## 2. Diagnosis

### 2.1 Where the e is lost

The derivation on pp. 93–95 evolves `ρ_e` (E.420, a number density built from `Ψ⁺Ψ`) and `j_e` (E.419, which contains e) with `H⁰_e` (E.421) and then uses the continuity equation `ρ̇_e = -∇j_e`. That equation is only true for a charge density. With a number density it must read

    e ρ̇_e = -∇·j_e .

Everything after that inherits the missing e: E.422–E.425, the display after E.427, and E.428.

### 2.2 Free-electron check

For `φ_k = e^{ik·r}/√Ω`: `ρ_{k,k+q} = 1`, `j_{k,k+q} = (eħ/2m_e)(2k+q)` (from E.384), `E_k - E_{k+q} = -(ħ²/2m_e)(2k·q + q²)`. Hence

    -ħ q·j_{k,k+q} / (E_k - E_{k+q}) = e ,  not 1.

The script reproduces this (`free_electron_ratio = -3.794680 = e` in (eV Å)^½).

### 2.3 Dimensional analysis (CGS-Gaussian, from the script)

| Quantity | Dimension | Required | |
|---|---|---|---|
| E.428 printed, left side `ρ` | 1 | esu | mismatch |
| E.428 right side `ħ q·j/ΔE` | esu | esu | |
| E.428′ left side `eρ` | esu | esu | OK |
| p.95 / p.98 printed `Z* q/(Ω_uc q²) √(ħ/mω) ρ` | g^½ cm^½ s^-1 = statvolt | erg | mismatch (a potential, not an energy) |
| p.95′ / p.98′ with e | erg | erg | OK |
| E.415 `J_el* D⁰ J_ion` | erg | erg | OK (charge is in j) |
| E.409 `4π\|j_ασ\|²/Ω_uc` | s^-2 | s^-2 | OK |
| p.101 printed `v_Q Π^{e,0}_ρρ` | esu^-2 | 1 | mismatch |
| p.101′ `v_Q e² Π^{e,0}_ρρ` | 1 | 1 | OK |
| E.312 Coulomb printed `ρ_e ρ_ion/\|r-r'\|` | statvolt | erg | mismatch |
| E.312′ `eρ_e ρ_ion/\|r-r'\|` | erg | erg | OK |
| E.237 A² term `e²/(mc²) ρ_e A²` (number density) | erg | erg | OK |

### 2.4 Where the charge is *not* lost (E.414–E.416, E.417–E.418, E.429–E.448)

- E.414 `J^μ = √(1/cΩ) (ħω/ΔE) j^μ` carries e through E.384. E.415, E.416 and their insertions E.417, E.418, E.429, E.430, E.447, E.448 are dimensionally correct (energy); the e is lost only when `q·j` is traded for `ρ` with E.428.
- E.431–E.432: `ε^e = δ - D⁰Π^{e,0}_jj` is dimensionless because `Π^{e,0}_jj` (E.383, E.321) carries e².
- E.435: `D^e_11 = v + v⊗Π^{e,0}_ρρ⊗v` with `v = 1/|r-r'|` requires `Π^{e,0}_ρρ` to be the charge–charge correlator. Via the connection formula E.350, `Π_ρρ = -(c/ω²)∂∂'Π_jj`, it automatically is (Π_jj carries e²). E.436–E.446 are then correct. The only inconsistency is the *definition* E.320 (for i = e) if `ρ^e = Ψ⁺Ψ` is inserted literally; the continuity equations E.335, E.340, E.342 and the transformation E.377 (`ρ^i → j^i`, `v → D^{0l}`) also only make sense with a charge density. Hence the one-sentence fix in Sec. 6 (item 6 of the .tex).
- New finding: the **explicit 2D bubble on p. 101** is written with the number vertex `|ρ_{nK,mK+Q}|²` and no e², while E.461, E.462, E.469, E.470 use it with `v_Q = 2π/Q`. It needs `2e²/S` (p101′). Not listed in the errata so far.
- The p. 96 transcription of VG Eq. 4 keeps VG's e² but uses the manuscript's charge-dimensioned Z*, so it has e³.

## 3. Choice of fix

| | Fix A (adopted): ρ stays a number density, write eρ where charge is meant | Fix B: `ρ_e ≡ eΨ⁺Ψ` | Fix C: also redefine the vertex `ρ_{nk,mk+q}` as a charge vertex |
|---|---|---|---|
| E.106, E.426 | unchanged | unchanged | `e⟨φ\|e^{-i(q+G)r}\|φ⟩`; changes E.105–E.107, E.130–E.131 (where §3 uses −e with e > 0) |
| E.239, E.420 | unchanged | `eΨ⁺Ψ`, factor e | factor e |
| E.234 text, E.235, E.237/E.312 Coulomb | `eρ_e` (errata row 28) | unchanged (text: `ρ_e = eΨ⁺Ψ`) | unchanged |
| A² terms of E.237, E.270, E.273, E.312; Π^{e,dia} of E.279, E.321 | unchanged | `e² → e` | `e² → e` |
| p. 94 continuity, E.422 (written with ρ_e) | e on the left | consistent as printed | consistent as printed |
| E.423, E.424, the two displays between E.424 and E.425, E.425 (E.420 substituted explicitly as `φ*φ c⁺c`) | e on the left | e on the left (inherited from E.420) | e on the left (inherited from E.420) |
| display after E.427 (written with ρ_{nk,mk+q}) | e on the left | e on the left (E.426 is still the number vertex) | consistent (ρ is the redefined charge vertex) |
| E.428 | `eρ` on the left | **still needs** `eρ` (E.426 is a pure wave-function integral) | unchanged |
| p.95, p.98 | gain e | gain e | unchanged |
| p.96 display | `e² → \|e\|` | `e² → \|e\|` | `ρ` contains e |
| p.101 bubble | `2e²/S` | `2e²/S` | unchanged |
| §6 (E.320–E.350) | one sentence: ρ^e there means eρ^e | consistent | consistent |

Fix A changes no definition and no result written with `j^μ`, `j_ασ` or `W^e`, keeps `ρ_{nk,mk+q}` the standard dimensionless pair density of RPA/BSE codes, and puts e only where charge conservation and the Coulomb coupling demand it. Fix B does not remove the e from E.428 (the vertex E.426 contains no ρ_e). Because E.420 is substituted explicitly in E.423–E.425 and in the displays between and after them, Fix B needs the same e there as Fix A. Within §7.2 it therefore edits everything Fix A edits except the continuity equation and E.422. On top of that it makes eight edits elsewhere: E.239, E.420, four A² terms (E.237, E.270, E.273, E.312) and two diamagnetic self-energies (E.279, E.321). It saves only the Coulomb terms of E.234–E.237 and E.312 and the §6 sentence. Fix C leaves E.428 and the boxes formally intact but redefines a §3 object that already appears with an explicit e in E.131, and it drags the §3/§4 sign convention of e into the vertex definition; it still needs e in E.423–E.425 and the displays between them, and it changes E.426.

## 4. Corrected derivation (Fix A)

### 4.1 Continuity and the charge–current relation (pp. 93–95)

With `c_nk(t) = e^{iH⁰_e t/ħ} c_nk e^{-iH⁰_e t/ħ}` (the draft omits ħ; needed for `iħρ̇ = [ρ, H]`, which the draft uses):

    e ρ̇_e = -∇·j_e                                     (p94′)
    e [ρ_e(r,t), H⁰_e]_- = -iħ ∇·j_e(r,t)               (422′)

Commutator: `[c⁺_i c_j, c⁺_l c_l] = (δ_jl - δ_il) c⁺_i c_l`, so `[ρ_e, H⁰_e] = Σ (E_mk' - E_nk) φ*_nk φ_mk' c⁺_nk c_mk'` (E.424 left side, now times e). Fourier transform with `∫_Ω dr e^{-iqr}`, integrate the divergence by parts (`∫ e^{-iqr} ∂_μ F = i q_μ ∫ e^{-iqr} F`):

    e Σ ∫ e^{-iqr} φ*_nk φ_mk' (E_mk' - E_nk) c⁺c = ħ Σ_μ q_μ Σ (eħ/2im_e) ∫ e^{-iqr} {…} c⁺c      (425′)

With E.426 and E.427 and comparing coefficients:

    e ρ_{nk,mk+q} (E_mk+q - E_nk) = ħ q·j_{nk,mk+q}
    ⇒  e ρ_{nk,mk+q} = -Σ_μ ħ q_μ j^μ_{nk,mk+q} / (E_nk - E_mk+q)                           (428′)

Sign and factor agree with the free-electron check. The identity holds for any local `H⁰_e = T + V(r)` (`[e^{-iqr}, V] = 0`); with nonlocal pseudopotentials or scissor shifts j must be the velocity-operator matrix element, as in optics codes. The script verifies (428′) for Bloch states of a non-centrosymmetric local potential to `|eρ + ħq·j/ΔE|/|e| ≤ 3×10⁻⁸` (10⁻¹¹ for the two larger q; the bound is set by eigenvector noise divided by small intraband ΔE), and shows that `(-ħq·j/ΔE)/ρ = -3.7946805 (eV Å)^½ = e` with spread 10⁻⁹.

### 4.2 Bare vertex (p. 95)

E.417 with E.375. The longitudinal part of E.375 equals E.374 (= E.442) for every c: `D⁰·q̂ = -(4πc/ω²) q̂` exactly. This projection is linear in c, but its factor c cancels the `1/c` of the E.417 prefactor. The transverse part is `(4π/c)(δ - q̂q̂)/(q² - ω²/c²)`. The prefactor of E.417 times `-4πc/ω²` gives `-4πiħ/(Ω_uc ΔE)`, which depends on neither ω nor c. Hence, for every c,

    g^{dyn,0}_L = (4πi/Ω_uc) Σ_a Σ_βμ √(ħ/2Nm_aω_σ) Z*_{aβμ} e_{aμσ} q_β/q² Σ_α (-ħ q_α [j^α]^*)/(E_nk - E_mk+q)

and with `Σ_α -ħ q_α [j^α]^*/ΔE = [eρ]^* = e ρ^*` (e real):

    lim_{c→∞} g^{dyn,0}_{nmσ}(k,q,ω) = (4πi e/Ω_uc) Σ_a Σ_μν √(ħ/2Nm_aω_σ(q)) q_μ Z*_{aμν} e_{aνσ}(q)/q² [ρ_{nk,mk+q}]^*      (p95′)

The transverse remainder is of relative order `(ω/cq)²` (script: 5.0×10⁻⁷, 5.0×10⁻⁵, 5.1×10⁻³ for `ω/cq = 10⁻³, 10⁻², 10⁻¹`; for GaAs at q = 0.29 Å⁻¹ it is 4×10⁻⁹). The longitudinal part alone equals (p95′) to 1.3×10⁻¹⁴ at every c tested (script §2e, `q̂·D⁰·q̂/(4πc/ω²) = −1` to 10⁻¹³).

### 4.3 Screened vertex (p. 98)

E.429 → E.430 (with E.295) → E.446 inserted → E.447 → E.448 (all correct). With (428′), `Σ_β -ħ q_β j^β/ΔE = e ρ`:

    lim_{c→∞} g^{dyn}_{nmσ}(k,q,ω) = (4πi e/Ω_uc) Σ_a Σ_μν √(ħ/2Nm_aω_σ(q)) q_μ Z*_{aμν} e*_{aνσ}(q)/(q·ε^e(q,ω)·q) ρ_{nk,mk+q}      (p98′)

For isotropic ε^e the longitudinal part of D^e is `-(4πc/ω²) q̂q̂/ε_L` at every c, so (p98′) is again the longitudinal part of E.430 for every c.

### 4.4 Independent route: the Coulomb term of (235′)/(312′)

`H_{e-ion} = e ∫ dr ρ_e(r) φ_ion(r)`, `φ_ion(r) = ∫ dr' ρ_ion(r')/|r - r'|`, `ρ_ion = -∇·P`, `P = Σ_i p_i δ(r - R⁰_i)`, `p_{iμ} = Σ Z*_{aμν} Q_{iaν}` (E.241), Q from E.242. Fourier: `φ_ion(r) = Σ_i ∫ dq/(2π)³ (-4πi q·p_i/q²) e^{iq(r - R_i)}` (E.121 with ε = 1, 1/ε0 → 4π). Using `Σ_i e^{i(q'-q)R_i} = Nδ` (G = 0) and `∫dq/(2π)³ → Ω⁻¹Σ_q`:

    H_{e-ion} = Σ_{kq nm σ} [ -(4πi e/Ω_uc) Σ_a Σ_μν √(ħ/2Nm_aω_σ) q_μ Z*_{aμν} e_{aνσ}(q)/q² ⟨φ_mk+q|e^{iqr}|φ_nk⟩ ] c⁺_{mk+q} c_{nk} (a_qσ + a⁺_{-qσ})

The bracket is E.131 at ε^∞ = 1, G = 0 in Gaussian units (with §3's e > 0 written as `-e` in §4's signed convention, and k, n, m relabelled), is exactly `g⁰_mnσ` of (p96′) below, and is `-(p95′)`. The e of `eρ_e` in the Coulomb term is the e that (p95′) gains; without it H_{e-ion} would not be an energy.

### 4.5 2D bubble (p. 101)

    Π^{e,0}_ρρ(Q,ω) = (2e²/S) Σ_{K∈SBZ} Σ_nm (f_nK - f_mK+Q)/(ħω + E_nK - E_mK+Q + iη) |ρ_{nK,mK+Q}|²     (p101′)

Check: for a 2D electron gas, `Π^{e,0}_ρρ(Q, 0) = -e² m*/(πħ²)` for `Q ≤ 2k_F` (T → 0), and E.470 gives `ε_2D = 1 + 2m*e²/(ħ²Q) = 1 + q_TF/Q` (Stern 1967; vacuum background, since E.470 has none). Script, numerical K sum (n = 10¹² cm⁻², m* = 0.067, T = 20 K, Q = 0.004 Å⁻¹): `k_F = 0.0251 Å⁻¹`, so `Q < 2k_F = 0.0501 Å⁻¹`, and `k_BT = 1.72 meV ≪ E_F = 35.7 meV` (both asserted). Bubble −2.79880×10⁻³ eV⁻¹Å⁻² (exact −2.79880×10⁻³), `ε_2D = 64.306` vs Stern 64.306 (`q_TF = 0.2532 Å⁻¹`). The printed form gives `1 + 4.396 (eV Å)⁻¹`, not a number.

## 5. Reduction to Verdi–Giustino (Giustino Eq. 180), units and normalization

Giustino Eq. 180 (p. 015003-35, SI):

    g^L_mnν(k,q) = i (4π/Ω)(e²/4πε0) Σ_κ (ħ/2N_p M_κ ω_qν)^½ Σ_{G≠-q} (q+G)·Z*_κ·e_κν(q) / ((q+G)·ε^∞·(q+G)) ⟨ψ_mk+q| e^{i(q+G)·(r-τ_κ)} |ψ_nk⟩_sc

Dictionary to EPP (all items checked on the PDFs):

| Giustino / VG | EPP | Note |
|---|---|---|
| `e²/4πε0` (SI) | `e²` (Gaussian) | p. 96 states this |
| `Z*_κ` dimensionless | `Z*_{aμν} = \|e\| Z*_κ` | E.110 "dimension of charge" |
| Ω = unit-cell volume | Ω_uc | |
| N_p (BvK supercell), ψ normalized in the supercell (`⟨⟩_sc`) | N, φ normalized in `Ω = NΩ_uc` (E.426 `∫_Ω`, E.291) | `[ρ_{nk,mk+q}]^* = ⟨φ_mk+q\|e^{iqr}\|φ_nk⟩` is dimensionless and → δ_nm as q → 0 |
| N_p⁻½ inside the square root | N⁻½ inside the square root (E.417 "Ω = NΩ_uc") | both correspond to `H_ep = Σ g c⁺c(a + a⁺)` without a separate N⁻½ (EPP E.108; VG convention). In Giustino's own Eq. 1/38 convention N_p must be dropped from Eq. 180 (errata row 1) |
| `M_κ`, `e_κν(q)` | `m_a`, `e_{aνσ}(q)` | |
| `τ_κ` | 0 (dipoles at `R⁰_i`, E.314) | Giustino footnote 7: VG Eq. 4 also lacks the τ phase |
| G sum | G = 0 | p. 95 "we have neglected to write G = 0" |

Written in EPP notation (G = 0, τ = 0, ε^∞ = 1, Gaussian):

    g⁰_mnσ(k,q) = (4πi|e|/Ω_uc) Σ_a Σ_μν √(ħ/2Nm_aω_σ(q)) q_μ Z*_{aμν} e_{aνσ}(q)/q² ⟨φ_mk+q|e^{iq·r}|φ_nk⟩          (p96′)

Units: `[|e| Z*] = esu² = erg·cm`; `[1/Ω_uc] = cm⁻³`; `[q/q²] = cm`; `[√(ħ/mω)] = cm`; so `[g] = erg`. In SI: `[e²/4πε0] = J·m`, the same counting gives J.

Results:

    lim g^{dyn,0}_nmσ = (e/|e|) g⁰_mnσ = -g⁰_mnσ                                              (428a)
    lim g^{dyn}_nmσ   = -(e/|e|) q²/(q·ε^e(q,ω)·q) [g⁰_mnσ]^*                                  (448a)

So `|lim g^{dyn,0}| = |g^VG(ε^∞ = 1)|` and `|lim g^{dyn}| = |g^VG|` with `ε^∞ → ε^e(q,ω)`; for real ε^e, (448a) is exactly `[g^VG(ε^∞ → ε^e)]^*` (EPP's "conjugated complex part"). The overall sign `e/|e|` of (428a) is a convention of the diagrammatic vertices E.415–E.416 (the factor i and the order of the vertices). It **cannot** be absorbed in `e_{aνσ}(q)`: (p95′) and (p96′) are linear in the same polarization vector, so `e_a → z e_a` (|z| = 1, including z = −1) leaves the ratio at −1. It also cannot be removed from (428a) and (448a) at once, because `p98′ = −[p95′]^* q²/(q·ε^e·q)` is fixed by the factor i in E.415/E.416. Only the product `g^{dyn} χ⁰ g^{dyn,0}` enters E.413, and it is fixed: it is invariant under `e_a → z e_a`, and with Ω⁰ it reproduces E.452. So the sign cancels in `ΔΩ^e` and `|g|²`. It would matter only if (p95′) were added to vertices in the Hamiltonian convention (e.g. a short-range `g^S` as in VG's `g = g^S + g^L`); there `g⁰_mnσ = −lim g^{dyn,0}_nmσ` must be used. Script §2f checks all of this for z = 1, −1 and e^{0.7i}. Together, `lim g^{dyn} = -[lim g^{dyn,0}]^* q²/(q·ε^e·q)`, the relation already derived in primer §B.8.

Limits: c → ∞ is built in (and the longitudinal part holds at every c, §4.2). For `ħω ≪ E_gap` and q → 0, `ε^e(q,ω)` → the clamped-ion electronic tensor (in EPP the head of the RPA dielectric matrix, no local fields, no xc), and (448a) coincides with VG. Intraband (`ρ → δ_nm`) gives the 1/q Fröhlich divergence for LO modes only (`q·Σ Z* e_TO = 0`).

Origin of VG's e²: one e is in `Z*_EPP`, the other is the e of (428′). D is innocent, which answers the margin note.

## 6. Consistency with the manuscript's own results

1. **E.413 = Σ g^{dyn} χ⁰ g^{dyn,0}.** E.413 is written with currents (charge-correct). Evaluated at c → ∞ with E.446, it equals `Σ g^{dyn}χ⁰g^{dyn,0}` with (p95′)/(p98′) to 10⁻¹³ relative in the model; with the printed vertices the ratio is `14.3996 = e²` in eV Å, i.e. not a number.
2. **Longitudinal ε from currents equals the charge bubble.** From E.432 and `D⁰·q̂ = -(4πc/ω²)q̂`: `q̂·ε^e·q̂ = 1 + (4πc/ω²) q̂·Π^{e,0}_jj·q̂`. With E.383 and (428′):

        q̂·ε^e·q̂ = 1 + (8πe²/Ωq²) Σ (ΔE/ħω)(f - f)|ρ|²/(ħω + ΔE + iη) = 1 - (4π/q²) e² Π_num(q,ω)

   where `Π_num = (2/Ω) Σ (f - f)|ρ|²/(ħω + ΔE + iη)` (3D analogue of p. 101) and the sum rule `Σ_{knm} (f_nk - f_mk+q)|ρ_{nk,mk+q}|² = 0` (completeness) was used. The e² of (p101′) is exactly this one. Script: 1.374487 (currents) vs 1.374487 (charge bubble) vs 1.026 (printed, dimensionful).
3. **Ω⁰ + ΔΩ^e = E.452 (screened Gonze–Lee).** With (p95′)/(p98′), `A ≡ Σ_a √(ħ/2m_aω_σ) q·Z*_a·e_aσ`, `χ⁰ = 2(f-f)[1/ħω - (1 + iη/ħω)/(ħω + ΔE + iη)]` and the sum rule:

        ΔΩ^e = (4π|A|²/Ω_uc q²) · 4πe²Π_num/(q·ε^e·q),    Ω⁰ = 4π|A|²/(Ω_uc q²)
        Ω⁰ + ΔΩ^e = 4π|A|²/(Ω_uc q·ε^e·q) = ħW^e|q·j_σ|²/(2Ω_uc ω_σ)        (E.452, W^e = 4π/(q·ε^e·q))

   For q → 0, ω → 0 this is the non-analytic term `4π|q̂·Z*·e/√m|²/(Ω_uc q̂·ε^∞·q̂)` of Gonze–Lee. Script: 11.455390 meV both ways (|q| = 0.37 Å⁻¹), 6.911846 (0.29 Å⁻¹), 13.396988 (0.065 Å⁻¹). With the printed vertices ΔΩ is smaller by e² and the identity fails.
4. **E.131.** The corrected vertices are E.131 (§3 electrostatics, one explicit e) at ε^∞ = 1 or ε^∞ → ε^e, G = 0, in Gaussian units, as the text before E.417 intends. This holds up to the overall sign e/|e| (§3 has e > 0 and writes the electron charge as −e in E.123; from §4 on, e is signed) and the relabelling n ↔ m, k → k + q (E.131 has `ρ_{nk,mk−q}(−G)`). The printed vertices do not reduce to E.131 even with these allowances.
5. **E.454.** Unaffected (written with `j_ασ` and W^e). GaAs: `ħω_LO = √(33.3² + 13.50²) = 35.93 meV` with Z* = 2.07|e| (measured 36.2).

## 7. Numerical checks (script output)

### 7.1 Model (3D plane waves, a = 4 Å simple cubic, 125 PWs, local potential with an odd (111) component so that ρ and j are complex; 216 k-points; band 1 occupied, gap ≥ 9.09 eV; ionic vertex GaAs-like, Z* = ±2.07|e|, general polarization; ħω = 36.2 meV, ħω_σ = 33.3 meV)

| q | max \|eρ + ħq·j/ΔE\|/\|e\| | ε_L currents | ε_L = 1 − v e²Π | 1 − vΠ (printed) | ΔΩ E.413 [meV] | ΔΩ p95′/p98′ [meV] | E.413 / printed | Ω⁰+ΔΩ [meV] | E.452 [meV] |
|---|---|---|---|---|---|---|---|---|---|
| (1,1,0)b/6, 0.370 Å⁻¹ | 2.6e-11 | 1.374487 | 1.374487 | 1.0260 | −4.289894 | −4.289894 | 14.3996 | 11.455390 | 11.455390 |
| (1,2,0)b/12, 0.293 Å⁻¹ | 4.4e-11 | 1.380288 | 1.380288 | 1.0264 | −2.628492 | −2.628492 | 14.3996 | 6.911846 | 6.911846 |
| (1,0,0)b/24, 0.065 Å⁻¹ | 2.6e-8 | 1.390870 | 1.390870 | 1.0271 | −5.236484 | −5.236484 | 14.3996 | 13.396988 | 13.396988 |

`max |g^{dyn,0}(E.418, j) − g(p95′)| / max |g| ≤ 1e-13` for all q. The model is a wide-gap toy insulator, not GaAs; it tests algebra, prefactors, conjugations and the sum rule.

### 7.2 Finite c (E.417 with the full E.375 vs p95′)

| ω/(cq) | rel. deviation | x²/(1−x²) | longitudinal part vs p95′ |
|---|---|---|---|
| 10⁻³ | 5.0e-7 | 1.0e-6 | 1.3e-14 |
| 10⁻² | 5.0e-5 | 1.0e-4 | 1.3e-14 |
| 0.1 | 5.1e-3 | 1.0e-2 | 1.3e-14 |
| 0.3 | 5.0e-2 | 9.9e-2 | 1.3e-14 |
| 0.6 | 0.28 | 0.56 | 1.3e-14 |

The deviation is the transverse part, ∝ (ω/cq)²; the longitudinal part equals (p95′) at every c (`q̂·D⁰·q̂ = −4πc/ω²` for every c; its c cancels the 1/c of E.417).

### 7.3 GaAs (literature values, illustrative: ħω_TO = 33.3 meV, ħω_LO = 36.2 meV, ε∞ = 10.9, Z* = 2.07|e|, Ω_uc = 45.2 Å³, μ = 36.1 amu; LO mode, G = 0, ρ = 1, ω_LO under the square root in all columns)

| q [Å⁻¹] | √N\|g\| p98′ (ε^e→ε∞) [meV] | Giustino Eq. 180 in SI [meV] | Fröhlich (LST) [meV] | p95′ unscreened [meV] | p98 as printed |
|---|---|---|---|---|---|
| 0.01 | 3040.51 | 3040.47 | 3197.57 | 33141.6 | 0.8013 (eV/Å)^½ |
| 0.05 | 608.10 | 608.09 | 639.51 | 6628.3 | 0.1603 (eV/Å)^½ |
| 0.10 | 304.05 | 304.05 | 319.76 | 3314.2 | 0.0801 (eV/Å)^½ |

- G.180 was evaluated with SI constants (C, ε0, ħ, amu, m); it equals the Gaussian (p98′) to 1.4×10⁻⁵ (rounding of the constants). This is the unit conversion check.
- Fröhlich (Giustino Eq. 7 form): `√N|g_F| = [4πe²/(Ω_uc q²) · ħω_LO/2 · (1/ε∞ − 1/ε0)]^½`, ε0 = 12.881 from LST. Ratio p98′/Fröhlich = 0.9509 = Z*/Z*_LST with `Z*_LST = [Ω_uc μ ε∞ (ω_LO² − ω_TO²)/(4πe²)]^½ = 2.177|e|`: the two agree identically when `ω_LO² − ω_TO² = 4πZ*²/(Ω_uc μ ε∞)` (E.454), and the 5 % gap is a property of the input set.
- Polaron constant (m* = 0.067): α = 0.0708 from LST, 0.0640 from the corrected Born-charge vertex (ratio (Z*/Z*_LST)²).
- The printed p98 evaluates to a number in (eV/Å)^½, which cannot be compared with an energy in any unit system.
- Aside (not part of this issue): with the bare ω_σ = ω_TO under the square root, as EPP writes it, p98′ is larger by √(36.2/33.3) = 1.043; see the quasiparticle-weight caveat in primer §B.8/§X.2.

### 7.4 2DEG (p101′) — see §4.5.

### 7.5 Sign convention of (428a)/(448a) (script §2f; model, q = (1,2,0)b/12, one k, n = 0, m = 3)

| rephasing `e_a → z e_a` | p95′ / g⁰_VG | p98′ / ([g⁰_VG]^*/ε_L) | g^{dyn} g^{dyn,0} [eV²] |
|---|---|---|---|
| z = 1 | −1 | +1 | −1.512308e-05 |
| z = −1 | −1 | +1 | −1.512308e-05 |
| z = e^{0.7i} | −1 | +1 | −1.512308e-05 |

The relative signs in (428a) and (448a) do not depend on the phase of the polarization vector, and the product that enters E.413 is invariant.

## 8. Downstream equations

**Change (factor e or e²):**
- p. 94 continuity display → (p94′); E.422 → (422′); left sides of E.423, E.424, the two unnumbered displays between E.424 and E.425, E.425, and the display after E.427 (×e); E.428 → (428′).
- p. 95 vertex → (p95′).
- p. 96 display → (p96′) (e² → |e| with charge-Z*), new (428a); margin note deleted; replacement paragraph.
- p. 98 box → (p98′); p. 99 sentence replaced, new (448a).
- p. 101 bubble → (p101′).
- Same convention (no displayed result changes): E.234 → (234′) and the sentence after it; E.235 → (235′); last terms of E.237 and E.312 (`ρ_e → eρ_e`); one sentence in §6 after E.315 (ρ^e in E.320–E.350 means eρ^e). These coincide with primer errata row 28.
- Prose in the same passage: p. 93 Heisenberg evolution (ħ); p. 94 "So the equation for charge density is" → "So the equation for the electron (number) density is"; "retarded limit" → "non-retarded (instantaneous) limit" before E.418; the claim before E.417 ("reduce to standard bare Fröhlich vertex Eq. (131)") → "reduces, up to the overall sign e/|e| [see (428a)] and the relabelling n ↔ m, k → k + q, to the unscreened (ε^∞ = 1), G = 0 form of (131) in Gaussian units".

**Do not change:** E.105–E.131 (E.131 already has its e); E.166, E.179–E.180 (read ρ_e as the charge density; one bracketed remark); E.239–E.240; E.321; E.383–E.421 (includes E.413–E.418: written with j); E.426–E.427; E.429–E.454 (E.447–E.448 written with j; E.449–E.454 with j_ασ and W^e); E.455–E.472 except the p. 101 bubble (E.461, E.462, E.469–E.472 are right once the bubble has e²); the diamagnetic A² terms of E.237, E.270, E.273 and E.312 and the diamagnetic self-energies Π^{e,dia} of E.279 and E.321.

## 9. Replacement prose

All replacement text is in `issue-efactor.tex`, items 1–6:
1. the p. 94 "charge density" sentence, continuity sentence and E.422–E.428 (plus the three small fixes in the same passage, the one before E.417 now with the sign/relabelling caveat);
2. "After using (428′) in (418) we get" + (p95′);
3. the p. 96 paragraph, (p96′), (428a), the status of the sign (convention of E.415–E.416, not absorbable in e_aνσ, irrelevant for E.413; Hamiltonian sign when combined with g^S), explanation of VG's e², deletion of the margin note;
4. "and using (428′) we get" + (p98′); the p. 99 sentence with (448a) and the precise statement of agreement (ε^e real, ħω ≪ gap, q → 0, no LFE);
5. p. 101 bubble with e² and the Stern check (Q ≤ 2k_F, T → 0, vacuum background);
6. the Hamiltonian convention (234′), (235′), E.237/E.312 last terms, one sentence in §6.

## 10. What the primer should change

1. **§B.8 `sections/B08-frohlich.html`**
   - "From currents to charges": state the corrected continuity `eρ̇_e = −∇·j_e` (E p94′ / E.422′) instead of quoting the printed `ρ̇_e = −∇·j_e (E.422)`, and tag the identity `E.428′` instead of `cf. E.428`.
   - Replace the two displays "as printed" (p. 95 and p. 98) by the corrected `(E p.95′)` and `(E p.98′)` with the factor e, and add the two comparison lines `lim g^{dyn,0} = (e/|e|) g⁰_VG`, `lim g^{dyn} = −(e/|e|) q²/(q·ε^e·q) [g⁰_VG]^*` (E.428a, E.448a). The sentence "The bare vertex is Verdi–Giustino with ε∞ = 1 and G = 0 …" can then cite them, adding "up to the overall sign e/|e| and the relabelling n ↔ m, k → k + q". Add one sentence saying that the sign is a convention of the factorization E.415–E.416 (the factor i). It cannot be absorbed in the polarization vector, and only `g^{dyn}χ⁰g^{dyn,0}` (E.413) is observable, which is fixed. (The primer does not currently claim that the sign is absorbable; it simply omits the sign.)
   - Source callout "the missing e²": say the correction restores e at its origin (continuity on p. 94), add the dimension argument (printed vertices have the dimension of a potential, erg/esu), note that the p. 96 transcription has e³ with charge-dimensioned Z*, that the p. 101 bubble needs e² (link §B.10), and that the margin note can be removed. The "three smaller slips" are now part of the correction.
   - Optional: one line of numbers, GaAs √N|g| = 304 meV at q = 0.1 Å⁻¹ from both (p98′) and G Eq. 180 (3314 meV unscreened), and "E.413 = Σ g χ⁰ g only with the e".
2. **§B.10 `sections/B10-2d.html`**: the p. 101 display of `Π^{e,0}_ρρ(Q,ω)` must carry `2e²/S` (tag `E p.101′`), with a short source note: without e², E.470 is not dimensionless; with it a 2DEG gives Stern's `1 + 2m*e²/(ħ²Q)` (for Q ≤ 2k_F, T → 0, no background dielectric constant). The figures use `r_0 Q` and are unaffected.
3. **§B.9 `sections/B09-electrostatic.html`**: where `ε_L = 1 − (4π/q²)Π^{e,0}_ρρ` is written (derivation of E.445), say that `Π^{e,0}_ρρ` is the charge–charge bubble, e² times the number bubble built from `ρ_{nk,mk+q}`, as E.350 guarantees. Optionally add the derivation `Ω⁰ + ΔΩ^e = E.452` with the corrected vertices (§6.3 above) as a `details` box: it shows the corrected vertices reproduce the screened Gonze–Lee term.
4. **§B.4 `sections/B04-fields.html` and §B.6 `sections/B06-propagator.html`**: the displays tagged E.234, E.235, E.237, E.312 reproduce the printed number-density form while the page convention redefines `ρ_e = eΨ⁺Ψ` in Coulomb terms. For consistency with the authors' correction, show the primed forms with explicit `eρ_e` (tags E.234′, E.235′, E.237′, E.312′) and keep `Ψ⁺Ψ` in the A² term; add in §B.6 that `ρ^e` in E.320–E.350 means `eρ^e`.
5. **§X.2 `sections/X02-dictionary.html`**: rows on E.428 and on the vertex should cite E.428′, E p.95′, E p.98′; the figure legend "EPP p. 98 as printed (bare ω_σ)" and the caption sentence "EPP's p. 98 vertex as printed" should read "p. 98′ (e restored)", because the code (`pref = 4π/Ω · E2 · Z/ε`) already includes the e.
6. **§Z.2 `sections/Z02-errata.html`**: row 54 — cite `corrections/issue-efactor.tex`, list (p94′), (422′), (428′), (p95′), (p96′), (428a), (p98′), (448a); add the p. 96 display (e³ with charge-Z*) and a new row (or an addition to row 54) for the p. 101 bubble missing e²; mark row 28 as the same convention (234′, 235′, 237′, 312′); state that the margin note can be removed.

## 11. Open points

- Fix A (number density, explicit e) is recommended; the authors might prefer Fix B/C (Section 3 lists the extra edits each needs).
- Item 6 (Hamiltonian convention) overlaps errata row 28; if that is issued separately, item 6 can be dropped from this correction without affecting items 1–5.
- The overall sign of `g^{dyn,0}` relative to the Hamiltonian coefficient (`p95′ = −g⁰_VG` with e < 0) is a convention of the diagrammatic factorization of E.411–E.416. It cannot be absorbed in `e_{aνσ}(q)` and it is unobservable inside EPP's self-energy. It matters only when (p95′) is combined with Hamiltonian-convention vertices (e.g. VG's `g^S + g^L`).
- Not addressed here: ω_σ (bare) vs ω_LO in the zero-point amplitude; missing τ_κ phases (as in VG Eq. 4); ε^e is the RPA head without local fields; (428′) assumes a local H⁰_e.
- GaAs inputs: Z* = 2.07|e| and the measured ω_TO, ω_LO, ε∞ are mutually inconsistent at the 5 % level in |g| (Z*_LST = 2.18|e|).

## 12. Review responses

Two skeptics raised 14 objections (several are duplicates). Each was checked against the PDF and the derivation. All are accepted; one is accepted with a refinement (O2/O7, Fix C entry). Unless noted, no number changed. The script was extended and re-run, and all assertions pass.

| # | Objection (short) | Check | Action |
|---|---|---|---|
| O1 | Fix B/C list omits E.270, E.273 (A² terms) and E.279 (Π^{e,dia}) | PDF p. 56: E.270, E.273 carry `e²/(2m_e c²)∫ρ^e A²`; p. 58: E.279 `−(e²/mc)ρ^e δ δ_αβ`. Under `ρ_e = eΨ⁺Ψ` each goes `e² → e`, like E.237, E.312, E.321 (p. 69). | Accepted. Added to the tex ("Choice of fix", both the Fix A justification and the alternative's list) and to the report (§1 table, §3 table, §8 "Do not change"). |
| O2, O7 | §3 table: under B/C, E.423–E.425 etc. are not "consistent as printed" | PDF pp. 94–95: only E.422 is written with ρ_e; the left sides of E.423, E.424, the two displays between E.424 and E.425, and E.425 are written with E.420 substituted (`φ*φ c⁺c`, no e). The display after E.427 is written with ρ_{nk,mk+q}. | Accepted, with a refinement for Fix C. Under B every one of these displays needs e. Under C, E.423–E.425 and the displays between them need e, but the display after E.427 is consistent, because there ρ is the redefined charge vertex (C changes E.426 instead). The table is split into three rows. The tex now lists these displays for the alternative and says that within §7.2 B edits everything A edits except the continuity equation and E.422, plus eight edits elsewhere. The report text after the table says the same. |
| O3, O12 | "constant phase that can be absorbed in e_aνσ(q)" is wrong | (p95′) and (p96′) are linear in the same e_aνσ, so rephasing leaves p95′/g⁰ = −1. p98′ = −[p95′]^* q²/(q·ε·q) is fixed by the factor i in E.415/E.416. Only g^{dyn}χ⁰g^{dyn,0} enters E.413 (p. 92). Script §2f confirms: the ratios are −1 and +1 and the product is invariant for z = 1, −1, e^{0.7i}. | Accepted. The tex item 3 sentence is replaced: the sign is a convention of E.415–E.416 and is not absorbable; the product in E.413 is fixed and with Ω⁰ reproduces E.452; the sign cancels in ΔΩ^e and \|g\|². Also added: when combining with Hamiltonian-convention vertices (VG's g^S + g^L), use g⁰_mnσ = −lim g^{dyn,0}_nmσ. Report §5, §11 and primer item §10.1 updated. |
| O4, O8 | Justification (ii): "longitudinal projection of (375) does not depend on c" is false as written | PDF p. 83: q̂·D⁰·q̂ = −4πc/ω², linear in c; it equals E.374 (p. 83) = E.442 (p. 98) for every c, and its c cancels the 1/c of the E.417 prefactor (p. 93). | Accepted. The tex (ii) is reworded accordingly, adding the corresponding statement for D^e (`−4πc/(ω²ε_L)`), and report §4.2 is reworded the same way. New script check: q̂·D⁰·q̂/(4πc/ω²) = −1 to 10⁻¹³ and longitudinal part = p95′ to 1.3×10⁻¹⁴ for ω/cq = 10⁻³ … 0.6 (§7.2 table). |
| O5, O9 | Replacement prose before (417) drops the sign and the relabelling | PDF p. 24: E.123 writes the electron charge as −e (e > 0); p. 25: E.131 has ρ_{nk,mk−q}(−G). Relabelled (n ↔ m, k → k + q), E.131 at ε = 1, G = 0 is g⁰_mnσ, and p95′ = (e/\|e\|) g⁰_mnσ. | Accepted. The replacement sentence now reads "reduces, up to the overall sign e/\|e\| [see (428a)] and the relabelling n ↔ m, k → k + q, to …", with a bracketed note on the §3/§4 sign of e. The Justification and report §6.4 now state the relabelling. |
| O6 | Stern formula quoted without its range | Static 2D Lindhard: χ⁰ = −m*/(πħ²) for Q ≤ 2k_F at T = 0. The check has k_F = 0.0251 Å⁻¹, so Q = 0.004 < 2k_F = 0.0501 Å⁻¹, and k_BT = 1.72 meV ≪ E_F = 35.7 meV. | Accepted. The tex now says "valid for Q ≤ 2k_F at T → 0 and without a background dielectric constant". Report §4.5 updated. The script now asserts Q < 2k_F and k_BT < 0.1 E_F. |
| O10 | p. 94 "So the equation for charge density is" left unchanged | PDF p. 94: the sentence directly precedes `iħρ̇_e = [ρ_e, H⁰_e]_−`. | Accepted. Tex item 1 now replaces it with "So the equation for the electron (number) density is". Report §8 and §9 updated. |
| O11 | "M = 36.1 amu" is the reduced mass | 69.723 × 74.922/144.645 = 36.11 amu. The script uses it as μ in `√(ħ/2μω_LO)`. | Accepted. The tex now reads "reduced mass μ = 36.1 amu", and a comment was added in the script. No number changes. |
| O13 | "it is dimensionless" reads as if Ψ⁺Ψ were dimensionless | Ψ⁺Ψ is in cm⁻³; ρ_{nk,mk+q} = ∫_Ω e^{−iqr}φ*φ with φ normalized in Ω is dimensionless. | Accepted. The tex now says "ρ_e has the dimension of an inverse volume and the matrix element ρ_{nk,mk+q} is dimensionless". |
| O14 | Report header says EPP pp. 60–104; the PDF has 102 pages | `pdfinfo`: 102 pages. | Accepted. The header now reads 55–58, 60–102 and lists the pages re-read in this round. |

No objection was rejected.

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
  - The vertices (p95′), (p96′), (p98′) now also carry the basis phase. τ_κ = **r**⁰_a replaces τ_κ = 0 (dipoles at
    **R**⁰_i, §5 table row "τ_κ") in the (p96′) paragraph, after (448a) and in check (iv). Both sides of (428a) and
    (448a) carry the same phase, so the relations are unchanged (check_merge M1, M3).
  - (448a) now shows both forms in one display. The second form, lim g^dyn = −[lim g^{dyn,0}]* q²/(**q**·ε·**q**),
    is traced to the +i of both E.415 and E.416, with D⁰ (442) real. The overall sign e/|e| of (428a) compares the
    diagrammatic vertex with the Hamiltonian coefficient, which is −(p95′) for either sign of e (merged §5.2;
    check_merge M2 builds both from E.415/E.416 directly).
  - Displayed in addition: (p93c′) (Heisenberg evolution with ħ), (237′) (with the basis phase) and (312′). The
    remark before E.417 now refers to (131′).
  - GaAs: 304.05 meV with q→0 LO eigenvectors on both sides. With finite-q eigenvectors and the phase (force-constant
    model of the phase report) the common value at q = 0.1 Å⁻¹ is 303.0 meV along [111] and unchanged along [100].
