# 2DMaterialsOpticalProps

*A fully ab‑initio toolkit for linear–response and spectroscopy of quasi two‑dimensional (2D) materials using Quantum‑Electrodynamic Many‑Body Theory.*

---

## ✨ What this code does

1. **Starts from a Bloch basis** obtained from periodic planewave Density Functional Theory (DFT) calculation (currently interfaced only to [**Quantum ESPRESSO**](https://www.quantum-espresso.org/)) 
2. **Builds many‑body propagators** and in reciprocal space.
3. **Evaluates response functions** at successive levels of approximation:
   - **RPA** current‑current tensor ( `Optabs.f90` )
   - **\(G^0W^0\)** self‑energy & screened Coulomb interaction ( `Sloss.f90` )
   - **Bethe–Salpeter Equation (BSE)** for ladder electron‑hole correlations ( `BSE.f90` )
   - **Dyson equation for the photon propagator** in layered media ( `Photon.f90` )
4. **Outputs spectra** such as optical absorption, electron energy‑loss, dispersion relations for the characterisation of plasmon/exciton‑polaritons in arbitrary van der Waals heterostructures

The formalism follows the Quantum‑Electrodynamic Many‑Body framework detailed in *“Nanophotonics of Thin‑film Materials: a Many‑body Approach”* (PhD thesis, 2023) — with the present implementation also extended to includes GW self-energy corrections and self‑consistent screening.

---

## Repository layout

```
├── Modules/        ! shared modules (I/O, numerics, BZ, response funcs)
├── Chi_RPA/        ! k‑resolved & folded current–current response (RPA)
├── W/              ! screened Coulomb interaction W(q, ω)
├── Chi_ladd_up/    ! ↑‑spin ladder polarizability & BSE kernels
├── Chi_ladd_down/  ! ↓‑spin counterpart
├── BSE/            ! exciton solver & optical spectra
├── doc/            ! FORD automatic documentation (needs to be compiled)
├── Makefile        ! configure compilers/libraries & build all targets
└── README.md       ! this file
```

---

## Quick start

```bash
# Clone
$ git clone https://github.com/Nevensky/2DMaterialsOpticalProps.git
$ cd 2DMaterialsOpticalProps

# Edit the top of Makefile to select your Fortran compiler & BLAS/LAPACK
FC = ifort          # or ifx / gfortran / nvfortran …

# Build all executables
$ make -j8            # produces Optabs.x  Sloss.x  BSE.x  Photon.x
```

### 1. Ground‑state DFT

```
$ pw.x  < scf.in > scf.out
$ pw.x  < nscf.in > nscf.out     # dense k‑grid
$ bands.x –postproc outputs
```

### 2. Reads QE data → Bloch basis 

The module ( `Modules/read_qe.f90`) collects

- Bloch wave‑functions (`.wfc*`)
- cell vectors, symmetries, k-point meshe (in the ireducible Broullouin zone), etc. (`*.xml`)

You can explore these outputs with the helper script in `doc/qe2bloch.py` 
### 3. Linear‑response workflow

```bash
# RPA current‑current tensor
$ Optabs.x < optabs.in > optabs.log   # writes Chi_RPA/*

# GW self‑energy & screened Coulomb W
$ Sloss.x  < sloss.in  > sloss.log    # writes W/*

# Bethe–Salpeter equation
$ BSE.x    < bse.in    > bse.log      # writes Chi_ladd_*/  & spectra.dat

# Photon Dyson equation – obtain plasmon–polaritons, s‑SNOM maps …
$ Photon.x < photon.in > photon.log   # writes photon_modes.dat
```


---

## Prerequisites

| Requirement             | Tested versions                         |
|-------------------------|-----------------------------------------|
| Fortran compiler        |  ifort ≥ 2008, gfortran ≥ 12            |
| BLAS / LAPACK           | OpenBLAS, Intel MKL (recommended)       |
| FFT library (optional)  | FFTW3 (included in Intel MKL)           |
| OpenMP                  | Intel OpenMP multithreading recommended |
| MPI                     | Not implemented                         |
| Python ≥ 3.8 (optional) | numpy, matplotlib for post‑processing   |

> **Note:** OpenMP multithreading is highly recommended especially for large number of atoms

---

## Input files at a glance

| File        | Key variables (examples)                       | Description                   |
| ----------- | ---------------------------------------------- | ----------------------------- |
| `optabs.in` | `nk = 64 64`, `ecut = 15 Ry`, `broad = 50 meV` | RPA χ(𝐪,ω) parameters        |
| `sloss.in`  | `nq = 100`, `nω = 700`, `GW_scheme = G0W0`     | GW self‑energy, W(q,ω) grid   |
| `bse.in`    | `nband_v = 4`, `nband_c = 4`, `eps_inf = 3.5`  | BSE kernel & exciton settings |
| `photon.in` | `substrate = Al2O3`, `polarization = p`        | Layered‑medium Dyson solver   |

All input keywords are documented inside the example `*.in` files.

---

## How the theory fits together

```
DFT  →  RPA  →  GW Σ & W  →  BSE ladder  →  Dyson equation for D(𝐪,ω)
```

- **RPA** calculate the current-current (photonic) response arising from electron-hole excitations
- **GW** calculates the the single‑particle self-energy corrections and yields the screened Coulomb W
- **BSE** captures excitonic effects (electron-hole binding) within the ladder approximation vital for insulating and semiconducting materials
- **Dyson’s equation** re‑embeds the electronic kernel into a fully dynamic screened retarded photon propagator which correctly predicts both the transverse and longitudinal response properties. Its real and imaginary parts can be related to observable quantities (optical spectra, dispersion relations, renormalization of energies, relaxation rates, etc.), predicting plasmon-polaritons, exciton‑polaritons, etc.

The derivation up to the BSE level is covered in Chapters 1 and 4 of the thesis and related publications Ref [] [] []. GW is implemented according to Ref [].

---

## Citing this work

If you use *2DMaterialsOpticalProps* in your research, please cite
1. Golenić, N., de Gironcoli, S., & Despoja, V. (2024). Optically driven plasmons in graphene/hBN van der Waals heterostructures: simulating s-SNOM measurements. Nanophotonics, 13(15), 2765-2780.
2. Golenić, N., & Despoja, V. (2023). Trapped photons: Transverse plasmons in layered semiconducting heterostructures. Physical Review B, 108(12), L121402.
3. Golenić, N., de Gironcoli, S., & Despoja, V. (2024). Tailored plasmon polariton landscape in graphene/boron nitride patterned heterostructures. npj 2D Materials and Applications, 8(1), 37.
4. Novko, D., Šunjić, M., & Despoja, V. (2016). Optical absorption and conductivity in quasi-two-dimensional crystals from first principles: Application to graphene. Physical Review B, 93(12), 125413.
5. Despoja, V., Novko, D., Lončarić, I., Golenić, N., Marušić, L., & Silkin, V. M. (2019). Strong acoustic plasmons in chemically doped graphene induced by a nearby metal surface. Physical Review B, 100(19), 195401.
6. N. Golenić, Nanophotonics of Thin‑film Materials: a Many‑body Approach, PhD thesis, Univ. of Trieste (2023)
7. **Link to this repository**

---

## Contributing & support

Pull requests are welcome! Feel free to open an Issue for bug reports or feature requests. For help with running the code, please contact `neven.golenic [at] gmail.com`.

---

## License

Distributed under the **MIT License**. See `LICENSE` for details.

