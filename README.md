# asymptoticRTE
asymptotic convergence of the diffusion equation to the radiative transfer equation in the time-domain and steady-state regimes

# asymptoticRTE

Reference implementations for photon **fluence** in an infinite, homogeneous medium using the **Radiative Transfer Equation (RTE)** and the **Diffusion Equation (DE)**. Includes both **time-domain** and **continuous-wave (CW)** regimes, with utilities to study asymptotic behavior and to test different expressions of the diffusion coefficient found in the literature.

The code was designed and tested in **GNU Octave**, aiming for MATLAB compatibility.

---

## Content

- **RTE (isotropic)**: `RTE_iso` (time-domain), `RTE_iso_CW` (CW).
- **RTE (anisotropic, HG)** via **P\_N** approximation: `RTE_aniso` (time-domain), `RTE_aniso_CW` (CW).
- **Diffusion equation** closed forms: `DE` (time-domain), `DE_CW` (CW).

---

## Quickstart

Running `setup.m` from the repo root adds all relevant folders to the path and generates (optionally) figures showing the time-domain and steady-state asymptotic convergence between the RTE ($\phi$) and the DE ($\psi$).
