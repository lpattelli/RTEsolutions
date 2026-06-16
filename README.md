# RTEsolutions

Reference implementations for photon fluence in an infinite, homogeneous
medium using the radiative transfer equation (RTE) and the diffusion equation
(DE). The code covers time-domain and continuous-wave (CW) regimes.

The code is designed for GNU Octave, with MATLAB-compatible syntax where
practical.

## Setup

From the repository root:

```matlab
setup
```

This adds `src/`, `utils/`, and `figures/` to the path. It does not generate
figures by default. To generate the production figure windows, run:

```matlab
Figure1
Figure2
Figure3
Figure4
```

## Source Functions

- `RTE_iso`: isotropic time-domain RTE fluence.
- `RTE_aniso`: time-domain PN fluence for Henyey-Greenstein anisotropy.
- `RTE_iso_CW`: isotropic CW RTE fluence.
- `RTE_aniso_CW`: CW PN fluence for Henyey-Greenstein anisotropy.
- `DE`: time-domain diffusion fluence.
- `DE_CW`: CW diffusion fluence.

## Asymptotic Helpers

- `cw_iso_asymptotics`: exact isotropic CW pole/residue data.
- `cw_aniso_asymptotics`: direct PN dominant-pole and residue data for HG
  anisotropy.
- `compute_a`: grid computation of the pole-matched parameter `a` using the
  direct PN dominant pole.
- `HG_high_albedo_coeffs`: leading high-albedo coefficients for HG scattering.
- `holte_hg_asymptotics`: legacy comparison helper for the approximate Holte HG
  expansion.

## Production Figures

- `Figure1`: time-domain pointwise RTE/DE convergence.
- `Figure2`: time-domain limits for alternative diffusion coefficients.
- `Figure3`: CW radial RTE/DE comparisons.
- `Figure4`: CW anisotropic pole, ratio, and residue maps.

## Quickstart

Plot the time-domain ratio between anisotropic RTE and DE for `g = 0.9`:

```matlab
setup

mm = 1e-3;
ps = 1e-12;
c = 0.299792458*mm/ps;

r    = 10*mm;
mua  = 0.01/mm;
g    = 0.9;
musp = 1/mm;
mus  = musp/(1 - g);
D    = (1/3)/musp;

t = logspace(0, 5, 501)*ps;
t = t(t > r/c);

phi = RTE_aniso(t, r, mua, mus, g, c);
psi = DE(t, r, D, mua, c);

figure;
semilogx(t/ps, phi ./ psi);
xlabel('t [ps]');
ylabel('\phi_{RTE}/\psi_{DE}');
grid on;
```
