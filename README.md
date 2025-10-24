# asymptoticRTE

Reference implementations for photon fluence in an infinite, homogeneous medium using the Radiative Transfer Equation (RTE) and the Diffusion Equation (DE). Includes both time-domain and continuous-wave (CW) regimes, with utilities to study asymptotic behavior and to test different expressions of the diffusion coefficient found in the literature.

The code was designed and tested in GNU Octave, aiming for MATLAB compatibility.


## Content

- RTE (isotropic): `RTE_iso` (time-domain), `RTE_iso_CW` (CW).
- RTE (anisotropic, HG) via P<sub>N</sub> approximation: `RTE_aniso` (time-domain), `RTE_aniso_CW` (CW).
- Diffusion equation closed forms: `DE` (time-domain), `DE_CW` (CW).


## Setup

Running `setup.m` from the repo root adds all relevant folders to the path and generates (optionally) figures showing the time-domain and steady-state asymptotic convergence between the RTE ($\phi$) and the DE ($\psi$).


## Quickstart

Plot the ratio between the RTE and the DE for anisotropic scattering ($g=0$)

``` matlab
% time-domain convergence of the RTE fluence for g = 0.9 and the DE fluence
mm = 1e-3;
ps = 1e-12;
c = 0.299792458*mm/ps;

r    = 10*mm;         % source–detector distance
mua  = 0.01/mm;       % absorption
g    = 0.9;           % anisotropy
musr = 1/mm;          % reduced scattering
mus  = musr/(1 - g);  % scattering used by RTE_aniso
D    = (1/3)/musr;    % diffusion coefficient for DE

% time axis > r/c
t = logspace(0, 5, 501)*ps;
t = t(t > r/c);

% fluence (time-domain)
phi = RTE_aniso(t, r, mua, mus, g, c);
psi = DE(t, r, D, mua, c);

% ratio
ratio = phi ./ psi;

figure; hold on
title('time-domain ratio: RTE\_aniso (g = 0.9) vs DE')
semilogx(t/ps, ratio)
line([t(1), t(end)]/ps, [1, 1], 'LineStyle', '--', 'Color', 'k')
xlabel('t [ps]'); ylabel('\phi_{RTE}/\psi_{DE}')
xlim([t(1), t(end)]/ps)
ylim([0, 1.2])
grid on

```
