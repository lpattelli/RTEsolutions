function [a_pole_matrix, ratio_matrix, g_vals, mua_vals, D_pole_matrix, D_res_matrix] = compute_a(g_vals, eta_vals, musp, N)
% COMPUTE_A - Compute exponent-matched CW parameter a from the dominant PN pole
%
% Syntax:
%   [a_pole_matrix, ratio_matrix, g_vals, mua_vals] = compute_a
%   [a_pole_matrix, ratio_matrix, g_vals, mua_vals, D_pole_matrix, D_res_matrix] = compute_a(g_vals, eta_vals, musp, N)
%
% Description:
%   Computes the pole-matched parameter a on a grid of absorption and
%   anisotropy values. The calculation uses the dominant PN eigenvalue of the
%   CW transport matrix, not a finite-radius best-fit search.
%
%   The diffusion coefficient is parameterized as
%      D = 1/[3(mus' + a*mua)].
%
% Inputs:
%   g_vals   [1xG double, -]          HG anisotropy values (default linspace(0,0.99,34))
%   eta_vals [1xM double, -]          mua/mus' values (default logspace(-5,0,51))
%   musp     [1x1 double, length^-1]  Reduced scattering coefficient (default 1)
%   N        [1x1 int]                PN order (default 101)
%
% Outputs:
%   a_pole_matrix [GxM double, -]      Pole-matched a values
%   ratio_matrix  [GxM double, -]      lim_{r->inf} Phi_RTE/Phi_DE at D=D_pole
%   g_vals        [1xG double, -]      Grid used along anisotropy
%   mua_vals      [1xM double, length^-1] Absorption grid
%   D_pole_matrix [GxM double, length] Dominant-pole matched D
%   D_res_matrix  [GxM double, length] Dominant-residue matched D

  if nargin < 1 || isempty(g_vals)
    g_vals = linspace(0, 0.99, 34);
  end
  if nargin < 2 || isempty(eta_vals)
    eta_vals = logspace(-5, 0, 51);
  end
  if nargin < 3 || isempty(musp)
    musp = 1;
  end
  if nargin < 4 || isempty(N)
    N = 101;
  end

  mua_vals = eta_vals * musp;

  a_pole_matrix = zeros(numel(g_vals), numel(mua_vals));
  ratio_matrix  = zeros(numel(g_vals), numel(mua_vals));
  D_pole_matrix = zeros(numel(g_vals), numel(mua_vals));
  D_res_matrix  = zeros(numel(g_vals), numel(mua_vals));

  for ig = 1:numel(g_vals)
    for im = 1:numel(mua_vals)
      [D_pole, D_res, a_pole, ~, ratio] = cw_aniso_asymptotics(mua_vals(im), musp, g_vals(ig), N);

      a_pole_matrix(ig, im) = a_pole;
      ratio_matrix(ig, im)  = ratio;
      D_pole_matrix(ig, im) = D_pole;
      D_res_matrix(ig, im)  = D_res;
    end
  end
end
