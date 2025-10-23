function [best_a_matrix, avg_ratio_matrix] = compute_a
% COMPUTE_A - Optimize parameter 'a' over (mua, g) grid and save database
%
% Syntax:
%   [best_a_matrix, avg_ratio_matrix] = compute_a
%
% Description:
%   For each (mua, g), solves for the value of a in D = (1/3)/(mus + a*mua)
%   that maximizes asymptotic agreement between DE_CW and RTE_aniso_CW.
%   Users can set a custom grid of (mua, g) values, and tune convergence settings.
%   Results are saved to data/a_values.mat.
%
% Outputs:
%   best_a_matrix  [GxM double]  Optimal a over grid
%   avg_ratio_matrix [GxM double] Asymptotic ratio (convergence metric)

  mm = 1e-3;
  r = logspace(-3, 4, 701) * mm;
  mus_reduced = 1 / mm;

  % build mesh by concatenating different mua ranges
  tmp1 = logspace(-5, -2, 16);
  tmp2 = logspace(-2, -1, 11);
  tmp3 = logspace(-1,  0, 51);
  mua_vals = [tmp1, tmp2(2:end), tmp3(2:end)] / mm;  % fine mesh

  g_vals = [0.00:0.01:0.99];

  best_a_matrix  = zeros(length(g_vals), length(mua_vals));
  avg_ratio_matrix = zeros(length(g_vals), length(mua_vals));

  % optimization settings
  options = optimset('TolX', 1e-15, 'TolFun', 1e-15, 'MaxIter', 1000);
  Npoints = 20;  % number of tail points used in the metric

  % loop g then mua (smoother evolution along g)
  for i = 1:length(g_vals)
    g   = g_vals(i);
    mus = mus_reduced / (1 - g);
    a_bnd = [0, 1];

    for j = 1:length(mua_vals)
      mua = mua_vals(j);

      % optimize a and compute its corresponding asymptotic ratio
      best_a = fminbnd(@(a) convergence_metric(a, r, mua, mus, g, Npoints), a_bnd(1), a_bnd(2), options);
      best_a_matrix(i,j) = best_a;

      [~, avg_ratio] = convergence_metric(best_a, r, mua, mus, g, Npoints);
      avg_ratio_matrix(i,j) = avg_ratio;

      % progress
      fprintf('\rProgress: %d / %d', ((i - 1) * length(mua_vals) + j), numel(best_a_matrix));

      % narrow bracket for the next mua at same g (improves convergence)
      a_bnd = best_a * [0.9, 1.1];
    end
  end
  fprintf('\n');

  % save raw output to the data/ folder
  baseDir = fileparts(mfilename('fullpath'));  % utils/
  rootDir = fileparts(baseDir);
  dataDir = fullfile(rootDir, 'data');

  save('-v7', fullfile(dataDir, 'a_values.mat'), ...
    'mm', 'r', 'mua_vals', 'mus_reduced', 'g_vals', 'best_a_matrix', 'avg_ratio_matrix');
end


function [score, avg_ratio] = convergence_metric(a, r, mua, mus, g, Npoints)
% CONVERGENCE_METRIC - Stability metric of the asymptotic RTE/DE CW ratio
%
% Syntax:
%   [score, avg_ratio] = convergence_metric(a, r, mua, mus, g, Npoints)
%
% Inputs:
%   a        [1x1 double, -]             Absorption weight in D = (1/3)/(mus + a*mua)
%   r        [1xR double, length]         Radii for CW evaluation (tail analyzed)
%   mua      [1x1 double, length^-1]      Absorption coefficient
%   mus      [1x1 double, length^-1]      Scattering coefficient
%   g        [1x1 double, -]              HG anisotropy factor
%   Npoints  [1x1 int]                    # of tail points for averaging
%
% Outputs:
%   score     [1x1 double, -]     Mean absolute adjacent diff over last Npoints
%   avg_ratio [1x1 double, -]     Mean ratio over last Npoints

  D = 1 / (3 * (mus * (1 - g) + a * mua));

  % RTE / DE (CW)
  ratio = RTE_aniso_CW(r, mua, mus, g) ./ DE_CW(r, D, mua);

  % remove NaNs as well as the last (typically unstable) point
  valid_idx = find(~isnan(ratio));
  ratio = ratio(valid_idx(1:end-1));

  % compute convergence score and average over last Npoints
  score = mean(abs(diff(ratio(end-Npoints:end))));
  avg_ratio = mean(ratio(end-Npoints:end));
end

