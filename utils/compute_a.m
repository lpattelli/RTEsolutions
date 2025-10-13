function [best_a_matrix, avg_ratio_matrix] = compute_a
  % Compute best "a" over a mesh of mua and g values, plot and save results

  mm = 1e-3;
  r = logspace(-3, 4, 701) * mm;
  mus_reduced = 1 / mm;
  % Define mesh
  mua_vals = [logspace(-5, -2, 16), ...
              logspace(-2, -1, 11)(2:end), ...
              logspace(-1, 0, 51)(2:end)] / mm;   % fine mesh, long runtime ...
  g_vals = [0.0:0.01:0.99];                       % fine mesh, long runtime ...

  best_a_matrix = zeros(length(g_vals), length(mua_vals));
  avg_ratio_matrix = zeros(length(g_vals), length(mua_vals));

  % optimization settings
  options = optimset('TolX', 1e-15, 'TolFun', 1e-15, 'MaxIter', 1000);
  Npoints = 20;  % number of (non-NaN) tail points to be considered

  % loop over g first, then mua (variation with constant mua does not show abrupt changes)
  for i = 1:length(g_vals)
      g = g_vals(i);
      mus = mus_reduced / (1 - g);
      a_bnd = [0, 1];

      for j = 1:length(mua_vals)
          mua = mua_vals(j);

          % apply last successful boundaries
          best_a = fminbnd(@(a) convergence_metric(a, r, mua, mus, g, Npoints), a_bnd(1), a_bnd(2), options);
          best_a_matrix(i,j) = best_a;
          [~, avg_ratio] = convergence_metric(best_a, r, mua, mus, g, Npoints);
          avg_ratio_matrix(i,j) = avg_ratio;  % call the function one more time to get the asymptotic ratio

          % print progress info
          fprintf('\rProgress: %d / %d', ((i - 1) * length(mua_vals) + j), (length(mua_vals) * length(g_vals)));

          % define narrow boundaries for the next iteration
          a_bnd = best_a * [0.9, 1.1];
      end
  end
  fprintf('\n'); % clear console

  % save raw output to the data/ folder
  baseDir = fileparts(mfilename('fullpath'));  % utils/
  rootDir = fileparts(baseDir);
  dataDir = fullfile(rootDir, 'data');

  save(fullfile(dataDir, 'a_values.mat'), ...
    'mm', 'r', 'mua_vals', 'mus_reduced', 'g_vals', 'best_a_matrix', 'avg_ratio_matrix');
end

function [score, avg_ratio] = convergence_metric(a, r, mua, mus, g, Npoints)
  % Evaluates how stable the asymptotic ratio tail is for a given "a"

  D = 1 / (3 * (mus * (1 - g) + a * mua));
  ratio = RTE_aniso_CW(r, mua, mus, g) ./ DE_CW(D, r, mua);

  % remove NaN values
  valid_indices = find(~isnan(ratio));
  ratio = ratio(valid_indices);
  r = r(valid_indices);

  % discard last valid point (numerically unstable)
  ratio = ratio(1:end-1);
  r = r(1:end-1);

  % compute convergence score and convergence value
  score = mean(abs(diff(ratio(end-Npoints:end))));
  avg_ratio = mean(ratio(end-Npoints:end));
end

