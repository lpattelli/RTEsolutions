function Figure2
%% Figure2 - CW regime analysis, optimal D expressions and asymptotics
%
% Authors: André Liemert (ILM-ULM)
%          Lorenzo Pattelli (INRIM)
%          Fabrizio Martelli (UNIFI)
% License: MIT

  mm = 1e-3;

  %% Figure 2a - Check different expressions of the diffusion coefficient
  r   = logspace(-1, 3.6, 2301)*mm;
  mua = 0.01/mm;
  mus = 1/mm;

  a = [0, 1/5, 0.200228571055617, 1/3, 1];
  labels = arrayfun(@(x) sprintf('a=%.3g', x), a, 'UniformOutput', false);

  ref = RTE_iso_CW(r, mua, mus);

  figure
  hold on
  set(gca, 'ColorOrder', cool(numel(a)));
  title(sprintf('Figure 2a: RTE/DE ratio (CW) (\\mu_s = %g m^{-1}, \\mu_a = %g m^{-1})', mus, mua))
  for idx = 1:numel(a)
    D = (1/3)/(mus + a(idx)*mua);
    ratio = ref ./ DE_CW(r, D, mua);
    semilogx(r/mm, ratio, 'DisplayName', labels{idx})
  end
  line([r(1),r(end)]/mm, [1,1], 'linestyle', '--', 'color', 'k', 'HandleVisibility', 'off')
  xlim([r(1),r(end)]/mm)
  ylim([0.98, 1.02])
  xlabel('r [mm]')
  ylabel('\phi/\psi')
  set(gca, 'xscale', 'log')
  legend('show')

  %% Figure 2b - Different absorption levels
  muas = logspace(-5, -1, 5)/mm;
  as   = [0.200000228687555, 0.200002285706599, 0.20002285714323, 0.200228571055617, 0.202285284334453];

  figure
  hold on
  set(gca, 'ColorOrder', cool(numel(muas)));
  title('Figure 2b: CW convergence for different \mu_a')
  for idx = 1:numel(muas)
    D = (1/3)/(mus + as(idx)*muas(idx));
    ratio = RTE_iso_CW(r, muas(idx), mus) ./ DE_CW(r, D, muas(idx));
    semilogx(r/mm, ratio, 'DisplayName', sprintf('\\mu_a/\\mu_s = %g, a = 0.2 + %g', muas(idx)/mus, as(idx)-0.2))
  end
  line([r(1), r(end)]/mm, [1,1], 'linestyle', '--', 'color', 'k', 'HandleVisibility', 'off')
  xlim([r(1), r(end)]/mm)
  ylim([0.98, 1.02])
  xlabel('r [mm]')
  ylabel('\phi/\psi')
  set(gca, 'xscale', 'log')
  legend('show')

  %% Figure 2c–d - Surfaces: optimal a and asymptotic value
  % load values from pre-computed database
  baseDir = fileparts(mfilename('fullpath'));
  rootDir = fileparts(baseDir);
  dataDir = fullfile(rootDir, 'data');
  load(fullfile(dataDir, 'a_values.mat'));  % provides: g_vals, mua_vals, mus_reduced, best_a_matrix, avg_ratio_matrix

  % grids
  [G, MUA_ratio_log10] = meshgrid(g_vals, log10(mua_vals/mus_reduced));
  Z1 = best_a_matrix.';
  Z2 = avg_ratio_matrix.';

  figure

  for panel = 1:2
    subplot(1,2,panel)

    if panel==1
      Z = Z1;
    else
      Z = Z2;
    end

    surf(G, MUA_ratio_log10, Z, 'EdgeColor','none','FaceAlpha',0.5);
    hold on
    xlabel('g')
    ylabel('log_{10}(\mu_a/\mu_s'')')

    if panel==1
      zlabel('a')
      levels = 0.21:0.09:0.57;
    else
      zlabel('lim_{r \rightarrow \infty} \phi/\psi')
      levels = [0.09, 0.9, 0.99, 0.999, 0.9999, 0.99999];
    end

    % get current colormap and color limits
    cmap = colormap;
    nColors = size(cmap,1);
    clim = caxis;

    % draw contours projected on base plane for context
    C = contourc(g_vals, log10(mua_vals/mus_reduced), Z, levels);
    idx = 1;
    while idx < size(C,2)
      lvl  = C(1,idx);
      npts = C(2,idx);
      xs   = C(1, idx+1:idx+npts);
      ys   = C(2, idx+1:idx+npts);

      % map contour level to colormap index
      t = (lvl - clim(1)) / (clim(2) - clim(1));
      t = max(0, min(1, t));  % clamp to [0,1]
      colorIdx = round(1 + t * (nColors - 1));
      lineColor = cmap(colorIdx, :);

      plot3(xs, ys, zeros(1,npts), '-', 'LineWidth', 1.2, 'Color', lineColor)
      idx = idx + npts + 1;
    end
    view(60,30)
    hold off
  end

  clim = caxis;
  cmap = colormap;

  S = axes('visible', 'off', 'title', 'Figure 2c-d: convergence parameter and convergence value');

  %% Figure 2e-f - Holte/Aronson expansion prediction on the same grid
  % Implements Aronson&Corngold Eq.(14)-(15) (Holte expansion) with HG: f_l=g^l.

  best_a_holte   = zeros(size(best_a_matrix));   % same orientation as loaded matrices: (g_idx, mua_idx)
  avg_ratio_holte = zeros(size(avg_ratio_matrix));

  for ig = 1:numel(g_vals)
    g = g_vals(ig);

    for im = 1:numel(mua_vals)
      mua  = mua_vals(im);
      musp = mus_reduced;            % this is mu_s'' in your notation
      [aH, ratioH] = holte_a_ratio(mua, musp, g);

      best_a_holte(ig, im)    = aH;
      avg_ratio_holte(ig, im) = ratioH;
    end
  end

  % grids already defined above:
  % [G, MUA_ratio_log10] = meshgrid(g_vals, log10(mua_vals/mus_reduced));

  Z1H = best_a_holte.';        % transpose to match meshgrid shape (numel(mua_vals) x numel(g_vals))
  Z2H = avg_ratio_holte.';

  figure
  for panel = 1:2
    subplot(1,2,panel)

    if panel==1
      Z = Z1H;
    else
      Z = Z2H;
    end

    surf(G, MUA_ratio_log10, Z, 'EdgeColor','none','FaceAlpha',0.5);
    hold on
    xlabel('g')
    ylabel('log_{10}(\mu_a/\mu_s'')')

    if panel==1
      zlabel('a (Holte Eq.14)')
      levels = 0.21:0.09:0.57;
    else
      zlabel('lim_{r \rightarrow \infty} \phi/\psi (Holte-based)')
      levels = [0.09, 0.9, 0.99, 0.999, 0.9999, 0.99999];
    end

    % get current colormap and color limits
    cmap = colormap;
    nColors = size(cmap,1);
    clim = caxis;

    % draw contours projected on base plane for context
    C = contourc(g_vals, log10(mua_vals/mus_reduced), Z, levels);
    idx = 1;
    while idx < size(C,2)
      lvl  = C(1,idx);
      npts = C(2,idx);
      xs   = C(1, idx+1:idx+npts);
      ys   = C(2, idx+1:idx+npts);

      % map contour level to colormap index
      t = (lvl - clim(1)) / (clim(2) - clim(1));
      t = max(0, min(1, t));
      colorIdx = round(1 + t * (nColors - 1));
      lineColor = cmap(colorIdx, :);

      plot3(xs, ys, zeros(1,npts), '-', 'LineWidth', 1.2, 'Color', lineColor)
      idx = idx + npts + 1;
    end

    view(60,30)
    hold off
  end

  S = axes('visible', 'off', 'title', 'Figure 2e-f: Holte/Aronson expansion prediction');


function [a, ratio] = holte_a_ratio(mua, musp, g)
  % Holte expansion as quoted by Aronson & Corngold (1999), Eq.(14)-(15)
  % musp = mu_s' (reduced scattering), mu_tr = musp + mua
  % HG phase function: f_l = g^l

  mu_tr = musp + mua;

  % convert mu_s' -> mu_s assuming HG with anisotropy g: mu_s' = mu_s(1-g)
  denom = max(1 - g, 1e-12);
  mus   = musp / denom;
  mu_t  = mus + mua;

  omega = mus / mu_t;    % = \varpi in Aronson&Corngold

  % HG Legendre moments
  f1 = g;
  f2 = g^2;
  f3 = g^3;
  f4 = g^4;

  % h_l = (2l+1)(1 - omega f_l); with f0=1
  h0 = 1 * (1 - omega);         % = mua/mu_t
  h1 = 3 * (1 - omega*f1);      % = 3 mu_tr/mu_t
  h2 = 5 * (1 - omega*f2);
  h3 = 7 * (1 - omega*f3);
  h4 = 9 * (1 - omega*f4);

  % protect divisions
  h2 = max(h2, 1e-15);
  h3 = max(h3, 1e-15);
  h4 = max(h4, 1e-15);

  x = 4*h0/h2;

  % Bracket in Eq.(14) (truncated at the order written in Aronson&Corngold)
  B = 1 ...
      - x ...
      + x^2 * (1 - (9*h1)/(4*h3)) ...
      - x^3 * ( 1 - (27*h1)/(4*h3) + (81/16)*(h1^2/(16*h3^2)) + 9*(h1^2/h3^2)*(h2/h4) );

  % They state: 1/D = 3 mu_tr * B => D = 1/(3 mu_tr B)
  % Mapping to our parameterization D = 1/(3(musp + a*mua)):
  % musp + a*mua = mu_tr * B  -> a = (mu_tr*B - musp)/mua
  a = (mu_tr*B - musp) / mua;

  % Asymptotic amplitude ratio (exact for isotropic when B is exact; here used as Holte-based predictor)
  % Derived from the dominant-mode prefactor vs diffusion prefactor after enforcing exponent matching:
  ratio = 2*(mu_tr - 3*mua*B) / (musp*(3*B - 1));

  % clean up numerics near the breakdown region
  if ~isfinite(a),     a = NaN;     end
  if ~isfinite(ratio), ratio = NaN; end
  ratio = max(0, min(1, ratio));
end

end

