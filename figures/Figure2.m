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

      plot3(xs, ys, zeros(1,npts), '-', 'LineWidth', 1.2)
      idx = idx + npts + 1;
    end
    view(60,30)
    hold off
  end

  clim = caxis;
  cmap = colormap;

  S = axes('visible', 'off', 'title', 'Figure 2c-d: convergence parameter and convergence value');
end

