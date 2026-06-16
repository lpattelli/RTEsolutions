function Figure1
% FIGURE1 - Time-domain pointwise RTE/DE convergence
%
% Description:
%   Generates panels 1a-b: the isotropic time-domain RTE/DE ratio as a
%   function of time and distance, and the anisotropic similarity check for a
%   fixed source-detector distance.
%
% Authors: Andre Liemert (ILM-ULM)
%          Lorenzo Pattelli (INRIM)
%          Fabrizio Martelli (UNIFI)
% License: MIT

  mm = 1e-3;
  ps = 1e-12;
  c  = 0.299792458*mm/ps;

  figure('Name', 'Figure 1: time-domain convergence');

  %% Figure 1a - RTE/DE ratio vs time and distance
  mua = 0.01/mm;
  mus = 1/mm;
  D   = (1/3)/mus;

  t  = logspace(1, 5, 801)*ps;
  rs = logspace(-1, 2, 10)*mm;
  zlimits = [0.9, 1.1];
  data = nan(numel(rs), numel(t));

  for ir = 1:numel(rs)
    postballistic = t > rs(ir)/c;
    data(ir, postballistic) = RTE_iso(t(postballistic), rs(ir), mua, mus, c) ./ ...
                              DE(t(postballistic), rs(ir), D, mua, c);
    invalid_tail = isnan(data(ir,:)) | ~isfinite(data(ir,:)) | data(ir,:) < 0.01;
    last_invalid = find(invalid_tail, 1, 'last');
    if ~isempty(last_invalid) && last_invalid < numel(t)
      data(ir, 1:last_invalid) = NaN;
    end
  end
  plotdata = data;
  plotdata(plotdata < zlimits(1) | plotdata > zlimits(2)) = NaN;

  subplot(1, 2, 1);
  hold on;
  set(gca, 'ColorOrder', cool(numel(rs)));
  for ir = numel(rs):-1:1
    plot3(t/ps, rs(ir)*ones(size(t))/mm, plotdata(ir, :));
  end
  zlim(zlimits);
  set(gca, 'xscale', 'log', 'yscale', 'log');
  xlabel('t [ps]');
  ylabel('r [mm]');
  zlabel('\phi/\psi');
  title('Figure 1a');
  view([30, 10]);
  grid on;

  %% Figure 1b - Similarity relations under anisotropic scattering
  musp = 1/mm;
  r    = 10*mm;
  t    = logspace(2, 5, 301)*ps;
  gs   = linspace(0, 0.9, 10);
  musv = musp ./ (1 - gs);

  D = (1/3)/musp;
  ref = DE(t, r, D, mua, c);

  subplot(1, 2, 2);
  hold on;
  set(gca, 'ColorOrder', cool(numel(gs)));
  for ig = 1:numel(gs)
    ratio = RTE_aniso(t, r, mua, musv(ig), gs(ig), c) ./ ref;
    semilogx(t/ps, ratio, 'DisplayName', sprintf('g = %.2g', gs(ig)));
  end
  line([t(1), t(end)]/ps, [1, 1], 'linestyle', '--', 'color', 'k', 'HandleVisibility', 'off');
  xlim([t(1), t(end)]/ps);
  ylim([0.97, 1.03]);
  xlabel('t [ps]');
  ylabel('\phi/\psi');
  title('Figure 1b');
  legend('show');
  legend('location', 'northeast');
  grid on;
end
