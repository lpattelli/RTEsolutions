function Figure3
% FIGURE3 - CW radial RTE/DE comparisons
%
% Description:
%   Generates panels 3a-b: the CW RTE/DE radial ratio for several scalar
%   diffusion coefficients and for several absorption levels using the
%   pole-matched coefficient.
%
% Authors: Andre Liemert (ILM-ULM)
%          Lorenzo Pattelli (INRIM)
%          Fabrizio Martelli (UNIFI)
% License: MIT

  mm = 1e-3;
  mus = 1/mm;
  r = logspace(-1, 3.6, 2301)*mm;

  figure('Name', 'Figure 3: CW radial comparison');

  %% Figure 3a - Diffusion-coefficient choices at fixed optical parameters
  mua = 0.01/mm;
  [~, a_pole] = cw_iso_asymptotics(mua, mus);
  a_values = [0, 1/5, a_pole, 1/3, 1];

  ref = RTE_iso_CW(r, mua, mus);

  subplot(1, 2, 1);
  hold on;
  set(gca, 'ColorOrder', cool(numel(a_values)));
  for ia = 1:numel(a_values)
    D = diffusion_from_a(a_values(ia), mua, mus);
    ratio = ref ./ DE_CW(r, D, mua);
    semilogx(r/mm, ratio, 'DisplayName', sprintf('a = %.4g', a_values(ia)));
  end
  line([r(1), r(end)]/mm, [1, 1], 'linestyle', '--', 'color', 'k', 'HandleVisibility', 'off');
  xlim([r(1), r(end)]/mm);
  ylim([0.98, 1.02]);
  xlabel('r [mm]');
  ylabel('\phi/\psi');
  title('Figure 3a');
  legend('show');
  legend('location', 'northeast');
  grid on;

  %% Figure 3b - Pole-matched coefficient for different absorption levels
  muas = logspace(-5, -1, 5)/mm;

  subplot(1, 2, 2);
  hold on;
  set(gca, 'ColorOrder', cool(numel(muas)));
  for im = 1:numel(muas)
    [~, a_pole] = cw_iso_asymptotics(muas(im), mus);
    D = diffusion_from_a(a_pole, muas(im), mus);
    ratio = RTE_iso_CW(r, muas(im), mus) ./ DE_CW(r, D, muas(im));
    semilogx(r/mm, ratio, 'DisplayName', sprintf('\\mu_a/\\mu_s = %.1e', muas(im)/mus));
  end
  line([r(1), r(end)]/mm, [1, 1], 'linestyle', '--', 'color', 'k', 'HandleVisibility', 'off');
  xlim([r(1), r(end)]/mm);
  ylim([0.98, 1.02]);
  xlabel('r [mm]');
  ylabel('\phi/\psi');
  title('Figure 3b');
  legend('show');
  legend('location', 'northeast');
  grid on;
end


function D = diffusion_from_a(a, mua, mus)
% DIFFUSION_FROM_A - Scalar diffusion coefficient parameterized by a

  D = 1 / (3 * (mus + a*mua));
end
