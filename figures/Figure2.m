function Figure2
% FIGURE2 - Time-domain limits for alternative diffusion coefficients
%
% Description:
%   Generates panels 2a-b: the finite-time RTE/DE ratio for several diffusion
%   coefficients and the corresponding analytical infinite-time limits.
%
% Authors: Andre Liemert (ILM-ULM)
%          Lorenzo Pattelli (INRIM)
%          Fabrizio Martelli (UNIFI)
% License: MIT

  mm = 1e-3;
  ps = 1e-12;
  c  = 0.299792458*mm/ps;

  figure('Name', 'Figure 2: time-domain modified diffusion');

  %% Figure 2a - RTE/DE ratio for different expressions of D
  mua = 0.01/mm;
  mus = 1/mm;
  r   = 10*mm;
  t   = logspace(2, 5, 301)*ps;

  a_values = linspace(0, 1, 11);
  ref = RTE_iso(t, r, mua, mus, c);

  subplot(1, 2, 1);
  hold on;
  set(gca, 'ColorOrder', cool(numel(a_values)));
  for ia = 1:numel(a_values)
    D = (1/3) / (mus + a_values(ia)*mua);
    ratio = ref ./ DE(t, r, D, mua, c);
    semilogx(t/ps, ratio, 'DisplayName', sprintf('a = %.2g', a_values(ia)));
  end
  line([t(1), t(end)]/ps, [1, 1], 'linestyle', '--', 'color', 'k', 'HandleVisibility', 'off');
  xlim([t(1), t(end)]/ps);
  ylim([0.97, 1.03]);
  xlabel('t [ps]');
  ylabel('\phi/\psi');
  title('Figure 2a');
  legend('location', 'northeast');
  grid on;

  %% Figure 2b - Analytical infinite-time limit
  muas = logspace(-4, 4, 81)/mm;
  a_lines = [0, 0.01, 0.1, 1];

  subplot(1, 2, 2);
  hold on;
  set(gca, 'ColorOrder', cool(numel(a_lines)));
  for ia = numel(a_lines):-1:1
    xaxis = muas / mus;
    yaxis = exact_limit(a_lines(ia), muas, mus);
    loglog(xaxis, yaxis, 'DisplayName', sprintf('a = %.2g', a_lines(ia)));
  end
  ylim([1e-4, 1e1]);
  xlabel('\mu_a/\mu_s');
  ylabel('lim_{t \rightarrow \infty} \phi/\psi');
  title('Figure 2b');
  legend('location', 'southwest');
  grid on;
end


function value = exact_limit(a, mua, mus)
% EXACT_LIMIT - Analytical asymptotic limit for the time-domain RTE/DE ratio
%
% Syntax:
%   value = exact_limit(a, mua, mus)
%
% Inputs:
%   a    [array double, -]            Parameter in D = 1/[3(mus+a*mua)]
%   mua  [array double, length^-1]    Absorption coefficient
%   mus  [array double, length^-1]    Scattering coefficient
%
% Output:
%   value [array double, -]           lim_{t->inf} Phi_RTE/Phi_DE

  value = (mus ./ (a.*mua + mus)).^(3/2);
end
