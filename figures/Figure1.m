function Figure1
%% Figure1 - Asymptotic analysis of RTE and DE in time domain
%
% Authors: André Liemert (ILM-ULM)
%          Lorenzo Pattelli (INRIM)
%          Fabrizio Martelli (UNIFI)
% License: MIT
%
% Description:
%   Generates panels (a–d) comparing RTE and DE solutions

  % units and constants
  mm = 1e-3;
  ps = 1e-12;
  c  = 0.299792458*mm/ps;

  %% Figure 1a - RTE/DE ratio vs time and distance (isotropic case)
  mua = 0.01/mm;
  mus = 1/mm;
  D   = (1/3)/mus;

  t  = logspace(1,5,1001)*ps;
  rs = logspace(-1, 2, 10)*mm;
  data = nan(numel(rs), numel(t));

  for idx = 1:numel(rs)
    postballistic = t > rs(idx)/c;
    data(idx,postballistic) = RTE_iso(t(postballistic), rs(idx), mua, mus, c) ./ ...
                              DE(t(postballistic), rs(idx), D, mua, c);
    starttime = find(diff(data(idx,:) < 0.01), 1, 'last');
    data(idx, 1:starttime) = NaN;
  end

  figure
  hold on
  title('Figure 1a: asymptotic RTE/DE ratio convergence vs distance')
  set(gca, 'ColorOrder', cool(numel(rs)));
  for idx = numel(rs):-1:1
    plot3(t/ps, rs(idx)*ones(size(t))/mm, data(idx, :))
  end
  zlim([0.9, 1.1])
  set(gca, 'xscale', 'log', 'yscale', 'log')
  xlabel('t [ps]')
  ylabel('r [mm]')
  zlabel('\phi/\psi')
  view([30, 10])
  data(isinf(data)) = NaN;

  %% Figure 1b - Similarity relations under anisotropic scattering
  musp = 1/mm;
  r    = 10*mm;
  t    = logspace(2,5,301)*ps;
  gs   = linspace(0., 0.9, 10);
  musv = musp ./ (1 - gs);

  D = (1/3)/musp;
  ref = DE(t, r, D, mua, c);

  figure
  hold on
  title(sprintf('Figure 1b: RTE/DE ratio (r = %g m, \\mu_s'' = %g m^{-1})', r, musp))
  set(gca, 'ColorOrder', cool(numel(gs)));
  for idx = 1:numel(gs)
    ratio = RTE_aniso(t, r, mua, musv(idx), gs(idx), c) ./ ref;
    semilogx(t/ps, ratio, 'DisplayName', sprintf('g = %g', gs(idx)))
  end
  line([t(1), t(end)]/ps, [1,1], 'linestyle', '--', 'color', 'k', 'HandleVisibility', 'off')
  xlim([t(1), t(end)]/ps)
  ylim([0.97, 1.03])
  xlabel('time [ps]')
  ylabel('\phi/\psi')
  legend('show')

  %% Figure 1c - RTE/DE ratio for different expressions of D
  mus = 1/mm;
  r   = 10*mm;
  t   = logspace(2,5,301)*ps;

  a  = linspace(0,1,11);
  Ds = (1/3)./(mus + a*mua);

  ref = RTE_iso(t, r, mua, mus, c);

  figure
  hold on
  title(sprintf('Figure 1c: RTE/DE ratio (r = %g m, \\mu_s = %g m^{-1})', r, mus))
  set(gca, 'ColorOrder', cool(numel(a)));
  for idx = 1:numel(Ds)
    ratio = ref ./ DE(t, r, Ds(idx), mua, c);
    semilogx(t/ps, ratio)
  end
  line([t(1), t(end)]/ps, [1,1], 'linestyle', '--', 'color', 'k', 'HandleVisibility', 'off')
  xlim([t(1), t(end)]/ps)
  ylim([0.97, 1.03])
  xlabel('time [ps]')
  ylabel('\phi/\psi')

  %% Figure 1d - Asymptotic limit
  mus   = 1/mm;
  muas  = logspace(-4, 4, 81)/mm;
  a = [0, 0.01, 0.1, 1];

  figure
  hold on
  title('Figure 1d: asymptotic limit of the \phi/\psi ratio')
  set(gca, 'ColorOrder', cool(numel(a)));

  for ia = numel(a):-1:1
    for muaidx = numel(muas):-1:1
      xaxis(ia, muaidx) = muas(muaidx)/mus;
      %% uncomment the two ia below for the explicit limit calculation
      % Ds = (1/3)./(mus + a(ia)*muas(muaidx));
      % curves(ia, muaidx) = find_limit(r, Ds, mus, c, 1e-15);
      curves(ia, muaidx) = exact_limit(a(ia), muas(muaidx), mus);
    end
    plot(xaxis(ia,:), curves(ia,:), 'DisplayName', sprintf('a = %g', a(ia)))
  end

  ylim([1e-4, 1e1])
  set(gca, 'xscale', 'log', 'yscale', 'log')
  xlabel('\mu_a/\mu_s');
  ylabel('asymptotic RTE/DE limit')
  legend('location', 'southwest')
end


function value = find_limit(r, D, mus, v, reltol)
% FIND_LIMIT - Iteratively approach the asymptotic RTE/DE time-domain limit
%
% Syntax:
%   value = find_limit(r, D, mus, v, reltol)
%
% Inputs:
%   r      [1x1 double, length]
%   D      [1x1 double, length]            Diffusion coefficient (DE model)
%   mus    [1x1 double, length^-1]
%   v      [1x1 double, length time^-1]
%   reltol [1x1 double, -]                 Relative tolerance target
%
% Output:
%   value  [1x1 double, -]                 Approximated asymptotic value

  targetreltol = reltol;
  reltol = Inf;
  t = 2*r/v;
  value0 = limit_function(r, t, D, mus, v);
  while reltol > targetreltol
    t = 2*t;
    value = limit_function(r, t, D, mus, v);
    reltol = abs(value - value0)/value;
    value0 = value;
  end
end


function value = limit_function(r, t, D, mus, v)
% LIMIT_FUNCTION - Limiting expression used for RTE/DE asymptotic ratio
%
% Syntax:
%   value = limit_function(r, t, D, mus, v)
%
% Inputs:
%   r   [1x1 double, length]
%   t   [1x1 double, time]
%   D   [1x1 double, length]
%   mus [1x1 double, length^-1]
%   v   [1x1 double, length time^-1]
%
% Output:
%   value [1x1 double, -]

  Dw = D;
  D0 = (1/3)/mus;

  alpha = 4;
  lim = 1/(v*t*D0^(alpha-1))^(1/alpha);
  E = sqrt(D0*v*t) * lim - 1i * r / sqrt(4*D0*v*t);
  try
    % Octave supports complex erf; keep fallback for MATLAB
    cErf = erf(E);
  catch
    % Faddeyeva fallback
    cErf = 1 - Faddeyeva(1i*E, 1e-14)*exp(-E^2);
  end

  value = (Dw/D0)^1.5 * real(cErf) * exp(r^2*(1/Dw-1/D0) / (4*v*t)) - ...
          (Dw/D0) * sin(lim * r) * sqrt(4*pi*Dw*v*t) / (pi*r) * ...
          exp(-D0*v*t*lim^2 + r^2 / (4*Dw*v*t));
end


function value = exact_limit(a, mua, mus)
% EXACT_LIMIT - Analytical asymptotic limit for the RTE/DE ratio
%
% Syntax:
%   value = exact_limit(a, mua, mus)
%
% Inputs:
%   a    [1x1 double, -]            Parameter in D = (1/3)/(mus + a*mua)
%   mua  [1x1 double, length^-1]    Absorption coefficient
%   mus  [1x1 double, length^-1]    Scattering coefficient
%
% Output:
%   value [1x1 double, -]           lim_{t -> inf} (RTE/DE)

  value = (mus/(a*mua + mus))^(3/2);
end

