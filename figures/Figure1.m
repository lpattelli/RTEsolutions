function Figure1
  %% Asymptotic analysis of the RTE and DE solutions in the time-domain
  %
  % Authors: André Liemert (ILM-ULM)
  %          Lorenzo Pattelli (INRIM)
  %          Fabrizio Martelli (UNIFI)
  %
  % License: MIT

  % define units and constants
  mm = 1e-3;
  ps = 1e-12;
  c = 0.299792458*mm/ps;

  %% Figure 1a: RTE/DE ratio at different distances from the source

  % define optical parameters
  mua = 0.01/mm;
  mus = 1/mm;
  D = (1/3)/mus;  % isotropic case

  % define time axis and distance values
  t = logspace(1,5,1001)*ps;
  rs = logspace(-1, 2, 10)*mm;
  data = nan(numel(rs), numel(t));

  for idx = 1:numel(rs)
    postballistic = t > rs(idx)/c;
    data(idx,postballistic) = RTE_iso(rs(idx), t(postballistic), mus, mua, c) ./ ...
                              DE(D, rs(idx), t(postballistic), mua, c);
    starttime = find(diff(data(idx,:) < 0.01), 1, 'last');
    data(idx, 1:starttime) = NaN;
  end

  figure
  hold on
  title('Figure 1a: asymptotic RTE/DE ratio convergence as a function of distance')
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


  %% Figure 1b: validity of the similarity relations in the anisotropic scattering case
  clearvars -except mm ps c mua

  % set illustrative distance value
  musp = 1/mm;

  r = 10*mm;
  t = logspace(2,5,301)*ps;

  gs = linspace(0., 0.9, 10);
  mus = musp ./ (1 - gs);

  D = (1/3)/musp;

  ref = DE(D, r, t, mua, c);

  figure
  hold on
  title(sprintf('Figure 1b: asymptotic RTE/DE ratio convergence (r = %g m, µ_s'' = %g m^{-1})', r, musp))
  set(gca, 'ColorOrder', cool(numel(gs)));
  ratio = zeros(numel(gs), numel(t));
  for idx = 1:numel(gs)
    ratio(idx,:) = RTE_aniso(t, r, mus(idx), mua, gs(idx), c) ./ ref;
    semilogx(t/ps, ratio(idx,:), 'DisplayName', sprintf('g = %g', gs(idx)))
  end
  line([t(1), t(end)]/ps, [1,1], 'linestyle', '--', 'color', 'k', 'HandleVisibility', 'off')
  xlim([t(1), t(end)]/ps)
  ylim([0.97, 1.03])
  xlabel('time [ps]')
  ylabel('\phi/\psi')


  %% Figure 1c: RTE/DE ratio using different expressions for the diffusion coefficient
  clearvars -except mm ps c mua

  % set illustrative distance value
  mus = 1/mm;

  r = 10*mm;
  t = logspace(2,5,301)*ps;

  % different D values based on the "a" parameter in [0, 1]
  a = linspace(0,1,11);
  Ds = (1/3)./(mus + a*mua);

  % reference solution using the asymptotic RTE solution
  ref = RTE_iso(r, t, mus, mua, c);

  figure
  hold on
  title(sprintf('Figure 1c: asymptotic RTE/DE ratio convergence (r = %g m, µ_s = %g m^{-1})', r, mus))
  set(gca, 'ColorOrder', cool(numel(a)));
  ratio = zeros(numel(Ds), numel(t));
  for idx = 1:numel(Ds)
    ratio(idx,:) = ref ./ DE(Ds(idx), r, t, mua, c);
    semilogx(t/ps, ratio(idx,:))
  end
  line([t(1), t(end)]/ps, [1,1], 'linestyle', '--', 'color', 'k', 'HandleVisibility', 'off')
  xlim([t(1), t(end)]/ps)
  ylim([0.97, 1.03])
  xlabel('time [ps]')
  ylabel('\phi/\psi')


  %% Figure 1d: asymptotic limit
  clearvars -except mm ps c r

  mus = 1/mm;
  muas = logspace(-4, 4, 81)/mm;

  a = [0, 0.01, 0.1, 1];

  figure
  hold on
  title('Figure 1d: asymptotic limit of the \phi/\psi ratio')
  set(gca, 'ColorOrder', cool(numel(a)));

  xaxis = zeros(numel(a), numel(muas));
  curves = zeros(numel(a), numel(muas));

  for Didx = 1:numel(a)
    for muaidx = 1:numel(muas)
      xaxis(Didx, muaidx) = muas(muaidx)/mus;
      %% uncomment the two lines below for the explicit limit calculation
      % Ds = (1/3)./(mus + a(Didx)*muas(muaidx));
      % curves(Didx, muaidx) = find_limit(r, Ds, mus, c, 1e-15);
      curves(Didx, muaidx) = exact_limit(a(Didx), mus, muas(muaidx));
    end
      plot(xaxis(Didx,:), curves(Didx,:), 'DisplayName', ['a =', num2str(a(Didx))])
  end

  ylim([1e-4, 1e1])
  set(gca, 'xscale', 'log', 'yscale', 'log')
  xlabel('µ_a/µ_s');
  ylabel('asymptotic RTE/DE limit')
  legend('location', 'southwest')
end

function value = find_limit(r, D, mus, v, reltol)
  % FIND_LIMIT find the asymptotic limit within target reltol, iteratively
  %
  % Input:
  % r [scalar, length]:             source-detector distance
  % D [scalar, length]:             diffusion coefficient
  % v [scalar, length time^-1]:     speed of light
  % reltol [scalar, -]:             relative tolerance for the limit convergence
  %
  % Output:
  % value [scalar, -]:              asymptotic value

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
  % LIMIT_FUNCTION limiting expression for the RTE/DE limit
  %
  % Input:
  % r [scalar, length]:             source-detector distance
  % t [scalar, time]:               time
  % D [scalar, length]:             diffusion coefficient
  % v [scalar, length time^-1]:     speed of light
  %
  % Output:
  % value [scalar, -]:              RTE/DE ratio calculated at time t

  Dw = D;
  D = (1/3)/mus;

  alpha = 4;
  lim = 1/(v*t*D^(alpha-1))^(1/alpha);
  E = sqrt(D*v*t) * lim - 1i * r / sqrt(4*D*v*t);
  try % Octave
    cErf = erf(E);
  catch % MATLAB (requires Faddeyeva function)
    cErf = 1 - Faddeyeva(1i*E, 1e-14)*exp(-E^2);
  end
  value = (Dw/D)^1.5 * real(cErf) * exp(r^2*(1/Dw-1/D) / (4*v*t)) - ...
          (Dw/D) * sin(lim * r) * sqrt(4*pi*Dw*v*t) / (pi*r) * ...
          exp(-D*v*t*lim^2 + r^2 / (4*Dw*v*t));
end

function value = exact_limit(a, mus, mua)
  % EXACT_LIMIT simple analytical equation for the asymptotic value
  %
  % Input:
  % a [scalar, ]:                 factor appearing in the absorption-dependent D
  % mus [scalar, length^-1]:      scattering coefficient
  % mua [scalar, length^-1]:      absorption coefficient
  %
  % Output:
  % value [scalar, -]:            lim t->oo RTE/DE ratio

  value = (mus/(a*mua + mus))^(3/2);
end

