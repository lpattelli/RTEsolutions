function out = RTE_iso(r, t, mus, mua, v)
  % RTE_ISO time-dependent fluence predicted by the radiative transfer equation
  %
  % Input:
  % r [scalar, length]:                 source-detector distance
  % t [1D array, time]:                 time
  % mus [scalar, length^-1]:            scattering coefficient
  % mua [scalar, length^-1]:            absorption coefficient
  % v [scalar, length time^-1]:         speed of light
  %
  % Output:
  % out [1D array, time^-1 length^-2]:  fluence

  out = zeros(size(t));
  valid_idx = (t > r/v);

  if ~any(valid_idx), return; end

  t_valid = t(valid_idx);

  f1 = @(k, t) ((v*t-k)./(v*t+k)).^(mus*(r-k)/2) .* ...
               ((pi^3 - 3*pi*(log((v*t-k)./(v*t+k))).^2) .* ...
               cos(mus*pi*(r-k)/2) + ...
               log((v*t-k)./(v*t+k)) .* ...
               (3*pi^2 - (log((v*t-k)./(v*t+k))).^2) .* ...
               sin(mus*pi*(r-k)/2));

  f2 = @(k, t) exp((k.*cot(k/mus) - mus - mua) .* v .* t)./(sin(k/mus)).^2 .* ...
               sin(k*r) .* k.^3;

  % Compute integrals using `arrayfun`
  F1 = arrayfun(@(t) quadgk(@(k) f1(k, t), 0, r, 'AbsTol', 1e-300), t_valid);
  F2 = arrayfun(@(t) quadgk(@(k) f2(k, t), 0, mus*pi/2, 'AbsTol', 1e-300), t_valid);

  % Compute fluence for valid time indices
  out(valid_idx) = (mus^2 ./ (16*pi^2*r*t_valid) .* exp(-(mus+mua) .* v .* t_valid) .* F1 + ...
                    v ./ (2*pi^2*mus^2*r) .* F2 + mus ./ (4*pi*r*t_valid) .* ...
                    exp(-(mus+mua) .* v .* t_valid) .* log((v*t_valid + r) ./ (v*t_valid - r)));
end

