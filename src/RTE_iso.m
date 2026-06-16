function fluence = RTE_iso(t, r, mua, mus, v)
% RTE_ISO - Time-domain fluence from the RTE with isotropic scattering
%
% Syntax:
%   fluence = RTE_iso(t, r, mua, mus, v)
%
% Inputs:
%   t    [1xT double, time]           Time points
%   r    [1x1 double, length]         Source-detector distance
%   mua  [1x1 double, length^-1]      Absorption coefficient
%   mus  [1x1 double, length^-1]      (Isotropic) scattering coefficient
%   v    [1x1 double, length time^-1] Speed of light in medium
%
% Outputs:
%   fluence [1xT double, time^-1 length^-2]  Fluence

  fluence = zeros(size(t));
  valid_idx = (t > r/v);
  if ~any(valid_idx), return; end

  t_valid = t(valid_idx);

  f1 = @(k, tt) ((v*tt - k)./(v*tt + k)).^(mus*(r - k)/2) .* ...
                ((pi^3 - 3*pi*(log((v*tt - k)./(v*tt + k))).^2) .* ...
                cos(mus*pi*(r - k)/2) + ...
                log((v*tt - k)./(v*tt + k)) .* ...
                (3*pi^2 - (log((v*tt - k)./(v*tt + k))).^2) .* ...
                sin(mus*pi*(r - k)/2));

  f2 = @(k, tt) exp((k.*cot(k/mus) - mus - mua) .* v .* tt)./(sin(k/mus)).^2 .* ...
                sin(k*r) .* k.^3;

  F1 = arrayfun(@(tt) quadgk(@(k) f1(k, tt), 0, r,        'AbsTol', 1e-300), t_valid);
  F2 = arrayfun(@(tt) quadgk(@(k) f2(k, tt), 0, mus*pi/2, 'AbsTol', 1e-300), t_valid);

  fluence(valid_idx) = ( mus^2 ./ (16*pi^2*r*t_valid) .* exp(-(mus+mua) .* v .* t_valid) .* F1 + ...
                         v ./ (  2*pi^2*mus^2*r) .* F2 + ...
                         mus ./ (  4*pi*r*t_valid) .* exp(-(mus+mua).*v.*t_valid) .* ...
                         log((v*t_valid + r) ./ (v*t_valid - r)) );
end

