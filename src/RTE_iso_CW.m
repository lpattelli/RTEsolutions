function fluence = RTE_iso_CW(r, mua, mus)
% RTE_ISO_CW - CW fluence from the RTE with isotropic scattering
%
% Syntax:
%   fluence = RTE_iso_CW(r, mua, mus)
%
% Inputs:
%   r    [1xR double, length]        Source-detector distance(s)
%   mua  [1x1 double, length^-1]     Absorption coefficient
%   mus  [1x1 double, length^-1]     Scattering coefficient
%
% Outputs:
%   fluence [1xR double, length^-2]  Fluence

  fluence = zeros(size(r));

  mut = mua + mus;
  h = @(x) x - tanh((mut)/mus*x);
  if mua == 0
    L = 0;
  else
    L = fzero(h, sqrt(3*mua/mus));
  end

  for u = 1:numel(fluence)
    if mua == 0
      D = (1/3)/mus;
      T1 = 1/(4*pi*D*r(u));   % Poisson-like solution
    else
      T1 = mut^3*L^2*(1-L^2)/(mut*L^2-mua)*exp(-mut*L*r(u))/(2*pi*mus*r(u));
    end
    f = @(w) exp(-w*r(u))./((1+mus./(2*w).*log(abs((w-mut)./(w+mut)))).^2 + (mus*pi./(2*w)).^2);
    T2 = quadgk(f, mut, inf, 'AbsTol', 1e-30)/(4*pi*r(u));
    fluence(u) = T1 + T2;
  end
end

