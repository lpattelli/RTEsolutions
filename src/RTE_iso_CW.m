function out = RTE_iso_CW(r, mus, mua)
  % RTE_iso_CW continuous wave fluence predicted by the radiative transfer
  % equation for the case of isotropic scattering
  %
  % Input:
  % r [scalar, length]:                 source-detector distance
  % mus [scalar, length^-1]:            scattering coefficient
  % mua [scalar, length^-1]:            absorption coefficient
  %
  % Output:
  % out [1D array, length^-2]:          fluence

  out = zeros(size(r));

  mut = mua + mus;
  h = @(x) x-tanh((mut)/mus*x);
  if mua == 0
    L = 0;
  else
    L = fzero(h, sqrt(3*mua/mus));
  end

  for u = 1:numel(out) % the loop could be removed using arrayfun?
    if mua == 0
      T1 = 1/(4*pi*D1*r(u));   % Poisson-like solution
    else
      T1 = mut^3*L^2*(1-L^2)/(mut*L^2-mua)*exp(-mut*L*r(u))/(2*pi*mus*r(u));
    end
    f = @(w) exp(-w*r(u))./((1+mus./(2*w).*log(abs((w-mut)./(w+mut)))).^2+(mus*pi./(2*w)).^2);
    T2 = quadgk(f, mut, inf, 'AbsTol', 1e-30)/(4*pi*r(u));
    out(u) = T1 + T2;
  end
end

