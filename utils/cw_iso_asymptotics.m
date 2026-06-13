function [lambda, a_star, ratio, D_star, A1, nu1] = cw_iso_asymptotics(mua, mus)
% CW_ISO_ASYMPTOTICS - Exact isotropic CW pole and residue data
%
% Syntax:
%   [lambda, a_star, ratio, D_star, A1, nu1] = cw_iso_asymptotics(mua, mus)
%
% Description:
%   Computes the dominant discrete-mode attenuation and asymptotic mismatch
%   for the isotropic infinite-medium CW RTE Green's function.
%
% Inputs:
%   mua [array double, length^-1] Absorption coefficient(s)
%   mus [array double, length^-1] Scattering coefficient(s), scalar or same size as mua
%
% Outputs:
%   lambda [array double, -]        Root of atanh(lambda)/lambda = 1 + mua/mus
%   a_star [array double, -]        Exponent-matched parameter in D = 1/[3(mus+a*mua)]
%   ratio  [array double, -]        lim_{r->inf} Phi_RTE/Phi_DE at D = D_star
%   D_star [array double, length]   Exponent-matched diffusion coefficient
%   A1     [array double, length^-1] Dominant RTE mode prefactor
%   nu1    [array double, length^-1] Dominant RTE attenuation constant

  if isscalar(mus)
    mus = mus + zeros(size(mua));
  end
  if isscalar(mua)
    mua = mua + zeros(size(mus));
  end

  lambda = arrayfun(@(ma, ms) solve_lambda(ma / ms), mua, mus);
  mut = mua + mus;
  nu1 = mut .* lambda;

  D_star = zeros(size(mua));
  A1     = zeros(size(mua));
  a_star = zeros(size(mua));
  ratio  = zeros(size(mua));

  nonzero = mua > 0;

  D_star(nonzero) = mua(nonzero) ./ (nu1(nonzero) .^ 2);
  A1(nonzero) = 2 .* mut(nonzero).^3 .* lambda(nonzero).^2 .* (1 - lambda(nonzero).^2) ./ ...
                (mus(nonzero) .* (mut(nonzero).*lambda(nonzero).^2 - mua(nonzero)));
  ratio(nonzero) = A1(nonzero) .* D_star(nonzero);
  a_star(nonzero) = (1 ./ (3 .* D_star(nonzero)) - mus(nonzero)) ./ mua(nonzero);

  % limiting values for the nonabsorbing case
  zero = ~nonzero;
  D_star(zero) = 1 ./ (3 .* mus(zero));
  A1(zero)     = 1 ./ D_star(zero);
  ratio(zero)  = 1;
  a_star(zero) = 1/5;
end


function lambda = solve_lambda(rho)
% SOLVE_LAMBDA - Root of atanh(lambda)/lambda = 1 + rho on (0,1)

  if rho <= 0
    lambda = 0;
    return;
  end

  f = @(L) atanh(L)./L - (1 + rho);

  % Bracketed solve is robust from high albedo to moderately absorbing cases.
  lo = sqrt(eps);
  hi = 1 - sqrt(eps);
  lambda = fzero(f, [lo, hi]);
end
