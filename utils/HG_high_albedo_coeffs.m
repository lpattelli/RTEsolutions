function [a0, a1, ratio_slope] = HG_high_albedo_coeffs(g)
% HG_HIGH_ALBEDO_COEFFS - Leading HG high-albedo CW asymptotic coefficients
%
% Syntax:
%   [a0, a1, ratio_slope] = HG_high_albedo_coeffs(g)
%
% Description:
%   Coefficients for eta = mua/mus' in the high-albedo expansions
%      a_star(g,eta) = a0(g) + a1(g)*eta + O(eta^2)
%      Phi_RTE/Phi_DE = 1 - ratio_slope(g)*eta + O(eta^2)
%   where D_star is the exponent-matched CW diffusion coefficient.
%
% Inputs:
%   g [array double, -] HG anisotropy factor
%
% Outputs:
%   a0          [array double, -] zero-absorption exponent-matched a
%   a1          [array double, -] first absorption correction for a_star
%   ratio_slope [array double, -] leading slope of the residue mismatch

  a0 = 1 - 4 ./ (5 .* (1 + g));
  a1 = -4 .* (35.*g.^3 + 7.*g.^2 + 7.*g - 1) ./ ...
       (175 .* (1 + g).^2 .* (1 + g + g.^2));
  ratio_slope = 4 ./ (5 .* (1 + g));
end
