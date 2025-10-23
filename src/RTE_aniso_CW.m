function fluence = RTE_aniso_CW(r, mua, mus, g, N)
% RTE_ANISO_CW - CW fluence via PN approximation with HG anisotropy
%
% Syntax:
%   fluence = RTE_aniso_CW(r, mua, mus, g)
%   fluence = RTE_aniso_CW(r, mua, mus, g, N)
%
% Inputs:
%   r    [1xR double, length]        Source-detector distance(s)
%   mua  [1x1 double, length^-1]     Absorption coefficient
%   mus  [1x1 double, length^-1]     Scattering coefficient
%   g    [1x1 double, -]             HG anisotropy factor
%   N    [1x1 int]                   PN order (default 101)
%
% Outputs:
%   fluence [1xR double, length^-2]  Fluence

  if nargin < 5
    N = 101;
  end

  sigma = mua + (1 - g .^ (0:N)) * mus;   % (1+N)
  beta  = (1:N) ./ sqrt((4*(1:N).^2 - 1) .* sigma(1:N) .* sigma(2:N+1));
  A     = diag(beta, 1) + diag(beta, -1);

  [V, D] = eig(A);
  valid_idx = find(diag(D) > 0);
  D = diag(D);
  D = D(valid_idx);
  V = V(1, valid_idx)';

  % vectorized over r via arrayfun
  fluence = arrayfun(@(R) sum(exp(-R ./ D) .* (V .^ 2) ./ (D .^ 2)) / (2 * pi * mua * R), r);
end

