function phi = RTE_aniso_CW(r, mua, mus, g, N)
  % RTE_aniso_CW continuous wave fluence predicted by the radiative transfer
  % equation in an infinite medium with anisotropic scattering
  %
  % Input:
  % r [1D array, length]:               source-detector distance
  % mua [scalar, length^-1]:            absorption coefficient
  % mus [scalar, length^-1]:            scattering coefficient
  % g [scalar]                          scattering anisotropy factor
  % N [integer]                         PN order (default 101)
  %
  % Output:
  % phi [1D array, length^-2]:          fluence

  if nargin < 5
      N = 101;
  end

  sigma = mua + (1 - g .^ (0:N)) * mus;
  beta = (1:N) ./ sqrt((4 * (1:N).^2 - 1) .* sigma(1:N) .* sigma(2:N+1));
  A = diag(beta, 1) + diag(beta, -1);
  [V, D] = eig(A);

  valid_indices = find(diag(D) > 0);
  D = diag(D);
  D = D(valid_indices);
  V = V(1, valid_indices)';

  phi = arrayfun(@(R) sum(exp(-R ./ D) .* (V .^ 2) ./ (D .^ 2)) / (2 * pi * mua * R), r);
end

