function fluence = RTE_aniso(t, r, mua, mus, g, v, N, Nk, R_s_mult)
% RTE_ANISO - Time-domain fluence via PN approximation with HG anisotropy
%
% Syntax:
%   fluence = RTE_aniso(t, r, mua, mus, g, v)
%   fluence = RTE_aniso(t, r, mua, mus, g, v, N, Nk, R_s_mult)
%
% Inputs:
%   t         [1xT double, time]           Time points
%   r         [1x1 double, length]         Source-detector distance
%   mua       [1x1 double, length^-1]      Absorption coefficient
%   mus       [1x1 double, length^-1]      Scattering coefficient
%   g         [1x1 double, -]              HG anisotropy factor
%   v         [1x1 double, length time^-1] Speed of light in medium
%   N         [1x1 int]                    PN order (default 151)
%   Nk        [1x1 int]                    Number of wave numbers (default 500)
%   R_s_mult  [1x1 double]                 Sphere radius multiplier (default 100)
%
% Outputs:
%   fluence   [1xT double, time^-1 length^-2]  Fluence

  if nargin < 9, R_s_mult = 100; end
  if nargin < 8, Nk       = 500; end
  if nargin < 7, N        = 151; end

  R_s   = R_s_mult * r;
  l     = (0:N).';
  sigma = mua + (1 - g.^l) * mus;       % (N+1x1)
  eig_k = zeros(N+1, Nk);
  Res   = zeros(N+1, Nk);

  for k = Nk:-1:1
    ek = (k*pi) / R_s;
    od = 1i*ek*(1:N) ./ sqrt((1:2:2*N-1) .* (3:2:2*N+1));
    A  = diag(od, +1) + diag(sigma) + diag(od, -1);

    [V, Dm] = eig(A);
    b = V \ eye(N+1, 1);
    eig_k(:, k) = diag(Dm);                           % N+1 x Nk
    Res(:, k) = ek * sin(r*ek) * (V(1,:) .* b.').';   % N+1 x Nk
  end

  ev = eig_k(:);     % (N+1)*Nk x 1
  rv = Res(:);       % (N+1)*Nk x 1

  E       = exp(-(ev*v) * t(:)');   % (N+1)*Nk x T
  sumE    = real(rv.' * E);         % 1 x T
  fluence = sumE * v/(pi*r*R_s) / 2;

  % zero out pre-ballistic times
  fluence(t < r/v) = 0;
end
