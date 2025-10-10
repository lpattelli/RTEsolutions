function out = RTE_aniso(t, r, mus, mua, g, v, N, Nk, R_s_mult)
% RTE_aniso
%   Compute the time-dependent fluence in an infinite homogeneous medium
%   using the PN approximation to the radiative transfer equation with
%   Henyey–Greenstein anisotropic scattering.
%
% Inputs
%   t           [1×T double]   Time array at which to evaluate fluence
%   r           [1×1 double]   Source–detector separation
%   mus         [1×1 double]   Scattering coefficient
%   mua         [1×1 double]   Absorption coefficient
%   g           [1×1 double]   HG asymmetry factor
%   v           [1×1 double]   Speed of light in the medium
%   N           [1×1 int]      PN expansion order (optional; default = 151)
%   Nk          [1×1 int]      Number of wave numbers (optional; default = 500)
%   R_s_mult    [1×1 double]   Sphere radius multiplier (optional; default = 100)
%
% Output
%   out         [1×T double]   Time-dependent fluence

  % set default arguments
  if nargin < 9, R_s_mult = 100; end
  if nargin < 8, Nk = 500; end
  if nargin < 7, N = 151; end

  % Henyey–Greenstein phase function
  R_s   = R_s_mult * r;
  l     = (0:N).';
  sigma = mua + (1 - g.^l) * mus;       % (N+1×1)

  for k = Nk:-1:1
    ek          = (k * pi)/R_s;
    od          = 1i*ek*(1:N)./sqrt((1:2:2*N-1).*(3:2:2*N+1));
    A           = diag(od, +1) + diag(sigma) + diag(od, -1);
    [V,D]       = eig(A);
    b           = V \ eye(N+1,1);
    eig_k(:,k)  = diag(D);               % N+1 × Nk
    Res(:,k)    = ek * sin(r*ek) * (V(1,:).*b.').';  % N+1 × Nk
  end

  ev = eig_k(:);     % (N+1)*Nk × 1
  rv = Res(:);       % same size

  E     = exp(- (ev*v) * t(:)' );        % (N+1)*Nk × T
  sumE  = real(rv.' * E);                % 1×T
  out   = sumE * v/(pi*r*R_s) / 2;       % 1×T

  % zero‐out pre-ballistic times
  out(t < r/v) = 0;
end

