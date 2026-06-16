function [D_pole, D_res, a_pole, a_res, ratio, d_star, v_star] = cw_aniso_asymptotics(mua, musp, g, N)
% CW_ANISO_ASYMPTOTICS - Dominant PN pole and residue data for HG scattering
%
% Syntax:
%   [D_pole, D_res, a_pole, a_res, ratio] = cw_aniso_asymptotics(mua, musp, g)
%   [D_pole, D_res, a_pole, a_res, ratio, d_star, v_star] = cw_aniso_asymptotics(mua, musp, g, N)
%
% Description:
%   Computes the scalar diffusion coefficients obtained by matching either the
%   dominant CW RTE pole or its leading residue. The transport calculation uses
%   the same PN matrix as RTE_aniso_CW for Henyey-Greenstein anisotropy.
%
%   The pole-matched coefficient is
%      D_pole = mua*d_star^2,
%   where d_star is the largest positive PN attenuation length. The leading
%   RTE residue gives the asymptotic pole-matched RTE/DE ratio
%      ratio = 2*v_star^2,
%   where v_star is the first component of the normalized eigenvector. The
%   residue-matched diffusion coefficient is D_res = D_pole/ratio.
%
% Inputs:
%   mua  [array double, length^-1] Absorption coefficient(s)
%   musp [array double, length^-1] Reduced scattering coefficient(s)
%   g    [array double, -]         HG anisotropy factor(s)
%   N    [1x1 int]                 PN order (default 101)
%
% Outputs:
%   D_pole [array double, length] Dominant-pole matched diffusion coefficient
%   D_res  [array double, length] Dominant-residue matched diffusion coefficient
%   a_pole [array double, -]      D_pole mapped to D=1/[3(mus'+a*mua)]
%   a_res  [array double, -]      D_res mapped to D=1/[3(mus'+a*mua)]
%   ratio  [array double, -]      lim_{r->inf} Phi_RTE/Phi_DE at D=D_pole
%   d_star [array double, length] Dominant PN attenuation length
%   v_star [array double, -]      First eigenvector component of dominant mode

  if nargin < 4
    N = 101;
  end

  [mua, musp, g] = common_size(mua, musp, g);

  D_pole = zeros(size(mua));
  D_res  = zeros(size(mua));
  a_pole = zeros(size(mua));
  a_res  = zeros(size(mua));
  ratio  = zeros(size(mua));
  d_star = zeros(size(mua));
  v_star = zeros(size(mua));

  for idx = 1:numel(mua)
    [D_pole(idx), D_res(idx), a_pole(idx), a_res(idx), ratio(idx), d_star(idx), v_star(idx)] = ...
      scalar_cw_aniso_asymptotics(mua(idx), musp(idx), g(idx), N);
  end
end


function [D_pole, D_res, a_pole, a_res, ratio, d_star, v_star] = scalar_cw_aniso_asymptotics(mua, musp, g, N)
  if mua == 0
    [a0] = HG_high_albedo_coeffs(g);
    D_pole = 1 / (3*musp);
    D_res  = D_pole;
    a_pole = a0;
    a_res  = a0;
    ratio  = 1;
    d_star = Inf;
    v_star = 1 / sqrt(2);
    return
  end

  denom = max(1 - g, 1e-12);
  mus = musp / denom;

  sigma = mua + (1 - g .^ (0:N)) * mus;
  beta  = (1:N) ./ sqrt((4*(1:N).^2 - 1) .* sigma(1:N) .* sigma(2:N+1));
  A     = diag(beta, 1) + diag(beta, -1);

  [V, Lambda] = eig(A);
  d = diag(Lambda);
  valid = find(d > 0);

  if isempty(valid)
    D_pole = NaN;
    D_res  = NaN;
    a_pole = NaN;
    a_res  = NaN;
    ratio  = NaN;
    d_star = NaN;
    v_star = NaN;
    return
  end

  d = d(valid);
  v0 = V(1, valid).';
  [d_star, imax] = max(d);
  v_star = v0(imax);

  D_pole = mua * d_star^2;
  ratio  = 2 * v_star^2;
  D_res  = D_pole / ratio;

  a_pole = (1/(3*D_pole) - musp) / mua;
  a_res  = (1/(3*D_res)  - musp) / mua;
end


function varargout = common_size(varargin)
  sz = [];
  for k = 1:nargin
    if ~isscalar(varargin{k})
      sz = size(varargin{k});
      break
    end
  end

  if isempty(sz)
    sz = [1, 1];
  end

  varargout = cell(1, nargin);
  for k = 1:nargin
    if isscalar(varargin{k})
      varargout{k} = varargin{k} + zeros(sz);
    elseif isequal(size(varargin{k}), sz)
      varargout{k} = varargin{k};
    else
      error('cw_aniso_asymptotics: inputs must be scalar or have matching sizes');
    end
  end
end
