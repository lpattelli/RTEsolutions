function [a_holte, ratio_holte] = holte_hg_asymptotics(mua, musp, g)
% HOLTE_HG_ASYMPTOTICS - Holte/Aronson HG expansion for CW exponent matching
%
% Syntax:
%   [a_holte, ratio_holte] = holte_hg_asymptotics(mua, musp, g)
%
% Description:
%   Evaluates the Holte expansion as quoted by Aronson and Corngold for a
%   Henyey-Greenstein phase function. This helper is retained for comparison
%   with the direct dominant-pole calculation. It is an asymptotic expansion
%   for the exponent-matched branch and is not the default production route.
%
% Inputs:
%   mua  [array double, length^-1] Absorption coefficient(s)
%   musp [array double, length^-1] Reduced scattering coefficient(s)
%   g    [array double, -]         HG anisotropy factor(s)
%
% Outputs:
%   a_holte    [array double, -] Holte exponent-matched parameter
%   ratio_holte [array double, -] Associated asymptotic ratio estimate

  [mua, musp, g] = common_size(mua, musp, g);

  a_holte = zeros(size(mua));
  ratio_holte = zeros(size(mua));

  for idx = 1:numel(mua)
    [a_holte(idx), ratio_holte(idx)] = scalar_holte_hg(mua(idx), musp(idx), g(idx));
  end
end


function [a, ratio] = scalar_holte_hg(mua, musp, g)
  if mua == 0
    [a] = HG_high_albedo_coeffs(g);
    ratio = 1;
    return
  end

  mu_tr = musp + mua;

  denom = max(1 - g, 1e-12);
  mus = musp / denom;
  mu_t = mus + mua;
  omega = mus / mu_t;

  f1 = g;
  f2 = g^2;
  f3 = g^3;
  f4 = g^4;

  h0 = 1 * (1 - omega);
  h1 = 3 * (1 - omega*f1);
  h2 = max(5 * (1 - omega*f2), 1e-15);
  h3 = max(7 * (1 - omega*f3), 1e-15);
  h4 = max(9 * (1 - omega*f4), 1e-15);

  x = 4*h0/h2;
  B = 1 ...
      - x ...
      + x^2 * (1 - (9*h1)/(4*h3)) ...
      - x^3 * (1 - (27*h1)/(4*h3) + (81/16)*(h1^2/(16*h3^2)) + 9*(h1^2/h3^2)*(h2/h4));

  a = (mu_tr*B - musp) / mua;
  ratio = 2*(mu_tr - 3*mua*B) / (musp*(3*B - 1));

  if ~isfinite(a)
    a = NaN;
  end
  if ~isfinite(ratio)
    ratio = NaN;
  else
    ratio = max(0, min(1, ratio));
  end
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
      error('holte_hg_asymptotics: inputs must be scalar or have matching sizes');
    end
  end
end
