function fluence = DE(t, r, D, mua, v)
% DE - Time-domain fluence from the diffusion equation in an infinite medium
%
% Syntax:
%   fluence = DE(t, r, D, mua, v)
%
% Description:
%   Closed-form Green's function of the diffusion equation in 3D.
%
% Inputs:
%   t   [1xT double, time]           Time points
%   r   [1x1 double, length]         Source-detector distance
%   D   [1x1 double, length]         Diffusion coefficient
%   mua [1x1 double, length^-1]      Absorption coefficient
%   v   [1x1 double, length time^-1] Speed of light in medium
%
% Outputs:
%   fluence [1xT double, time^-1 length^-2]  Fluence

  fluence = v ./ (4*pi*D*v*t).^1.5 .* exp(-mua*v*t) .* exp(-r^2./(4*D*v*t));
end

