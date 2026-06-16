function fluence = DE_CW(r, D, mua)
% DE_CW - Continuous-wave fluence from the diffusion equation in 3D
%
% Syntax:
%   fluence = DE_CW(r, D, mua)
%
% Inputs:
%   r   [1xR double, length]         Source-detector distance(s)
%   D   [1x1 double, length]         Diffusion coefficient
%   mua [1x1 double, length^-1]      Absorption coefficient
%
% Outputs:
%   fluence [1xR double, length^-2]  Fluence

  fluence = exp(-sqrt(mua/D).*r) ./ (4*pi*D.*r);
end

