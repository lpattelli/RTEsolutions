function out = DE(D, r, t, mua, v)
  % DE time-dependent fluence predicted by the diffusion equation
  %
  % Input:
  % D [scalar, length]:                 diffusion coefficient
  % r [scalar, length]:                 source-detector distance
  % t [1D array, time]:                 time
  % mua [scalar, length^-1]:            absorption coefficient
  % v [scalar, length time^-1]:         speed of light
  %
  % Output:
  % out [1D array, time^-1 length^-2]:  fluence

  out = v ./ (4*pi*D*v*t).^1.5 .* exp(-mua*v*t) .* exp(-r^2./(4*D*v*t));
end

