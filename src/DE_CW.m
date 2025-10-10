function out = DE_CW(D, r, mua)
  % DE_cw continuous wave fluence predicted by the diffusion equation
  %
  % Input:
  % D [scalar, length]:                 diffusion coefficient
  % r [1D array, length]:               source-detector distance
  % mua [scalar, length^-1]:            absorption coefficient
  %
  % Output:
  % out [1D array, length^-2]:          fluence

  out = exp(-sqrt(mua/D)*r)./(4*pi*D*r);
end

