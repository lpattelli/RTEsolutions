% setup.m
%
% Description:
%   Add RTEsolutions source folders to the Octave/MATLAB path.
%
% Authors: Andre Liemert (ILM-ULM)
%          Lorenzo Pattelli (INRIM)
%          Fabrizio Martelli (UNIFI)
%
% License: MIT

rootDir = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(rootDir, 'src')));
addpath(genpath(fullfile(rootDir, 'utils')));
addpath(genpath(fullfile(rootDir, 'figures')));

fprintf('RTEsolutions paths added:\n');
fprintf('  %s\n', fullfile(rootDir, 'src'));
fprintf('  %s\n', fullfile(rootDir, 'utils'));
fprintf('  %s\n', fullfile(rootDir, 'figures'));

% Set this to true to generate all production figure windows after setup.
RUN_FIGURES = false;

if RUN_FIGURES
  Figure1;
  Figure2;
  Figure3;
  Figure4;
end
