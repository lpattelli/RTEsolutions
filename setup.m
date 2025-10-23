% setup.m
%
% Description:
%   Setup script adding necessary folders to the path and (optionally) generates figures
%
% Authors: André Liemert (ILM-ULM)
%          Lorenzo Pattelli (INRIM)
%          Fabrizio Martelli (UNIFI)
%
% License: MIT

% Add directories to the path
addpath(genpath(fullfile(pwd, 'src')));
addpath(genpath(fullfile(pwd, 'utils')));
addpath(genpath(fullfile(pwd, 'data')));
addpath(genpath(fullfile(pwd, 'figures')));

fprintf('Paths added:\n  %s\n  %s\n  %s\n  %s\n', ...
  fullfile(pwd,'src'), fullfile(pwd,'utils'), fullfile(pwd,'data'), fullfile(pwd,'figures'));

% Optional: plot all figures
RUN_FIGURES = true;  % set false to skip plotting when setting up
if RUN_FIGURES
  disp('Generating Figure 1'); Figure1;
  disp('Generating Figure 2'); Figure2;
end

