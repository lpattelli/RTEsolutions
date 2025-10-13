% setup.m
%
% Description:
%   Setup script adding necessary folders to the path
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

% Optional: plot all figures
disp('Generating Figure 1');
Figure1;

disp('Generating Figure 2');
Figure2;

