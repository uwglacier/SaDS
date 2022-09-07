% Script to run once per session to set up all the Matlab paths

% If you've installed the 'palettes' package
% (https://github.com/timghill/palettes)
addpath(genpath('~/MATLAB Add-Ons/palettes/'));

% CMOCEAN package
addpath(genpath('~/MATLAB Add-Ons/cmocean/'))

% Meshing
addpath(genpath('meshes/'))
addpath(genpath('mattri/'))

% Add model code directories
addpath(genpath('model/'))

% Add helper function directories
addpath(genpath('functions/'))
