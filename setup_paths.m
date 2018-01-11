% Set up environments for Allan Hills scripts.
%
% Laura Kehrl, UW, 3/14/2017

CODE_HOME = '~/Dropbox/AllanHills2016/Code/';
REPO_HOME = fullfile(CODE_HOME,'allanhills/');
DATA_HOME = '~/Data/';

% Add other toolboxes
addpath(fullfile(CODE_HOME,'spicker'));
addpath(fullfile(CODE_HOME,'export_fig'));
addpath(fullfile(CODE_HOME,'cubehelix'));
addpath(fullfile(CODE_HOME,'AntarcticMappingTools'));
addpath(fullfile(CODE_HOME,'bedmap2_toolbox_v4.6'));
addpath(fullfile(CODE_HOME,'measures'));

% Add paths for this repository
addpath(genpath(REPO_HOME));

% Add data paths
addpath(fullfile(DATA_HOME,'Velocity/AntarcticaMeasures'));
addpath(fullfile(DATA_HOME,'Bed/BedMap2'));
addpath(fullfile(CODE_HOME,'Antarctic Basins'));

% Add paths for radar data processed by Howard Conway
addpath('~/Dropbox/Allan Hills');
addpath('~/Dropbox/Allan Hills/BIT58');
