%load in spatial autocorrelation file

% paths
addpath(genpath('F:\cdrive-samira_nejrgeco'));

% load saline data (all_corr (62x12), all_dist (62x12), and partition)
load('F:\nejrgeco_code_for_publication\figure 2\saline_spatialcorr.mat');

% get the order of runs
clear;
smooth_kernel = 15;
px_to_micron = [1.47 1.47 0.74 0.74 0.74 0.74 0.74 0.74 0.74 0.74 0.74 0.74];
times = [5 5 5 5 5 5 10 5 5 5 5 10];
mice = {'PrL2.1', 'PrL2.1', 'PrL3.2', 'PrL3.2', 'PrL3.2', 'PrL3.2', 'PrL3.2', 'PrL3.4', 'PrL3.4', 'PrL3.4', 'PrL3.4', 'PrL3.4'};
runs = {'2', '4', '1', '2', '3', '4', '5', '1', '2', '3', '4', '5'}
condition = 'saline'


%%
smooth_kernel = 15;
px_to_micron = [1.47 1.47 0.74 0.74 0.74 0.74 0.74 0.74 0.74 0.74 0.74 0.74];
times = [5 5 5 5 5 5 5 5 5 5 10 10];
mice = {'PrL2.1', 'PrL2.1', 'PrL3.2', 'PrL3.2', 'PrL3.2', 'PrL3.2', 'PrL3.4', 'PrL3.4', 'PrL3.4', 'PrL3.4', 'PrL3.4', 'PrL3.4'};
runs = {'1', '2', '1', '2', '3', '4', '1', '2', '3', '4', '5', '6'}
condition = 'desipramine'

for z = 1:12
    mouse = mice{z};
    run = runs{z}
    addpath(genpath(['F:\' mouse]));
    
    % get cell names
    files = dir(fullfile(['F:\' mouse '\' mouse ' GLM predictors\reg2_' mouse '_' condition(1:3) run], '*.csv'));
    names = cell(size(files,1), 1);
    for i = 1:length(names)
        names{i} = [files(i).folder '\' files(i).name];
    end

    % get the low and high synchrony time periods
    T = load(names{1});
    idx = T.mcor_all < mean(T.mcor_all);


end


