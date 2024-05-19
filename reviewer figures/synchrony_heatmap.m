% heatmaps NE synchrony
clear;
D = 'H:\PrL*\*GLM predictors*\*reg2_*sal*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('H:\PrL2.1')); addpath(genpath('H:\PrL3.2')); addpath(genpath('H:\PrL3.4')); addpath(genpath('H:\samira_nejrgeco'));
all_names1 = cell(size(files,1), 1);
all_synchrony = zeros(length(all_names1), 4200);
for i = 1:length(all_names1)
    all_names1{i} = [files(i).folder '\' files(i).name];
    T = readtable(all_names1{i});
    mcor_all = T.mcor_all; mcor_all = mcor_all';
    all_synchrony(i, :) = mcor_all(1:4200);
end

% plot
figure;
imagesc(all_synchrony);
colorbar;
title('NE Synchrony for All Runs');