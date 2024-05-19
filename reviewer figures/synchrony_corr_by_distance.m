clear;
addpath(genpath('H:\nejrgeco_code_for_publication\'));

% set global variables
distance_steps = 1:200;
all_corr = [];



% run identifiers
runs = [{'PrL2.1' 'sal' 'saline' '2' '5min'};...
        {'PrL2.1' 'sal' 'saline' '4' '5min'};...
        {'PrL3.2' 'sal' 'saline' '1' '5min'};...
        {'PrL3.2' 'sal' 'saline' '2' '5min'};...
        {'PrL3.2' 'sal' 'saline' '3' '5min'};...
        {'PrL3.2' 'sal' 'saline' '4' '5min'};...
        {'PrL3.2' 'sal' 'saline' '5' '10min'};...
        {'PrL3.4' 'sal' 'saline' '1' '5min'};...
        {'PrL3.4' 'sal' 'saline' '2' '5min'};...
        {'PrL3.4' 'sal' 'saline' '3' '5min'};...
        {'PrL3.4' 'sal' 'saline' '4' '5min'};...
        {'PrL3.4' 'sal' 'saline' '5' '10min'}];

       
for j = 1:12

    % get the run identifier
    run = runs(j, :)
    
    % get all centroids for the run
    addpath(genpath(['H:\' run{1}])); 
    filename = [run{1} '_roi_struct_' run{3} '_run' run{4} '_' run{5} '.mat'];
    load(filename)
    channel_origin = {roi_data.channel_origin};
    roi_data2 = roi_data(strcmp(channel_origin, 'jrgeco'));
    center = {roi_data2.centroid};
    center = cell2mat(center');
    
    % get all synchrony time courses for the run
    D = ['H:\*' run{1} '*\*GLM predictors*\*reg2*' run{2} run{4}];
    files = dir(fullfile(D, '*.csv'));
    folder = {files.folder};
    names = natsortfiles({files.name});
    names = cellfun(@(x,y)[x '\' y], folder, names, 'UniformOutput', false);
    mcors = zeros(4200, length(names));
    for i = 1:length(names)
        T = readtable(names{i});
        mcors(:,i) = T.mcor_bglg(1:4200);
    end
    
    % get distance and correlation matrix, linearize, and bin
    mcor_corr = corrcoef(mcors); mcor_corr = triu(mcor_corr, 1); mcor_corr = mcor_corr(find(mcor_corr)); %mcor_corr = reshape(triu(mcor_corr), [], 1);
    mcor_dist = squareform(pdist(center, 'Euclidean')); mcor_dist = triu(mcor_dist, 1); mcor_dist = round(mcor_dist(find(mcor_dist))); %mcor_dist = reshape(triu(mcor_dist), [], 1); mcor_dist = round(mcor_dist);
    mcor_corr_steps = zeros(size(distance_steps));
    for i = 1:length(distance_steps)
        mcor_corr_steps(i) = mean(mcor_corr(find(mcor_dist == distance_steps(i))));
    end
    all_corr = [all_corr; mcor_corr_steps];
end

%% plot

% determine downsample factor
x1 = distance_steps;
dwnsmp = 10;

% binning the correlations
new_corr = mat2cell(all_corr, 12, repmat(dwnsmp, [1, length(x1)/dwnsmp]));
new_corr = cellfun(@(x)reshape(x,[],1), new_corr, 'UniformOutput', false);
new_corr_cell = new_corr;
new_corr = cell2mat(new_corr);

% binning the x distances
new_x = mat2cell(x1, 1, repmat(dwnsmp, [1, length(x1)/dwnsmp]));
new_x = cellfun(@mean, new_x, 'UniformOutput', true);

% remove final two values because they are empty
new_corr = new_corr(:, 1:18);
new_x = new_x(1:18);
new_corr_cell = new_corr_cell(1:18);

% Get the 95% CI
new_corr_nonan = cellfun(@(y)y(~isnan(y)), new_corr_cell, ...
    'UniformOutput', false);
bci = [];
for i = 1:length(new_corr_nonan)
    bci_temp =  bootci(120, @mean, new_corr_nonan{i});
    bci = [bci bci_temp];
end

% plot
figure;
plot(new_x, nanmean(new_corr), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
hold on;
patch([new_x fliplr(new_x)], [bci(1, :) fliplr(bci(2, :))], 'k', 'FaceAlpha',0.2, 'EdgeColor','none');

% figure properties
ylim([0, 1]); yticks([0, 0.5, 1]);
xlim([5, 150]); xticks([0, 25, 50, 75, 100, 125, 150]);
set(gca, 'TickDir', 'out');
title('Synchrony Similarity by Distance');
ylabel('Synchrony Similarity'); xlabel('Distance');
axis square;
box off;



