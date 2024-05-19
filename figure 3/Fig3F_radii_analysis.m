clear; 

% load GLM files ( single radius )
D = 'F:\*PrL*\*GLM predictors*\*reg2*sal*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('F:\PrL2.1')); addpath(genpath('F:\PrL3.2')); addpath(genpath('F:\PrL3.4')); addpath('F:\nejrgeco_code_for_publication');
names = cell(size(files,1), 1);
for i = 1:length(names)
    names{i} = [files(i).folder '\' files(i).name];
end

% load local and bulk green for different radii
D = 'F:\radii\*PrL*\*GLM predictors*\*reg*sal*'; 
files = dir(fullfile(D, '*.mat'));
addpath(genpath('F:\radii')); 
names2 = {files.name}';

% sanity check for data loading
% generate unique IDs for GLM files
[C,~] = cellfun(@(x) strsplit(x,{'\', '_GLM_predictors_', '.csv'}, 'CollapseDelimiters',true), names, 'UniformOutput', false);
UID = cellfun(@(x) [x{5} '_' x{6}], C, 'UniformOutput', false);
% generate unique IDs for radii files
[C2,~] = cellfun(@(x) strsplit(x,{'_radii_', '.mat'}, 'CollapseDelimiters',true), names2, 'UniformOutput', false);
UID2 = cellfun(@(x) [x{1} '_' x{2}], C2, 'UniformOutput', false);
sum(cellfun(@(x,y) strcmp(x,y), UID, UID2)) % output should equal 326


%get radii steps from one file
load(names2{1});
radii_steps = S.radii_steps;

% set up empty radii matrix
all_radii = zeros(length(names2), length(radii_steps)); % num cells X radii steps

% colors
rgb_ne_global   = [40, 181, 156]/255;
rgb_ne_local    = [1, 88, 97]/255; 
rgb_ca_global   = [253, 59, 34]/255;
rgb_ca_local    = [148, 36, 36]/255;

%% fit model across radii (low synchrony)

for z = 1:length(radii_steps)-1 % iterate over radii_steps vector 
    
    coeffNames = [];
    coeffVals = [];
    pvals = [];
    
    for i = 1:length(names2) % iterate over all cells (length(names2))
        name = names{i}; % load GLM predictors 
        T = readtable(names{i});
        name2 = names2{i}; % load radii for that cell
        load(name2);
        lg_all = S.lg_all;
        bg_all = S.bg_all;
        green_donut = lg_all(z, :)'; % get local NE for specific radius
        bulk_green = bg_all(z, :)'; % get global NE for specific radius
        mcor = movcorr(green_donut, bulk_green, 450); % generate new NE synchrony time series
        idx = mcor < mean(mcor); % we want low synchrony condition
        X = [green_donut, T.bulk_red];
        %X = [green_donut, bulk_green, green_donut.*bulk_green, T.bulk_red];
        [b, dev, stats] = glmfit(X(idx, :), T.red_tc(idx), 'normal');
        Varnames = {'intercept', 'green_donut', 'bulk_red'};
        coeffNames = [coeffNames; Varnames'];
        pvals = [pvals; stats.p];
        coeffVals = [coeffVals; b];
    end
    
    gdo = coeffVals(strcmp(coeffNames, 'green_donut'));
    pval_gdo = pvals(strcmp(coeffNames, 'green_donut'));
    gdo(find(pval_gdo < 0.01)) = NaN;
    all_radii(:, z) = gdo;

end

%% plotting

figure('Position', [0, 0, 2000, 800]);

% select the radii to plot
new_radii = [5:5:95];
[~,idx] = ismember(new_radii,radii_steps);
all_radii_new = all_radii(:, idx);

% get means
y = squeeze(nanmean(all_radii_new)); 
y = y'; 

% bootstrapping
bootk = 1000;
err = bootci(bootk, @nanmean, all_radii_new(:, isfinite(y))); 
err = err'; 
vals = bootstrp(bootk, @nanmean, all_radii_new(:, isfinite(y))); 

% scatter of points (bootstrapped means)
sz = size(vals);
all_points_y = reshape(vals, [sz(1) * sz(2), 1]); 
all_points_x = sort(repmat([new_radii]', bootk, 1));
scatter((all_points_x + rand(size(all_points_x)) - .5), all_points_y, 3, [.75, .75, .75], 'filled');
hold on;

% % scatter of points (true measured values)
% all_points_y = []; all_points_x = []; 
% sz = size(all_radii_new);
% all_points_y = reshape(all_radii_new, [sz(1) * sz(2), 1]);
% all_points_nans = isnan(all_points_y);
% all_points_y = all_points_y(~all_points_nans);
% all_points_x = sort(repmat([new_radii]', sz(1), 1));
% all_points_x = all_points_x(~all_points_nans);
% scatter((all_points_x + rand(size(all_points_x)) - .5), all_points_y, 3, [.75, .75, .75]);
% hold on;

% error bars and line
errorbar(radii_steps(idx), y, y - err(:, 1), err(:, 2) - y, '.k', 'Capsize', 10);
hold on;
plot(radii_steps(idx), y, 'Color', 'k', 'LineWidth', 1);

% plot significance levels for each radii
p_val = zeros(1, length(new_radii));
for i = 1:length(new_radii)
    test_data = all_radii_new(:,i);
    [~,p,~,~] = ttest(test_data);
    p_val(i) = p;
    B = rmoutliers(test_data);
    if (p < 0.05 && p >= .01)
        text(radii_steps(i), .019, '*', 'FontSize', 40', 'HorizontalAlignment', 'center');
    elseif (p < 0.01 && p >= .001)
        text(radii_steps(i), .019, '**', 'FontSize', 40', 'HorizontalAlignment', 'center');
    elseif (p < 0.001)
        text(radii_steps(i), .019, '***', 'FontSize', 40', 'HorizontalAlignment', 'center');
    end
end

% figure properties
yline(0);
title('Local NE Beta Weight by Radius, All Time Points', 'FontSize', 14);
ax = gca;
ax.FontSize = 35;
ylabel('Beta Weight');
xlabel('Radius (um)');
ylim([-0.02, 0.025]);
yticks([-.02:.01:.03]);
box off;
