%% DESIPRAMINE VS. SALINE RAW DATA EXAMPLES
clear; 

% get saline example for local and global NE
D = 'H:\*PrL*\*GLM predictors*\*reg2*sal*';
files = dir(fullfile(D, '*.csv'));
names = natsortfiles({files.name});
T = readtable(names{50}); %315 %60

% add paths
addpath(genpath('H:\nejrgeco_code_for_publication')); 
addpath(genpath('H:\PrL2.1')); addpath(genpath('H:\PrL3.2')); addpath(genpath('H:\PrL3.4')); 

% set colors
rgb_ne_global   = [40, 181, 156]/255;
rgb_ne_local    = [1, 88, 97]/255; 
rgb_ca_global   = [253, 59, 34]/255;
rgb_ca_local    = [148, 36, 36]/255;

% get desipramine examples
% raw
D2 = 'H:\*PrL*\*data*\';
addpath(genpath(D2));
files2 = dir(fullfile(D2, '*roi_struct*des*.mat'));
names2 = natsortfiles({files2.name});
load(names2{1});
green = reshape(cell2mat({roi_data.tc_g_raw}), 4200, [])';
% detrended
D3 = 'H:\*PrL*\*GLM predictors*\*reg2*des*';
files3 = dir(fullfile(D3, '*.csv'));
names3 = natsortfiles({files3.name});
T2 = readtable(names3{14});

% plot
f = figure;
f.Position = [100, 100, 1800, 900];
% saline
h1 = plot(T.green_donut + 4, 'Color', rgb_ne_local, 'LineWidth', 2); hold on; 
plot(T.bulk_green + 8, 'Color', rgb_ne_local, 'LineWidth', 2); hold on;
% desipramine
plot(T2.green_donut + 12, 'Color', rgb_ca_global, 'LineWidth', 2); hold on;
plot(T2.bulk_green + 16, 'Color', rgb_ca_global, 'LineWidth', 2); ylim([0, 20]); hold on;
plot(zscore(smooth(mean(green), 15/4200)) + 20, 'Color', rgb_ca_global, 'LineWidth', 2); hold on;

% figure properties
yticks([4, 8, 12, 16, 20]); yticklabels({'Local NE', 'Global NE', 'Local NE, detrended', 'Global NE, detrended', 'Global NE'});
ylim([0, 25]);
h = gca; 
h.XAxis.Visible = 'off';
box off
xlim([50, 4150]);
xticks(50:900:4150);
xticklabels({'', '10', '20', '30', '40', '50', '60', '70', '80', '90'});
xlabel('Time (seconds)');
ax = gca;
sgtitle('Raw Data Inputs to GLM', 'FontSize', 30);

%% COMPARE SALINE AND DESIPRAMINE GLM (note this is for all time points)
% get data (saline)
clear;
D = 'H:\*PrL*\*GLM predictors*\*reg2*sal*';
files = dir(fullfile(D, '*.csv'));
names = cell(size(files,1), 1);
for i = 1:length(names)
    names{i} = [files(i).folder '\' files(i).name];
end

% fit models (saline)
coeffNames = [];
coeffVals = [];
pvals = [];
for i = 1:length(names)
    name = names{i};
    T = readtable(names{i});
    X = [T.green_donut, T.bulk_green, T.bulk_red, T.green_donut.*T.bulk_green];
    [b, dev, stats] = glmfit(X, T.red_tc, 'normal');
    Varnames = {'intercept', 'green_donut', 'bulk_green', 'bulk_red', 'green_donut:bulk_green'};
    coeffNames = [coeffNames; Varnames'];
    pvals = [pvals; stats.p];
    coeffVals = [coeffVals; b];
end

% get predictors of interest (saline) and set up Xcell and Xmean
gdo = coeffVals(strcmp(coeffNames, 'green_donut'));
pval_gdo = pvals(strcmp(coeffNames, 'green_donut'));
gdo = {gdo};
bg = coeffVals(strcmp(coeffNames, 'bulk_green'));
pval_bg = pvals(strcmp(coeffNames, 'bulk_green'));
bg = {bg};
Xcell = [gdo bg];
Xmean = [mean(gdo{:}) mean(bg{:})];


% get data (desipramine)
D = 'H:\*PrL*\*GLM predictors*\*reg2*des*';
files = dir(fullfile(D, '*.csv'));
names = cell(size(files,1), 1);
for i = 1:length(names)
    names{i} = [files(i).folder '\' files(i).name];
end

% fit models (desipramine)
coeffNames = [];
coeffVals = [];
pvals = [];
for i = 1:length(names)
    name = names{i};
    T = readtable(names{i});
    X = [T.green_donut, T.bulk_green, T.bulk_red, T.green_donut.*T.bulk_green];
    [b, dev, stats] = glmfit(X, T.red_tc, 'normal');
    Varnames = {'intercept', 'green_donut', 'bulk_green', 'bulk_red', 'green_donut:bulk_green'};
    coeffNames = [coeffNames; Varnames'];
    pvals = [pvals; stats.p];
    coeffVals = [coeffVals; b];
end

% get predictors of interest (desipramine) and add to Xcell and Xmean
gdo = coeffVals(strcmp(coeffNames, 'green_donut'));
pval_gdo = pvals(strcmp(coeffNames, 'green_donut'));
gdo = {gdo};
bg = coeffVals(strcmp(coeffNames, 'bulk_green'));
pval_bg = pvals(strcmp(coeffNames, 'bulk_green'));
bg = {bg};
Xcell = [Xcell gdo bg];
Xmean = [Xmean mean(gdo{:}) mean(bg{:})];

% set up data for plotting
Xcell = Xcell([1 3 2 4]); % reorder Xcell so that local green and bulk green from both conditions are side-by-side
E = []; % get error
for i = 1:length(Xcell)
    [~,~,ci,~] = ttest(Xcell{i});
    E = [E ci];
end
error =  E(2,:)-E(1,:);
P = NaN(4, 4); % get p matrix
pairs =  [1  2; 3 4];
for i = 1:size(pairs,1)
[~,p, ~,~] = ttest2(Xcell{pairs(i, 1)}, Xcell{pairs(i, 2)});
P(pairs(i, 1), pairs(i, 2)) = p; P(pairs(i, 2), pairs(i, 1)) = p; 
end


% plot
superbar(Xmean([1 3 2 4]), 'E', error([1 3 2 4]), 'P', P, 'BarFaceColor', [0.2 0.9 0.4; 0.1 0.6 0.1; 0.2 0.9 0.4; 0.1 0.6 0.1]);

% figure properties
xticks([1.5, 3.5]);
xticklabels({'Local NE', 'global NE'});
yticks([-0.02, 0, 0.06]);
xtickangle(45);
legend({'saline', 'desipramine'});
title('Saline vs Desipramine Predictor Distributions');

%% get full stats for t-tests

pairs =  [1  2; 3 4];
for i = 1:size(pairs,1)
[h,p,ci,stats] = ttest2(Xcell{pairs(i, 1)}, Xcell{pairs(i, 2)})
end
