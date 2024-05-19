%% COMPARE REGRESSED AND NON-REGRESSED

% get data (non-regressed)
D = 'F:\*PrL*\*GLM predictors*\*noreg*sal*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('F:\PrL2.1')); addpath(genpath('F:\PrL3.2')); addpath(genpath('F:\PrL3.4')); addpath(genpath('F:\nejrgeco_code_for_publication'));
names = cell(size(files,1), 1);
for i = 1:length(names)
    names{i} = [files(i).folder '\' files(i).name];
end
setup_figprop;

% fit models (non-regressed)
pvals = [];
coeffNames = [];
coeffVals = [];
for i = 1:length(names)
    name = names{i};
    T = readtable(names{i});
    condition = T.mcor_all < mean(T.mcor_all);
    X = [T.green_donut(condition), T.bulk_green(condition), T.bulk_red(condition), T.green_donut(condition).*T.bulk_green(condition)];
    [b, dev, stats] = glmfit(X, T.red_tc(condition), 'normal');
    Varnames = {'intercept', 'green_donut', 'bulk_green', 'bulk_red', 'green_donut:bulk_green'};
    coeffNames = [coeffNames; Varnames'];
    pvals = [pvals; stats.p];
    coeffVals = [coeffVals; b];
end

% separate results by predictor (non-regressed)
violin_mat_reg = zeros(length(coeffNames)/5, 4);
predictor_names = [{'green_donut'} {'bulk_green'} {'bulk_red'} {'green_donut:bulk_green'}];
labels = {'local NE', 'global NE', 'global Ca', 'local NE:global NE'};
for i = 1:4
    violin_mat_reg(:,i) = coeffVals(strcmp(coeffNames, predictor_names{i}));
end

% get data (regressed)
D = 'F:\*PrL*\*GLM predictors*\*reg2*sal*';
files = dir(fullfile(D, '*.csv'));
names = cell(size(files,1), 1);
for i = 1:length(names)
    names{i} = [files(i).folder '\' files(i).name];
end

% fit models (regressed)
pvals = [];
coeffNames = [];
coeffVals = [];
for i = 1:length(names)
    name = names{i};
    T = readtable(names{i});
    condition = T.mcor_all < mean(T.mcor_all);
    X = [T.green_donut(condition), T.bulk_green(condition), T.bulk_red(condition), T.green_donut(condition).*T.bulk_green(condition)];
    [b, dev, stats] = glmfit(X, T.red_tc(condition), 'normal');
    Varnames = {'intercept', 'green_donut', 'bulk_green', 'bulk_red', 'green_donut:bulk_green'};
    coeffNames = [coeffNames; Varnames'];
    pvals = [pvals; stats.p];
    coeffVals = [coeffVals; b];
end

% separate results by predictor (regressed)
violin_mat = zeros(length(coeffNames)/5, 4);
predictor_names = [{'green_donut'} {'bulk_green'} {'bulk_red'} {'green_donut:bulk_green'}];
labels = {'local NE', 'global NE', 'global Ca', 'local NE:global NE'};
for i = 1:4
    violin_mat(:,i) = coeffVals(strcmp(coeffNames, predictor_names{i}));
end

% prep data for plotting
[h,p] = ttest2(violin_mat_reg(:,1), violin_mat(:,1))
X = [violin_mat_reg(:,1) violin_mat(:,1) violin_mat_reg(:,4) violin_mat(:,4)]; % input matrix for superbar function
E = std(X)/sqrt(size(X,1)); % standard error
xlabels = {'local NE', 'global NE'};
P = NaN(4, 4); % make P matrix for superbar function
pairs =  [1  2; 3 4];
for i = 1:size(pairs,1)
[~,p, ~,~] = ttest2(X(:, pairs(i, 1)), X(:, pairs(i,2)));
P(pairs(i, 1), pairs(i, 2)) = p; P(pairs(i, 2), pairs(i, 1)) = p; 
end

%plot 
f = figure;
f.Position = [100, 100, 500 800]
superbar(mean(X), 'E', E, 'P', P, 'BarFaceColor', [0.2 0.9 0.4; 0.1 0.6 0.1; 0.2 0.9 0.4; 0.1 0.6 0.1]);

% figure properties
xticklabels(xlabels);
xticks([1.5, 3.5]);
xticklabels({'Local NE', 'global NE'});
yticks([-0.02, 0, 0.06]);
xtickangle(0);
legend({'non-regressed', 'regressed'});
title('Effect of Regression on Results', 'FontSize', 20);
ax = gca;
xlim([0, 2.5]);