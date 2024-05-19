%% PLOT RAW INPUTS TO GLM
clear;

% load data
addpath(genpath('F:\PrL2.1')); addpath(genpath('F:\PrL3.2')); addpath(genpath('F:\PrL3.4')); addpath(genpath('F:\nejrgeco_code_forpublication\'));
D = 'F:\*PrL*\*GLM predictors*\*reg2*sal*';
files = dir(fullfile(D, '*.csv'));
names = cell(size(files,1), 1);
for i = 1:length(names)
    names{i} = [files(i).folder '\' files(i).name];
end
T = readtable(names{315});
%setup_figprop;

% colors
rgb_ne_global   = [40, 181, 156]/255;
rgb_ne_local    = [1, 88, 97]/255; 
rgb_ca_global   = [253, 59, 34]/255;
rgb_ca_local    = [148, 36, 36]/255;

% plot local NE, global NE, cell Ca2+, and global Ca2+
f = figure;
f.Position = [100, 100, 1800, 900];
subplot(6,1,[1 2 3 4]);
plot(T.green_donut + 4, 'Color', rgb_ne_local, 'LineWidth', 2); hold on; 
plot(T.bulk_green + 8, 'Color', rgb_ne_global, 'LineWidth', 2); hold on;
plot(T.bulk_red + 12, 'Color', rgb_ca_global, 'LineWidth', 2); hold on;
plot(T.red_tc + 16, 'Color', rgb_ca_local, 'LineWidth', 2); ylim([0, 20]);

% figure properties
yticks([4, 8, 12, 16]); yticklabels({'Local NE', 'Global NE', 'Global Ca2+', 'Cell Ca2+'});
h = gca; 
h.XAxis.Visible = 'off';
box off
xlim([50, 4150]);

% plot NE synchrony
subplot(6, 1, [5 6]);
plot(T.mcor_bglg, 'k', 'LineWidth', 2);

% figure properties
yticks(0.75); yticklabels('NE Synchrony');
box off
xlim([50, 4150]);
xticks(50:900:4150);
xticklabels({'', '10', '20', '30', '40', '50', '60', '70', '80', '90'});
xlabel('Time (seconds)');
sgtitle('Raw Data Inputs to GLM', 'FontSize', 30);

%% forced entry glm (all times)
% load data
clear;
D = 'F:\*PrL*\*GLM predictors*\*reg2*sal*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('F:\PrL2.1')); addpath(genpath('F:\PrL3.2')); addpath(genpath('F:\PrL3.4')); addpath(genpath('F:\samira_nejrgeco'));
names = cell(size(files,1), 1);
for i = 1:length(names)
    names{i} = [files(i).folder '\' files(i).name];
end


% fit model for all cells
pvals = [];
coeffNames = [];
coeffVals = [];
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

% plot
plot_GLM(coeffNames, coeffVals, 'all time points')
%% forced entry glm by NE synchrony (LOW NE Synchrony)
% load data
clear;
D = 'H:\*PrL*\*GLM predictors*\*reg2*sal*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('H:\PrL2.1')); addpath(genpath('H:\PrL3.2')); addpath(genpath('H:\PrL3.4')); addpath(genpath('H:\nejrgeco_code_for_publication'));
names = cell(size(files,1), 1);
for i = 1:length(names)
    names{i} = [files(i).folder '\' files(i).name];
end

% fit model for all cells
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

% plot
plot_GLM(coeffNames, coeffVals, 'low synchrony');

%%
% get means and run 1-sample t-tests
gdo = coeffVals(strcmp('green_donut', coeffNames));
disp("mean beta local NE")
mean(gdo)
disp("SD local NE")
std(gdo)
disp("SE local NE")
std(gdo)./sqrt(length(gdo))
disp("1-sample t-test local NE")
[h,p,ci,stats] = ttest(gdo)

bulk_red = coeffVals(strcmp('bulk_red', coeffNames));
disp("mean beta global Ca2+")
mean(bulk_red)
disp("SD local NE")
std(bulk_red)
disp("SE local NE")
std(bulk_red)./sqrt(length(bulk_red))
disp("1-sample t-test local NE")
[h,p,ci,stats] = ttest(bulk_red)


%% forced entry glm by NE synchrony (HIGH NE Synchrony)
% load data
clear;
D = 'F:\*PrL*\*GLM predictors*\*reg2*sal*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('F:\PrL2.1')); addpath(genpath('F:\PrL3.2')); addpath(genpath('F:\PrL3.4')); addpath(genpath('F:\nejrgeco_code_for_publication'));
names = cell(size(files,1), 1);
for i = 1:length(names)
    names{i} = [files(i).folder '\' files(i).name];
end

% fit model for all cells
pvals = [];
coeffNames = [];
coeffVals = [];
for i = 1:length(names)
    name = names{i};
    T = readtable(names{i});
    condition = T.mcor_all > mean(T.mcor_all);
    X = [T.green_donut(condition), T.bulk_green(condition), T.bulk_red(condition), T.green_donut(condition).*T.bulk_green(condition)];
    [b, dev, stats] = glmfit(X, T.red_tc(condition), 'normal');
    Varnames = {'intercept', 'green_donut', 'bulk_green', 'bulk_red', 'green_donut:bulk_green'};
    coeffNames = [coeffNames; Varnames'];
    pvals = [pvals; stats.p];
    coeffVals = [coeffVals; b];
end

% plot
plot_GLM(coeffNames, coeffVals, 'high synchrony');

%% plotting function

% inputs: coeffNames, coeffVals, condition
function plot_GLM(coeffNames, coeffVals, condition)
    f = figure;
    violin_mat = zeros(length(coeffNames)/5, 4);
    predictor_names = [{'green_donut'} {'bulk_green'} {'bulk_red'} {'green_donut:bulk_green'}];
    labels = {'local NE', 'global NE', 'global Ca', 'local NE:global NE'};
    for i = 1:4
        violin_mat(:,i) = coeffVals(strcmp(coeffNames, predictor_names{i}));
    end
    b = violin(violin_mat, 'medc', []); 
    
    p_val = zeros(1, 4);
    for i = 1:4
        [~,p,~,~] = ttest(coeffVals(strcmp(coeffNames, predictor_names{i})));
        p_val(i) = p;
        p
        if (p < 0.05 && p >= 0.01)
            text(i, .8, '', 'FontSize', 40, 'HorizontalAlignment', 'center');
        elseif (p < 0.01 && p >= 0.001)
            text(i, .8, '', 'FontSize', 40, 'HorizontalAlignment', 'center');
        elseif (p < 0.001)
            text(i, .8, '***', 'FontSize', 40, 'HorizontalAlignment', 'center');
        end
    end
    
    yline(0);
    ylim([-.5,1]);
    xticks(1:6);
    xtickangle(30);
    ax = gca;
    set(gca, 'TickLabelInterpreter', 'none');
    xticklabels(labels);
    yticks([-1, -0.5, 0, 0.5, 1]);
    title(['Beta Weight of GLM Predictors, ' condition]);
    ylabel('Beta Weight');
    box off
    axis square
    legend off
    set(b, {'linew'}, {2})
end
