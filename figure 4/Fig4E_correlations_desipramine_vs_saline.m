%% FIGURE 4E local:local green desipramine vs. saline
% set file paths
addpath(genpath('F:\PrL2.1')); addpath(genpath('F:\PrL3.2')); addpath(genpath('F:\PrL3.4'));

% load in desipramine file names
D = 'F:\*PrL*\*GLM predictors*\*reg2_*des*';
files = dir(fullfile(D, '*.csv'));
all_names = cell(size(files,1), 1);
for i = 1:length(all_names)
    all_names{i} = [files(i).folder '\' files(i).name];
end

% get local:local correlations (desipramine, lg_lg_des)
within_group_correlations_des = [];
run_names = cellfun(@(x)x(38:48), all_names, 'UniformOutput', false);
unique_run_names = unique(run_names);
lg_lg_des = [];
for i = 1:length(unique_run_names)
    file_names = all_names(contains(all_names, unique_run_names{i}));
    all_lg = [];
    for k = 1:length(file_names)
        T = readtable(file_names{k});
        all_lg = [all_lg T.green_donut];
    end
    cor_pairs = corr(all_lg);
    lg_lg_des = [lg_lg_des cor_pairs(find(tril(cor_pairs, -1) ~= 0))'];
end

% load in saline file names
D = 'F:\*PrL*\*GLM predictors*\*reg2_*sal*';
files = dir(fullfile(D, '*.csv'));
all_names = cell(size(files,1), 1);
for i = 1:length(all_names)
    all_names{i} = [files(i).folder '\' files(i).name];
end

% get local:local correlations (saline, lg_lg_sal)
within_group_correlations_sal = [];
run_names = cellfun(@(x)x(38:48), all_names, 'UniformOutput', false);
unique_run_names = unique(run_names);
lg_lg_sal = []
for i = 1:length(unique_run_names)
    file_names = all_names(contains(all_names, unique_run_names{i}));
    all_lg = [];
    for k = 1:length(file_names)
        T = readtable(file_names{k});
        all_lg = [all_lg T.green_donut];
    end
    cor_pairs = corr(all_lg);
    lg_lg_sal = [lg_lg_sal cor_pairs(find(tril(cor_pairs, -1) ~= 0))'];
end

% get cumulative probability densities and run KS-test
[f_des,x_des] = ecdf(lg_lg_des);
[f_sal,x_sal] = ecdf(lg_lg_sal);
[~,p] = kstest2(lg_lg_des, lg_lg_sal)

% plot
f = figure;
plot(x_des, f_des, 'LineWidth', 2, 'Color', 'm'); 
hold on; 
plot(x_sal, f_sal, 'LineWidth', 2, 'Color', 'b');
%if (p < 0.001)
    %text(0.2, 0.7, 'KS test: ***', 'FontSize', 20);
%end

% figure properties
f.Position = [100, 100, 1500, 750];
title('CDF of Local/Local NE Correlations');
legend({'desipramine', 'saline'}, 'FontSize', 20, 'Location', 'northwest');
ylabel('Cumulative fraction');
xlabel('Correlation');
yticks([0, 0.25, 0.50, 0.75, 1.00]);
xticks([0, 0.25, 0.50, 0.75, 1.00]);
xlim([0, 1]);
box off;
legend boxoff;

%% FIGURE 4E local:global green desipramine vs. saline

% set file paths
addpath(genpath('F:\PrL2.1')); addpath(genpath('F:\PrL3.2')); addpath(genpath('F:\PrL3.4'));

% load in desipramine file names
D = 'F:\*PrL*\*GLM predictors*\*reg2_*des*';
files = dir(fullfile(D, '*.csv'));
all_names = cell(size(files,1), 1);
for i = 1:length(all_names)
    all_names{i} = [files(i).folder '\' files(i).name];
end

% get local:global green desipramine correlations (lg_bg_des)
lg_bg_des = [];
for i = 1:length(all_names)
    T = readtable(all_names{i});
    cor = corr(T.green_donut, T.bulk_green);
    lg_bg_des = [lg_bg_des cor];
end

% load in saline file names
D = 'F:\*PrL*\*GLM predictors*\*reg2_*sal*';
files = dir(fullfile(D, '*.csv'));
all_names = cell(size(files,1), 1);
for i = 1:length(all_names)
    all_names{i} = [files(i).folder '\' files(i).name];
end

% get local:global green saline correlations (lg_bg_sal)
lg_bg_sal = [];
for i = 1:length(all_names)
    T = readtable(all_names{i});
    cor = corr(T.green_donut, T.bulk_green);
    lg_bg_sal = [lg_bg_sal cor];
end

% get cumulative probability distribution for correlations
% also run KS-test
[f_des,x_des] = ecdf(lg_bg_des);
[f_sal,x_sal] = ecdf(lg_bg_sal);
[~,p] = kstest2(lg_bg_des, lg_bg_sal)

% plot
f = figure;
plot(x_des, f_des, 'LineWidth', 2, 'Color', 'm'); 
hold on; 
plot(x_sal, f_sal, 'LineWidth', 2, 'Color', 'b');
%if (p < 0.001)
    %text(0.2, 0.7, 'KS test: ***', 'FontSize', 20);
%end

% figure properties
f.Position = [100, 100, 1500, 750]
title('Local/Global NE Correlations');
legend({'desipramine', 'saline'}, 'FontSize', 20, 'Location', 'northwest');
ylabel('Cumulative Fraction');
xlabel('Correlation');
yticks([0, 0.25, 0.50, 0.75, 1.00]);
xticks([0, 0.25, 0.50, 0.75, 1.00]);
xlim([0, 1]);
box off;
legend boxoff;
