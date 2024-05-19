%% LOCAL/LOCAL vs LOCAL/GLOBAL

D = 'F:\PrL*\*GLM predictors*\*reg_*sal*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('F:\nejrgeco_code_for_publication'));
addpath(genpath('F:\PrL2.1')); addpath(genpath('F:\PrL3.2')); addpath(genpath('F:\PrL3.4')); 
all_names = natsortfiles({files.name});

within_group_correlations_sal = [];
run_names = cellfun(@(x)x(1:11), all_names, 'UniformOutput', false);
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

% LOCAL/GLOBAL CORRELATIONS
lg_bg_sal = [];
for i = 1:length(all_names)
    T = readtable(all_names{i});
    cor = corr(T.green_donut, T.bulk_green);
    lg_bg_sal = [lg_bg_sal cor];
end

%%
clear;

% paths
addpath('C:\samira_nejrgeco');
setup_figprop;

% colors
rgb_ne_global   = [40, 181, 156]/255;
rgb_ne_local    = [1, 88, 97]/255; 
rgb_ca_global   = [253, 59, 34]/255;
rgb_ca_local    = [148, 36, 36]/255;

D = 'C:\*PrL*\*GLM predictors*\*reg_*sal*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('C:\PrL2.1')); addpath(genpath('C:\PrL3.2')); addpath(genpath('C:\PrL3.4')); addpath(genpath('C:\samira_nejrgeco'));
all_names = natsortfiles({files.name});

%% RAW SYNCHRONY TRACES FOR FIG 2E

D = 'F:\*PrL*\*GLM predictors*\*reg3_*sal*';
files = dir(fullfile(D, '*PrL2.1_sal2*.csv'));
addpath(genpath('F:\PrL2.1')); addpath(genpath('F:\PrL3.2')); addpath(genpath('F:\PrL3.4')); addpath(genpath('F:\nejrgeco_code_for_publication'));
all_names = natsortfiles({files.name});

f = figure;
f.Position = [100, 100, 1500, 500];
avg_synchrony = zeros(4200,1);
for i = 1:length(all_names)
    if i == 6
%         T = readtable(all_names{i});
%         plot(T.mcor_bglg, 'LineWidth', 3, 'color', (rgb_ne_global + rgb_ne_local)/2); 
%         hold on;
    else
        T = readtable(all_names{i});
        plot(T.mcor_bglg, 'LineWidth', 1, 'color', [0.75 0.75 0.75 1]); 
        hold on;
    end
    avg_synchrony = avg_synchrony+T.mcor_bglg;
end
T = readtable(all_names{6});
plot(T.mcor_bglg, 'LineWidth', 3, 'color', (rgb_ne_global + rgb_ne_local)/2); 
avg_synchrony = avg_synchrony/length(all_names); 
plot(avg_synchrony, 'LineWidth', 5, 'color', [0 0 0]); 

ylim([-0.4, 1]);
xlim([0, 4200]);
yticks([-0.25:0.25:1]);
ylabel('Correlation (r)'); 
xticks(0:900:4200); 
xticklabels({'', '1', '2', '3', '4', '5'}); 
xlabel('Time (minutes)')
box off; 

%% RAW DATA EXAMPLE (global ne - local ne - correlation)
D = 'C:\*PrL*\*GLM predictors*\*reg3_*sal*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('C:\PrL2.1')); addpath(genpath('C:\PrL3.2')); addpath(genpath('C:\PrL3.4')); addpath(genpath('C:\samira_nejrgeco'));
all_names = natsortfiles({files.name});

f = figure;
f.Position = [100, 100, 1500, 500];
i = 6; %5
T = readtable(all_names{i});

% global ne
subplot(3, 1, 1);
plot(T.bulk_green, 'LineWidth', 1.5, 'Color', rgb_ne_global); 
ylim([-2.5, 2.5]); yticks([0]); yticklabels('global NE');
box off; ax = gca; ax.XAxis.Visible = 'off'; xlim([1000, 2500]); 
%ax = gca; ax.LineWidth = 2; ax.FontSize = 15;

% local ne
subplot(3, 1, 2);
plot(T.green_donut, 'LineWidth', 1.5, 'Color', rgb_ne_local); 
ylim([-2.5, 2.5]); yticks([0]); yticklabels('local NE');
box off; ax = gca; ax.XAxis.Visible = 'off'; xlim([1000, 2500]); 
%ax = gca; ax.LineWidth = 2; ax.FontSize = 15;

% correlation 
subplot(3,1,3);
plot(T.mcor_bglg, 'LineWidth', 2, 'Color', [.5, .5, .5]); 
xticks(0:450:4200); 
xticklabels({'', '0.5' '1', '1.5', '2', '2.5', '3', '3.5', '4', '4.5', '5'}); 
ylabel('correlation'); 
xlim([1000, 2500]);
xlabel('Time (minutes)')
%ax = gca; ax.LineWidth = 2; ax.FontSize = 15;
box off;

sgtitle('NE Synchrony Example', 'FontSize', 20);

%% DESIPRAMINE VS. SALINE RAW TIME COURSES EXAMPLES
D = 'C:\*PrL*\*GLM predictors*\*reg2_*sal*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('C:\PrL2.1')); addpath(genpath('C:\PrL3.2')); addpath(genpath('C:\PrL3.4')); addpath(genpath('C:\samira_nejrgeco'));
sal_names = natsortfiles({files.name});

D = 'C:\*PrL*\*GLM predictors*\*reg2_*des*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('C:\PrL2.1')); addpath(genpath('C:\PrL3.2')); addpath(genpath('C:\PrL3.4')); addpath(genpath('C:\samira_nejrgeco'));
des_names = natsortfiles({files.name});

f = figure;
f.Position = [0, 0, 2000, 1500];

subplot(2,1,1);
for i = [10, 40, 80, 120, 160, 200] 
sal_names{i}
T = readtable(sal_names{i});
plot(T.bulk_green+i/10, 'g', 'LineWidth', 2, 'Color', rgb_ne_global);  xticks([]); xlim([10, 4100]); ylabel('saline');  yticks([]);
%ax = gca; ax.LineWidth = 2; ax.FontSize = 15; xlabel('');
hold on;
end
box off;

subplot(2,1,2);
for i = [10, 40, 80, 120, 160, 200] 
T = readtable(des_names{i});
plot(T.bulk_green+i/10, 'b', 'LineWidth', 2, 'Color', rgb_ne_global);  xticks(0:900:4200); xlim([10, 4100]); ylabel('desipramine'); yticks([]);
hold on;
end
xticklabels({'', '1', '2', '3', '4', '5'}); 
ax = gca; ax.LineWidth = 2; ax.FontSize = 15; xlabel('Time (minutes)');
box off;

%ax = gca; ax.LineWidth = 2; ax.FontSize = 15;
sgtitle('Desipramine vs. Saline Example', 'FontSize', 20);

%% DESIPRAMINE VS. SALINE RAW TIME COURSES EXAMPLES 
D = 'C:\*PrL*\*GLM predictors*\*reg2_*sal*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('C:\PrL2.1')); addpath(genpath('C:\PrL3.2')); addpath(genpath('C:\PrL3.4')); addpath(genpath('C:\samira_nejrgeco'));
sal_names = natsortfiles({files.name});

D = 'C:\*PrL*\*GLM predictors*\*reg2_*des*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('C:\PrL2.1')); addpath(genpath('C:\PrL3.2')); addpath(genpath('C:\PrL3.4')); addpath(genpath('C:\samira_nejrgeco'));
des_names = natsortfiles({files.name});

D = 'C:\*PrL*\*data*\';
files = dir(fullfile(D, '*roi_struct*des*.mat'));
addpath(genpath('C:\PrL2.1')); addpath(genpath('C:\PrL3.2')); addpath(genpath('C:\PrL3.4')); addpath(genpath('C:\samira_nejrgeco'));
all_names = natsortfiles({files.name});

f = figure;
f.Position = [0, 0, 2000, 1500];

subplot(3,1,1);
for i = [10, 40, 80, 120, 160, 200] 
sal_names{i}
T = readtable(sal_names{i});
plot(T.bulk_green+i/10, 'g', 'LineWidth', 2, 'Color', rgb_ne_global);  xticks([]); xlim([10, 4100]); ylabel('saline');  yticks([]);
%ax = gca; ax.LineWidth = 2; ax.FontSize = 15; xlabel('');
hold on;
end
box off;

subplot(3,1,2);
y = [1, 4, 8, 12, 16, 20];
x = [1, 7, 9, 11, 12];
for i = 1:5
    load(all_names{x(i)});
    green = reshape(cell2mat({roi_data.tc_g_raw}), 4200, [])';
    plot(zscore(smooth(mean(green), 15/4200)) + y(i));
    hold on;
end
xticks([]); xlim([10, 4100]); ylabel({'global NE desipramine', '(non-detrended)'});  yticks([]); box off;

subplot(3,1,3);
for i = [10, 40, 80, 120, 160, 200] 
T = readtable(des_names{i});
plot(T.bulk_green+i/10, 'b', 'LineWidth', 2, 'Color', rgb_ne_global);  xticks(0:900:4200); xlim([10, 4100]); ylabel({'desipramine', '(detrended)'}); yticks([]);
hold on;
end
xticklabels({'', '1', '2', '3', '4', '5'}); 
ax = gca; ax.LineWidth = 2; ax.FontSize = 15; xlabel('Time (minutes)');
box off;

%ax = gca; ax.LineWidth = 2; ax.FontSize = 15;
sgtitle('Desipramine vs. Saline Example', 'FontSize', 20);

%% DESIPRAMINE VS. SALINE RAW DATA EXAMPLES UPDATED
D = 'F:\*PrL*\*GLM predictors*\*reg2*sal*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('F:\PrL2.1')); addpath(genpath('F:\PrL3.2')); addpath(genpath('F:\PrL3.4')); addpath(genpath('F:\samira_nejrgeco'));
names = natsortfiles({files.name});
T = readtable(names{50}); %50

D = 'F:\*PrL*\*data*\';
files = dir(fullfile(D, '*roi_struct*des*.mat'));
addpath(genpath('F:\PrL2.1')); addpath(genpath('F:\PrL3.2')); addpath(genpath('F:\PrL3.4')); addpath(genpath('F:\samira_nejrgeco'));
all_names = natsortfiles({files.name});
load(all_names{1});
green = reshape(cell2mat({roi_data.tc_g_raw}), 4200, [])';

rgb_ne_global   = [40, 181, 156]/255;
rgb_ne_local    = [1, 88, 97]/255; 
rgb_ca_global   = [253, 59, 34]/255;
rgb_ca_local    = [148, 36, 36]/255;

f = figure;
f.Position = [100, 100, 1800, 900];

h1 = plot(T.green_donut + 4, 'Color', rgb_ne_local, 'LineWidth', 2); hold on; 
plot(T.bulk_green + 8, 'Color', rgb_ne_local, 'LineWidth', 2); hold on;

D = 'F:\*PrL*\*GLM predictors*\*reg2*des*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('F:\PrL2.1')); addpath(genpath('F:\PrL3.2')); addpath(genpath('F:\PrL3.4')); addpath(genpath('F:\samira_nejrgeco'));
names = natsortfiles({files.name});
T = readtable(names{90}); %14 %40 #90 %180 %260

h2 = plot(T.green_donut + 12, 'Color', rgb_ca_global, 'LineWidth', 2); hold on;
plot(T.bulk_green + 16, 'Color', rgb_ca_global, 'LineWidth', 2); ylim([0, 20]); hold on;
plot(zscore(smooth(mean(green), 15/4200)) + 20, 'Color', rgb_ca_global, 'LineWidth', 2); hold on;

yticks([4, 8, 12, 16, 20]); yticklabels({'Local NE', 'Global NE', 'Local NE, detrended', 'Global NE, detrended', 'Global NE'});
ylim([0, 25]);
h = gca; 
h.XAxis.Visible = 'off';
% h.FontSize = 25;
% h.LineWidth = 6;
box off
xlim([50, 4150]);
xticks(50:900:4150);
xticklabels({'', '10', '20', '30', '40', '50', '60', '70', '80', '90'});
xlabel('Time (seconds)');
legend([h1, h2], {'saline', 'desipramine'});
ax = gca;
% ax.LineWidth = 6;
% ax.FontSize = 25;
sgtitle('Raw Data Inputs to GLM', 'FontSize', 30);

%%
figure;
y = [1, 4, 8, 12, 16, 20];
x = [1, 7, 9, 11, 12];
for i = 1
    load(all_names{x(i)});
    green = reshape(cell2mat({roi_data.tc_g_raw}), 4200, []);
    for i = 1:length(roi_data)
        green(:,i) = smooth(zscore(green(:,i)'), 15/4200);
    end
    for k = 1:length(roi_data)
        plot(green(:,k)+2*k); hold on;
        xlim([1000, 2000]);
    end
end

xticks([]); xlim([1000, 3000]); ylabel({'global NE desipramine', '(non-detrended)'});  yticks([]); box off;


%% LOCAL/RED CORRELATIONS
lg_lr_sal = [];
for i = 1:length(all_names)
    T = readtable(all_names{i});
    cor = corr(T.green_donut, T.red_tc);
    lg_lr_sal = [lg_lr_sal cor];
end
%%
D = 'C:\PrL*\*GLM predictors*\*reg_*des*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('C:\PrL2.1')); addpath(genpath('C:\PrL3.2')); addpath(genpath('C:\PrL3.4')); addpath(genpath('C:\samira_nejrgeco'));
all_names = natsortfiles({files.name});

%%
within_group_correlations_sal = [];
run_names = cellfun(@(x)x(1:11), all_names, 'UniformOutput', false);
unique_run_names = unique(run_names);
lg_lg_des = []
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

%% LOCAL/GLOBAL CORRELATIONS
lg_bg_des = [];
for i = 1:length(all_names)
    T = readtable(all_names{i});
    cor = corr(T.green_donut, T.bulk_green);
    lg_bg_des = [lg_bg_des cor];
end
%% LOCAL/RED CORRELATIONS
lg_lr_des = [];
for i = 1:length(all_names)
    T = readtable(all_names{i});
    cor = corr(T.green_donut, T.red_tc);
    lg_lr_des = [lg_lr_des cor];
end
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

%% ALL CELLS CORRELATION MATRIX
run_names = cellfun(@(x)x(1:11), all_names, 'UniformOutput', false);
run_cutoffs = [];
unique_run_names = unique(run_names);
for i = 1:length(unique_run_names)
    tf = strcmp(run_names, unique_run_names{i});
    run_cutoffs = [run_cutoffs max(find(tf == 1))];
end
run_cutoffs = run_cutoffs(find(run_cutoffs < max(run_cutoffs)));

mcor1 = ones(4200, length(all_names));
mcor2 = ones(4200, length(all_names));
for i = 1:length(all_names)
    T = readtable(all_names{i});
    mcor1(:, i) = movcorr(T.bulk_green(1:4200), T.green_donut(1:4200), 450);
    mcor2(:, i) = movcorr(T.bulk_green(1:4200), T.red_tc(1:4200), 450);
end

figure;
ax = axes;
imagesc(corr(mcor1)); hold on;
for i = 1:length(run_cutoffs)
    yline(run_cutoffs(i)+0.5, 'LineWidth', 3, 'Color', [0 0 0]);
    xline(run_cutoffs(i)+0.5, 'LineWidth', 3, 'Color', [0 0 0]);
    hold on;
end
yticks(12:28:320);
yticklabels(unique_run_names);
ax.TickLabelInterpreter = 'none';
title("Correlation Matrix of all runs Local/Bulk Green Sliding Correlation"); colorbar;
%% COMPARE DESIPRAMINE VS. SALINE WITHIN RUN AND WITHIN CELL CORRELATIONS
D = 'C:\*PrL*\*GLM predictors*\*reg_*sal*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('C:\PrL2.1')); addpath(genpath('C:\PrL3.2')); addpath(genpath('C:\PrL3.4')); addpath(genpath('C:\samira_nejrgeco'));
all_names = natsortfiles({files.name});

within_group_correlations_sal = [];
run_names = cellfun(@(x)x(1:11), all_names, 'UniformOutput', false);
unique_run_names = unique(run_names);
for i = 1:length(unique_run_names)
    file_names = all_names(contains(all_names, unique_run_names{i}));
    T = readtable(file_names{1});
    mcor = ones(length(T.red_tc), length(file_names));
    for k = 1:length(file_names)
        T = readtable(file_names{k});
        mcor(:, k) = movcorr(T.bulk_green, T.green_donut, 450);
    end    
    cor_pairs = corr(mcor);
    within_group_correlations_sal = [within_group_correlations_sal cor_pairs(find(tril(cor_pairs, -1) ~= 0))'];
end

D = 'C:\*PrL*\*GLM predictors*\*reg_*des*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('C:\PrL2.1')); addpath(genpath('C:\PrL3.2')); addpath(genpath('C:\PrL3.4')); addpath(genpath('C:\samira_nejrgeco'));
all_names = natsortfiles({files.name});

within_group_correlations_des = [];
run_names = cellfun(@(x)x(1:11), all_names, 'UniformOutput', false);
unique_run_names = unique(run_names);
for i = 1:length(unique_run_names)
    file_names = all_names(contains(all_names, unique_run_names{i}));
    T = readtable(file_names{1});
    mcor = ones(length(T.red_tc), length(file_names));
    for k = 1:length(file_names)
        T = readtable(file_names{k});
        mcor(:, k) = movcorr(T.bulk_green, T.green_donut, 450);
    end    
    cor_pairs = corr(mcor);
    within_group_correlations_des = [within_group_correlations_des cor_pairs(find(tril(cor_pairs, -1) ~= 0))'];
end

% LOOK AT BETWEEN CELL CORRELATION
D = 'C:\*PrL*\*GLM predictors*\*reg_*sal*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('C:\PrL2.1')); addpath(genpath('C:\PrL3.2')); addpath(genpath('C:\PrL3.4')); addpath(genpath('C:\samira_nejrgeco'));
all_names = natsortfiles({files.name});

within_cell_correlations_sal = [];
run_names = cellfun(@(x)x(1:11), all_names, 'UniformOutput', false);
unique_run_names = unique(run_names);
mice = {'PrL2.1', 'PrL3.2', 'PrL3.4'};


for k = 1:length(mice)
    mouse_runs = unique_run_names(contains(unique_run_names, mice{k}));
    mouse_files = all_names(contains(all_names, mouse_runs));
    for i = 1:length(all_names(contains(all_names, mouse_runs{1})))
        cellnum = num2str(i);
        cell_files = mouse_files(contains(mouse_files, ['cell' cellnum]));
        mcor = zeros(4200, length(cell_files));
        for z = 1:length(cell_files)
            T = readtable(cell_files{z});
            mcor(:,z) = T.mcor_bglg(1:4200);
        end
        cor_pairs = corr(mcor);
        within_cell_correlations_sal = [within_cell_correlations_sal cor_pairs(find(tril(cor_pairs, -1) ~= 0))'];
    end
end

D = 'C:\*PrL*\*GLM predictors*\*reg_*des*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('C:\PrL2.1')); addpath(genpath('C:\PrL3.2')); addpath(genpath('C:\PrL3.4')); addpath(genpath('C:\samira_nejrgeco'));
all_names = natsortfiles({files.name});

within_cell_correlations_des = [];
run_names = cellfun(@(x)x(1:11), all_names, 'UniformOutput', false);
unique_run_names = unique(run_names);
mice = {'PrL2.1', 'PrL3.2', 'PrL3.4'};


for k = 1:length(mice)
    mouse_runs = unique_run_names(contains(unique_run_names, mice{k}));
    mouse_files = all_names(contains(all_names, mouse_runs));
    for i = 1:length(all_names(contains(all_names, mouse_runs{1})))
        cellnum = num2str(i);
        cell_files = mouse_files(contains(mouse_files, ['cell' cellnum]));
        mcor = zeros(4200, length(cell_files));
        for z = 1:length(cell_files)
            T = readtable(cell_files{z});
            mcor(:,z) = T.mcor_bglg(1:4200);
        end
        cor_pairs = corr(mcor);
        within_cell_correlations_des = [within_cell_correlations_des cor_pairs(find(tril(cor_pairs, -1) ~= 0))'];
    end
end

%%
figure;
[f_des,x_des] = ecdf(within_group_correlations_des);
[f_sal,x_sal] = ecdf(within_group_correlations_sal);
plot(x_des, f_des, 'LineWidth', 2, 'Color', 'k'); 
hold on; 
plot(x_sal, f_sal, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '--');
[~,p] = kstest2(within_group_correlations_des, within_group_correlations_sal);
% if (p < 0.001)
%     text(-0.1, 0.7, 'KS test: ***', 'FontSize', 20);
% end
title('CDF of NE Synchrony Correlations');
legend({'desipramine', 'saline'}, 'FontSize', 20, 'Location', 'northwest');
legend boxoff
ax = gca;
% ax.LineWidth = 3;
yticks([0, 0.25, 0.50, 0.75, 1]);
ylabel('Cumulative fraction');
xlabel('Correlation');
xticks([-0.5, 0, 0.5, 1]);
xlim([-.5, 1])
% ax.FontSize = 20;
box off;
axis square;

%%
figure;
h1 = histogram(within_group_correlations_sal, 'FaceColor', [1 0 0], 'FaceAlpha', 0.2); hold on;
h2 = histogram(within_group_correlations_des, 'FaceColor', [0 0 1], 'FaceAlpha', 0.2); hold on;
h3 = histogram(within_cell_correlations_sal, 'FaceColor', [0 1 0], 'FaceAlpha', 0.2); hold on;
h4 = histogram(within_cell_correlations_des, 'FaceColor', [1 0 1], 'FaceAlpha', 0.2);
legend({'saline, within run', 'desipramine, within run', 'saline, within cell', 'desipramine, within cell'});
title({'Saline vs. Desipramine within run and within cell correlations'; 'local to bulk green correlations'});
xlabel('correlation')

%% LOOP THROUGH ALL CORRELATION PLOTS
clear;
D = 'C:\*PrL*\*GLM predictors*\*reg_*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('C:\PrL2.1')); addpath(genpath('C:\PrL3.2')); addpath(genpath('C:\PrL3.4')); addpath(genpath('C:\samira_nejrgeco'));
all_names = natsortfiles({files.name});

% f = figure;
% plot(50, 50, 'bx');
% export_fig C:\samira_nejrgeco\lg_bg_relationship\mcor_by_run_condition.pdf
% close(f);

run_names = cellfun(@(x)x(1:11), all_names, 'UniformOutput', false);
run_cutoffs = [];
unique_run_names = unique(run_names);

for k = 6 %1:length(unique_run_names)
run_names = all_names(contains(all_names, unique_run_names{k}));
mcor1 = ones(4200, length(run_names));
mcor2 = ones(4200, length(run_names));
for i = 1:length(run_names)
    T = readtable(run_names{i});
    mcor1(:, i) = movcorr(T.bulk_green(1:4200), T.green_donut(1:4200), 450);
    mcor2(:, i) = movcorr(T.bulk_green(1:4200), T.red_tc(1:4200), 450);
end
f = figure;
f.Position = [100, 100, 1500, 500];
% if contains(unique_run_names{k}, 'sal')
    plot(mcor1, 'LineWidth', 0.2, 'Color', [0.5, 0.5, 0.5, 0.5]);
    hold on; 
    plot(mean(mcor1, 2), 'LineWidth', 3, 'Color', [0 0 0]);
    ylabel('Correlation');
% elseif contains(unique_run_names{k}, 'des')
%     plot(mcor1, 'LineWidth', 0.5, 'Color', [255,192,203]/256);
%     hold on; plot(mean(mcor1, 2), 'LineWidth', 3, 'Color', [1 0 0]);
%     ylabel('correlation');
% end
title(unique_run_names{k}, 'Interpreter', 'none');
ax = gca(f);
% ax.LineWidth = 3;
% ax.FontSize = 25;
xticks(0:900:4200);
xticklabels({'', '1', '2', '3', '4', '5'});xlabel('Time (minutes)');
yticks([0, 0.25, 0.50, 0.75, 1.00]);
box off
xlim([0, 4200]);
title('Sliding Correlation between local and global NE');
% export_fig C:\samira_nejrgeco\lg_bg_relationship\mcor_by_run_condition1.pdf
% append_pdfs('C:\samira_nejrgeco\lg_bg_relationship\mcor_by_run_condition.pdf', 'C:\samira_nejrgeco\lg_bg_relationship\mcor_by_run_condition1.pdf');
% close(f);
end