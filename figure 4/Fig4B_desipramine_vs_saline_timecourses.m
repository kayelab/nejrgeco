% DESIPRAMINE VS. SALINE RAW DATA EXAMPLE

% load in file paths
addpath(genpath('F:\PrL2.1')); addpath(genpath('F:\PrL3.2')); addpath(genpath('F:\PrL3.4')); addpath(genpath('F:\samira_nejrgeco'));

% load in colors
rgb_ne_global   = [40, 181, 156]/255;
rgb_ne_local    = [1, 88, 97]/255; 
rgb_ca_global   = [253, 59, 34]/255;
rgb_ca_local    = [148, 36, 36]/255;

% load in saline data
D = 'F:\*PrL*\*GLM predictors*\*reg2*sal*';
files = dir(fullfile(D, '*.csv'));
names = natsortfiles({files.name});
T = readtable(names{50}); %50

% load in 1 non-detrended desipramine global NE timecourse
D = 'F:\*PrL*\*data*\';
files = dir(fullfile(D, '*roi_struct*des*.mat'));
all_names = natsortfiles({files.name});
load(all_names{1});
green = reshape(cell2mat({roi_data.tc_g_raw}), 4200, [])';

% load in desipramine data
D = 'F:\*PrL*\*GLM predictors*\*reg2*des*';
files = dir(fullfile(D, '*.csv'));
names = natsortfiles({files.name});
T2 = readtable(names{90}); %14 %40 #90 %180 %260

% plot 
f = figure;
h1 = plot(T.green_donut + 4, 'Color', rgb_ne_local, 'LineWidth', 2); hold on; 
plot(T.bulk_green + 8, 'Color', rgb_ne_local, 'LineWidth', 2); hold on;
h2 = plot(T2.green_donut + 12, 'Color', rgb_ca_global, 'LineWidth', 2); hold on;
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
legend([h1, h2], {'saline', 'desipramine'});
sgtitle('Example GRAB_NE Timecourses', 'FontSize', 30);