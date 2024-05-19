%% load in data
clear;

% paths
addpath(genpath('H:\nejrgeco_code_for_publication\'));
addpath(genpath('H:\PrL2.1')); addpath(genpath('H:\PrL3.2')); addpath(genpath('H:\PrL3.4'));

%%
% set x vector
x1 = 1:1:100;
[all_y,all_b_sal,all_p_sal] = logistic_regression(x1,'sal','Blue');
figure;
plot_log(x1, all_y, 'saline', 'Release'); hold on;
[all_y,all_b_des,all_p_des] = logistic_regression(x1,'des','Blue');
plot_log(x1, all_y, 'desipramine', 'Release');

%% get stats
disp("proportion significant saline, reuptake")
(sum(all_p_sal<0.05))/length(all_p_sal)
disp("proportion significant desipramine, reuptake")
(sum(all_p_des<0.05))/length(all_p_des)
%%
% set colors
h = get(gca,'Children');
h1 = h(1); set(h1, 'FaceColor', 'm');
h2 = h(2); set(h2, 'Color', 'm');
h3 = h(3); set(h3, 'FaceColor', 'b');
h4 = h(4); set(h4, 'Color', 'b');
set(gca,'TickDir','out');
%% compare saline and desipramine
% set x vector
x1 = 1:1:100;
Xb = [];
Xp = [];

% saline, release
[all_y,all_b,all_p] = logistic_regression(x1,'sal','Yellow');
Xb = [Xb {all_b'}]; Xp = [Xp {all_p'}];
plot_log(x1, all_y, 'saline', 'Release');

% saline, reuptake
[all_y,all_b,all_p] = logistic_regression(x1,'sal','Blue');
Xb = [Xb {all_b'}]; Xp = [Xp {all_p'}];
plot_log(x1, all_y, 'saline', 'Reuptake');

% desipramine, release
[all_y,all_b,all_p] = logistic_regression(x1,'des','Yellow');
Xb = [Xb {all_b'}]; Xp = [Xp {all_p'}];
plot_log(x1, all_y, 'desipramine', 'Release');

% desipramine, reuptake
[all_y,all_b,all_p] = logistic_regression(x1,'des','Blue');
Xb = [Xb {all_b'}]; Xp = [Xp {all_p'}];
plot_log(x1, all_y, 'desipramine', 'Reuptake');

%%
% reorder Xb
Xb_new = Xb([1 3 2 4])
figure;
% set up data for comparison
E = []; % get error
for i = 1:size(Xb_new, 2)
    [~,~,ci,~] = ttest(Xb_new{i});
    E = [E ci];
end
error =  E(2,:)-E(1,:);
P = NaN(4, 4); % get p matrix
pairs =  [1  2; 3 4];
for i = 1:size(pairs,1)
    [~,p, ~,~] = ttest2(Xb_new{pairs(i, 1)}, Xb_new{pairs(i, 2)});
    P(pairs(i, 1), pairs(i, 2)) = p; P(pairs(i, 2), pairs(i, 1)) = p; 
end
Xmean = cellfun(@mean, Xb_new);

% plot
superbar(Xmean, 'E', error([1 3 2 4]), 'P', P,...
    'BarFaceColor', [0 0 1; 1 0 1; 0 0 1; 1 0 1]);

% figure properties
xticks([1.5, 3.5]);
xticklabels({'Release beta', 'Reuptake beta'});
yticks([-0.02, 0, 0.3]);
xtickangle(45);
legend({'saline', 'desipramine'});
title('Saline vs Desipramine Predictor Distributions');


%% calculate how often b is significant
sum(all_p<0.05)/length(all_p)
mean(all_b(all_p<0.05))

%% PLOTTING FUNCTION
% drug = 'saline' or 'desipramine'
% condition = 'Release' or 'Reuptake'

function plot_log(x1, all_y, drug, condition)
    plot(x1, mean(all_y), 'k', 'LineWidth', 2);
    bci = bootci(300, @mean, all_y);
    hold on;
    patch([x1 fliplr(x1)], [bci(1, :) fliplr(bci(2, :))], 'k', 'FaceAlpha',0.2, 'EdgeColor','none')
    
    % figure properties
    xlim([1,max(x1)]);
    ylim([0, 1])
    title(['NE synchrony vs. Distance to ' condition]);
    subtitle(drug);
    xlabel('Distance to Release (microns)');
    ylabel('NE synchrony');
    xticks([0, 25, 50, 75, 100]);
    yticks([0.25, 0.5, 0.75, 1]);
    ax = gca;
    ax.FontSize = 25;
    box off;
    axis square;
end
%% LOGISTIC REGRESSION FUNCTION
% drug = 'sal' or 'des'
% condition = 'Yellow' or 'Blue'

function [all_y,all_b,all_p] = logistic_regression(x1,drug,condition)
    
    % load data
    D = ['H:\*PrL*\*GLM predictors*\*reg2*' drug '*'];
    files = dir(fullfile(D, '*.csv'));
    names = cell(size(files,1), 1);
    for i = 1:length(names)
        names{i} = [files(i).folder '\' files(i).name];
    end

    % set empty vectors to fill
    all_y = [];
    all_b = [];
    all_p = [];

    % loop through all local fields
    for i = 1:length(names)
        T = readtable(names{i});
        eval(['x = T.dSpreading' condition ';']);
        y = T.mcor_bglg;
        
        dwnsmp = 1:15:length(x);
        mdl = fitlm(log(x(dwnsmp)), y(dwnsmp));
        coeff = mdl.Coefficients;
        b = coeff.Estimate;
        p = coeff.pValue;
        y1 = b(1)+ b(2)*log(x1);
        all_p = [all_p p(2)];
        all_b = [all_b b(2)];
        all_y = [all_y; y1];
    end
end
