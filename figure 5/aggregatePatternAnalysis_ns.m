%% SPATIAL DISTRIBUTION AND SPATIAL CORRELATION MATRIX
%%BUILT IN FUNCTION FOR PLOTTING SPATIAL DISTRIBUTION
%% setup

% paths
addpath('F:\cdrive-samira_nejrgeco\');
setup_figprop;

% colors
rgb_ne_global   = [40, 181, 156]/255;
rgb_ne_local    = [1, 88, 97]/255; 
rgb_ca_global   = [253, 59, 34]/255;
rgb_ca_local    = [148, 36, 36]/255;

%% PLOT RAW DATA (TEMPORAL)
condition = 'saline'
mouse = 'PrL3.2'
run = '2'
addpath(genpath(['F:\' mouse]));
addpath(genpath('F:\samira_nejrgeco'));
flow_file = dir(fullfile(['F:\' mouse '\' mouse ' optic flow data\'], ['*' condition(1:3) run  '*' ]));
flow_file = flow_file.name;
load(flow_file);

%%
allPatternLocs = patternDetection.allPatternLocs;
Vy = patternDetection.Vy;
Vx = patternDetection.Vx;
pattTypes = patternDetection.pattTypes;
activeArray = zeros(9, size(Vx,3));
for i = 1:length(pattTypes)
    if ~isempty(allPatternLocs{i})
        time = allPatternLocs{i}(:,3);
        for z = 1:size(Vx, 3)
        activeArray(i,z) = sum(time == z);
        end
        activeArray(i,:) = smooth(activeArray(i,:), 75/4200);
    end
end

D = ['F:\*' mouse '*\*GLM predictors*\*reg3_*' condition(1:3) run '*'];
files = dir(fullfile(D, '*.csv'));
addpath(genpath('C:\PrL2.1')); addpath(genpath('C:\PrL3.2')); addpath(genpath('C:\PrL3.4')); addpath(genpath('C:\samira_nejrgeco'));
names = natsortfiles({files.name});
T = readtable(names{20});

f = figure;
f.Position = [100, 100, 1700, 900];
subplot(2,1,1); 
plot(activeArray(5,:), 'LineWidth', 2, 'Color', [255/256, 215/256 0]); 
hold on; 
plot(activeArray(6,:), 'LineWidth', 2, 'Color', [0 0 1]); 
legend({'release', 'reuptake'}, 'LineWidth', 3); 
xlim([1000, 3000]);
ax = gca;
ax.XAxis.Visible = 'off'; 
box off;
% ax.LineWidth = 3;
yticks(15); yticklabels('total patterns');
% ax.FontSize = 25;

dSpreadingYellow = smooth(T.dSpreadingYellowRaw, 0.05);
dSpreadingBlue = smooth(T.dSpreadingBlueRaw, 0.05);
b1 = activeArray(5, :) > median(activeArray(5, :));
b2 = activeArray(6, :) > median(activeArray(6, :));
% b1 = dSpreadingYellow >= 25;
% b2 = dSpreadingBlue >= 25;
subplot(2,1,2); 
colormap(flipud(hot));
imagesc([b1; b2]);
xlim([1000, 3000]);
box off
ax = gca; 
%ax.FontSize = 25;
xticks(1000:450:3000); 
xticklabels({'', '1', '2', '3', '4', ''});
xlabel('Time (minutes)');
%ax.LineWidth = 3;
yticks([1, 2]); yticklabels({'Release', 'Reuptake'});
sgtitle('Temporal Distribution of Release and Reuptake', 'FontSize', 30);
%%
%load in control structure
% controls
letter_codes = {'R', 'P', 'A'};
load('F:/cdrive-samira_nejrgeco/Optical Flow and Pattern Analysis/controlstruct.mat');
eval(['control_spatial_corr = {controlstruct.' letter_codes{1} 'sCorrMat};'])
eval(['control_temp_corr = {controlstruct.' letter_codes{1} 'tCorrMat};'])
CsCorrMat = control_spatial_corr{1};
CtCorrMat = control_temp_corr{1};
for i = 2:length(controlstruct)
    CsCorrMat = cat(3, CsCorrMat, control_spatial_corr{i});
    CtCorrMat = cat(3, CtCorrMat, control_temp_corr{i});
end

%set up runs to loop through
conditions = [repmat({'saline'}, 1, 12) repmat({'desipramine'}, 1, 12)];
mice = [repmat({'PrL2.1'}, 1, 2) repmat({'PrL3.2'}, 1, 5) repmat({'PrL3.4'}, 1, 5) repmat({'PrL2.1'}, 1, 2) repmat({'PrL3.2'}, 1, 4) repmat({'PrL3.4'}, 1, 6) ];
runs = {'2', '4', '1', '2', '3', '4', '5', '1', '2', '3', '4', '5', '1', '2', '1', '2', '3', '4', '1', '2', '3', '4', '5', '6'};

%get names of GRABNE fluorescence recordings
addpath(genpath('F:\PrL2.1')); addpath(genpath('F:\PrL3.2')); addpath(genpath('F:\PrL3.4'));
D1 = 'F:\*PrL*\*PrL* tif images\';
files = dir(fullfile(D1,'*GRABNE*.tif'));
names = {files.name};
names = natsortfiles(names);
names(contains(names, '10minB')) = [];

% set empty arrays to fill
spatial_cormat = zeros(9,9); % the temporal correlation matrix of all 9 patterns to each other
temporal_cormat = zeros(9,9); % the spatial correlation matrix of all 9 patterns to each other

for z = 1:24
condition = conditions{z};
mouse = mice{z};
run = runs{z}
addpath(genpath(['F:\' mouse]));
addpath(genpath('F:\cdrive-samira_nejrgeco'));

% flow file
flow_file = dir(fullfile(['F:\' mouse '\' mouse ' optic flow data\'], ['*' condition(1:3) run  '*' ]));
flow_file = flow_file.name;
load(flow_file);

% GRABNE file
fluorescence_file = dir(fullfile(['F:\' mouse '\' mouse ' tif images\'], ['*' condition '*run' run  '*' ])); 
fluorescence_file = fluorescence_file.name;
sframe = 301;
[gr, ~, ~, ~] = loadrdz(fluorescence_file, sframe);
filtered_green = double(gr);
filtered_green = imgaussfilt(mean(filtered_green,3), 4);
Z = reshape(filtered_green, [], 1);

% get the distribution of patterns over time 
allPatternLocs = patternDetection.allPatternLocs;
Vy = patternDetection.Vy;
Vx = patternDetection.Vx;
pattTypes = patternDetection.pattTypes;
activeArray = zeros(9, size(Vx,3));
for i = 1:length(pattTypes)
    if ~isempty(allPatternLocs{i})
        time = allPatternLocs{i}(:,3);
        for w = 1:size(Vx, 3)
        activeArray(i,w) = sum(time == w);
        end
    end
end

% get the distribution of patterns over space
ind = 1:length(allPatternLocs);
distribution = cell(size(allPatternLocs));
spatial_correlation = zeros(size(Vx,1)*size(Vx,2), length(allPatternLocs));
for k = ind(~cellfun('isempty',allPatternLocs))
    spatial_distribution = zeros(size(Vx,1), size(Vx,2));
    x = round(allPatternLocs{k}(:,1));
    y = round(allPatternLocs{k}(:,2));
    for i = 1:length(x)
    spatial_distribution(x(i),y(i)) = spatial_distribution(x(i),y(i)) + 1;
    end
    spatial_distribution = imgaussfilt(spatial_distribution, 4);
    distribution{k} = spatial_distribution;
    spatial_correlation(:,k) = reshape(spatial_distribution, size(Vx,1)*size(Vx,2), 1); %linearize the spatial distribution of pattens
end

%SPATIAL CORRELATION
% take the partial correlation to account for the patterns' correlation to
% differences in GRABNE expression
corrmat = partialcorr(spatial_correlation, Z);
spatial_cormat = cat(3, spatial_cormat, corrmat);

% TEMPORAL CORRELATION
% temporally smooth the array of patterns over time (5 s kernel)
activeArraySmooth = ones(size(activeArray));
for i = 1:length(pattTypes)
activeArraySmooth(i,:) = smooth(activeArray(i,:), 75/4200);
end
corrmat = corrcoef(activeArraySmooth');
temporal_cormat = cat(3, temporal_cormat, corrmat);

end

%%
sal_tcor = temporal_cormat(3:9, 3:9, 2:size(temporal_cormat,3));
sal_scor = spatial_cormat(3:9, 3:9, 2:size(spatial_cormat,3));

%%
des_tcor = temporal_cormat(3:6, 3:6, 2:size(temporal_cormat,3));
des_scor = spatial_cormat(3:6, 3:6, 2:size(spatial_cormat,3));

%%
sal_tcor = cat(3, sal_tcor, NaN(4,4));
sal_scor = cat(3, sal_scor, NaN(4,4));

%%
[r, c] = ind2sub([4, 4], find(tril(ones(4,4), -1) ~= 0));
X = [];
xlabels = [];
pattCompare = pattTypes(3:6);
for i = 1:length(r)
    X = [X squeeze(sal_tcor(r(i), c(i), :)) squeeze(des_tcor(r(i), c(i), :)) squeeze(sal_scor(r(i), c(i), :)) squeeze(des_scor(r(i), c(i), :))];
    xlabels = [xlabels {[pattCompare{r(i)} ':' pattCompare{c(i)}]}];
end
Xmean = nanmean(X,1);
E = [];
for i = 1:length(Xmean)
    [~,~,ci,~] = ttest(X(:,i));
    E = [E ci];
end
P = NaN(24, 24);
pairs =  [1  2; 3 4; 5 6; 7 8; 9 10; 11 12; 13 14; 15 16; 17 18; 19 20; 21 22; 23 24];
for i = 1:size(pairs,1)
[~,p, ~,~] = ttest(X(:,pairs(i, 1)), X(:,pairs(i, 2)));
P(pairs(i, 1), pairs(i, 2)) = p; P(pairs(i, 2), pairs(i, 1)) = p; 
end
superbar(Xmean, 'E', E(2,:)-E(1,:), 'P', P, 'BarFaceColor', repmat([.8 .2 .8; 0.2 0.7 0.2; 0.3 0.4 0.3; 0.7 0.1 0], 6, 1));
xticks(2:4:22);
xticklabels(xlabels);
xtickangle(45);
legend({'sal tcor', 'des tcor', 'sal scor', 'des scor'});
title('saline vs. desipramine spatial and temporal coupling of different patterns, PrL3.4');

%% compare temporal and spatial cor FIGURE 5C
sal_tcor = temporal_cormat(5:8, 5:8, 2:size(temporal_cormat,3));
sal_scor = spatial_cormat(5:8, 5:8, 2:size(spatial_cormat,3));
f = figure;
hold on;
[r, c] = ind2sub([4, 4], find(tril(ones(4,4), -1) ~= 0));
X = [];
xlabels = [];
% pattCompare = pattTypes(3:6);
pattCompare = {'release', 'reuptake', 'spiral-in', 'spiral-out'};
for i = 1:length(r)
    X = [X squeeze(sal_tcor(r(i), c(i), :)) squeeze(sal_scor(r(i), c(i), :))];
    xlabels = [xlabels {[pattCompare{r(i)} ':' pattCompare{c(i)}]}];
end
Xmean = nanmean(X,1);
Xmean = Xmean(:, [1 2 11 12]);
X = X(:, [1 2 11 12]);
E = [];
for i = 1:length(Xmean)
    [~,~,ci,~] = ttest(X(:,i));
    E = [E ci];
end
P = NaN(12, 12);
pairs =  [1  2; 3 4];
% for i = 1:size(pairs,1)
% [~,p, ~,~] = ttest2(X(:,pairs(i, 1)), X(:,pairs(i, 2)));
% P(pairs(i, 1), pairs(i, 2)) = p; P(pairs(i, 2), pairs(i, 1)) = p; 
% end
P = NaN(size(Xmean));
for i = 1:length(Xmean)
[~,p, ~,~] = ttest(X(:,i));
P(i) = p; 
end
% superbar(Xmean, 'E', E(2,:)-E(1,:), 'P', P, 'BarFaceColor', repmat([0 1 1; 0.9 0.5 1], 3, 1));
superbar(Xmean, 'E', E(2,:)-E(1,:), 'BarWidth', .75, 'BarFaceColor', repmat([0.33 0.33 0.33; 0.66 0.66 0.66], 3, 1));
scatter([1+((randi([-10, 10], size(X, 1), 1))/100); 2+((randi([-10, 10], size(X, 1), 1))/100)], [X(:, 1); X(:, 2)], 'k', 'filled');
xticks([1, 2]);
xlim([.5, 2.5]);
xticklabels({'Temporal', 'Spatial'});
% xtickangle(30);
ylabel('Correlation');
% legend({'Temporal correlation', 'Spatial correlation'}, 'FontSize', 20, 'Location', 'NorthWest');
% legend boxoff
title('Release and Reuptake Spatial & Temporal Correlation', 'FontSize', 24);
ax = gca(f);
% ax.FontSize = 25;
ylim([-1, 1]); yticks([-1:.5:1]); 
% set(gca, 'linewidth', 6);
axis square;

%%
figure;
err = bootci(12, @mean, [squeeze(sal_tcor(5,6,1:12)) squeeze(temporal_cormat(5,6,1:12)) squeeze(spatial_cormat(5,6,13:24)) squeeze(temporal_cormat(5,6,13:24))]);
data = mean([squeeze(spatial_cormat(5,6,1:12)) squeeze(temporal_cormat(5,6,1:12)) squeeze(spatial_cormat(5,6,13:24)) squeeze(temporal_cormat(5,6,13:24))], 1);
bar(1:4, data); hold on;
er = errorbar(1:4, data, err(1,:),err(2,:), 'ko', 'LineWidth', 3); 
xticklabels({'spatial cor saline', 'temporal cor saline', 'spatial cor desipramine', 'temporal cor desipramine'});
title('spreading Yellow/ spreading Blue Temporal and Spatial Correlations');
ylabel('spreading Yellow/ spreading Blue correlation');
xtickangle(45);


%% SPATIAL CORMAT (ALL PATTERNS)
%spatial cormat
figure;
imagesc(mean(spatial_cormat(:,:,13:size(spatial_cormat,3)),3)); colorbar; colormap(cool);
% make matrix of p-values and cross out values that aren't significant
p_values = nan(9,9);
v = find(~isnan(corrcoef(spatial_correlation)));
for i = v'
    [r, c] = ind2sub([9,9], i);
    distribution = squeeze(CsCorrMat(r, c, :));
    if(corrmat(r,c) > 0)
        nmore = sum(distribution > corrmat(r,c));
    end
    if(corrmat(r,c) < 0)
            nmore = sum(distribution < corrmat(r,c));
    end
    nequal = sum(distribution == corrmat(r,c));
    pval = (nmore + 0.5*nequal)/length(distribution);
    p_values(r,c) = pval;
end
[r, c] = ind2sub([9,9], find(p_values >= 0.05));
for i = 1:length(r)
    x = [r(i) - 0.5; r(i) + 0.5];
    y = [c(i) + 0.5; c(i) - 0.5];
    line(x,y, 'Color', 'black', 'LineWidth', 2);
    x = [r(i) + 0.5; r(i) - 0.5];
    y = [c(i) + 0.5; c(i) - 0.5];
    line(x,y, 'Color', 'black', 'LineWidth', 2);
    clear x y
end


N = length(pattTypes);
xticks(0:1:N); yticks(0:1:N);
labels = pattTypes;
labels{5} = 'NE release';
labels{6} = 'NE reuptake';
xticklabels(cat(2, {''}, labels)); yticklabels(cat(2, {''}, labels));
xtickangle(-45);
colormap(cool);
hold on;
title('Spatial Correlation of Patterns');
xlim([4.5,8.5]);
ylim([4.5,8.5]);
plot([4.5; 6.5], [6.5; 6.5], 'k', 'LineWidth', 4); plot([4.5, 4.5], [4.5; 6.5], 'k', 'LineWidth', 4);
plot([6.5; 6.5], [4.5; 6.5], 'k', 'LineWidth', 4); plot([4.5; 6.5], [4.5, 4.5], 'k', 'LineWidth', 4);
ax = gca;
ax.FontSize = 15;
axis square;
colorbar;
caxis([-0.0, 1]);

%% TEMPORAL CORMAT (ALL PATTERNS)
%temporal cormat
figure;
imagesc(mean(temporal_cormat(:,:,13:size(temporal_cormat,3)),3)); colorbar; colormap(cool);
% make matrix of p-values and cross out values that aren't significant
p_values = nan(9,9);
v = find(~isnan(corrcoef(spatial_correlation)));
for i = v'
    [r, c] = ind2sub([9,9], i);
    distribution = squeeze(CtCorrMat(r, c, :));
    if(corrmat(r,c) > 0)
        nmore = sum(distribution > corrmat(r,c));
    end
    if(corrmat(r,c) < 0)
            nmore = sum(distribution < corrmat(r,c));
    end
    nequal = sum(distribution == corrmat(r,c));
    pval = (nmore + 0.5*nequal)/length(distribution);
    p_values(r,c) = pval;
end
[r, c] = ind2sub([9,9], find(p_values >= 0.05));
for i = 1:length(r)
    x = [r(i) - 0.5; r(i) + 0.5];
    y = [c(i) + 0.5; c(i) - 0.5];
    line(x,y, 'Color', 'black', 'LineWidth', 2);
    x = [r(i) + 0.5; r(i) - 0.5];
    y = [c(i) + 0.5; c(i) - 0.5];
    line(x,y, 'Color', 'black', 'LineWidth', 2);
    clear x y
end



N = length(pattTypes);
xticks(0:1:N); yticks(0:1:N);
labels = pattTypes;
labels{5} = 'NE release';
labels{6} = 'NE reuptake';
xticklabels(cat(2, {''}, labels)); yticklabels(cat(2, {''}, labels));
xtickangle(-45);
colormap(cool);
hold on;
title('Temporal Correlation of Patterns');
xlim([4.5,8.5]);
ylim([4.5,8.5]);
plot([4.5; 6.5], [6.5; 6.5], 'k', 'LineWidth', 4); plot([4.5, 4.5], [4.5; 6.5], 'k', 'LineWidth', 4);
plot([6.5; 6.5], [4.5; 6.5], 'k', 'LineWidth', 4); plot([4.5; 6.5], [4.5, 4.5], 'k', 'LineWidth', 4);
ax = gca;
ax.FontSize = 15;
axis square;
colorbar;
caxis([-0.0, 1]);


%%
%temporal cormat
figure;
ax1 = axes;
imagesc(ax1, mean(temporal_cormat(:,:,13:size(temporal_cormat,3)),3)); 
colormap(ax1, cool); caxis(ax1, [-0.2, 0.7]); 
axis square; hold all; 

% make matrix of p-values and cross out values that aren't significant
p_values = nan(9,9);
v = find(~isnan(corrcoef(activeArraySmooth')));
for i = v'
    [r, c] = ind2sub([9,9], i);
    distribution = squeeze(CtCorrMat(r, c, :));
    if(corrmat(r,c) > 0)
        nmore = sum(distribution > corrmat(r,c));
    end
    if(corrmat(r,c) < 0)
            nmore = sum(distribution < corrmat(r,c));
    end
    nequal = sum(distribution == corrmat(r,c));
    pval = (nmore + 0.5*nequal)/length(distribution);
    p_values(r,c) = pval;
end
[r, c] = ind2sub([9,9], find(p_values >= 0.05));
for i = 1:length(r)
    x = [r(i) - 0.5; r(i) + 0.5];
    y = [c(i) + 0.5; c(i) - 0.5];
    line(ax1, x,y, 'Color', 'black', 'LineWidth', 4);
    x = [r(i) + 0.5; r(i) - 0.5];
    y = [c(i) + 0.5; c(i) - 0.5];
    line(ax1, x,y, 'Color', 'black', 'LineWidth', 4);
    clear x y
end

pattTypes{3} = 'shrinkingRelease';
pattTypes{4} = 'shrinkingReuptake';
pattTypes{5} = 'growingRelease';
pattTypes{6} = 'growingReuptake'
title('ALL SALINE RUNS Temporal Correlation Matrix');
N = length(pattTypes);
xticks(ax1, 0:1:N); yticks(ax1, 0:1:N);
xticklabels(ax1, cat(2, {''}, pattTypes)); yticklabels(ax1, cat(2, {''}, pattTypes));
xtickangle(ax1, -45);
set(ax1, 'FontSize', 30)

tri = ones(9,9);
tri = triu(tri, 1);
ax2 = axes;
imagesc(ax2, tri, 'alphadata', tri); axis square; linkaxes([ax1, ax2]); ax2.Visible = 'off';
caxis(ax2, [0, 1]); colormap(ax2, gray); 

title('ALL DESIPRAMINE RUNS Temporal Correlation Matrix');
N = length(pattTypes);
xticks(ax2, 0:1:N); yticks(ax2, 0:1:N);
xticklabels(ax2, cat(2, {''}, pattTypes)); yticklabels(ax2, cat(2, {''}, pattTypes));
xtickangle(ax2, -45);
set(ax2, 'FontSize', 30)

hold on;
for i = 0.5:1:9.5
    yline(i, 'LineWidth', 4);
    xline(i, 'LineWidth', 4);
end

idx = find(tri == 0);
for i = 1:length(idx)
    [r,c] = ind2sub([9,9], idx(i));
    if(p_values(r,c) > 0.05)
            text(r-0.25,c, 'n.s.', 'FontSize', 30);
    elseif(p_values(r,c) < 0.05 && p_values(r,c) >= 0.01)
            text(r-0.1,c, '**', 'FontSize', 30);
    elseif(p_values(r,c) < 0.01 && p_values(r,c) >= 0.001)
            text(r-0.1,c, '***', 'FontSize', 30);
    elseif(p_values(r,c) < 0.001)
            text(r-0.3,c, '>***', 'FontSize', 30);
    end

end
xlim([2.5, 6.5]);
ylim([2.5, 6.5]);






