% COMPARE TEMPORAL AND SPATIAL CORRELATION OF RELEASE AND REUPTAKE

% set file paths
addpath(genpath('H:\nejrgeco_code_for_publication\'));
addpath(genpath('H:\nejrgeco_code_for_publication\figure 5\NeuroPattToolbox'));

% load in control structure
letter_codes = {'R', 'P', 'A'};
load('H:\nejrgeco_code_for_publication\figure 5\controlstruct.mat');
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
addpath(genpath('H:\PrL2.1')); addpath(genpath('H:\PrL3.2')); addpath(genpath('H:\PrL3.4'));
D1 = 'H:\*PrL*\*PrL* tif images\';
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
addpath(genpath(['H:\' mouse]));

% flow file
flow_file = dir(fullfile(['H:\' mouse '\' mouse ' optic flow data\'], ['*' condition(1:3) run  '*' ]));
flow_file = flow_file.name;
load(flow_file);

% GRABNE file
fluorescence_file = dir(fullfile(['H:\' mouse '\' mouse ' tif images\'], ['*' condition '*run' run  '*' ])); 
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

%% compare temporal and spatial cor FIGURE 5C
% sal_tcor = temporal_cormat(5:8, 5:8, 2:size(temporal_cormat,3));
% sal_scor = spatial_cormat(5:8, 5:8, 2:size(spatial_cormat,3));

sal_tcor = temporal_cormat(5:8, 5:8, 1:25);
sal_scor = spatial_cormat(5:8, 5:8, 1:25);

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
title('Desipramine Spatial & Temporal Correlation', 'FontSize', 24);
ax = gca(f);
% ax.FontSize = 25;
ylim([-1, 1]); yticks([-1:.5:1]); 
% set(gca, 'linewidth', 6);
axis square;

%% get statistics
%Xmean = nanmean(X,1);
disp("mean temporal correlations")
mean(X(:,1))
disp("SE temporal correlations")
std(X(:,1))./sqrt(length(X(:,1)))
disp("1-sample t-test temporal correlations")
[h,p,ci,stats] = ttest(X(:,1))
disp("mean spatial correlations")
mean(X(:,2))
disp("SE spatial correlations")
std(X(:,2))./sqrt(length(X(:,2)))
disp("1-sample t-test spatial correlations")
[h,p,ci,stats] = ttest(X(:,2))