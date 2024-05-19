clear;

% get data
addpath(genpath('F:\PrL2.1')); addpath(genpath('F:\PrL3.2')); addpath(genpath('F:\PrL3.4')); addpath(genpath('F:\nejrgeco_code_forpublication\'));
D = 'F:\*PrL*\*GLM predictors*\*reg2*sal*';
files = dir(fullfile(D, '*.csv'));
names = cell(size(files,1), 1);
for i = 1:length(names)
    names{i} = [files(i).folder '\' files(i).name];
end

% define key variables
timelags = -15:15:15
timelagsBlue = cell(length(timelags), 1);
timelagsYellow = cell(length(timelags),1);
pBlue = cell(length(timelags), 1);
pYellow = cell(length(timelags), 1);

% fit model across timelags ( code fits two models)
% model 1: dSpreadingBlue ~ synchrony 
% model 2: dSpreadingYellow ~ synchrony
for k = 1:length(timelags) % iterate across time lags
    k
    b_fit1 = [];
    b_fit2 = [];
    p_fit1 = [];
    p_fit2 = [];
    
    for i = 1:length(names) % iterate across cells
        % load cell data
        T = readtable(names{i});
        
        % fit model 1 ( dSpreadingBlue ~ synchrony )
        if (timelags(k) < 0)
            X = [T.dSpreadingBlue(abs(timelags(k))+1:length(T.dSpreadingBlue))];
            [b, dev, stats] = glmfit( X, T.mcor_bglg(1:(length(T.mcor_bglg)-abs(timelags(k)))) );
        elseif (timelags(k) >= 0)
            X = [T.dSpreadingBlue(1:(length(T.dSpreadingBlue)-abs(timelags(k))))];
            [b, dev, stats] = glmfit( X, T.mcor_bglg(abs(timelags(k))+1:length(T.mcor_bglg)) );
        end
        b_fit1(i) = b(2);
        p_val = stats.p;
        p_fit1(i) = p_val(2);
        
        % fit model 1 ( dSpreadingYellow ~ synchrony )
        if (timelags(k) < 0)
            X = [T.dSpreadingYellow(abs(timelags(k))+1:length(T.dSpreadingBlue))];
            [b, dev, stats] = glmfit( X, T.mcor_bglg(1:(length(T.mcor_bglg)-abs(timelags(k)))) );
        elseif (timelags(k) >= 0)
            X = [T.dSpreadingYellow(1:(length(T.dSpreadingYellow)-abs(timelags(k))))];
            [b, dev, stats] = glmfit( X, T.mcor_bglg(abs(timelags(k))+1:length(T.mcor_bglg)) );
        end
        b_fit2(i) = b(2);
        p_val = stats.p;
        p_fit2(i) = p_val(2);
    end

    timelagsBlue{k} = b_fit1;
    timelagsYellow{k} = b_fit2;
    pBlue{k} = p_fit1;
    pYellow{k} = p_fit2;

end

%% TIMELAGS PROCESSING
% remove rows that aren't significant at zero timelag (46)
timelagsBlueMat = cell2mat(timelagsBlue)';
timelagsYellowMat = cell2mat(timelagsYellow)';
pBlueMat = cell2mat(pBlue)';
pYellowMat = cell2mat(pYellow)';

% remove values that aren't significant
timelagsBlueMat(find(pBlueMat >= 0.05)) = NaN;
timelagsYellowMat(find(pYellowMat >= 0.05)) = NaN;

timelagsBlueMat(find(pBlueMat(:,round(length(timelags/2)) ) < 0.05), :) = NaN;
timelagsYellowMat(find(pYellowMat(:,round(length(timelags)/2)) < 0.05), :) = NaN;

timelagsYellow2 = mat2cell(timelagsYellowMat, size(timelagsYellowMat,1), ones(size(timelagsYellowMat,2), 1) )';
timelagsBlue2 = mat2cell(timelagsBlueMat, size(timelagsBlueMat,1), ones(size(timelagsBlueMat,2), 1) )';

%% MAKE SURE YOU ARE USING SALINE DATA
timelagsBlueMat = cell2mat(timelagsBlue)';
timelagsYellowMat = cell2mat(timelagsYellow)';

f = figure;
X = [];
for i = [ 3 2 1] % these are in reverse order so that the sign of the timelags is reverse (current is: positive timelag = preceding)
    X = [X timelagsBlueMat(:,i) timelagsYellowMat(:,i)];
end
Xmean = nanmean(X,1);
E = [];
for i = 1:length(Xmean)
    [~,~,ci,~] = ttest(X(:,i));
    E = [E ci];
end
P = NaN(size(Xmean));
% for i = 1:length(Xmean)
% [~,p, ~,~] = ttest(X(:,i));
% P(i) = p; 
% end
% superbar(Xmean, 'E', E(2,:)-E(1,:), 'P', P, 'BarFaceColor', repmat([0 0 1; 255/256 215/256 0], 3, 1));
superbar(Xmean, 'E', E(2,:)-E(1,:), 'BarFaceColor', repmat([0, 0, 0; .66, .66, .66], 3, 1));
xticks(1.5:2:5.5);
xticklabels({'1s preceding', 'no timelag', '1s lagging'});
xtickangle(30);
ylabel('Beta Weight');
legend({'growing reuptake', 'growing release'}, 'FontSize', 24);
legend boxoff
title('Beta Weights Across Timelags', 'FontSize', 40);
ylim([0, 0.012]);
% yticks([0, 0.012]);
% ax = gca(f);
% ax.FontSize = 25;
% set(gca, 'linewidth', 6);
axis square

%% PLOTTING RESULTS
f = figure;
f.Position = [100 100 700 900];
subplot(3,1,1);
meanBlue = cellfun(@nanmean, timelagsBlue2);
meanYellow = cellfun(@nanmean, timelagsYellow2);
[~,~,ciBlue, ~] = cellfun(@ttest, timelagsBlue2, 'UniformOutput', false); ciBlue = cell2mat(ciBlue')';
[~,~,ciYellow, ~] = cellfun(@ttest, timelagsYellow2, 'UniformOutput', false); ciYellow = cell2mat(ciYellow')';
errorbar(timelags, meanBlue, ciBlue(:,2), ciBlue(:,1), 'b-o', 'LineWidth', 3, 'CapSize', 12); hold on;
errorbar(timelags, meanYellow, ciYellow(:,2), ciYellow(:,1), '-o', 'Color', [255, 215, 0]/255, 'LineWidth', 3, 'CapSize', 12);
legend({'dSpreadingBlue', 'dSpreadingYellow'});
yline(0);
ylim([-0.005, 0.015]);
title('beta of dSpreadingBlue and dSpreadingYellow Across Timelags');
ylabel('beta');
xlabel('frames');

subplot(3,1,2);
errorbar(timelags, meanBlue, ciBlue(:,2), ciBlue(:,1), 'b-o', 'LineWidth', 3, 'CapSize', 12); 
yline(0);
ylim([-0.005, 0.015]);
title('beta of dSpreadingBlue Across Timelags');
ylabel('beta');
xlabel('frames');

subplot(3,1,3);
errorbar(timelags, meanYellow, ciYellow(:,2), ciYellow(:,1), '-o', 'Color', [255, 215, 0]/255, 'LineWidth', 3, 'CapSize', 12);
yline(0);
ylim([-0.005, 0.015]);
title('beta of dSpreadingYellow Across Timelags');
ylabel('beta');
xlabel('frames');


%% T-TESTING
ttestArray = zeros(length(timelags), length(timelags));
for i = 1:length(timelags)
    for k = 1:length(timelags)
        ttestArray(i,k) = ttest2(timelagsBlue{i}, timelagsBlue{k});
    end
end
imagesc(ttestArray);
%%
figure;
b_fit1_sig = b_fit1(find(p_fit1 < 0.05 & p_fit2 < 0.05));
b_fit2_sig = b_fit2(find(p_fit1 < 0.05 & p_fit2 < 0.05));
names_sig = names(find(p_fit1 < 0.05 & p_fit2 < 0.05));
run_names_sig = cellfun(@(x)x(1:11), names_sig, 'UniformOutput', false);
g = gscatter(b_fit1_sig', b_fit2_sig', run_names_sig',[], '.', 30);
yline(0); xline(0);
hold on;
fit = fitlm(table(b_fit1_sig', b_fit2_sig'));
coef = fit.Coefficients;
int = coef(1,1); int = int.Estimate;
slope = coef(2,1); slope = slope.Estimate;
text(0.01, 0.015, sprintf('slope = %f', slope));
x = -0.03: 0.01:0.03;
plot(x, x*slope + int, 'LineWidth', 3, 'Color', 'k');
xlabel('Beta dSpreadingBlue'); ylabel('Beta dSpreadingYellow');
title({'Relationship Between Beta dSpreadingBlue and Beta dSpreadingYellow'; ...
    'all regressed desipramine runs'});
group = g(3);
group.Color = [0.7 0.7 0.5];
group = g(4);
group.Color = [34 139 34]/256;
group = g(8);
group.Color = [ 135 206 250]/256;
group = g(12);
group.Color = [219 112 147]/256;