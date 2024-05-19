clear;
f = figure;
f.Position = [100, 100, 900, 900]
T = readtable('H:\nejrgeco_code_for_publication\figure 1\Innervation_density_data.xlsx');
PFC = T.ratio;
V1 = T.ratioV1;
X = [PFC V1];
X(find(PFC == max(X)), :) = [];
E = std(X)/sqrt(size(X,1));
Xmean = mean(X);
superbar(Xmean, 'E', E, 'BarFaceColor', [0.5 0.5 0.5; 0.7 0.7 0.7]);
hold on;
scatter([1+((randi([-10, 10], size(X, 1), 1))/100); 2+((randi([-10, 10], size(X, 1), 1))/100)], [X(:, 1); X(:, 2)], 'k', 'filled');
xticks([1, 2]);
xticklabels({'PFC', 'V1'});
ax = gca;
ax.FontSize = 30;
ax.LineWidth = 6;
ylabel('total axon area (%)'); title('LC Innervation Density in PFC vs. V1');

%% get stats
[h, p, ci, stats] = ttest2(X(:,1), X(:,2))
