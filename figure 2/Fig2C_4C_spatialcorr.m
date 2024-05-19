%% PLOT SALINE AND DESIPRAMINE SPATIAL CORR

% setup
addpath('H:\nejrgeco_code_for_publication\');

% load saline data (all_corr, all_dist, and partition)
load('H:\nejrgeco_code_for_publication\figure 2\desipramine_spatialcorr.mat');

% plot 
figure;
x1 = mean(all_dist, 2)*partition*0.886;
mean_dropoff = mean(all_corr, 2);
desipramine = plot(x1, mean(all_corr, 2), 'm', 'LineWidth', 1);
bci = bootci(200, @mean, all_corr');
hold on;
patch([x1' fliplr(x1')], [bci(1, :) fliplr(bci(2, :))], 'm', 'FaceAlpha',0.2, 'EdgeColor','none');
hold on;

% load desipramine data
load('H:\nejrgeco_code_for_publication\figure 2\saline_spatialcorr.mat');

% plot 
x1 = mean(all_dist, 2)*partition*0.886;
mean_dropoff = mean(all_corr, 2);
saline = plot(x1, mean(all_corr, 2), 'b', 'LineWidth', 1);
bci = bootci(200, @mean, all_corr');
hold on;
patch([x1' fliplr(x1')], [bci(1, :) fliplr(bci(2, :))], 'b', 'FaceAlpha',0.2, 'EdgeColor','none')


% figure properties
xticks([0, 25, 50, 75, 100]);
yticks([0, round(max(bci(2,:))/2, 2), round(max(bci(2,:)), 2)]);
legend([saline desipramine],'saline','desipramine');
xlim([partition,100]);
ylim([0, round(max(bci(2,:)), 2]);
title('Spatial Autocorrelation');
xlabel('Distance (microns)');
ylabel('Autocorrelation');
box off;
axis square;
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 30;
%% PLOT ONLY SALINE
clear;

% load saline data (all_corr, all_dist, partition)
load('F:\nejrgeco_code_for_publication\figure 2\saline_spatialcorr.mat');

% plot
figure;
x1 = mean(all_dist, 2)*partition*0.886;
mean_dropoff = mean(all_corr, 2);
saline = plot(x1, mean(all_corr, 2), 'b', 'LineWidth', 1);
bci = bootci(200, @mean, all_corr');
hold on;
patch([x1' fliplr(x1')], [bci(1, :) fliplr(bci(2, :))], 'b', 'FaceAlpha', 0.2, 'EdgeColor','none');

% figure properties
xticks([0, 25, 50, 75, 100]);
yticks(0:0.03:max(bci(2,:)));
xlim([partition,100]);
ylim([0, max(bci(2,:))]);
title('Spatial Autocorrelation');
xlabel('Distance (microns)');
ylabel('Autocorrelation');
box off;
axis square;
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 30;

%% COMPARE LARGE GRID AND SMALL GRID
clear;

% load saline small grid data (all_corr, all_dist, and partition)
load('H:\nejrgeco_code_for_publication\figure 2\saline_spatialcorr.mat');

% plot 
figure;
x1 = mean(all_dist, 2)*partition*0.886;
mean_dropoff = mean(all_corr, 2);
smallgrid = plot(x1, mean(all_corr, 2), 'b', 'LineWidth', 1);
bci = bootci(200, @mean, all_corr');
hold on;
patch([x1' fliplr(x1')], [bci(1, :) fliplr(bci(2, :))], 'b', 'FaceAlpha',0.2, 'EdgeColor','none');
hold on;

% load saline large grid data (all_corr, all_dist, and partition)
load('H:\nejrgeco_code_for_publication\figure 2\saline_spatialcorr_largegrid.mat');

% plot 
x1 = mean(all_dist, 2)*partition*0.886;
mean_dropoff = mean(all_corr, 2);
largegrid = plot(x1, mean(all_corr, 2), 'm', 'LineWidth', 1);
bci = bootci(200, @mean, all_corr');
hold on;
patch([x1' fliplr(x1')], [bci(1, :) fliplr(bci(2, :))], 'm', 'FaceAlpha',0.2, 'EdgeColor','none')

% figure properties
xticks([0, 25, 50, 75, 100]);
yticks(0:0.03:max(bci(2,:)));
xlim([partition,100]);
ylim([0, max(bci(2,:))]);
title('Spatial Autocorrelation');
xlabel('Distance (microns)');
ylabel('Autocorrelation');
legend([smallgrid largegrid],'10 micron patch','20 micron patch');
box off;
axis square;
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 30;

%% CURVE FITTING

% load saline data
load('F:\nejrgeco_code_for_publication\figure 2\saline_spatialcorr.mat');

% plot
figure;
x1 = mean(all_dist, 2)*4*0.886;
mean_dropoff = mean(all_corr, 2);
f = fit(x1(2:12), mean_dropoff(2:12), 'exp2');
hold on;
saline = plot(x1, mean(all_corr, 2), 'b', 'LineWidth', 3);
bci = bootci(200, @mean, all_corr');
hold on;
patch([x1' fliplr(x1')], [bci(1, :) fliplr(bci(2, :))], 'b', 'FaceAlpha', 0.2, 'EdgeColor','none');
hold on;
fittedcurve = plot(f, 'k');
print('saline coefficients')
f

% figure properties
xticks([0, 25, 50, 75, 100]);
yticks([0, 0.03, 0.06]);
legend([saline fittedcurve],'saline','fitted curve');
xlim([4,100]);
ylim([0, 0.06]);
title('Spatial Autocorrelation');
xlabel('Distance (microns)');
ylabel('Autocorrelation');
box off;
axis square;
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 30;


