%load in red and green raw movie
clear;

% paths
addpath(genpath('F:\nejrgeco_code_for_publication'));
addpath(genpath('F:\cdrive-samira_nejrgeco')); % for loadrdz function

% set variables to iterate through 
smooth_kernel = 15;
% saline
px_to_micron = [1.47 1.47 0.74 0.74 0.74 0.74 0.74 0.74 0.74 0.74 0.74 0.74];
times = [5 5 5 5 5 5 10 5 5 5 5 10];
mice = {'PrL2.1', 'PrL2.1', 'PrL3.2', 'PrL3.2', 'PrL3.2', 'PrL3.2', 'PrL3.2', 'PrL3.4', 'PrL3.4', 'PrL3.4', 'PrL3.4', 'PrL3.4'};
runs = {'2', '4', '1', '2', '3', '4', '5', '1', '2', '3', '4', '5'}
condition = 'saline'
% desipramine
% px_to_micron = [1.47 1.47 0.74 0.74 0.74 0.74 0.74 0.74 0.74 0.74 0.74 0.74];
% times = [5 5 5 5 5 5 5 5 5 5 10 10];
% mice = [repmat({'PrL2.1'}, 1, 2) repmat({'PrL3.2'}, 1, 4) repmat({'PrL3.4'}, 1, 6)];
% runs = {'1', '2', '1', '2', '3', '4', '1', '2', '3', '4', '5', '6'};
% condition = 'desipramine'

% set up empty matrices
spatial_corr = [];
all_corr_low = [];
all_corr_high = [];
% all_corr = [];
all_dist = []; %distances should be the same, but check just to be sure

%%
for z = 1:12
    %% LOAD 3D RAW GREEN MOVIE (gr, 256x256x4200)
    % load in mouse, run, and time, and pixel conversion
    tic
    z
    time = times(z);
    mouse = mice{z}
    run = runs{z}
    px_mc = px_to_micron(z);

    % get a list of files
    addpath(genpath(['F:\' mouse]));
    jRGECO_files = dir(fullfile(['F:\' mouse '\' mouse ' tif images\'], ['*' condition '*run' run '*jRGECO*' ]));
    jRGECO_files = {jRGECO_files.name};
    GRABNE_files =  dir(fullfile(['F:\' mouse '\' mouse ' tif images\'], ['*' condition '*run' run '*GRABNE*' ]));
    GRABNE_files = {GRABNE_files.name};
    
    % load in raw green and red movie (gr, rr)
    sframe = 301;
            switch time
                case 5
                    tiffile_g  = GRABNE_files{1};
                    tiffile_r = jRGECO_files{1};
                    [gr, ~, long_gr, ~] = loadrdz(tiffile_g, sframe);
                    [rr, ~, long_rr, ~] = loadrdz(tiffile_r, sframe);
                    filtered_green = double(gr);
                case 10
                    tiffile_g  = GRABNE_files{1}
                    tiffile_g2  = GRABNE_files{2}
                    tiffile_r  = jRGECO_files{1}
                    tiffile_r2  = jRGECO_files{2}
                    [gr1, ~, long_gr1, ~] = loadrdz(tiffile_g, sframe);
                    [gr2, ~, long_gr2, ~] = loadrdz(tiffile_g2, sframe);
                    gr = cat(3, gr1, gr2); clear gr1 gr2
                    gr = double(gr);
                    long_gr = cat(2, long_gr1, long_gr2); clear long_gr1 long_gr2
    
                    [rr1, ~, long_rr1, ~] = loadrdz(tiffile_r, sframe);
                    [rr2, ~, long_rr2, ~] = loadrdz(tiffile_r2, sframe);
                    rr = cat(3, rr1, rr2); clear rr1 rr2
                    rr = double(rr);
                    long_rr = cat(2, long_rr1, long_rr2); clear long_rr1 long_rr2
            end
    
    %% PARTITION RAW GREEN MOVIE INTO GRID
    partition = 8;
    % rr_cell = mat2cell(rr, 4*ones(1,64), 4*ones(1,64), size(rr, 3));
    % gr_cell = mat2cell(gr, 4*ones(1,64), 4*ones(1,64), size(gr, 3));
    % rr_cell = cellfun(@(x)squeeze(mean(x, [1 2])), rr_cell, 'UniformOutput', false);
    % gr_cell = cellfun(@(x)squeeze(mean(x, [1 2])), gr_cell, 'UniformOutput', false);
    
    rr_cell = mat2cell(rr, partition*ones(1,256/partition), partition*ones(1,256/partition), size(rr, 3));
    gr_cell = mat2cell(gr, partition*ones(1,256/partition), partition*ones(1,256/partition), size(gr, 3));
    rr_cell = cellfun(@(x)squeeze(mean(x, [1 2])), rr_cell, 'UniformOutput', false);
    gr_cell = cellfun(@(x)squeeze(mean(x, [1 2])), gr_cell, 'UniformOutput', false);
    
    %% LOOP THROUGH CORRELATIONS FOR EVERY PATCH

    % get the synchrony time series
    files = dir(fullfile(['F:\' mouse '\' mouse ' GLM predictors\reg2_' mouse '_' condition(1:3) run], '*.csv')); % get files
    names = cell(size(files,1), 1);
    for i = 1:length(names)
        names{i} = [files(i).folder '\' files(i).name];
    end
    T = readtable(names{1});
    idx_high = T.mcor_all > mean(T.mcor_all);
    idx_low = T.mcor_all < mean(T.mcor_all); % get the low and high synchrony time periods in the form of idx
    
    % set empty correlation and distance variables to fill
    corr_low = [];
    corr_high = [];
    % corr = [];
    distance = [];

    % loop through patch
    for i = 5:25 % for partition = 4, use 20:40; for partition = 8, use 5:25
        i
        for j = 5:25 % for partition = 4, use 20:40; for partition = 8, use 5:25

            % process the target patch
            green = gr_cell{i,j};
            red = rr_cell{i,j};
            [~,~, target_px,~,~] = regress(reshape((green-mean(green)), [length(green),1]), [ ones(length(red), 1) reshape((red-mean(red)), [length(red), 1]) ]);
            target_px = smooth(zscore(detrend(target_px)), smooth_kernel/length(target_px)); 
            target_px_low = target_px(idx_low); target_px_high = target_px(idx_high); % split into high and low synchrony
            
            % set empty variables to fill
            corr_1px_low = [];
            corr_1px_high = [];
            % corr_1px = [];
            distance_1px = [];

            %correlate the target pixels to all other pixels
            for k = 1:256/partition
                for l = 1:256/partition

                    %process the comparison pixel
                    green = gr_cell{k,l};
                    red = rr_cell{i,j};
                    [~,~, comparison_px,~,~] = regress(reshape((green-mean(green)), [length(green),1]), [ ones(length(red), 1) reshape((red-mean(red)), [length(red), 1]) ]);
                    comparison_px = smooth(zscore(detrend(comparison_px)), smooth_kernel/length(comparison_px)); 
                    comparison_px_low = comparison_px(idx_low); comparison_px_high = comparison_px(idx_high); % split into high and low synchrony
                    
                    %get the correlation and the distance
                    c_low = corrcoef(target_px_low,comparison_px_low);
                    corr_1px_low = [corr_1px_low c_low(1,2)];
                    c_high = corrcoef(target_px_high,comparison_px_high);
                    corr_1px_high = [corr_1px_high c_high(1,2)];
                    % c = corrcoef(target_px,comparison_px);
                    % corr_1px = [corr_1px c(1,2)];
                    
                    %get the distance
                    distance_1px = [distance_1px pdist2([i,j], [k,l], 'euclidean')]; % distance is the same for high and low
                end
            end
            %end of correlating target patch to all other patches

            % add to corr and distance
            corr_low = [corr_low; corr_1px_low];
            corr_high = [corr_high; corr_1px_high];
            % corr = [corr; corr_1px];
            distance = [distance; distance_1px];
        end 
    end
    %end of looping through target pixels

    %save mean
    y_low = reshape(corr_low, [], 1);
    y_high = reshape(corr_high, [], 1);
    % y = reshape(corr, [], 1);
    x = reshape(distance, [], 1);
    distance2 = round(x);
    x_new = 1:max(distance2);
    x_new = x_new';
    y_mean_low = zeros(max(distance2),1);
    y_mean_high = zeros(max(distance2),1);
    % y_mean = zeros(max(distance2),1);
    for i = 1:max(distance2)
        y_mean_low(i) = mean(y_low(find(distance2 == i)));
        y_mean_high(i) = mean(y_high(find(distance2 == i)));
        % y_mean(i) = mean(y(find(distance2 == i)));
    end

    % save mean
    all_corr_low = [all_corr_low y_mean_low];
    all_corr_high = [all_corr_high y_mean_high];
    % all_corr = [all_corr y_mean];
    all_dist = [all_dist x_new];
    toc
end
save("saline_spatialcorr_largegrid_lowsynchrony.mat","all_corr_low","all_dist", "partition");
save("saline_spatialcorr_largegrid_highsynchrony.mat","all_corr_high","all_dist", "partition");

% %% plotting 2
% figure;
% x1 = mean(all_dist, 2)*4*0.886;
% plot(x1, mean(all_corr, 2), 'k', 'LineWidth', 2);
% bci = bootci(200, @mean, all_corr');
% hold on;
% plot(x1, bci(2, :), 'k--', 'LineWidth', 1);
% hold on;
% plot(x1, bci(1, :), 'k--', 'LineWidth', 1);
% xlim([4,100]);
% ylim([0, .06]);
% title('Spatial Autocorrelation');
% xlabel('Distance (microns)');
% ylabel('Autocorrelation');
% box off;
% axis square;
% %ax = gca; ax.LineWidth = 3; ax.FontSize = 14;
% %% plot
% y = reshape(corr, [], 1);
% x = reshape(distance, [], 1);
% scatter(x,y);
% %%
% distance2 = round(x);
% y_mean = zeros(max(distance2),1);
% for i = 1:max(distance2)
%     y_mean(i) = mean(y(find(distance2 == i)));
% end
% %%
% scatter(distance2*4*1.47, y, '*r', 'MarkerFaceAlpha', 0.3); hold on;
% x_new = 1:max(distance2);
% plot(x_new*4*1.47, y_mean, 'ok-');
% ylim([-0.5, 0.5]);
% yline(0);
% box off;
% ax = gca;
% ax.LineWidth = 3;
% xlabel('Distance (microns)');
% ylabel('correlation');
% title('Correlation vs. Distance');
% axis square;
% xlim([0, 350]);
% ax.FontSize = 14;
% %%
% plot(x_new*4*1.47, y_mean, 'ok-');
% ylim([-0.02, 0.04]);
% yline(0);
% box off;
% % ax = gca;
% % ax.LineWidth = 3;
% xlabel('Distance (microns)');
% ylabel('correlation');
% title('Correlation vs. Distance');
% axis square;
% xlim([0, 350]);
% % ax.FontSize = 14;


