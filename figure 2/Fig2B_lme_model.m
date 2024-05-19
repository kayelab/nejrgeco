%% RUN FIRST TO LOAD DATA
clear;

% paths
addpath('H:\nejrgeco_code_for_publication\');
D = 'H:\*PrL*\*GLM predictors*\*reg2*';
files = dir(fullfile(D, '*.csv'));
names = cell(size(files,1), 1);
for i = 1:length(names)
    names{i} = [files(i).folder '\' files(i).name];
end
all_names = names;
%%
% Get correlations and hierarchical levels as inputs for model

% LOCAL/GLOBAL CORRELATIONS
lg_bg = [];
mouse2 = [];
drug2 = [];
run2 = [];
id2 = [];
corr_type2 = [];

for i = 1:length(all_names)
    T = readtable(all_names{i});
    cor = corr(T.green_donut, T.bulk_green);
    
    % collect additional variables
    name = all_names{i};
    mouse2 = [mouse2 {name(50:55)}];
    drug2 = [drug2 {name(57:59)}];
    run2 = [run2 {str2num(name(60))}];
    corr_type2 = [corr_type2 {'lg_bg'}];
    [C,~] = strsplit(name,{'cell', '.csv'}, 'DelimiterType','RegularExpression');
    if strcmp(name(50:55), 'PrL2.1')
        id2 = [id2 str2num(C{2})];
    elseif strcmp(name(50:55), 'PrL3.2')
        id2 = [id2 str2num(C{2})+19];
    elseif strcmp(name(50:55), 'PrL3.4')
        id2 = [id2 str2num(C{2})+19+28];
    end
    
    lg_bg = [lg_bg cor];
end

%% fit model

tbl2 = table(lg_bg' , id2', drug2', mouse2', 'VariableNames', {'correlations', 'cell_id', 'drug', 'mouse'});
lme = fitlme(tbl2,'correlations~drug+(drug|mouse)+(drug|cell_id)');
lme

%%
% Get correlations and hierarchical levels as inputs for model

% LOCAL/LOCAL CORRELATIONS

mouse = [];
drug = [];
run = [];
id = [];
corr_type = [];

run_names = cellfun(@(x)x(50:60), all_names, 'UniformOutput', false);
unique_run_names = unique(run_names);
lg_lg = [];
for i = 1:length(unique_run_names)
    file_names = all_names(contains(all_names, unique_run_names{i}));
    all_lg = [];
    for k = 1:length(file_names)
        T = readtable(file_names{k});
        all_lg = [all_lg T.green_donut];
    end
    name = file_names{1};
    cor_pairs = corr(all_lg);
    lg_lg = [lg_lg cor_pairs(find(tril(cor_pairs, -1) ~= 0))'];
    mouse = [mouse repmat({name(50:55)}, 1, length(cor_pairs(find(tril(cor_pairs, -1) ~= 0))))];
    drug = [drug repmat({name(57:59)}, 1, length(cor_pairs(find(tril(cor_pairs, -1) ~= 0))))];
    run = [run repmat({str2num(name(60))}, 1, length(cor_pairs(find(tril(cor_pairs, -1) ~= 0))))];
    corr_type = [corr_type repmat({'lg_lg'}, 1, length(cor_pairs(find(tril(cor_pairs, -1) ~= 0))))];
    
    lengths = fliplr(1:length(file_names)-1);
    indices = [];
    for z = 1:length(lengths)
        indices = [indices repmat(z, 1, lengths(z))];
    end
    
    if strcmp(name(50:55), 'PrL2.1')
        id = [id indices];
    elseif strcmp(name(50:55), 'PrL3.2')
        id = [id indices+19];
    elseif strcmp(name(50:55), 'PrL3.4')
        id = [id indices+19+28];
    end
end

%% fit model

tbl = table(lg_lg' , id', drug', mouse', 'VariableNames', {'correlations', 'cell_id', 'drug', 'mouse'});
lme = fitlme(tbl2,'correlations~drug+(drug|mouse)+(drug|cell_id)');
lme
%% EXTRA CODE

% tbl3 = table(lg_bg' , id2', drug2', mouse2', run2', corr_type2', 'VariableNames', {'correlations', 'cell_id', 'drug', 'mouse', 'run', 'corr_type'});
% lme = fitlme(tbl3,'correlations~drug+(drug|mouse)+(drug|cell_id)');
% lme

% cell_correlation ~ cell_id + (cell_id | mouse) + (cell_id | drug)
% tbl = table([lg_lg'; lg_bg'] , [id'; id2'], [drug'; drug2'], [mouse'; mouse2'], [run'; run2'], [corr_type'; corr_type2'], 'VariableNames', {'correlations', 'cell_id', 'drug', 'mouse', 'run', 'corr_type'});
% lme = fitlme(tbl,'correlations~drug+corr_type+(drug|mouse)+(drug|cell_id)+(corr_type|mouse) + (corr_type|cell_id)');
% lme

D = 'F:\PrL*\*GLM predictors*\*reg2_*sal*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('F:\PrL2.1')); addpath(genpath('F:\PrL3.2')); addpath(genpath('F:\PrL3.4')); addpath(genpath('F:\samira_nejrgeco'));
all_names = cell(size(files,1), 1);
for i = 1:length(all_names)
    all_names{i} = [files(i).folder '\' files(i).name];
end

run_names = cellfun(@(x)x(50:60), all_names, 'UniformOutput', false);
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
D = 'F:\PrL*\*GLM predictors*\*reg2_*des*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('F:\PrL2.1')); addpath(genpath('F:\PrL3.2')); addpath(genpath('F:\PrL3.4')); addpath(genpath('F:\samira_nejrgeco'));
all_names = cell(size(files,1), 1);
for i = 1:length(all_names)
    all_names{i} = [files(i).folder '\' files(i).name];
end

D = 'F:\PrL*\*GLM predictors*\*reg2_*sal*';
files = dir(fullfile(D, '*.csv'));
addpath(genpath('F:\PrL2.1')); addpath(genpath('F:\PrL3.2')); addpath(genpath('F:\PrL3.4')); addpath(genpath('F:\samira_nejrgeco'));
all_names1 = cell(size(files,1), 1);
for i = 1:length(all_names1)
    all_names1{i} = [files(i).folder '\' files(i).name];
end

all_names = [all_names; all_names1];

% LOCAL/LOCAL CORRELATIONS

mouse = [];
drug = [];
run = [];
id = [];
corr_type = [];

run_names = cellfun(@(x)x(50:60), all_names, 'UniformOutput', false);
unique_run_names = unique(run_names);
lg_lg = [];
for i = 1:length(unique_run_names)
    file_names = all_names(contains(all_names, unique_run_names{i}));
    all_lg = [];
    for k = 1:length(file_names)
        T = readtable(file_names{k});
        all_lg = [all_lg T.green_donut];
    end
    name = file_names{1};
    cor_pairs = corr(all_lg);
    lg_lg = [lg_lg cor_pairs(find(tril(cor_pairs, -1) ~= 0))'];
    mouse = [mouse repmat({name(50:55)}, 1, length(cor_pairs(find(tril(cor_pairs, -1) ~= 0))))];
    drug = [drug repmat({name(57:59)}, 1, length(cor_pairs(find(tril(cor_pairs, -1) ~= 0))))];
    run = [run repmat({str2num(name(60))}, 1, length(cor_pairs(find(tril(cor_pairs, -1) ~= 0))))];
    corr_type = [corr_type repmat({'lg_lg'}, 1, length(cor_pairs(find(tril(cor_pairs, -1) ~= 0))))];
    
    lengths = fliplr(1:length(file_names)-1);
    indices = [];
    for z = 1:length(lengths)
        indices = [indices repmat(z, 1, lengths(z))];
    end
    
    if strcmp(name(50:55), 'PrL2.1')
        id = [id indices];
    elseif strcmp(name(50:55), 'PrL3.2')
        id = [id indices+19];
    elseif strcmp(name(50:55), 'PrL3.4')
        id = [id indices+19+28];
    end
end
