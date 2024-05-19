clear;
all_events = [];

D = 'H:\*PrL*\*GLM predictors*\*reg2*sal*';
files = dir(fullfile(D, '*.csv'));
folder = {files.folder};
names = natsortfiles({files.name});
names = cellfun(@(x,y)[x '\' y], folder, names, 'UniformOutput', false);


for i = 1:length(names)
    % get synchrony
    T = readtable(names{i});
    mcor_bglg = T.mcor_bglg;

    % get desynchrony event lengths
    x = mcor_bglg<mean(mcor_bglg); x = double(x)'; % desynchrony segments are 1
    idx = find([1,diff(x),1]);
    lengths = idx(2:end)-idx(1:end-1);
    segment_ids = x([true, diff(x) ~= 0]); % collapse the segment
    desynchrony_lengths = lengths(find(segment_ids));
    all_events = [all_events desynchrony_lengths];
end

%%
% set histogram edges
edges = 0:1:60;

% plot
figure;
histogram(all_events/15, edges, 'FaceColor', 'k', 'EdgeColor','none', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2, Normalization="percentage");

% plot formatting
set(gca, 'TickDir', 'out');
ytickformat("percentage")
yticks([0, 5, 10, 15, 20]);
title('Event Lengths');
xlabel('time in seconds');
ylabel('frequency');
axis square;
box off;





