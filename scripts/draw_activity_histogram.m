function draw_activity_histogram(raw_data,activity)

activities = cell2mat(activity{:});

figure
histogram2(raw_data(:,5), raw_data(:,8), 100, ...
    'DisplayStyle','tile', 'ShowEmptyBins','off', 'EdgeAlpha',0)
grid off
