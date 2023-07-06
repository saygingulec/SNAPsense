function draw_activity_histogram(raw_data,activity,cell_name)

warning off
mkdir figures\activity_histogram
warning on

activities = cell2mat(activity);

figure
hold on
histogram2(raw_data(:,5), raw_data(:,8), 100, ...
    'DisplayStyle','tile', 'ShowEmptyBins','off', 'EdgeAlpha',0)
daspect([1 1 1])
active_points = raw_data(activities,[5 8]);
active_bound = convhull(active_points);
plot(active_points(active_bound,1), active_points(active_bound,2),'g',LineWidth=1)
legend(["" "Active"])
hold off

savefig(f,"figures/activity_histogram/" + cell_name + "_activity_histogram.fig");
print(f,'-dtiff',"figures/activity_histogram/" + cell_name + "_activity_histogram.tif",'-r300');

end
