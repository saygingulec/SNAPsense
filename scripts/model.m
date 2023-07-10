function model(min_frame,nlayers,layer_width,exclude_clusters,cell_name,exposure,shape,all_tracks,all_centroids,all_activities,all_clusters,all_edge_dists,clusfile)
%%
disp("Modeling the dynamics of " + cell_name)

fig_save_dir = "figures/model/";
warning off
mkdir(fig_save_dir)
warning on


% Variables
exposure = exposure * 0.001; % exposure in seconds
min_lt = min_frame*exposure;
if ~exist('exclude_clusters','var')
    exclude_clusters = 0;
end
segmented_clusters = load(clusfile).segmented_clusters;

rates_mat = [];  % in case you want to summarise all conditions

% conditions: [clustering distance-to-edge]
cluster_conditions = ["", " In Cluster", " Not In Cluster"];
layer_conditions = ["", " Layer 1", " Layer 2", " Layer 3", " Layer 4", " Layer Int"];
if ~exclude_clusters
    conds = combvec([1 0 -1], [0:nlayers Inf])';
else
    conds = combvec([1], [0:nlayers Inf])';
    all_clusters = ones(size(tracks));
end

% is_long    = cellfun(@height, all_tracks) >= min_frame;
is_long_lt = cellfun(@(x) x(end,10) - x(1,10) + 1, all_tracks).*exposure >= min_lt;

for conditon = 1:height(conds)
%% Filter tracks

% Track selection
if conds(conditon,2) == Inf  % Internal tracks
    selected = is_long_lt & logical(all_clusters) + conds(conditon,1) & all_edge_dists > nlayers*layer_width;
elseif conds(conditon,2)  % Layer tracks
    selected = is_long_lt & logical(all_clusters) + conds(conditon,1) & all_edge_dists > (conds(conditon,2)-1)*layer_width & all_edge_dists <= (conds(conditon,2))*layer_width;
else  % All tracks
    selected = is_long_lt & logical(all_clusters) + conds(conditon,1);
end

% Label
if conds(conditon,2) == Inf
    label = cluster_conditions(2-conds(conditon,1)) + layer_conditions(end);
else
    label = cluster_conditions(2-conds(conditon,1)) + layer_conditions(conds(conditon,2)+1);
end

disp(label)

tracks     = all_tracks(selected);
centroids  = all_centroids(selected,:);
activities = all_activities(selected);
clusters   = all_clusters(selected);
edge_dists = all_edge_dists(selected);

% Names
suptitle_suffix = cell_name + label;
save_suffix = replace(lower(suptitle_suffix), ' ', '_');

% Check selected tracks
figure
hold on
scatter(centroids(:,1),centroids(:,2),'.')
plot(shape,FaceAlpha=0)
cellfun(@(x) plot(x, 'FaceAlpha', 0), segmented_clusters)
title(suptitle_suffix + newline + "Selected Tracks", Interpreter="none")
hold off
savefig('figures/model/' + save_suffix + '_selected_tracks.fig');
print('figures/model/' + save_suffix + '_selected_tracks.tiff', '-dtiff','-r300');


%% Calculate rates and save plots

[koff,kdi,kda,koi,ka,koa,ki,pona,poni,psucs,psuc,inact_pop,act_pop] = calculate_rates(exposure, suptitle_suffix, save_suffix, fig_save_dir, tracks, activities);
rates_mat = [rates_mat;koff,kdi,kda,koi,ka,koa,ki,pona,poni,psucs,psuc,inact_pop,act_pop];

save("data/" + cell_name + "/" + save_suffix + "_rates.mat", ...
     "koff", "pona", "poni", "koa", "ki", "koi", "ka", "psuc", "psucs", ...
     "inact_pop", "act_pop", "tracks", "shape", "save_suffix", ...
     "cell_name", "activities", "all_clusters", "edge_dists", "min_frame", "label") %"B", 


end
