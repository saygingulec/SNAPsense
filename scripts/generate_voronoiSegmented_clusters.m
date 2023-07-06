% GENERATE_VORONOISEGMENTED_CLUSTERS
% Generates Voronoi segmentations from experimental data,
% then uses the 50x simulation-based threshold to delete irrelevant
% polygons. Merges remaining touching polygons to obtain clusters (polyshapes)
%
% Input:
%
% Saves:
%
% Part of the cluster_segmentation.m pipeline.


function generate_voronoiSegmented_clusters(threshold_data, exp_distrs, centroids, cell_name, shape, pad_size)

warning off
mkdir('figures/raw_clusters')
warning on

disp("Generating clusters of " + cell_name)

threshold_data = load(threshold_data);
thresh = threshold_data.cutoffs;
exp_distrs = load(exp_distrs);
polys = exp_distrs.polys;
poly_dens = exp_distrs.density_distrs;

% Find centroids in ROI
centroids = centroids( isinterior(shape,centroids(:,1),centroids(:,2)) , : );

% Find polys denser than simulated threshold
dense_polys = polys(poly_dens >= thresh);
npolys = numel(dense_polys);

% Use graph theoretic method to get connected polygons.
adjacency_mat = false(npolys,npolys);
for j=1:npolys
    for k=j:npolys
        v1 = dense_polys{j}.Vertices;
        v2 = dense_polys{k}.Vertices;
        if share_vertex(v1,v2)
            adjacency_mat(j,k) = true;
            adjacency_mat(k,j) = true;
        end
    end
end
% Get the IDs of base polygons composing unique clusters
polygraph = graph(adjacency_mat);
clusters=conncomp(polygraph); % Each polygon is assigned to a cluster ID
cluster_list = unique(clusters);
nclus = numel(cluster_list);
segmented_clusters = cell(nclus,1); % we'll remove empty entries after
for j=1:nclus
    curr_clus = cluster_list(j);
    curr_poly_in_clus = find(curr_clus==clusters);

    % Want to have at least three sufficiently dense polygons to start a cluster
    if numel(curr_poly_in_clus) < 3
        continue;
    end
    % Merge together polygons that are attached
    curr_seg_clus = union(dense_polys{curr_poly_in_clus(1)},...
                          dense_polys{curr_poly_in_clus(2)});
    for k=3:numel(curr_poly_in_clus)
        curr_seg_clus = union(curr_seg_clus,...
                              dense_polys{curr_poly_in_clus(k)});
    end
    segmented_clusters{j} = curr_seg_clus;
end
segmented_clusters = segmented_clusters(cellfun(@(x) ~isempty(x),segmented_clusters));

% Save output
save("data/" + cell_name + "/" + cell_name + '_voronoi_clus.mat','segmented_clusters','pad_size','threshold_data');


% Plot
figure(Visible="on")
hold on
plot(shape, 'FaceColor', 'none')
scatter(centroids(:,1),centroids(:,2),10,'.')
cellfun(@(x) plot(x), segmented_clusters)
title([cell_name ' Raw Clusters'], Interpreter="none")
daspect([1 1 1])
hold off
savefig("figures/raw_clusters/" + cell_name + "_raw_clusters.fig");
print("figures/raw_clusters/" + cell_name + "_raw_clusters.tiff", '-dtiff','-r300');
close

end

function dotouch = share_vertex(v1,v2)
nv1 = size(v1,1);
dotouch=false;
for i=1:nv1
    if any(v1(i,:)==v2)
        dotouch = true;
        break;
    end
end
end
