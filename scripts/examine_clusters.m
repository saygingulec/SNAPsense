function [RefI_ClusID, RefI_ClusID_minN, RefII_ClusID, RefII_ClusID_minN] = examine_clusters(minN, clusfile, refIclusfile, refIIclusfile, tracks, cell_name, shape, centroids, fig_save_dir)

set(groot,'defaultAxesTickLabelInterpreter','none');

segmented_clusters = load(clusfile).segmented_clusters;
refI_trackIDs = load(refIclusfile).cluster_track_IDs;
refII_trackIDs = load(refIIclusfile).refined_cluster_track_IDs;


%% Assign cluster IDs to tracks

% RefI clusters
RefI_ClusID = zeros(height(tracks), 1);
RefI_ClusID_minN = zeros(height(tracks), 1);  % clusters with at least N tracks
for clusterID = 1:height(refI_trackIDs)
    RefI_ClusID(refI_trackIDs{clusterID}) = clusterID;
    if numel(refI_trackIDs{clusterID}) >= minN
        RefI_ClusID_minN(refI_trackIDs{clusterID}) = clusterID;
    end
end

% RefII clusters
RefII_ClusID = zeros(height(tracks), 1);
RefII_ClusID_minN = zeros(height(tracks), 1);  % clusters with at least N tracks
for clusterID = 1:height(refII_trackIDs)
    RefII_ClusID(refII_trackIDs{clusterID}) = clusterID;
    if numel(refII_trackIDs{clusterID}) >= minN
        RefII_ClusID_minN(refII_trackIDs{clusterID}) = clusterID;
    end
end
    

%% Figure RefII_minN clusters

warning off
mkdir('figures/RefII_minN_clusters')
warning on

figure(Visible="on")
hold on
plot(shape, FaceColor='none')
scatter(centroids(:,1), centroids(:,2), 10, [.7 .7 .7], '.')
scatter(centroids(logical(RefII_ClusID_minN),1), centroids(logical(RefII_ClusID_minN),2), 10, '.')
cellfun(@(x) plot(x, 'FaceAlpha', 0), segmented_clusters)
hold off
daspect([1 1 1])
title([cell_name ' Refined Clusters'], Interpreter="none")
savefig(fig_save_dir + "/" + cell_name + '_refined_clusters.fig');
print(fig_save_dir + "/" + cell_name + '_refined_clusters.tiff', '-dtiff','-r300');


end