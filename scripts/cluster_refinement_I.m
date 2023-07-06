% CLUSTER_REFINEMENT_I
% Searches through tracks not assigned to clusters in the initial Voronoi
% segmentation step, then assigns those tracks to clusters if they diffuse
% into a cluster's Voronoi polygon.
%
% If a track happens to be found within more than one Voronoi polygon,
% make a warning and but assign it to the cluster it spends more time in.
% Random assignment in case of a tie.
%
% Transforms the Voronoi polygons into track IDs identifying the cluster.
%
% Input:
%
% Saves:
%
% Part of the cluster_segmentation.m pipeline.
% 2018 May 8 / Mike Pablo
function cluster_refinement_I(clusfile, tracks, cell_name,centroids,shape,pixel_size)
% also returns the total # of tracks [that were unassigned to clusters] and that visited multiple voronoi polygons,
% and the # of tracks that were unassigned to clusters.

disp("Cluster Refinement I: Identifying tracks that diffuse in to the clusters for " + cell_name)

vorSeg = load(clusfile).segmented_clusters;

% Find centroids outside the ROI
toReject = ~isinterior(shape,centroids(:,1),centroids(:,2));

% Find centroids used for cluster detection (they in cluster by definition)
n_clus = numel(vorSeg);
cluster_track_IDs = cell(n_clus,1);
for j=1:numel(tracks)
    if toReject(j) % If rejected already, don't bother checking
        continue;
    end
    for k=1:n_clus % Find which cluster the centroid belongs to
        if isinterior(vorSeg{k},centroids(j,:))
            cluster_track_IDs{k} = [cluster_track_IDs{k},j];
            break;
        end
    end
end

% Tracks that aren't in clusters, and are within the cell ROI are examined further
diffuseToCluster_candidateIDs = false(numel(tracks),1);
diffuseToCluster_candidateIDs(~toReject) = true; % want tracks inside cell ROI
diffuseToCluster_candidateIDs(cell2mat(cluster_track_IDs')) = false; % don't want tracks already in clusters

diffuseToCluster_candidateIDs = find(diffuseToCluster_candidateIDs); % want numerical indexing

nMultVisTrackFromUnassigned = 0;
nUnassigned = numel(diffuseToCluster_candidateIDs);
for j=1:nUnassigned
    coords = tracks{diffuseToCluster_candidateIDs(j)}(:,1:2)*pixel_size;

    cluster_assigned = false;

    % Each position is tested to see whether it visited a cluster
    cluster_visit_matrix = false(n_clus,size(coords,1));
    for k=1:n_clus
        cluster_visit_matrix(k,:) = isinterior(vorSeg{k},coords);
        if sum(cluster_visit_matrix(k,:)) > (size(coords,1)/2) % more than half the positions are inside a cluster
            cluster_track_IDs{k} = [cluster_track_IDs{k},diffuseToCluster_candidateIDs(j)];
            cluster_assigned = true;
        end
    end

    clusters_visited = logical(sum(cluster_visit_matrix,2));
    if sum(clusters_visited)>1
        nMultVisTrackFromUnassigned = nMultVisTrackFromUnassigned+1;
    end

    % If the track wasn't immediately assigned, need to determine the maximum entry along the visit matrix
    if ~cluster_assigned
        if sum(cluster_visit_matrix(:)) == 0
            continue; % no cluster assignment for this track
        else
            cluster_tallys = sum(cluster_visit_matrix,2);
            [~,favorite_cluster] = max(cluster_tallys);
            cluster_track_IDs{favorite_cluster} = [cluster_track_IDs{favorite_cluster},diffuseToCluster_candidateIDs(j)];
        end
    end
end

save("data/" + cell_name + "/" + cell_name + '_refinementI_clusTrackIDs.mat', ...
     'cluster_track_IDs','nMultVisTrackFromUnassigned','nUnassigned');

end
