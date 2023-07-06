% CLUSTER_REFINEMENT_II
% Loads the refined cluster track IDs produced by cluster_refinement_I and
% the source datafiles to calculate recruitment intervals (time between
% addition of each new track to the cluster). We produce a plot of the current
% recruitment intervals and report the % of intervals below that value.
% Empirically, we refine the assignment based on the observation that
% ~90-95% of cluster recruitment % intervals are < 4 seconds, and the
% distribution strongly tailed off.
%
% Input:
%   ...
%   RIcutoff: numeric value for temporal cluster segmentation [seconds]
%   ...
%
% Saves:
%
% Part of the cluster_segmentation.m pipeline.
% 2018 May 9 / Mike Pablo

function cluster_refinement_II(refIclusfile, RIcutoff, exposure, tracks, cell_name)

tagclusID = load(refIclusfile);

figSubDir = "RIs_pre_refinementII/";
parentFigDir = 'figures/' + figSubDir;
warning('off');
mkdir(parentFigDir);
warning('on');

refined_cluster_track_IDs = refine_cluster_results(tagclusID.cluster_track_IDs, ...
                                                   exposure, tracks, RIcutoff,...
                                                   parentFigDir,cell_name);

save("data/" + cell_name + "/" + cell_name + '_refinementII_clusTrackIDs.mat','refined_cluster_track_IDs');

end

function cluster_track_IDs = refine_cluster_results(clusTrackIDs,exposure,curr_cell,RIcutoff,parentFigDir,figName)
% Calculates recruitment interval data and reassigns the clusters
% If recruitment interval is greater th an RIcutoff, then consider it a new cluster

expdata.smLinked = curr_cell;
initial_nclus = numel(clusTrackIDs);
cluster_track_IDs = {};

original_RI_distribution = cell(numel(clusTrackIDs),1);

%% Reassign tracks to clusters.
for j=1:initial_nclus
    trackIDs = clusTrackIDs{j};
    track_starts = cellfun(@(x) x(1,end-1),expdata.smLinked(trackIDs));

    % Determine track starts and corresponding IDs for the cluster.
    trackStartandID = [track_starts(:), trackIDs(:)]; % trackIDs is 1xN vector, starts is Nx1
    trackStartandID=sortrows(trackStartandID,1); % Sort by start point
    all_recruitment_intervals = (trackStartandID(2:end,1)-trackStartandID(1:end-1,1))*(exposure*0.001); % in seconds(?)
    original_RI_distribution{j} = all_recruitment_intervals;
    % Does this cluster need temporal separation?
    if any(all_recruitment_intervals>RIcutoff)
        temp_RIs = [all_recruitment_intervals;nan];
        subcluster_assignments = nan(numel(trackIDs),1);
        curr_subcluster = 1;

        % Assign all the tracks into subclusters.
        for k=1:(numel(temp_RIs)-1) % last entry in tempRIs is edge case
           if temp_RIs(k) < RIcutoff
               subcluster_assignments(k) = curr_subcluster;
           else
               subcluster_assignments(k) = curr_subcluster;
               curr_subcluster = curr_subcluster+1;
           end
        end
        subcluster_assignments(end) = curr_subcluster;
        nsubclusters = curr_subcluster;

        for k=1:nsubclusters
            subcluster_idxs  = subcluster_assignments==k;
            cluster_track_IDs = [cluster_track_IDs,trackStartandID(subcluster_idxs,2)'];
        end
    else % the whole cluster was 'normal', we don't even need to use re-sorted information.
        cluster_track_IDs = [cluster_track_IDs,trackIDs];
    end
end
cluster_track_IDs = cluster_track_IDs';

original_RI_distribution = cell2mat(original_RI_distribution);

% Plot the original RI distribution
f=figure;
try
    histogram(original_RI_distribution,'binedges',0:0.02:max(original_RI_distribution))
catch
    warning(['Could not graph RI distribution for ' figName])
end
xlabel('recruitment interval (s)')
ylabel('count')
percent_below_threshold = (sum(original_RI_distribution<RIcutoff)/numel(original_RI_distribution))*100;
title(sprintf('%.2f %% of raw RIs < %.2f',percent_below_threshold,RIcutoff));
set(gca,'fontsize',14,'fontname','arial','tickdir','out')

savename = parentFigDir + figName;
savefig(f,savename + '.fig');
print(f,'-dtiff',savename + '.tif','-r300');

set(gca,'xlim',[0 RIcutoff]);
savename = parentFigDir + figName;
savefig(f,savename + '_zoom.fig');
print(f,'-dtiff',savename + '_zoom.tif','-r300');

% close(f);
fprintf('For %s, %.2f %% of input RIs < %.2f\n',figName,percent_below_threshold,RIcutoff);
end
