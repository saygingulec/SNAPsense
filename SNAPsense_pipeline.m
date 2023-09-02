function SNAPsense_pipeline(cell_name)

% The cell_name follows the format: Date_CellLine_Strain_Replicate_CellID_Exposure
% cell_name leads to cell_name_cell_info.mat
% cell_name_cell_info.mat contains the struct 'cell' with at least the following attributes:
%     cell.Path: Path to unprocessed data.
%     cell.Date: Date of the experiment.
%     cell.Cell_Line
%     cell.Strain
%     cell.Replicate
%     cell.Cell_ID: Name given by the experimenter
%     cell.Exposure: Float in ms units. Default is 20 ms.
%     cell.Cell_Name: Name used in the pipeline
%     cell.Raw_Data: Matrix for each localization with columns:
%           [yy (pixels), xx (pixels), i_ch1, bkg_ch1, i_corrected_ch1, i_ch2, bkg_ch2, i_corrected_ch2, i_corrected_ch2/i_corrected_ch1, t, trackID];
%     cell.Tracks: Cell array where each cell contains the same columns for one track.
%     cell.Active: Cell array where each cell contains the activity (0:inactive, 1:active) for one track.
%     cell.ROI: ImageJ ROI object
%     cell.Shape: ImageJ ROI as a polyshape object in um
%     cell.Edge_Dists: Distance to cell edge for each localization

% Parameters
pixel_size = 0.1067;  % pixel size in um
min_frame = 4;  % Minimum lifetime in frames for tracks to be used in 
nlayers = 4;  % Number of layers by the edge
layer_width = 1;  % Layer width in um

% Cluster analysis parameters
RIcutoff = 2.5;  % Recruitment interval cutoff for the temporal separation of clusters
minN = 10;  % Minimum number of particles/tracks(?) for a cluster to count
pad_size = 0.1; % Ignore voronoi polygons x um near the edge.
nsims = 5; % Number of random distributions to simulate.

% Paths
addpath(genpath('scripts/'))
fig_save_dir = "results/figures/";
rate_save_dir = "results/rates/";
warning off
mkdir(fig_save_dir)
mkdir(rate_save_dir)
warning on
cell_folder = "data/" + cell_name + "/";
cell_info_mat = cell_folder + cell_name + "_cell_info.mat";
exp_distrs = cell_folder + cell_name + '_exp_den_distr.mat';
sim_distrs = cell_folder + cell_name + '_sim_den_distr.mat';
threshold_data = cell_folder + cell_name + '_voronoi_thresh.mat';
clusfile = cell_folder + cell_name + '_voronoi_clus.mat';
refIclusfile = cell_folder + cell_name + '_refinementI_clusTrackIDs.mat';
refIIclusfile = cell_folder + cell_name + '_refinementII_clusTrackIDs.mat';

% Load cell info
load(cell_info_mat,'cell_info');


%% Check data
draw_centroids(cell_info.Centroids,cell_info.Shape,cell_info.Cell_Name,fig_save_dir+"/draw_centroids/")
draw_activity_histogram(cell_info.Raw_Data, cell_info.Active, cell_info.Cell_Name,fig_save_dir+"/activity_histogram/")

%% Cluster analysis

% Simulate 50 cells with the same ROI and number of centroids as the cell
% in the experiment and outputs all_cells_sim_den_distr_50x.mat with 
% density_distrs and pad_size.
simulate_uniform_density_distributions(cell_info.Centroids,cell_info.Cell_Name,cell_info.Shape,nsims,pad_size,sim_distrs);

% Calculate the experimentally observed density distributions. Outputs
% all_cells_exp_den_distr.mat with density_distrs and pad_size.
generate_density_distributions(cell_info.Centroids,cell_info.Cell_Name,cell_info.Shape,pad_size,exp_distrs);

% Compare the densities of experimental voronoi segments to the simulated
% ones and decides on a cutoff point for the cluster density. Saves
% 'cutoffs', 'bin_cutoffs', 'binwidth' to 'voronoi_thresh.mat'.
compare_voronoicell_size_distributions(exp_distrs, sim_distrs, cell_info.Cell_Name, threshold_data, fig_save_dir+"/vor_seg_clus_comp/"); 

% Pick the voronoi segments that dense enough to be clusters and merge
% adjacent ones. Saves 'segmented_clusters', 'pad_size', 'threshold_data'
% to 'voronoi_clus.mat'.
generate_voronoiSegmented_clusters(threshold_data, exp_distrs, cell_info.Centroids,cell_info.Cell_Name, cell_info.Shape, pad_size, clusfile, fig_save_dir+"/raw_clusters/");

% Searches through tracks not assigned to clusters in the initial Voronoi
% segmentation step, then assigns those tracks to clusters if they diffuse
% into a cluster's Voronoi polygon. Saves 'cluster_track_IDs',
% 'nMultVisTrackFromUnassigned' (number of tracks that visited multiple
% clusters) and 'nUnassigned' (number of unassigned tracks) in
% 'refinementI_clusTrackIDs'.
cluster_refinement_I(clusfile, cell_info.Tracks, cell_info.Cell_Name, cell_info.Centroids,cell_info.Shape,pixel_size,refIclusfile);

% Temporally separate clusters
cluster_refinement_II(refIclusfile, RIcutoff, cell_info.Exposure, cell_info.Tracks, cell_info.Cell_Name,refIIclusfile,fig_save_dir+"/RIs_pre_refinementII/"); 

% Add clusterIDs and make RefII_minN cluster figures
[cell_info.RefI_ClusID, cell_info.RefI_ClusID_minN, cell_info.RefII_ClusID, cell_info.RefII_ClusID_minN] ...
    = examine_clusters(minN, clusfile, refIclusfile, refIIclusfile, cell_info.Tracks, cell_info.Cell_Name, cell_info.Shape, cell_info.Centroids, fig_save_dir+"/RefII_minN_clusters/");

% Save results
cell_info.minN = minN;
save(cell_info_mat, 'cell_info')


%% Calculate rates
model(min_frame, nlayers, layer_width, 0, cell_info.Cell_Name, cell_info.Exposure, cell_info.Shape, cell_info.Tracks, ...
    cell_info.Centroids, cell_info.Active, cell_info.RefI_ClusID, cell_info.Edge_Dists,clusfile,fig_save_dir+"/model/",rate_save_dir)

end