function simulate_uniform_density_distributions(exp_centroids, cell_name, shape, nsims, pad_size)
% SIMULATE_UNIFORM_DENSITY_DISTRIBUTIONS_50X
% Simulates negative controls for the Voronoi-based cluster segmentation.
% Randomly places points within the cell ROI, corresponding to actual observed
% (and usable) trajectories.
%
% 'pad_size' : um; for ignoring voronoi polygons near the edge.
%              edge polygon sizes can be artificially small.
%
% Input:
%   
%
% Saves:
%
%
% This is for thresholding purposes, and must be compared vs. the
% experimental density distributions.
%
% Part of the cluster_segmentation.m pipeline.

rng('default'); % for reproducibility
disp("Simulating density distributions for " + cell_name)


%% Simulate

% Uses Voronoi-based approach to estimate local track density.
%	density_distrs is a vector is the density for each Voronoi polygon,
%       in localizations/nm2
% We reject polygons from analysis if their centroids are within 2*pad_size
%   of the cell edge, since they might have artificially high density.
% Returns 50x simulations of the density distributions.

% Find centroids in ROI
exp_centroids = isinterior(shape,exp_centroids(:,1),exp_centroids(:,2));

centroids = nan(height(exp_centroids),2); % x,y
window_min = min(shape.Vertices);
window_max = max(shape.Vertices);

density_distrs = cell(nsims,1);

for sim_rep_idx = 1:nsims
    disp("Simulation " + sim_rep_idx + " of " + nsims)

    % Sample random centroids
    for i=1:height(exp_centroids)
        % Rejection sampling: place a point randomly for every track.
        % If outside of the cell boundaries, try again.
        currcentroid = [unifrnd(window_min(1),window_max(1)),unifrnd(window_min(2),window_max(2))];
        while ~isinterior(shape,currcentroid(1),currcentroid(2))
            currcentroid = [unifrnd(window_min(1),window_max(1)),unifrnd(window_min(2),window_max(2))];
        end
        centroids(i,:) = currcentroid;
    end

    % Calculate Voronoi cells and local densities, constrained by cell ROI; um units
    [v,c,~]=VoronoiLimit(centroids(:,1),centroids(:,2),...
                         'bs_ext',shape.Vertices,...
                         'figure','off');
    npolys = size(c,1);
    xyden = nan(npolys,3);
    edgetol = 2*pad_size; % um; Exclude voronoi cell centroids near edge
                          %     they can cause artifacts in the analysis
    for j=1:npolys
        warning('off','all');
        ps = polyshape(v(c{j},:));
        warning('on','all');
        [centx,centy] = centroid(ps);
        cent_to_edge = abs(p_poly_dist(centx,centy,...
                          shape.Vertices(:,1),...
                          shape.Vertices(:,2)));
        if cent_to_edge < edgetol
           continue;
        else
           den = 1/area(ps);  % Density is now in localizations/um2
           xyden(j,:) = [centx,centy,den];
        end
    end
    to_keep = ~isnan(xyden(:,1));
    xyden = xyden(to_keep,:);

    xyden(:,3) = xyden(:,3)/(1000^2);  % Density is now in localizations/nm2
    density_distrs{sim_rep_idx} = xyden(:,3);
end % sim reps

save(sim_distrs,'density_distrs','pad_size', 'cell_name');

