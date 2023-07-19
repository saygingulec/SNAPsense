% COMPARE_VORONOICELL_SIZE_DISTRIBUTIONS_SIM50X
% Determines threshold densities for cluster segmentation.
%
% Input:
%   threshSuffix: string for file + figure naming, e.g.
%       'voronoi_thresholds'
%   parentFigDir: string for parent figure directory, e.g.
%       'figures_20180508/'
%
% Saves a .mat file to the working directory:
%   [dateprefix threshSuffix '.mat']
%
% This sets thresholds for the segmentation step.
%
% Part of the cluster_segmentation.m pipeline.

function compare_voronoicell_size_distributions(exp_distrs, sim_distrs, cell_name, threshold_data)
%%
warning('off');
mkdir('figures');
mkdir('figures/vor_seg_clus_comp')
warning('on');

exp = load(exp_distrs);
sim = load(sim_distrs);

% Get the localization areas in nm^2
exp_areas = 1./exp.density_distrs;

% the simulations are a nested cell array. There are 50 "cell array" elements
% per biological cell (ncells)
sim_areas = cellfun(@(densities) 1./densities, sim.density_distrs,'uniformoutput',false);
tag_color = [0,0.4470,0.7410];
sim_color = [0,0,0];

% Get cutoff
% It's useful to check the single cell comparisons to make sure no problems
% occurred during threshold determination.
disp('Calculating cutoff for ' + cell_name)

binwidth = 100; % nm^2; corresponds to 10 nm pixel
cutoffs = get_cutoff(exp_areas,sim_areas,binwidth);

% Save cutoff
cutoffs = 1./cutoffs; % convert area to density
bin_cutoffs = nan(size(cutoffs));
save(threshold_data,'cutoffs','bin_cutoffs','binwidth');

% Plot and save
plot_single_cell_smoothed_intersected(exp_areas, sim_areas, tag_color, ...
    sim_color, 'Exp', 'Sim', binwidth, [0 2000000], "figures/vor_seg_clus_comp/" + cell_name);
end


function cutoff = get_cutoff(expdata,simdata50x,binwidth)
%%
% Find the intesection of the probability distributions between the
% simulated mean and the experimental data. Apply smoothing prior to
% intersection checking; empirically chose 5pt Gaussian smooth.
% Check the single cell traces to verify that the cutoffs are consistent with
% the distributions' appearances.
[exp_distr,exp_edges]=histcounts(expdata,'binwidth',binwidth);
nsimreps = height(simdata50x);
sim_distr=cell(nsimreps,1);
sim_edges=cell(nsimreps,1);

% Figure out the appropriate bin limits across the whole sim data
maxbin = min(cellfun(@max,simdata50x));
minbin = min(cellfun(@min,simdata50x));

for i=1:nsimreps
    [sim_distr{i},sim_edges]=histcounts(simdata50x{i},'binwidth',binwidth,'binlimits',[minbin maxbin]);
end
mean_sim_distr = mean(cell2mat(sim_distr));
exp_cents = get_bin_centers(exp_edges);
sim_cents = get_bin_centers(sim_edges);
exp_distr = smoothdata(exp_distr,'gaussian',5);
mean_sim_distr = smoothdata(mean_sim_distr,'gaussian',5);
[x_ints,~,~,~] = intersections(exp_cents,exp_distr,sim_cents,mean_sim_distr);
cutoff = x_ints(1);
end


function plot_single_cell_smoothed_intersected(expdata,simdata50x,expcolor, ...
    simcolor,explabel,simlabel,binwidth,xplotwindow,savename)
%%
f= figure('visible','on'); hold on
[exp_distr,exp_edges]=histcounts(expdata,'binwidth',binwidth);
nsimreps = height(simdata50x);
sim_distr=cell(nsimreps,1);
sim_edges=cell(nsimreps,1);

% Figure out the appropriate bin limits across the whole sim data
maxbin = min(cellfun(@max,simdata50x));
minbin = min(cellfun(@min,simdata50x));

for i=1:nsimreps
    [sim_distr{i},sim_edges]=histcounts(simdata50x{i},'binwidth',binwidth,'binlimits',[minbin maxbin]);
end
mean_sim_distr = mean(cell2mat(sim_distr));

exp_cents = get_bin_centers(exp_edges);
sim_cents = get_bin_centers(sim_edges);
exp_distr = smoothdata(exp_distr,'gaussian',5);
sim_distr = smoothdata(mean_sim_distr,'gaussian',5);
[x_ints,y_ints,~,~] = intersections(exp_cents,exp_distr,sim_cents,sim_distr);
hExp=plot(exp_cents,exp_distr,'linewidth',4,'color',expcolor);
hSim=plot(sim_cents,sim_distr,'linewidth',4,'color',simcolor);
plot(x_ints(1),y_ints(1),'mo','markerfacecolor','m','markersize',8)
legend([hExp,hSim],explabel,simlabel,'location','northeast');
xlabel('Area of Voronoi cell (nm^2)')
ylabel('Counts');
set(gca,'fontsize',16,'fontname','arial','xlim',xplotwindow);

savefig(f,savename + '.fig');
print(f,'-dtiff',savename + '.tif','-r300');

% Zoomed figure
xlim([0 max(sim_cents)])
savefig(f,savename + '_zoomed.fig');
print(f,'-dtiff',savename + '_zoomed.tif','-r300');

close(f)
end

function centers = get_bin_centers(binedges)
centers = mean([binedges(1:end-1);binedges(2:end)]);
end
