% GENERATE_DENSITY_DISTRIBUTIONS
% Saves the densities of Voronoi cells obtained from Voronoi tessellation
% of experimental data.
%
% 'pad_size' : um; for ignoring voronoi polygons near the edge.
%              edge polygon sizes can be artificially small.
%
% Input:
%
% Saves:
%
% This is for thresholding purposes, and must be compared vs. the
% simulated density distributions.

function generate_density_distributions(centroids, cell_name, shape, pad_size)

disp("Calculating density distributions of " + cell_name)
    
% Uses Voronoi-based approach to estimate local track density.
%	density_distr is a vector is the density for each Voronoi polygon,
%       in localizations/nm2
% We reject polygons from analysis if their centroids are within 2*pad_size
%   of the cell edge, since they might have artificially high density.
% Returns 50x simulations of the density distributions.


% Find centroids in ROI
centroids = centroids( isinterior(shape,centroids(:,1),centroids(:,2)) , : );


% Calculate Voronoi cells and local densities, constrained by cell ROI; um units
[v,c,~]=VoronoiLimit(centroids(:,1),centroids(:,2),...
                     'bs_ext',shape.Vertices,...
                     'figure','off');
npolys = size(c,1);
polys = cell(npolys,1);
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
       polys{j,:} = ps;
    end
end
to_keep = ~isnan(xyden(:,1));
xyden = xyden(to_keep,:);
polys = polys(to_keep,:);

xyden(:,3) = xyden(:,3)/(1000^2);  % Density is now in localizations/nm2
density_distrs = xyden(:,3);

% Save output
save(exp_distrs,'density_distrs','pad_size','polys');


% %% Plot (takes too long to plot)
% warning off
% mkdir figures\vor_densities
% warning on
% cmap = vals2cmap(log(poly_dens));
% figure(Visible="off")
% hold on
% plot(shape, 'FaceColor', 'none')
% for poly_idx = 1:height(polys)
%     plot(polys{poly_idx}, 'FaceColor',cmap(poly_idx,:),'FaceAlpha',1)
% %     drawnow
% end
% title(cell_name + " Density (log(N/nm^2))", Interpreter="none")
% daspect([1 1 1])
% hold off
% savefig("figures/vor_densities/" + cell_name + "_raw_clusters.fig");
% print("figures/vor_densities/" + cell_name + "_raw_clusters.tiff", '-dtiff','-r300');
% close

end
