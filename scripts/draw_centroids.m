function draw_centroids(centroids, cell_ps, cell_name)

warning('off','MATLAB:MKDIR:DirectoryExists')

mkdir figures/draw_centroids

f = figure;
hold on
plot(cell_ps)
scatter(centroids(:,1), centroids(:,2), 1, '.')
xlim([min(cell_ps.Vertices(:,1)) max(cell_ps.Vertices(:,1))])
ylim([min(cell_ps.Vertices(:,2)) max(cell_ps.Vertices(:,2))])
title(cell_name, 'interpreter', 'none')
daspect([1 1 1])
hold off

savefig(f,"figures/draw_centroids/" + cell_name + "_roi_test.fig");
print(f,'-dtiff',"figures/draw_centroids/" + cell_name + "_roi_test.tif",'-r300');

end