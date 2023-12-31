function draw_centroids(centroids, cell_ps, cell_name, fig_save_dir)

warning('off','MATLAB:MKDIR:DirectoryExists')

if exist('fig_save_dir','var')
mkdir(fig_save_dir)
f = figure;
end

hold on
plot(cell_ps)
scatter(centroids(:,1), centroids(:,2), 1, '.')
xlim([min(cell_ps.Vertices(:,1)) max(cell_ps.Vertices(:,1))])
ylim([min(cell_ps.Vertices(:,2)) max(cell_ps.Vertices(:,2))])
title(cell_name, 'interpreter', 'none')
daspect([1 1 1])
hold off

if exist('fig_save_dir','var')
savefig(f,fig_save_dir + "/" + cell_name + "_roi_test.fig");
print(f,'-dtiff',fig_save_dir + "/" + cell_name + "_roi_test.tif",'-r300');
end

end