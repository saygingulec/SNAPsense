function rate_figure(pona,poni,koa,ki,koi,ka,active_marker_size,inactive_marker_size)

s = ["Start" "Start" "Active" "Active" "Inactive" "Inactive"];
t = ["Active" "Inactive" "End" "Inactive" "End" "Active"];
weights = [10*pona 10*poni koa ki koi ka];
G = digraph(s,t,weights);
LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
plot(G,'LineWidth',LWidths, 'ArrowSize', 15, ...
     'EdgeLabel',[compose("%.4g%%",G.Edges.Weight(1:2)*10); string(round(G.Edges.Weight(3:end), 2))], ...
     NodeColor=[0 0 0; 0 1 0; 0 0 1; 1 0 0], EdgeFontSize=10, NodeFontSize=10, ...
     MarkerSize=[5 20*active_marker_size 20*inactive_marker_size 5])

end