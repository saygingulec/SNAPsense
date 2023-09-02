function rate_figure(pona,poni,koa,ki,koi,ka,active_marker_size,inactive_marker_size,error_type)
%%
s = ["Start" "Start" "Active" "Active" "Inactive" "Inactive"];
t = ["Active" "Inactive" "End" "Inactive" "End" "Active"];
weights = mean([10*pona 10*poni koa ki koi ka],1);
if exist("error_type","var")
    if error_type == "std"
        errors = std([pona poni koa ki koi ka],1);
    elseif error_type == "sem"
        errors = std([pona poni koa ki koi ka],1) ./ sqrt(height([pona poni koa ki koi ka]));
    end
    labels = [compose("%.4g%%",weights(1:2)*10) string(round(weights(3:end), 2))] + repmat("\pm",1,6) + compose("%.2f", errors);
else
    labels = [compose("%.4g%%",weights(1:2)*10) string(round(weights(3:end), 2))];
end
G = digraph(table([s' t'], weights', labels', 'VariableNames', {'EndNodes', 'weights', 'labels'}));
LWidths = 5*G.Edges.weights/max(G.Edges.weights);
plot(G,'LineWidth',LWidths, 'EdgeLabel',G.Edges.labels, 'ArrowSize', 15, ...
     NodeColor=[0 0 0; 0 1 0; 0 0 1; 1 0 0], EdgeFontSize=10, NodeFontSize=10, ...
     MarkerSize=[5 20*mean(active_marker_size) 20*mean(inactive_marker_size) 5])

end