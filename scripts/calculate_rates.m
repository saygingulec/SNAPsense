function [koff,kdi,kda,koi,ka,koa,ki,pona,poni,psucs,psuc,inact_pop,act_pop] = calculate_rates(exposure, suptitle_suffix, save_suffix, save_dir, tracks, activities)

%% Lifetimes and Decay Rates
interp_acts = cellfun( ...
    @(track, act) ...
    interp1(track(:,10), single(act), track(1,10):track(end,10),'nearest'), ...
    tracks, activities, UniformOutput=false);

% Selection
is_all_inact = cellfun(@(x) all(~x),interp_acts);
is_all_act = cellfun(@(x) all(x),interp_acts);
is_switch = ~is_all_act & ~is_all_inact;


% TOTAL
total_lt = cellfun(@width,interp_acts);
min_length = min(total_lt);
total_lt = (total_lt - min_length)*exposure;

koff = exp_fit_figure(total_lt, "Total Lifetimes" + newline + suptitle_suffix);

savefig(save_dir + save_suffix + "_koff" + ".fig");
print('-dtiff',save_dir + save_suffix + "_koff"  + ".tif",'-r300');


% INACTIVE
inact_lt = [cellfun(@width,interp_acts(is_all_inact)) - min_length
            cell2mat(cellfun(@(x) (strfind([x==0,0],[1 0]) - strfind([0,x==0],[0 1]))',interp_acts(is_switch),UniformOutput=false))] ...
           * exposure;
inact_pop = (height(inact_lt));

kdi = exp_fit_figure(inact_lt, "Inactive Lifetimes" + newline + suptitle_suffix);

savefig(save_dir + save_suffix + "_kdi" + ".fig");
print('-dtiff',save_dir + save_suffix + "_kdi" + ".tif",'-r300');


% ACTIVE
act_lt = [cellfun(@width,interp_acts(is_all_act)) - min_length
          cell2mat(cellfun(@(x) (strfind([x==1,0],[1 0]) - strfind([0,x==1],[0 1]))',interp_acts(is_switch),UniformOutput=false))] ...
         * exposure;
act_pop = (height(act_lt));

kda = exp_fit_figure(act_lt, "Active Lifetimes" + newline + suptitle_suffix);

savefig(save_dir + save_suffix  + "_kda" + ".fig");
print('-dtiff',save_dir + save_suffix + "_kda" + ".tif",'-r300');


%% Calculate rates

n_oi = sum(cellfun(@(x) x(end)==0, activities));
n_oa = sum(cellfun(@(x) x(end)==1, activities));

koi = kdi * n_oi / inact_pop;
ka = kdi - koi;

koa = kda * n_oa / act_pop;
ki = kda - koa;

% Plot
f = figure;

s = ["Start" "Start" "Active" "Active" "Inactive" "Inactive"];
t = ["Active" "Inactive" "End" "Inactive" "End" "Active"];
pona = sum(cellfun(@(x) x(1), activities))/height(activities);
poni = sum(~cellfun(@(x) x(1), activities))/height(activities);
weights = [10*pona 10*poni koa ki koi ka];
G = digraph(s,t,weights);
LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
plot(G,'LineWidth',LWidths, 'ArrowSize', 15, ...
     'EdgeLabel',[compose("%.4g%%",G.Edges.Weight(1:2)*10); string(round(G.Edges.Weight(3:end), 2))], ...
     NodeColor=[0 0 0; 0 1 0; 0 0 1; 1 0 0], EdgeFontSize=10, NodeFontSize=10, ...
     MarkerSize=[5 20*sum(cellfun(@(x) sum(x),activities)) / sum(cellfun(@height,activities)) ...
                 20*sum(cellfun(@(x) sum(~x),activities)) / sum(cellfun(@height,activities)) 5])
title("Rate Constants" + newline + suptitle_suffix + " (n = " + height(activities) + ")", Interpreter="none")

savefig(f,save_dir + save_suffix + "_rates" + ".fig");
print(f,'-dtiff',save_dir + save_suffix + "_rates" + ".tif",'-r300');


%% Verify

psucs = (ka.*ki) ./ ( (ki+koa).*(ka+koi) );

acts = activities;
acts = acts(cellfun(@(x) x(1)==0,acts));
niai = cellfun(@(x) sum(diff(x)==-1), acts);

[y, x] = histcounts(categorical(niai)); % exposure*min_lifetime:exposure:1
x = str2double(x');
y = y';

[curve, gof] = fit(x,log(y),'poly1','Exclude', x>4);
psuc = exp(curve.p1);
rsquare = gof.rsquare;

f = figure;
hold on
bar(x,log(y))
plot(curve)
xticks(x)
text(x, log(y), num2str(y), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
text(0.9,0.8, ...
     "p_{suc}^{model} = " + round(psucs,2) + newline + ...
     "p_{suc}^{ver} = " + round(psuc,2) + newline + ...
     "R^2 = " + round(rsquare,2), ...
     'Units','normalized','HorizontalAlignment','right')
legend off
title("Verification" + newline + suptitle_suffix + " (n = " + height(activities) + ")", Interpreter="none")
xlabel('n')
ylabel('log(Frequency)')
xlim([-inf 7])
hold off

savefig(f,save_dir + save_suffix + "_verification" + ".fig");
print(f,'-dtiff',save_dir + save_suffix + "_verification" + ".tif",'-r300');

end




