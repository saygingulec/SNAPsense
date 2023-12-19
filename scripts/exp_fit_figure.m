function kd = exp_fit_figure(lt, tit)
%%

warning('off','curvefit:fit:noStartPoint')
    
% Survival funciton
[y, x] =  ecdf(lt,'function','survivor');
x = x(2:end);
y = y(1:end-1);

% Exponential fit
[curve, gof] = fit(x,y,'exp(b*x)','Exclude',x>0.4);
kd = -curve.b;
rsquare = gof.rsquare;

% % 2 Exp
% [curve, gof] = fit(x,y,'f*exp(b*x)+(1-f)*exp(d*x)','Exclude',x>0.6);
% kda1 = -curve.b;
% kda2 = -curve.d;
% rsquare = gof.rsquare;

% Plot
f = figure;
hold on
ecdf(lt,'function','survivor','Bounds','on')
plot(curve,'predfun')
title(tit + " (n = " + height(lt) + ")", Interpreter="none")
xlabel('Time (s)')
ylabel('Survival')
text(0.9,0.8,"k_d = " + round(kd,2) + newline + "R^2 = " + round(rsquare,2), ...
    'Units','normalized','HorizontalAlignment','right')
% text(0.9,0.8,"k_{da1} = " + round(kda1,2) + newline + "k_{da2} = " + round(kda2,2) + newline + "R^2 = " + round(rsquare,2), ...
%     'Units','normalized','HorizontalAlignment','right')
legend off
xlim([x(1) 0.4])
ylim([0 1])
pbaspect([1 1 1])
hold off
end