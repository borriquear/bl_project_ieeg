function [metricsall] = plotnetworkmetrics(hdl, BM, delta, contrast)
%% plotnetworkmetrics
%INPUT: hdl figure handler, Badj{1,:} Binary matrices, one for each threshold
%Badj{2,:} network metrics, one for each threshold, delta is the vector of
%threshold values
%OUTPUT:
if nargin > 3
    linespec = '--';
else
    linespec = '-';
end
[fi, col] = size(BM);
netmetlist = keys(BM{2,1});
for i=1:col
    mapobj = BM{2,i};
    Betti_v(i) = mapobj(netmetlist{1});
    clustering_v(i) = mapobj(netmetlist{2});
    density_v(i) = mapobj(netmetlist{3});
    pathlength_v(i) = mapobj(netmetlist{4});
%     transitivity_v(i) = mapobj(netmetlist{5});
%     transitivity_v(isnan(transitivity_v)) = 0 ;
end

figure(hdl)
set(gca, 'XTick',[min(delta) mean(delta) max(delta)])
%ylim([0 1]);
xlabel('delta'), ylabel('network metrics')
msgt = sprintf('Network metrics for BN threshold');
%legend(conditionslist{1}, conditionslist{2}, 'Interpreter', 'None')
%title(msgt,'interpreter', 'none');

plot(delta,Betti_v,'Color','m','LineStyle',linespec);hold on;
plot(delta,clustering_v,'Color','r','LineStyle',linespec);hold on;
plot(delta,density_v,'Color','b','LineStyle',linespec);hold on;
plot(delta,pathlength_v,'Color','g','LineStyle',linespec);hold on;
% plot(delta,transitivity_v);
legend('Betti 0','clustering', 'density','path length')
metricsall = {Betti_v, clustering_v,density_v,pathlength_v};
end