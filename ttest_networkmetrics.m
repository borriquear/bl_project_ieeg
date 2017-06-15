function [ttest_resuls] = ttest_networkmetrics(Badjall)
%%ttest_networkmetrics ttest for the network metrics. 4 ttest one for each
%%network metric. each vector is values of network metric for each
%%different threshold
%INPUT: Badjall cell(patients, conditions, frequencies,2)
%OUTPUT: ttest results for each
[nbpats,nbconds, nbfreqs,two] = size(Badjall);
condloop = 1:2;
for pati=1:nbpats
    for condj=condloop
        for freqk=1:nbfreqs
            Badj = Badjall{pati, condj, freqk,1};
            thresholdv = Badjall{pati, condj, freqk,2};
            fprintf('Calling to gettopolicalvaluesfromBadj pat %d/%d\n', pati, nbpats);
            netmetlist_el = gettopolicalvaluesfromBadj(Badj,thresholdv);
            netmetlist{pati, condj, freqk} = netmetlist_el;
        end
    end
end
ttest_resuls = getthettest(netmetlist)
end

function [metricsall] = gettopolicalvaluesfromBadj(Badj, thresholdv)
%%gettopolicalvaluesfromBadj from the list of adjacency matrices created
%%one for each threshold, return the vector containing the values of the
%%metrics, B0(min(thresiold)) ..... B(max(theshold));
%%Clustering(min(thresiold))...Clustering(max(threshold))
% fi number of patients, col number of metrics
%OUTPUT: metricsall = Betti_v, clustering_v, wiringcost_v, pathlength_v
[fi, col] = size(Badj);
netmetlist = keys(Badj{2,1});
for i=1:col
    %'B0'    'clustering'    'density'    'pathlength'
    mapobj = Badj{2,i};
    Betti_v(i) = mapobj(netmetlist{1});
    clustering_v(i) = mapobj(netmetlist{2});
    wiringcost_v(i) = mapobj(netmetlist{3});
    pathlength_v(i) = mapobj(netmetlist{4});
end
metricsall = {Betti_v, clustering_v, wiringcost_v, pathlength_v};
end

function [ttest_resuls] = getthettest(netmetlist)
%getthettest calculate the ttest for a matrix patient, cond, freq containin
%{Betti_v, clustering_v, wiringcost_v, pathlength_v}
%MBec = [];MBeo=[];MCec = [];MCeo=[];MWec = [];MWeo=[];MPec = [];MPeo=[];
% i the alpha band
foi = 3; % (3. 5. 9. 16. 28. 50.) alpha =3
%{1x4 cell}
%condition_ec{1} [1x1295 double]    [1x1295 double]    [1x1295 double]    [1x1295 double]
%condition_ec{end}
condition_ec = netmetlist(:,1,foi);
condition_eo = netmetlist(:,2,foi);
%build network metrics containing all the patients 

for i=1:length(condition_ec)
    x_ec = condition_ec{i};
    x_eo = condition_eo{i};
    MBec{i,:} = x_ec{1};
    MBeo{i,:} = x_eo{1};
    MCec{i,:} = x_ec{2};
    MCeo{i,:} = x_eo{2};
    MWec{i,:} = x_ec{3};
    MWeo{i,:} = x_eo{3};
    MPec{i,:} = x_ec{4};
    MPeo{i,:} = x_eo{4};    
end
%build one vector column for all patients
B_ec_all = cat(1,MBec{1,:}',MBec{2,:}',MBec{3,:}',MBec{4,:}',MBec{5,:}',MBec{6,:}',MBec{7,:}',MBec{8,:}',MBec{9,:}',MBec{10,:}',MBec{11,:}');
B_eo_all = cat(1,MBeo{1,:}',MBeo{2,:}',MBeo{3,:}',MBeo{4,:}',MBeo{5,:}',MBeo{6,:}',MBeo{7,:}',MBeo{8,:}',MBeo{9,:}',MBeo{10,:}',MBeo{11,:}');

C_ec_all = cat(1,MCec{1,:}',MCec{2,:}',MCec{3,:}',MCec{4,:}',MCec{5,:}',MCec{6,:}',MCec{7,:}',MCec{8,:}',MCec{9,:}',MCec{10,:}',MCec{11,:}');
C_eo_all = cat(1,MCeo{1,:}',MCeo{2,:}',MCeo{3,:}',MCeo{4,:}',MCeo{5,:}',MCeo{6,:}',MCeo{7,:}',MCeo{8,:}',MCeo{9,:}',MCeo{10,:}',MCeo{11,:}');


W_ec_all = cat(1,MWec{1,:}',MWec{2,:}',MWec{3,:}',MWec{4,:}',MWec{5,:}',MWec{6,:}',MWec{7,:}',MWec{8,:}',MWec{9,:}',MWec{10,:}',MWec{11,:}');
W_eo_all = cat(1,MWeo{1,:}',MWeo{2,:}',MWeo{3,:}',MWeo{4,:}',MWeo{5,:}',MWeo{6,:}',MWeo{7,:}',MWeo{8,:}',MWeo{9,:}',MWeo{10,:}',MWeo{11,:}');

P_ec_all = cat(1,MPec{1,:}',MPec{2,:}',MPec{3,:}',MPec{4,:}',MPec{5,:}',MPec{6,:}',MPec{7,:}',MPec{8,:}',MPec{9,:}',MPec{10,:}',MPec{11,:}');
P_eo_all = cat(1,MPeo{1,:}',MPeo{2,:}',MPeo{3,:}',MPeo{4,:}',MPeo{5,:}',MPeo{6,:}',MPeo{7,:}',MPeo{8,:}',MPeo{9,:}',MPeo{10,:}',MPeo{11,:}');

%do ttest
[Bh,Bp,Bci,Bstats] = ttest2(B_ec_all, B_eo_all);
[Ch,Cp,Cci,Cstats] = ttest2(C_ec_all, C_eo_all);
[Wh,Wp,Wci,Wstats] = ttest2(W_ec_all, W_eo_all);
[Ph,Pp,Pci,Pstats] = ttest2(P_ec_all, P_eo_all);
ttest_resuls = [{Bh,Bp,Bci,Bstats} {Ch,Cp,Cci,Cstats} {Wh,Wp,Wci,Wstats} {Ph,Pp,Pci,Pstats}]

disp('Calculating kruskalwallis MCP for Betti numbers... ')
matrix = [[B_ec_all, B_eo_all ]];
[mcp_p_b,mcp_c_b,mcp_m_b]  = comparemultiplecomparison(matrix);
matrix = [[C_ec_all, C_eo_all ]];
disp('Calculating kruskalwallis MCP for Clustering ... ')
[mcp_p_clus,mcp_c_clus,mcp_m_clus]  = comparemultiplecomparison(matrix);
disp('Calculating kruskalwallis MCP for W ... ')
matrix = [[W_ec_all, W_eo_all ]];
[mcp_p_w,mcp_c_w,mcp_m_w]  = comparemultiplecomparison(matrix);
disp('Calculating kruskalwallis MCP for P ... ')
matrix = [[P_ec_all, P_eo_all ]];
[mcp_p_p,mcp_c_p,mcp_m_p]  = comparemultiplecomparison(matrix);


%do bonferroni correctin for multiple comparison problem
%https://es.mathworks.com/matlabcentral/fileexchange/4114-t-test-with-bonferroni-correction/content/ttest_bonf.m

%[Bh_bf,Bp_bf,BsigPairs_bf] = ttest_bonf(MBec,[1 2],0.05,0)

end
function [p,c,m] = comparemultiplecomparison(matrixB)
    
    [p,tbl,stats] = kruskalwallis(matrixB)
    disp('computing  MCP... ')
    figure
    [c,m] = multcompare(stats, 'ctype','bonferroni', 'display','on');
end