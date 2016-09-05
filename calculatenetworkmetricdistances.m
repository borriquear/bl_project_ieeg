function [] = calculatenetworkmetricdistances(label)
%calculatenetworkmetricdistances calculate difference between network
%metrics for pair of conditions
global globalFsDir;
%global label;
globalFsDir = loadglobalFsDir();
global metricsnetwstr;
if strcmp(label, 'power') > 0
    fprintf('Calling to power calculatenetworkmetricdistances....\n');
    label = 'power_netmetrics.mat';
elseif strcmp(label, 'phase') > 0
    fprintf('Calling to phase calculatenetworkmetricdistances....\n');
    label = 'phase_netmetrics.mat';
end
filetoopen = fullfile(globalFsDir, label);
S = load(filetoopen);
calculatenetworkmetricdistances_aux(S, label);
end

function [] = calculatenetworkmetricdistances_aux(S, label)
%open file
%global globalFsDir;
%S = load(fullfile(globalFsDir, 'power_netmetrics.mat'));
%global label;
%[Lia1,Locb1] = ismember(A,B)
listoffreqs = [];
listofpats = [];
S = S.metricsnetwstr;
Srows = length(S);
%take first file to see how many conditions
nconds = 1;
i = 1;
while strcmp(S{1}{1}(1),S{i+1}{1}(1)) < 1
    nconds = nconds + 1;
    i = i+1;
end
totalcombinations = (Srows*Srows - Srows)/2;
totalrowofdifferences = {};
for i=1:Srows
    %S{i,1} {1x4 cell}    {1x16 cell}
    patient=  S{i}{1}(2);
    condition=  S{i}{1}(1);
    freq = S{i}{1}{3};
    channel_list = S{i}{1}{4}; channel_list = channel_list(2:end);
    corrmat = S{i, 4};
    thresh = S{i, 3};
    histo = S{i,2};
    nmetrics = S{i}{2};
    [Lia1] = ismember(freq,listoffreqs);
    tot_nb_edges = (pow2(length(corrmat),2) - length(corrmat))/2;
    %build list of frequencies
    if sum(Lia1(:)) == 0
        listoffreqs = [listoffreqs freq];
    end
    %build list of patients
    [rn,cn]=  find(strcmp(patient,listofpats));
    if isempty(rn) == 1
        listofpats = [listofpats patient];
    end
    for j=i+1:Srows
        patient_ot=  S{j}{1}(2);
        condition_ot=  S{j}{1}(1);
        freq_ot = S{j}{1}{3};
        channel_list_ot = S{j}{1}{4};channel_list = channel_list(2:end);
        corrmat_ot = S{j, 4};
        thresh_ot = S{j, 3};
        histo_ot = S{j,2};
        nmetrics_ot = S{j}{2};
        %nmetrics_ot = {degree_v, gtommatrix,matchingmatrix,density_coeff,clusteringcoeff_v,transcoeff,components_v,componentsizes_v,communityLovaffi_v,assort_coeff,richclub_v,coreperip_v,Dmatrix,charpathlength_coeff,BCvector_normalized,pagerank_vector}
        if strcmp(patient, patient_ot)==1 && freq == freq_ot && strcmp(condition, condition_ot) < 1
            conditiondoff= sprintf('%s-%s', condition{1}, condition_ot{1});
            kd = KLDiv(histo,histo_ot); jd = JSDiv(histo,histo_ot);
            nb_edges_norm_diff = (sum(corrmat(:)) - sum(corrmat_ot(:)))/tot_nb_edges;
            density_coeff_diff = nmetrics{4} - nmetrics_ot{4};
            transitivity_coeff_diff = nmetrics{6} - nmetrics_ot{6};
            assort_coeff_diff = nmetrics{10} - nmetrics_ot{10};
            charpathlength_coeff_diff = nmetrics{14} - nmetrics_ot{14};
            clustering_coeff_diff = mean(nmetrics{5}) - mean(nmetrics_ot{5});
            degree_coeff_diff = mean(nmetrics{1}) - mean(nmetrics_ot{1});
            pagerank_diff = mean(nmetrics{end}) - mean(nmetrics_ot{end});
            BC_norm_diff = mean(nmetrics{end-1}) - mean(nmetrics_ot{end-1});
            % add file in difference object
            rowofdifferences= [patient, freq, conditiondoff, kd, jd, nb_edges_norm_diff, degree_coeff_diff,density_coeff_diff, clustering_coeff_diff, transitivity_coeff_diff, assort_coeff_diff, charpathlength_coeff_diff,BC_norm_diff, pagerank_diff];
            totalrowofdifferences = [totalrowofdifferences;rowofdifferences];
        end
    end
end
fprintf('Finished calculating totalrowofdifferences\n')
disp(totalrowofdifferences)
%save totalrowofdifferences file
[resvector] = calculatedistancresults(totalrowofdifferences, listofpats, listoffreqs, label);
%resvector 1x3 (diffconds) each cell has freqsxpatsxvalues
%plot for specfici frequency and value matric
% f=1 delta.. f= length(listoffreqs) = gamma
%k=1=  kd, k=2= jd, 3 nb_edges_norm_diff, 4 degree_coeff_diff, 5 density_coeff_diff, 6 clustering_coeff_diff, 7 transitivity_coeff_diff, 8 assort_coeff_diff, 9 charpathlength_coeff_diff, 10 BC_norm_diff, 11 pagerank_diff]
%resvector, 1,f,k,label
for fq=1:length(listoffreqs)
    for k=1:11
        fprintf('Calling to plotdistanceresults fq=%s measure= %s \n', getgreeksymbolfreq(listoffreqs(fq)), getmetriclabelfromindex(k));
        plotdistanceresults(resvector, listofpats, listoffreqs, fq, k, label);
    end
    %correlation between pair of measures
    plotcorrelationpairs(resvector, listofpats, listoffreqs, fq, label)
end
end

function [resvector] = calculatedistancresults(totalrowofdifferences,listofpats, listoffreqs, label)
%freqs
%pats
%build result vector
lencols = length(totalrowofdifferences(1,:));
resultsperpatfq_ECEO = {}; resultsperpatfq_ECHYP = {}; resultsperpatfq_EOHYP = {};
for f=1:length(listoffreqs)
    for p=1:length(listofpats)
        resultsperpatfq = findrowfrompatfrq(totalrowofdifferences,listofpats,p,listoffreqs,f);
        for ifp = 1:size(resultsperpatfq,1)
            for k=4:lencols
                if strcmp(resultsperpatfq(ifp,3),'ECPRE-EOPRE') == 1
                    resultsperpatfq_ECEO(f,p,k-3)= resultsperpatfq(ifp,k);  %ECPRE-EOPRE #
                elseif strcmp(resultsperpatfq(ifp,3),'ECPRE-HYP') == 1
                    resultsperpatfq_ECHYP(f,p,k-3)= resultsperpatfq(ifp,k);  %ECPRE-EOPRE #
                elseif  strcmp(resultsperpatfq(ifp,3),'EOPRE-HYP') == 1
                    resultsperpatfq_EOHYP(f,p,k-3)= resultsperpatfq(ifp,k);  %ECPRE-EOPRE #
                else
                end
            end
        end
    end
end
resvector{1} = resultsperpatfq_ECEO;
resvector{2} = resultsperpatfq_ECHYP;
resvector{3} = resultsperpatfq_EOHYP;
end

function resultsperpatfq = findrowfrompatfrq(totalrowofdifferences,listofpats,p,listoffreqs,f)

lencols = length(totalrowofdifferences(1,:));
resultsperpatfq = [];
for i=1:length(totalrowofdifferences)
    row = totalrowofdifferences(i,:);
    if sum(strcmp(listofpats{p},row)) > 0 && isequal(listoffreqs(f), row{2})== 1
        %buil results for patient anf frequency
        %[patient, freq, conditiondoff, kd, jd, nb_edges_norm_diff, degree_coeff_diff,density_coeff_diff, clustering_coeff_diff, transitivity_coeff_diff, assort_coeff_diff, charpathlength_coeff_diff,BC_norm_diff, pagerank_diff]
        resultsperpatfq = [resultsperpatfq; row];
    end
end
end

function [] = plotdistanceresults(resvector, listofpats, listoffreqs, fq, k, label)

figure;
label  = strsplit(label, '_'); label = label{1};
band = listoffreqs(fq);
metric = getmetriclabelfromindex(k);
nbconds = length(resvector);
condl1 = resvector{1};
kvalues1 = condl1(fq,:,k);
kvaluesv1 = kvalues1{1:end};
plot(1:1:11, [kvalues1{1:end}], 'r+');
hold on;
condl2 = resvector{2};
kvalues2 = condl2(fq,:,k);
kvaluesv2 = kvalues2{1:end};
plot(1:1:11, [kvalues2{1:end}], 'gd');
hold on;
condl3 = resvector{3};
kvalues3 = condl3(fq,:,k);
kvaluesv3 = kvalues3{1:end};
plot(1:1:11, [kvalues3{1:end}], 'bs');
%end
legend('EC-EO','EC-HYP','EO-HYP');
ylabellab = sprintf('%s difference between 2 conditions', metric);
xlabel('Patients'), ylabel(ylabellab), grid on
%phase/powr corr band metric
sp = sprintf('%s corr fq=%.2f %s',label, band, metric);
title(sp);
end

function [] = plotcorrelationpairs(resvector, listofpats, listoffreqs, fq, label)
%k=1=  kd, k=2= jd, 3 nb_edges_norm_diff, 4 degree_coeff_diff, 5 density_coeff_diff, 6 clustering_coeff_diff, 7 transitivity_coeff_diff, 8 assort_coeff_diff, 9 charpathlength_coeff_diff, 10 BC_norm_diff, 11 pagerank_diff]
%initialize
kkd1=[]; jjsd1 = []; eedges1= []; ddegree1 = []; ddensity1= []; cclustering1 = [];  ttransitivity1 =[];  aassort1 = []; ppathl1 = []; ppager1 = [];cond1multiv=[];
kkd2=[]; jjsd2 = []; eedges2= []; ddegree2 = []; ddensity2= []; cclustering2 = [];  ttransitivity2 =[];  aassort2 = []; ppathl2 = []; ppager2 = [];cond2multiv=[];
kkd3=[]; jjsd3 = []; eedges3= []; ddegree3 = []; ddensity3= []; cclustering3 = [];  ttransitivity3 =[];  aassort3 = []; ppathl3 = []; ppager3 = [];cond3multiv=[];

%cond1 = ECPRE-EOPRE
condl1 = resvector{1};
kvalues1 = condl1(fq,:,:);
kd1= kvalues1(1,:,1)';
jsd1 = kvalues1(1,:,2)';
edges1 = kvalues1(1,:,3)';
degree1 = kvalues1(1,:,4)';
density1 = kvalues1(1,:,5)';
clustering1 = kvalues1(1,:,6)';
transitivity1 = kvalues1(1,:,7)';
assort1= kvalues1(1,:,8)';
pathl1= kvalues1(1,:,9)';
pager1= kvalues1(1,:,11)';
for i=1:length(listofpats)
    kkd1(i)= kd1{i}; jjsd1(i)= jsd1{i};eedges1(i) = edges1{i} ;ddegree1(i)= degree1{i}; ddensity1(i) = density1{i}; cclustering1(i)= clustering1{i};  ttransitivity1(i)= transitivity1{i};  aassort1(i)= assort1{i}; ppathl1(i)= pathl1{i}; ppager1(i) = pager1{i};
end
cond1multiv = [kkd1'  jjsd1' ddegree1' ddensity1'  cclustering1' ppathl1' ];
[R1,PValue1] = corrplot(cond1multiv,'tail','right');
tt= sprintf('Freq = %.2f, cond =%s',listoffreqs(fq), label );
title(tt);

condl2 = resvector{2};
kvalues2 = condl2(fq,:,:);
kd2= kvalues2(1,:,2)';
jsd2 = kvalues2(1,:,2)';
edges2 = kvalues2(1,:,3)';
degree2 = kvalues2(1,:,4)';
density2 = kvalues2(1,:,5)';
clustering2 = kvalues2(1,:,6)';
transitivity2 = kvalues2(1,:,7)';
assort2= kvalues2(1,:,8)';
pathl2= kvalues2(1,:,9)';
pager2= kvalues2(1,:,9)';
for i=1:length(listofpats)
    kkd2(i)= kd2{i}; jjsd2(i)= jsd2{i};eedges2(i) = edges2{i} ;ddegree2(i)= degree2{i}; ddensity2(i) = density2{i}; cclustering2(i)= clustering2{i};  ttransitivity2(i)= transitivity2{i};  aassort2(i)= assort2{i}; ppathl2(i)= pathl2{i}; ppager2(i) = pager2{i};
end
cond2multiv = [kkd2'  jjsd2'  ddegree2' ddensity2'  cclustering2' ppathl2' ];
[R2,PValue2] = corrplot(cond2multiv,'tail','right');

condl3 = resvector{3};
kvalues3 = condl3(fq,:,:);
kd3= kvalues3(1,:,3)';
jsd3 = kvalues3(1,:,3)';
edges3 = kvalues3(1,:,3)';
degree3 = kvalues3(1,:,4)';
density3 = kvalues3(1,:,5)';
clustering3 = kvalues3(1,:,6)';
transitivity3 = kvalues3(1,:,7)';
assort3= kvalues3(1,:,8)';
pathl3= kvalues3(1,:,9)';
pager3= kvalues3(1,:,9)';
for i=1:length(listofpats)
    kkd3(i)= kd3{i}; jjsd3(i)= jsd3{i};eedges3(i) = edges3{i} ;ddegree3(i)= degree3{i}; ddensity3(i) = density3{i}; cclustering3(i)= clustering3{i};  ttransitivity3(i)= transitivity3{i};  aassort3(i)= assort3{i}; ppathl3(i)= pathl3{i}; ppager3(i) = pager3{i};
end
cond3multiv = [kkd3'  jjsd3'  ddegree3' ddensity3'  cclustering3' ppathl3' ];
[R3,PValue3] = corrplot(cond3multiv,'tail','right');

end

function metric = getmetriclabelfromindex(k)
%k=1=  kd, k=2= jd, 3 nb_edges_norm_diff, 4 degree_coeff_diff, 5 density_coeff_diff, 6 clustering_coeff_diff, 7 transitivity_coeff_diff, 8 assort_coeff_diff, 9 charpathlength_coeff_diff, 10 BC_norm_diff, 11 pagerank_diff]
switch(k)
    case 1
        metric = 'KL';
    case 2
        metric = 'JS';
    case 3
        metric = '#edges';
    case 4
        metric = 'degree coeff';
    case 5
        metric = 'density coeff';
    case 6
        metric = 'clustering coeff';
    case 7
        metric = 'transitivity coeff';
    case 8
        metric = 'assort coeff';
    case 9
        metric = 'path length';
    case 10
        metric =  'BC norm';
    case 11
        metric =  'page rank';
    otherwise
        fprintf('ERROR')
end
end