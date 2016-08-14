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
            nb_edges_norm_diff = (sum(corrmat(:)) - sum(corrmat_ot(:)))/length(corrmat);
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
plotdistanceresults(resvector, listofpats, listoffreqs,length(listoffreqs), 6, label);
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

label  = strsplit(label, '_'); label = label{1};
band = listoffreqs(fq);
metric = getmetriclabelfromindex(k);
figure;
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