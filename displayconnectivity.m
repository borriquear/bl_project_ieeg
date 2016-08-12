function [] = displayconnectivity(label)
%dispplayconnectivity: displyas adj matrix and network 
% IN: label = power|phase

if strcmp(label, 'power') > 0 
    fprintf('Calling to displaypowerconnectivity....\n');
    displaypowerconnectivity();
elseif strcmp(label, 'phase') > 0 
    fprintf('Calling to displayphaseconnectivity....\n');
    displayphaseconnectivity
end
end

function displaypowerconnectivity()
% displaypowerconnectivity shows the correlation matrices for power
% from power_conn.mat
global globalFsDir;
globalFsDir = loadglobalFsDir();
pfile = fullfile(globalFsDir, 'powerconn_matrices');
matf = matfile(pfile);
powercorrelatiomatrices = matf.powerconn_matrix;
centerfrequencies = powercorrelatiomatrices.freqsl;
conditionslistl = powercorrelatiomatrices.conditionsl;
eegpatientl = powercorrelatiomatrices.patientsl;
power_mtxs = powercorrelatiomatrices.power_matrix;

%h = zeros(length(eegpatientl),length(conditionslistl),length(centerfrequencies));
h2 = zeros(length(eegpatientl),length(centerfrequencies));
patientslist = {'TWH030','TWH031','TWH033','TWH037','TWH038','TWH042','TWH043','TWH045','TWH047', 'TWH048','TWH049'};
powermats = {};
%for indpat=length(eegpatientl):length(eegpatientl)
for indpat=1:1
    curpatient = eegpatientl{indpat};
    powermats = power_mtxs(indpat,:,:);
    [myfullname, EEG, channel_labels, eegdate, eegsession] = initialize_EEG_variables(curpatient,'HYP');
    for j= 3:length(centerfrequencies)
        h2(indpat, j) = figure;
        displaypowercorrelationmatrix_allconds(powermats, curpatient,conditionslistl, centerfrequencies(j),j,h2(indpat, j), channel_labels)%,channel_labels);
        [allmetrics] = graphtheoryanalysis_power(powermats, curpatient,conditionslistl, centerfrequencies(j),j, channel_labels);
    end
end
end

function [] = displaypowercorrelationmatrix_allconds(powm, curpatient,conditionslistl, curfreq,indfreq, hfig, channel_labels)
% display one figure with all conditions, per patient and frequency
figure(hfig);
nbconds = length(conditionslistl);
vecconds = {'ECPRE', 'EOPRE', 'HYP', 'ECPOST', 'EOPOST'};
for i=1:nbconds
    subplot(1,nbconds,i)
    powpercf = powm{1,i,indfreq}
    powpercf = triu(powpercf,1); %get triangular superior
    powpercf = eye(length(powpercf))+ powpercf; %ones in the diagonal
    imagesc(powpercf);
    colormap('jet');
    colorbar;
    caxis([0 1])
    msgtitle = sprintf('Power Corr. matrix %s, frq=%.2f %s', curpatient,curfreq, vecconds{i});
    title(msgtitle);
end
end

function [allmetrics] = graphtheoryanalysis_power(corrMatrix, currentpat, conditionslistl, currfreq,indfreq, channel_labels)

legendofmatrices = {};
typeofgraph =  {'undirectedu','undirectedw','directed'};
indextype = 1; %1 undirected unweighted(0,1), undirected weighted, directed
typeofgraph = typeofgraph{indextype};
nbconds = length(conditionslistl);
vecconds = {'ECPRE', 'EOPRE', 'HYP', 'ECPOST', 'EOPOST'};
strNames = channel_labels(2:end); %delete 'Event' channel
if strcmp(typeofgraph,'undirectedw') ==1
    fprintf('Showing Graph for Undirected Weighted correlation matrix\n');
elseif strcmp(typeofgraph,'undirectedu')
    fprintf('Showing Graph for Undirected Binary correlation matrix\n');
end
nstds = 1; % number of std to calculate the threshold
for i=1:nbconds
    figure;
    powpercf = corrMatrix{1,i,indfreq}
    powpercf = triu(powpercf,1); %get triangular superior
    powpercf = eye(length(powpercf))+ powpercf; %ones in the diagonal
    if i == 1
        %fix the threshld for ECPRE
        thrd = mean(powpercf(:)) + nstds*(std(powpercf(:)));
    end
    [threshold, corrmatthres] = calculatethresholdmatrix(powpercf, nstds, thrd);
    fprintf('The threshold is mean + %d * std = %.2f\n', nstds,threshold);
    myColorMap = lines(length(corrmatthres));
    circularGraph(corrmatthres,'Colormap',myColorMap,'Label',strNames);
    %title(['Condition: ',currcond,', Frequency:', getgreeksymbolfreq(currfreq), ', Patient:',currentpat, ', Thres. = m+ n*std=', num2str(threshold)]);
    cond = vecconds(i); cond = cond{1};
    TIT = sprintf('Power Spearman corr. %s in fq=%.2f(%s) %s Thr= mean + %d*std= %.3f',cond, currfreq,getgreeksymbolfreq(currfreq), currentpat, nstds, threshold);
    %title([vecconds(i),',in', getgreeksymbolfreq(currfreq), ', Patient:',currentpat, ', Thresh. = m+ ',num2str(nstds),'*std']);
    title(TIT);
    drawnow update
    legendofmatrices = {vecconds(i), currentpat, currfreq, channel_labels};
    %schemaball(strNames, corrMatrix);
    [allmetrics] = calculategraphmatrixmetrics(corrmatthres, legendofmatrices);
    phasemetrics = {legendofmatrices, allmetrics};
    %save mat file
end
end

function displayphaseconnectivity()
% displayphaseconnectivity shows the correlation matrices for phase
% measures from phase_conn.mat
global globalFsDir;
globalFsDir = loadglobalFsDir();
pfile = fullfile(globalFsDir, 'phaseconn_matrices');
matf = matfile(pfile);
phasecorrlatiomatrices = matf.phaseconn_matrix;
centerfrequencies = phasecorrlatiomatrices.freqsl;
conditionslistl = phasecorrlatiomatrices.conditionsl;
eegpatientl = phasecorrlatiomatrices.patientsl;
ispc_mtxs = phasecorrlatiomatrices.ispc_matrix;
pli_mtxs = phasecorrlatiomatrices.pli_matrix;
icoh_mtxs = phasecorrlatiomatrices.icoh_matrix;
h = zeros(length(eegpatientl),length(conditionslistl),length(centerfrequencies));
h2 = zeros(length(eegpatientl),length(centerfrequencies));
patientslist = {'TWH030','TWH031','TWH033','TWH037','TWH038','TWH042','TWH043','TWH045','TWH047', 'TWH048','TWH049'};
plis = {};rs = {};pliscond = [];pliscond= [];pliscondfreq = []; rscondfreq = [];
%for indpat=length(eegpatientl):length(eegpatientl)
for indpat=1:1
    curpatient = eegpatientl{indpat};
    plis = pli_mtxs(indpat,:);
    rs = ispc_mtxs(indpat,:);
    [myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(curpatient,'HYP');
    %     for indcond=1:length(conditionslistl)
    %         curcondition = conditionslistl{indcond};
    %         pliscond = plis(indcond);
    %         rscond = rs(indcond);
    %         [myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(curpatient,curcondition);
    for j= 3:length(centerfrequencies)
        %display correlation matrix
        %one figure for each freq band for corrplot
        %h(indpat,indcond, j) = figure;
        h2(indpat, j) = figure;
        %            pliscondfreq = pliscond{1}; pliscondfreq = pliscondfreq(:,:,j);
        %            rscondfreq = rscond{1}; rscondfreq = rscondfreq(:,:,j);
        %displayphasecorrelationmatrix(pliscondfreq, rscondfreq, curpatient,curcondition, centerfrequencies(j),h(indpat,indcond, j),channel_labels);
        displayRphasecorrelationmatrix_allconds( rs, curpatient,conditionslistl, centerfrequencies(j),j,h2(indpat, j), channel_labels)%,channel_labels);
        
        % generate the network
        [allmetrics] = graphtheoryanalysis_phase(rs, curpatient,conditionslistl, centerfrequencies(j),j, channel_labels);
    end
    % end
end
end


function [] = displayphasecorrelationmatrix(rs,curpatient,curcondition,curfreq,hfig,channels_patl)
% display one figure with allmeasures, per patient and frequency and
% condition
figure(hfig);
subplot(1,2,1)
imagesc(rs(1));
colormap('jet');
colorbar;
caxis([0 1])
msgtitle = sprintf('PLI Phase Correlation matrix %s, %s, frq=%.2f', curpatient,curcondition,curfreq);
title(msgtitle);
subplot(1,2,2)
imagesc(rs(2));
colormap('jet');
colorbar;
caxis([0 1])
msgtitle = sprintf('R Phase Correlation matrix %s, %s, frq=%.2f', curpatient,curcondition,curfreq);
title(msgtitle);

%fprintf('Patient %s, Cond %s, Freq %d. The mean in || of the correlation amtrix= %.3f \n',curpatient,curcondition,curfreq,  mean(abs(corr_mat(:))))
end

function [] = displayRphasecorrelationmatrix_allconds(rs, curpatient,conditionslistl, curfreq,indfreq, hfig, channel_labels)
% display one figure with all conditions, per patient and frequency
figure(hfig);
nbconds = length(conditionslistl);
vecconds = {'ECPRE', 'EOPRE', 'HYP', 'ECPOST', 'EOPOST'};
for i=1:nbconds
    
    subplot(1,nbconds,i)
    rsvec= rs(i);
    rsvec = rsvec{1};
    vecrs = rsvec(:,:,indfreq);
    imagesc(vecrs);
    colormap('jet');
    colorbar;
    caxis([0 1])
    msgtitle = sprintf('R Corr. matrix %s, frq=%.2f %s', curpatient,curfreq, vecconds{i});
    title(msgtitle);
end
end


function [allmetrics] = graphtheoryanalysis_phase(corrMatrix, currentpat, conditionslistl, currfreq,indfreq, channel_labels)

global globalFsDir;
[globalFsDir] = loadglobalFsDir();
legendofmatrices = {};
typeofgraph =  {'undirectedu','undirectedw','directed'};
indextype = 1; %1 undirected unweighted(0,1), undirected weighted, directed
typeofgraph = typeofgraph{indextype};
nbconds = length(conditionslistl);
vecconds = {'ECPRE', 'EOPRE', 'HYP', 'ECPOST', 'EOPOST'};
strNames = channel_labels(2:end); %delete 'Event' channel
if strcmp(typeofgraph,'undirectedw') ==1
    fprintf('Showing Graph for Undirected Weighted correlation matrix\n');
elseif strcmp(typeofgraph,'undirectedu')
    fprintf('Showing Graph for Undirected Binary correlation matrix\n');
end
nstds = 1; % number of std to calculate the threshold
for i=1:nbconds
    figure;
    rsvec= corrMatrix(i);
    rsvec = rsvec{1};
    vecrs = rsvec(:,:,indfreq);
    if i == 1
        %fix the threshld for ECPRE
        thrd = mean(vecrs(:)) + nstds*(std(vecrs(:)));
    end
    [threshold, corrmatthres] = calculatethresholdmatrix(vecrs, nstds, thrd);
    fprintf('The threshold is mean + %d * std = %.2f\n', nstds,threshold);
    myColorMap = lines(length(corrmatthres));
    circularGraph(corrmatthres,'Colormap',myColorMap,'Label',strNames);
    %title(['Condition: ',currcond,', Frequency:', getgreeksymbolfreq(currfreq), ', Patient:',currentpat, ', Thres. = m+ n*std=', num2str(threshold)]);
    cond = vecconds(i); cond = cond{1};
    TIT = sprintf('Phase corr. %s in fq=%.2f(%s) %s Thr= mean + %d*std= %.3f',cond, currfreq,getgreeksymbolfreq(currfreq), currentpat, nstds, threshold);
    %title([vecconds(i),',in', getgreeksymbolfreq(currfreq), ', Patient:',currentpat, ', Thresh. = m+ ',num2str(nstds),'*std']);
    title(TIT);
    drawnow update
    legendofmatrices = {vecconds(i), currentpat, currfreq, channel_labels};
    %schemaball(strNames, corrMatrix);
    [allmetrics] = calculategraphmatrixmetrics(corrmatthres, legendofmatrices);
    phasemetrics = {legendofmatrices, allmetrics};
    %save mat file
end
end