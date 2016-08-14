function [] = displayconnectivity(label)
%dispplayconnectivity: displyas adj matrix and network
% IN: label = power|phase
%clear all;
global globalFsDir;
globalFsDir = loadglobalFsDir();
%global metricsnetwstr;
clear metricsnetwstr;
if strcmp(label, 'power') > 0
    fprintf('Calling to displaypowerconnectivity....\n');
    displaypowerconnectivity();
elseif strcmp(label, 'phase') > 0
    fprintf('Calling to displayphaseconnectivity....\n');
    displayphaseconnectivity();
end
end

function displaypowerconnectivity()
% displaypowerconnectivity shows the correlation matrices for power
% from power_conn.mat
global globalFsDir;
%global metricsnetwstr;
metricsnetwstr_ij = {}; metricsnetwstr = {};
allmetrics_ij = {}; allmetrics = {};
pfile = fullfile(globalFsDir, 'powerconn_matrices');
matf = matfile(pfile);
powercorrelatiomatrices = matf.powerconn_matrix;
centerfrequencies = powercorrelatiomatrices.freqsl;
conditionslistl = powercorrelatiomatrices.conditionsl;
eegpatientl = powercorrelatiomatrices.patientsl;
power_mtxs = powercorrelatiomatrices.power_matrix;

%h = zeros(length(eegpatientl),length(conditionslistl),length(centerfrequencies));
h2 = zeros(length(eegpatientl),length(centerfrequencies));
%patientslist = {'TWH030','TWH031','TWH033','TWH037','TWH038','TWH042','TWH043','TWH045','TWH047', 'TWH048','TWH049'};
powermats = {};
%for indpat=length(eegpatientl):length(eegpatientl)
for indpat=1:length(eegpatientl)
    curpatient = eegpatientl{indpat};
    powermats = power_mtxs(indpat,:,:);
    [myfullname, EEG, channel_labels, eegdate, eegsession] = initialize_EEG_variables(curpatient,'HYP');
    for j= 3:length(centerfrequencies)
        h2(indpat, j) = figure;
        displaypowercorrelationmatrix_allconds(powermats, curpatient,conditionslistl, centerfrequencies(j),j,h2(indpat, j), channel_labels)%,channel_labels);
        [allmetrics_ij, metricsnetwstr_ij] = graphtheoryanalysis_power(powermats, curpatient,conditionslistl, centerfrequencies(j),j, channel_labels);
        allmetrics = [allmetrics; allmetrics_ij];
        metricsnetwstr = [metricsnetwstr; metricsnetwstr_ij];
        %update object of network metrics
        %pause(1);
        %close(h2(indpat, j));
    end
end
%metricsnetwstr
%creatematfileofnetworkmetrics(powermats, curpatient,conditionslistl, centerfrequencies(j),j, channel_labels, allmetrics)
%nbofrows = length(metricsnetwstr); %pats*conds*freqs
% metricsnetwstr(1) {1x2 cell}
% metricsnetwstr{1} {1x4 cell}    {1x16 cell} %legend metric
% metricsnetwstr{i,2} = [ 0's 1's]  metricsnetwstr{i,3} = [thr]  metricsnetwstr{i,3} = [matrix]
matfnetm = fullfile(globalFsDir,'power_netmetrics.mat' );
fprintf('Saving power based network metrics in file  .....\n');
save(matfnetm, 'metricsnetwstr');
disp(metricsnetwstr)
fprintf('DONE! power based network metrics in file %s \n', matfnetm);

end

function [] = displaypowercorrelationmatrix_allconds(powm, curpatient,conditionslistl, curfreq,indfreq, hfig, channel_labels)
% display one figure with all conditions, per patient and frequency
figure(hfig);
nbconds = length(conditionslistl);
%global metricsnetwstr;
vecconds = {'ECPRE', 'EOPRE', 'HYP', 'ECPOST', 'EOPOST'};
for i=1:nbconds
    subplot(1,nbconds,i)
    powpercf = powm{1,i,indfreq};
    powpercf = triu(powpercf,1); %get triangular superior
    powpercf = eye(length(powpercf))+ powpercf; %ones in the diagonal
    imagesc(powpercf);
    colormap('jet');
    colorbar;
    caxis([0 1])
    msgtitle = sprintf('Power Corr. matrix %s, frq=%.2f %s', curpatient,curfreq, vecconds{i});
    title(msgtitle);
end
pause(1);
close(hfig);
end

function [allmetrics,metricsnetwstr] = graphtheoryanalysis_power(corrMatrix, currentpat, conditionslistl, currfreq,indfreq, channel_labels)
%graphtheoryanalysis_power calculate network quantities (clustering etc)
%for a corr mat, note that corrMatrix is avectorof matrices one for each
%condition
legendofmatrices = {};
metricsnetwstr = {};
%global metricsnetwstr;
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
hfig = figure;
for i=1:nbconds
    
    subplot(1,nbconds,i);
    hold on;
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
    %
    legendofmatrices = {vecconds{i}, currentpat, currfreq, channel_labels};
    %schemaball(strNames, corrMatrix);
    [allmetrics] = calculategraphmatrixmetrics(corrmatthres, legendofmatrices);
    legendofmatrices = {legendofmatrices, allmetrics}
    %save mat file
    % calculate the histograme for the distances
    %nbins = 3; [hp, ep] = calculatehistogram(corrmatthres, nbins);
    [histog, ep] = calculatehistogram(corrmatthres);
    phasemetrics = {legendofmatrices, histog, threshold, corrmatthres };
    metricsnetwstr = [metricsnetwstr ; phasemetrics];
end
drawnow update
pause(1);
close(hfig);
end

function displayphaseconnectivity()
% displayphaseconnectivity shows the correlation matrices for phase
% measures from phase_conn.mat
global globalFsDir;
globalFsDir = loadglobalFsDir();
pfile = fullfile(globalFsDir, 'phaseconn_matrices');
matf = matfile(pfile);
metricsnetwstr = {};
phasecorrlatiomatrices = matf.phaseconn_matrix;
centerfrequencies = phasecorrlatiomatrices.freqsl;
conditionslistl = phasecorrlatiomatrices.conditionsl;
eegpatientl = phasecorrlatiomatrices.patientsl;
ispc_mtxs = phasecorrlatiomatrices.ispc_matrix;
pli_mtxs = phasecorrlatiomatrices.pli_matrix;
icoh_mtxs = phasecorrlatiomatrices.icoh_matrix;
h = zeros(length(eegpatientl),length(conditionslistl),length(centerfrequencies));
h2 = zeros(length(eegpatientl),length(centerfrequencies));
%patientslist = {'TWH030','TWH031','TWH033','TWH037','TWH038','TWH042','TWH043','TWH045','TWH047', 'TWH048','TWH049'};
plis = {};
rs = {}; pliscond = [];pliscond= [];pliscondfreq = []; rscondfreq = [];
allmetrics_ij = {}; allmetrics = {}; metricsnetwstr = {}; metricsnetwstr_ij = {};
%for indpat=length(eegpatientl):length(eegpatientl)
for indpat=1:length(eegpatientl)
    curpatient = eegpatientl{indpat};
    plis = pli_mtxs(indpat,:);
    rs = ispc_mtxs(indpat,:);
    [myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(curpatient,'HYP');
    for j= 3:length(centerfrequencies)
        h2(indpat, j) = figure;
        displayRphasecorrelationmatrix_allconds( rs, curpatient,conditionslistl, centerfrequencies(j),j,h2(indpat, j), channel_labels)%,channel_labels);
        % generate the network
        [allmetrics_ij, metricsnetwstr_ij] = graphtheoryanalysis_phase(rs, curpatient,conditionslistl, centerfrequencies(j),j, channel_labels);
         allmetrics = [allmetrics; allmetrics_ij];
         metricsnetwstr = [metricsnetwstr; metricsnetwstr_ij];
    % end
end
matfnetm = fullfile(globalFsDir,'phase_netmetrics.mat' );
disp(metricsnetwstr);
fprintf('Saving phase based network metrics in file  .....\n');
save(matfnetm, 'metricsnetwstr');
fprintf('DONE! phase based network metrics in file %s \n', matfnetm);
end
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
pause(1);
close(hfig);
end


function [allmetrics, metricsnetwstr] = graphtheoryanalysis_phase(corrMatrix, currentpat, conditionslistl, currfreq,indfreq, channel_labels)

%global globalFsDir;
%global metricsnetwstr;
legendofmatrices = {};
phasemetrics = {};
metricsnetwstr = {};
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
hh= figure;
for i=1:nbconds
    subplot(1, nbconds,i);
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
    legendofmatrices = {vecconds{i}, currentpat, currfreq, channel_labels};
    %schemaball(strNames, corrMatrix);
    [allmetrics] = calculategraphmatrixmetrics(corrmatthres, legendofmatrices);
    legendofmatrices = {legendofmatrices, allmetrics}
    %save mat file
    % calculate the histograme for the distances
    %nbins = 3; [hp, ep] = calculatehistogram(corrmatthres, nbins);
    [histog, ep] = calculatehistogram(corrmatthres);
    phasemetrics = {legendofmatrices, histog, threshold, corrmatthres };
    metricsnetwstr = [metricsnetwstr ; phasemetrics];
end
pause(1);
close(hh);
end
