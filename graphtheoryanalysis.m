%THIS FUNCTION IS DEPRECATED. _power
% use graphtheoryanalysis_phase and graphtheoryanalysis
function [allmetrics] = graphtheoryanalysis(corrMatrix, currentpat, currcond, currfreq, channel_labels)
%graphtheoaryanalysis  build a network from corr_matrix in powerconnectivity_freq_*.mat and calculates network 
% IN: optional argument. corrMatrix == []  get the correlation matrix
% : corrMatrix != [] correlation matrix from which to obtain the
% network , from the mat file powerconnectivity
% set currfreq
global globalFsDir;
[globalFsDir] = loadglobalFsDir();
% global currfreq;
% global currcond;
legendofmatrices = {};
allmetricpercfp = [];
% patientsList = {'TWH024','TWH027','TWH028','TWH030','TWH031','TWH033','TWH034'};
% currentpat = patientsList{4};
typeofgraph =  {'undirectedu','undirectedw','directed'};
indextype = 1; %1 undirected unweighted(0,1), undirected weighted, directed
typeofgraph = typeofgraph{indextype};
mattoload = strcat('networkmetrics_freq_',num2str(currfreq),'_',currcond,'.mat');
%end
figure;
strNames = channel_labels(2:end); %delete 'Event' channel
%strNames= channel_labels;
if strcmp(typeofgraph,'undirectedw') ==1
    fprintf('Showing Graph for Undirected Weighted correlation matrix\n');
elseif strcmp(typeofgraph,'undirectedu')
    fprintf('Showing Graph for Undirected Binary correlation matrix\n');
    nstds =2;
    [threshold, corrMatrix] = calculatethresholdmatrix(corrMatrix, nstds);
    %     listofcorrmatrices(1) = corrMatrix;
    legendofmatrices = {currcond, currentpat, currfreq, channel_labels};
end
myColorMap = lines(length(corrMatrix));
circularGraph(corrMatrix,'Colormap',myColorMap,'Label',strNames);
%title(['Condition: ',currcond,', Frequency:', getgreeksymbolfreq(currfreq), ', Patient:',currentpat, ', Thres. = m+ n*std=', num2str(threshold)]);
title([currcond,',in', getgreeksymbolfreq(currfreq), ', Patient:',currentpat, ', Thresh. = m+ ',num2str(nstds),'*std']);
drawnow update
%schemaball(strNames, corrMatrix);
[allmetrics] = calculategraphmatrixmetrics(corrMatrix, legendofmatrices);
allmetrics = {legendofmatrices, allmetrics};
%allmetricpercfp = [allmetrics; allmetricpercfp]; %for more than 1 patient
%save mat file
[globalFsDir] =loadglobalFsDir();
patpath = strcat(globalFsDir,currentpat);
%mattoload = strcat('networkmetrics_freq_',num2str(currfreq),'_',currcond,'_', currentpat,'_',patdate,'_',patsession,'.mat');
fftfile = fullfile(patpath,'data','figures', mattoload);
save(fftfile,'allmetrics')
end