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

% centerfrequencies = {2, 6 , 10, 23.5, 40};
% currfreq = centerfrequencies{1};
% currcond = {'HYP', 'EO PRE', 'EC POST'};
% currcond = currcond{2};

%if nargin < 1
%     disp('Finding the corrrelation matrix prior to Load it...')
%     %correlation matrix for single patient
%     fprintf('Calculating Network metrics for Patient %s:\n\n',currentpat);
%     h = figure;
%     matfileid = 'powerconnectivity_freq_';
%     patpath = strcat(globalFsDir,currentpat);
%     [myfullname, EEG, channel_labels, patdate, patsession] = initialize_EEG_variables(currentpat,currcond);
%     mattoload = strcat(matfileid,num2str(currfreq),'_',currcond,'_', currentpat,'_',patdate,'_',patsession,'.mat');
%     fftfile = fullfile(patpath,'data','figures', mattoload);
%     fprintf('Opening correlation matrix....%s\n',fftfile);
%     matf= matfile(fftfile);
%     corrmatpersubband = matf.corr_matrix;
%     corrMatrix = corrmatpersubband; 
    %mattoload = strcat('networkmetrics_freq_',num2str(currfreq),'_',currcond,'_', currentpat,'_',patdate,'_',patsession,'.mat');

    %else
    mattoload = strcat('networkmetrics_freq_',num2str(currfreq),'_',currcond,'.mat');
%end
figure;
strNames = channel_labels(2:end); %delete 'Event' channel
%strNames= channel_labels;
if strcmp(typeofgraph,'undirectedw') ==1
    fprintf('Showing Graph for Undirected Weighted correlation matrix\n');
elseif strcmp(typeofgraph,'undirectedu')
    fprintf('Showing Graph for Undirected Binary correlation matrix\n');
    [threshold,nstds, corrMatrix] = calculatethresholdmatrix(corrMatrix);
    %     listofcorrmatrices(1) = corrMatrix;
    legendofmatrices = {currcond, currentpat, currfreq, channel_labels};
end
myColorMap = lines(length(corrMatrix));
circularGraph(corrMatrix,'Colormap',myColorMap,'Label',strNames);
%title(['Condition: ',currcond,', Frequency:', getgreeksymbolfreq(currfreq), ', Patient:',currentpat, ', Thres. = m+ n*std=', num2str(threshold)]);
title([currcond,',in', getgreeksymbolfreq(currfreq), ', Patient:',currentpat, ', Thresh. = m+ 1*std']);
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

function [threshold, nstds, corrMatrix] = calculatethresholdmatrix(corrMatrix)
%calculatethresholdmatrix(corrMatrix) returns the threshold and the
%threshold matrix given a correlation matrix
%IN: corrMatrix [0,1] if [-1,1] the function calculates abs(corrMatrix)
%OUT: threshold, new threshold matrix

meanmatrix = mean2(abs(corrMatrix)); stdmatrix = std2(abs(corrMatrix));
nstds = 1; %number of standard deviations
threshold = meanmatrix + nstds*stdmatrix;
fprintf('Corr matrix mean=%2.4f +(n)%d*(std)%2.4f = %2.4f\n', meanmatrix, nstds,stdmatrix, threshold);
corrMatrix(corrMatrix >  threshold) = 1;
corrMatrix(corrMatrix <= threshold) = 0;
end