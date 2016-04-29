function [] = comparecorrelationmatrix
% comparecorrelationmatrix distances and metric of the
% correlation matrix created per subject and band
eegcondition = 'HYP';
centerfrequencies = {2, 6 , 10, 23.5, 40};
eegpatientl = { 'TWH027','TWH024','TWH028','TWH030', 'TWH031','TWH033','TWH034'};

corrmetricsperb =zeros(length(eegpatientl),length(centerfrequencies))
% make sure that globalFsDir is assigned
if ~exist('globalFsDir','var')
    fprintf('globalFsDir not found, loading it...')
    eval('global globalFsDir');
    myp = 'D:\BIAL PROJECT\patients\';
    eval(['globalFsDir=' 'myp']);
end
%file where to write the report
metrtoload = 'corrmat_metrics.txt';
corrmetricfile = fullfile(globalFsDir, metrtoload);
fileID = fopen(corrmetricfile,'w');
arrayofcorrmatrix = cell(length(eegpatientl),length(centerfrequencies))
for indpat=1:length(eegpatientl)
    [myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(eegpatientl{indpat},eegcondition);
    patpath = strcat(globalFsDir,eegpatientl{indpat});
    fprintf(fileID,'Patient= %s\n',eegpatientl{indpat});
    for indexfreq = 1:length(centerfrequencies) % delta, theta, alpha, beta, gamma
        centerfreq = centerfrequencies{indexfreq};
        fprintf(fileID,'Frequency= %s\n',num2str(centerfreq));
        fprintf('Calculating metric of the correlation matrix for patient %s, cond %s, band %s \n',eegpatientl{indpat},eegcondition,num2str(centerfreq));
        mattoload = strcat('powerconnectivity_freq_',num2str(centerfreq),'_',eegcondition,'_', eegpatientl{indpat},'_',eegdate,'_',eegsession,'.mat');
        fftfile = fullfile(patpath,'data','figures', mattoload);
        fprintf('\tOpening correlation matrix....%s\n',fftfile);
        matf= matfile(fftfile);
        corrmatpersubband = matf.corr_matrix;
        meantwo = mean2(corrmatpersubband);
        %detertwo = det(corrmatpersubband);
        fprintf('\tMean of correlation matrix =%s\n',num2str(meantwo));
        fprintf(fileID,'The mean of the corr matrix= %.4f\n',meantwo);
        corrmetricsperb(indpat,indexfreq) = meantwo
        %calculate distance between corr matrices
        arrayofcorrmatrix{indpat,indexfreq}= corrmatpersubband;
    end
    
end
calculatecorrmatricesdistances(arrayofcorrmatrix)
fclose(fileID);
end
function [] = calculatecorrmatricesdistances(arrayofcorrmatrix)
% calculatecorrmatricesdistances calculates the euclidean distance among
% all the correlation matrices of the array

end