%% Script that goes step by step processing the data
%0. the EEG objects for each patients need to be already created (EEGLab)

%% Create mat file with power FFT 
%1.Define list of patients and conditions to analyze
patientslist = {'TWH030','TWH031','TWH033', 'TWH034','TWH037','TWH038','TWH042','TWH043'};
conditionslist = {'EC_PRE', 'HYP'};%,'EC_POST'};
% Generate the mat filewith power fft 
ip = 1; ic =1;
for i =ip:length(patientslist)
    eegpatient = patientslist{i};
    for j = ic:length(conditionslist)
        eegcond = conditionslist{j};
        fprintf('Calculating mat file with power for Patient %s, Condition %s\n',eegpatient,eegcond );
        %plotsinglechannel = 0 to do not plot each channel
        matfilename = fftperperpatient(eegpatient, eegcond);
        fprintf('DONE: mat file with power Patient %s, Condition %s in %s\n',eegpatient,eegcond, matfilename)
    end
end
%% 2. Display the power spectra for the patients
ip = 1; ic =1;
patientslist = {'TWH030','TWH031','TWH033', 'TWH034','TWH037','TWH038','TWH042','TWH043'};
conditionslist = {'EC_PRE', 'HYP'};
powerspecmatrix = {};
powerfreqsindexes = {};
powerspecmatrix_freqbands = {};
xlimtodisplay = 30;
% plot Power Spectra (mean , std) ALL channels per patient/condition
for i =ip:length(patientslist)
    eegpatient = patientslist{i};
    for j = ic:length(conditionslist)
        eegcond = conditionslist{j};
        fprintf('Calling to loadpowerspectrum, Patient %s, Condition %s\n', eegpatient,eegcond );
        [frequencies,frexidx1, frexidx2, powerx, powerx_mean_perbands] = loadpowerspectrumperpatient(eegpatient,eegcond);
        frex = frequencies(frexidx1:frexidx2);
        powerfreqsindexes{i,j} = frex;
        powerspecmatrix{i,j} = powerx;
        powerspecmatrix_freqbands{i,j} = powerx_mean_perbands;
    end
end
% plot power spectra per ROIs per areas ALL patients condition
rois = {'HD','T','F', 'IH'};
plotpowerspectrumallpatientsROI(patientslist, conditionslist,powerspecmatrix,powerfreqsindexes, powerspecmatrix_freqbands, rois);
% Statistical significance between conditions in power spectra
plotstatisticalsignificance_powerspec(patientslist, conditionslist, powerspecmatrix, powerfreqsindexes, powerspecmatrix_freqbands, rois);
%% Spectrogram
fprintf('Calling to plotspectrogramperpatient, Patient %s, Condition %s\n', eegpatient,eegcond)
plotspectrogram(patientslist, conditionslist, rois)
%% create power based correlation matrix
% creates the correlation matrix with Spearman, powerbased 
patientslist = { 'TWH031','TWH033', 'TWH034','TWH037','TWH038','TWH042','TWH043'};
patientslist = { 'TWH034','TWH037','TWH038','TWH042','TWH043'};
conditionslist = {'EC_PRE', 'HYP'}; 
centerfrequencies = {2};%, 6 , 10, 23.5, 40};
createpowerbcorrmatrix(patientslist, conditionslist, centerfrequencies)
%% network analysis
% displaypowerconnectivity (display corr.matrx 
% needs the mat with the corr_matrix and the network (graphtheoryanalysis.m)
%displays the correlation matrix and the undirected network 
displaypowerconnectivity(patientslist, conditionslist, centerfrequencies);


%% Phase-based Analysis
%Calculate the Inter Siete Phase Clustering
% R = ||1/n \sum_t=1,n e^i(\phi_ch1,t - \phi_ch2,t) 
%Get the  bivariate phase difference for every 2 channels
% append the result to the fft file
calculatephasedifferences(patientslist, conditionslist)
