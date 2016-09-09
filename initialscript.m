
%% Script that goes step by step processing the data
%0. The EEG objects for each patients need to be already created (EEGLab)
% cuteoneepochNEW.m, read overleaf documenttation        
% Patients conditions
% TWH024 =               HYP
% TWH027 =               HYP
% TWH028 =               HYP EC_POST
% TWH030 = EC_PRE EO_PRE HYP EC_POST  
% TWH031 = EC_PRE EO_PRE HYP EC_POST  
% TWH033 = ALL
% TWH034 = EC_PRE        HYP EC_POST EO_POST
% TWH035 = EC_PRE        HYP EC_POST EO_POST
% TWH037 = ALL
% TWH038 = ALL
% TWH042 = ALL
% TWH043 = ALL
% TWH045 = EC_PRE EO_PRE HYP
% TWH047 = ALL
% TWH048 = ALL
% TWH049 = ALL
%To save figures 
%>>print -f3 -djpeg 'D:\BIAL PROJECT\patients\figure_results\ECEO-Allroisfigi1'
%>>savefig('D:\BIAL PROJECT\patients\figure_results\ECEO-Allroisfigi1')
%% Create mat file with power/amplitud calculated via the FFT. Power analysis, identify channels with most power and the frequency bands that pick up maximum power
%OUTPUT: % globalFsDir\eegpatient\data\figures\fft_%pat_%cond_%date_%se.mat -> 'ampli_fft','power_fft','power_fft_perband', 'power_fft_mean_perband','channel_labels'
global globalFsDir;
globalFsDir = loadglobalFsDir();
patientslist = {'TWH030','TWH031','TWH033','TWH037','TWH038','TWH042','TWH043','TWH045','TWH047', 'TWH048','TWH049'};
%patientslist = {'TWH047', 'TWH048'};
conditionslist = {'EC_PRE', 'EO_PRE', 'HYP', 'EC_POST','EO_POST'};%,'EC_POST'};
conditionslist = {'EC_PRE', 'EO_PRE', 'HYP'};
% Generate the mat file (fft_%s_%s_%s_%s.mat',eegcond, eegpatient, eegdate, eegsession)
%power_fft_perband (1x50) with the power for each frequency from 1 to 51
%power_fft_mean_perband (nbchanneslx5) for each channel the % of power
%relative to the other bands. The sum of each row is 1.
%power_fft (nbchanneslxtimepoints)
% that calculates amplitude and power from the
% fft to find the frequency components of the signal (patient, condition)
% ranges of frequency f = 0:50 (bins);hz = linspace(0,nyquistfreq,floor(n/2)+1);
%For a given signal, the power spectrum gives the portion of a
% signal's power (energy per unit time) falling within given frequency bins.
%FFT fft(X) computes the discrete Fourier transform (DFT) of X using a fast Fourier transform (FFT) algorithm.
%     signalXF = fft(signal)/n;
%     ampli_fft(irow-1,:) = 2*abs(signalXF(1:length(hz)));
%     power_fft(irow-1,:) = abs(signalXF(1:length(hz))).^2;
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
%%  1.2. Display the power spectra for the patients 
% Load the (fft_%s_%s_%s_%s.mat',eegcond, eegpatient, eegdate, eegsession)
% for each patient,condition
ip = 1; ic =1;
patientslist = {'TWH030','TWH031','TWH033','TWH037','TWH038','TWH042','TWH043','TWH045','TWH047', 'TWH048','TWH049'};
%conditionslist = {'EC_PRE', 'EO_PRE', 'HYP', 'EC_POST','EO_POST'};
conditionslist = {'EC_PRE', 'EO_PRE'}%, 'HYP'};

powerspecmatrix = {};
powerfreqsindexes = {};
powerspecmatrix_freqbands = {};
%Set xlimtodisplay = 50 in plotpowerspectrumallpatientsROI.m
%plot Power Spectra (mean, std) ALL channels per patient/condition
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
%rois = 'All' plot all areas, rois ~= All, plot all Areas and Region of
%interest, eg HD
roislist = {'All','T','F','FP','IH','Grid', 'HD','NOHD', 'D', 'BiTemp'};
rois = roislist(8);
plotpowerspectrumallpatientsROI(patientslist, conditionslist, powerspecmatrix, powerfreqsindexes, powerspecmatrix_freqbands, rois);

% Statistical significance between conditions in power spectra.
% Calculate ttest to compare
plotstatisticalsignificance_powerspec(patientslist, conditionslist, powerspecmatrix, powerfreqsindexes, powerspecmatrix_freqbands, rois);
%% 1.3 Spectrogram
fprintf('Calling to plotspectrogramperpatient, Patient %s, Condition %s\n', eegpatient,eegcond)
plotspectrogram(patientslist, conditionslist, rois)
%% 1.4 create power based correlation matrix
% creates the correlation matrix with Spearman, powerbased 
patientslist = {'TWH030','TWH031','TWH033','TWH037','TWH038','TWH042','TWH043','TWH045','TWH047', 'TWH048','TWH049'};
conditionslist = {'EC_PRE', 'EO_PRE', 'HYP'};%,'EC_POST'};
centerfrequencies  = logspace(log10(1),log10(50),8);
%create powerconn matrices, first for segments of temporalwindow secs in powerconn_matrices.mat
%1. Create powerconn_matrices_tw.mat for temporalwindow >0 (cut the signal into segments)
temporalwindow = 4; % 4 seconds powerconn_matrices_tw.mat
createpowerbcorrmatrix(patientslist, conditionslist, centerfrequencies, temporalwindow);
powerconnmatf = fullfile(globalFsDir, 'powerconn_matrices_tw.mat');
%temporalwindow: 0 means correlation calculated with entire epoch length
fprintf('Created file of power connectivity in %s\n', powerconnmatf);
%2.Create powerconn_matrices.mat now for for temporalwindow = 0 , entire epoch
createpowerbcorrmatrix(patientslist, conditionslist, centerfrequencies);
powerconnmatf = fullfile(globalFsDir, 'powerconn_matrices.mat');
fprintf('Created file of power connectivity in %s\n', powerconnmatf);
%% 1.5 network analysis
%displays the correlation matrix and the undirected network 
%[srate, min_freq, max_freq, num_frex, time, n_wavelet, half_wavelet, freqs2use, s, wavelet_cycles]= initialize_wavelet();
%centerfrequenciesl  = logspace(log10(min_freq),log10(max_freq),8)
% centerfrequencies  = logspace(log10(1),log10(50),8);
% centerfrequencies  = centerfrequencies(3:end);
%display power connectivity (adj matrix and network ) from "power_conn.mat"

%INPUT: correlation matrix mat file
%OUTPUT: power_netmetrics.mat, phase_netmetrics.mat
displayconnectivity('power');
%Creates metric distance file in D:\BIAL PROJECT\patients\power_netmetrics.mat
%calculate network metric distances between pair of conditions
%INPUT:power_netmetrics.mat
%OUTPUT:
calculatenetworkmetricdistances('power');
%%%%%%%%%%%%%%%%%%%%%%%%%%% Phase-based Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Phase-based Analysis
% Calculate the Inter Site Phase Clustering
% R = ||1/n \sum_t=1,n e^i(\phi_ch1,t - \phi_ch2,t) 
% Get the  bivariate phase difference for every 2 channels
% append the result to the fft file
patientslist = {'TWH030','TWH031','TWH033','TWH037','TWH038','TWH042','TWH043','TWH045','TWH047', 'TWH048','TWH049'};
conditionslist = {'EC_PRE', 'EO_PRE', 'HYP'};%,'EC_POST'};
centerfrequencies  = logspace(log10(1),log10(50),8);
patientslist = {'TWH030'}; conditionslist = {'HYP'}; centerfrequencies = 1;
%temporalwindow in miliseconds
temporalwindow = 5;% temporalwindow: 5 argumentin phase_conn_tw.mat
calculatephasedifferences(patientslist, conditionslist, temporalwindow*1000);
matfilename = fullfile(globalFsDir, 'phaseconn_matrices.mat');
fprintf('Created file of phase connectivity in %s\n', matfilename);
%display phase connectivity , adj matrix and network, from "phase_conn.mat"
%displayphaseconnectivity();
displayconnectivity('phase');
%Creates metric distance file in D:\BIAL PROJECT\patients\phase_netmetrics.mat 
calculatenetworkmetricdistances('phase');
%% Wiring distance
%IN: phaseconn_matrices_tw.mat, 'powerconn_matrices_tw.mat
%OUTPUT: wiringcost_matrices.mat
%calls to calculateEuclideandistancematrix
wiring_matrices = calculatewiringcostmatrices();
%disp(wiring_matrices);
%plot the wiring cost matrices
%Plot ALL electrodes or a subset
plotwiringcost(wiring_matrices);%ALL electrodes

typeobject = {'phaseconn_matrices', 'powerconn_matrices', 'wiring_matrices'};
electrodelist = {'T','F','FP','IH','Grid', 'HD','NOHD', 'D', 'BiTemp'};
plotwiringcost(wiring_matrices, typeobject(3),electrodelist(9));
%plot (P*F*frequency)^2, distance*Hz is distance/s^2 energy to move a
%punctual mass between 2 electrodes.
plotkineticenergy(wiring_matrices);
%scatter PhysicalXfunctional, electrode to electrode
fprintf('Plotting the scatter of physical and functional distance\n');
scatterPF(wiring_matrices);
