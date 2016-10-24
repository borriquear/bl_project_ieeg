
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

%% 1. Power an dphase Analysis
% 1.1 Create mat file with power/amplitud calculated via the FFT. Power analysis, identify channels with most power and the frequency bands that pick up maximum power
%OUTPUT: % globalFsDir\eegpatient\data\figures\fft_%pat_%cond_%date_%se.mat -> 'ampli_fft','power_fft','power_fft_perband', 'power_fft_mean_perband','channel_labels'
global globalFsDir;
globalFsDir = loadglobalFsDir();
% patientslist = {'TWH030','TWH031','TWH033','TWH037','TWH038','TWH042','TWH043','TWH045','TWH047', 'TWH048','TWH049'};
% %patientslist = {'TWH047', 'TWH048'};
% conditionslist = {'EC_PRE', 'EO_PRE', 'HYP', 'EC_POST','EO_POST'};%,'EC_POST'};
% conditionslist = {'EC_PRE', 'EO_PRE'}%, 'HYP'};
patientslist = {'TWH030','TWH031','TWH033','TWH037','TWH038','TWH042','TWH043','TWH045','TWH047', 'TWH048','TWH049'};
conditionslist = {'EC_PRE', 'EO_PRE'}
temporalw = 5; % [mat file, and object with power spectra contained in mat file]
[powerx_perband, averagemeanpower] = waveletpowerspectraperperpatient(patientslist, conditionslist, temporalw);
%fprintf('DONE: mat file with power Patient %s, Condition %s in %s\n',eegpatient,eegcond, matfilename)
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
rois = roislist(9);
%ympb all (conds x5) % of power inthat band. condmeanrois for roi
[ympb,condmeanrois] = plotpowerspectrumallpatientsROI(patientslist, conditionslist, powerspecmatrix, powerfreqsindexes, powerspecmatrix_freqbands, rois);

% Statistical significance between conditions in power spectra.
% Calculate ttest to compare
plotstatisticalsignificance_powerspec(patientslist, conditionslist, powerspecmatrix, powerfreqsindexes, powerspecmatrix_freqbands, rois);
%% 1.3 Spectrogram
fprintf('Calling to plotspectrogramperpatient, Patient %s, Condition %s\n', eegpatient,eegcond)
plotspectrogram(patientslist, conditionslist, rois)

%% 2. NETWORK ANALYSIS: Create power_conn.mat and phase_conn.mat
% 2.1 Power adjacency matrices, create correlation matrix with Spearman based correlation 
%in power_conn.mat
%%%%%%%%%%%%%%%%%%%% POWER ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%% PHASE ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.2 Phase-based Analysis
% Calculate the Inter Site Phase Clustering (R = ||1/n \sum_t=1,n e^i(\phi_ch1,t - \phi_ch2,t)) and PLI.  
% Create correlation matrix with coherence and PLI based correlation in phaseconn_matrices.mat
% Get the  bivariate phase difference for every 2 channels append the result to the fft file
patientslist = {'TWH030','TWH031','TWH033','TWH037','TWH038','TWH042','TWH043','TWH045','TWH047', 'TWH048','TWH049'};
conditionslist = {'EC_PRE', 'EO_PRE', 'HYP'};%,'EC_POST'};
centerfrequencies  = logspace(log10(1),log10(50),8);
patientslist = {'TWH030'}; conditionslist = {'HYP'}; centerfrequencies = 1;
%temporalwindow in miliseconds
temporalwindow = 5;% temporalwindow: 5 argumentin phase_conn_tw.mat
calculatephasedifferences(patientslist, conditionslist, temporalwindow*1000);
matfilename = fullfile(globalFsDir, 'phaseconn_matrices.mat');
fprintf('Created file of phase connectivity in %s\n', matfilename);
%% 3. Display Networks calling displayconnectivity('po'|'pa') and calculatenetworkmetricdistances('po'|'pa')
%%%%%% POWER NETWORK DISPLAY %%%%%%%%%
% 3.1 Display Power-based correlation matrix and the undirected associated network 
%[srate, min_freq, max_freq, num_frex, time, n_wavelet, half_wavelet, freqs2use, s, wavelet_cycles]= initialize_wavelet();
%centerfrequenciesl  = logspace(log10(min_freq),log10(max_freq),8)
% centerfrequencies  = logspace(log10(1),log10(50),8);
% centerfrequencies  = centerfrequencies(3:end);
%display power connectivity (adj matrix and network ) from "power_conn.mat"
%displayconnectivity('power'|'phase')
%calculatenetworkmetricdistances('power'|'phase')
%INPUT: correlation matrix mat file
%OUTPUT: power_netmetrics.mat, phase_netmetrics.mat
displayconnectivity('power');
%INPUT: power_netmetrics.mat (metric distances between pair of conditions)
%OUTPUT:
calculatenetworkmetricdistances('power');
%%%%%% PHASE NETWORK DISPLAY %%%%%%%%%
%3.2 Display Phase-based correlation matrix
%display phase connectivity , adj matrix and network, from "phase_conn.mat"
displayconnectivity('phase');
%Creates metric distance file in D:\BIAL PROJECT\patients\phase_netmetrics.mat 
calculatenetworkmetricdistances('phase');
%% 4. Wiring distance (calls to calculateEuclideandistancematrix)
% 4.1. Create the wiring_matrices with the wiring cost for both power and
% phase based correlation matrices
% IN: phaseconn_matrices_tw.mat, 'powerconn_matrices_tw.mat (tw = when conn 
% matrix calculated using a moving window)
% OUTPUT: wiringcost_matrices.mat
wiring_matrices = calculatewiringcostmatrices();
%disp(wiring_matrices);
% 4.2 Plot the wiring cost matrices (ALL electrodes | subset)
plot_all_electrodes = 1;
if plot_all_electrodes == 1
    plotwiringcost(wiring_matrices); %ALL electrodes
else
    typeobject = {'phaseconn_matrices', 'powerconn_matrices', 'wiring_matrices'};
    electrodelist = {'T','F','FP','IH','Grid','HD','NOHD', 'D', 'BiTemp'};
    plotwiringcost(wiring_matrices, typeobject(3),electrodelist(6));
end
plot_scatter = 1;
if plot_scatter == 1
    fprintf('Plotting the scatter of physical and functional distance\n');
    scatterPF(wiring_matrices);
end

%% 5. TOPOLOGY  

%% 5.1 Initial Conditions 
%Calculates the object Badjall containing the binary matrix and the network
%metrics associated for that metric, for each patient, condition and
%frequency
%Load the wiring_matrices 
%clear all
% CREATE wiring_matrices  COMPLETE, select the subset of pats, conds and
% freqs for Badjall

% if exist('wiring_matrices') > 0
%     fprintf('wiring_matrices object already exist phaseconn and power_conn matrices\n')
% else
%     fprintf('Loading wiring matrices from phaseconn and power_conn matrices\n')
wiring_matrices = calculatewiringcostmatrices();
% end
%Initial conditions: select the patients, conditions and frequencies from wiring_matrices we
%want to setudy
initpat = 1; 
nbpats = length(wiring_matrices.patientslist);
listofselectedpats = wiring_matrices.patientslist(initpat:nbpats);
% wiring_matrices.frequencylist=
% 1.0000    1.7487    3.0579    5.3472    9.3506   16.3512   28.5930   50.0000
freqlist = wiring_matrices.frequencylist; 
initfreq = 3;%3 = 3.079, 5 = alpha
freq_gap = initfreq +5;%how many different frequencies (other than initfreq) to study (1 till theta, 2alpha)
freqlist = freqlist(initfreq:freq_gap);

nbfreqs = length(freqlist);
%condition EC EO HYP
initcond=1; 
nbconds = size(wiring_matrices.conditionslist, 2); 

connmatrixrequired = [1,0,0]; 
%get the wiring_ispc|pli|power oject from the wiring_matrices object
if find(connmatrixrequired) ==1
    wcconn = wiring_matrices.wiring_ispc;
    typeofconnectivity = 'coherence';
elseif find(connmatrixrequired) == 2
    wcconn = wiring_matrices.wiring_pli;
    typeofconnectivity = 'pli';
elseif find(connmatrixrequired) == 3
    wcconn = wiring_matrices.wiring_power;
    typeofconnectivity = 'power';
end
%Choose between building binary matrices or weighted matrices upon 
%all possible thresholds. If binary = 0 then Weighted
binary = 0;
%Set of Threshold MAtrices Binary Matrices depending on binary =1|0

%loop to calculate Badjall object containing the binary matrices and the
%network metrics associated with each BM
%% Calculates Badjall
nbpats = 1 + 1;% initpat + additional pats 
nbconds =3; % 2 ec-eo, 3 ec-hyp
nbfreqs = 6; % 1 till 3.0579, 2 till 5.3472 3 till 9.3506,  4- 16.3512,   5- 28.5930,  6 50.0000
Badjall = cell(nbpats,nbconds,nbfreqs,2 );
%array with all the initialconditions ie_patcondfreq initpat, nbpats, initcond, nbconds, initfreq, nbfreqs
ie_patcondfreq = [initpat,nbpats,initcond,nbconds,initfreq,freq_gap];
for pati=initpat:nbpats
    cupat = listofselectedpats(pati);
    patindex = pati-initpat + 1;
    for condj=1:nbconds
        cucond =  wiring_matrices.conditionslist(condj);
        for freqk=1:nbfreqs
            %Badjall{patindex,condj,freqk,:} = [];
            cufreq = freqlist(freqk);
            Badj = [];thresholdv = [];
            wcmatrix = []; wcmatrix_norm = [];
            fprintf('Calculating wcmatrix for pat=%s cond=%s freq=%s \n ',cupat{1}, cucond{1}, num2str(cufreq))
            wcmatrix= wcconn{initpat+pati-1,condj,freqk+initfreq-1};
            %normalize the matrix between 0 and 1
            if strcmp(typeofconnectivity, 'power') ==1
                %normalize between -1 and 1 (Power correlation is Spearman)
                wcmatrix_norm = -1 + 2.*(wcmatrix(:) - min(wcmatrix(:)))./(max(wcmatrix(:)) - min(wcmatrix(:)));
                wcmatrix_norm = reshape(wcmatrix_norm, size(wcmatrix,1),size(wcmatrix,2));
            else
                %normalize between 0 and 1 (coherence and PLI are [0,1])
                wcmatrix_norm = (wcmatrix - min(wcmatrix(:)))/(max(wcmatrix(:)) - min(wcmatrix(:)));
            end
            [Badj, thresholdv] = calc_threshold_wmatrix(wcmatrix_norm, binary);
            fprintf('Saving Badj and thresholdv in Badjall \n');
            
            Badjall{patindex,condj,freqk,1} = Badj;
            Badjall{patindex,condj,freqk,2} = thresholdv;
        end
    end
end

%% 5.2
fprintf('Badjall built, click OK to plot the network metric for each delta\n');
f = warndlg('Badjall built, Click OK to Continue with the Ploting ', 'Program interruption');
drawnow 
waitfor(f);
disp('Plotting Badjall is an array Badj for eachpat, cond anf freq: Badj{1,:}Binary matrices, one for each threshold and Badj{2,:} network metrics, one for each threshold');
%ie_patcondfreq = [initpat,nbpats,initcond,2,initfreq,freq_gap];
[netmetlist] = plotnetworkmetrics_all(Badjall, wiring_matrices, ie_patcondfreq, binary, typeofconnectivity)
%% 5.3 Calculates the statistical effect between the conditions for the binary 
% networks
% H0= \mu condition 1 = \mu condition 2
%ttest for eachpatiend and frequency between two conditions
%cond1 =1; cond2=2; %cond3=3 %hypnois
fprintf('Calculating the statistical significance for the two conditions\n');
f = warndlg('Calculating the statistical significance for the two conditions. Click OK', 'Program interruption');
drawnow 
waitfor(f);
%netmetlist = {Betti_v, clustering_v,density_v,pathlength_v};
alphav = 0.001;
for i=1:nbpats
    cupat = listofselectedpats(pati + initpat - 1);
    for f=1:nfreqs
        vcond_a = {};vcond_a_betti = {}; vcond_a_clustering ={};vcond_a_path ={};
        vcond_b = {};vcond_b_betti = {};vcond_a_clustering ={};vcond_a_path ={};
        vcond_a = netmetlist{i,cond1,f};
        vcond_b = netmetlist{i,cond2,f};
        vcond_a_betti = vcond_a{1}; vcond_b_betti = vcond_b{1};
        [h,p] = ttest2(vcond_a_betti,vcond_b_betti, 'Alpha',alphav);
        fprintf('Betti ttest for Pat= %s, freq=%.2f: h = %d p=%.5f alpha =%.4f\n',listofselectedpats{i + initpat - 1},freqlist(f),h,p, alphav);
        vcond_a_clustering = vcond_a{1}; vcond_b_clustering = vcond_b{1};
        [h,p] = ttest2(vcond_a_clustering,vcond_b_clustering, 'Alpha',alphav);
        fprintf('Clustering ttest for Pat= %s, freq=%.2f: h = %d p=%.5f alpha =%.4f\n',listofselectedpats{i + initpat - 1},freqlist(f),h,p, alphav);
        vcond_a_path = vcond_a{1}; vcond_b_path = vcond_b{1};
        [h,p] = ttest2(vcond_a_clustering,vcond_b_clustering, 'Alpha',alphav);
        fprintf('Path ttest for Pat= %s, freq=%.2f: h = %d p=%.5f alpha =%.4f \n',listofselectedpats{i + initpat - 1},freqlist(f),h,p, alphav);
    end
end
% The result h is 1 if the test rejects the null hypothesis (data in vectors x and y comes from independent random samples from normal distributions with equal means and equal but unknown variances) at the 5% significance level, and 0 otherwise.
% Save the ttest results in a mat file















%%%%%% WORK ON DECONSTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Episurg to plot bonitos cerebros
%Load condition [power|other]data and display coloured electrodes
%clear all
patientid = 'TWH030';
patientcond = 'HYP';
nbpatient = 1; %TWH030
[myfullname, EEG, channel_labels, patdate, patsession] =  initialize_EEG_variables(patientid,patientcond)
cfg=[];
cfg.view='r';
cfg.ignoreDepthElec='n';
cfg.opaqueness=.1;
cfg.figId=1;
%cfg.overlayParcellation='DK';
cfg.showLabels='y';
cfg.title=[];
cfgOut=plotPialSurf(patientid,cfg);

%% Plot bars between electrodes, color coded to reflect bipolar reference correlation value
connmat = wiring_matrices.wiring_ispc{nbpatient};
nbelecs = length(connmat);
pairs=[];
totpairs = (nbelecs^2 - nbelecs)/2;
 ct=0;
 eleclabels = {};
for a=1:8,
    ct=ct+1;
    idxch = a+1;
    pairs{ct,1}= channel_labels{idxch};
   	pairs{ct,2}= channel_labels{idxch+1};
    pairs{ct,3}=rand(1,3); % RGB val
    pairs{ct,4}='l';
    pairs{ct,5}=channel_labels{idxch};
    pairs{ct,6}='lineWidth'  ;
end

cfg.view='l';
cfg.figId=2;
cfg.pairs=pairs;
cfg.showLabels='y';
cfg.elecUnits='l';
cfg.title='TWH030: Stimulus Correlations';
cfg_out=plotPialSurf('TWH037',cfg);
 
 %% Kinetic Energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot data points in the K axis and the percentage of power relative to
% other frequency bands. Need to have calculated wiring_matrices and
% powerspecmatrix_freqbands(fft section 1)
plotkineticenergy(wiring_matrices, powerspecmatrix_freqbands);
%scatter PhysicalXfunctional, electrode to electrode