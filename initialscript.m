
%% Script that goes step by step processing the data
%0. The EEG objects for each patients need to be already created (EEGLab)
% cuteoneepochNEW.m, read overleaf documentation           
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
%>>print -f3 -djpeg '/Users/jaime/BIALPROJECT/patients/MCP-local-PLI-nbedges'
%>>savefig('D:\BIAL PROJECT\patients\figure_results\ECEO-Allroisfigi1')

%% 1. Power and phase Analysis  

%% 1.1 Display the MEAN power spectra for the patients for the entire brain
% Power analysis, identify channels with most power and the frequency bands that pick up maximum power
% Create mat file with power/amplitud calculated via the FFT. 
% OUTPUT: % globalFsDir\eegpatient\data\figures\fft_%pat_%cond_%date_%se.mat -> 'ampli_fft','power_fft','power_fft_perband', 'power_fft_mean_perband','channel_labels'
global globalFsDir;
globalFsDir = loadglobalFsDir();
patientslist = {'TWH030','TWH031','TWH033','TWH037','TWH038','TWH042','TWH043','TWH045','TWH047', 'TWH048','TWH049'};
conditionslist = {'EC_PRE','EO_PRE'}
temporalw = 5; % [mat file, and object with power spectra contained in mat file]
[powerx_perband, averagemeanpower] = waveletpowerspectraperperpatient(patientslist, conditionslist);

%%  1.2. Display the power spectra for the patients per ROI
% Requirement: Needs the fft_%pat_%cond_%date_%se.mat files created in 1.1
% from that will load the objects that contain the amplitude 
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
rois = roislist(1);
%ympb all (conds x5) % of power in that band. condmeanrois for roi
[ympb,condmeanrois] = plotpowerspectrumallpatientsROI(patientslist, conditionslist, powerspecmatrix, powerfreqsindexes, powerspecmatrix_freqbands, rois);

% ttest Statistical significance between conditions in power spectra.
plotstatisticalsignificance_powerspec(patientslist, conditionslist, powerspecmatrix, powerfreqsindexes, powerspecmatrix_freqbands, rois);
%% 1.3 Spectrogram
%plot the spectogram (time x frequency), 2 charts per patients, entire time and some short window in the middle 
fprintf('Calling to plotspectrogramperpatient, Patient %s, Condition %s\n', eegpatient,eegcond)
patientslist = {'TWH033'}
plotspectrogram(patientslist, conditionslist, rois)

%% 2. NETWORK ANALYSIS: Create power_conn.mat and phase_conn.mat
% 2.1 Power adjacency matrices, create correlation matrix with Spearman based correlation 
%%%%%%%%%%%%%%%%%%%% POWER ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
patientslist = {'TWH030','TWH031','TWH033','TWH037','TWH038','TWH042','TWH043','TWH045','TWH047', 'TWH048','TWH049'};
conditionslist = {'EC_PRE', 'EO_PRE'} %, 'HYP','EC_POST'};

centerfrequencies  = logspace(log10(1),log10(50),8);
%create powerconn matrices, first for segments of temporalwindow secs in powerconn_matrices.mat
%1. Create powerconn_matrices_tw.mat for temporalwindow >0 (cut the signal into segments)
temporalwindow = 5; % 4 seconds powerconn_matrices_tw.mat
powerconn_matrix = createpowerbcorrmatrix(patientslist, conditionslist, centerfrequencies, temporalwindow);
powerconnmatf = fullfile(globalFsDir, 'powerconn_matrices_tw.mat');

%temporalwindow: 0 means correlation calculated with entire epoch length
% generate powerconn_matrices.mat for no temporal window
temporalwindow = 0;
if temporalwindow > 0 %skip this if dont want to generate power_connfor tw=0
    fprintf('Created file of power connectivity in %s\n', powerconnmatf);
    %2.Create powerconn_matrices.mat now for for temporalwindow = 0 , entire epoch
    createpowerbcorrmatrix(patientslist, conditionslist, centerfrequencies);
    powerconnmatf = fullfile(globalFsDir, 'powerconn_matrices.mat');
    fprintf('Created file of power connectivity in %s\n', powerconnmatf);
end
%%%%%%%%%%%%%%%%%%%% PHASE ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% 2.2 Phase-based Analysis
% Calculate the Inter Site Phase Clustering (R = ||1/n \sum_t=1,n e^i(\phi_ch1,t - \phi_ch2,t)) and PLI.  
% Create correlation matrix with coherence and PLI based correlation in phaseconn_matrices.mat
% Get the  bivariate phase difference for every 2 channels append the result to the fft file
patientslist = {'TWH030','TWH031','TWH033','TWH037','TWH038','TWH042','TWH043','TWH045','TWH047', 'TWH048','TWH049'};
conditionslist = {'EC_PRE', 'EO_PRE', 'HYP'};%,'EC_POST'};
centerfrequencies  = logspace(log10(1),log10(50),8);
%patientslist = {'TWH030'}; conditionslist = {'HYP'}; centerfrequencies = 1;
%temporalwindow in miliseconds
temporalwindow = 5;% temporalwindow: 5 argumentin phase_conn_tw.mat
phaseconn_matrix = calculatephasedifferences(patientslist, conditionslist, temporalwindow*1000);
matfilename = fullfile(globalFsDir, 'phaseconn_matrices.mat');
fprintf('Created file of phase connectivity in %s\n', matfilename);
%% 3. Display correlation matrix and network calling displayconnectivity('po'|'pa') and calculatenetworkmetricdistances('po'|'pa')
% 3.1 Display Power|Phase-based correlation matrix and the undirected associated network 
%display power|phase connectivity (adj matrix and network ) from "power|phase conn_matrices_tw.mat"
type_analysis= {'power','phase'} %displayconnectivity('power'|'phase')
for i =1:length(type_analysis)
    fprintf('Displaying connectivity for %s\n', type_analysis{i});
    %INPUT: correlation matrix: power|phase conn_matrices_tw
    %OUTPUT: power_netmetrics.mat, phase_netmetrics.mat
    displayconnectivity(type_analysis(i))
    fprintf('Calcultating network metrics for %s\n', type_analysis{i});
    %INPUT: power_netmetrics.mat (metric distances between pair of conditions)
    %OUTPUT:
    calculatenetworkmetricdistances(type_analysis(i));
end    

%% 4. Wiring distance (calls to calculateEuclideandistancematrix)
% 4.1. Create the wiring_matrices with the wiring cost for both power and
% phase based correlation matrices
% IN: scope, phaseconn_matrices_tw.mat, 'powerconn_matrices_tw.mat (tw = when conn 
% matrix calculated using a moving window)
% OUTPUT: wiringcost_matrices.mat

globalFsDir = loadglobalFsDir();
scope = {'local','meso'}; %W = D.*F, W= D./LH
for i =1:length(scope)
    fprintf('Calculating wiring cost for for %s\n', scope{i})
    %INPUT:power|phase_conn_matrices_tw.mat
    %OUT: The Wiring Matrices D*F
    wiring_matrices = calculatewiringcostmatrices(scope{i}); 
    if i == 1 
        wiring_matrices_local = wiring_matrices;
    end
    if i == 2
        wiring_matrices_meso = wiring_matrices;
    end
    fprintf('Created file wiringcost_matrices.mat_%s', scope{i});
end

%% 5. ttest for power, phase connectivity mattix and for wiring_local and wiringcost_messo

%ttest and report for power based connectivity
filename='/Users/jaime/BIALPROJECT/patients/powerconn_matrices_tw.mat';
powhdl = load(filename);
performttestandwilcox_connectivity(powhdl, 'power_conn');
%ttest and report for phase based connectivity
filename='/Users/jaime/BIALPROJECT/patients/phaseconn_matrices_tw.mat';
phasehdl = load(filename);
performttestandwilcox_connectivity(phasehdl, 'pli');
performttestandwilcox_connectivity(phasehdl, 'i');
%ttest for wiring cost local and mesos
performttestandwilcox(wiring_matrices_local, 'wc_local');
performttestandwilcox(wiring_matrices_meso, 'wc_meso');

%%
%Plot wiring cost and do ttest and Man whithney brain all
plotwiringcost(wiring_matrices_local);
plotwiringcost(wiring_matrices_meso);
%Plot wiring cost and do ttest and Man whithney for regions
electrodelist = {'T','F','FP','IH','Grid','HD','NOHD', 'D', 'BiTemp'};
typeobject = {'phaseconn_matrices', 'powerconn_matrices', 'wiring_matrices'};
for i=1:length(electrodelist)
    plotwiringcost(wiring_matrices_local, typeobject{3},electrodelist(i));
    savefiguresloop(i);
    plotwiringcost(wiring_matrices_meso, typeobject{3},electrodelist(i));
    savefiguresloop(i);
end
% Plot Differences in connectivity matrices

plotwiringcost(typeobject{1});
plotwiringcost(typeobject{2});

% Plot Differences in wiring cost matrices
%Plot per ROIs
typeobject = {'phaseconn_matrices', 'powerconn_matrices', 'wiring_matrices'};
electrodelist = {'T','F','FP','IH','Grid','HD','NOHD', 'D', 'BiTemp'};
for i=1:length(electrodelist)
    plotwiringcost(wiring_matrices, typeobject(3),electrodelist(i));
    savefiguresloop(i);
end


% 4.2 Plot the wiring cost matrices (ALL electrodes | subset)
plot_all_electrodes = 0;
if plot_all_electrodes == 1
    plotwiringcost(wiring_matrices); %ALL electrodes
    savefiguresloop(0);
else
    typeobject = {'phaseconn_matrices', 'powerconn_matrices', 'wiring_matrices'};
    electrodelist = {'T','F','FP','IH','Grid','HD','NOHD', 'D', 'BiTemp'};
    for i=1:length(electrodelist)
        plotwiringcost(wiring_matrices, typeobject(3),electrodelist(i));
        savefiguresloop(i);
    end
end

%plot_scatter = 0;
%if plot_scatter == 1
%    fprintf('Plotting the scatter of physical and functional distance\n');
%    scatterPF(wiring_matrices);
%end

%% 5. TOPOLOGY  

%% 5.1 Initial Conditions 
%Calculates the object Badjall containing the binary matrix and the network
%metrics associated for that metric, for each patient, condition and
%frequency
%Initial conditions: select the patients, conditions and frequencies from wiring_matrices we
%want to setudy
initpat = 1; 
nbpats = length(wiring_matrices.patientslist);
listofselectedpats = wiring_matrices.patientslist(initpat:nbpats);
% wiring_matrices.frequencylist=
% 1.0000    1.7487    3.0579    5.3472    9.3506   16.3512   28.5930   50.0000
freqlist = wiring_matrices.frequencylist; 
initfreq = 3; %3 = 3.079, 5 = alpha
freq_gap = initfreq + 5; %how many different frequencies (other than initfreq) to study (1 till theta, 2alpha)
freqlist = freqlist(initfreq:freq_gap);

nbfreqs = length(freqlist);
%condition EC EO HYP
initcond=1; 
nbconds = size(wiring_matrices.conditionslist, 2); 

connmatrixrequired = [0,1,0]; %100 010 001 
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
%% Calculates Badjall
%Choose between building binary matrices or weighted matrices upon 
%all possible thresholds. If binary = 0 then Weighted
binary = 1;
%Set of Threshold MAtrices Binary Matrices depending on binary =1|0
nbpats = 1 + 10;% initpat + additional pats 
nbconds = 2; % 2 ec-eo, 3 ec-hyp
%nbfreqs = 3; % 1 till 3.0579, 2 till 5.3472 3 till 9.3506,  4- 16.3512,   5- 28.5930,  6 50.0000
%Bdjall contains all networks but for only one measure (R, PLI or Power)
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
            [Badj,thresholdv] = calc_threshold_wmatrix(wcmatrix_norm, binary);
            fprintf('Saving Badj and thresholdv in Badjall \n');
            Badjall{patindex,condj,freqk,1} = Badj;
            Badjall{patindex,condj,freqk,2} = thresholdv;
        end
    end
end
fprintf('Done! Badjall is being created \n');
fna = sprintf('Bdjall_local_%s.mat',typeofconnectivity);
fprintf('Saving Bdjall in\n');
%save(fna, 'Badjall', '-v7.3')
%% t-test of B0, clustering, wc and path length are differen t for the two conditions with Kruskal Wallis and MCP
%frequency fixed inside
[ttest_resuls] = ttest_networkmetrics(Badjall);
%ttest_resuls = {[Bh,Bp,Bci,Bstats],[Ch,Cp,Cci,Cstats],[Wh,Wp,Wci,Wstats] ,[Ph,Pp,Pci,Pstats]  }


%% 5.2.1 Plot Badjall matrix results (overall metrics for all networks)

disp('Plotting Badjall is an array Badj for each pat, cond and freq: Badj{1,:}Binary matrices, one for each threshold and Badj{2,:} network metrics, one for each threshold');
%ie_patcondfreq = [initpat,nbpats,initcond,2,initfreq,freq_gap];
[netmetlist] = plotnetworkmetrics_all(Badjall, wiring_matrices, ie_patcondfreq, binary, typeofconnectivity)
%[netmetlist] = plotnetworkmetrics_all(Badjall, wiring_matrices, ie_patcondfreq, binary, 'power')
%% 5.2.2 Plot selected matrices (overall metrics for all networks)
if exist('Badjall') > 0
    fprintf('wiring_matrices object already exist phaseconn and power_conn matrices\n');
    %[nbpats, nbconds, nbfreqs,2] =size(Badjall);
else
    warning('Badjall NOT FOUND\n')
    return
end

figure
mypat = 6; mypatlabel = 'TWH042';
mycond = 1; %EC EO HYP
myfreq = 3;  % 1 is 3.05
myconnmat = Badjall(mypat,mycond,myfreq,1);
mythreshold_v = Badjall(mypat,mycond,myfreq,2);
mythreshold_v = mythreshold_v{1};
mycurrmat = myconnmat{1};
[myfullname, EEG, channel_labels, eegdate, eegsession] = initialize_EEG_variables(mypatlabel,'HYP');
strNames = channel_labels(2:end);

subplot(1,3,1)
threscut = round((length(mythreshold_v) -length(mythreshold_v)+500 ));
corrMatrix = mycurrmat{mycond,threscut};
myColorMap = lines(length(corrMatrix));
circularGraph(corrMatrix,'Colormap',myColorMap,'Label',strNames);
subplot(1,3,2)
threscut = round((length(mythreshold_v) -length(mythreshold_v)+800 ));
corrMatrix = mycurrmat{mycond,threscut};
myColorMap = lines(length(corrMatrix));
circularGraph(corrMatrix,'Colormap',myColorMap,'Label',strNames);
subplot(1,3,3)
threscut = round((length(mythreshold_v) - 1 ));
corrMatrix = mycurrmat{mycond,threscut};
myColorMap = lines(length(corrMatrix));
circularGraph(corrMatrix,'Colormap',myColorMap,'Label',strNames);

msgtitle = sprintf('Wiring cost network TWH042 alpha band threshold=%d , pat=%s', mypatlabel,threscut )
title(msgtitle)

%% 5.3 Calculates the statistical effect between the conditions for the binary 
% networks
% H0= \mu condition 1 = \mu condition 2
%ttest for eachpatiend and frequency between two conditions
%cond1 =1; cond2=2; %cond3=3 %hypnosis
fprintf('Calculating the statistical significance for the two conditions\n');
f = warndlg('Calculating the statistical significance for the two conditions. Click OK', 'Program interruption');
drawnow 
waitfor(f);
%netmetlist = {Betti_v, clustering_v,density_v,pathlength_v};
alphav = 0.01;
cond1 = 1; cond2 = 2;
p_metrics = cell(nbpats,nbfreqs,4);
for i=1:nbpats
    cupat = listofselectedpats(pati + initpat - 1);
    for f=1:nbfreqs
        vcond_a = {};vcond_a_betti = {}; vcond_a_clustering ={};vcond_a_path ={};
        vcond_b = {};vcond_b_betti = {};vcond_b_clustering ={};vcond_b_path ={};
        vcond_a = netmetlist{i,cond1,f};
        vcond_b = netmetlist{i,cond2,f};
        vcond_a_betti = vcond_a{1}; vcond_b_betti = vcond_b{1};
        [h,p] = ttest2(vcond_a_betti,vcond_b_betti, 'Alpha',alphav);
        p_metrics{i,f,1} = p;
        fprintf('Betti ttest for Pat= %s, freq=%.2f: h = %d p=%.5f alpha =%.4f\n',listofselectedpats{i + initpat - 1},freqlist(f),h,p, alphav);
        vcond_a_clustering = vcond_a{2}; vcond_b_clustering = vcond_b{2};
        [h,p] = ttest2(vcond_a_clustering,vcond_b_clustering, 'Alpha',alphav);
        p_metrics{i,f,2} = p;
        fprintf('Clustering ttest for Pat= %s, freq=%.2f: h = %d p=%.5f alpha =%.4f\n',listofselectedpats{i + initpat - 1},freqlist(f),h,p, alphav);
        vcond_a_wiring = vcond_a{3}; vcond_b_wiring = vcond_b{3};
        [h,p] = ttest2(vcond_a_wiring, vcond_b_wiring, 'Alpha',alphav);
        p_metrics{i,f,3} = p;
        fprintf('Nb of Edges ttest for Pat= %s, freq=%.2f: h = %d p=%.5f alpha =%.4f \n',listofselectedpats{i + initpat - 1},freqlist(f),h,p, alphav);
        vcond_a_path = vcond_a{4}; vcond_b_path = vcond_b{4};
        [h,p] = ttest2(vcond_a_path, vcond_b_path, 'Alpha',alphav);
        p_metrics{i,f,4} = p;
        fprintf('Path ttest for Pat= %s, freq=%.2f: h = %d p=%.5f alpha =%.4f \n',listofselectedpats{i + initpat - 1},freqlist(f),h,p, alphav);
        fprintf('FREQ =%.2f\n', f)
    end
    fprintf('=======PATIENT======= =%d\n\n', i)
end
% The result h is 1 if the test rejects the null hypothesis (data in vectors 
% x and y comes from independent random samples from normal distributions with equal means and equal but unknown variances) at the 1% significance level, and 0 otherwise.
% Save the ttest results in a mat file
folderderst = 'D:\BIALPROJECT\patients\figure_results\';
fileorig = sprintf('p-test-01-BCP');
filedest =fullfile(folderderst,fileorig );
save('filedest','p_metrics');
%% Mutual Information
patientlist = {'TWH030'}%,'TWH031','TWH034','TWH033','TWH038','TWH042'};%,'TWH043','TWH037'
%Patient 37 and 43  have NOT HD,  only Deep
conditionslist = {'EC_PRE', 'EO_PRE'}%,'HYP','EC_POST','EO_POST'};
[mi_matrix] = mutualinformation_matrix(patientlist, conditionslist)









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



