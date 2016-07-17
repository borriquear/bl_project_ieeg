function [frex,frexidx1,frexidx2, powerx, powerx_mean_perbands] = loadpowerspectrumperpatient(eegpatient, eegcondition)
%plotpowerspectrum displays power spectrum for a fft mat file. the power
%spectrum is already being generated, the function ONLY displays results
% IN eegpatient, eegcondition
% OUT frex, powerx: vector of frequencies i'm interested in eg 0:50 and
% power values in that frequency domain
% REQUIRES fft_*.matfiles, created at fourieranalysis.mat
% For a given signal, the power spectrum gives a plot of the portion of a
% signal's power (energy per unit time) falling within given frequency bins.

global globalFsDir;
[globalFsDir] = loadglobalFsDir();
nyquist = 500; % for 1000 sampling rate
fprintf('Obtaining the matfile for patient %s and condition %s\n',eegpatient,eegcondition );
[matfile] = getfftmatfilefrompatcond(eegpatient,eegcondition);
%load mat file and save the vector
%[myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(eegpatient,eegcond);
patdir = fullfile(globalFsDir, eegpatient, 'data\figures');
matfile = fullfile(patdir,matfile);
fprintf('Opening mat file : %s \n',matfile);
fhd = load(matfile); %'power_fft','power_fft_perband', 'channel_labels'
%from 0 to nyquist, size(power_fft,2) linearly spaced points(Hzs)
frequencies = linspace(0,nyquist,size(fhd.power_fft,2));
%to plot only significative frequencies, let us say between 0 and 50
freqiwant1=0;freqiwant2=50;
[~,frexidx1] = min(abs(frequencies-freqiwant1));
[~,frexidx2] = min(abs(frequencies-freqiwant2));
%plot(frequencies(frexidx1:frexidx2),mean(power_fft(:,(frexidx1:frexidx2))))
frex = frequencies(frexidx1:frexidx2);
powerx = fhd.power_fft(:,(frexidx1:frexidx2));
powerx_mean_perbands = fhd.power_fft_mean_perband;
end


function [pre] = getfftmatfilefrompatcond(patientid, pre_label)
% getfftmatfilefrompatcond return the mat file that contains the fft results
% for power
if strcmp(patientid, 'TWH043') == 1
    if strcmp(pre_label, 'EC_PRE') == 1
        pre = 'fft_EC_PRE_TWH043_05042016_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') == 1
        pre = 'fft_EO_PRE_TWH043_05042016_s1.mat';
    elseif strcmp(pre_label, 'HYP') == 1
        pre = 'fft_HYP_TWH043_05042016_s1.mat';
    elseif strcmp(pre_label, 'EC_POST') == 1
        pre = 'fft_EC_POST_TWH043_05042016_s1.mat';
    elseif strcmp(pre_label, 'EO_POST') == 1
        pre = 'fft_EO_POST_TWH043_05042016_s1.mat';
    end
elseif strcmp(patientid, 'TWH042') ==1
    if strcmp(pre_label, 'EC_PRE')  == 1
        pre = 'fft_EC_PRE_TWH042_05042016_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') == 1
        pre = 'fft_EO_PRE_TWH042_05042016_s1.mat';
    elseif strcmp(pre_label, 'HYP') == 1
        pre = 'fft_HYP_TWH042_05042016_s1.mat';
    elseif  strcmp(pre_label, 'EC_POST') == 1
        pre = 'fft_EC_POST_TWH042_05042016_s1.mat';
    elseif  strcmp(pre_label, 'EO_POST') == 1
        pre = 'fft_EO_POST_TWH042_05042016_s1.mat';
    end
elseif strcmp(patientid, 'TWH038') ==1
    if strcmp(pre_label, 'EC_PRE') == 1
        pre = 'fft_EC_PRE_TWH038_03082016_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') == 1
        pre = 'fft_EO_PRE_TWH038_03082016_s1.mat';
    elseif strcmp(pre_label, 'HYP') == 1
        pre = 'fft_HYP_TWH038_03082016_s1.mat';
    elseif strcmp(pre_label, 'EC_POST') == 1
        pre = 'fft_EC_POST_TWH038_03082016_s1.mat';
    elseif strcmp(pre_label, 'EO_POST') == 1
        pre = 'fft_EO_POST_TWH038_03082016_s1.mat';
    end
elseif strcmp(patientid, 'TWH037') ==1
    if strcmp(pre_label, 'EC_PRE') == 1
        pre = 'fft_EC_PRE_TWH037_03142016_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') == 1
        pre = 'fft_EO_PRE_TWH037_03142016_s1.mat';
    elseif strcmp(pre_label, 'HYP') == 1
        pre = 'fft_HYP_TWH037_03142016_s1.mat';
    elseif  strcmp(pre_label, 'EC_POST') == 1
        pre = 'fft_EC_POST_TWH037_03142016_s1.mat';
    elseif  strcmp(pre_label, 'EO_POST') == 1
        pre = 'fft_EO_POST_TWH037_03142016_s1.mat';
    end
elseif strcmp(patientid, 'TWH034') ==1
    if strcmp(pre_label, 'EC_PRE') == 1
        pre = 'fft_EC_PRE_TWH034_02092016_s2.mat';
    elseif  strcmp(pre_label, 'HYP') == 1
        pre  = 'fft_HYP_TWH034_02092016_s2.mat';
    elseif  strcmp(pre_label, 'EC_POST') == 1
        pre  = 'fft_EC_POST_TWH034_02092016_s2.mat';
    elseif  strcmp(pre_label, 'EO_POST') == 1
        pre  = 'fft_EO_POST_TWH034_02092016_s2.mat';
    end
elseif strcmp(patientid, 'TWH033') ==1
    if strcmp(pre_label, 'EC_PRE') == 1
        pre = 'fft_EC_PRE_TWH033_02032016_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') == 1
        pre = 'fft_EO_PRE_TWH033_02032016_s1.mat';
    elseif strcmp(pre_label, 'HYP') == 1
        pre = 'fft_HYP_TWH033_02032016_s1.mat';
    elseif strcmp(pre_label, 'EC_POST') == 1
        pre = 'fft_EC_POST_TWH033_02032016_s1.mat';
    elseif strcmp(pre_label, 'EO_POST') == 1
        pre = 'fft_EO_POST_TWH033_02032016_s1.mat';
    end
elseif strcmp(patientid, 'TWH031') ==1
    if strcmp(pre_label, 'EC_PRE') == 1
        pre = 'fft_EC_PRE_TWH031_12012015_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') == 1
        pre = 'fft_EO_PRE_TWH031_12012015_s1.mat';
    elseif strcmp(pre_label, 'HYP') == 1
        pre = 'fft_HYP_TWH031_12012015_s1.mat';
    elseif strcmp(pre_label, 'EC_POST') == 1
        pre = 'fft_EC_POST_TWH031_12012015_s1.mat';
    end
elseif strcmp(patientid, 'TWH030') ==1
    if strcmp(pre_label, 'EC_PRE') == 1
        pre = 'fft_EC_PRE_TWH030_11172015_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') == 1
        pre = 'fft_EO_PRE_TWH030_11172015_s1.mat';
    elseif  strcmp(pre_label, 'HYP') == 1
        pre = 'fft_HYP_TWH030_11172015_s1.mat';
    elseif strcmp(pre_label, 'EC_POST') == 1
        pre ='fft_EC_POST_TWH030_11172015_s1.mat';
    end
end
end
%>>print -f1 -djpeg twh045condA % make a jpg of Figure 1
%>>savefig('twh045condA') % save current figure to disk as a MATLAB fig file
