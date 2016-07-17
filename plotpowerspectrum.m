function [ ] = plotpowerspectrumperpatient(eegpatient, eegcondition)
%plotpowerspectrum displays power spectrum for a fft mat file. the power
%spectrum is already being generated, the function ONLY displays results
% IN []
% OUT []
% REQUIRES fft_*.matfiles, created at fourieranalysis.mat
% For a given signal, the power spectrum gives a plot of the portion of a
% signal's power (energy per unit time) falling within given frequency bins.
% The most common way of generating a power spectrum is by using a discrete Fourier transform, but other techniques such as the maximum entropy method can also be used.
global globalFsDir;
[globalFsDir] = loadglobalFsDir();
% patientlist = {'TWH030','TWH031','TWH034','TWH033','TWH037','TWH038','TWH042','TWH043'};%,'TWH043','TWH037'
%Patient 37 and 43  have NOT HD,  only Deep
% conditionslist = {'EC_PRE', 'HYP','EC_POST',}; %'EO_PRE' 'EO_POST'
nyquist = 500; % for 1000 sampling rate
powerspecmatrix = {};
powerfreqsindexes = {};

        fprintf('Obtaining the matfile for patient %s and condition %s\n',eegpatient,eegcondition );
        [matfile] = getfftmatfilefrompatcond(eegpatient,eegcondition);
        %load mat file and save the vector
        %[myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(eegpatient,eegcond);
        patdir = fullfile(globalFsDir,eegpatient, 'data\figures');
        matfile = fullfile(patdir,matfile);
        fprintf('Opening mat file : %s \n',matfile);
        fhd = load(matfile); %channel_labels, ampli_fft, power_fft,frqperband,percentlistoffrqperband
        frequencies = linspace(0,nyquist,size(fhd.power_fft,2)); %from 0 to nyquist, size(power_fft,2) linearly spaced points(Hzs)
        %to plotonly significative frequencies, let us say between 0 and 50
        freqiwant1=0;freqiwant2=50;
        [~,frexidx1] = min(abs(frequencies-freqiwant1));
        [~,frexidx2] = min(abs(frequencies-freqiwant2));
        %plot(frequencies(frexidx1:frexidx2),mean(power_fft(:,(frexidx1:frexidx2))))
        frex = frequencies(frexidx1:frexidx2);
        powerx = fhd.power_fft(:,(frexidx1:frexidx2));
        powerfreqsindexes{ip,ic} = frex;
        powerspecmatrix{ip,ic} = powerx;
        xlimtodisplay = 30;
        % plot the powerspecmatrix per patient

            fighdl(ip,ic,1) = figure;
            maxvalue = max(powerx(:));
            meanvx = mean(powerx);
            stdx= std(powerx);
            maxofmean =max(meanvx);
            % to plot the the mean values, one figure per patient
            %plot(frequencies(frexidx1:frexidx2),power_fft(:,(frexidx1:frexidx2)));
            %hold on
            plot(frequencies(frexidx1:frexidx2),meanvx); hold on;
            plot(frequencies(frexidx1:frexidx2),stdx);
            set(gca, 'YLim', [0, 1.1*maxofmean], 'XLim', [0 xlimtodisplay]);
            ylabel('Power')
            msgt = sprintf('Power spectra from FFT %s %s',eegpatient, eegcondition);
            title(msgt);
            legend('mean', 'std');
            % uncomment this to plot all channels, one figure per patient
            %             fighdl(ip,ic,2) = figure;
            %             plot(frequencies(frexidx1:frexidx2),powerx);
            %             set(gca, 'YLim', [0, 1.1*maxvalue], 'XLim', [0 xlimtodisplay]);
            %             ylabel('Power')
            %             msgt = sprintf('Power spectra from FFT %s %s',eegpatient, eegcondition);
            %             title(msgt);
            %             %legend('Power all channels');
            %             legend(fhd.channel_labels(2:end));

%plot all patients mean of all channels in one figure
%figure1=figure('Position', [100, 100, 1024, 1200]);
figm = figure;
suptitle(['Mean Power spectra All areas']);
figstd = figure;
suptitle(['Mean and std Power spectra All areas']);
nbp = length(patientlist);
nbc = length(conditionslist);
rois = {'HD','T','F', 'IH'};
for ip=1:length(patientlist)
    eegpatient = patientlist{ip};
    for ic=1:length(conditionslist)
        eegcondition = conditionslist{ic};
        %ha = tight_subplot(nbp,nbc,3*(ip-1) + ic,[.01 .01])
        powerx = powerspecmatrix{ip,ic};
        meanvx = mean(powerx);
        stdx = std(powerx);
        maxofmean =max(meanvx);
        freqsindexl = powerfreqsindexes{ip,ic};
        figure(figm);
        subplot(nbp,nbc,3*(ip-1) + ic);
        % plot mean as such
        plot(freqsindexl,meanvx);
        set(gca, 'YLim', [0, 1.1*maxofmean], 'XLim', [0 xlimtodisplay]);
        ylabel('Power ')
        msgt = sprintf('Mean Power spectra from FFT %s %s', eegpatient, eegcondition);
        title(msgt);
        % plot mean and error bars
        figure(figstd);
        subplot(nbp,nbc,3*(ip-1) + ic);
        errorbar(freqsindexl,meanvx,stdx);
        maxofmeanandstd = max(maxofmean, max(stdx));
        set(gca, 'YLim', [0, 1.1*maxofmeanandstd], 'XLim', [0 xlimtodisplay]);
        ylabel('Power')
        msgt = sprintf('Mean and Std Power spectra from FFT %s %s', eegpatient, eegcondition);
        title(msgt);
    end
end
if isempty(rois) == 0
    for ia=1:length(rois)
        curoi = rois{ia};
        figureareas(ia) = figure;
        hst = suptitle(['Mean Power spectra ROI =',curoi]);
        fprintf('Printing Power spectra for ROI =',curoi );
        for ip=1:length(patientlist)
            eegpatient = patientlist{ip};
            [channelclassindexes] = getindexesfromlabel(eegpatient, curoi); %got the ids of the channels we want to display
            channelclassindexes = channelclassindexes -1;
            if isempty(channelclassindexes) < 1
                fprintf('Printing Power spectra for ROI = %s in Patient= %s', curoi, eegpatient);
                for ic=1:length(conditionslist)
                    eegcondition = conditionslist{ic};
                    %ha = tight_subplot(nbp,nbc,3*(ip-1) + ic,[.01 .01])
                    powerx = powerspecmatrix{ip,ic};
                    powerx = powerx([channelclassindexes],:)
                    meanvx = mean(powerx);
                    stdx = std(powerx);
                    maxofmean =max(meanvx);
                    freqsindexl = powerfreqsindexes{ip,ic};
                    figure(figureareas(ia));
                    subplot(nbp,nbc,3*(ip-1) + ic);
                    % plot mean as such
                    plot(freqsindexl,meanvx);
                    set(gca, 'YLim', [0, 1.1*maxofmean], 'XLim', [0 xlimtodisplay]);
                    ylabel('Power ')
                    msgt = sprintf('Mean Power spectra ROI= %s, %s, %s', curoi, eegpatient, eegcondition);
                    title(msgt);
                end
            else
                fprintf('Skipping Patient %s, she has not area %s\n',eegpatient,curoi);
            end
        end
    end
end

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
