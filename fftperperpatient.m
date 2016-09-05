function filematname = fftperperpatient(eegpatient, eegcond)
%IN patient, condition
%OUT matfile handle
%
global t;
global tot_seconds;
global hz;
global scaledlog;
global f; %frequencies vector
global freq_bands;
global  globalFsDir;
global plotsinglechannel;
globalFsDir = loadglobalFsDir();
fprintf('Calling to fftperperpatient for patient: %s and condition: %s\n', eegpatient,eegcond);
[myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(eegpatient,eegcond);
fprintf('EEG object for patient: %s and condition: %s initialized \n', eegpatient,eegcond);
% nb of trials is 1 in Bial database
trial2plot = EEG.trials;
% set the channels we need to work on
chani = 2; %initial channel (chani = 1 is reserved for the 'Event' channel)
chanend = EEG.nbchan; %end channel
tot_channels = chanend - chani+1;
fprintf('Total number of channels =%s\n', num2str(EEG.nbchan))
srate = EEG.srate;
tot_seconds = floor(EEG.pnts/srate);
fprintf('Total seconds= %s for %s %s\n',num2str(tot_seconds), eegpatient, eegcond);
% remaining data points beyond tot_seconds (1 second has srate data points)
remainingpoints = EEG.pnts - tot_seconds*srate;
%time vector normalized
t = 0:1/srate:(tot_seconds + (remainingpoints/srate));
% end -1 to have samne length as signal
t = t(1:end-1);
% t =horzcat(t,restsecondslist);
n = length(t); % = EEG.pnts
nyquistfreq = srate/2;
% There is one half only of frequencies of fourier coefficients
hz = linspace(0,nyquistfreq,floor(n/2)+1);
f = 0:50;
frex_idx = sort(dsearchn(hz',f'));
freq_bands = {'\delta'; '\theta'; '\alpha'; '\beta'; '\gamma'};
% initialize vecgor with power coefficients that result form FFT of the signal
ampli_fft = zeros(tot_channels, length(hz));
power_fft = zeros(tot_channels, length(hz));
power_fft_mean_perband = zeros(tot_channels, length(freq_bands));
scaledlog = 0; %scale the amplitude and the power
plotsinglechannel = ones(chanend); %plot every single channel
plotsinglechannel = zeros(chanend);
for irow =chani:chanend
    %irow =1 is for 'Event'
    chan2use = channel_labels{irow};
    % check if channel has interictal spikes, if yes still use it
    isinterictalspikes = interictalspikeschannels(eegpatient, irow-1);
    if isinterictalspikes == 1
        fprintf('Fund interictal in patient:, condition%s, channel=%s!!\n', eegpatient, eegcond, chan2use);
    end
    % Build the signal vector from EEG.data object
    signal = EEG.data(irow,:,trial2plot);
    signalX = zeros(size(signal));
    % fourier time is not the time in seconds but the time normalized
    fouriertime = (0:n)/n;
    fouriertime = fouriertime(1:end-1);
    fprintf('Calculating FFT for Patient %s, Cond %s, channel = %s, channel_id=%d / channel_total=%d ...\n', eegpatient, eegcond, chan2use, irow-1, chanend-1);
    %FFT fft(X) computes the discrete Fourier transform (DFT) of X using a fast Fourier transform (FFT) algorithm.
    signalXF = fft(signal)/n;
    msgtitle = sprintf('Cond=%s Channel=%s Patient=%s \n',eegcond, chan2use,eegpatient);
    ampli_fft(irow-1,:) = 2*abs(signalXF(1:length(hz)));
    power_fft(irow-1,:) = abs(signalXF(1:length(hz))).^2;
    %extract information about specific frequencies
    power_fft_perband = abs(signalXF(frex_idx)).^2; %frex_idx = sort(dsearchn(hz',f'));
    if scaledlog == 1
        ampli_fft = scalebylogarithm(ampli_fft);
        power_fft = scalebylogarithm(power_fft);
        power_fft_perband = scalebylogarithm(power_fft_perband);
    end
    % plot signal and power per patient
    bandvaluesavg = plotsignalpowerperchannel(eegpatient,eegcond,chan2use, irow, signal, power_fft,power_fft_perband);
    power_fft_mean_perband(irow-1,:) = bandvaluesavg;
end
% saving the matfile
patdir = fullfile(globalFsDir, eegpatient, 'data\figures');
if ~exist(patdir, 'dir')
    mkdir(patdir);
end
fprintf('Saving the mat file with power vector in %s \n',patdir);
filematname = sprintf('fft_%s_%s_%s_%s.mat',eegcond, eegpatient, eegdate, eegsession);
filematname = fullfile(patdir,filematname);
save(filematname, 'ampli_fft','power_fft','power_fft_perband', 'power_fft_mean_perband','channel_labels');
fprintf('DONE mat file with power vector in %s \n', filematname);
end

function [bandvaluesavg] = plotsignalpowerperchannel(eegpatient,eegcond,chan2use, irow, signal, power_fft, power_fft_perband)
%plot signal and power extracted from FFT for a single channel
global t;
global tot_seconds;
global hz;
global scaledlog;
global freq_bands;
global plotsinglechannel;
%indexes of frequency bands in the hz vector
%db = dsearchn(hz',[0:3.5]'); tb = dsearchn(hz',[4:7]'); ab  = dsearchn(hz',[8:12]');  bb = dsearchn(hz',[13:35]');gb  = dsearchn(hz',[36:50]');
db = [1:4];tb = [5:8]; ab = [9:12]; bb = [13:35]; gb = [36:51];
dbv = mean(power_fft_perband(1,db)); tbv = mean(power_fft_perband(1,tb));
abv = mean(power_fft_perband(1,ab)); bbv = mean(power_fft_perband(1,bb)); gbv = mean(power_fft_perband(1,gb));
bandvaluesavg = [dbv tbv abv bbv gbv];
if plotsinglechannel(irow) == 1
    fprintf('Showing chart for patient:%s, condition:%s, channel:%s\n', eegpatient, eegcond, chan2use);
    msgtitle = sprintf('Patient=%s Cond=%s Channel=%s  ',eegpatient, eegcond, chan2use);
    figure;
    subplot(3,1,1)
    plot(t, signal, 'b')
    xlabel('Time (s)'), ylabel('Signal Amplitude')
    set(gca,'xlim',[0 tot_seconds])
    title(msgtitle)
    subplot(3,1,2)
    plot(hz,power_fft(irow-1,:),'b')
    xlabel('Frequencies (Hz)'), ylabel('Power')
    if  scaledlog == 1
        ylabel('log Power (dB)')
    end
    set(gca,'xlim',[0 50])
    legend({'FFT, Power'})
    subplot(3,2,5)
    bar(power_fft_perband(1,:))
    xlabel('Frequency Bands'), ylabel('Power')
    set(gca, 'XLim', [0 30], 'XTick',0:4:30,'XTickLabel',0:4:30);
    % aggregate per frequency band, mean value of the histogramain the
    % frequency band range
    subplot(3,2,6);
    bar(bandvaluesavg);
    xlabel('Frequency Bands'), ylabel('Power per Band')
    set(gca, 'XTickLabel',freq_bands);
end
end

function [scaledsignal] = scalebylogarithm(signal)
%% scalebylogarithm scaledsignal = 10*log10(signal)
scaledsignal = 10*log(signal);
a = 0;
end