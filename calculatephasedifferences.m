function [] = calculatephasedifferences(patientslist, conditionslist)
global min_freq;
global max_freq;
global num_frex;
global srate;
global time;
%global n_wavelet;
global n_convolution;
global half_wavelet;
global freqs2use;
global s;
global wavelet_cycles;
global num_cycles; %these 2 vars are the same, use one or the other
global timewindow;
global globalFsDir;
[globalFsDir] = loadglobalFsDir();
[srate, min_freq, max_freq, num_frex, time, n_wavelet, half_wavelet, freqs2use, s, wavelet_cycles]= initialize_wavelet();
freqs2use  = logspace(log10(min_freq),log10(max_freq),8);
num_cycles    = logspace(log10(4),log10(8),length(freqs2use));
%uncomment this and comment previous line for fixednumber of wavelet cycles
%wavelet_cycles = 7
fprintf('Wavelet parameters loaded\n')
timewindow = linspace(1.5,3,length(freqs2use)); % number of cycles on either end of the center point (1.5 means a total of 3 cycles))
ispc_matrix{length(patientslist), length(conditionslist)} = [];
phase_matrix{length(patientslist), length(conditionslist),length(freqs2use)} = [];
for ip=1:length(patientslist)
    eegpatient = patientslist{ip};
    for ic=1:length(conditionslist)
        eegcondition = conditionslist{ic};
        [myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(eegpatient,eegcondition);
        % initialize
        %times2save = EEG.pnts;
        n_data               = EEG.pnts*EEG.trials;
        n_convolution        = n_wavelet+n_data-1;
        %ispc_matrix{ip, ic} = zeros(length(freqs2use),length(times2save));
        %ps_matrix{ip, ic}  = zeros(length(freqs2use),length(times2save));
        chani = 2; %initial channel (chani = 1 is reserved for the 'Event' channel)
        chanend = EEG.nbchan; %end channel
        phase_matrix = {};
        for irow =chani:chanend
            ichan2use = channel_labels{irow};
            for jcol =chani:chanend
                jchan2use = channel_labels{jcol};
                %calculate phase difference for two channels
                data_fft1 = fft(reshape(EEG.data(irow,:,:),1,n_data),n_convolution);
                data_fft2 = fft(reshape(EEG.data(jcol,:,:),1,n_data),n_convolution);
                for fi=1:length(freqs2use)
                    % create wavelet and take FFT
                    phaseij = calculatephasedifferences2channels(EEG, data_fft1,data_fft2, fi );
                    phase_matrix{irow, jcol, fi} = phaseij;
                    fprintf(' ISPC=%.3f . %s:%s, channels %s-%s frq=%d\n', phaseij, eegpatient, eegcondition, ichan2use, jchan2use, freqs2use(fi));
                end % end frequency loop
            end
        end   
        ispc_matrix{ip,ic} = phase_matrix;
    end
end

end

function [phaseij] = calculatephasedifferences2channels(EEG, data_fft1,data_fft2,fi)
%
global freqs2use;
global n_convolution;
global half_wavelet;
global num_cycles;
global time;
global timewindow;
s = num_cycles(fi)/(2*pi*freqs2use(fi));
wavelet_fft = fft( exp(2*1i*pi*freqs2use(fi).*time) .* exp(-time.^2./(2*(s^2))) ,n_convolution);
% phase angles from channel 1 via convolution
convolution_result_fft = ifft(wavelet_fft.*data_fft1,n_convolution);
convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
phase_sig1 = angle(reshape(convolution_result_fft,EEG.pnts,EEG.trials));
% phase angles from channel 2 via convolution
convolution_result_fft = ifft(wavelet_fft.*data_fft2,n_convolution);
convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
phase_sig2 = angle(reshape(convolution_result_fft,EEG.pnts,EEG.trials));
% phase angle differences
phase_diffs = phase_sig1-phase_sig2;

% compute ICPS over trials
%ps(fi,:) = abs(mean(exp(1i*phase_diffs(times2saveidx,:)),2));
ps(fi,:) = abs(mean(exp(1i*phase_diffs(:)),2));
% compute time window in indices for this frequency
time_window_idx = round((1000/freqs2use(fi))*timewindow(fi)/(1000/EEG.srate));
%     time_window_idx = round(300/(1000/EEG.srate)); % set 300 to 100 for figure 3c/d
for ti=time_window_idx+1:EEG.pnts-time_window_idx
    % compute phase synchronization
    phasesynch = abs(mean(exp(1i*phase_diffs(ti-time_window_idx:ti+time_window_idx,:)),1));
    % average over trials
    ispc(ti) = mean(phasesynch);
end
phaseij = mean(ispc(:));
end