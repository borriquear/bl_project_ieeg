% eegpatient = 'TWH034'; eegcond = 'EC_PRE' ;chanidx(1) = 2; chanidx(2) = 33
%[myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(eegpatient,eegcond);

function [] = spectralcoherence2channels(EEG, chanidx)
%spectralcoherence2channels
%IN EEG, chanidx(1), chanidx(2)
%%
[srate, min_freq, max_freq, num_frex, time, n_wavelet, half_wavelet, frex, s, num_cycles]= initialize_wavelet();
%overwrite num_cyces
freqs2use  = logspace(log10(min_freq),log10(max_freq),15); % 4-30 Hz in 15 steps
%num_cycles    = logspace(log10(4),log10(8),length(freqs2use));

stepsize = 100; %ms
epochtime = length(EEG.data);


% wavelet and FFT parameters
% time          = -1:1/EEG.srate:1;
half_wavelet  = (length(time)-1)/2;
%num_cycles    = logspace(log10(4),log10(8),length(freqs2use));
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;
% time in indices
%times2saveidx = dsearchn(EEG.times',times2save');
times2saveidx = EEG.times;
%baselineidx   = dsearchn(times2save',baselinetm');

% chanidx(1) = find(strcmpi(channel1,{EEG.chanlocs.labels}));
% chanidx(2) = find(strcmpi(channel2,{EEG.chanlocs.labels}));
%idxchannel1, idxchannel2
% data FFTs
data_fft1 = fft(reshape(EEG.data(chanidx(1),:,:),1,n_data),n_convolution);
data_fft2 = fft(reshape(EEG.data(chanidx(2),:,:),1,n_data),n_convolution);
% initialize

%spectcoher = zeros(length(freqs2use),length(times2save));
spectcoher = zeros(length(freqs2use),epochtime);
for fi=1:length(freqs2use)
    wavelet_fft = fft( exp(2*1i*pi*freqs2use(fi).*time) .* exp(-time.^2./(2*(num_cycles^2))) ,n_convolution);
    
    % phase angles from channel 1 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft1,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    sig1 = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
    
    % phase angles from channel 2 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft2,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    sig2 = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
    
    % compute power and cross-spectral power
%     spec1 = mean(sig1.*conj(sig1),2);
%     spec2 = mean(sig2.*conj(sig2),2);
%     specX = abs(mean(sig1.*conj(sig2),2)).^2;
%     spectcoher(fi,:) = specX./(spec1.*spec2);
%     %alternative notation for the same procedure, using the Euler-like expression: Me^ik
%     spec1 = mean(abs(sig1).^2,2);
%     spec2 = mean(abs(sig2).^2,2);
%     specX = abs(mean( abs(sig1).*abs(sig2) .* exp(1i*(angle(sig1)-angle(sig2))) ,2)).^2;
%     
%     %compute spectral coherence, using only requested time points
%     spectcoher(fi,:) = specX./(spec1.*spec2);
%     
%     %yet another equivalent notation, just FYI
%     spec1 = sum(sig1.*conj(sig1),2);
%     spec2 = sum(sig2.*conj(sig2),2);
%     specX = sum(sig1.*conj(sig2),2);
%     %spectcoher(fi,:) = abs(specX(times2saveidx)./sqrt(spec1(times2saveidx).*spec2(times2saveidx))).^2;
%     spectcoher(fi,:) = abs(specX ./ sqrt(spec1 .* spec2)).^2;
    
    
    %imaginary coherence
    spec1 = sum(sig1.*conj(sig1),2);
    spec2 = sum(sig2.*conj(sig2),2);
    specX = sum(sig1.*conj(sig2),2);
    %spectcoher(fi,:) = abs(imag(specX(times2saveidx)./sqrt(spec1(times2saveidx).*spec2(times2saveidx))));
    spectcoher(fi,:) = abs(imag(specX./sqrt(spec1.*spec2)));
end
end