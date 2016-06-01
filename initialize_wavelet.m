function [half_of_wavelet_size, n_wavelet, n_data, n_convolution, wavelet_cycles, wtime] =   initialize_wavelet(EEG_study)
% initialize_wavelet  initialize complex wavelet parameters
% wavelet parameters
wtime = -1:1/EEG_study.srate:1; % window width
wavelet_cycles= 7; %time-frequency prcision trade-off 
% 
%w = 2(wavelet_cycles/(2 pi freq)^2)
half_of_wavelet_size = (length(wtime)-1)/2;
% convolution parameters
n_wavelet     = length(wtime);
n_data        = EEG_study.pnts*EEG_study.trials;
n_convolution = n_wavelet+n_data-1;

end