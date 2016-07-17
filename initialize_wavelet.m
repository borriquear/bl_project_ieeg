function [srate, min_freq, max_freq, num_frex, time, n_wavelet, half_of_wavelet_size, frex, s, wavelet_cycles] =   initialize_wavelet()
min_freq =  1;
max_freq = 50;
num_frex = 30;
srate = 1000;
time = -1:1/srate:1;
n_wavelet = length(time);
half_of_wavelet_size = (n_wavelet-1)/2;
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
% number of wavelet cycles changes in function of the frequency
s = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
% fixed number of wavelet cycles between 7 and 10 (long epochs)
wavelet_cycles = 7; %Abs min max 4-14
end