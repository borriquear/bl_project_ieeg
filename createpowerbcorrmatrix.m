function powerconn_matrix = createpowerbcorrmatrix(eegpatientl, eegconditionl, centerfrequenciesl)
%createpowerbcorrmatrix calls to powerbasedconnectivity(eegpatient, eegcondition, centerfreq)
% to calculate correlation analysis per patient and frequency band
%OUT: mat file, powerconnectivity_freq_pat_cond  -- corr_matrix', 'channel_labels

%eegconditionl = {'EC_PRE', 'EO_PRE'}; 
%eegpatientl =  {'TWH030', 'TWH031', 'TWH033','TWH037','TWH038','TWH042'};
%[srate, min_freq, max_freq, num_frex, time, n_wavelet, half_wavelet, freqs2use, s, wavelet_cycles]= initialize_wavelet();
%centerfrequenciesl  = logspace(log10(min_freq),log10(max_freq),8);
global globalFsDir;
globalFsDir = loadglobalFsDir();
powerconn_matrix = struct;
%open file with power conn matrices all patients, conds and freqs
powerconnmatf = fullfile(globalFsDir, 'powerconn_matrices.mat');
%fh = load(powerconnmatf);
%corr_matrix_list = zeros(length(eegpatientl),length(eegconditionl),length(centerfrequenciesl));
corr_matrix_list = {};
if exist(powerconnmatf, 'file') ~= 2
    %if file doesnt exist, create file from scratch
    fprintf('Creating %s from scratch \n', powerconnmatf)
    for indpat=1:length(eegpatientl)
        eegpatient = eegpatientl{indpat};
        for indcond=1:length(eegconditionl)
            eegcondition = eegconditionl{indcond};
            for indexfreq = 1:length(centerfrequenciesl) % delta, theta, alpha, beta, gamma
                centerfreq = centerfrequenciesl(indexfreq);
                fprintf('Calling to powerbasedconnectivityperpatient %s %s %s',eegpatient, eegcondition, num2str(centerfreq) )
                corr_matrix = createpowerbcorrmatrixperpatient(eegpatient, eegcondition, centerfreq);
                corr_matrix_list{indpat,indcond,indexfreq} = corr_matrix;
                %corr_matrix_list(indpat,indcond,indexfreq) = corr_matrix;
            end
        end
    end
    powerconn_matrix.power_matrix = corr_matrix_list;
    powerconn_matrix.patientsl = eegpatientl;
    powerconn_matrix.conditionsl = eegconditionl;
    powerconn_matrix.freqsl = centerfrequenciesl;
    save(powerconnmatf,'powerconn_matrix');
else
    %added matrices from new patient or condition
end
end


function [ corr_matrix ] = createpowerbcorrmatrixperpatient(eegpatient, eegcondition, centerfreq)
% powerbasedconnectivityperpatient Returns the correlation matrix for the time frequency power correlation
% between all pair of electrodes for patient and condition for each
% frequency band
% IN: patient, condition,centerfreq e.g. patient='TWH020', condition =
% 'HYP', centerfreq = 2
% OUT: corr_matrix nxn correlation matrix, where n is the number of
% electrodes of the patient (powerconnectivity_freq_)This file is needed in
% displaypowerconnectivity
%% 1. Load epoch and Initialize data
powerconnperpatient = 0;
fprintf('Loading EEG for patient:%s and condition%s\n',eegpatient,eegcondition);
[myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(eegpatient,eegcondition);
[eegpathname eegfilename eegextname]= fileparts(myfullname);
%times2plot =  dsearchn(EEG.times',[ EEG.times(1) EEG.times(end)]');
trial2plot = EEG.trials;
chani = 2;
chanend = EEG.nbchan;
tot_channels = chanend - chani + 1;
%scaledlog = 1 10*log(abs(fft(signal))
%scaledlog = 0;
disp(['Total number of channels=' num2str(EEG.nbchan)]);
srate = EEG.srate;
disp(['Sampling rate=', num2str(srate)]);
%total seconds of the epoch
tot_seconds = floor(EEG.pnts/srate);
disp(['Total seconds of the session=' num2str(tot_seconds)])
% remaining points beyond tot_seconds
remainingpoints = EEG.pnts - tot_seconds*srate;
%time vector normalized
t = 0:1/srate:(tot_seconds + ((remainingpoints/srate)) );
% end -1 to have samne length as signal
t = t(1:end-1);
% t =horzcat(t,restsecondslist);
%nyquistfreq = srate/2;
%baseidx = dsearchn(EEG.times',[1000 2000]');

%initialize the correlation between any two pairs
corr_matrix =zeros(tot_channels,tot_channels);
corrtype = 'Spearman';

%[half_of_wavelet_size, n_wavelet, n_data n_convolution, wavelet_cycles,wtime] = initialize_wavelet(EEG);
[srate, min_freq, max_freq, num_frex, wtime, n_wavelet, half_of_wavelet_size, frex, s, wavelet_cycles]= initialize_wavelet();
n_data               = EEG.pnts*EEG.trials;
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
%n = length(t); % = EEG.pnts  = n_wavelet
for irow =chani:chanend
    ichan2use = channel_labels(irow);
    fft_data1 = fft(reshape(EEG.data(strcmpi(ichan2use,{EEG.chanlocs.labels}),:,:),1,n_data),n_convolution);
    for jcol =chani:chanend
        jchan2use = channel_labels(jcol);
        fprintf('Calculating in patient %s Freq %s the %s correlation between channels %s and %s \n',eegpatient,num2str(centerfreq), corrtype,  ichan2use{1}, jchan2use{1});
        r = 0;
        %for windowsc=1:tot_seconds
        %time = -1:1/EEG.srate:1;
        %time = windowsc-1:1/EEG.srate:windowsc;
        %n_convolution = length(time);
        
        fft_data2 = fft(reshape(EEG.data(strcmpi(jchan2use,{EEG.chanlocs.labels}),:,:),1,n_data),n_convolution);
        %We use only one wavelet because we are dpoing the t-f decomposition for one frequency band (centrerfrequency) if we have no hypothesis about the frequecnie swe need to use a vector of frequencies, freq_min, freq_max
        fft_wavelet = fft(exp(2*1i*pi*centerfreq.*wtime) .* exp(-wtime.^2./(2*( wavelet_cycles /(2*pi*centerfreq))^2)),n_convolution);
        convolution_result_fft  = ifft(fft_wavelet.*fft_data1,n_convolution) * sqrt(wavelet_cycles /(2*pi*centerfreq));
        convolution_result_fft  = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
        convolution_result_fft1 = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
        %convolution_result_fft1 = reshape(convolution_result_fft,length(time),EEG.trials);
        %fft_wavelet             = fft(exp(2*1i*pi*centerfreq.*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*centerfreq))^2)),n_convolution);
        convolution_result_fft  = ifft(fft_wavelet.*fft_data2,n_convolution) * sqrt(wavelet_cycles /(2*pi*centerfreq));
        convolution_result_fft  = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
        convolution_result_fft2 = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
        %convolution_result_fft2 = reshape(convolution_result_fft,length(time),EEG.trials);
        r =  corr(abs(convolution_result_fft1(:,trial2plot)).^2,abs(convolution_result_fft2(:,trial2plot)).^2,'type','s','rows','pairwise');
        %rmean = r;
        %disp(rmean);
        fprintf('Spearman mean between channels %s and %s  is %.4f\n',ichan2use{1},jchan2use{1},r);
        corr_matrix(irow-1,jcol-1) = r;
    end
end
%entropym = sum(corr_matrix(corr_matrix~=0).*log(corr_matrix(corr_matrix~=0)));
if powerconnperpatient == 1
    %save mat file for each patient in the patient directory, this is
    %redundant since we have all in globalFsDir/powerconn_matrices.mat
    [globalFsDir] = loadglobalFsDir();
    patpath = strcat(globalFsDir,eegpatient);
    mattoload = strcat('powerconnectivity_freq_',num2str(centerfreq),'_',eegcondition,'_', eegpatient,'_',eegdate,'_',eegsession,'.mat');
    fftfile = fullfile(patpath,'data','figures', mattoload);
    save(fftfile,'corr_matrix', 'channel_labels')
end
end