function powerbasedconnectivityall()
% powerbasedconnectivityall Calls to powerbasedconnectivity(eegpatient, 
%eegcondition, centerfreq)
%Connectivity based on power correlation analysis per patient
eegcondition = 'HYP';
centerfrequencies = {2, 6 , 10, 23.5, 40};

%eegpatientl = { 'TWH027','TWH024','TWH028','TWH030', 'TWH031','TWH033','TWH034'};
%eegpatientl = {'TWH031','TWH033','TWH034'};
eegpatientl = {'TWH027' ,'TWH024' , 'TWH033','TWH034'};
%sigue por aqui
centerfrequencies = {10, 23.5, 40};
eegpatientl = {'TWH034'};


for indpat=1:length(eegpatientl)
    for indexfreq = 1:length(centerfrequencies) % delta, theta, alpha, beta, gamma
        centerfreq = centerfrequencies{indexfreq};
        powerbasedconnectivity( eegpatientl{indpat}, eegcondition,centerfreq );
    end
end
end
function [ corr_matrix ] = powerbasedconnectivity( eegpatient, eegcondition, centerfreq)
%powerbasedconnectivity Power based connectivity analysis
% Returns the correlation matrix for the time frequency power correlation
% between all pair of electrodes for patient and condition
% IN: patient, condition e.g. patient='TWH020', condition = 'HYP'
% OUT: corr_matrix nxn correlation matrix, where n is the number of
% electrodes of the patient (powerconnectivity_freq_)This file is need in
% displaypowerconnectivity
%% 1. Load epoch and Initialize data
disp('Loading the EEG data....\n')
%eegpatient= 'TWH030';
%eegcond = 'HYP';
fprintf('Loading EEG for patient:%s and condition%s\n',eegpatient,eegcondition);
[myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(eegpatient,eegcondition);
% get condition, sesssion, patient and date
[eegpathname eegfilename eegextname]= fileparts(myfullname);
times2plot =  dsearchn(EEG.times',[ EEG.times(1) EEG.times(end)]');
trial2plot = EEG.trials;
%[pathstr, subjandcond, matext] = fileparts(myfullname);
% Choose the channel(s) we want to analyze
chani = 2;
chanend = EEG.nbchan;
tot_channels = chanend - chani+1;
figcounter = 1;
hfigures = [];
%scaledlog = 1 10*log(abs(fft(signal))
scaledlog = 0;
disp('Data Loaded OK!');

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
n = length(t); % = EEG.pnts
nyquistfreq = srate/2;
%clculate the correlation between any two pairs
corr_matrix =zeros(tot_channels,tot_channels);
corrtype = 'Spearman';

[half_of_wavelet_size, n_wavelet, n_data n_convolution, wavelet_cycles] = initialize_wavelet(EEG);
for irow =chani:chanend
    ichan2use = channel_labels(irow);
    for jcol =chani:chanend
        jchan2use = channel_labels(jcol);
        fprintf('Calculating in patient %s Freq %s the %s correlation between channels %s and %s \n',eegpatient,num2str(centerfreq), corrtype,  ichan2use{1}, jchan2use{1});
        r = 0;
        %for windowsc=1:tot_seconds
        time = -1:1/EEG.srate:1;
        %time = windowsc-1:1/EEG.srate:windowsc;
        %n_convolution = length(time);
        fft_data1 = fft(reshape(EEG.data(strcmpi(ichan2use,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_convolution);
        fft_data2 = fft(reshape(EEG.data(strcmpi(jchan2use,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_convolution);
        fft_wavelet = fft(exp(2*1i*pi*centerfreq.*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*centerfreq))^2)),n_convolution);
        %for windowsc=1:tot_seconds
        % create wavelet and run convolution
        %time = windowsc-1:1/EEG.srate:windowsc;
        %             inittw = (length(time) -1)*(windowsc -1) + 1;
        %             endtw = inittw + length(time) -1;
        %             fft_data1 = fft_data1_orig(inittw:endtw);
        %             fft_data2 = fft_data2_orig(inittw:endtw);
        % create wavelet and run convolution
        fft_wavelet             = fft(exp(2*1i*pi*centerfreq.*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*centerfreq))^2)),n_convolution);
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
        %disp(rforwindow);
        %             if isnan(rforwindow)
        %                 rforwindow =1;
        %             end
        %             r= r + rforwindow;
        %r = corr(abs(convolution_result_fft1(:,trial2plot)).^2,abs(convolution_result_fft2(:,trial2plot)).^2,'type','s');
        %fprintf(); legend([ 'Spearman Rho = ' num2str(r) ]);
        %         end
        %         rmean = r/tot_seconds;
        rmean = r;
        %disp(rmean);
        fprintf('Spearman mean between channels %s and %s  is %.4f\n',ichan2use{1},jchan2use{1},rmean);
        corr_matrix(irow-1,jcol-1) = rmean;
    end
end


%entropym = sum(corr_matrix(corr_matrix~=0).*log(corr_matrix(corr_matrix~=0)));
[globalFsDir] = loadglobalFsDir();
patpath = strcat(globalFsDir,eegpatient);
mattoload = strcat('powerconnectivity_freq_',num2str(centerfreq),'_',eegcondition,'_', eegpatient,'_',eegdate,'_',eegsession,'.mat');
fftfile = fullfile(patpath,'data','figures', mattoload);
fprintf('Appending percentlistoffrqperband at file %s\n',fftfile);
save(fftfile,'corr_matrix', 'channel_labels')
end

