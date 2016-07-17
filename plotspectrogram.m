function [] = plotspectrogram(patientlist,conditionslist, rois )
% plotspectrogram plot powerxfreqxtime
global min_freq;
global max_freq;
global num_frex;
global srate;
global time;
global n_conv_pow2;
%global n_wavelet;
global n_convolution;
global half_of_wavelet_size;
global frex;
global s;
global baseidx;
global wavelet_cycles;
global globalFsDir;
[globalFsDir] = loadglobalFsDir();
[srate, min_freq, max_freq, num_frex, time, n_wavelet, half_of_wavelet_size, frex, s, wavelet_cycles]= initialize_wavelet();
fprintf('Wavelet parameters loaded\n')
% patientlist = {'TWH030'}%,'TWH031','TWH034','TWH033','TWH037','TWH038','TWH042','TWH043'};%,'TWH043','TWH037'
% %Patient 37 and 43  have NOT HD,  only Deep
% conditionslist = {'EC_PRE', 'HYP'}%,'EC_POST',}; %'EO_PRE' 'EO_POST'
% rois = {'F','HD'}%,'T', 'IH'}; %comment this for single channel, otherwise spectrogram for ROIs
singlechannel = isempty(rois);
for ip=1:length(patientlist)
    eegpatient = patientlist{ip};
    for ic=1:length(conditionslist)
        eegcondition = conditionslist{ic};
        [myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(eegpatient,eegcondition);
        % definte convolution parameters
        %         n_wavelet            = length(time);
        %srate = EEG.srate;
        n_data               = EEG.pnts*EEG.trials;
        n_convolution        = n_wavelet+n_data-1;
        n_conv_pow2          = pow2(nextpow2(n_convolution));
        %         half_of_wavelet_size = (n_wavelet-1)/2;
        % get FFT of data
        for iroi=1:length(rois)
            curoi = rois{iroi};
            % get channel indexes for curoi eg 'HD'
            [channelclassindexes] = getindexesfromlabel(eegpatient, curoi); %got the ids of the channels we want to display
            if isempty(channelclassindexes) == 0
                %         channelclassindexes = channelclassindexes -1;
                %chan2use = channel_labels{2};
                eegpowerperch = {};
                for iidx=1:length(channelclassindexes)
                    %eegfft = fft(reshape(EEG.data(strcmpi(chan2use,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);
                    indexc = channelclassindexes(iidx);
                    eegfft = fft(reshape(EEG.data(indexc,:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);
                    % call for one patient
                    eegpowerperch{iidx} = spectrogram_onechannel(EEG, eegfft);
                end
                %get mean of eegpowerperch
                %mean(A,2) is a column vector containing the mean of each row
                eepowersum = zeros(num_frex, EEG.pnts*EEG.trials );
                for k=1:iidx
                    eepowersum = eepowersum + eegpowerperch{iidx} ;
                end
                eegpower = eepowersum/k;
                eepowerallpatsconds{ip,ic,iroi}= eegpower;
                %plot spectrogams
                plotspectrogramperroi(eegpower,EEG,curoi,eegcondition,eegpatient);
            else
                fprintf('Skipping ROI =%s in Patient=%s', curoi,eegpatient);
            end
        end% rois
    end
end
end

function  eegpower = spectrogram_onechannel(EEG, eegfft)
%global min_freq;
%global max_freq;
global num_frex;
global srate;
global time;
%global n_wavelet;
global half_of_wavelet_size;
global frex;
global s;
global n_conv_pow2;
global n_convolution;
global baseidx;
global wavelet_cycles;
%global globalFsDir;
eegpower = zeros(num_frex,EEG.pnts); % frequencies X time X trials
%         baseidx = dsearchn(EEG.times',[-500 -200]');
% tlag = 10; %30 seconds of the last 2 minutes
% a = EEG.pnts -EEG.srate*tlag-1000; b = EEG.pnts -EEG.srate*tlag;
a = 1;b = a + 1000; %the BL is 2 seconds
baseidx = dsearchn(EEG.times',[a b]');

% loop through frequencies and compute synchronization
for fi=1:num_frex
    %variable wavelet cycles
    %wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
    % fixed wavelet cycles
    wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(wavelet_cycles^2))) , n_conv_pow2 );
    % convolution
    eegconv = ifft(wavelet.*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
    % Average power over trials (this code performs baseline transform,
    % which you will learn about in chapter 18)
    temppower = mean(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2);
    eegpower(fi,:) = 10*log10(temppower./mean(temppower(baseidx(1):baseidx(2))));
end
end

function [] = plotspectrogramperroi(eegpower,EEG,curoi, eegc,eegp)
global frex;
global min_freq;
global max_freq;
global baseidx;
figure;
subplot(121);
contourf(EEG.times,frex,eegpower,40,'linecolor','none');
%set(gca,'clim',[min(min(eegpower))  max(max(eegpower))],'xlim',[EEG.times(1) EEG.times(end)],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
%set(gca,'clim',[min(eegpower(:))  max(eegpower(:))],'xlim',[EEG.times(end -1000) EEG.times(end)],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
%set(gca,'clim',[min(eegpower(:))  max(eegpower(:))],'xlim', baseidx,'xtick',baseidx ,'xticklabel',round(baseidx)/1000,'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
set(gca,'clim',[min(eegpower(:))  0.5*max(eegpower(:))],'xlim', [EEG.times(1) EEG.times(end)],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10);

xlabel('Time(s)'),ylabel('Log Frequency(Hz)');
title(['ROI=' curoi,' ',eegc,' ', eegp])
colorbar
subplot(122)
contourf(EEG.times,frex,eegpower,40,'linecolor','none')
%set(gca,'clim',[-3 3],'xlim',[EEG.times(1) EEG.times(end)])
set(gca,'clim',[min(eegpower(:))  max(eegpower(:))],'xlim', [90000 92000] ,'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10);
xlabel('Time(ms)'),ylabel('Log Frequency(Hz)');
title(['ROI=' curoi,' ',eegc,' ', eegp, ' at 90s 200ms'])
colorbar
end