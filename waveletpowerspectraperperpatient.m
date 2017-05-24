function [powerx_perband, averagemeanpower] = waveletpowerspectraperperpatient(eegpatientl, eegconditionl, temporalw)
% waveletpowerspectraperperpatient Calculates the powerspectra for a list
% of patients and conditions. It creates a .mat object that contains the
% poweer spectra for num_frex frequencies. It calls the initialize_EEG_variables to create the EEG object and the 
% initialize_wavelet functions. The rer 
%
% Input:
% eegpatientl list of patients
% eegconditionl list of conditions (EO, EC, HYP)
% Output: The return arguments are created in the plotpowerspectra function
% powerx_perband  power spectrum for each frequency band
% averagemeanpower
powerspectra_matrix_list = {};
for indpat=1:length(eegpatientl)
    eegpatient = eegpatientl{indpat};
    for indcond=1:length(eegconditionl)
        eegcondition = eegconditionl{indcond};
        [myfullname, EEG, channel_labels, patdate, patsession] =  initialize_EEG_variables(eegpatient,eegcondition)
        [srate, min_freq, max_freq, num_frex, time, n_wavelet, half_of_wavelet_size, frex, s, wavelet_cycles] =   initialize_wavelet();
        %baseidx = dsearchn(EEG.times',[0 length(EEG.times)]');
        baseidx = dsearchn(EEG.times',[0 3*srate]');
        trial2plot = EEG.trials; chani = 2; chanend = EEG.nbchan; tot_channels = chanend - chani + 1;
        n_data               = EEG.pnts*EEG.trials;
        n_convolution        = n_wavelet+n_data-1;
        n_conv_pow2          = pow2(nextpow2(n_convolution));
        % get FFT of data for eachchannel
        eegpowermean = zeros(tot_channels, num_frex, 1);
        for irow =chani:chanend
            chan2use = channel_labels(irow);
            %eegfft = fft(reshape(EEG.data(strcmpi(chan2use,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);
            eegfft = fft(reshape(EEG.data(irow,:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);
            
            % initialize
            %eegpower = zeros(tot_channels, num_frex, EEG.pnts); % frequencies X time X trials
            
            for fi=1:num_frex
                fprintf('Calculating wavelet for pat %s cond %s channel %s at frq = %s\n', eegpatient, eegcondition, chan2use{1}, num2str(fi))
                %wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
                wavelet = fft( sqrt(1/(7*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(7^2))) , n_conv_pow2 );
                % convolution
                eegconv = ifft(wavelet.*eegfft);
                eegconv = eegconv(1:n_convolution);
                eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
                
                % Average power over trials
                temppower = mean(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2);
                %eegpower(irow-1, fi,:) = 10*log10(temppower./mean(temppower(baseidx(1):baseidx(2))));
                %decibel conversion The decibel (dB) is a ratio between the strength of one signal (frequency-band-specific power)
                % and the strength of another signal (a baseline level of power in that same frequency band).
                % The decibel scale is convenient because it is robust to many limitations of raw power, such
                % as 1/ f power scaling and subject-specific and electrode-specific idiosyncratic characteristics.
                % The base unit is called a bel and is the logarithm of a ratio of numbers. Typically, tens of bels
                % are used, hence decibels.
                %eegpowermean(irow-1, fi) = mean(10*log10(temppower./mean(temppower(baseidx(1):baseidx(2)))));
                %eegpowermean(irow-1, fi) = mean(10*log10(temppower));
                eegpowermean(irow-1, fi) = mean(temppower); %Raw power
            end
        end
        fprintf('Updating powerspectra_matrix_list{indpat, indcond, fi} ')
        powerspectra_matrix_list{indpat, indcond} = eegpowermean;
    end
end
[powerx_perband, averagemeanpower] = plotpowerspectra(powerspectra_matrix_list, eegpatientl, eegconditionl, frex);
end

function [powerx_perband, averagemeanpower] = plotpowerspectra(powerspectra_matrix_list, patientlist, conditionslist, frex)
% plotpowerspectra Plots the power spectrum for different frequency bands 
% 
% Input:
% powerspectra_matrix_list matrix pats x conds
% frex is a list of frequencies
% Output: 
% powerx_perband, averagemeanpower

initfreq=3; %delta band from 3.5
ympb= []; condmeanrois = [];
%conditionslist = {'EC', 'EO'}
figm = figure;%powerx (plot)
nbp = length(patientlist);
nbc = length(conditionslist);
xlimtodisplay = 50;
matrixtest = powerspectra_matrix_list{1,1};
nbcols = size(matrixtest,2);
freq_bands = {'\delta'; '\theta'; '\alpha'; '\beta 1'; '\beta 2';'\gamma'};
nbfreqs = length(frex) -2;
%Charts ALL ROIs
ylistpat = zeros(length(patientlist),length(conditionslist), nbfreqs);
freqsindexl = frex(3:end);
%powerx_perband ={};
powerx_perband = zeros(length(patientlist),length(conditionslist),length(freqsindexl));
stdx_perband = zeros(length(patientlist),length(conditionslist),length(freqsindexl));
for ip=1:length(patientlist)
    eegpatient = patientlist{ip};
    for ic =1:length(conditionslist)
        %yel = [];
        eegcondition = conditionslist{ic};
        %ha = tight_subplot(nbp,nbc,3*(ip-1) + ic,[.01 .01])
        powerx = powerspectra_matrix_list{ip,ic}; meanvx = mean(powerx); stdx = std(powerx);
        powerx = powerx(:,initfreq:end); meanvx = meanvx(initfreq:end); stdx = stdx(initfreq:end);
        maxofmean =max(meanvx); minofmean = min(meanvx);
        %meanperfreq(fq) = zeros(1,length(freq_bands));
        for fq=1:length(freq_bands)
            powerx_perband(ip,ic,fq) = mean(powerx(:,fq));
            stdx_perband(ip,ic,fq) = std(powerx(:,fq));
            %meanperfreq(fq) = mean(powerx(:,fq));
        end
        figure(figm);
        subplot(nbp,nbc,nbc*(ip-1) + ic);
        % Figuer 1: plot mean as such
        plot(freqsindexl,meanvx);
        set(gca, 'YLim', [ minofmean 1.01*maxofmean], 'XLim', [0 xlimtodisplay]);
        ystr = sprintf('Power %s',eegpatient);
        ylabel(' Raw Power \mu V^2');
        msgt = sprintf('Power spectrum all ROIS %s', eegcondition);
        if ip == 1
            title(msgt, 'interpreter', 'none');
        end
    end
    %         yel = mean(powerx_mean_perbands{ip,ic})/normfqb;
    %         ylistpat(ip,ic,:) = yel;
end
%obtain  powerx_perband with the  mean power per band the mean all patients each band
powerspectoplot = zeros(length(conditionslist),length(freq_bands));
error_bars= zeros(length(conditionslist),length(freq_bands));
for ic = 1:length(conditionslist)
    for fq=1:length(freq_bands)
        powerspectoplot(ic,fq) = mean(powerx_perband(:,ic,fq));
        error_bars(ic,fq) = std(stdx_perband(:,ic,fq));
    end
    %normalize   
end
% plot the mean over subjects
figure
X = [powerspectoplot(1,:)/sum(powerspectoplot(1,:));powerspectoplot(2,:)/sum(powerspectoplot(2,:))]
Y = [error_bars(1,:)/sum(error_bars(1,:));error_bars(2,:)/sum(error_bars(2,:))]
bar(X')
set(gca, 'XTick', 1:length(freq_bands), 'XTickLabel',freq_bands);
ylim([0 1]);
xlabel('Frequency Bands'), ylabel('Normalized Power \muV^2 per Band')
msgt = sprintf('Mean Power spectra All patients-ROIS');
legend(conditionslist{1}, conditionslist{2}, 'interpreter', 'None')
% legend('EC','EO')
title(msgt,'interpreter', 'none');
%plot the imagesc per patient
% plot error bars
fig_eb = figure;   %powerxperband (bars)
errorbar(X', Y');
set(gca, 'XTick', 1:length(freq_bands), 'XTickLabel',freq_bands);
msgt = sprintf('Error bars of Power spectra All patients-ROIS');
title(msgt,'interpreter', 'none');


figure;
% eyes closed
cond1 = powerx_perband(:,1,:);
cond1 = reshape(cond1, [length(patientlist) length(freq_bands)]);
% eyes open
cond2 = powerx_perband(:,2,:);
cond2 = reshape(cond2, [length(patientlist) length(freq_bands)]);
subplot(1,2,1);
c1 = zeros(length(patientlist),length(freq_bands)); c2 = zeros(length(patientlist),length(freq_bands));
for i=1:length(patientlist)
    for j=1:length(freq_bands)
        c1(i,j) = cond1(i,j)/sum(cond1(i,:));
        c2(i,j) =  cond2(i,j)/sum(cond2(i,:));
    end
end
imagesc(c1);
set(gca,'XTick', [1:length(freq_bands)], 'YTick',[]);
set(gca,'XTickLabel', freq_bands);
ylabel('Patients'), xlabel('Frequency bands'); %cd= condicslabels(cc);
msgt= sprintf('Power normalized all ROIs in %s',conditionslist{1});
title(msgt,'interpreter', 'none')
caxis([0 1]);
colorbar;
subplot(1,2,2)
imagesc(c2);
set(gca,'XTick', [1:length(freq_bands)], 'YTick',[]);
set(gca,'XTickLabel', freq_bands);
ylabel('Patients'), xlabel('Frequency bands'); 
msgt= sprintf('Power normalized all ROIs in %s',conditionslist{2});
title(msgt,'interpreter', 'none')
caxis([0 1]);
colorbar;
averagemeanpower = X;
end