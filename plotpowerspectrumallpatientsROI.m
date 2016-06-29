function [] = plotpowerspectrumallpatientsROI(patientlist, conditionslist, powerspecmatrix, powerfreqsindexes, powerx_mean_perbands, rois )
%plotpowerspectrumallpatientsROI plot all patients mean of all channels in one figure
%figure1=figure('Position', [100, 100, 1024, 1200]);
figm = figure;%powerx (plot)
figp = figure;   %powerxperband (bars)
% figstd = figure;
% suptitle(['Mean and std Power spectra All areas']);
nbp = length(patientlist);
nbc = length(conditionslist);
xlimtodisplay = 30;
freq_bands = {'\delta'; '\theta'; '\alpha'; '\beta'; '\gamma'};
% % rois = {'HD','T','F', 'IH'};
%
for ip=1:length(patientlist)
    eegpatient = patientlist{ip};
    for ic =1:length(conditionslist)
        eegcondition = conditionslist{ic};
        %ha = tight_subplot(nbp,nbc,3*(ip-1) + ic,[.01 .01])
        powerx = powerspecmatrix{ip,ic};meanvx = mean(powerx);stdx = std(powerx);
        maxofmean =max(meanvx);
        freqsindexl = powerfreqsindexes{ip,ic};
        figure(figm);
        subplot(nbp,nbc,nbc*(ip-1) + ic);
        % plot mean as such
        plot(freqsindexl,meanvx);
        set(gca, 'YLim', [0, 1.1*maxofmean], 'XLim', [0 xlimtodisplay]);
        ylabel('Power ')
        msgt = sprintf('Mean Power spectra All ROIS %s %s', eegpatient, eegcondition);
        title(msgt);
        %plot bars of plot per band
        figure(figp);
        subplot(nbp,nbc,nbc*(ip-1) + ic);
        % plot mean as such
        %remove the outliers
        bar(mean(powerx_mean_perbands{ip,ic}));
        %normalized
        normfqb = sum(mean(powerx_mean_perbands{ip,ic}));
        bar(mean(powerx_mean_perbands{ip,ic})/normfqb);
        set(gca, 'XTickLabel',freq_bands);
        xlabel('Frequency Bands'), ylabel('Normalized Power per Band')
        msgt = sprintf('Mean Power spectra per Band All ROIS %s %s', eegpatient, eegcondition);
        title(msgt);
    end
end

if isempty(rois) == 0
    for ia=1:length(rois)
        figureareas(ia) = figure;
        figureareas_bars(ia) = figure;
        curoi = rois{ia};
        %hst = suptitle(['Mean Power spectra ROI =',curoi]);
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
                    powerx = powerx([channelclassindexes],:);
                    meanvx = mean(powerx);
                    stdx = std(powerx);
                    maxofmean =max(meanvx);
                    freqsindexl = powerfreqsindexes{ip,ic};
                    figure(figureareas(ia));
                    subplot(nbp,nbc,nbc*(ip-1) + ic);
                    % plot mean as such
                    plot(freqsindexl,meanvx);
                    set(gca, 'YLim', [0, 1.1*maxofmean], 'XLim', [0 xlimtodisplay]);
                    ylabel('Power ')
                    msgt = sprintf('Mean Power spectra ROI= %s, %s, %s', curoi, eegpatient, eegcondition);
                    title(msgt);
                    %plot bars avg frequency
                    figure(figureareas_bars(ia));
                    power_avg_bands = powerx_mean_perbands{ip,ic};
                    subplot(nbp,nbc,nbc*(ip-1) + ic);
                    %absolute value
                    %bar(mean(power_avg_bands));
                    %normalized per band
                    % plot mean as such
                    normfqb = sum(mean(power_avg_bands));
                    bar(mean(power_avg_bands)/normfqb);
                    set(gca, 'XTickLabel',freq_bands);
                    xlabel('Frequency Bands'), ylabel('Power per Band')
                    msgt = sprintf('Mean Power spectra per Band ROI= %s %s %s', curoi, eegpatient, eegcondition);
                    title(msgt);
                end
            else
                %close(figureareas(ia));
                fprintf('Skipping Patient %s, she has not area %s\n',eegpatient,curoi);
            end
        end
    end
end
%     end
% end
end