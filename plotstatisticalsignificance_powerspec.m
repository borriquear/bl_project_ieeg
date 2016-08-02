function [] = plotstatisticalsignificance_powerspec(patientlist, conditionslist, powerspecmatrix, powerfreqsindexes, powerx_mean_perbands, rois )
%plotstatisticalsignificance_powerspec calculate p-values

% suptitle(['Mean and std Power spectra All areas']);
nbp = length(patientlist);
nbc = length(conditionslist);
nbc = 1; %compare conditions in one chart
xlimtodisplay = 30;
freq_bands = {'\delta'; '\theta'; '\alpha'; '\beta'; '\gamma'};
% % rois = {'HD','T','F', 'IH'};
%
%powerforcond = zeros(length(conditionslist),length(conditionslist), size(powerx,2));
powerforcond = {};
allrois = figure;
maxofmeanall = 0;
for ip=1:length(patientlist)
    eegpatient = patientlist{ip};
    for ic =1:length(conditionslist)
        eegcondition = conditionslist(ic);
        %ha = tight_subplot(nbp,nbc,3*(ip-1) + ic,[.01 .01])
        powerx = powerspecmatrix{ip,ic};
        meanvx = mean(powerx); stdx = std(powerx);maxofmean =max(meanvx);
        freqsindexl = powerfreqsindexes{ip,ic};
        powerforcond{ip,ic}= meanvx;
        figure(allrois);
        subplot(nbp,1,ip);
        hold on
        plot(freqsindexl, meanvx);
        ylabel('Power (dB)');
        if maxofmeanall < maxofmean
            maxofmeanall = maxofmean;
            set(gca, 'YLim', [0, 1.1*maxofmean], 'XLim', [0 xlimtodisplay]);
        end
        msgt = sprintf('Power spectra All ROIs Patient=%s',eegpatient);
        title(msgt);
    end
    legend(conditionslist);
end
%calculate p-values for powerforcond
% x = powerx_mean_perbands(:,1);
% xx1 = cat(1, x{:});
% y = powerx_mean_perbands(:,2);
% yy2= cat(1,y{:});
% %[h,p] = ttest2(xx1, yy2);
% [h p] = ttest2(mean(xx1), mean(yy2));
% fprintf('p-value for Power spectra in conditions % and % =%s \n',conditionslist{1},conditionslist{2}, num2str(p) )
if length(conditionslist) > 1
    for i=1:length(patientlist)
        fprintf('Calculating ttest2 for patient=%s, %s-%s\n',  patientlist{i}, conditionslist{1}, conditionslist{2});
        pause(1.0);
        [h, p]= ttest2(powerx_mean_perbands{i,1},powerx_mean_perbands{i,2});
        disp([h p])
    end
end
%Calculating ttest2 for patient=TWH030, EC_PRE-HYP
%          0         0         0         0         0    0.4009    0.7562    0.0655    0.1730    0.9251
%
% Calculating ttest2 for patient=TWH031, EC_PRE-HYP
%          0         0         0         0         0    0.7964    0.5792    0.7634    0.6787    0.7822
%
% Calculating ttest2 for patient=TWH033, EC_PRE-HYP
%          0         0         0         0         0    0.6853    0.5187    0.9559    0.1529    0.3895
%
% Calculating ttest2 for patient=TWH034, EC_PRE-HYP
%          0         0         0         0         0    0.3152    0.3146    0.3381    0.3731    0.3108
%
% Calculating ttest2 for patient=TWH037, EC_PRE-HYP
%     1.0000         0         0         0         0    0.0000    0.1627    0.4738    0.2731    0.3007
%
% Calculating ttest2 for patient=TWH038, EC_PRE-HYP
%          0         0         0         0         0    0.1974    0.0555    0.9106    0.9883    0.4745
%
% Calculating ttest2 for patient=TWH042, EC_PRE-HYP
%          0    1.0000    1.0000         0         0    0.3457    0.0020    0.0002    0.8920    0.7899
%
% Calculating ttest2 for patient=TWH043, EC_PRE-HYP
%          0    1.0000         0         0         0    0.5678    0.0153    0.9536    0.6819    0.4217

meanavgbandroipatcond = {};
if isempty(rois) == 0
    for ia=1:length(rois)
        curoi = rois{ia};
        %hst = suptitle(['Mean Power spectra ROI =',curoi]);
        for ip=1:length(patientlist)
            eegpatient = patientlist{ip};
            [channelclassindexes] = getindexesfromlabel(eegpatient, curoi); %got the ids of the channels we want to display
            channelclassindexes = channelclassindexes -1;
            if isempty(channelclassindexes) < 1
                skippatient = 0;
                for ic=1:length(conditionslist)
                    eegcondition = conditionslist{ic};
                    %ha = tight_subplot(nbp,nbc,3*(ip-1) + ic,[.01 .01])
                    powerx = powerspecmatrix{ip,ic};powerx = powerx([channelclassindexes],:); meanvx = mean(powerx);
                    stdx = std(powerx); maxofmean = max(meanvx);
                    freqsindexl = powerfreqsindexes{ip,ic};
                    
                    % plot mean as such
                    powerforcondroi{ip,ic}= meanvx;
                    power_avg_bands = powerx_mean_perbands{ip,ic};
                    power_avg_bands = power_avg_bands([channelclassindexes],:);
                    %meanavgbandroi = mean(power_avg_bands);
                    meanavgbandroipatcond{ip,ic} = power_avg_bands;
                    %normfqb = sum(mean(power_avg_bands));
                    %bar(mean(power_avg_bands)/normfqb);
                end
            else
                %close(figureareas(ia));
                skippatient = 1;
                fprintf('Skipping Patient %s, she has not area %s\n',eegpatient,curoi);
            end
            if skippatient < 1
                if length(conditionslist) > 1
                    fprintf('Calculating ttest2 for patient=%s, %s-%s in ROI=%s\n',  patientlist{ip}, conditionslist{1}, conditionslist{2}, curoi);
                    pause(1.0);
                    [h, p]= ttest2(meanavgbandroipatcond{ip,1},meanavgbandroipatcond{ip,2});
                    disp([h p]);
                end
            end
        end
    end
end
%     end
% end
end