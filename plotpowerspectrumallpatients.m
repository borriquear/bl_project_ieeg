function [] = plotpowerspectrumallpatients(patientlist,conditionslist,powerspecmatrix,powerfreqsindexes )
%plot all patients mean of all channels in one figure
%figure1=figure('Position', [100, 100, 1024, 1200]);
figm = figure;
suptitle(['Mean Power spectra All areas']);
figstd = figure;
suptitle(['Mean and std Power spectra All areas']);
nbp = length(patientlist);
nbc = length(conditionslist);
xlimtodisplay = 30;
rois = {'HD','T','F', 'IH'};

for ip=1:length(patientlist)
    eegpatient = patientlist{ip}
    for ic =1:length(conditionslist)
        eegcondition = conditionslist{ic}
        %ha = tight_subplot(nbp,nbc,3*(ip-1) + ic,[.01 .01])
        powerx = powerspecmatrix{ip,ic};
        meanvx = mean(powerx);
        stdx = std(powerx);
        maxofmean =max(meanvx);
        freqsindexl = powerfreqsindexes{ip,ic};
        figure(figm);
        subplot(nbp,nbc,3*(ip-1) + ic);
        % plot mean as such
        plot(freqsindexl,meanvx);
        set(gca, 'YLim', [0, 1.1*maxofmean], 'XLim', [0 xlimtodisplay]);
        ylabel('Power ')
        msgt = sprintf('Mean Power spectra from FFT %s %s', eegpatient, eegcondition);
        title(msgt);
        % plot mean and error bars
        figure(figstd);
        subplot(nbp,nbc,3*(ip-1) + ic);
        errorbar(freqsindexl,meanvx,stdx);
        maxofmeanandstd = max(maxofmean, max(stdx));
        set(gca, 'YLim', [0, 1.1*maxofmeanandstd], 'XLim', [0 xlimtodisplay]);
        ylabel('Power')
        msgt = sprintf('Mean and Std Power spectra from FFT %s %s', eegpatient, eegcondition);
        title(msgt);
        
        if isempty(rois) == 0
            for ia=1:length(rois)
                curoi = rois{ia};
                figureareas(ia) = figure;
                hst = suptitle(['Mean Power spectra ROI =',curoi]);
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
                            powerx = powerx([channelclassindexes],:)
                            meanvx = mean(powerx);
                            stdx = std(powerx);
                            maxofmean =max(meanvx);
                            freqsindexl = powerfreqsindexes{ip,ic};
                            figure(figureareas(ia));
                            subplot(nbp,nbc,3*(ip-1) + ic);
                            % plot mean as such
                            plot(freqsindexl,meanvx);
                            set(gca, 'YLim', [0, 1.1*maxofmean], 'XLim', [0 xlimtodisplay]);
                            ylabel('Power ')
                            msgt = sprintf('Mean Power spectra ROI= %s, %s, %s', curoi, eegpatient, eegcondition);
                            title(msgt);
                        end
                    else
                        fprintf('Skipping Patient %s, she has not area %s\n',eegpatient,curoi);
                    end
                end
            end
        end
    end
end
end