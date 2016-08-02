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
nbfreqs = length(freq_bands);
% % rois = {'HD','T','F', 'IH'};
%Charts ALL ROIs
ylistpat = zeros(length(patientlist),length(conditionslist),nbfreqs );

for ip=1:length(patientlist)
    eegpatient = patientlist{ip};
    for ic =1:length(conditionslist)
        yel = [];
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
        yel = mean(powerx_mean_perbands{ip,ic})/normfqb;
        ylistpat(ip,ic,:) = yel;
    end
end
%Show imagesc chart p[atientsxbands for pat andcond
figure;
for cc=1:length(conditionslist)
    subplot(1, length(conditionslist),cc);
    C = reshape(ylistpat(:,cc,:),[],length(freq_bands));
    new_C = C;%delete rows all zeros
    new_C(~any(C,2),:) = [];
    imagesc(new_C);
    set(gca,'XTick', [1:nbfreqs], 'YTick',[]);
    set(gca,'XTickLabel', freq_bands);
    ylabel('Patients'), xlabel('Frequency bands');
    condl = conditionslist(cc);
    msgt= sprintf('Power normalized ALL regions in %s' ,condl{1});
    title(msgt)
    colorbar;
    caxis([0 1]);
end
%Charts per ROIs
displayperRois = 1; %0
if displayperRois > 0
    ycondp = zeros(length(patientlist),length(conditionslist),nbfreqs );
    yrois = {};
    if isempty(rois) == 0
        for ia=1:length(rois)
            ycondp = zeros(length(patientlist),length(conditionslist),nbfreqs );
            figureareas(ia) = figure;
            figureareas_bars(ia) = figure;
            figuresc(ia) = figure;
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
                        yel = [];
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
                        powepp = powerx_mean_perbands{ip,ic};
                        power_avg_bands = powepp(channelclassindexes, :);
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
                        yel = mean(power_avg_bands)/normfqb;
                        ycondp(ip,ic,:) = yel;
                        
                        
                    end
                else
                    %close(figureareas(ia));
                    fprintf('Skipping Patient %s, she has not area %s\n',eegpatient,curoi);
                end
                
            end
            yrois{ia} = ycondp;
        end
        for ir=1:length(rois)
            curoi = rois{ir};
            figure(figuresc(ir));
            for cc=1:length(conditionslist)
                subplot(1, length(conditionslist),cc);
                actualmatrix = yrois{ir}
                C = reshape(actualmatrix(:,cc,:),[],length(freq_bands));
                new_C = C;%delete rows all zeros
                new_C(~any(C,2),:) = [];
                imagesc(new_C);
                set(gca,'XTick', [1:nbfreqs], 'YTick',[]);
                set(gca,'XTickLabel', freq_bands);
                ylabel('Patients'), xlabel('Frequency bands'); cd= conditionslist(cc);
                msgt= sprintf('Power normalized in ROI=%s in %s',curoi ,cd{1});
                title(msgt)
                colorbar;
                caxis([0 1]);
            end
        end
    end
end %end displayperROIs
end