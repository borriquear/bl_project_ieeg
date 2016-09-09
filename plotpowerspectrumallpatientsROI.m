function [] = plotpowerspectrumallpatientsROI(patientlist, conditionslist, powerspecmatrix, powerfreqsindexes, powerx_mean_perbands, rois )
%plotpowerspectrumallpatientsROI plot all patients mean of all channels in one figure
%figure1=figure('Position', [100, 100, 1024, 1200]);
if strcmpi(rois, 'All') == 1
    displayperRois = 0; 
    fprintf('Plotting FFT for all Areas \n');
else
    displayperRois = 1; 
    fprintf('Plotting FFT for all Areas and ROIs\n')
end

figm = figure;%powerx (plot)
figp = figure;   %powerxperband (bars)
% figstd = figure;
% suptitle(['Mean and std Power spectra All areas']);
nbp = length(patientlist);
nbc = length(conditionslist);
xlimtodisplay = 50;
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
        ystr = sprintf('Power %s',eegpatient);
        ylabel(ystr);
        msgt = sprintf('Mean Power spectra All ROIS %s', eegcondition);
        if ip == 1
            title(msgt, 'interpreter', 'none');
        end
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
        ystr2 = sprintf('Normalized Power %s', eegpatient);
        xlabel('Frequency Bands'), ylabel(ystr2);
        msgt = sprintf('Mean Power spectra  All ROIS %s', eegcondition);
        if ip == 1
            title(msgt, 'interpreter', 'none');
        end
        yel = mean(powerx_mean_perbands{ip,ic})/normfqb;
        ylistpat(ip,ic,:) = yel;
    end
end
%mean all patients per frequency band. from figure 2, one row and = number
%of conditions
figure;
ympb = zeros(nbc,nbfreqs);
for idxc =1:nbc
    for idxf=1:nbfreqs
        ympb(idxc, idxf)= mean(ylistpat(:,idxc,idxf))
    end
%     subplot(1, nbc,idxc);
%     bar(ympb(idxc,:));
%     set(gca, 'XTick', 1:nbfreqs, 'XTickLabel',freq_bands);
%     ylim([0 1]);
%     xlabel('Frequency Bands'), ylabel('Normalized Power per Band')
%     msgt = sprintf('Mean Power spectra All patients-ROIS %s',conditionslist{idxc});
%     title(msgt,'interpreter', 'none');
    %%%%
    bar_handle(idxc) = bar(ympb(idxc,:));
    set(bar_handle(idxc),'FaceColor',[mod(nbc-idxc+1,2),mod(nbc-idxc+1,2),mod(nbc-idxc+1,2)])
    hold on
    set(gca, 'XTick', 1:nbfreqs, 'XTickLabel',freq_bands);      
end
%show power per frequency bar in the same plot fro all conditions 
%ympb = (conditions, freq bands)
ylim([0 1]);
xlabel('Frequency Bands'), ylabel('Normalized Power per Band')
msgt = sprintf('Mean Power spectra All patients-ROIS');
%legend(conditionslist{1}, conditionslist{2}, 'Interpreter', 'None')
legend('EC PRE', 'EO PRE', 'Interpreter', 'None')
title(msgt,'interpreter', 'none');


%Show imagesc chart p[atientsxbands for pat andcond
figure;
for cc=1:length(conditionslist)
    subplot(1, length(conditionslist),cc);
    C = reshape(ylistpat(:,cc,:),[],length(freq_bands));
    new_C = C;%delete rows all zeros
    new_C(~any(C,2),:) = [];
    imagesc(new_C);
    set(gca,'XTick', [1:nbfreqs], 'YTickLabel',patientlist,'ytick', 1:length(patientlist));
    set(gca,'XTickLabel', freq_bands);
    ylabel('Patients'), xlabel('Frequency bands');
    condl = conditionslist(cc);
    msgt= sprintf('Power normalized ALL regions in %s' ,condl{1});
    title(msgt,'interpreter', 'none')
    colorbar;
    caxis([0 1]);
end
%Charts per ROIs

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
                        title(msgt,'interpreter', 'none');
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
                        title(msgt,'interpreter', 'none');
                        yel = mean(power_avg_bands)/normfqb;
                        ycondp(ip,ic,:) = yel;
                    end
                else
                    %close(figureareas(ia));
                    fprintf('Skipping Patient %s, she has not area %s\n',eegpatient,curoi);
                end
            end
            yrois{ia} = ycondp;
            
            
%         %%%%% plot relative power for ROIs    
%             figure;
%             ympb = zeros(nbc,nbfreqs);
%             for idxc =1:nbc
%                 for idxf=1:nbfreqs
%                     ympb(idxc, idxf)= mean(ylistpat(:,idxc,idxf))
%                 end
%                 subplot(1, nbc,idxc);
%                 bar(ympb(idxc,:));
%                 set(gca, 'XTick', 1:nbfreqs, 'XTickLabel',freq_bands);
%                 ylim([0 1]);
%                 xlabel('Frequency Bands'), ylabel('Normalized Power per Band')
%                 msgt = sprintf('Mean Power spectra All patients-ROIS %s',conditionslist{idxc});
%                 title(msgt,'interpreter', 'none');
%             end
            
        end
        
        for ir=1:length(rois)
            curoi = rois{ir};
            figure(figuresc(ir));
            actualmatrix = yrois{ir};
            for cc=1:length(conditionslist)
                subplot(1, length(conditionslist),cc);
                C = reshape(actualmatrix(:,cc,:),[],length(freq_bands));
                new_C = C;%delete rows all zeros
                new_C(~any(C,2),:) = [];
                imagesc(new_C);
                set(gca,'XTick', [1:nbfreqs], 'YTick',[]);
                set(gca,'XTickLabel', freq_bands);
                ylabel('Patients'), xlabel('Frequency bands'); cd= conditionslist(cc);
                msgt= sprintf('Power normalized in ROI=%s in %s',curoi ,cd{1});
                title(msgt,'interpreter', 'none')
                colorbar;
                caxis([0 1]);
            end
        end
        %plot mean per patient for ROIS
        figure;
        [np nc nf] = size(actualmatrix);
        condmeanrois = zeros(nc,nf);
        for cc=1:nc
            for cf=1:nf
                condmeanrois(cc,cf ) = mean(actualmatrix(:,cc,cf));
            end
        end 
        %plot mean per ROI
        for cc=1:nc
            bar_handle(cc) = bar(condmeanrois(cc,:));
            hold on
            set(bar_handle(cc),'FaceColor',[mod(nc-cc+1,2),mod(nc-cc+1,2),mod(nc-cc+1,2)])
            hold on
            set(gca, 'XTick', 1:nbfreqs, 'XTickLabel',freq_bands);     
        end
        ylim([0 1]);
        xlabel('Frequency Bands'), ylabel('Normalized Power per Band')
        msgt = sprintf('Mean Power spectra All patients-ROI=%s',curoi);
        %legend(conditionslist{1}, conditionslist{2}, 'Interpreter', 'None')
        legend('EC PRE', 'EO PRE', 'Interpreter', 'None')
        title(msgt,'interpreter', 'none');

    end
end %end displayperROIs
end