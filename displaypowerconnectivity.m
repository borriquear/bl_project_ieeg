function displaypowerconnectivity(eegpatientl, conditionslistl, centerfrequencies)
% displaypowerconnectivity shows the correlation matix calculated in powerbasedconnectivityall.mat
% The file "powerconnectivity_freq_[10]_HYP_TWH[034_02092016_s2].mat" must
% exist (Connectivity is always about correlation)

allmatrixpowcon = cell(length(eegpatientl),length(centerfrequencies));
channels_pat = [];
%load allmatrixpowconmatrix
for indpat=1:length(eegpatientl)
    curpatient = eegpatientl{indpat};
    for indcond=1:length(conditionslistl)
        curcondition = conditionslistl{indcond};
        [myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(curpatient,curcondition);
        channels_pat{indpat} = channel_labels;
        fprintf('Calculating corr. matrix for patient %s', eegpatientl{indpat});
        for indexfreq = 1:length(centerfrequencies) % delta, theta, alpha, beta, gamma
            centerfreq = centerfrequencies(indexfreq);
            initelem = zeros(length(channel_labels)-1,length(channel_labels)-1);
            allmatrixpowcon{indpat,indcond, indexfreq} = initelem;
            allmatrixpowcon{indpat,indcond, indexfreq} = opencorrelationmatrix(curpatient, curcondition,centerfreq,eegdate, eegsession);
        end
        %drawnow limitrate nocallbacks
    end
end
%display charts from allmatrixpowcon

for indpat=1:length(eegpatientl)
    curpatient = eegpatientl{indpat};
    for indcond=1:length(conditionslistl)
        curcondition = conditionslistl{indcond};
        [myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(curpatient,curcondition);
        for j=1:length(centerfrequencies)
            %display correlation matrix
            %one figure for each freq band for corrplot
           h(j) = figure;
           displaycorrelationmatrix(allmatrixpowcon{indpat,indcond, j},curpatient,curcondition, centerfrequencies(j),h(j),channels_pat{indpat});
           % generate the network
           [allmetrics] = graphtheoryanalysis(allmatrixpowcon{indpat,indcond, j}, curpatient,curcondition, centerfrequencies(j), channels_pat{indpat});
        end
   end
end
end

function [corrmatpersubband] = opencorrelationmatrix( eegpatient, eegcondition,centerfreq,eegdate, eegsession )
%opencorrelationmatrix( eegpatient, eegcondition,centerfreq,eegdate,
%eegsession ) Open the correlation matrix

fprintf('Open mat file with the correlation matrix for patient %s at frequency %s in %s\n',eegpatient, num2str(centerfreq),eegcondition)
[globalFsDir] = loadglobalFsDir();
patpath = strcat(globalFsDir,eegpatient);
mattoload = strcat('powerconnectivity_freq_',num2str(centerfreq),'_',eegcondition,'_', eegpatient,'_',eegdate,'_',eegsession,'.mat');
fftfile = fullfile(patpath,'data','figures', mattoload);
fprintf('Opening correlation matrix....%s\n',fftfile);
matf= matfile(fftfile);
corrmatpersubband = matf.corr_matrix;
end

function [] = displaycorrelationmatrix(corr_mat,curpatient,curcondition,curfreq,hfig,channels_patl);
figure(hfig);
imagesc(corr_mat);
colormap('jet');
colorbar;
msgtitle = sprintf('Correlation matrix Spearman power %s, %s, frq=%d', curpatient,curcondition,curfreq);
title(msgtitle);
fprintf('Patient %s, Cond %s, Freq %d. The mean in || of the correlation amtrix= %.3f \n',curpatient,curcondition,curfreq,  mean(abs(corr_mat(:))))
end



% function  displaymeancorrelationperlectrode(corrmat, patientid, cond, freq, h, channels_pat, nrows)
% %displaycorrelationmatrix(corrmat) display BARS mean correlation per electrode
% %find how many bands to display
% figure(h);
% fprintf('Displaying correlation matrix per patient %s at frequency %s in %s\n',patientid, num2str(freq),cond);
% %patientid = 'L'
% msgtitle = sprintf('Sum of power correlation per channel, Patient=%s, Cond=%s',patientid,cond);
% switch patientid
%     case {'TWH027','TWH030' }
%         xlabeljump=5; %36 channels
%     case {'TWH028','TWH033','TWH034', 'TWH031', 'TWH037'}
%         xlabeljump=15; %88 channels
%     otherwise
%         xlabeljump=5; %36 channels
% end
% % adjust the xticklabel with the nb of channels
% axes( 'Position', [0, 0.95, 1, 0.05] ) ;
% set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
% text( 0.5, 0.3, msgtitle, 'FontSize', 12', 'FontWeight', 'Bold', ...
%     'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' )
% xchannel_limit = length(channels_pat)-1;
% switch freq %2, 6 , 10, 23.5, 40};
%     case 2
%         freq_band = '\delta';
%         indexchart = 1;
%     case 6
%         freq_band = '\theta';
%         indexchart = 2;
%     case 10
%         freq_band = '\alpha';
%         indexchart = 3;
%     case 23.5
%         freq_band = '\beta';
%         indexchart = 4;
%     case 40
%         freq_band = '\gamma';
%         indexchart = 5;
%     otherwise
%         fprintf('ERROR frequency band missing!!');
%         return;
% end
% %calculate the mean per row or col
% [fil,col] = size(corrmat);
% powerperelec = zeros(fil);
% for indm = 1:fil
%     powerperelec(indm) = mean(corrmat(indm,:));
% end
% if nrows == 1
%     %to display only one band
%     subplot(1,1,indexchart)
% else
%     subplot(2,3,indexchart)
% end
% bar(powerperelec);
% % patientid = 'H';
% % msgtitleband = sprintf('corrplot, Patient=%s, Cond=%s,Band=%',patientid,cond,freq_band);
% % title(msgtitleband);
% set(gca, 'XLim', [1 xchannel_limit],'YLim',[0 1],'XTick',[1:xchannel_limit],'XTickLabel',[channels_pat(2:end)]);
% xticklabel_rotate([],45,[],'Fontsize',6);
% xlabel('Channels'), ylabel('Power-based Spearman corr.');
% grid on ;
% end