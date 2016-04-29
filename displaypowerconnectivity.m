
function displaypowerconnectivity()
% Display power connectivity analysis from the correlation matrix
% calculated in powerbasedconnectivityall.mat
% The file "powerconnectivity_freq_[10]_HYP_TWH[034_02092016_s2].mat" must
% exist
eegcondition = 'HYP';
centerfrequencies = {2, 6 , 10, 23.5, 40};
eegpatientl = { 'TWH027','TWH024','TWH028','TWH030', 'TWH031','TWH033','TWH034'};
allmatrixpowcon = cell(length(eegpatientl),length(centerfrequencies));
channels_pat = [];
for indpat=1:length(eegpatientl)
    [myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(eegpatientl{indpat},eegcondition);
    channels_pat{indpat} = channel_labels;
    fprintf('Calculating corr. matrix for patient %s', eegpatientl{indpat});
    for indexfreq = 1:length(centerfrequencies) % delta, theta, alpha, beta, gamma
        centerfreq = centerfrequencies{indexfreq};
        initelem = zeros(length(channel_labels)-1,length(channel_labels)-1);
        allmatrixpowcon{indpat,indexfreq} = initelem;
        allmatrixpowcon{indpat,indexfreq} = opencorrelationmatrix( eegpatientl{indpat}, eegcondition,centerfreq,eegdate, eegsession);
    end
    drawnow limitrate nocallbacks
end
%display charts from allmatrixpowcon
for indpat=1:length(eegpatientl)
    [myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(eegpatientl{indpat},eegcondition);
    h = figure;
    for j=1:length(centerfrequencies)
        %[matcorrtodisp]= allmatrixpowcon{indpat,j};
        displaymeancorrelationperlectrode(allmatrixpowcon{indpat,j},eegpatientl{indpat},eegcondition, centerfrequencies{j},h,channels_pat{indpat} );
    end
end
end

function [corrmatpersubband] = opencorrelationmatrix( eegpatient, eegcondition,centerfreq,eegdate, eegsession )
%opencorrelationmatrix( eegpatient, eegcondition,centerfreq,eegdate,
%eegsession ) Open the correlation matrix
%IN:patient
%OUT:correlation matrix for the specified input
fprintf('Open mat file with the correlation matrix for patient %s at frequency %s in %s\n',eegpatient, num2str(centerfreq),eegcondition)
if ~exist('globalFsDir','var')
    fprintf('globalFsDir not found, loading it...')
    eval('global globalFsDir');
    myp = 'D:\BIAL PROJECT\patients\'
    eval(['globalFsDir=' 'myp']);
end
patpath = strcat(globalFsDir,eegpatient);
mattoload = strcat('powerconnectivity_freq_',num2str(centerfreq),'_',eegcondition,'_', eegpatient,'_',eegdate,'_',eegsession,'.mat');
fftfile = fullfile(patpath,'data','figures', mattoload);
fprintf('Opening correlation matrix....%s\n',fftfile);
matf= matfile(fftfile);
corrmatpersubband = matf.corr_matrix;
end

function  displaymeancorrelationperlectrode(corrmat, patientid, cond, freq, h, channels_pat)
%displaycorrelationmatrix(corrmat) display mean correlation per electrode
%find how many bands to display
figure(h);
fprintf('Displaying correlation matrix per patient %s at frequency %s in %s\n',patientid, num2str(freq),cond)
msgtitle = sprintf('Power correlation per Frequency Band, Patient=%s, Cond=%s',patientid,cond);
switch patientid
    case {'TWH027','TWH030' }
        xlabeljump=5; %36 channels
    case {'TWH028','TWH033','TWH034', 'TWH031'}
        xlabeljump=15; %88 channels
    otherwise
        xlabeljump=5; %36 channels
end
% adjust the xticklabel with the nb of channels
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0.3, msgtitle, 'FontSize', 12', 'FontWeight', 'Bold', ...
    'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' )
xchannel_limit = length(channels_pat)-1;
switch freq %2, 6 , 10, 23.5, 40};
    case 2
        freq_band = '\delta';
        indexchart = 1;
    case 6
        freq_band = '\theta';
        indexchart = 2;
    case 10
        freq_band = '\alpha';
        indexchart = 3;
    case 23.5
        freq_band = '\beta';
        indexchart = 4;
    case 40
        freq_band = '\gamma';
        indexchart = 5;
    otherwise
        fprintf('ERROR frequency band missing!!');
        return;
end
%calculate the mean per row or col
[fil,col] = size(corrmat);
powerperelec = zeros(fil);
for indm = 1:fil
    powerperelec(indm) = mean(corrmat(indm,:));
end
subplot(2,3,indexchart)
bar(powerperelec)
title(freq_band); %max(max(powerperelec))+0.4
set(gca, 'XLim', [1 xchannel_limit],'ylim', [0 1],'XTick',[1:xlabeljump:length(channels_pat)-1]);
xticklabel_rotate([],45,[],'Fontsize',6);
%set(gca, 'XLim', [1 xchannel_limit],'ylim', [0 1],'XTick',[1:xlabeljump:length(channels_pat)-1]);
xlabel('Channels'), ylabel('Power-based Spearman corr. ');
grid on ;
end