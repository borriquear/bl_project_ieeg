function [] = powerpercentchannelandband_all(list_of_patientid, list_of_patientcond)
%% powerpercentchannelandband_all compares power-frequency across conditions per same patient or accross patients per one condition
%IN: list_of_patientid = {pat1} list_of_patientcond={cond1, cond2,...condN}
% or list_of_patientid = {pat1...patN} list_of_patientcond={cond1}
% example 1: powerpercentchannelandband_all({'TWH030'}, {'EC_PRE', 'HYP', 'EC_POST'})
% example 2: powerpercentchannelandband_all({'TWH030', 'TWH027', 'TWH024'}, {'HYP'})

hold on;
sizeconds = size(list_of_patientcond);
sizeconds = sizeconds(end);
sizepats = size(list_of_patientid)
sizepats = sizepats(end);
labelchannelclass ='All';
if sizepats == 1 
    %for one patient all conditions
    vectorstoplot = sizeconds; %if 1 patient we plot the conditions
    patientid = list_of_patientid{1};
    for i=1:sizeconds
        hold on;
        h = figure;
        patientcond = list_of_patientcond{i};
        [myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(patientid,patientcond);
        powerpercentchannelandband(patientid, patientcond,EEG,labelchannelclass,h); 
        %load the vector
        conditions{i}= loadconditionsvector(patientid,patientcond, EEG);
    end
else
    if size(list_of_patientcond) == 1 
        % compare different patients for the same condition
        vectorstoplot = sizepats; %if 1 condition we plot all the patients
        patientcond = list_of_patientcond{1};
        for i=1:sizepats
            h = figure;
            patientid = list_of_patientid{i};
            [myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(patientid,patientcond);
            powerpercentchannelandband(patientid, patientcond,EEG,labelchannelclass,h );
            %load the vector
            conditions{i}= loadconditionsvector(patientid,patientcond, EEG);
        end   
    else
        masgerror = sprintf('Violated Condition: 1 patient, N conditions OR N patients,1 condition pats=%s conds=%s', num2str(size(list_of_patientid)), num2str(size(list_of_patientcond)));
        error(masgerror);
    end
end
% plot conditions vector
%%%%%%%%%%%%%%%%%
freq_bands = {'\delta'; '\theta'; '\alpha'; '\beta'; '\gamma'};
index_of_channels_of_interest = getindexesfromlabel(patientid, labelchannelclass);
xchannel_limit = length(index_of_channels_of_interest);
channel_labels = {EEG.chanlocs.labels}; %channels +1 (Event);
channel_labels = channel_labels(index_of_channels_of_interest); 

hsubband = subplot(2,3,1);
%hsubband = findobj('Tag','delta');
for icond=1:vectorstoplot
    plot(hsubband, conditions{icond}.percentlistoffrqperband(:,:,1));
    hold on;
end
title(freq_bands{1});
%set(gca, 'XLim', [1 xchannel_limit],'ylim', [0 1],'XTick',[1:xlabeljump:EEG.nbchan]);
set(gca, 'XLim', [1 xchannel_limit],'Ylim', [0 1],'XTick',[1:xchannel_limit],'XTickLabel',[channel_labels]);
xticklabel_rotate([],45,[],'Fontsize',6);
xlabel('Channels'), ylabel(' \mu V^2/Hz ');
grid on ;
hsubband = subplot(2,3,2);
%hsubband = findobj('Tag','theta');
for icond=1:vectorstoplot
    plot(hsubband, conditions{icond}.percentlistoffrqperband(:,:,2));
    hold on;
end
title(freq_bands{2})
set(gca, 'XLim', [1 xchannel_limit],'Ylim', [0 1],'XTick',[1:xchannel_limit],'XTickLabel',[channel_labels]);
xticklabel_rotate([],45,[],'Fontsize',6);xlabel('Channels'), ylabel(' \mu V^2/Hz ');
grid on ;
hold on;
hsubband = subplot(2,3,3);
%hsubband = findobj('Tag','alpha');
for icond=1:vectorstoplot
    plot(hsubband, conditions{icond}.percentlistoffrqperband(:,:,3));
    hold on;
end
xlabel('Channels'), ylabel(' \mu V^2/Hz')
set(gca, 'XLim', [1 xchannel_limit],'Ylim', [0 1],'XTick',[1:xchannel_limit],'XTickLabel',[channel_labels]);
xticklabel_rotate([],45,[],'Fontsize',6);title(freq_bands{3});
grid on ;
hsubband = subplot(2,3,4)
% hsubband = findobj('Tag','beta');
for icond=1:vectorstoplot
    plot(hsubband, conditions{icond}.percentlistoffrqperband(:,:,4));
    hold on;
end

xlabel('Channels'), ylabel(' \mu V^2/Hz ')
set(gca, 'XLim', [1 xchannel_limit],'Ylim', [0 1],'XTick',[1:xchannel_limit],'XTickLabel',[channel_labels]);
xticklabel_rotate([],45,[],'Fontsize',6);title(freq_bands{4});
grid on ;
hsubband = subplot(2,3,5)
% hsubband = findobj('Tag','gamma');
for icond=1:vectorstoplot
    plot(hsubband, conditions{icond}.percentlistoffrqperband(:,:,5));
    hold on;
end
if sizepats  == 1
    legend(hsubband, list_of_patientcond);
    hst = suptitle(['Patient:', patientid, ', Relative \mu V^2/Hz per band in ', labelchannelclass])
elseif sizeconds ==1
    legend(hsubband, list_of_patientid);
    hst = suptitle(['Relative \mu V^2/Hz per band in Channels=',labelchannelclass, ', Condition=',patientcond])
end

xlabel('Channels'), ylabel(' \mu V^2/Hz ')
set(gca, 'XLim', [1 xchannel_limit],'Ylim', [0 1],'XTick',[1:xchannel_limit],'XTickLabel',[channel_labels]);
xticklabel_rotate([],45,[],'Fontsize',6);title(freq_bands{5});
grid on ;hold on;
%suptitle('Relative \mu V^2/Hz per band Patient=H, condition=all')
set(hst,'Position', [0.5 0 0], 'FontSize',11,'FontWeight','normal','HorizontalAlignment' ,'center','VerticalAlignment', 'top')
end
function conditions = loadconditionsvector(patientid,patientcond, EEG)
        [patdir, patfile,patdate,patsession] = getfullpathfrompatient(patientid,patientcond);
        myfullname = fullfile(patdir, patfile);
        mattoload = strcat('fft_',patientcond,'_', patientid,'_',patdate,'_',patsession,'.mat');
        channel_labels = {EEG.chanlocs.labels}; %chan2use = channel_labels(irow);
        total_chan = size(channel_labels,2);
        fftfile = fullfile(patdir,'figures', mattoload);
        conditions = load(fftfile,'percentlistoffrqperband');
end