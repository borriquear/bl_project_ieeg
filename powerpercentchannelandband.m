function []  = powerpercentchannelandband(patientid, patientcond,EEG,label, h)
%% powerpercentchannelandband relative power per band for one patient and one condition
%INPUT:  patientid, eegcond, EEG = 'HYP', 'ECPRE' ...
%E.g. (TWH028) patientid= 'TWH030'; eegcond = 'HYP', EEG
%OUTPUT: Figure with power relative per band and quantities in mat file 
%(appended at the end of the mat file calculated in fouriertransformanalysis)
% NEEDS fouriertransformanalysis function has created the fft_ mat file
if nargin ==5
    hold on;
    figure(h);
else
    h = figure;
    %msgtitle = sprintf('Relative P/Hz per bands, Cond=%s Patient=%s',patientcond,patientid);
end
freq_bands = {'\delta'; '\theta'; '\alpha'; '\beta'; '\gamma'};
num_bands = length(freq_bands);

listoffrqperband = zeros(1,[],num_bands);
percentlistoffrqperband  = zeros(1,[],num_bands);
% label = 'All'; %by default
% label = 'BiTemp'; %All, T, AT etc. in creatematfilewithchannelclassperpatient.mat
index_of_channels_of_interest = getindexesfromlabel(patientid, label);
fprintf('The list fo channel of interest is %s \n', label);
disp(index_of_channels_of_interest);

%conditiotostring = strcat('BL',patientcond);
[patdir, patfile,patdate,patsession] = getfullpathfrompatient(patientid,patientcond);
myfullname = fullfile(patdir, patfile);
mattoload = strcat('fft_',patientcond,'_', patientid,'_',patdate,'_',patsession,'.mat');
fprintf('Mat file to load and save percentlistoffrqperband is%s\n ',mattoload );

fftfile = fullfile(patdir,'figures', mattoload);
load(fftfile);
% load channel_labels after load file to overload the channel_labels
% variable of the file
channel_labels = {EEG.chanlocs.labels}; %channels +1 (Event);
channel_labels = channel_labels(index_of_channels_of_interest); %channels labels subset

%listoflabelsfound = [];
for bandidx=1:num_bands
    fprintf('Calculating for %s band relative power vector\n', freq_bands{bandidx})
    %load eeg to get chan_labels
    bitmpcounter = 0; %counts how many electrodes are labeled as 'label' (bi temp, All ...)
    %for chidx=2:EEG.nbchan %channel_labels2
    for chidx =1:length(index_of_channels_of_interest)
        indexlabel = index_of_channels_of_interest(chidx);
        listoffrqperband(1,chidx,bandidx) = requested_frequences_power_bnds(indexlabel-1,bandidx);
        %listoflabelsfound =[listoflabelsfound channel_labels(chidx)];
        bitmpcounter = bitmpcounter + 1;
        %indexlabel = getnameidx(total_channel_list_labels,channel_labels(chidx));
%         if indexlabel > 0
%             listoffrqperband(1,chidx-1,bandidx) = requested_frequences_power_bnds(indexlabel,bandidx);
%             bitmpcounter = bitmpcounter + 1;
%             listoflabelsfound =[listoflabelsfound channel_labels(chidx)];
%         else
%             fprintf('That is too bad, patient:%s has not typical bi temp config, we show what she has ch:%d\n',patientid, chidx);
%             listoffrqperband(1,chidx-1,bandidx) = requested_frequences_power_bnds(chidx-1,bandidx);
%         end
        
    end
    fprintf('Number of bi temp labeled channels is :%d, patient%s\n',bitmpcounter,patientid);
    %disp(listoflabelsfound)
end

%for chidx=2:EEG.nbchan
for chidx =1:length(index_of_channels_of_interest)
    for bandidx=1:num_bands
        sumparallbands =sum(listoffrqperband(1,chidx,1:end));
        elementvalue = listoffrqperband(1,chidx,bandidx);
        percactualband = elementvalue/sumparallbands;
        percentlistoffrqperband(1,chidx,bandidx) = percactualband;
        fprintf('Patient:%s Channel:%d, Band:%d, percentlistoffrqperband:%.2f\n',patientid,chidx,bandidx,percactualband)
    end
end

% adjust the xticklabel with the nb of channels
xchannel_limit = length(index_of_channels_of_interest);
xlabeljump = fix(xchannel_limit/6);
%xchannel_limit = EEG.nbchan-1;
%suptitle(msgtitle);
% - Build title axes and title.
hold on;
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
%text( 0.5, 0.3, 'Relative \mu V^2/Hz per band Patient=H','FontSize', 11', 'FontWeight', 'Bold', ...
%   'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' )
%text( 0.5, 0.3,patientcond)
hst = suptitle(['Relative \mu V^2/Hz per band. Patient=',patientid,' ,condition=' ,patientcond, ', ChannelClass=',label ])
%hst = suptitle('Relative \mu V^2/Hz per band Patient=H, all conditions')
set(hst,'Position', [0.5 0.0 0], 'FontSize',11,'FontWeight','normal','HorizontalAlignment' ,'center','VerticalAlignment', 'top')

hsubband = subplot(2,3,1);
%hsubband = findobj('Tag','delta');
plot(hsubband, percentlistoffrqperband(1,:,1));
title(freq_bands{1});
%set(gca, 'XLim', [1 xchannel_limit],'ylim', [0 1],'XTick',[1:xlabeljump:EEG.nbchan]);
set(gca, 'XLim', [1 xchannel_limit],'Ylim', [0 1],'XTick',[1:xchannel_limit],'XTickLabel',[channel_labels]);
xticklabel_rotate([],45,[],'Fontsize',6);
xlabel('Channels'), ylabel('\mu V^2/Hz');
grid on ;
hold on;
hsubband = subplot(2,3,2);
%hsubband = findobj('Tag','theta');
plot(hsubband,percentlistoffrqperband(1,:,2))

title(freq_bands{2})
set(gca, 'XLim', [1 xchannel_limit],'Ylim', [0 1],'XTick',[1:xchannel_limit],'XTickLabel',[channel_labels]);
xticklabel_rotate([],45,[],'Fontsize',6);xlabel('Channels'), ylabel('\mu V^2/Hz');
grid on ;
hold on;
hsubband = subplot(2,3,3);
%hsubband = findobj('Tag','alpha');
plot(hsubband,percentlistoffrqperband(1,:,3))

xlabel('Channels'), ylabel('\mu V^2/Hz')
set(gca, 'XLim', [1 xchannel_limit],'Ylim', [0 1],'XTick',[1:xchannel_limit],'XTickLabel',[channel_labels]);
xticklabel_rotate([],45,[],'Fontsize',6);title(freq_bands{3});
grid on ;
hold on;
hsubband = subplot(2,3,4)
% hsubband = findobj('Tag','beta');
plot(hsubband,percentlistoffrqperband(1,:,4))

xlabel('Channels'), ylabel('\mu V^2/Hz')
set(gca, 'XLim', [1 xchannel_limit],'Ylim', [0 1],'XTick',[1:xchannel_limit],'XTickLabel',[channel_labels]);
xticklabel_rotate([],45,[],'Fontsize',6);title(freq_bands{4});
grid on ;
hold on;
hsubband = subplot(2,3,5)
% hsubband = findobj('Tag','gamma');
plot(hsubband,percentlistoffrqperband(1,:,5));

xlabel('Channels'), ylabel('\mu V^2/Hz')
set(gca, 'XLim', [1 xchannel_limit],'Ylim', [0 1],'XTick',[1:xchannel_limit],'XTickLabel',[channel_labels]);
xticklabel_rotate([],45,[],'Fontsize',6);title(freq_bands{5});
grid on ;hold on;
savefigure(myfullname, h,1,'powerrelativeperband');
%
% make sure that globalFsDir is assigned
[globalFsDir] = loadglobalFsDir();
patpath = strcat(globalFsDir,patientid);
fftfile = fullfile(patpath,'data','figures', mattoload);
fprintf('Appending percentlistoffrqperband at file %s\n',fftfile);
save(fftfile,'percentlistoffrqperband', '-append')
end