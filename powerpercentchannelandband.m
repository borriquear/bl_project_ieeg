function []  = powerpercentchannelandband(patientid, patientcond,EEG,h)
%% powerpercentchannelandband relative power per band for one patient and one condition
%INPUT:  patientid, eegcond, EEG = 'HYP', 'ECPRE' ...
%E.g. (TWH028) patientid= 'TWH030'; eegcond = 'HYP', EEG
%OUTPUT: Figure with power relative per band and quantities in mat file 
%(appended at the end of the mat file calculated in fouriertransformanalysis)
% NEEDS fouriertransformanalysis function has created the fft_ mat file
if nargin == 4
    hold on;
    figure(h);
else
    h = figure;
    %msgtitle = sprintf('Relative P/Hz per bands, Cond=%s Patient=%s',patientcond,patientid);
end
%msgtitle = text('Relative \mu V^2/Hz per band, Cond=%s Patient=%s',patientcond,patientid);
%title('Relative \mu V^2/Hz per band Patient=H, condtion=HYP') 
brainregion = 'temphd';
total_channel_list_labels = getlistofchannelsofinterest(brainregion)
freq_bands = {'\delta'; '\theta'; '\alpha'; '\beta'; '\gamma'};
num_bands = 5;

listoffrqperband = zeros(1,[],num_bands);
percentlistoffrqperband  = zeros(1,[],num_bands);

%patidx= 1;
conditiotostring = strcat('BL',patientcond);
[patdir, patfile,patdate,patsession] = getfullpathfrompatient(patientid,patientcond);
myfullname = fullfile(patdir, patfile);
mattoload = strcat('fft_',patientcond,'_', patientid,'_',patdate,'_',patsession,'.mat');
channel_labels = {EEG.chanlocs.labels}; %chan2use = channel_labels(irow);
total_chan = size(channel_labels,2);

fftfile = fullfile(patdir,'figures', mattoload);
load(fftfile);
listoflabelsfound = [];
for bandidx=1:num_bands
    fprintf('Calculating for %s band relative power vector\n', freq_bands{bandidx})
    %load eeg to get chan_labels
    bitmpcounter = 0; %counts how many electrodes arelabeled as bi temp
    for chidx=2:EEG.nbchan 
        indexlabel = getnameidx(total_channel_list_labels,channel_labels(chidx));
        if indexlabel > 0
            listoffrqperband(1,chidx-1,bandidx) = requested_frequences_power_bnds(indexlabel,bandidx);
            bitmpcounter = bitmpcounter + 1;
            listoflabelsfound =[listoflabelsfound channel_labels(chidx)];
        else
            fprintf('That is too bad, patient:%s has not typical bi temp config, we show what she has ch:%d\n',patientid, chidx);
            listoffrqperband(1,chidx-1,bandidx) = requested_frequences_power_bnds(chidx-1,bandidx);
        end
    end
    fprintf('Number of bi temp labeled channels is :%d, patient%s\n',bitmpcounter,patientid);
    disp(listoflabelsfound)
end

for chidx=2:EEG.nbchan
    for bandidx=1:5
        sumparallbands =sum(listoffrqperband(1,chidx-1,1:end));
        elementvalue = listoffrqperband(1,chidx-1,bandidx);
        percactualband = elementvalue/sumparallbands;
        percentlistoffrqperband(1,chidx-1,bandidx) = percactualband;
        fprintf('Patient:%s Band:%d, Channel:%d is percentlistoffrqperband:%.2f\n',patientid,bandidx,chidx-1,percactualband)
    end
end

% adjust the xticklabel with the nb of channels
switch patientid
    case {'TWH024','TWH027','TWH030', 'TWH038','TWH042'}
        xlabeljump=5; %36 channels
    case {'TWH028','TWH033','TWH034', 'TWH031','TWH037','TWH043'}
        xlabeljump=15; %88 channels
    otherwise
        xlabeljump=5; %36 channels
end
xchannel_limit = EEG.nbchan-1;
%suptitle(msgtitle);
% - Build title axes and title.
hold on;
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
%text( 0.5, 0.3, 'Relative \mu V^2/Hz per band Patient=H','FontSize', 11', 'FontWeight', 'Bold', ...
%   'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' )
%text( 0.5, 0.3,patientcond)
hst = suptitle(['Relative \mu V^2/Hz per band Patient=H, condition= ' patientcond ])
%hst = suptitle('Relative \mu V^2/Hz per band Patient=H, all conditions')
set(hst,'Position', [0.5 0.0 0], 'FontSize',11,'FontWeight','normal','HorizontalAlignment' ,'center','VerticalAlignment', 'top')

% tf = ishold;
% if tf < 1 
%     hold on;
% end
% hsub1 = subplot(2,3,1);set(hsub1,'Tag','delta');
% hsub2 = subplot(2,3,2);set(hsub2,'Tag','theta');
% hsub3 = subplot(2,3,3);set(hsub3,'Tag','alpha');
% hsub4 = subplot(2,3,4);set(hsub4,'Tag','beta');
% hsub5 = subplot(2,3,5);set(hsub5,'Tag','gamma');


hsubband = subplot(2,3,1);
%hsubband = findobj('Tag','delta');
plot(hsubband, percentlistoffrqperband(1,:,1));
title(freq_bands{1});
%set(gca, 'XLim', [1 xchannel_limit],'ylim', [0 1],'XTick',[1:xlabeljump:EEG.nbchan]);
set(gca, 'XLim', [1 xchannel_limit],'Ylim', [0 1],'XTick',[1:xchannel_limit],'XTickLabel',[total_channel_list_labels]);
xticklabel_rotate([],45,[],'Fontsize',6);
xlabel('Channels'), ylabel('\mu V^2/Hz');
grid on ;
hold on;
hsubband = subplot(2,3,2);
%hsubband = findobj('Tag','theta');
plot(hsubband,percentlistoffrqperband(1,:,2))

title(freq_bands{2})
set(gca, 'XLim', [1 xchannel_limit],'Ylim', [0 1],'XTick',[1:xchannel_limit],'XTickLabel',[total_channel_list_labels]);
xticklabel_rotate([],45,[],'Fontsize',6);xlabel('Channels'), ylabel('\mu V^2/Hz');
grid on ;
hold on;
hsubband = subplot(2,3,3);
%hsubband = findobj('Tag','alpha');
plot(hsubband,percentlistoffrqperband(1,:,3))

xlabel('Channels'), ylabel('\mu V^2/Hz')
set(gca, 'XLim', [1 xchannel_limit],'Ylim', [0 1],'XTick',[1:xchannel_limit],'XTickLabel',[total_channel_list_labels]);
xticklabel_rotate([],45,[],'Fontsize',6);title(freq_bands{3});
grid on ;
hold on;
hsubband = subplot(2,3,4)
% hsubband = findobj('Tag','beta');
plot(hsubband,percentlistoffrqperband(1,:,4))

xlabel('Channels'), ylabel('\mu V^2/Hz')
set(gca, 'XLim', [1 xchannel_limit],'Ylim', [0 1],'XTick',[1:xchannel_limit],'XTickLabel',[total_channel_list_labels]);
xticklabel_rotate([],45,[],'Fontsize',6);title(freq_bands{4});
grid on ;
hold on;
hsubband = subplot(2,3,5)
% hsubband = findobj('Tag','gamma');
plot(hsubband,percentlistoffrqperband(1,:,5));

xlabel('Channels'), ylabel('\mu V^2/Hz')
set(gca, 'XLim', [1 xchannel_limit],'Ylim', [0 1],'XTick',[1:xchannel_limit],'XTickLabel',[total_channel_list_labels]);
xticklabel_rotate([],45,[],'Fontsize',6);title(freq_bands{5});
grid on ;hold on;
savefigure(myfullname, h,1,'powerrelativeperband');
%
% make sure that globalFsDir is assigned
if ~exist('globalFsDir','var') 
   fprintf('globalFsDir not found, loading it...')
   eval('global globalFsDir');
   myp = 'D:\BIAL PROJECT\patients\';
   eval(['globalFsDir=' 'myp']); 
end
patpath = strcat(globalFsDir,patientid);
mattoload = strcat('fft_',patientcond,'_', patientid,'_',patdate,'_',patsession,'.mat');
fftfile = fullfile(patpath,'data','figures', mattoload);
fprintf('Appending percentlistoffrqperband at file %s\n',fftfile);
save(fftfile,'percentlistoffrqperband', '-append')
end