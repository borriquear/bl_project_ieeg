function []  = anovatesthypnosis()
%% anovatesthypnosis relative power per band for one patient and one condition
% at a time
%INPUT: patidx= 1...5, eegcond = 'HYP', 'ECPRE' ...
%E.g. (TWH028) patidx = 2; eegcond = 'HYP' anovatesthypnosis(patidx, eegcond)
% Labels concident with TMPBMP
%       patient 24 = 36/ 36  pero si nos fijamos en la labesl (4/36) 
%       patient 27 = 36/36
%       patient 28 = 0/88
%       patient 30 = 34/36
%       patient 31 = 36/104
%       patient 33 = 36/82
%       patient 34 = 16/108
total_channel_list_labels = {'LHD1','LHD2','LHD3', 'LHD4', 'LAT1','LAT2'...
    ,'LAT3','LAT4', 'LMT1','LMT2', 'LMT3','LMT4','LPT1' ,'LPT2','LPT3','LPT4','LPT5','LPT6','RHD1'...
    ,'RHD2','RHD3','RHD4','RAT1','RAT2','RAT3','RAT4','RMT1','RMT2','RMT3', 'RMT4','RPT1','RPT2'...
    ,'RPT3','RPT4','RPT5','RPT6'};
freq_bands = {'\delta'; '\theta'; '\alpha'; '\beta'; '\gamma'};
num_bands = 5;

listoffrqperband = zeros(1,[],num_bands);
percentlistoffrqperband  = zeros(1,[],num_bands);

patientid= 'TWH024';
patientcond = 'HYP';
%patidx= 1;
conditiotostring = strcat('BL',patientcond);
[patdir, patfile,patdate,patsession,patname] = getfullpathfrompatient(patientid,patientcond);
myfullname = fullfile(patdir, patfile);
EEG = load(myfullname);
if isfield(EEG, 'EEG_cut_BL_HYP')
    EEG = EEG.EEG_cut_BL_HYP;
elseif isfield(EEG,'EEG')
    EEG = EEG.EEG;
end
mattoload = strcat('fft_',conditiotostring,'_', patname,'_',patdate,'_',patsession,'.mat');

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
            fprintf('Thats too bad, patient:%s has not typical bi temp config, we show what she has ch:%d\n',patientid, chidx);
            listoffrqperband(1,chidx-1,bandidx) = requested_frequences_power_bnds(chidx-1,bandidx);
        end
    end
    fprintf('Number of bitemp labeled channels is :%d, patient%s\n',bitmpcounter,patientid);
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

msgtitle = sprintf('Relative P/Hz per bands, Cond=%s Patient=%s',patientcond,patientid);
% adjust the xticklabel with the nb of channels
switch patientid
    case {'TWH027','TWH030' }
        xlabeljump=5; %36 channels
    case {'TWH028','TWH033','TWH034', 'TWH031'}
        xlabeljump=15; %88 channels
    otherwise
        xlabeljump=5; %36 channels
end
h = figure;
%suptitle(msgtitle);
% - Build title axes and title.
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
text( 0.5, 0.3, msgtitle, 'FontSize', 12', 'FontWeight', 'Bold', ...
    'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' )
xchannel_limit = EEG.nbchan;
subplot(2,3,1)
plot(percentlistoffrqperband(1,:,1))
title(freq_bands{1});
set(gca, 'XLim', [1 xchannel_limit],'ylim', [0 1],'XTick',[1:xlabeljump:EEG.nbchan]);
xlabel('Channels'), ylabel('Power/Hz');
grid on ;
subplot(2,3,2)
plot(percentlistoffrqperband(1,:,2))
title(freq_bands{2})
set(gca, 'xlim', [1 xchannel_limit],'ylim', [0 1],'XTick',[1:xlabeljump:EEG.nbchan]);
xlabel('Channels'), ylabel('Power/Hz');
grid on ;
subplot(2,3,3)
plot(percentlistoffrqperband(1,:,3))
xlabel('Channels'), ylabel('Power/Hz')
set(gca, 'xlim', [1 xchannel_limit],'ylim', [0 1],'XTick',[1:xlabeljump:EEG.nbchan]);
title(freq_bands{3});
grid on ;
subplot(2,3,4)
plot(percentlistoffrqperband(1,:,4))
xlabel('Channels'), ylabel('Power/Hz')
set(gca, 'xlim', [1 xchannel_limit],'ylim', [0 0.5],'XTick',[1:xlabeljump:EEG.nbchan] );
title(freq_bands{4});
grid on ;
subplot(2,3,5)
plot(percentlistoffrqperband(1,:,5))
xlabel('Channels'), ylabel('Power/Hz')
set(gca, 'xlim', [1 xchannel_limit],'ylim', [0 0.5],'XTick',[1:xlabeljump:EEG.nbchan] );
title(freq_bands{5});
grid on ;
savefigure(myfullname, h,1,'powerrelativeperband')
end