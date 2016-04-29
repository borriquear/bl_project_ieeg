function []   = powerperspecificelectrodes(patientHigh, patientLow)
%%powerperspecificelectrodes(patientHigh, patientLow)
% IN:patientid, patientcondition
% OUT: display charts for the common channels between the two groups
fprintf('Verify whether the electrodes are present in all subjects\n');
freq_bands = {'\delta'; '\theta'; '\alpha'; '\beta'; '\gamma'};
patientcond = 'HYP'
% patientHigh = {'TWH030','TWH028'};
% patientLow = {'TWH034', 'TWH027','TWH033','TWH024'};
%common channels for a list of patients
[commonchannelslistH] = commonchannelsfromlistofpatients(patientHigh,patientcond);
[commonchannelslistL] = commonchannelsfromlistofpatients(patientLow,patientcond);
listofChannels = intersect(commonchannelslistH,commonchannelslistL);
tot_channels = size(listofChannels,2);
% the channels in common are
fprintf('The channels that patients have in common is\n:')
disp(listofChannels)
tot_patientsH = size(patientHigh,2);
tot_patientsL = size(patientLow,2);
Highmatrices = zeros(tot_patientsH,tot_channels,5);
Lowmatrices = zeros(tot_patientsL,tot_channels,5);
meanHighmatrices=zeros(tot_patientsH,tot_channels,5);
meanLowmatrices= zeros(tot_patientsH,tot_channels,5);
vecofchvaluesL = zeros(tot_channels,5);
vecofchvaluesH = zeros(tot_channels,5);
for i=1:size(patientHigh,2)
    patientid = patientHigh(i);
    patientHvalscurrent = getelectrodenamesfromfile(patientid{1}, patientcond);
    %labelspat = {};valspat={};bandaspat={};
    labelspat = patientHvalscurrent{1};
    valspat =  patientHvalscurrent{2};
    bandaspat =  patientHvalscurrent{3};
    % find the channel index that corresponds with labelspat
    for c =1:tot_channels
        labelelectrode = listofChannels(c);
        found = sum(strcmp({labelspat{:}},labelelectrode));
        if found == 1
            indexh = find(strcmp({labelspat{:}},labelelectrode),1);
            fprintf('Found label:%s in index:%d\n',labelelectrode{1},indexh);
        else
            fprintf('Label %s NOT FOND!!!\n',labelelectrode{1});
        end
        for j=1:5
            Highmatrices(i,c,j)= bandaspat(1,indexh,j);
        end
    end
end
% low patients
for i=1:size(patientLow,2)
    patientid = patientLow(i);
    patientHvalscurrent = getelectrodenamesfromfile(patientid{1}, patientcond);
    %labelspat = {};valspat={};bandaspat={};
    labelspat = patientHvalscurrent{1};
    valspat =  patientHvalscurrent{2};
    bandaspat =  patientHvalscurrent{3};
    % find the channel index that corresponds with labelspat
    for c =1:tot_channels
        labelelectrode = listofChannels(c);
        found = sum(strcmp({labelspat{:}},labelelectrode));
        if found == 1
            indexh = find(strcmp({labelspat{:}},labelelectrode),1);
            fprintf('Found label:%s in index:%d\n',labelelectrode{1},indexh);
        else
            fprintf('Label %s NOT FOND!!!\n',labelelectrode{1});
        end
        for j=1:5
            Lowmatrices(i,c,j)= bandaspat(1,indexh,j);
        end
    end
end
%build matrices for mean
for i=1:tot_channels
    for j=1:5
        meanHighmatrices(1,i,j) = sum(Highmatrices(:,i,j))/tot_patientsH;
        meanLowmatrices(1,i,j) = sum(Lowmatrices(:,i,j))/tot_patientsL;
    end
end
%plot chart one per each band for low and high values
for i=1:5
    ch_low{i} = meanLowmatrices(1,:,i); % 
    canalbandaL =  ch_low{i};
    ch_high{i} = meanHighmatrices(1,:,i);
    canalbandaH =  ch_high{i};
    for j=1:tot_channels
        vecofchvaluesL(j,i) = canalbandaL(j);
        vecofchvaluesH(j,i) = canalbandaH(j);
    end
end
fprintf('Printing the chart with the average values for Highs and Lows...\n');
    
h = figure;
msgtitle = sprintf('Relative P/Hz per bands, Highs vs Lows \n');
%suptitle(msgtitle);
% - Build title axes and title.
%axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
xchannel_limit = tot_channels;
xlabeljump = xchannel_limit/4;
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, -0.5, msgtitle, 'FontSize', 10', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
subplot(2,3,1)
plot(vecofchvaluesH(:,1))
hold on
plot(vecofchvaluesL(:,1))
title(freq_bands{1});
set(gca, 'XLim', [1 xchannel_limit],'ylim', [0 1],'XTickLabel',listofChannels,'XTick',1:xchannel_limit);
xticklabel_rotate([],45,[],'Fontsize',6);
xlabel('Channels'), ylabel('Power/Hz');
legend('High hypnotizability patients','Low hypnotizability patients')
grid on ;

subplot(2,3,2)
plot(vecofchvaluesH(:,2))
hold on
plot(vecofchvaluesL(:,2))

title(freq_bands{2});
set(gca, 'XLim', [1 xchannel_limit],'ylim', [0 1],'XTickLabel',listofChannels,'XTick',1:xchannel_limit);
xticklabel_rotate([],45,[],'Fontsize',6);
xlabel('Channels'), ylabel('Power/Hz');
grid on ;

subplot(2,3,3)
plot(vecofchvaluesH(:,3))
hold on
plot(vecofchvaluesL(:,3))
title(freq_bands{3});
set(gca, 'XLim', [1 xchannel_limit],'ylim', [0 1],'XTickLabel',listofChannels,'XTick',1:xchannel_limit);
xticklabel_rotate([],45,[],'Fontsize',6);
xlabel('Channels'), ylabel('Power/Hz');
grid on ;
subplot(2,3,4)
plot(vecofchvaluesH(:,4))
hold on
plot(vecofchvaluesL(:,4))
title(freq_bands{4});
set(gca, 'XLim', [1 xchannel_limit],'ylim', [0 0.5],'XTickLabel',listofChannels,'XTick',1:xchannel_limit); 
xticklabel_rotate([],45,[],'Fontsize',6);
xlabel('Channels'), ylabel('Power/Hz');
grid on ;
subplot(2,3,5)
plot(vecofchvaluesH(:,5))
hold on
plot(vecofchvaluesL(:,5))
title(freq_bands{5});
set(gca, 'XLim', [1 xchannel_limit],'ylim', [0 0.25],'XTickLabel',listofChannels,'XTick',1:xchannel_limit); 
xticklabel_rotate([],45,[],'Fontsize',6);
xlabel('Channels'), ylabel('Power/Hz');
grid on ;

% make sure that globalFsDir is assigned
if ~exist('globalFsDir','var') 
   fprintf('globalFsDir not found, loading it...')
   eval('global globalFsDir');
   myp = 'D:\BIAL PROJECT\patients\';
   eval(['globalFsDir=' 'myp']); 
end
%patpath = strcat(globalFsDir,patientid);
% patpath = globalFsDir;
% figtosave = strcat('power_relative_Highs_vs_Lows',patientcond,'_','.mat');
% fftfile = fullfile(patpath,'data','figures', figtosave);
% savefig(h,fftfile);

end
function [commonchannelslist] = commonchannelsfromlistofpatients(listpatients,condition)
tot_channels = size(listpatients,2);
%commonchannelslist = {};
commonchannelslist = [];

%channel_labels = {};
channel_labels = [];
for i=1:tot_channels
    patientid = listpatients{i};
    [myfullname, EEG, channel_labels] =  initialize_EEG_variables(patientid,condition);
    %delete 'Event' channel
    channel_labels = channel_labels(2:end);
    %commonchannelslist = {commonchannelslist channel_labels};
    commonchannelslist = [commonchannelslist channel_labels];
    %channel_labels = {};
end
end
