%% Load condition [power|other]data and display coloured electrodes
clear all
patientid = 'TWH028';
patientid = 'TWH030';
patientid = 'TWH037';

condition = 'ECPRE'
condition = 'EOPRE'
condition = 'ECPOST'
condition = 'HYP'

[myfullname, EEG, channel_labels] =  initialize_EEG_variables(patientid,condition);
%load data of electrodes to display
%%%%????? check how loads this!!!!
[quantitytomeasure, frqperband] = loadconditionpowerdata(myfullname,patientid,condition);
[celloflabelswithq] = getelectrodenamesfromfile(patientid, condition);
% make sure that globalFsDir exists
[globalFsDir] =loadglobalFsDir()
% Processing Electrodes with value coloured
tot_channels = size(celloflabelswithq{1},1);
elecNames=cell(tot_channels,1);
initchan =1;
listofelecnames = celloflabelswithq{1}
for i=initchan:tot_channels
    %chan2use = channel_labels(i)
    chan2use = listofelecnames(i)
    %%elecNames{i-1}=sprintf('%s',chan2use{1});
    elecNames{i}=sprintf('%s',chan2use{1});  
end
cfg=[];
cfg.view='l';
cfg.elecShape='sphere';
%cfg.elecColors=rand(8,1);
cfg.elecColors= quantitytomeasure;
%cfg.elecColors = celloflabelswithq{2};
cfg.elecColorScale='minmax';
cfg.showLabels='y';
cfg.elecUnits='r';
cfg.elecNames=elecNames;
cfg.elecSize=3;
%cfg.title='PT001: Stimulus Correlations';
cfg.title=sprintf('Power across freqs: %s, %s',patientid,condition)
cfg.title=sprintf('Power: %s, %s',patientid,condition)
cfgOut=plotPialSurf(patientid,cfg);
fprintf('Done');
%call savefigure
%% Processing Pair correlations coloured
pairs=[];
ct=0;
for a=1:7,
  ct=ct+1;
  pairs{ct,1}=sprintf('LGd%d',a+8);
  pairs{ct,2}=sprintf('LGd%d',a+1+8);
  pairs{ct,3}=rand(1,3); % RGB val
  pairs{ct,4}='L';
end
cfg=[];
cfg.view='l';
cfg.figId=2;
cfg.pairs=pairs;
cfg.showLabels='n';
cfg.elecUnits='r';
% if patient has no grid, abovelabels cant work
%patientid = 'PT001'
cfg.title=sprintf('%s: Pair Correlationse',patientid)
cfg_out=plotPialSurf(patientid,cfg);
%% Mapping Electrode Locations to Average Brains
% The surface based mapping is derived by trying to 
% align the gyrification patterns of the individual to that of the average 
% and has been shown to be superior to volume-based mapping 
% (i.e., warping the entire brain to that of an average brain) when doing 
% fMRI group analyses (Fischl et al., 1999).
% To use this Matlab code you first need to make sure you've downloaded 
% the FreeSurfer average brain and placed it in your FreeSurfer subject
%folder along with all your patient folders.

%collect the coordinates of each individual's brain relative to the average brain 
%patientid1 = 'TWH030';

patientid2 = 'TWH036';
groupAvgCoords=[];
groupLabels=[];
groupIsLeft=[];
%subs={patientid1,patientid2};
subs={patientid2};
for a=1:length(subs),
    
    fprintf('Working on Participant %s\n',subs{a});
    fprintf('Showing Average on the Left and patient: %s on the Right\n',subs{a});
    [avgCoords, elecNames, isLeft]=pial2AvgBrain(subs{a},[]);
    groupAvgCoords=[groupAvgCoords; avgCoords];
    groupLabels=[groupLabels; elecNames];
    groupIsLeft=[groupIsLeft; isLeft];
end
cfg=[];
cfg.view='l';
cfg.elecCoord=[groupAvgCoords groupIsLeft];
cfg.elecNames=groupLabels;
cfg.showLabels='n';
cfg.title=sprintf('%s: on Avg. Brain',patientid2)
cfgOut=plotPialSurf('fsaverage',cfg);

%% Mapping Electrodes to Atlases
% The Desikan-Killiany Atlas (Desikan et al., 2006) is a 35 area cortical atlas that is based on gyral morphology.
%DK Atlas here: http://episurg.pbworks.com/w/file/104506753/desikanKillianyAreas.txt
patientid = 'TWH033'
fprintf('Showing Pial morphology according to DK Atlas (35 areas) for: %s\n',patientid);
cfg=[];
cfg.view='l';
cfg.overlayParcellation='DK';
cfg.title=sprintf('%s:  DK Atlas (35 areas)',patientid)
cfgOut=plotPialSurf(patientid,cfg);
%To get the Desikan-Killiany area that is closest to each subdural electrode
patientid = 'TWH036'
parcOut=elec2Parc(patientid,'DK');

% Destrieux Atlas is also a gyral morphology based atlas
% it divides the cortical surface into 75 areas. 
%To view an individual's pial surface color coded according to the atlas
% The atlas here: http://episurg.pbworks.com/w/file/104599657/destrieuxTable.pdf
cfg=[];
cfg.view='l';
cfg.overlayParcellation='D';
cfg.title=sprintf('%s:  Destrieux Atlas (75 areas)',patientid)
cfgOut=plotPialSurf(patientid,cfg);
% To get the Destrieux area that is closest to each subdural electrode
%parcOut=elec2Parc(patientid,'D');
% Intrinic Functional Connectvity Atalas
% https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation_Yeo2011
% Yeo et al. resting-state fMRI Atlas
createIndivYeoMapping('TWH028');
cfg=[];
cfg.view='l';
cfg.overlayParcellation='Y17';
cfg.title='PT001: Yeo 17-Area'; 
cfgOut=plotPialSurf('TWH036',cfg);
%% Viewing the anatomical location of electrodes
%Destrieux Atlas of electrodes
cfg=[];
cfg.view='omni';
cfg.overlayParcellation='D';
cfg.showLabels='y';
cfg.title=sprintf('%s:Destrieux Atlas',patientid)
%'TWH033: Destrieux Atlas';
cfgOut=plotPialSurf(patientid,cfg);

%% View the location of depth electrodes with semi-transparent pial surface
cfg=[];
cfg.view='li';
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.5;
%cfg.onlyShow={'Dp1','Dp2','Dp3','Dp4','Dp5','Dp6','Dp7','Dp8','Da1','Da2','Da3','Da4','Da5','Da6','Da7','Da8'};
cfg.onlyShow={'rpt2','tp1'};
%cfg.title='TWH033';
cfg.title=sprintf('%s:semi-transparent pial surface',patientid)
cfgOut=plotPialSurf(patientid,cfg);