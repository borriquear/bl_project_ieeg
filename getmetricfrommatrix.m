function  [mean_value] = getmetricfrommatrix()
%getmetricfrommatrix returns the mean value for a matrix contained in the matfile
%The patients and conditions is specified in
%patientlist and conditionslist
global globalFsDir;
[globalFsDir] = loadglobalFsDir();
patientlist = {'TWH030','TWH031','TWH034','TWH033','TWH038','TWH042'};
% NOTE: 37 and 43 DO NOT HAVE HD
%34 do not have eo pre
conditionslist = {'EC_PRE', 'EO_PRE','HYP','EC_POST','EO_POST'};
variablename = matrix_MI; % power-based correlation matrix 
label1= 'HD';
label2 = 'NOHD';
indexp = 1;
indexc = 1;
for ip=indexp:indexp %length(patientlist)
    for ic=indexc:indexc
        eegpatient = patientlist{ip};
        eegcond = conditionslist{ic};
        fprintf('%s, %s, Calculating the mean MI for the Block Corr Matrix defined as %s and %s',eegpatient,eegcond ,label1 ,label2 )
        classlabel1 =  getindexesfromlabel(eegpatient, label1);
        classlabel2 =  getindexesfromlabel(eegpatient, label2);
        [myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(eegpatient,eegcond);
        patdir = fullfile(globalFsDir,eegpatient, 'data\figures');
        fprintf('Opening the mat file with the MI matrix in %s \n',patdir);
        filematname = sprintf('mutualinformation_%s_%s_%s_%s.mat',eegcond, eegpatient,eegdate,eegsession);
        % to open another file just start here
        fh = load(filematname);
        
        matrixmi = fh.variablename;
        blockmatrixmi = matrixmi(classlabel1, classlabel2);
        meanblockmatrixmi = mean2(blockmatrixmi);
        fprintf('The mean MI between %s and %s is:%.2f \n',label1,label2,meanblockmatrixmi);
    end
end

end

function  [mean_mi,mean_entropy_x,mean_entropy_y,mean_entropy_xy, mean_nbins] = mutualinformation2channels(electrodes4mi, EEG)
%mutualinformation mutual information between two channels
%sliding segment

% eegpatient = 'TWH033'
% eegcond = 'HYP'
% [myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(eegpatient,eegcond);
%electrodes4mi = {'RFP1','LHD2'}
timewindow = 400; % in ms
stepsize = 100; %ms
epochtime = length(EEG.data);
times2save = 1:stepsize:epochtime;
times2save =  timewindow:stepsize:length(EEG.times) -1000;

% convert ms to indices
timewindowidx = round(timewindow/(1000/EEG.srate)/2);
times2saveidx = zeros(size(times2save));
for i=1:length(times2save)
    [junk,times2saveidx(i)] = min(abs(EEG.times-times2save(i)));
end
electrodesidx(1) = find(strcmpi(electrodes4mi{1},{EEG.chanlocs.labels}));
electrodesidx(2) = find(strcmpi(electrodes4mi{2},{EEG.chanlocs.labels}));
% initialize outputs
entropy = zeros(3,length(times2save));
mi      = zeros(2,length(times2save));
nbins   = zeros(1,length(times2save));
for timei = 1:length(times2save)
    datax = EEG.data(electrodesidx(1),times2saveidx(timei)-timewindowidx:times2saveidx(timei)+timewindowidx);
    datay = EEG.data(electrodesidx(2),times2saveidx(timei)-timewindowidx:times2saveidx(timei)+timewindowidx);
    
    [mi(1,timei),entropy(:,timei),nbins(timei)] = mutualinformationx(datax,datay); %variable number of bins
    %[mi(2,timei),entropy(:,timei)             ] = mutualinformationx(datax,datay,20); %fixed umber of bins
end
mean_mi = mean(mi(1,:)); %mean mutual information in bits
mean_entropy_x = mean(entropy(1,:)); %mean entropy channel x
mean_entropy_y = mean(entropy(2,:)); %mean entropy channel y
mean_entropy_xy = mean(entropy(3,:)); %mean joint entropy channels x,y
mean_nbins = mean(nbins); % number of bins used for discretization
%                  (based on Freedman-Diaconis rule)
plotfigures = 0;
if plotfigures > 0
    figure
    set(gcf,'name',[ 'Mutual information between ' electrodes4mi{1} ' and ' electrodes4mi{2} ])
    
    subplot(221)
    plot(times2save,mi(1,:))
    xlabel('Time (ms)'), ylabel('MI (bits)')
    title('Variable bin length')
    set(gca,'xlim',[times2save(1)-50 times2save(end)+50],'ylim',[min(mi(:))-.01 max(mi(:))+.01])
    
    subplot(222)
    plot(nbins,mi(1,:),'.')
    xlabel('bin length'), ylabel('MI (bits)')
    title('Bin length vs. MI')
    set(gca,'ylim',[min(mi(:))-.01 max(mi(:))+.01])
    
    subplot(223)
    plot(times2save,mi(2,:))
    xlabel('Time (ms)'), ylabel('MI (bits)')
    title('Constant bin length')
    set(gca,'xlim',[times2save(1)-50 times2save(end)+50],'ylim',[min(mi(:))-.01 max(mi(:))+.01])
end
end
%%
