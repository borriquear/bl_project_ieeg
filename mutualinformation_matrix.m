function  [mi_matrix] = mutualinformation_matrix()
%mutualinformation_matrix returns the MI matrix a symmetric matrix with the
%MI for each pair of channels. The matrix is saved as a mat file in the patient/data/figures folder
%The patients and conditions is specified in patientlist and conditionslist
% Call this once to create the mat file, compitationally EXPENSIVE!
global globalFsDir;
[globalFsDir] = loadglobalFsDir();
patientlist = {'TWH030','TWH031','TWH034','TWH033','TWH038','TWH042'};%,'TWH043','TWH037'
%Patient 37 and 43  have NOT HD,  only Deep
conditionslist = {'EC_PRE', 'EO_PRE','HYP','EC_POST','EO_POST'};
indexp = 1;
indexc = 1;
isinlist = 0; % isinlist = 1 To calculate MI ALL Pairs
label1 = 'All'; label2 = 'All';
if isinlist < 1 % isinlist = 0 To calculate among subgroups of channels
    label1 = 'HD'; label2 = 'NOHD';
end
meanMImatrix = zeros(length(patientlist), length(conditionslist));
for ip=1:length(patientlist)
    eegpatient = patientlist{ip};
    label1indexes = getindexesfromlabel(eegpatient, label1);
    label2indexes = getindexesfromlabel(eegpatient, label2);
    for ic=1:length(conditionslist)
        eegcond = conditionslist{ic};
        %patient 34 missing EO PRE
        if strcmp(eegpatient, 'TWH034') +  strcmp(eegcond, 'EO_PRE') > 1
            continue
        end
        if strcmp(eegpatient, 'TWH030') +  strcmp(eegcond, 'EO_POST') > 1
            continue
        end
        if strcmp(eegpatient, 'TWH031') +  strcmp(eegcond, 'EO_POST') > 1
            continue
        end
        [myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(eegpatient,eegcond);
        patdir = fullfile(globalFsDir,eegpatient, 'data\figures');
        fprintf('Saving the mat file with the MI matrix in %s \n',patdir);
        filematname = sprintf('mutualinformation_%s_%s_%s_%s_%s_%s.mat',label1,label2,eegcond, eegpatient,eegdate,eegsession);
        %get all the combination of pairs
        channel_labels= {EEG.chanlocs.labels};
        channel_pairs = combntns(channel_labels(2:end),2);
        total_pairs= length(channel_pairs);
        matrix_MI = zeros(length(channel_labels)-1);
        %total_pairs = 3;
        for i=1:total_pairs
            electrodes4mi = channel_pairs(i,:);
            electrodesidx(1) = find(strcmpi(electrodes4mi{1},{EEG.chanlocs.labels}));
            electrodesidx(2) = find(strcmpi(electrodes4mi{2},{EEG.chanlocs.labels}));
            % mi is non positive and symmetric, thus we do the triangular superior
            % isinlist = 1 to calculatethe MI for all pair
            isinlist = ismember(electrodesidx(1),label1indexes) && ismember(electrodesidx(2),label2indexes);
            if (electrodesidx(1) <= electrodesidx(2)) && (isinlist)
                fprintf('%s, %s : Calculating MI between channels %s and %s\n',eegpatient, eegcond, electrodes4mi{1} , electrodes4mi{2});
                [mi, hx, hy, hxy, nb] = mutualinformation2channels(electrodes4mi, EEG);
                channel_entropy(electrodesidx(1)) = hx; channel_entropy(electrodesidx(2)) = hy;
                fprintf(' MI = %s , H(x)=%s H(y)=%s H(xy)=%s Nbins= %s\n',num2str(mi) , num2str(hx),num2str(hy),num2str(hxy),num2str(nb));
                matrix_MI(electrodesidx(1)-1,electrodesidx(2)-1) = mi;
                matrix_MI(electrodesidx(2)-1,electrodesidx(1)-1) = matrix_MI(electrodesidx(1)-1,electrodesidx(2)-1);
            end
            
        end
        %disp(size(matrix_MI))
        fprintf('Printing out the MI matrix\n')
        meanMI = sum(sum(matrix_MI))/nnz(matrix_MI);
        fprintf('The mean MI of the metrix is:%.2f\n', meanMI);
        meanMImatrix(ip,ic) = meanMI;
        %save matrix in mat file
        filematname = fullfile(patdir,filematname);
        fprintf('Saving the MI matrix in %s\n',filematname);
        save(filematname, 'eegpatient','eegcond','channel_labels', 'matrix_MI', 'channel_entropy');
    end
end
disp('The mean MI matrix values patientsxconditions :')
disp(meanMImatrix)
disp('Calling displayMIcharts to display MIcharts')
displayMIcharts()
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
times2save =  timewindow:stepsize:length(EEG.times)- 1000;

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
mi      = zeros(1,length(times2save)); %zeros(2,...
nbins   = zeros(1,length(times2save));
for timei = 1:length(times2save)
    try
        datax = EEG.data(electrodesidx(1),times2saveidx(timei)-timewindowidx:times2saveidx(timei)+timewindowidx);
        datay = EEG.data(electrodesidx(2),times2saveidx(timei)-timewindowidx:times2saveidx(timei)+timewindowidx);
        [mi(1,timei),entropy(:,timei),nbins(timei)] = mutualinformationx(datax,datay); %variable number of bins
        %[mi(2,timei),entropy(:,timei)             ] = mutualinformationx(datax,datay,20); %fixed umber of bins
    catch err
        if strcmp(err.identifier, 'MATLAB:badsubscript')
            fprintf('CAUTION: (Line 108 mutualinformation2channels) exceed EEG.data size\n')
            continue
        else
            rethrow(err); %some other unexpected error. Better stop
        end
    end
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

function [] = displayMIcharts()
%%displayMIcharts display the MI among labels .
% Requirement: mutualinformation_%s_%s_%s_%s_%s_%s.mat already created
global globalFsDir;
[globalFsDir] = loadglobalFsDir();
patientlist = {'TWH030','TWH031','TWH033','TWH034','TWH038','TWH042'}%,'TWH043'};%,'TWH037'
%Patient 37 and 43  have NOT HD,  only Deep
conditionslist = {'EC_PRE', 'EO_PRE','HYP','EC_POST','EO_POST'};
indexp = 1;
indexc = 1;
MIValues = zeros(size(patientlist,2), size(conditionslist,2),2);
figuresl = zeros(length(patientlist));
maxhd = max(max(MIValues(:,:,1)));
maxhdno = max(max(MIValues(:,:,2)));
maxt = max(maxhd,maxhdno);
for ip=1:length(patientlist)
    eegpatient = patientlist{ip};
    for ic=1:length(conditionslist)
        eegcond = conditionslist{ic};
        [isvalid, hashd] = ifisvalidpatcondcombination(eegpatient,eegcond);
        legend_text{ic}= eegcond;
        if (isvalid + hashd < 2)
            fprintf('CAUTION: Omitting this patient %s condition %s',eegpatient, eegcond)
            continue
        end
        [myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(eegpatient,eegcond);
        %HD - HD
        label1 = 'HD'; label2 = label1;
        patdir = fullfile(globalFsDir,eegpatient, 'data\figures');
        fprintf('Opening the the mat file with the MI matrix HD-HD %s \n',patdir);
        filematnamehdhd = sprintf('mutualinformation_%s_%s_%s_%s_%s_%s.mat',label1,label2, eegcond, eegpatient,eegdate,eegsession);
        filematnamehdhd = fullfile(patdir,filematnamehdhd);
        fhd2 = load(filematnamehdhd);
        meanMI = sum(sum(fhd2.matrix_MI))/nnz(fhd2.matrix_MI);
        MIValues(ip,ic,1) = meanMI ;
        %HD NO HD
        fprintf('Opening the the mat file with the MI matrix HD-NOHD %s \n',patdir);
        label2 = 'NOHD';
        filematnamehdno = sprintf('mutualinformation_%s_%s_%s_%s_%s_%s.mat',label1,label2, eegcond, eegpatient,eegdate,eegsession);
        filematnamehdno = fullfile(patdir,filematnamehdno);
        fhdno = load(filematnamehdno);
        meanMI = sum(sum(fhdno.matrix_MI))/nnz(fhdno.matrix_MI);
        MIValues(ip,ic,2) = meanMI ;
    end   
end
% plot MIValues in figuresl figures
for ip=1:length(patientlist)
   eegpatient = patientlist{ip};
   legend_patients{ip}= eegpatient;
    if strcmp(eegpatient,'TWH037') + strcmp(eegpatient,'TWH043') > 0
        continue
    end
    [legend_text, veccondinexes] = legendtextforpatient(eegpatient);
    difference_hdinout(ip,[veccondinexes]) =  MIValues(ip,[veccondinexes],1) - MIValues(ip,[veccondinexes],2);  %difference between MI(HD,HD) - MI(HD,NOHD)
    figuresl(ip) = figure;
    figure(figuresl(ip));
    plot(MIValues(ip,[veccondinexes],1));
    hold on
    plot(MIValues(ip,[veccondinexes],2));
    set(gca, 'YLim', [0, 2],'xticklabel', legend_text,'XTick', [1:size(veccondinexes,2)])
    ylabel('MI (bits)')
    msgt = sprintf('MI(HD,HD) and MI(HD,NOHD) electrodes in patient %s',eegpatient);
    title(msgt);
    legend('HD-HD', 'HD-NOHD');
end
% plot MIValues all encompased in one figure
fh_totdif = figure; % bar MI(HDHD) - MI(HDNOHD) per condition and patient
bar(difference_hdinout);set(gca, 'YLim', [0, 1.5*max(max(difference_hdinout))],'xticklabel', legend_patients,'XTick', [1:size(legend_patients,2)])
xlabel('Patients'),ylabel('MI(HD,HD)- MI(HD,NOHD) (bits)')
msgt = sprintf('MI(HD,HD)- MI(HD,NOHD) all patients');
title(msgt);
legend(legend_text);

%plot both MI(HD,HD) and MI(HD,NOHD) for eachpatient
fh_totpats = figure; % bar MI(HDHD) - MI(HDNOHD) per condition and patient
%set(gcf,'Color',[1,0.4,0.6]) % bakground color
plot(MIValues(:,:,1)');
hold all
plot(MIValues(:,:,2)', '-.');
set(gca, 'YLim', [0, 1.2*max(MIValues(:))],'xticklabel', legend_text,'XTick', [1:size(legend_text,2)]);
set(gca, 'ColorOrderIndex', 1)
xlabel('Patients'),ylabel('MI(HD,HD)- MI(HD,NOHD) (bits)')
msgt = sprintf('MI(HD,HD)- MI(HD,NOHD) all patients');
title(msgt);
legend(legend_patients);
end
