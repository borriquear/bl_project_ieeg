function [electrodesincommon] = calculateintersectionchannels(patientsL)
% calculateintersectionchannels  given a list of patients returns the
% electrodes in common
%   electrodesincommon = calculateintersectionchannels(A) e.g. A={'TWH01',....,'TWH11' }
%   electrodesincommon = {'LHD1', ..., 'RHD4'}
%
%patientsL = { 'TWH027','TWH024','TWH030', 'TWH031','TWH033','TWH034'};
%all patients empty set, no 28 LHD1-4
electrodesincommon = {};
disp('Loading globalFsDir...')
if ~exist('globalFsDir','var')
    fprintf('globalFsDir not found, loading it...')
    eval('global globalFsDir');
    myp = 'D:\BIAL PROJECT\patients\'
    eval(['globalFsDir=' 'myp']);
end

for indpat=1:length(patientsL)
        currentpat = patientsL{indpat};
        %fprintf('The current patient is:%s\n',currentpat);
        patpath = strcat(globalFsDir,currentpat);
        fprintf('The current path is %s\n',patpath);
        fprintf('Loading the EEG object for patient  %s\n',currentpat);
        [myfullname, EEG, channel_labels, patdate, patsession] = initialize_EEG_variables(currentpat);
        if size(electrodesincommon,1) == 0
            electrodesincommon = {channel_labels{2:end}};
        else
            electrodesincommon = intersect(electrodesincommon, {channel_labels{2:end}} );
        end
        %disp('Elctrodes in common so far:');
        %disp(electrodesincommon);
end
end
