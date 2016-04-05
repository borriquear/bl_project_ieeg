function [myfullname, EEG, channel_labels, patdate, patsession] =  initialize_EEG_variables(patientid,patientcond)
%% initialize_EEG_variables  load mat file and initialize EEG lab variables .
%   [myfullname, EEG_study, channel_labels] = initialize_EEG_variables()
%IN: patientid,patientcond
%Output: myfullname:full path of the mat file with the signal
%EEG: the struct created by EEGLAB, channel_labels:label of the channels,
%patientcond, patientid, patdate, patsession
% patientid= 'TWH024';
% patientcond = 'HYP';
[mydir, myfile,patdate,patsession] = getfullpathfrompatient(patientid,patientcond);
myfullname = fullfile(mydir, myfile);
disp('Loading mat file...')
EEG = load(myfullname);
if isstruct(EEG) == 1
    if (strcmp(mydir, 'D:\BIAL PROJECT\patients\TWH033\data') + strcmp(mydir, 'D:\BIAL PROJECT\patients\TWH034\data') + strcmp(mydir, 'D:\BIAL PROJECT\patients\TWH024\data')) == 1
        EEG = EEG.EEG
    elseif (strcmp(mydir, 'D:\BIAL PROJECT\patients\TWH030\data') + strcmp(mydir, 'D:\BIAL PROJECT\patients\TWH027\data') + strcmp(mydir, 'D:\BIAL PROJECT\patients\TWH028\data')+ strcmp(mydir, 'D:\BIAL PROJECT\patients\TWH031\data')) == 1
        EEG = EEG.EEG_cut_BL_HYP
        %EEG = EEG.EEG_cut_BL_EC_PRE
        %EEG = EEG.EEG_cut_BL_EO_PRE
        %EEG = EEG.EEG_cut_BL_EC_POST
    end
end
disp(['File ' myfile ' loaded!' ])
channel_labels = {EEG.chanlocs.labels};
disp([ 'Displaying the label of all the channels....' ])
initchan = 2 % channel 1 is the null Event channel
for i=initchan:EEG.nbchan
    disp(['Index ' i ' is the channel' channel_labels(i)])
end
end