function [myfullname, EEG, channel_labels, patdate, patsession] =  initialize_EEG_variables(patientid,patientcond)
%% initialize_EEG_variables  load mat file and initialize EEG lab variables .
%   [myfullname, EEG_study, channel_labels] = initialize_EEG_variables()e
%IN: patientid,patientcond
%Output: myfullname:full path of the mat file with the signal
%EEG: the struct created by EEGLAB, channel_labels:label of the channels,
%patientcond, patientid, patdate, patsession
% patientid= 'TWH024';
% patientcond = 'HYP';
if nargin < 2
    %if this function is called ith no condition we just assume one
    patientcond = 'HYP';
end
[mydir, myfile,patdate,patsession] = getfullpathfrompatient(patientid,patientcond);
myfullname = fullfile(mydir, myfile);
disp('Loading mat file...')
EEG = load(myfullname);
if isstruct(EEG) == 1
    if (strcmp(patientid, 'TWH033') + strcmp(patientid, 'TWH034') + strcmp(patientid, 'TWH024')) == 1
        EEG = EEG.EEG;
    elseif (strcmp(patientid,'TWH037')+ strcmp(patientid,'TWH038')+ strcmp(patientid,'TWH042')+ strcmp(patientid,'TWH043')+ strcmp(patientid,'TWH045')+ strcmp(patientid,'TWH049')+ strcmp(patientid,'TWH047')+ strcmp(patientid,'TWH048')) == 1
        EEG = EEG.EEGepocht10t20;
        if strcmp(patientid,'TWH042') == 1
            %42 has 48 channels
            EEG.nbchan = 48;
        end
    elseif (strcmp(patientid, 'TWH030') + strcmp(patientid, 'TWH027') + strcmp(patientid, 'TWH028')+ strcmp(patientid, 'TWH031')) == 1
        if isfield(EEG, 'EEG_cut_BL_HYP') EEG = EEG.EEG_cut_BL_HYP;
        elseif isfield(EEG, 'EEG_cut_BL_EC_PRE') EEG = EEG.EEG_cut_BL_EC_PRE;
        elseif isfield(EEG, 'EEG_cut_BL_EO_PRE') 
            EEG = EEG.EEG_cut_BL_EO_PRE;
        elseif isfield(EEG, 'EEG_cut_BL_EC_POST') EEG = EEG.EEG_cut_BL_EC_POST;
        end
    end
end
% correct wrong labeling in the EEG object compared with the electrodeNames
% files
EEG = fixwronglabeling(patientid, EEG);
disp(['File ' myfile ' loaded!' ])
channel_labels = {EEG.chanlocs.labels};
%initchan = 2; % channel 1 is the null Event channel
%disp([ 'Displaying the label of all the channels....' ])
%for i=initchan:EEG.nbchan
%    disp(['Index ' i ' is the channel' channel_labels(i)])
%end
end