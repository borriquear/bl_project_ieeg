function [patdir, patfile,patdate,patsession] = getfullpathfrompatient(patientid, patientcond, entireorcut)
%% [patdir, patfile] = getfullpathfrompatient(patientid, patientcond)
% returns the path and the file name of the mat file for that patient and
% condition.
%This function assigns the session and date for the patient and condition
%The mat file needs to start with EEG_cut_BL_ or EEG_entire_session_BL_
%Input: patientid 'TWH030' patientcond 'HYP', entireorcut optional if
%exists file we look for is the entire session
%Output: patdir, patfile,patdate,patsession

%patdir1 = 'D:\BIAL PROJECT\patients\';
patdir1 =  loadglobalFsDir();
patdir = strcat(patientid, '\data');
patdir = fullfile(patdir1, patdir);
patsession = 's1'; %by default session is s1
switch patientid
    case 'TWH024'
        patdate = '09192015';
    case 'TWH027'
        patdate = '10222015';
        patsession = 's2';
    case 'TWH028'
        %patname = 'cj28';
        patdate = '10202015';
    case 'TWH030'
        %patname = 'was30';
        patdate = '11172015';
    case 'TWH031'
        %patname = 'sm31';
        patdate = '12012015';
    case 'TWH033'
        patdate = '02032016';
    case 'TWH034'
        patdate = '02092016';
        patsession ='s2';
    case 'TWH037'
        patdate = '03142016';
        patsession ='s1';
    case 'TWH038'
        patdate = '03082016';
        patsession ='s1';
    case 'TWH042'
        patdate = '05042016';
        patsession ='s1';
    case 'TWH043'
        patdate = '05042016';
        patsession ='s1';
    case 'TWH047'
        patdate = '06282016';
        patsession ='s1';
    case 'TWH048'
        patdate = '06292016';
        patsession ='s1';
    case 'TWH045'
        patdate = '06022016';
        patsession ='s1';
    case 'TWH049'
        patdate = '06292016';
        patsession ='s1';
    otherwise
        warningMessage = sprintf('Error: patient %s does not exist:\n%s', patientid);
        uiwait(msgbox(warningMessage));
end
patfile1 = 'EEG_cut_';
if nargin > 2
    patfile1 = 'EEG_entire_session_BL_';
else
    if (  strcmp(patientid, 'TWH030') + strcmp(patientid, 'TWH031')+  strcmp(patientid, 'TWH033') + strcmp(patientid, 'TWH034')) == 1
        patfile1 =  'EEG_cut_BL_';
    end
end
patfile2 = patientcond;
patfile3 = strcat('_', patientid, '_', patdate,'_',patsession, '.mat');
patfile = strcat(patfile1,patfile2,patfile3);

fullFileName =  fullfile(patdir, patfile);

if exist(fullFileName, 'file')
    fprintf('Cool! I have found the mat file:%s, now I will load the EEG struct...',fullFileName);
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end
end