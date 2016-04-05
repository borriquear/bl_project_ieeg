function [patdir, patfile,patdate,patsession] = getfullpathfrompatient(patientid, patientcond, entireorcut)
%% [patdir, patfile] = getfullpathfrompatient(patientid, patientcond)
% returns the path and the file name of the mat file for that patient and
% condition
%Input: patientid 'TWH030' patientcond 'HYP', entireorcut optional if
%exists file we look for is the entire session
%Output: patdir, patfile

patdir1 = 'D:\BIAL PROJECT\patients\';
patdir = strcat(patientid, '\data');
patdir = fullfile(patdir1, patdir);
patsession = 's1'; %by default session is s1
switch patientid
    case 'TWH024'
        if strcmp(patientcond,'HYP') == 1
            %patname = 'fo24';
            patdate = '09192015';
        else
        end
    case 'TWH027'
        if strcmp(patientcond,'HYP') == 1
            %patname = 'bs27';
            patdate = '10222015';
            patsession = 's2';
        else
        end
    case 'TWH028'
        if strcmp(patientcond,'HYP') == 1
            %patname = 'cj28';
            patdate = '10202015';
        else
        end
    case 'TWH030'
        if strcmp(patientcond,'HYP') == 1
            %patname = 'was30';
            patdate = '11172015';
        else
        end
    case 'TWH031'
        if strcmp(patientcond,'HYP') == 1
            %patname = 'sm31';
            patdate = '12012015';
        else
        end
    case 'TWH033'
        if strcmp(patientcond,'HYP') == 1
            %patname = 'nk33';
            patdate = '02032016';
        else
        end
    case 'TWH034'
        if strcmp(patientcond,'HYP') == 1
            %patname = 'mj34';
            patdate = '02092016';
            patsession ='s2';
        else
        end
    otherwise
        warningMessage = sprintf('Error: patient%s does not exist:\n%s', patientid);
        uiwait(msgbox(warningMessage));
end
if nargin > 2
    patfile1 = 'EEG_entire_session_BL_';
else
    patfile1 = 'EEG_cut_BL_';
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