%% Preprocesssing iEEG data from Natus acquisition system
% Author: Jaime Gomez Ramirez, Ph.D., mathematical neuroscience
% The Hospital for Sick Chldren, Toronto, ON
% email address: jd.gomezramirez@gmail.com  
% Website: aslab.org/~gomez
% November 2015; Last revision: 16-February-2015

function [] = cutoneepoch()
%% cutoneepoch generates one epoch segment from a mat file of the entire session 
% DEPENDENCIES: It needs the .mat file of the entire session filtered in
% localdir/localmat
%The .mat file of the entire session has to be already been generated from the corresponding EDF file and
% band and notch filtered. The mat file has been already generated using the GUI of EEGLAB. 
%OUTPUT: it creates a mat file for the specified epoch times

%% Load the .mat file with the entire session. This section can be skipped if we load the file form the GUI
% localmat = 'dummy1chan-filternadnotch.mat'
% %patient , date and session must correspond with the patient of the entire sessiondefined in localmat
%patientlist = {'fo24' 'sb27' 'jc28'  'nk33' 'mj34' 'ms31' 'TWH32' 'TWH33' 'TWH34' 'TWH35'};
patientid= 'TWH024';
conditionlist = {'HYP' 'EC_POST' 'EO_POST' 'EC_PRE' 'EO_PRE' };
conditionindex = 1;
patientcond = conditionlist{conditionindex};

[mydir, myfile,patdate,patsession,patname] = getfullpathfrompatient(patientid,patientcond,'entire_session');
myfullname = fullfile(mydir, myfile);
disp('Loading mat file to get the EEG object...')
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
%% Specify the times to cut inthe session file initially load
% check excel file to see the initial and final time of the epoch
switch patientid
    case 'TWH024'
        if strcmp(patientcond,'HYP') == 1
            hmsdate_startsession = '11:55:29';
            hmsdate_initepoch = '12:45:57' ;
            hmsdate_endepoch = '12:48:40';
            se_patient = 'fo24';
        end
    case 'TWH025'
        if strcmp(patientcond,'HYP') == 1
        end
    case 'TWH037'
        if strcmp(patientcond,'HYP') == 1 
        else
        end
    case 'TWH038'
        if strcmp(patientcond,'HYP') == 1  
        else
        end
    otherwise
        warningMessage = sprintf('Error: patient%s does not exist:\n%s', patientid);
        uiwait(msgbox(warningMessage));
end

% get the seconds from the actualtime
% t0...t10.....t20, t21 is the duration in seconds 
%EEG = EEGentire;
disp(sprintf('From time to seconds, init session= %s , initial epoch= %s respectively', hmsdate_startsession, hmsdate_initepoch ))
[t10] = fromtime_to_seconds(hmsdate_startsession, hmsdate_initepoch);
disp(['seconds, t10 = ', num2str(t10)])
[t21] = fromtime_to_seconds(hmsdate_initepoch, hmsdate_endepoch);
disp(['seconds, t21 = ', num2str(t21)])
t20 = t21 + t10;
timelimits = [t10 t20];
disp(['seconds t20=', num2str(t20)])
%scale timelimits with the samploing rate 
datapointtimelimits = timelimits*EEG.srate;
disp(['Creating epoch between [t10, t20]s=' num2str(timelimits)])
% EEG = eeg_eegrej returbns the reamining set
rightend = [t20*EEG.srate EEG.pnts];
leftend = [0 t10*EEG.srate];
disp('getting the left side of the epoch')
EEGepoch0t20 =  eeg_eegrej(EEG, rightend);
disp('getting the final epoch')
EEGepocht10t20 =  eeg_eegrej(EEGepoch0t20, leftend);
%EEGepoch =  eeg_eegrej(EEG, datapointtimelimits);
disp(['EEGepoch created '])
[localdir, localfile, localextension] = fileparts(myfullname);
%cuttosave = 'EEG_cut_BL_HYP_fo24_s1.mat'
%cutfiletosave = fullfile(localdir, cuttosave);
%save cutfiletosave -struct EEGepoch;
% save .mat file
disp(['Saving EEGepoch as a mat file '])

savematfileforepoch(localdir, EEGepocht10t20, se_patient, patdate, patsession,patientcond)
disp(['Created mat file for :' se_patient, patdate, patsession])

end

function [] = savematfileforepoch(localdir, EEG, se_patient, se_date, nu_session,patientcond )
%% savematfile save mat file for the patient, condition etc specified in the function
% INPUT: directory where mat file is saved
% EEG object
% se_patient, se_date, nu_session to write the destination mat file name
%disp(['Patient, date and session are:' num2str(se_patient) , num2str(se_date), num2str(nu_session) ])
%EEG = EEGepoch
matfh = 'EEG_cut_';
TF = strncmpi(patientcond,'BL',2)
if TF == 1
    fprintf('savematfileforepoch for patient:%s and condition:%s', patientid,patientcond );
else
    patientcond = strcat('BL_', patientcond);
end
disp(['Condition is:' patientcond])
%matfh = sprintf('EEG_cut_%s_%s_%s_%s', conditionlist{conditionindex},se_patient{1}, se_date, nu_session)
matfh = sprintf('EEG_cut_%s_%s_%s_%s.mat',  patientcond, se_patient, se_date, nu_session)
disp(['mat file is:' matfh])
fullmatfile = fullfile(localdir, matfh);
oldFolder = cd(localdir)
disp(['Saving EEG epoch in file:', matfh])
save(matfh,'EEG') 
end

function [total_seconds_difference] = fromtime_to_seconds(hmsdate_first,hmsdate_last)
%FUNCTION_NAME - function that takes two strings with format hour:minute:second and returns the seconds, the first is the initial time and the second is the last time 
%and calculates the number of seconds that passed between the two instants 
%
% Inputs:
%    hmsdate_first - literal with the form 'h:m:s' for initial time
%    hmsdate_last -  literal with the form 'h:m:s' for final time
% Outputs:
%    total_seconds_difference  = sesoncdsinit- seconds final
% Example: 
%   [output1] = fromtime_to_seconds('12:34:59', '12:39:52')
%
[Y, M, D, H, MN, S] = datevec(hmsdate_first);
total_seconds_init = H*3600+MN*60+S;
[Y, M, D, H, MN, S] = datevec(hmsdate_last);
total_seconds_final = H*3600+MN*60+S;
total_seconds_difference = total_seconds_final - total_seconds_init;
end