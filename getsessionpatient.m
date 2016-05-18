function [eegcond, eegpatient,eegdate,eegsession ] = getsessionpatient(eegfilename)
%% [eegcond, eegpatient,eegdate,eegsession ] = getsessionpatient(eegfilename)  get patient name etc from file name
%% returns the condition pateint name date and session given the complete path of the signal(mat file)
[eegbl remain] = strtok(eegfilename, 'EEG_cut_');
[eegblcond remain]= strtok(remain, '_');
[eegblcond2 remain]= strtok(remain, '_');
if (strcmpi(eegblcond2,'PRE') == 1) || (strcmpi(eegblcond2,'POST') == 1)
    eegcond = strcat(eegbl, eegblcond);
    eegcond = strcat(eegcond, eegblcond2);
    [eegpatient remain]= strtok(remain, '_');
    [eegdate remain]= strtok(remain, '_');
    eegsession = strtok(remain, '_');
else
    eegcond = strcat(eegbl, eegblcond);
    eegpatient = eegblcond2;
    [eegdate eegsession]= strtok(remain, '_');
    eegsession = strtok(eegsession, '_');
end
end