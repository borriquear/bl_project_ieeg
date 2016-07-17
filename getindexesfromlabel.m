function [channelclassindexes] = getindexesfromlabel(patientid, label)
% getindexesfromlabel returns the field indexes from the structure patientchannelClass
% for patientid and label
%example: getindexesfromlabel('TWH024', 'MT') ->[10 11 12 13 28 29 30 31]
%REQUIREMENT: mat file globalFsDir/channelclasses.mat has to exist
global globalFsDir;
[globalFsDir] = loadglobalFsDir();
filechannel = 'channelclasses.mat';
filechannel = fullfile(globalFsDir,filechannel);
fh = load(filechannel,'patientchannelClass');
if isequal(patientid, 'TWH024') == 1
    idpatnb = 1;
elseif isequal(patientid, 'TWH027') == 1
    idpatnb = 2;
elseif isequal(patientid, 'TWH028') == 1
    idpatnb = 3;
elseif isequal(patientid, 'TWH030') == 1
    idpatnb = 4;
elseif isequal(patientid, 'TWH031') == 1
    idpatnb = 5;
elseif isequal(patientid, 'TWH033') == 1
    idpatnb = 6;
elseif isequal(patientid, 'TWH034') == 1
    idpatnb = 7;
elseif isequal(patientid, 'TWH037') == 1
    idpatnb = 8;
elseif isequal(patientid, 'TWH038') == 1
    idpatnb = 9;
elseif isequal(patientid, 'TWH042') == 1
    idpatnb = 10;
elseif isequal(patientid, 'TWH043') == 1
    idpatnb = 11;
else
    error('ERROR: Patientif is incorrect, patient doesnt exist')
end
[truefalse, indexchannels] = ismember(label, fh.patientchannelClass(idpatnb).labels);
if truefalse <1
    channelclassindexes = [];
    fprintf('ERROR: Label doesnt exist\n')
else
    channelclassindexes = fh.patientchannelClass(idpatnb).indexes{indexchannels};
    %channelclassindexes = channelclassindexes-1;
end
end