function [dataobject_cast] = castimplantwithlabel(dataobject, typeobject, label)
if strcmp(typeobject, 'phaseconn_matrices') == 1
    dataobject_cast = castimplantwithlabel_phaseconn(dataobject, label)
elseif strcmp(typeobject, 'powerconn_matrices') == 1
elseif strcmp(typeobject, 'wiring_matrices') == 1
    dataobject_cast= castimplantwithlabel_wiring(dataobject, label)
end
end
function [dataobject_cast] = castimplantwithlabel_wiring(dataobject, label)
% castimplantwithlabel_wiring returns the data structure casted according to the
% label which refer to elctrode location for wiring_matrcies
%example: castimplantwithlabal(dataobject, 'T') -> dataobject_casted
%REQUIREMENT: mat file globalFsDir/channelclasses.mat has to exist
%IN: dataobject wiring_matrices, label 'H', 'T','F','IH','FP', 'Grid'
global globalFsDir;
[globalFsDir] = loadglobalFsDir();
% filechannel = 'channelclasses.mat';
% filechannel = fullfile(globalFsDir,filechannel);
% fh = load(filechannel,'patientchannelClass');
dataobject_cast = dataobject;
nbpats = size(dataobject.patientslist,2);
freqlist = dataobject.frequencylist;
nbfreqs = length(freqlist);
nbconds = size(dataobject.conditionslist, 2);
listtodelete = [];
for i=1:nbpats
    patientid = dataobject.patientslist{i};
    channelclassindexes = getindexesfromlabel(patientid, label);
    channelclassindexes = channelclassindexes -1;
    eucldistmatrix = dataobject.distMatrixcell{i};
    eucldistmatrix = eucldistmatrix([channelclassindexes],[channelclassindexes]);
    dataobject_cast.distMatrixcell{i}= eucldistmatrix;
    if length(channelclassindexes) > 0
        %casting patient
        fprintf('Casting patient %s for Label %s \n',patientid, label{1})
        for j=1:nbconds
%             matricesfqs_ispc = dataobject.wiring_ispc{i,j};
%             matricesfqs_pli = dataobject.wiring_pli{i,j};
% %             matricesfqs_icoh = dataobject.pli_matrix{i,j};
%             matricesfqs_power = dataobject.wiring_power{i,j};
            for k=1:nbfreqs
                matrix_ispc = dataobject.wiring_ispc{i,j,k};%matricesfqs_ispc(:,:,k);
                matrix_ispc_cast = matrix_ispc([channelclassindexes],[channelclassindexes]);
                dataobject_cast.wiring_ispc{i,j,k}=matrix_ispc_cast;
                matrix_pli = dataobject.wiring_pli{i,j,k};
                matrix_pli_cast = matrix_pli([channelclassindexes],[channelclassindexes]);
                dataobject_cast.wiring_pli{i,j,k}=matrix_pli_cast;
                matrix_power =dataobject.wiring_power{i,j,k};
                matrix_power_cast = matrix_power([channelclassindexes],[channelclassindexes]);
                dataobject_cast.wiring_power{i,j,k}=matrix_power_cast;
%                 matrix_icoh = matricesfqs_ispc(:,:,k);
%                 matrix_icoh_cast = matrix_icoh([channelclassindexes],[channelclassindexes]);
%                 dataobject_cast.icoh_matrix{i,j,k}=matrix_icoh_cast;
            end
        end
    else
        %delete patient
        listtodelete= [listtodelete i];
        fprintf('Patient %s has not label %s, deleting from cast object\n', patientid, label{1});
    end
end
%delete patients that dont have label if any
if length(listtodelete)> 0
    fprintf('Deleting patients missing the label ... \n');
    dataobject_cast.wiring_ispc([listtodelete],:,:) = [];
    dataobject_cast.wiring_pli([listtodelete],:,:) = [];
    dataobject_cast.wiring_power([listtodelete],:,:) = [];
%     dataobject_cast.icoh_matrix([listtodelete],:) = [];
    dataobject_cast.patientslist([listtodelete]) = [];
    dataobject_cast.distMatrixcell([listtodelete]) = [];
end
end


function [dataobject_cast] = castimplantwithlabel_phaseconn(dataobject, label)
% castimplantwithlabel_phaseconn returns the data structure casted according to the
% label which refer to elctrode location for phaseconn_matrcies
%example: castimplantwithlabal(dataobject, 'T') -> dataobject_casted
%REQUIREMENT: mat file globalFsDir/channelclasses.mat has to exist
%IN: dataobject wiring_matrices, label 'H', 'T','F','IH','FP', 'Grid'
global globalFsDir;
[globalFsDir] = loadglobalFsDir();
% filechannel = 'channelclasses.mat';
% filechannel = fullfile(globalFsDir,filechannel);
% fh = load(filechannel,'patientchannelClass');
dataobject_cast = dataobject;
nbpats = size(dataobject.patientsl,2);
freqlist = dataobject.freqsl;
nbfreqs = length(freqlist);
nbconds = size(dataobject.conditionsl, 2);
listtodelete = [];
for i=1:nbpats
    patientid = dataobject.patientsl{i};
    channelclassindexes = getindexesfromlabel(patientid, label);
    channelclassindexes = channelclassindexes -1;
    if length(channelclassindexes) > 0
        %casting patient
        for j=1:nbconds
            matricesfqs_ispc = dataobject.ispc_matrix{i,j};
            matricesfqs_pli = dataobject.pli_matrix{i,j};
            matricesfqs_icoh = dataobject.pli_matrix{i,j};
            for k=1:nbfreqs
                matrix_ispc = matricesfqs_ispc(:,:,k);
                matrix_ispc_cast = matrix_ispc([channelclassindexes],[channelclassindexes]);
                dataobject_cast.ispc_matrix{i,j,k}=matrix_ispc_cast;
                matrix_pli = matricesfqs_pli(:,:,k);
                matrix_pli_cast = matrix_pli([channelclassindexes],[channelclassindexes]);
                dataobject_cast.pli_matrix{i,j,k}=matrix_pli_cast;
                matrix_icoh = matricesfqs_ispc(:,:,k);
                matrix_icoh_cast = matrix_icoh([channelclassindexes],[channelclassindexes]);
                dataobject_cast.icoh_matrix{i,j,k}=matrix_icoh_cast;
            end
        end
    else
        %delete patient
        listtodelete= [listtodelete i];
        fprintf('Patient %s has not label %s, deleting from cast object\n', patientid, label);
    end
end
%delete patients that dont have label if any
if length(listtodelete)> 0
    fprintf('Deleting patients missing the label ... \n')
    dataobject_cast.ispc_matrix([listtodelete],:) = [];
    dataobject_cast.pli_matrix([listtodelete],:) = [];
    dataobject_cast.icoh_matrix([listtodelete],:) = [];
    dataobject_cast.patientsl([listtodelete]) = [];
end
end