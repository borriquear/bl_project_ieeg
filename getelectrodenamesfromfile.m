function celloflabelswithq = getelectrodenamesfromfile(patientid, patientcond)
%%  getelectrodenamesfromfile(patientid, patientcondition) DEPRECATED FUNCTION DELETE
% IN:patientid, patientcondition
% OUT: celloflabelswithq : a cell array containing: electrodevector(row) = 'label' from elec_recon\electrodeNames
%     listtoplot  listtoplot(5) = quantitytomeasure(5), for the 5th element
%     in the elctrodeNodes file
[mydir, myfile,patdate,patsession] = getfullpathfrompatient(patientid,patientcond);

% make sure that globalFsDir is assigned
if ~exist('globalFsDir','var') 
   fprintf('globalFsDir not found, loading it...')
   eval('global globalFsDir');
   myp = 'D:\BIAL PROJECT\patients\';
   eval(['globalFsDir=' 'myp']); 
end
patpath = strcat(globalFsDir,patientid);
patpath = fullfile(patpath, 'elec_recon');
elecnodesname = strcat(patientid, '.electrodeNames');
elecnodesname = fullfile(patpath, elecnodesname);
fprintf('Loading  %s.electrodeNames file\n', patientid);
%load(elecnodesname);
%fprintf('Loading  %s.electrodeNames file\n', patientid)
fid = fopen(elecnodesname,'rt');
if fid == -1
    error('File :%s could not been found',elecnodesname);
end
tmp = textscan(fid, '%s %s %s', 'Headerlines', 2);
electrodevector = tmp{1};
nbchannels = size(tmp{1},1);

fprintf('Array for %s.electrodeNames file has been created\n', patientid);

mattoload = strcat('fft_',patientcond,'_', patientid,'_',patdate,'_',patsession,'.mat');
fprintf('Loading variables quantitytomeasure and channel_labels from %s file %s.electrodeNames file has been created\n', mattoload);
fftfile = fullfile(mydir,'figures', mattoload);
%load(fftfile);
load(fftfile,'quantitytomeasure','channel_labels','percentlistoffrqperband');
%build hash vector containing the quantitytomeasure to plot the
%electrodecolors
for i=1:size(quantitytomeasure,1)
    %index of i_th element in the electrodeNode file
    labelelectrode = electrodevector(i);
    found = sum(strcmpi({channel_labels{2:end}},labelelectrode));
    if found == 1
        indexh = find(strcmpi({channel_labels{2:end}},labelelectrode),1);
        %fprintf('Found label:%s in index:%d\n',labelelectrode{1},indexh);
        listtoplot(i) = quantitytomeasure(indexh);
    else
        fprintf('Label %s NOT FOND!!!\n',labelelectrode{1});
        listtoplot(i) = fixwronglabels(labelelectrode, patientid,quantitytomeasure);
    end
end
celloflabelswithq = {electrodevector,listtoplot',percentlistoffrqperband};
lablist = celloflabelswithq{1};
vallab = celloflabelswithq{2};
for i=1:size(celloflabelswithq{1},1)
    fprintf('Label:%s, Value:%.4f\n', lablist{i},vallab(i));
end
fclose(fid);
end
function [quvalueforlabel] = fixwronglabels(labelelectrode, patientid,quantitytomeasure)
if strcmpi(patientid, 'TWH030')== 1
    if strcmpi(labelelectrode, 'RPT6') == 1
        quvalueforlabel = quantitytomeasure(end);
    elseif strcmpi(labelelectrode, 'RPT5') == 1
        quvalueforlabel = quantitytomeasure(end-1);
    end
else
    %put something if i do not find the electrode name in the elctrodename
    %file
    quvalueforlabel= quantitytomeasure(1);
end    
end