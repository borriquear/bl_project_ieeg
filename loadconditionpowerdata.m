function [quantitytomeasure,frqperband ] = loadconditionpowerdata(myfullname,patientid,condition)
%% load power data , mat file, containing the variable for a specific condition
% a patient
fprintf('Calling loadconditionpowerdata function, patient%s, %condition..\n',patientid,condition)
fprintf('Cleaning workspace..\n')
% make sure that globalFsDir is assigned
[globalFsDir] = loadglobalFsDir();
patpath = strcat(globalFsDir,patientid);
patpath = fullfile(patpath, 'data','figures');
[mydir, myfile,patdate,patsession] = getfullpathfrompatient(patientid,condition);
matfileh = strcat('fft_',condition,'_', patientid,'_',patdate,'_',patsession, '.mat' );
fileh = fullfile(patpath,matfileh);

fprintf('Loading in memory file %s:\n',fileh);
load(fileh,'quantitytomeasure','frqperband');
fprintf('Done');

end