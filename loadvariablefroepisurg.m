function [quantitytomeasure, channel_labels ] = loadvariablefroepisurg(myfullname)
%% loadvariablefroepisurg load mat file (fft_condition_patient) containing the variable for a specific condition
% a patient and returns the variables we are interested in
%[globalFsDir] = loadglobalFsDir();
fileh = myfullname;

fprintf('Loading in memory file %s:\n',fileh);
load(fileh,'percentlistoffrqperband','channel_labels');
fprintf('percentlistoffrqperband is \n');
quantitytomeasure = percentlistoffrqperband;
disp(quantitytomeasure);

end