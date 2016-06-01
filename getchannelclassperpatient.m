function [channelclass] = getchannelclassperpatient(patientid)
% getchannelclassperpatient returns object channelclass with the classses
% of channel per patientid, mat file channelclasses_patient.mat needs to
% exist, is created manually
%IN: patientid eg 'TWH037'
%OUT channelclass = {'frontal', [inex of frontal channels], [labels of frontal channels]; 'lateral', [],[] ;}
[myfullname, EEG, channel_labels,eegdate,eegsession] = initialize_EEG_variables(eegpatientl{indpat},eegcondition);
%open myfullname/channelclasses_patient.mat

end