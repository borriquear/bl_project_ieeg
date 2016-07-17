function [isinterictal] = interictalspikeschannels(patient, irow)
%% return 1 is the channel irow of patient is a channel with interictal spikes
isinterictal = 0;
if (strcmp(patient, 'TWH034') ==1 && (irow == 68)) ||  (strcmp(patient, 'TWH030') == 1 && (irow == 8)) || (strcmp(patient, 'TWH031') ==1 && (irow == 17) )
    %LHD4$ , irow = 68 (twh34)
    %LAT4, irow = 8  (twh30)
    %LMF5, irow = 17 (twh031)
    fprintf('Channel number %d, patient: %s has INTERICTAL spikes\n', irow, patient);
    isinterictal = 1;   
end
end