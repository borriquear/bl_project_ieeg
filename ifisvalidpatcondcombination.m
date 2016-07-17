function [existcond,hashd] =  ifisvalidpatcondcombination(eegpatient,eegcond)
%ifisvalidpatcondcombination returns 1 is the cobination patient and
%conditon is valid, tjat is, there is an epoch created
existcond = 1;
hashd = 1;
if (strcmp(eegpatient,'TWH030') && strcmp(eegcond, 'EO_POST')) || (strcmp(eegpatient,'TWH031') && strcmp(eegcond, 'EO_POST')) || (strcmp(eegpatient,'TWH034') && strcmp(eegcond, 'EO_PRE'))
    existcond = 0;
end
if (strcmp(eegpatient,'TWH037') + strcmp(eegpatient,'TWH043') > 0)
    hashd = 0;
end