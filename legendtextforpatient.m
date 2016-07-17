function [legend_text, veccondinexes] = legendtextforpatient(eegpatient)
legend_text{1}= 'EC\_PRE';legend_text{2}='EO\_PRE' ;legend_text{3}= 'HYP' ;legend_text{4}= 'EC\_POST';legend_text{5}= 'EO\_POST';
veccondinexes = [1 2 3 4 5];
if strcmp(eegpatient, 'TWH030') || strcmp(eegpatient, 'TWH031')
    legend_text = [];legend_text{1}= 'EC\_PRE';legend_text{2}='EO\_PRE' ;legend_text{3}= 'HYP' ;legend_text{4}= 'EC\_POST';
    veccondinexes = []; veccondinexes = [1 2 3 4];
elseif strcmp(eegpatient,'TWH034')
    legend_text = [];legend_text{1}= 'EC\_PRE' ;legend_text{2}= 'HYP' ;legend_text{3}= 'EC\_POST';legend_text{4}= 'EO\_POST';
    veccondinexes = [];veccondinexes = [1 3 4 5];
end
end