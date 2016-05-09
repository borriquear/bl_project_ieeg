function [freq_band] = getgreeksymbolfreq(freq)
switch freq %2, 6 , 10, 23.5, 40};
    case 2
        freq_band = '\delta';
    case 6
        freq_band = '\theta';
    case 10
        freq_band = '\alpha';
    case 23.5
        freq_band = '\beta';
    case 40
        freq_band = '\gamma';
    otherwise
        fprintf('ERROR frequency band missing!!');
        return;
end
end