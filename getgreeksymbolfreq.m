function [freq_band] = getgreeksymbolfreq(freq)
switch freq %2, 6 , 10, 23.5, 40};
    case freq < 4 %2
        freq_band = '\delta';
    case freq > 4 && freq < 8
        freq_band = '\theta';
    case freq > 8 && freq < 12
        freq_band = '\alpha';
    case freq > 12 && freq < 30
        freq_band = '\beta';
    case freq > 30
        freq_band = '\gamma';
    otherwise
        fprintf('ERROR frequency band missing!!');
        return;
end
end