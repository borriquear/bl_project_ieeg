function [freq_band] = getgreeksymbolfreq(freq)
freq = double(freq);
if  freq < 4.0 %2
        freq_band = '\delta';
elseif freq > 4 && freq < 8
        freq_band = '\theta';
elseif freq > 8 && freq < 12
        freq_band = '\alpha';
elseif freq > 12 && freq < 30
        freq_band = '\beta';
elseif freq > 30
        freq_band = '\gamma';
else
        fprintf('ERROR frequency band missing!!\n');
        return;
end
end