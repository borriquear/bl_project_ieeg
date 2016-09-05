
matf = 'EEG_cut_EC_POST_TWH048_06292016_s1.mat';
fh = fullfile('D:\BIAL PROJECT\patients\TWH048\data\', matf);
S = load(fh);
count = 0;
EEGepocht10t20 = S.EEGepocht10t20;
for i=1:EEGepocht10t20.nbchan
    ilabel = EEGepocht10t20.chanlocs(i).labels;
    if  strcmp(ilabel,'VentD1') == 1 EEGepocht10t20.chanlocs(i).labels = 'LVD1'
    elseif strcmp(ilabel,'VentD2') == 1 EEGepocht10t20.chanlocs(i).labels = 'LVD2'
    elseif strcmp(ilabel,'VentD3') == 1 EEGepocht10t20.chanlocs(i).labels = 'LVD3'
    elseif strcmp(ilabel,'VentD4') == 1 EEGepocht10t20.chanlocs(i).labels = 'LVD4'
    elseif strcmp(ilabel,'MidD1') == 1 EEGepocht10t20.chanlocs(i).labels = 'LMD1'
    elseif strcmp(ilabel,'MidD2') == 1 EEGepocht10t20.chanlocs(i).labels = 'LMD2'
    elseif strcmp(ilabel,'MidD3') == 1 EEGepocht10t20.chanlocs(i).labels = 'LMD3'
    elseif strcmp(ilabel,'MidD4') == 1 EEGepocht10t20.chanlocs(i).labels = 'LMD4'
    elseif strcmp(ilabel,'DorsD1') == 1 EEGepocht10t20.chanlocs(i).labels = 'LDD1'
    elseif strcmp(ilabel,'DorsD2') == 1 EEGepocht10t20.chanlocs(i).labels = 'LDD2'
    elseif strcmp(ilabel,'DorsD3') == 1 EEGepocht10t20.chanlocs(i).labels = 'LDD3'
    elseif strcmp(ilabel,'DorsD4') == 1 EEGepocht10t20.chanlocs(i).labels = 'LDD4'
        
    else
        count = count + 1;
    end

end
    disp([count, EEGepocht10t20.nbchan] )
    save(fh, 'EEGepocht10t20');
%%
matf = 'EEG_cut_EC_POST_TWH037_03142016_s1.mat';
fh = fullfile('D:\BIAL PROJECT\patients\TWH037\data\', matf);
S = load(fh);
count = 0;
EEGepocht10t20 = S.EEGepocht10t20;

for i=1:EEGepocht10t20.nbchan
    ilabel = EEGepocht10t20.chanlocs(i).labels;
    if  strcmp(ilabel,'LAT1') == 1 EEGepocht10t20.chanlocs(i).labels = 'LAF1'
    elseif strcmp(ilabel,'LAT2') == 1  S.EEGepocht10t20.chanlocs(i).labels = 'LAF2'
    elseif strcmp(ilabel,'LAT3') == 1  EEGepocht10t20.chanlocs(i).labels = 'LAF3'
    elseif strcmp(ilabel,'LAT4') == 1  EEGepocht10t20.chanlocs(i).labels = 'LAF4'
    elseif strcmp(ilabel,'LAT5') == 1  EEGepocht10t20.chanlocs(i).labels = 'LAF5'
    elseif strcmp(ilabel,'LAT6') == 1  EEGepocht10t20.chanlocs(i).labels = 'LAF6'
    elseif strcmp(ilabel,'LAT7') == 1  EEGepocht10t20.chanlocs(i).labels = 'LAF7'
    elseif strcmp(ilabel,'LAT8') == 1  EEGepocht10t20.chanlocs(i).labels = 'LAF8'
    elseif strcmp(ilabel,'LMT1') == 1  EEGepocht10t20.chanlocs(i).labels = 'LMF1'
    elseif strcmp(ilabel,'LMT2') == 1  EEGepocht10t20.chanlocs(i).labels = 'LMF2'
    elseif strcmp(ilabel,'LMT3') == 1  EEGepocht10t20.chanlocs(i).labels = 'LMF3'
    elseif strcmp(ilabel,'LMT4') == 1  EEGepocht10t20.chanlocs(i).labels = 'LMF4'
    elseif strcmp(ilabel,'LMT5') == 1  EEGepocht10t20.chanlocs(i).labels = 'LMF5'
    elseif strcmp(ilabel,'LMT6') == 1  EEGepocht10t20.chanlocs(i).labels = 'LMF6'
    elseif strcmp(ilabel,'LMT7') == 1  EEGepocht10t20.chanlocs(i).labels = 'LMF7'
    elseif strcmp(ilabel,'LMT8') == 1  EEGepocht10t20.chanlocs(i).labels = 'LMF8'
    elseif strcmp(ilabel,'LPT1') == 1  EEGepocht10t20.chanlocs(i).labels = 'LPF1'
    elseif strcmp(ilabel,'LPT2') == 1  EEGepocht10t20.chanlocs(i).labels = 'LPF2'
    elseif strcmp(ilabel,'LPT3') == 1  EEGepocht10t20.chanlocs(i).labels = 'LPF3'
    elseif strcmp(ilabel,'LPT4') == 1  EEGepocht10t20.chanlocs(i).labels = 'LPF4'
    elseif strcmp(ilabel,'LPT5') == 1  EEGepocht10t20.chanlocs(i).labels = 'LPF5'
    elseif strcmp(ilabel,'LPT6') == 1  EEGepocht10t20.chanlocs(i).labels = 'LPF6'
    elseif strcmp(ilabel,'LPT7') == 1  EEGepocht10t20.chanlocs(i).labels = 'LPF7'
    elseif strcmp(ilabel,'LPT8') == 1  EEGepocht10t20.chanlocs(i).labels = 'LPF8'
        
        
    elseif  strcmp(ilabel,'RAT1') == 1 EEGepocht10t20.chanlocs(i).labels = 'RAF1'
    elseif strcmp(ilabel,'RAT2') == 1  EEGepocht10t20.chanlocs(i).labels = 'RAF2'
    elseif strcmp(ilabel,'RAT3') == 1  EEGepocht10t20.chanlocs(i).labels = 'RAF3'
    elseif strcmp(ilabel,'RAT4') == 1  EEGepocht10t20.chanlocs(i).labels = 'RAF4'
    elseif strcmp(ilabel,'RAT5') == 1  EEGepocht10t20.chanlocs(i).labels = 'RAF5'
    elseif strcmp(ilabel,'RAT6') == 1  EEGepocht10t20.chanlocs(i).labels = 'RAF6'
    elseif strcmp(ilabel,'RAT7') == 1  EEGepocht10t20.chanlocs(i).labels = 'RAF7'
    elseif strcmp(ilabel,'RAT8') == 1  EEGepocht10t20.chanlocs(i).labels = 'RAF8'
    elseif strcmp(ilabel,'RMT1') == 1  EEGepocht10t20.chanlocs(i).labels = 'RPF1'
    elseif strcmp(ilabel,'RMT2') == 1  EEGepocht10t20.chanlocs(i).labels = 'RPF2'
    elseif strcmp(ilabel,'RMT3') == 1  EEGepocht10t20.chanlocs(i).labels = 'RPF3'
    elseif strcmp(ilabel,'RMT4') == 1  EEGepocht10t20.chanlocs(i).labels = 'RPF4'
    elseif strcmp(ilabel,'RMT5') == 1  EEGepocht10t20.chanlocs(i).labels = 'RPF5'
    elseif strcmp(ilabel,'RMT6') == 1  EEGepocht10t20.chanlocs(i).labels = 'RPF6'
    elseif strcmp(ilabel,'RMT7') == 1  EEGepocht10t20.chanlocs(i).labels = 'RPF7'
    elseif strcmp(ilabel,'RMT8') == 1  EEGepocht10t20.chanlocs(i).labels = 'RPF8'
    elseif strcmp(ilabel,'RPT1') == 1  EEGepocht10t20.chanlocs(i).labels = 'RPF1'
    elseif strcmp(ilabel,'RPT2') == 1  EEGepocht10t20.chanlocs(i).labels = 'RPF2'
    elseif strcmp(ilabel,'RPT3') == 1  EEGepocht10t20.chanlocs(i).labels = 'RPF3'
    elseif strcmp(ilabel,'RPT4') == 1  EEGepocht10t20.chanlocs(i).labels = 'RPF4'
    elseif strcmp(ilabel,'RPT5') == 1  EEGepocht10t20.chanlocs(i).labels = 'RPF5'
    elseif strcmp(ilabel,'RPT6') == 1  EEGepocht10t20.chanlocs(i).labels = 'RPF6'
    elseif strcmp(ilabel,'RPT7') == 1  EEGepocht10t20.chanlocs(i).labels = 'RPF7'
    elseif strcmp(ilabel,'RPT8') == 1  EEGepocht10t20.chanlocs(i).labels = 'RPF8'
    else
        count = count + 1;
        
        %printf('not temp = cunt %d / %d \n', count, EEGepocht10t20.nbchan)
    end
end
disp([count, EEGepocht10t20.nbchan] )
save(fh, 'EEGepocht10t20');
