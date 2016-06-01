function listofchannels = getlistofchannelsofinterest(brainregion)
% getlistofchannelsofinterest return the channel labels we are interested
% in
if strcmp(brainregion, 'temphd') == 1
    listofchannels = {'LHD1','LHD2','LHD3', 'LHD4', 'LAT1','LAT2'...
    ,'LAT3','LAT4', 'LMT1','LMT2', 'LMT3','LMT4','LPT1' ,'LPT2','LPT3','LPT4','LPT5','LPT6','RHD1'...
    ,'RHD2','RHD3','RHD4','RAT1','RAT2','RAT3','RAT4','RMT1','RMT2','RMT3', 'RMT4','RPT1','RPT2'...
    ,'RPT3','RPT4','RPT5','RPT6'};
elseif strcmp(brainregion, 'hd') == 1
        listofchannels = {'LHD1','LHD2','LHD3', 'LHD4','RHD1','RHD2','RHD3','RHD4'};
end
end