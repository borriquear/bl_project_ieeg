function EEGcorr = fixwronglabeling(patientid, EEGorig)
% change wrong labels in EEG 
EEGcorr = EEGorig;
switch patientid
    case 'TWH030'
            EEGcorr.chanlocs(end-1).labels = 'RPT5'
            EEGcorr.chanlocs(end).labels = 'RPT6'
    case 'TWH024'
            EEGcorr.chanlocs(2).labels = 'LHD1'; EEGcorr.chanlocs(3).labels = 'LHD2';EEGcorr.chanlocs(4).labels = 'LHD3';EEGcorr.chanlocs(5).labels = 'LHD4';
            EEGcorr.chanlocs(6).labels = 'LAT1'; EEGcorr.chanlocs(7).labels = 'LAT2';EEGcorr.chanlocs(8).labels = 'LAT3';EEGcorr.chanlocs(9).labels = 'LAT4';
            EEGcorr.chanlocs(10).labels = 'LMT1'; EEGcorr.chanlocs(11).labels = 'LMT2';EEGcorr.chanlocs(12).labels = 'LMT3';EEGcorr.chanlocs(13).labels = 'LMT4';
            EEGcorr.chanlocs(14).labels = 'LPT1'; EEGcorr.chanlocs(15).labels = 'LPT2';EEGcorr.chanlocs(16).labels = 'LPT3';EEGcorr.chanlocs(17).labels = 'LPT4';EEGcorr.chanlocs(18).labels = 'LPT5';EEGcorr.chanlocs(19).labels = 'LPT6';
            EEGcorr.chanlocs(20).labels = 'RHD1'; EEGcorr.chanlocs(21).labels = 'RHD2';EEGcorr.chanlocs(22).labels = 'RHD3';EEGcorr.chanlocs(23).labels = 'RHD4';
            EEGcorr.chanlocs(24).labels = 'RAT1'; EEGcorr.chanlocs(25).labels = 'RAT2';EEGcorr.chanlocs(26).labels = 'RAT3';EEGcorr.chanlocs(27).labels = 'RAT4';
            EEGcorr.chanlocs(28).labels = 'RMT1'; EEGcorr.chanlocs(29).labels = 'RMT2';EEGcorr.chanlocs(30).labels = 'RMT3';EEGcorr.chanlocs(31).labels = 'RMT4';
            EEGcorr.chanlocs(32).labels = 'RPT1'; EEGcorr.chanlocs(33).labels = 'RPT2';EEGcorr.chanlocs(34).labels = 'RPT3';EEGcorr.chanlocs(35).labels = 'RPT4';EEGcorr.chanlocs(36).labels = 'RPT5';EEGcorr.chanlocs(37).labels = 'RPT6';
end
end