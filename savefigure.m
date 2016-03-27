function [] = savefigure(myfullname, hfigure,figcounter,typeofchart)
%Save figure atfter showing it for 5 seconds
[eegpathname eegfilename eegextname]= fileparts(myfullname);
[eegcond, eegpatient,eegdate,eegsession ] = getsessionpatient(eegfilename);
olderfolder =  cd(eegpathname);
figurespath = fullfile(eegpathname,'figures');
if ~exist(figurespath, 'dir')
    mkdir('figures');
end
donotsave = 0; % 
olderfolder2 = cd('figures');
currentFolder = pwd;
%fprintf('Saving figure in %s \n',currentFolder);
closewindow = 1; %by default close window
if nargin > 3
    if strcmp(typeofchart,'barfreqallconditions')  == 1
        cd ..
        filefigname = sprintf('power_all_conditions_%s_%s.fig', eegdate,eegsession);
        fprintf('Saving %s ...\n',filefigname);
    elseif strcmp(typeofchart,'barfreqonecondition')  == 1
        filefigname = sprintf('power_per_condition_%s_%s_%s_%s.fig',eegcond, eegpatient,eegdate,eegsession);
        fprintf('Saving %s ...\n',filefigname);
    elseif strcmp(typeofchart,'powerperchannel')  == 1
        %figcounter assigned to chan number
        filefigname = sprintf('power_per_channel_%d_%s_%s_%s_%s.fig',figcounter, eegcond, eegpatient,eegdate,eegsession);
        fprintf('Saving %s ...',filefigname);     
    elseif strcmp(typeofchart,'barfreqonechannel') == 1
        filefigname = sprintf('freqbands_per_channel_%d_%s_%s_%s_%s.fig',figcounter, eegcond, eegpatient,eegdate,eegsession);
        fprintf('Saving %s ...\n',filefigname); 
    elseif strcmp(typeofchart,'powerrelativeperband') == 1
        filefigname = sprintf('power_relative_per_band_%s_%s_%s_%s.fig', eegcond, eegpatient,eegdate,eegsession);
        closewindow = 0;
    end
else
    filefigname = sprintf('fft_%d_%s_%s_%s_%s.fig',figcounter,eegcond, eegpatient,eegdate,eegsession);
    %fprintf('Saving %s ...',filefigname);
    donotsave =1; %do not save the 3 charts per channel
end
if donotsave == 0
    tic;
    savefig(hfigure,filefigname);
    toc;
    %imshow(hfigure);
    %pause(1.2);
end
if closewindow == 1
    close(hfigure);
end
end
