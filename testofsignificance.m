%
function [tvalue,pvalue,t2value,Fvalue] = testofsignificance(HpatientsL, LpatientsL, conditionL, centerfrequencies)
% testofsignificance  given a list of patients group 1 , patients of group 2, electrodes and type(mono/multi) returns the
% t-test result to assess if we can separate the population into these two
% groups based on the mean of the electodesL
%   electrodesincommon = ttestofsignificance(H,L,ele,condL, centerfreq, t) e.g.
%   H={'TWH01',....,'TWH11' },L={'TWH21','TWH22' } conditionL = {'HYP'}, centerfrequencies = {2, 6 , 10, 23.5, 40};, type = 'mono/multi'
%   mono:eachvariable at a time, multi all the variables (B elements) at
%   the same time
%
%initialize output values
global globalFsDir;
disp('Loading globalFsDir...')
if ~exist('globalFsDir','var')
    fprintf('globalFsDir not found, loading it...')
    eval('global globalFsDir');
    myp = 'D:\BIAL PROJECT\patients\'
    eval(['globalFsDir=' 'myp']);
end
tvalue = 0;
pvalue = 0;
t2value = 0;
Fvalue = t2value;
% if strcmp(type, 'mono') == 1
% end
Allpatients = [HpatientsL,LpatientsL];
% if ~exist('electrodesL','var')
%   %electrodesL not invoked calculate the intersection of the elctrodes for
%   %all the patients
fprintf('Calculating the list electrode in common for all the patients...');
[electrodesL] = calculateintersectionchannels(Allpatients);
if isempty(electrodesL)
    msgerror = 'The patients have not any electrodes in common!!!END ';
    error(msgerror);
    %return
else
    fprintf('The list of electrodes in common are=');
    disp(electrodesL)
end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HpatientsL = {'TWH030'};
% LpatientsL = {'TWH031'};
% conditionL = {'HYP'};
% centerfrequencies = {2, 6 , 10, 23.5, 40};
% electrodesL = {'LHD1','LHD2','LHD3','LHD4'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


matfileid = 'powerconnectivity_freq_'; %header of matfile that contains the corr matrix
%Open the mat file that contains the variable for each condition
disp('Loading variables for condition H....')

% populationsample = {'pat1','pat2'; [1xnumel double],[1xnumel double]}
populationsample = containers.Map;
% populationsampleperfreq = containers.Map;
% populationsampleperfreqandcond = containers.Map;
%anidated data structures, populationsampleperfreqandcond('condition') =
%'HYP, populationsampleperfreqandcond('data') = populationsampleperfreq
for condind = 1:length(conditionL)
    currcond = conditionL{condind};
    %populationsampleperfreqandcond('condition') = currcond;
    hlist = {};
    plist = {};
    figh = figure;
    for freqind= 1:length(centerfrequencies)
        currfreq = centerfrequencies{freqind};
        %populationsampleperfreq('frequency')= currfreq;
        for indpat=1:length(Allpatients)
            currentpat = Allpatients{indpat};
            %fprintf('The current patient is:%s\n',currentpat);
            patpath = strcat(globalFsDir,currentpat);
            fprintf('The current path is %s\n', patpath);
            fprintf('Loading the EEG object for patient  %s\n',currentpat);
            [myfullname, EEG, channel_labels, patdate, patsession] = initialize_EEG_variables(currentpat,currcond);
            %open mat file and select the electrodesL in that mat file
            %powerconnectivity_freq_40_HYP_TWH030_11172015_s1
            mattoload = strcat(matfileid,num2str(currfreq),'_',currcond,'_', currentpat,'_',patdate,'_',patsession,'.mat');
            fftfile = fullfile(patpath,'data','figures', mattoload);
            fprintf('Opening correlation matrix....%s\n',fftfile);
            matf= matfile(fftfile);
            corrmatpersubband = matf.corr_matrix;
            chan_labels = matf.channel_labels;
            %get the row numbersof the electrodes i am interested in
            counterw = 1;
            fprintf('Finding the subset of the corr matrix for the electrodes list...:\n');
            for indchan =1:length(electrodesL)
                [elfound indelf] = find(ismember(chan_labels, electrodesL{indchan}));
                if elfound == 1
                    % vector of indices of electrodes in electrodesL
                    %e.g. is electrodesL = LHD1 and this is corresponds to
                    %the 36 electrode in the corr matrix,
                    %indexofcommonrows(1) =36
                    indexofcommonrows(counterw) = indelf-1;
                    counterw = counterw +1;
                end
            end
            disp(indexofcommonrows);
            %colapse the corr matrix in a vector, vecofcorr(i) (1xEEG.nbchan)= sum in
            %absvalue of the corr of electrode i with all the other
            %electrodes
            vecofcorr = sum(abs(corrmatpersubband));
            %vecofcorr = sum(corrmatpersubband);
            subsetofcorrmat = vecofcorr(indexofcommonrows);
            populationsample(currentpat) = subsetofcorrmat;
        end
        %populationsampleperfreq('data')= populationsample;
        %save object populationsample a matfile, for that freq and
        %condition
        newmatfileid = 'subsetofcorr_specchan_frq_';
        newmattoload = strcat(newmatfileid,num2str(currfreq),'_',currcond,'.mat');
        dirsaveallpatdata = 'figures_all_patients';
        destdir = fullfile(globalFsDir,dirsaveallpatdata);
        if exist(destdir, 'dir')
            fprintf('Saving mat fil at %s\n', destdir);
        else
            fprintf('mkdir ...%s\n', destdir);
            mkdir(destdir);
        end
        newfftfile = fullfile(destdir, newmattoload);
        fprintf('Saving subset of correlation matrix in %s\n',newfftfile);
        save(newfftfile,'populationsample', 'electrodesL')
        [h,p] = ttest2HvsL(newfftfile, HpatientsL, LpatientsL);
        hlist{freqind} = h;
        plist{freqind} = p;
    end
    %populationsampleperfreqandcond('data')=populationsampleperfreq;
    %print chart with results for all frequency bands
    drawttestresults(figh, HpatientsL, LpatientsL,electrodesL,hlist,plist,centerfrequencies);
end
fprintf('The final variables for the population and electrodes list is..:\n');
disp(keys(populationsample));
disp(values(populationsample));
[selectedpatients] = keys(populationsample);
[selectedcorrvecs] = values(populationsample);
fprintf('List of electrodes in common=');disp(electrodesL);
for i=1:length(selectedpatients)
    currpat = selectedpatients{i};
    currvalues = selectedcorrvecs{i};
    fprintf('For patient %s abs corr values =',currpat);
    disp(currvalues); fprintf('\n');
end
end

function [h,p] = ttest2HvsL(newfftfile, HpatientsL, LpatientsL)
%load mat file that contain the population variables
load(newfftfile,'-mat','populationsample', 'electrodesL')
%separate sample1 from sample 2
[sampleH, sampleL] = getpopulationsample(populationsample,HpatientsL, LpatientsL)
%do the t-test
fprintf('Calculating the t test...\n');
%remove empty cells
emptyCells = cellfun(@isempty,sampleL);
sampleL(emptyCells) = [];
emptyCells = cellfun(@isempty,sampleH);
sampleH(emptyCells) = [];
[h,p,ci,stats] = ttest2(horzcat(sampleH{1:end}),horzcat(sampleL{1:end}));
fprintf('h = %d , p= %.2f', h, p);
fprintf('ci and stats :');
disp(ci); disp(stats);
if h == 1
    fprintf('T-test2 does reject the null hypothesis at the default 5\% significance level.\n');
    fprintf('The variables are significantly different in High vs Lows\n');
elseif h == 0
    fprintf('T-test2 does not reject the null hypothesis at the default 5\% significance level.\n');
    fprintf('There is no statiscal justifcation to separate H vs L\n');
end

end


function [sampleH, sampleL] = getpopulationsample(populationsample, HpatientsL, LpatientsL)
sampleH = {};
sampleL = {};
[selectedpatients] = keys(populationsample);
[selectedcorrvecs] = values(populationsample);
for i=1:length(selectedpatients)
    currpat = selectedpatients{i};
    currvalues = selectedcorrvecs{i};
    [elfoundH indelfH] = find(ismember(HpatientsL, currpat));
    [elfoundL indelfL] = find(ismember(LpatientsL, currpat));
    if  elfoundH == 1
        sampleH{i} = currvalues;
    end
    if elfoundL == 1 %this allows having patients H and L , for exclusive H or L elseif
        sampleL{i} = currvalues;
    end
    if elfoundH + elfoundL == 0
        error('Patient do not found!!!');
    end
end
end

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

function [] = drawttestresults(figh, HpatientsL, LpatientsL,electrodesL,hlist,plist, centerfrequencies)

%freq name from currfreq
global globalFsDir;
figure(figh);
currcond = 'HYP';
xchannel_limit = length(electrodesL);
xlabeljump = 1;
paintedlegalready = 0;
for indexchart=1:length(centerfrequencies)
    %openfile
    currfreq = centerfrequencies{indexchart};
    newmatfileid = 'subsetofcorr_specchan_frq_';
    newmattoload = strcat(newmatfileid,num2str(currfreq),'_',currcond,'.mat');
    dirsaveallpatdata = 'figures_all_patients';
    newmattoload = fullfile(globalFsDir,dirsaveallpatdata,newmattoload);
    load(newmattoload,'-mat','populationsample', 'electrodesL');
    [sampleH, sampleL] = getpopulationsample(populationsample,HpatientsL, LpatientsL)
    subplot(2,3,indexchart);
    %plot
    hold on;
    %set(gca,'LineWidth',2, 'Color', 'r')
    cellfun(@plot,sampleH);
    
    %set(p,'Color','red')
    cellfun(@plot,sampleL);
    [freq_band] = getgreeksymbolfreq(currfreq);
    title(freq_band); %max(max(powerperelec))+0.4
    set(gca, 'XLim', [1 xchannel_limit],'XTick',[1:xlabeljump:xchannel_limit],'XTickLabel',[electrodesL]);
    xticklabel_rotate([],45,[],'Fontsize',6);
    %set(gca, 'XLim', [1 xchannel_limit],'ylim', [0 1],'XTick',[1:xlabeljump:length(channels_pat)-1]);
    xlabel('Channels'), ylabel('Power-based Spearman corr.');
    grid on ;
    if paintedlegalready == 0
        legend([HpatientsL,LpatientsL]);
        paintedlegalready = 1;
    end
    pandhvals=  sprintf('h=%s, p=%s',num2str(hlist{indexchart}),num2str(hlist{indexchart}) );
    text(indexchart*2, indexchart*2,pandhvals);
    
    
end
if indexchart == length(centerfrequencies)
    highspat=  sprintf(' %s ', HpatientsL{1:end});
    highspat= strcat('Highs:',highspat);
    highspat= strcat(highspat, '  ');
    lowspat=  sprintf(' %s ', LpatientsL{1:end});
    lowspat= strcat('Lows:',lowspat);
    allpatleg = strcat(highspat,lowspat);
    xlabel(allpatleg)
end
% MyBox = uicontrol('style','text')
% set(MyBox,'Hight Patients: Low Patients', HpatientsL, LpatientsL)
% xpos= 400; ypos=50; xsize= 100; ysize = 100;
% set(MyBox,'Position',[xpos,ypos,xsize,ysize])
end