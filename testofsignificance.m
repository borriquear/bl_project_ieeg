%
function [pcaresults] = testofsignificance(HpatientsL, LpatientsL, conditionL, centerfrequencies)
% testofsignificance  given a list of patients group 1, patients of group 2, 
%electrodes and type(mono/multi) returns the
% principal components and the t-hotelling test per condition and frequency
% band for all patients (highs and lows)
% IN: HpatientsL, LpatientsL, conditionL, centerfrequencies
% OUT: pcaresults ia a (cond,freq) cell .
% pcaresults{cond,freq}= {1x4 cell} coeff,score,latent,tsquare. coeff is the eigenvector matrix Z = Coeff X. The scores are 
% the data formed by transforming the original data into the 
% space of the principal components. The values of the vector latent are the
% variance of the columns of SCORE. 
% Hotelling's T2 is a measure of the multivariate distance of each observation 
% from the center of the data set.

%initialize output values
global globalFsDir;
global currfreq;
disp('Loading globalFsDir...')
if ~exist('globalFsDir','var')
    fprintf('globalFsDir not found, loading it...')
    eval('global globalFsDir');
    eval(['globalFsDir=' 'myp']);
end

% if strcmp(type, 'mono') == 1
% end
Allpatients = [HpatientsL,LpatientsL];
% if ~exist('electrodesL','var')
%   %electrodesL not invoked calculate the intersection of the elctrodes for
%   %all the patients
fprintf('Calculating the list of electrodes in common for all the patients...');
[electrodesL] = calculateintersectionchannels(Allpatients);
if isempty(electrodesL)
    msgerror = 'The patients have not any electrodes in common!!!END ';
    error(msgerror);
else
    fprintf('The list of electrodes in common are=');
    disp(electrodesL);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HpatientsL = {'TWH030'};
% LpatientsL = {'TWH031'};
% conditionL = {'HYP'};
% centerfrequencies = {2, 6 , 10, 23.5, 40};
% electrodesL = {'LHD1','LHD2','LHD3','LHD4'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%header of matfile that contains the corr matrix to analyze
%powerconnectivity_freq is the Spearman correlation matrix of power time
%series for every two nodes. aij = [-1 1], aii = 1
matfileid = 'powerconnectivity_freq_'; 

%Open the mat file that contains the variable for each condition
disp('Loading variables for condition H....')
populationsample = containers.Map;
pcaresults = {};
dirsolo = fullfile(globalFsDir,'figures_all_patients');
dirandname = fullfile(dirsolo, 'pcaresults.txt');
fileID = fopen(dirandname,'w');
for condind = 1:length(conditionL)
    currcond = conditionL{condind};
    figh = figure;
    for freqind= 1:length(centerfrequencies)
        currfreq = centerfrequencies{freqind};
        %populationsampleperfreq('frequency')= currfreq;
        for indpat=1:length(Allpatients)
            currentpat = Allpatients{indpat};
            fprintf('Tresting statistical relevance for cond=%s, freq=%s, pat=%s\n',conditionL{1}, num2str(currfreq),currentpat)
            %fprintf('The current patient is:%s\n',currentpat);
            patpath = strcat(globalFsDir,currentpat);
            %fprintf('The current path is %s\n', patpath);
            %fprintf('Loading the EEG object for patient  %s\n',currentpat);
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
%         disp('Calculating a t test...');
%         [h,p] = ttest2HvsL(newfftfile, HpatientsL, LpatientsL);
%         hlist{freqind} = h;plist{freqind} = p;
        disp('Calculating PCA for the p electrodes in common of the n patients...');
        caseinfo= {currcond};
        [coeff,score,latent,tsquare] = principalcomponents(newfftfile, HpatientsL, LpatientsL, caseinfo);
        %saving pca results in a cell array 
        pcaresults{condind, freqind} = {coeff,score,latent,tsquare};
        %save results in a file and print them out in screen
        writepcaresultsinfile(fileID, pcaresults,HpatientsL, LpatientsL, caseinfo, freqind);
    end
    %populationsampleperfreqandcond('data')=populationsampleperfreq;
    %print chart with results for all frequency bands
    drawttestresults(figh, HpatientsL, LpatientsL,electrodesL,centerfrequencies,currcond );
end
fclose(fileID);
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

% function [h,p] = ttest2HvsL(newfftfile, HpatientsL, LpatientsL)
% %load mat file that contain the population variables
% load(newfftfile,'-mat','populationsample', 'electrodesL')
% %separate sample1 from sample 2
% [sampleH, sampleL] = getpopulationsample(populationsample,HpatientsL, LpatientsL)
% %do the t-test
% fprintf('Calculating the t test...\n');
% %remove empty cells
% emptyCells = cellfun(@isempty,sampleL);
% sampleL(emptyCells) = [];
% emptyCells = cellfun(@isempty,sampleH);
% sampleH(emptyCells) = [];
% [h,p,ci,stats] = ttest2(horzcat(sampleH{1:end}),horzcat(sampleL{1:end}));
% fprintf('h = %d , p= %.2f', h, p);
% fprintf('ci and stats :');
% disp(ci); disp(stats);
% if h == 1
%     fprintf('T-test2 does reject the null hypothesis at the default 5\% significance level.\n');
%     fprintf('The variables are significantly different in High vs Lows\n');
% elseif h == 0
%     fprintf('T-test2 does not reject the null hypothesis at the default 5\% significance level.\n');
%     fprintf('There is no statiscal justifcation to separate H vs L\n');
% end
% 
% end


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

function [freq_band] = getfreqbandindfromind(freq)
switch freq %1 2 3 4 5 =2, 6 , 10, 23.5, 40};
    case 1
        freq_band = 2;
    case 2
        freq_band = 6;
    case 3
        freq_band = 10;
    case 4
        freq_band = 23.5;
    case 5
        freq_band = 40;
    otherwise
        errormsg ='ERROR index does not correspond with any frequency band!!';
        error(errormsg);
end
end

function [] = drawttestresults(figh, HpatientsL, LpatientsL,electrodesL, centerfrequencies, currcond)

global globalFsDir;
figure(figh);
xchannel_limit = length(electrodesL);
xlabeljump = 1;
paintedlegalready = 0;
ylabel_corrtype = 'Power-based Sp. corr.';
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
    xlabel('Channels'), ylabel(ylabel_corrtype);
    grid on ;
    if paintedlegalready == 0
        legend([HpatientsL,LpatientsL]);
        paintedlegalready = 1;
    end
    %pandhvals=  sprintf('h=%s, p=%s',num2str(hlist{indexchart}),num2str(hlist{indexchart}) );
    %text(indexchart*2, indexchart*2,pandhvals);    
end

annotation('textbox',...
    [0.6 0.05 0.5 0.10],...
    'String',{'Condition:', [currcond], 'Patients:',['H =' HpatientsL{1:end} '; L=' LpatientsL{1:end} ], [' Electrodes:' electrodesL{1:end} ]},...
    'FontSize',8,...
    'FontName','Arial',...
    'LineStyle','--',...
    'EdgeColor',[1 1 1],...
    'LineWidth',1,...
    'BackgroundColor',[0.9  0.9 0.9],...
    'Color',[0.84 0.16 0]);
% if indexchart == length(centerfrequencies)
%     highspat=  sprintf(' %s ', HpatientsL{1:end});
%     highspat= strcat('Highs:',highspat);
%     highspat= strcat(highspat, '  ');
%     lowspat=  sprintf(' %s ', LpatientsL{1:end});
%     lowspat= strcat('Lows:',lowspat);
%     allpatleg = strcat(highspat,lowspat);
%     xlabel(allpatleg)
% end
% MyBox = uicontrol('style','text')
% set(MyBox,'Hight Patients: Low Patients', HpatientsL, LpatientsL)
% xpos= 400; ypos=50; xsize= 100; ysize = 100;
% set(MyBox,'Position',[xpos,ypos,xsize,ysize])
end

function [coeff,score,latent,tsquare] = principalcomponents(newfftfile, HpatientsL, LpatientsL, caseinfo)
%principalcomponents does a PCA for corr_mat (Nxp) subjectsxvariables
%IN: corr_mat(Nxp), caseinfo(cell arryof string)
%OUT: See http://www.mathworks.com/help/stats/princomp.html?refresh=true

global globalFsDir;
if nargin < 2
    %only corr_mat as input arg.
else
    %patient, condition information included
end
%load corr matrix
disp('Loading the correlation matrix for')
disp(newfftfile)
load(newfftfile,'-mat','populationsample', 'electrodesL');
[selectedpatients] = keys(populationsample);
[selectedcorrvecs] = values(populationsample);
disp('p Patients for PCA:')
disp(selectedpatients)
disp('p Electrodes for PCA')
disp(electrodesL)
[nada n] = size(selectedcorrvecs); %n number of samples or patients
p = size(selectedcorrvecs{1}); %p number of variables
allpatients = reshape(cell2mat(selectedcorrvecs), [n,p(2)])
fprintf('PCA for n=%s samples and p=%s variables',num2str(n),num2str(p(end)));

disp('Performing PCA with standarized variables');
dataspacematrix = zscore(allpatients); %normalized data
%dataspacematrix = allpatients;
[coeff,score,latent,tsquare] = princomp(dataspacematrix);
disp('Plotting PCA Analysis');
%firstlowpatientindex firt row in the matrix with a low patient high1 ,
%high2.... low1...lown
plottingPCA(coeff,score,latent,tsquare,dataspacematrix,HpatientsL, LpatientsL);
end

function [] = plottingPCA(coeff,pc,latent,tsquare, dataspacematrix, HpatientsL, LpatientsL)
%plottingPCA plotting the PCA with the princomp results
global currfreq;
disp('The projected data on the PC componets (data in PCspace)=')
pc = pc'
disp(pc)
disp('Hotelling distance: Distance from each observation to the center of the data set')
disp(tsquare);
disp('Accumulated variance taken into account by each component, from highest to lowest')
eigenaccum = cumsum(latent)./sum(latent);
disp(eigenaccum);
firstlowpatientindex = size(HpatientsL) +1;
firstlowpatientindex = firstlowpatientindex(2);
h = figure;
%// Group 1
group1 = pc(:,1:firstlowpatientindex-1);
%// Group 2
group2 = pc(:,firstlowpatientindex:end);
%// Plot as separate colours
plot(group1(1,:), group1(2,:),'b*', group2(1,:), group2(2,:), 'r*'); 
title(['Projected power correlation per ',getgreeksymbolfreq(currfreq), ' band onto PC space'] );
legend(['Group H:' HpatientsL{1}], ['Group L:' LpatientsL{1:end}]);
xlabel('PC1') % x-axis label
ylabel('PC2') % y-axis label

% h2 = figure;
% gscatter(pc(:,1), pc(:,2), dataspacematrix', [], [], [], 'off', 'PC1', 'PC2')
% title('Projected Power correlation data onto PC space'), grid on
% disp('Save results in file')
% fileID = fopen('resultsPCA.txt','w');
% fprintf(fileID,'Condtion:%s\tFrequency:%s',conditionL{1},currfreq);
% fprintf(fileID,'%6s %12s\n','x','exp(x)');
% fprintf(fileID,'%6.2f %12.8f\n',A);
% fclose(fileID);
end
function [] =  writepcaresultsinfile(fileID, pcaresults,HpatientsL, LpatientsL, caseinfo, freqind)

formatSpec = '%1.4f\t%1.4f\t%1.4f\t%1.4f%1.4f\t%1.4f\t\n';

%if more than 1 condition {cond#,currfreq} 
[nrows,ncols] = size(pcaresults{1,freqind});
nrows = size(pcaresults{1,freqind}{4});
fprintf(fileID, ' High patients:\n');
fprintf(fileID,' %s ', HpatientsL{1,:});
fprintf(fileID, ' Low patients:\n');
fprintf(fileID,' %s ', LpatientsL{1,:});
fprintf(fileID, 'Condition(1=HYP) %s, Freq Band:%s\n',caseinfo{1},num2str(getgreeksymbolfreq(getfreqbandindfromind(freqind))));
fprintf(fileID,'\n\n');
for col = 1:ncols
    fprintf(fileID,'\n');
    if col == 1, pcar = 'coeff';
    elseif col == 2, pcar = 'score';
    elseif col == 3, pcar = 'latent';
    elseif col ==4, pcar = 'tsquare';
    end
    fprintf(fileID,'\t** %s **\n',pcar );
    fprintf(fileID,formatSpec,pcaresults{1,freqind}{col}');
    fprintf(fileID,'\n');
    fprintf('==Band= %s, PCA results for  %%s ==\n',num2str(getgreeksymbolfreq(getfreqbandindfromind(freqind))), pcar); 
    disp(pcaresults{1,freqind}{col}')
    fprintf('=================================\n\n');
end
%fclose(fileID);

end