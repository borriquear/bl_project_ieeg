function plotwiringcost(wiring_matrices, typeobject, electrodelabel)
%plotwiringcost plot bars of wiring cost per patient in all frequencies
%IN wiring_matrices

global nbpats;
global freqlist;
global nbfreqs;
global nbconds;
if nargin > 1
    fprintf('Casting wiring_matrices...\n');
    wiring_matrices = castimplantwithlabel(wiring_matrices, typeobject, electrodelabel);
    initpat = 1;
    nbpats = size(wiring_matrices.patientslist,2);
else
    electrodelabel = {'All'}; 
    initpat = 1;
    %nbpats = 8; %nbpats is last index of patient, nbpats >= initpat <= totalpats=size(wiring_matrices.patientslist,2);
    nbpats = size(wiring_matrices.patientslist,2);
end

freqlist = wiring_matrices.frequencylist;
%freqlist = [freqlist(2),freqlist(4), freqlist(5), freqlist(7), freqlist(7), freqlist(8)]
nbfreqs = length(freqlist);
nbconds = size( wiring_matrices.conditionslist, 2); nbconds = 2;
initfreq = 3;
initcond =1;
listofselectedpats = wiring_matrices.patientslist(initpat:nbpats);

%nbfreqs = length( wiring_matrices.frequencylist);
barmatrix_ispc = zeros(nbpats, nbconds*nbfreqs);
barmatrix_pli = zeros(nbpats, nbconds*nbfreqs);
barmatrix_power = zeros(nbpats, nbconds*nbfreqs);
barmatrix_chg_ispc = zeros(nbpats, nbconds*nbfreqs);
barmatrix_chg_pli = zeros(nbpats, nbconds*nbfreqs);
barmatrix_chg_power = zeros(nbpats, nbconds*nbfreqs);
%cols = cond*freq

barmatrix_ispc = zeros(nbpats-initpat+1, (nbconds - initcond +1)*(nbfreqs-initfreq+1));
barmatrix_pli = zeros(nbpats-initpat+1, (nbconds - initcond +1)*(nbfreqs-initfreq+1));
barmatrix_power = zeros(nbpats-initpat+1, (nbconds - initcond +1)*(nbfreqs-initfreq+1));
for p=initpat:nbpats
    for ic=initcond:nbconds
        for f=initfreq:nbfreqs
            barmatrix_ispc(p-initpat+1,(ic-1)*nbfreqs + f)=  mean2(wiring_matrices.wiring_ispc{p,ic,f});%closed
            barmatrix_pli(p-initpat+1,(ic-1)*nbfreqs + f)=  mean2(wiring_matrices.wiring_pli{p,ic,f});%closed
            barmatrix_power(p-initpat+1,(ic-1)*nbfreqs + f)=  mean2(wiring_matrices.wiring_power{p,ic,f});%closed
        end
    end
end

%plot 1 chart with change in R, PLI Power P- P*F for cond,frequency
sti = sprintf('Wiring Cost (P*.F) F=R|PLI|Power %d in conditions (%s ...%s) ROIS=%s', nbconds, wiring_matrices.conditionslist{1}, wiring_matrices.conditionslist{end}, electrodelabel{1} );
barf = figure;
subplot(1,3,1)
%bar(barmatrix_chg_ispc);
bar(barmatrix_ispc);
xlabel('Frequencies x Conds Patients '), ylabel('R wiring cost P.*F)')
ax = set(gca,'XTick',1:length(listofselectedpats),'XTickLabel', listofselectedpats,'Fontsize',7,'XTickLabelRotation', 45);%xticklabel_rotate([],45,[],'Fontsize',6);
%ax.XTickLabelRotation=45;
subplot(1,3,2)
%bar(barmatrix_chg_pli);
bar(barmatrix_pli);
xlabel('Frequencies x Conds Patients'), ylabel('PLI wiring cost  P.*F')
ax = set(gca,'XTick',1:length(listofselectedpats),'XTickLabel', listofselectedpats,'Fontsize',7,'XTickLabelRotation', 45);%xticklabel_rotate([],45,[],'Fontsize',6);
title(sti, 'Interpreter', 'none');
subplot(1,3,3)
%bar(barmatrix_chg_power);
bar(barmatrix_power);
ax = set(gca,'XTick',1:length(listofselectedpats),'XTickLabel', listofselectedpats,'Fontsize',7,'XTickLabelRotation', 45);%xticklabel_rotate([],45,[],'Fontsize',6);
xlabel('Frequencies x Conds Patients '), ylabel('Power wiring cost P.*F')
%plot the mean of all patients
ispc_diff = zeros(1,nbfreqs ); pli_diff = zeros(1, nbfreqs ); power_diff = zeros(1,nbfreqs );
for fqs=initfreq:nbfreqs
    ispc_diff(fqs) = mean(barmatrix_ispc(:, fqs)) - mean(barmatrix_ispc(:, nbfreqs + fqs ));
    pli_diff(fqs) = mean(barmatrix_pli(:, fqs)) - mean(barmatrix_pli(:, nbfreqs + fqs ));
    power_diff(fqs) = mean(barmatrix_power(:, fqs)) - mean(barmatrix_power(:, nbfreqs + fqs ));    
end
figure;
%R measure
subplot(1,3,1)
for fqs=initfreq:nbfreqs
    if fqs ==5 
        %delta band
        bar(fqs, ispc_diff(fqs), 'EdgeColor', [0 .9 .9], 'FaceColor', [0 .5 .5]);
    else
       bar(fqs, ispc_diff(fqs));
    end
    hold on
end
set(gca,'XTick', [initfreq:nbfreqs], 'XTickLabel', round(freqlist(initfreq:end),1,'significant'));
ylabel('Patients'), xlabel('Frequency bands')
ts = sprintf('Mean P*R_C-P*R_O ROI=%s', electrodelabel{1});
%legend('delta', 'theta', 'alpha', 'beta1','beta2','gamma');
title(ts);

%PLI measure
subplot(1,3,2)
for fqs=initfreq:nbfreqs
    if fqs ==5 
        %delta band
        bar(fqs, pli_diff(fqs), 'EdgeColor', [0 .9 .9], 'FaceColor', [0 .5 .5]);
    else
       bar(fqs, pli_diff(fqs));
    end
    hold on
end
set(gca,'XTick', [initfreq:nbfreqs], 'XTickLabel', round(freqlist(initfreq:end),1,'significant'));
ylabel('Patients'), xlabel('Frequency bands')
ts = sprintf('Mean P*PLI_C-P*PLI_O ROI=%s', electrodelabel{1});
%legend('delta', 'theta', 'alpha', 'beta1','beta2','gamma');
title(ts);

%Power measure
subplot(1,3,3)
for fqs=initfreq:nbfreqs
    if fqs ==5 
        %delta band
        bar(fqs, power_diff(fqs), 'EdgeColor', [0 .9 .9], 'FaceColor', [0 .5 .5]);
    else
       bar(fqs, power_diff(fqs));
    end
    hold on
end
set(gca,'XTick', [initfreq:nbfreqs], 'XTickLabel', round(freqlist(initfreq:end),1,'significant'));
ylabel('Patients'), xlabel('Frequency bands')
ts = sprintf('Mean P*Pw_C-P*Pw_O ROI=%s', electrodelabel{1});
legend('delta', 'theta', 'alpha', 'beta1','beta2','gamma');
title(ts);

%calculate cond1-cond2, condn-1- condn
if nbconds == 2
    pairlabellist= {'ECPRE-EOPRE'};
elseif nbconds == 3
    pairlabellist= {'ECPRE-EOPRE', 'ECPRE-HYP','EOPRE-HYP'};
end

for indexp=1:length(pairlabellist)
    pairlabel = pairlabellist{indexp};
    fprintf('Calculating distance for %s\n',pairlabel );
    [barmatrix_pli_df,barmatrix_ispc_df,barmatrix_power_df ]= calcPFdistance2conds(nbpats,initpat,nbfreqs,initfreq,barmatrix_pli,barmatrix_ispc,barmatrix_power,pairlabel,electrodelabel );
    fprintf('Ploting distance for %s\n',pairlabel );
    plotPFdistance2conds(barmatrix_pli_df,barmatrix_ispc_df,barmatrix_power_df,pairlabel,listofselectedpats,electrodelabel);
end
% fprintf('Plotting the scatter of physical and functional distance\n');
% scatterPF(wiring_matrices);
end

function [barmatrix_pli_df,barmatrix_ispc_df,barmatrix_power_df ]= calcPFdistance2conds(nbpats,initpat,nbfreqs,initfreq,barmatrix_pli,barmatrix_ispc,barmatrix_power, pairlabel,electrodelabel )

barmatrix_ispc_df = zeros(nbpats-initpat+1, 1*nbfreqs);
barmatrix_pli_df = zeros(nbpats-initpat+1, 1*nbfreqs);
barmatrix_power_df = zeros(nbpats-initpat+1, 1*nbfreqs);
if strcmp('ECPRE-EOPRE',pairlabel) == 1
    pos = 0; neg = nbfreqs ; %ec - eo 1 ...- 9...(for 8 greqs)
elseif strcmp('ECPRE-HYP',pairlabel) == 1 %ec - hyp 19....-1....
    pos = 0; neg = 2*nbfreqs ;
elseif strcmp('EOPRE-HYP',pairlabel) == 1 %eo - hyp 9 - 19
    pos =nbfreqs ; neg = 2*nbfreqs ;
else
    fprintf('ERROR: label %s NOT found', pairlabel);
end
for p=initpat:nbpats
    for f=initfreq:nbfreqs
        %EC - EO
        %difference in wiring cost respect to physical distance
        %wiring cost respect
        barmatrix_pli_df(p-initpat+1,f) = barmatrix_pli(p-initpat+1,pos+f) - barmatrix_pli(p-initpat+1,neg+f);
        barmatrix_ispc_df(p-initpat+1,f) = barmatrix_ispc(p-initpat+1,pos+f)- barmatrix_ispc(p-initpat+1,neg+f);
        barmatrix_power_df(p-initpat+1,f) = barmatrix_power(p-initpat+1,pos+f)- barmatrix_power(p-initpat+1,neg+f);
    end
end
end

function [] = plotPFdistance2conds(barmatrix_pli_df,barmatrix_ispc_df,barmatrix_power_df,pairlabel,patientslist,electrodelabel)
global freqlist;
sti = sprintf('Difference in Wiring Cost for %s ROIS=%s', pairlabel,electrodelabel{1});
barfdf = figure;
subplot(1,3,1)
%bar(barmatrix_chg_pli_df);
bar(barmatrix_pli_df(:,3:end));
xlabel('Frequencies Patients'), ylabel('Change in PLI wiring cost P*(F1-F2)')
ax = set(gca,'xtick',1:length(patientslist),'XTickLabel', patientslist,'Fontsize',7,'XTickLabelRotation', 45);%xticklabel_rotate([],45,[],'Fontsize',6);

subplot(1,3,2)
%bar(barmatrix_chg_ispc_df);
bar(barmatrix_ispc_df(:,3:end));
xlabel('Frequencies Patients '), ylabel('Change in R wiring cost  P*(F1-F2)')
ax = set(gca,'xtick',1:length(patientslist),'XTickLabel', patientslist,'Fontsize',7,'XTickLabelRotation', 45);%xticklabel_rotate([],45,[],'Fontsize',6);
title(sti, 'Interpreter', 'none');
subplot(1,3,3)
%bar(barmatrix_chg_power_df);
bar(barmatrix_power_df(:,3:end));
xlabel('Frequency-Patients EC-EO'), ylabel('Change in Power wiring cost  P*(F1-F2)')
ax = set(gca,'xtick',1:length(patientslist),'XTickLabel', patientslist,'Fontsize',7,'XTickLabelRotation', 45);%xticklabel_rotate([],45,[],'Fontsize',6);
legend('delta', 'theta', 'alpha', 'beta1','beta2','gamma');
% plot imagesc
initfreq = 3;
nbfreqs = size(barmatrix_ispc_df(:,3:end),2);
freqlist  = logspace(log10(1),log10(50),8);
freqlist = freqlist(initfreq:end);
imgf = figure;
subplot(1,3,1)
imagesc(barmatrix_pli_df(:,initfreq:end));
[hStrings,textColors ] = plotvaluesinimageesc(barmatrix_pli_df(:,initfreq:end));
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
colorbar;
set(gca,'XTick', [1:nbfreqs], 'XTickLabel', round(freqlist(1:end),1,'significant'), 'YTickLabel',patientslist,'ytick', 1:length(patientslist) );
%set(gca,'XTickLabel', round(freqlist(initfreq:end),1,'significant'));
ylabel('Patients'), xlabel('Frequency bands')
ts = sprintf('Wiring Cost %s, PLI based ROIS=%s', pairlabel, electrodelabel{1});
title(ts);

subplot(1,3,2);
imagesc(barmatrix_ispc_df(:,initfreq:end));
[hStrings,textColors ] = plotvaluesinimageesc(barmatrix_ispc_df(:,initfreq:end));
colorbar
set(gca,'XTick', [1:nbfreqs], 'XTickLabel', round(freqlist(1:end),1,'significant'), 'YTickLabel',patientslist,'ytick', 1:length(patientslist));
ylabel('Patients'), xlabel('Frequency bands')
ts = sprintf('Wiring Cost %s, R based', pairlabel);
title(ts);
subplot(1,3,3)
imagesc(barmatrix_power_df(:,initfreq:end));
[hStrings,textColors ] = plotvaluesinimageesc(barmatrix_power_df(:,initfreq:end));
colorbar;
set(gca,'XTick', [1:nbfreqs], 'XTickLabel', round(freqlist(1:end),1,'significant'), 'YTickLabel',patientslist,'ytick', 1:length(patientslist) );
ylabel('Patients'), xlabel('Frequency bands')
ts = sprintf('Wiring Cost %s, Power based ROIS=%s', pairlabel,electrodelabel{1});
title(ts);

%plot binary matrices +/-
bin = im2bw(barmatrix_pli_df(:,initfreq:end), 0.0);
bin2 = im2bw(barmatrix_ispc_df(:,initfreq:end), 0.0);
bin3 = im2bw(barmatrix_power_df(:,initfreq:end), 0.0);
imgfbin = figure;
subplot(1,3,1)
imagesc(bin);
set(gca,'XTick', [1:nbfreqs], 'XTickLabel', round(freqlist(1:end),1,'significant'), 'YTickLabel',patientslist,'ytick', 1:length(patientslist) );
ylabel('Patients'), xlabel('Frequency bands')
ts = sprintf('+/- wiring Cost %s, PLI based ROIS=%s',pairlabel,electrodelabel{1});
title(ts)

subplot(1,3,2)
imagesc(bin2);
set(gca,'XTick', [1:nbfreqs], 'XTickLabel', round(freqlist(1:end),1,'significant'), 'YTickLabel',patientslist,'ytick', 1:length(patientslist) );
ylabel('Patients'), xlabel('Frequency bands')
ts = sprintf('+/- wiring Cost %s, R  based ROIS=%s',pairlabel,electrodelabel{1});
title(ts)

subplot(1,3,3)
imagesc(bin3);
set(gca,'XTick', [1:nbfreqs], 'XTickLabel', round(freqlist(1:end),1,'significant'), 'YTickLabel',patientslist,'ytick', 1:length(patientslist) );
ylabel('Patients'), xlabel('Frequency bands')
ts = sprintf('+/- wiring Cost %s, Power based ROIS=%s',pairlabel, electrodelabel{1});
title(ts)
end

function [hStrings,textColors ] = plotvaluesinimageesc(mat)
%plotvaluesinimageesc plot imagesc with the numeric values
nfs = size(mat,2); npats =size(mat,1);
textStrings = num2str(mat(:),'%0.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding

%[x,y] = meshgrid(1:size(mat,2),1:11); %# Create x and y coordinates for the strings
[x,y] = meshgrid(1:size(mat,2),1:size(mat,1));
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
    'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
textColors = repmat(mat(:) > midValue,1,3);  %# Choose white or black for the
%#   text color of the strings so
%#   they can be easily seen over
%#   the background color
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
end
