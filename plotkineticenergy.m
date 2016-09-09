function kinetic_matrices = plotkineticenergy(wiring_matrices)
%plotkineticenergy plot bars of kinetic energy (P*F*frequency)^2
%the wise multiplication of the physical distance and the functional
%connectivity multiplied by the Hertz, gives us distance/s which is the
%velocity.
% 1/2v^2 is the kinetic energy to transport a punctual mass m=1 between two
% points, electrode positions.
%IN wiring_matrices

global nbpats;
global freqlist;
global nbfreqs;
global nbconds;
globalFsDir = loadglobalFsDir();
kinetic_matrices = struct ;
nbpats = size(wiring_matrices.patientslist,2);
freqlist = wiring_matrices.frequencylist;
nbfreqs = length(freqlist);
nbconds = size(wiring_matrices.conditionslist, 2);
%nbfreqs = length( wiring_matrices.frequencylist);
barmatrix_ispc = zeros(nbpats, nbconds*nbfreqs);
barmatrix_pli = zeros(nbpats, nbconds*nbfreqs);
barmatrix_power = zeros(nbpats, nbconds*nbfreqs);
barmatrix_chg_ispc = zeros(nbpats, nbconds*nbfreqs);
barmatrix_chg_pli = zeros(nbpats, nbconds*nbfreqs);
barmatrix_chg_power = zeros(nbpats, nbconds*nbfreqs);
%calculate the kinetic matrices
initfreq= 1;
wiring_ispc ={};
%wiring_ispc = p*f
for p=1:nbpats
    for ic=1:nbconds
        for f=initfreq:nbfreqs
            wiring_ispc{p,ic,f} = wiring_matrices.wiring_ispc{p,ic,f}*wiring_matrices.frequencylist(f);
            wiring_ispc{p,ic,f}= 0.5*(wiring_ispc{p,ic,f})^2;
            wiring_pli{p,ic,f} = wiring_matrices.wiring_pli{p,ic,f}*wiring_matrices.frequencylist(f);
            wiring_pli{p,ic,f}= 0.5*(wiring_pli{p,ic,f})^2
            wiring_power{p,ic,f} = wiring_matrices.wiring_power{p,ic,f}*wiring_matrices.frequencylist(f);
            wiring_power{p,ic,f}= 0.5*(wiring_power{p,ic,f})^2
        end
    end
end
kinetic_matrices.wiring_ispc = wiring_ispc;
kinetic_matrices.wiring_pli = wiring_pli;
kinetic_matrices.wiring_power = wiring_power;
kinetic_matrices.patientslist = wiring_matrices.patientslist;
kinetic_matrices.conditionslist = wiring_matrices.conditionslist;
kinetic_matrices.frequencylist = wiring_matrices.frequencylist;
kinetic_matrices.distMatrixcell = wiring_matrices.distMatrixcell;
%savethe kinetic_matrices structure
matfilename = 'wiringcost_matrices.mat';
matfilename = fullfile(globalFsDir, matfilename);
disp(kinetic_matrices)
save(matfilename,'kinetic_matrices');

for p=1:nbpats
    for ic=1:nbconds
        for f=initfreq:nbfreqs
            barmatrix_ispc(p,(ic-1)*nbfreqs + f)=  mean2(kinetic_matrices.wiring_ispc{p,ic,f});%closed
            barmatrix_pli(p,(ic-1)*nbfreqs + f)=  mean2(kinetic_matrices.wiring_pli{p,ic,f});%closed
            barmatrix_power(p,(ic-1)*nbfreqs + f)=  mean2(kinetic_matrices.wiring_power{p,ic,f});%closed
            %             barmatrix_chg_ispc(p,f*nbconds +(ic - nbconds)) = mean2(kinetic_matrices.distMatrixcell{p}) - mean2(kinetic_matrices.wiring_ispc{p,ic,f});
            %             barmatrix_chg_pli(p,f*nbconds +(ic - nbconds)) = mean2(kinetic_matrices.distMatrixcell{p}) - mean2(kinetic_matrices.wiring_pli{p,ic,f});
            %             barmatrix_chg_power(p,f*nbconds +(ic - nbconds))= mean2(kinetic_matrices.distMatrixcell{p}) - mean2(kinetic_matrices.wiring_power{p,ic,f});
        end
    end
end
%plot 1 chart with K = 0.5 (P.*F)^2
sti = sprintf('Kinetic Energy(P*.F*frequency) 1/2 (distance*Hz)^2 F=R|PLI|Power %d in conditions (%s ...%s)', nbconds, kinetic_matrices.conditionslist{1},  kinetic_matrices.conditionslist{end} );
barf = figure;
subplot(1,3,1)
%bar(barmatrix_chg_ispc);
bar(barmatrix_ispc);
xlabel('Conds x Frequencies Patients '), ylabel('Kinetic Energy for w=R. K = 1/2(w*P*f)^2 ')
ax = set(gca,'XTickLabel', kinetic_matrices.patientslist,'Fontsize',7,'XTickLabelRotation', 45);%xticklabel_rotate([],45,[],'Fontsize',6);
%ax.XTickLabelRotation=45;
subplot(1,3,2)
%bar(barmatrix_chg_pli);
bar(barmatrix_pli);
xlabel('Conds x Frequencies Patients '), ylabel('Kinetic Energy for w=PLI. K = 1/2(w*P*f)^2 ')
ax = set(gca,'XTickLabel', kinetic_matrices.patientslist,'Fontsize',7,'XTickLabelRotation', 45);%xticklabel_rotate([],45,[],'Fontsize',6);
title(sti, 'Interpreter', 'none');
subplot(1,3,3)
%bar(barmatrix_chg_power);
bar(barmatrix_power);
ax = set(gca,'XTickLabel', kinetic_matrices.patientslist,'Fontsize',7,'XTickLabelRotation', 45);%xticklabel_rotate([],45,[],'Fontsize',6);
xlabel('Conds xFrequencies Patients '), ylabel('Kinetic Energy for w=Power corr. K = 1/2(pow*P*f)^2 ')

%calculate cond1-cond2, condn-1- condn
pairlabellist= {'ECPRE-EOPRE', 'ECPRE-HYP','EOPRE-HYP'};
for indexp=1:length(pairlabellist)
    pairlabel = pairlabellist{indexp};
    fprintf('Calculating distance for %s\n',pairlabel );
    [barmatrix_pli_df,barmatrix_ispc_df,barmatrix_power_df ]= calcPFdistance2conds(nbpats,nbfreqs,barmatrix_pli,barmatrix_ispc,barmatrix_power,pairlabel );
    fprintf('Ploting distance for %s\n',pairlabel );
    plotPFdistance2conds(barmatrix_pli_df,barmatrix_ispc_df,barmatrix_power_df,pairlabel,kinetic_matrices.patientslist );
end
% fprintf('Plotting the scatter of physical and functional distance\n');
% scatterPF(kinetic_matrices);
end

function [barmatrix_pli_df,barmatrix_ispc_df,barmatrix_power_df ]= calcPFdistance2conds(nbpats,nbfreqs,barmatrix_pli,barmatrix_ispc,barmatrix_power, pairlabel )

barmatrix_ispc_df = zeros(nbpats, 1*nbfreqs);
barmatrix_pli_df = zeros(nbpats, 1*nbfreqs);
barmatrix_power_df = zeros(nbpats, 1*nbfreqs);
if strcmp('ECPRE-EOPRE',pairlabel) == 1
    pos = 0; neg = nbfreqs ; %ec - eo 1 ...- 9...(for 8 greqs)
elseif strcmp('ECPRE-HYP',pairlabel) == 1 %ec - hyp 19....-1....
    pos = 0; neg = 2*nbfreqs ;
elseif strcmp('EOPRE-HYP',pairlabel) == 1 %eo - hyp 9 - 19
    pos =nbfreqs ; neg = 2*nbfreqs ;
else
    fprintf('ERROR: label %s NOT found', pairlabel);
end
for p=1:nbpats
    for f=1:nbfreqs
        %EC - EO
        %difference in wiring cost respect to physical distance
        %wiring cost respect
        barmatrix_pli_df(p,f) = barmatrix_pli(p,pos+f) - barmatrix_pli(p,neg+f);
        barmatrix_ispc_df(p,f) = barmatrix_ispc(p,pos+f)- barmatrix_ispc(p,neg+f);
        barmatrix_power_df(p,f) = barmatrix_power(p,pos+f)- barmatrix_power(p,neg+f);
    end
end
end

function [] = plotPFdistance2conds(barmatrix_pli_df,barmatrix_ispc_df,barmatrix_power_df,pairlabel,patientslist )

sti = sprintf('Difference in Kinetic Energy for %s', pairlabel);
barfdf = figure;
subplot(1,3,1)
%bar(barmatrix_chg_pli_df);
bar(barmatrix_pli_df(:,3:end));
xlabel('CondsxFrequencies Patients'), ylabel('Change in K P*(PLI1-PLI2)')
ax = set(gca,'XTickLabel', patientslist,'Fontsize',7,'XTickLabelRotation', 45);%xticklabel_rotate([],45,[],'Fontsize',6);

subplot(1,3,2)
%bar(barmatrix_chg_ispc_df);
bar(barmatrix_ispc_df(:,3:end));
xlabel('CondsxFrequencies Patients '), ylabel('Change in K P*(R1-R2)')
ax = set(gca,'XTickLabel', patientslist,'Fontsize',7,'XTickLabelRotation', 45);%xticklabel_rotate([],45,[],'Fontsize',6);
title(sti, 'Interpreter', 'none');
subplot(1,3,3)
%bar(barmatrix_chg_power_df);
bar(barmatrix_power_df(:,3:end));
xlabel('Frequency-Patients EC-EO'), ylabel('Change in K P*(Pw1-Pw2)')
ax = set(gca,'XTickLabel', patientslist,'Fontsize',7,'XTickLabelRotation', 45);%xticklabel_rotate([],45,[],'Fontsize',6);

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
set(gca,'XTick', [1:nbfreqs], 'XTickLabel', round(freqlist(1:end),1,'significant'), 'YTickLabel',patientslist );
%set(gca,'XTickLabel', round(freqlist(initfreq:end),1,'significant'));
ylabel('Patients'), xlabel('Frequency bands')
ts = sprintf('K %s, PLI based', pairlabel);
title(ts);
subplot(1,3,2);
imagesc(barmatrix_ispc_df(:,initfreq:end));
[hStrings,textColors ] = plotvaluesinimageesc(barmatrix_ispc_df(:,initfreq:end));
colorbar
set(gca,'XTick', [1:nbfreqs], 'XTickLabel', round(freqlist(1:end),1,'significant'), 'YTickLabel',patientslist );
ylabel('Patients'), xlabel('Frequency bands')
ts = sprintf('K %s, R based', pairlabel);
title(ts);
subplot(1,3,3)
imagesc(barmatrix_power_df(:,initfreq:end));
[hStrings,textColors ] = plotvaluesinimageesc(barmatrix_power_df(:,initfreq:end));
colorbar;
set(gca,'XTick', [1:nbfreqs], 'XTickLabel', round(freqlist(1:end),1,'significant'), 'YTickLabel',patientslist );
ylabel('Patients'), xlabel('Frequency bands')
ts = sprintf('K %s, Power based', pairlabel);
title(ts);
%plot relationship between  and percentage of power. 2 gropus, 10 data
%points one for eaach gruoup.



%plot binary matrices +/-
bin = im2bw(barmatrix_pli_df(:,initfreq:end), 0.0);
bin2 = im2bw(barmatrix_ispc_df(:,initfreq:end), 0.0);
bin3 = im2bw(barmatrix_power_df(:,initfreq:end), 0.0);
imgfbin = figure;
subplot(1,3,1)
imagesc(bin);
set(gca,'XTick', [1:nbfreqs], 'XTickLabel', round(freqlist(1:end),1,'significant'), 'YTickLabel',patientslist );
ylabel('Patients'), xlabel('Frequency bands')
ts = sprintf('+/- K %s, PLI based',pairlabel );
title(ts)
subplot(1,3,2)
imagesc(bin2);
set(gca,'XTick', [1:nbfreqs], 'XTickLabel', round(freqlist(1:end),1,'significant'), 'YTickLabel',patientslist );
ylabel('Patients'), xlabel('Frequency bands')
ts = sprintf('+/- K %s, R  based',pairlabel );
title(ts)
subplot(1,3,3)
imagesc(bin3);
set(gca,'XTick', [1:nbfreqs], 'XTickLabel', round(freqlist(1:end),1,'significant'), 'YTickLabel',patientslist );
ylabel('Patients'), xlabel('Frequency bands')
ts = sprintf('+/- K %s, Power based',pairlabel );
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