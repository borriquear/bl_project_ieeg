function [] = displayconnectivitymatrix()
%mat file with the correlation matrix to open
%centerfrequencies = {2, 6 , 10, 23.5, 40};

matfile = 'D:\BIAL PROJECT\patients\TWH030\data\figures\powerconnectivity_freq_2_EC_POST_TWH030_11172015_s1.mat';
patientid = 'TWH030';
cond = 'EC POST';
fh = load(matfile, 'corr_matrix', 'channel_labels')
fprintf('Displaying correlation matrix in: %s\n',matfile);
[fil,col] = size(fh.corr_matrix);
% h = figure;
% % 
% corrplot(fh.corr_matrix);
% %plot variable correlation for a subset of channels(HPC)
% corrplot(corrmat([1,2,3,4,19,20,21,22],[1,2,3,4,19,20,21,22]));

h2 = figure;
imagesc(fh.corr_matrix);
colormap('jet');
colorbar;

set(gca, 'XLim', [1 fil],'XTick', [1:fil],'XTickLabel',fh.channel_labels(2:end),'YTick',[1:fil],'YTickLabel',fh.channel_labels(2:end));
xticklabel_rotate([],45,[],'Fontsize',6);
xlabel('Channels');
msgtitle = sprintf('Power-based Spearman correlation in Delta, Patient=%s, Cond=%s',patientid,cond);
title(msgtitle);

end