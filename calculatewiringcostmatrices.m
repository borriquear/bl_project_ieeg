function [wiring_matrices] = calculatewiringcostmatrices()
%calculatewiringcostmatrices calculate the wiring cost matrix and represents the
%results w = pHYS*fUNC
%conditionsL = {'EC_PRE', 'EC_POST'};
%[srate, min_freq, max_freq, num_frex, time, n_wavelet, half_of_wavelet_size, frex, s, wavelet_cycles] =   initialize_wavelet()
%freqs2use  = logspace(log10(min_freq),log10(max_freq),8);
globalFsDir = loadglobalFsDir();
%matfilename = fullfile(globalFsDir, 'phaseconn_matrices.mat');
matfilename = fullfile(globalFsDir, 'phaseconn_matrices_tw.mat');
fh = load(matfilename);
patientslist = fh.phaseconn_matrix.patientsl;
conditionslist = fh.phaseconn_matrix.conditionsl;
frequencylist = fh.phaseconn_matrix.freqsl;
%open power matrix mat file
%fh2 = load(fullfile(globalFsDir, 'powerconn_matrices.mat'));
fh2 = load(fullfile(globalFsDir, 'powerconn_matrices_tw.mat'));
wiring_matrices = struct ;
label = 'All';
%global initfreqindex;initfreqindex = 3;
%the distMatrix is casted with the vector of electrodes indices [2 3 4 ...66]
%patientsL  = { 'TWH030','TWH033','TWH037','TWH038','TWH042'};
%patientslist = {'TWH047'};
ispc_matrix_freqs = []; pli_matrix_freqs = []; power_matrix_freqs= []; power_matrix= [];
conditionsize= size(conditionslist,2);%-1; % EC EO HYP
for i=1:size(patientslist,2)
    channelclassindexes = getindexesfromlabel(patientslist{i}, label);
    channelclassindexes = channelclassindexes -1;
    %calculate euclidean matrix
    distMatrix = calculateEuclideandistancematrix(patientslist{i}, label);
    distMatrixcell{i} = distMatrix;
    fprintf('Calculated the (P) Euclidean matrix for %s  %s\n',patientslist{i}, label);
    for j=1:conditionsize
        ispc_matrix_freqs = fh.phaseconn_matrix.ispc_matrix{i,j};% [ 36x36x8  double]    [ 36x36x8  double]
        pli_matrix_freqs = fh.phaseconn_matrix.pli_matrix{i,j};
        %power_matrix_freqs = fh2.powerconn_matrix.power_matrix{i,j};
        icoh_matrix_freqs = fh.phaseconn_matrix.icoh_matrix{i,j};
        %initfreqindex = 3; %dont take 1 and 1.7  freqs
        for k=1:size(frequencylist,2)
            %matrix [channelsxchannels] per each frequency
            fprintf('Calculating (F) Functional matrices for %s, chanels= %s cond=%s Freq = %d...\n',patientslist{i}, label, conditionslist{j}, k);
            power_matrix_freqs = fh2.powerconn_matrix.power_matrix(:,:,k);
            power_matrix = power_matrix_freqs{i,j};
            power_matrix = power_matrix(channelclassindexes,channelclassindexes);
            %power_matrix_freqs  = power_matrix_freqs(:,:,k);
            ispc_matrix = ispc_matrix_freqs(:,:,k);
            % cast for Label type electrodes
            ispc_matrix= ispc_matrix(channelclassindexes,channelclassindexes);
            pli_matrix = pli_matrix_freqs(:,:,k);
            pli_matrix = pli_matrix(channelclassindexes,channelclassindexes);
            %icoh_matrix = icoh_matrix_freqs(:,:,k);
            wiring_ispc{i,j,k} = distMatrix.*ispc_matrix;
            wiring_pli{i,j,k} = distMatrix.*pli_matrix;
            %wiring_icoh{i,j,k} = distMatrix.*icoh_matrix;
            wiring_power{i,j,k} = distMatrix.*power_matrix;
        end
            fprintf('Mean wiring cost of R matrix  for Patient %s %s =%.2f\n',patientslist{i},conditionslist{j},mean(nonzeros(ispc_matrix(:))));
            fprintf('Mean wiring cost of PLI matrix  for Patient %s %s =%.2f\n',patientslist{i},conditionslist{j}, mean(nonzeros(pli_matrix(:))));
            fprintf('Mean wiring cost of Power matrix  for Patient %s %s =%.2f\n',patientslist{i},conditionslist{j}, mean(nonzeros(power_matrix(:))));
    end
end
wiring_matrices.wiring_ispc = wiring_ispc;
wiring_matrices.wiring_pli = wiring_pli;
wiring_matrices.wiring_power = wiring_power;
wiring_matrices.patientslist = patientslist;
wiring_matrices.conditionslist = conditionslist;
wiring_matrices.frequencylist = frequencylist;
wiring_matrices.distMatrixcell = distMatrixcell;
%savethe wiring_matrices structure
matfilename = 'wiringcost_matrices.mat';
matfilename = fullfile(globalFsDir, matfilename);
disp(wiring_matrices)
save(matfilename,'wiring_matrices');
end


% 
% function [allmetricslist] = plotwiringcostmatrices(wiring_matrices)
% % chose patient cond and frequency to plot
% patl = wiring_matrices.patientslist;
% %patl = {'TWH030'}
% condl= wiring_matrices.conditionslist;
% freql = wiring_matrices.frequencylist;
% allmetricslist = {};
% % i = 1; %pat
% % j = 1; %cond
% % k = 1; %freq
% for i=1:size(patl,2)
%     p = patl{i};
%     for j=1:size(condl,2)
%         c = condl{j};
%         for k=1:size(freql,2)
%             f = freql(k);
%             wiringcost_ispc = wiring_matrices.wiring_ispc{i,j,k};
%             wiringcost_pli = wiring_matrices.wiring_pli{i,j,k};
%             wiringcost_power = wiring_matrices.wiring_power{i,j,k};
%             phys_distance = wiring_matrices.distMatrixcell{i};
%             f1 = figure;
%             imagesc(phys_distance);
%             msgtitle = sprintf('Euclidean Distance BiTemp electrodes for Patient = %s',p);
%             title(msgtitle);
%             f2 = figure;
%             imagesc(wiringcost_ispc);
%             msgtitle = sprintf('Wiring Cost= P.*R for Patient = %s Cond=%s Freq=%s ',p,c,num2str(f));
%             title(msgtitle);
%             f3 = figure;
%             imagesc(wiringcost_pli);
%             msgtitle = sprintf('Wiring Cost= P.*PLI for Patient = %s Cond=%s Freq=%s ',p,c,num2str(f));
%             title(msgtitle);
%             colorbar;
%             f4 = figure;
%             %wiringchange due to connectivity;
%             wiring_change = phys_distance - wiringcost_ispc;
%             imagesc(wiring_change);
%             msgtitle = sprintf('Change in Wiring Cost due to phase connectivity (R) for Patient = %s Cond=%s Freq=%s ',p,c,num2str(f));
%             title(msgtitle);
%             colorbar;
%             fprintf('Patient = %s, Cond=%s, Freq=%s. Mean Wiring cost D=%.2f, D.*R=%.2f, D.*PLI=%.2f \n',p,c,num2str(f),sum(phys_distance)/sum(phys_distance~=0), sum(wiringcost_ispc)/sum(wiringcost_ispc~=0),sum(wiringcost_pli)/sum(wiringcost_pli~=0))
%             
%             
%             %wiringchange due to power;
%             f5 = figure; 
%             imagesc(wiringcost_power);
%             msgtitle = sprintf('Wiring Cost due to Power connectivity for Patient = %s Cond=%s Freq=%s ',p,c,num2str(f));
%             title(msgtitle);
%             colorbar;
%             fprintf('Patient = %s, Cond=%s, Freq=%s. Mean Wiring cost D=%.2f, D.*R=%.2f, D.*PLI=%.2f \n',p,c,num2str(f),sum(phys_distance)/sum(phys_distance~=0), sum(wiringcost_power)/sum(wiringcost_power~=0),sum(wiringcost_power)/sum(wiringcost_power~=0))
%             
%             
%             f6 = figure; %wiringchange due to power;
%             wiring_change = phys_distance - wiringcost_power;
%             imagesc(wiring_change);
%             msgtitle = sprintf('Change in Wiring Cost due to Power connectivity for Patient = %s Cond=%s Freq=%s ',p,c,num2str(f));
%             title(msgtitle);
%             colorbar;
%             fprintf('Patient = %s, Cond=%s, Freq=%s. Mean Wiring cost D=%.2f, D.*R=%.2f, D.*PLI=%.2f \n',p,c,num2str(f),sum(phys_distance)/sum(phys_distance~=0), sum(wiringcost_ispc)/sum(wiringcost_ispc~=0),sum(wiringcost_pli)/sum(wiringcost_pli~=0))
%             
%             
%             
%             %R=16.84, PLI=14.03
%             %plot the network
%             corrMatrix = wiringcost_ispc;
%             nstds = 2; %number of stds, threshold = mean +sntds*std
%             [threshold, corrMatrix] = calculatethresholdmatrix(corrMatrix, nstds);
%             [myfullname, EEG, channel_labels, patdate, patsession] =  initialize_EEG_variables(p,c)
%             strNames = channel_labels(2:end);
%             legendofmatrices = {c, p, f, channel_labels};
%             f7 = figure;
%             myColorMap = lines(length(corrMatrix));
%             circularGraph(corrMatrix,'Colormap',myColorMap,'Label',strNames);
%             %title(['Condition: ',currcond,', Frequency:', getgreeksymbolfreq(currfreq), ', Patient:',currentpat, ', Thres. = m+ n*std=', num2str(threshold)]);
%             title([c,',in', num2str(f), ', Patient:',p, ', Thresh. = m+ 1*std']);
%             drawnow update
%             %calculate network metrics for corrMatrix
%             [allmetrics] = calculategraphmatrixmetrics(corrMatrix, legendofmatrices);
%             allmetricslist{i,j,k} = {legendofmatrices, allmetrics, wiringcost_ispc, wiring_change, wiringcost_pli, phys_distance, wiringcost_power};
%             %close figures
%             figure(f1);pause(1);close(f1);
%             figure(f2); pause(1);close(f2);
%             figure(f3);pause(1);close(f3);figure(f4); pause(1);close(f4);figure(f5);pause(1);close(f5);figure(f6);pause(1);close(f6);
%         end
%     end
% end
% end
% 
% function plotcorrelationswiringnetwork(allmetricslist)
% %plotcorrelationswiringnetwork plot correlation between wiring cost and
% %network metrics defined in the cell allmetricslist
% %IN: allmetricslist{patients, conditions, frequencies}
% %fixed frequency
% [pats, conds, freqs] = size(allmetricslist);
% freqsinit = 1; % fixed one frequency
% freqsend = size(allmetricslist,3);
% for f=freqsinit:freqsend
%     fixedfreq = allmetricslist(:,:,f);
%     for c=1:conds
%         fixedcond = allmetricslist(:,c,f);
%         for p=1:pats
%             fixedpat = fixedcond{p};
%             legend = fixedpat{1}; % ' PRE' 'Pat' 'fq' ' channels'
%             wc_pli = fixedpat{3};
%             wc_r = fixedpat{4};
%             wc_change = fixedpat{5};
%             wc_phys = fixedpat{6};
%             metrics = fixedpat{2}; %list f metric function calculategraphmatrixmetrics
%             wc_power = fixedpat{end};
%             %allmetrics = {degree_v, gtommatrix,matchingmatrix,density_coeff,clusteringcoeff_v,transcoeff,components_v,componentsizes_v,communityLovaffi_v,assort_coeff, richclub_v, coreperip_v, Dmatrix, charpathlength_coeff, BCvector_normalized, pagerank_vector}
%             clusteringv = metrics{5};
%             clustering = mean(clusteringv);
%             clusteringmatrix(f,c,p) = clustering;
%             degreev = metrics{1};
%             degree = mean(degreev);
%             degreematrix(f,c,p) = degree;
%             pagerankv = metrics{end};
%             pagerank = mean(pagerankv);
%             pagerankmatrix(f,c,p) = pagerank;
%             charpathv = metrics{end-2};
%             charpath = mean(charpathv);
%             charpathmatrix(f,c,p) = charpath;
%             %         richclubv = metrics{end-5};
%             %         richclub = mean(richclubv);
%             %         richclubmatrix(f,c,p) = richclub;
%             assort_coeffv = metrics{end-6};
%             assort_coeff = mean(assort_coeffv);
%             assort_coeffvmatrix(f,c,p) = assort_coeff;
%             %wiringcost
%             wc_pli_mean = mean2(wc_pli);
%             wc_r_mean = mean2(wc_r);
%             wc_change_mean = mean2(wc_change);
%             wc_phys_mean = mean2(wc_phys);
%             wc_power_mean = mean2(wc_power);
%             wc_pli_x(f,c,p) = wc_pli_mean;
%             wc_r_x(f,c,p) = wc_r_mean;
%             wc_change_x(f,c,p) = wc_change_mean;
%             wc_phys_x(f,c,p) = wc_phys_mean;
%             wc_power_x(f,c,p)= wc_power_mean;
%         end
%     end
% end
% %plot the correlations between wc_* and degreematrix
% %cond 1
% figcorr= figure;
% cond =1; freq = 1;
% for freq=1:freqsend
%     subplot(freqsend,1,freq)
%     y= wc_r_x(freq,cond,:);
%     x = clusteringmatrix(freq,cond,:);
%     s1 = scatter(x,y,'MarkerFaceColor','red');
%     hold on;
%     cond = 2;
%     y2= wc_r_x(freq,cond,:);
%     x2 = clusteringmatrix(freq,cond,:);
%     %lsline
%     s2 =scatter(x2,y2,'MarkerFaceColor','blue');
%     xlabel('clusteringmatrix'), ylabel('Wiring Cost');
%     titls = sprintf('Frequency =%s' , num2str(freq))
%     title(titls);
%     %lsline;
%     xv = []; yv = [];xv2 = []; yv2 = [];
%     for pi=1:pats
%         xv = [xv x(pi)];
%         xv2 = [xv2 x2(pi)];
%         yv = [yv y(pi)];
%         yv2 = [yv2 y2(pi)];        
%     end
%     corr_closed = corr2(xv,yv);
%     corr_open = corr2(xv2,yv2);
%     corr1t = sprintf('EC: Wiring, clustering corr=%.2f', corr_closed);
%     corr2t = sprintf('EO: Wiring, clustering corr=%.2f', corr_open);
%     text(mean(xv),mean(yv),corr1t);
%     text(mean(xv2),mean(yv2),corr2t);
%     %legend('C', 'O');
% end
% for freq=1:freqsend
%     subplot(freqsend,1,freq)
%     y= wc_r_x(freq,cond,:);
%     x = charpathmatrix(freq,cond,:);
%     s1 = scatter(x,y,'MarkerFaceColor','red');
%     hold on;
%     cond = 2;
%     y2= wc_r_x(freq,cond,:);
%     x2 = charpathmatrix(freq,cond,:);
%     lsline;
%     s2 =scatter(x2,y2,'MarkerFaceColor','blue');
%     xlabel('charpathmatrix'), ylabel('Wiring Cost');%legend('C', 'O');
%     ts = sprintf('Frequency =%s', num2str(freq));
%     title(ts);
%     lsline;
%     xv = []; yv = [];xv2 = []; yv2 = [];
%     for pi=1:pats
%         xv = [xv x(pi)];
%         xv2 = [xv2 x2(pi)];
%         yv = [yv y(pi)];
%         yv2 = [yv2 y2(pi)];        
%     end
%     corr_closed = corr2(xv,yv);
%     corr_open = corr2(xv2,yv2);
%     corr1t = sprintf('EC: Wiring, charpathmatrix corr=%.2f', corr_closed);
%     corr2t = sprintf('EO: Wiring, charpathmatrix corr=%.2f', corr_open);
%     text(mean(x),mean(y),corr1t);
%     text(mean(x2),mean(y2),corr2t);
% end
% % Fit = polyfit(x,y,1);
% % plot(polyval(Fit,x));
% % Fit = polyfit(x2,y2,1);
% % plot(polyval(Fit,x2));
% 
% 
% end