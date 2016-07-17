function [] = calculatewiringcostmatrices()
%calculatewiringcostmatrices calculate the wiring cost matric and represents the
% results
%conditionsL = {'EC_PRE', 'EC_POST'};
%patientsL  = { 'TWH030', 'TWH031','TWH033','TWH034'};
%[srate, min_freq, max_freq, num_frex, time, n_wavelet, half_of_wavelet_size, frex, s, wavelet_cycles] =   initialize_wavelet()
%freqs2use  = logspace(log10(min_freq),log10(max_freq),8);
globalFsDir = loadglobalFsDir();
matfilename = fullfile(globalFsDir, 'phaseconn_matrices.mat');
fh = load(matfilename);
patientslist = fh.phaseconn_matrix.patientsl;
conditionslist = fh.phaseconn_matrix.conditionsl;
frequencylist = fh.phaseconn_matrix.freqsl;
wiring_matrices = struct ;
label = 'BiTemp';
%the distMatrix is casted with the vector of electrodes indices [2 3 4 ...66] 

for i=1:size(patientslist,2)
    if strcmp(patientslist{i}, 'TWH031') == 1
        %31 has surf folder missing
        continue
    end
    channelclassindexes = getindexesfromlabel(patientslist{i}, label);
    channelclassindexes = channelclassindexes -1;
    distMatrix = calculateEuclideandistancematrix(patientslist{i});
    distMatrixcell{i} = distMatrix;
    for j=1:size(conditionslist,2)
        ispc_matrix_freqs = fh.phaseconn_matrix.ispc_matrix{i,j};% [ 36x36x8  double]    [ 36x36x8  double]
        pli_matrix_freqs = fh.phaseconn_matrix.pli_matrix{i,j};
        %icoh_matrix_freqs = fh.phaseconn_matrix.icoh_matrix{i,j};
        for k=1:size(frequencylist,2)
            %matrix [channelsxchannels] per each frequency
            ispc_matrix = ispc_matrix_freqs(:,:,k);
            % cast for Label type electrodes
            ispc_matrix= ispc_matrix(channelclassindexes,channelclassindexes);
            pli_matrix = pli_matrix_freqs(:,:,k);
            pli_matrix= pli_matrix(channelclassindexes,channelclassindexes);
            %icoh_matrix = icoh_matrix_freqs(:,:,k);
            wiring_ispc{i,j,k} = distMatrix.*ispc_matrix;
            wiring_pli{i,j,k} = distMatrix.*pli_matrix;
            %wiring_icoh{i,j,k} = ditsMatrix.*icoh_matrix;
        end
    end
end
wiring_matrices.wiring_ispc = wiring_ispc;
wiring_matrices.wiring_pli = wiring_pli;
wiring_matrices.patientslist = patientslist;
wiring_matrices.conditionslist = conditionslist;
wiring_matrices.frequencylist = frequencylist;
wiring_matrices.distMatrixcell = distMatrixcell;
%savethe wiring_matrices structure
matfilename = 'wiringcost_matrices.mat';
matfilename = fullfile(globalFsDir, matfilename);
disp(wiring_matrices)
save(matfilename,'wiring_matrices');
%plot the wiring cost matrices
plotwiringcostmatrices(wiring_matrices);

end

function [] = plotwiringcostmatrices(wiring_matrices)
% chose patient cond and frequency to pot
patl = wiring_matrices.patientslist;
condl= wiring_matrices.conditionslist;
freql = wiring_matrices.frequencylist;
i = 6;
j = 1;
k = 1;
p = patl{i};
c = condl{j};
f = freql(k); 
wiringcost_ispc = wiring_matrices.wiring_ispc{i,j,k};
wiringcost_pli = wiring_matrices.wiring_pli{i,j,k};
phys_distance = wiring_matrices.distMatrixcell{i};
figure
imagesc(phys_distance);
msgtitle = sprintf('Euclidean Distance BiTemp electrodesfor Patient = %s',p);
title(msgtitle);
figure
imagesc(wiringcost_ispc);
msgtitle = sprintf('Wiring Cost= P.*R for Patient = %s ',p);
title(msgtitle);
figure
imagesc(wiringcost_pli);
msgtitle = sprintf('Wiring Cost= P.*PLI for Patient = %s ',p);
title(msgtitle);
colorbar;
figure
%wiringchange due to connectivity;
wiring_change = phys_distance - wiringcost_ispc
imagesc(wiring_change);
msgtitle = sprintf('Change in Wiring Cost due to phase connectivity (R) for Patient = %s ',p);
title(msgtitle);
colorbar;
fprintf('Patient = %s, Cond=%s, Freq=%s. Mean Wiring cost D=%.2f, D.*R=%.2f, D.*PLI=%.2f \n',p,c,num2str(f),sum(phys_distance)/sum(phys_distance~=0), sum(wiringcost_ispc)/sum(wiringcost_ispc~=0),sum(wiringcost_pli)/sum(wiringcost_pli~=0))
%R=16.84, PLI=14.03
%plot the network 
corrMatrix = wiringcost_ispc;
[threshold, nstds, corrMatrix] = calculatethresholdmatrix(corrMatrix);
[myfullname, EEG, channel_labels, patdate, patsession] =  initialize_EEG_variables(p,c)
strNames = channel_labels(2:end); 
legendofmatrices = {c, p, f, channel_labels};
figure;
myColorMap = lines(length(corrMatrix));
circularGraph(corrMatrix,'Colormap',myColorMap,'Label',strNames);
%title(['Condition: ',currcond,', Frequency:', getgreeksymbolfreq(currfreq), ', Patient:',currentpat, ', Thres. = m+ n*std=', num2str(threshold)]);
title([c,',in', getgreeksymbolfreq(f), ', Patient:',p, ', Thresh. = m+ 1*std']);
drawnow update
%calculate network metrics for corrMatrix
[allmetrics] = calculategraphmatrixmetrics(corrMatrix, legendofmatrices);
allmetrics = {legendofmatrices, allmetrics};
end