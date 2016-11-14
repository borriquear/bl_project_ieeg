function [wiring_matrices] = calculatewiringcostmatrices(scope)
%calculatewiringcostmatrices calculate the wiring cost matrix.
%PAteitns, frequency and conditions list is load from powerconn_matrices
%and phaseconn_matrices
% INPUT: scope = 'local' , 'meso'.
%                 scope == local -> W = P.*F (the wiring cost is between pairs)
%                 scope == meso ->  W = P./H (the wiring cost is taking into
%                 account all the neighbors of the node)
%OUTPUTS wiring_matrices = D*F  or D/H (depending on scope) D is
%                          the Euclidean Distance
%
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
    %from m to cm.
    distMatrix = distMatrix*10^2;
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
            ispc_matrix = ispc_matrix(channelclassindexes,channelclassindexes);
            pli_matrix = pli_matrix_freqs(:,:,k);
            pli_matrix = pli_matrix(channelclassindexes,channelclassindexes);
            icoh_matrix = icoh_matrix_freqs(:,:,k);
            icoh_matrix = icoh_matrix(channelclassindexes,channelclassindexes);
            if strcmp(scope, 'local') == 1
                %Calculate the wiring cost W = distMatrix.*F
                wiring_ispc{i,j,k} = distMatrix.*ispc_matrix;
                wiring_pli{i,j,k} = distMatrix.*pli_matrix;
                wiring_icoh{i,j,k} = distMatrix.*icoh_matrix;
                wiring_power{i,j,k} = distMatrix.*power_matrix;
            elseif strcmp(scope, 'meso') == 1
                % Calculate the wiring cost W = distMatrix./LH (Likelihood matrix of functional connectivity)
                % Transform [ispc, pli, power]_matrix into LH_[ispc, pli, power]_matrix
                % LH = pli(i,j)/ sum( pli(i,:) + pli(:,j))
                LH = wiringcostlikelihood(ispc_matrix);
                distperlh = distMatrix./LH; distperlh(isnan(distperlh))=0;distperlh(isinf(distperlh))=0;
                wiring_ispc{i,j,k} = distperlh;
                LH = wiringcostlikelihood(pli_matrix);
                distperlh = distMatrix./LH; distperlh(isnan(distperlh))=0;distperlh(isinf(distperlh))=0;
                wiring_pli{i,j,k} = distperlh;
                LH = wiringcostlikelihood(power_matrix);
                distperlh = distMatrix./LH; distperlh(isnan(distperlh))=0;distperlh(isinf(distperlh))=0;
                wiring_power{i,j,k} = distperlh;
                LH = wiringcostlikelihood(icoh_matrix);
                distperlh = distMatrix./LH; distperlh(isnan(distperlh))=0;distperlh(isinf(distperlh))=0;
                wiring_icoh{i,j,k} = distperlh;
            end    
        end
        %             fprintf('Mean wiring cost of R matrix  for Patient %s %s =%.2f\n',patientslist{i},conditionslist{j},mean(nonzeros(ispc_matrix(:))));
        %             fprintf('Mean wiring cost of PLI matrix  for Patient %s %s =%.2f\n',patientslist{i},conditionslist{j}, mean(nonzeros(pli_matrix(:))));
        %             fprintf('Mean wiring cost of Power matrix  for Patient %s %s =%.2f\n',patientslist{i},conditionslist{j}, mean(nonzeros(power_matrix(:))));
        %             fprintf('Mean wiring cost of icoh matrix  for Patient %s %s =%.2f\n',patientslist{i},conditionslist{j}, mean(nonzeros(icoh_matrix(:))));
    end
end
wiring_matrices.wiring_ispc = wiring_ispc;
wiring_matrices.wiring_pli = wiring_pli;
wiring_matrices.wiring_icoh = wiring_icoh;
wiring_matrices.wiring_power = wiring_power;
wiring_matrices.patientslist = patientslist;
wiring_matrices.conditionslist = conditionslist;
wiring_matrices.frequencylist = frequencylist;
wiring_matrices.distMatrixcell = distMatrixcell;
%savethe wiring_matrices structure
matfilename = sprintf('wiringcost_matrices.mat_%s', scope);
matfilename = fullfile(globalFsDir, matfilename);
disp(wiring_matrices)
save(matfilename,'wiring_matrices');
end

function LH = wiringcostlikelihood(conn_matrix)
%Transform the connectivity matrix with the likelihood matrix that the
%nodes i j are connected.
[nbr, nbc] = size(conn_matrix);
for i=1:nbr
    for j=1:nbc
        LH(i,j)= conn_matrix(i,j)/ (sum(conn_matrix(i,:)) + sum(conn_matrix(:,j)));
        LH(isnan(LH))=0;
    end
end
end