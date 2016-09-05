function [] = scatterPF(wiring_matrices)
%scatterPF scatter physical and functional distance in 1 axis, physical
%distance and functional distance
%plotting the scatter of physical and functional distance

%pairlist = {'EC_PRE', 'EO_PRE', 'HYP'};
phys_data = {};
functional_data = {};
globalFsDir = loadglobalFsDir();
%matfilename = fullfile(globalFsDir, 'phaseconn_matrices.mat');
matfilename = fullfile(globalFsDir, 'phaseconn_matrices_tw.mat');
fh = load(matfilename);
% matfilename = fullfile(globalFsDir, 'powerconn_matrices_tw.mat');
% fh2 =  load(matfilename);
 
nbpats = size(fh.phaseconn_matrix.patientsl,2);
pairlist = fh.phaseconn_matrix.conditionsl;
nbconds = length(fh.phaseconn_matrix.conditionsl);
nbfreqs = length(fh.phaseconn_matrix.freqsl);
initpats = 1;endpats = nbpats;
initconds = 1; endconds = nbconds; %EC EO HYP
initfreqs = 3; endfreqs = nbfreqs; %3.7...50
for i=initpats:endpats
    xlist= [];
    x = [];
    x = wiring_matrices.distMatrixcell{i}; % physical distance for each patient
    for iaux=1:size(x,1)
        for jaux=iaux+1:size(x,2)
            xlist = [xlist x(iaux,jaux)];
        end
    end
    phys_data{i} = xlist;
    for j=initconds:endconds
        ispc_matrix_freqs = fh.phaseconn_matrix.ispc_matrix{i,j};% [ 36x36x8  double]    [ 36x36x8  double]
        pli_matrix_freqs = fh.phaseconn_matrix.pli_matrix{i,j};
        %power_matrix_freqs = fh2.powerconn_matrix.power_matrix{i,j};
        %icoh_matrix_freqs = fh.phaseconn_matrix.icoh_matrix{i,j};
        for k=initfreqs:endfreqs 
            y_r=[]; y_pli = [];
            ispc_matrix = ispc_matrix_freqs(:,:,k);
            pli_matrix = pli_matrix_freqs(:,:,k);
            %functional distance for patient,cond,freq
            for iaux=1:size(x,1)
                for jaux=iaux+1:size(x,2)
                    y_r = [y_r ispc_matrix(iaux, jaux)];
                    y_pli = [y_pli pli_matrix(iaux, jaux)]; 
                end
            end
            functional_data{i,j,k,1} = y_r; functional_data{i,j,k,2} = y_pli;
        end
    end
end
plotscatterperpatient(wiring_matrices, phys_data, functional_data, initpats, endpats, initconds, endconds, initfreqs, endfreqs);
end

function plotscatterperpatient(wiring_matrices, phys_data, functional_data, initpats, endpats, initconds, endconds, initfreqs, endfreqs)
% plotscatterperpatient for all patients in onecond and freq

for c=initconds:endconds
    for f=initfreqs:endfreqs
        phys_list = []; %vector with physical distance between pair of electrodesfor all patients
        funct_list_1 = []; funct_list_2 = [];
        for p=initpats:endpats
            phys_pat = phys_data{p};
            phys_list = [phys_list phys_pat];
            phys_pat= [];
            funct_pat_1 = functional_data{p,c,f,1};
            funct_pat_2 = functional_data{p,c,f,2};
            funct_list_1 = [funct_list_1 funct_pat_1 ];
            funct_list_2 = [funct_list_2 funct_pat_2 ];
            funct_pat_1 = {}; funct_pat_2= {};
        end
        %plot all patients for that f and c
        if length(phys_list) ~=  length(funct_list_1)
            error('Error. \n Number of pairs of electrodes for Physical and Functional must be equal, %d %d',length(phys_list), length(funct_list_1) );
        else
            figure; %R value
            scatter(phys_list, funct_list_1);
            hold on
            [ne,dh] = hist3([phys_list', funct_list_1']);
            contour(dh{1},dh{2},ne)
            [rho, pval] = corr(phys_list', funct_list_1');
            stitle = sprintf('PxR correlation, Cond =%s, Frq=  %.2f, r=%.3f, p-value=%.3f', wiring_matrices.conditionslist{c},wiring_matrices.frequencylist(f), rho, pval);
            title(stitle, 'Interpreter', 'none');
            xlabel('Physical distance electrode pairs'), ylabel('Functional correlation')
            figure; %PLI value
            scatter(phys_list, funct_list_2);
            hold on
            [ne,dh] = hist3([phys_list', funct_list_2']);
            contour(dh{1},dh{2},ne)
            [rho, pval] = corr(phys_list', funct_list_2');
            stitle = sprintf('PxPLI correlation, Cond =%s, Frq=  %.2f, r=%.3f, p-value=%.3f', wiring_matrices.conditionslist{c},wiring_matrices.frequencylist(f), rho, pval);
            title(stitle, 'Interpreter', 'none');
            xlabel('Physical distance electrode pairs'), ylabel('Functional correlation')
            ylim([0 1])
        end
    end
end
end
