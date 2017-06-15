function [] = performttestandwilcox_connectivity(wiring_matrix_structure, type_of_matrix)
  
alpha = 0.05;
globalFsDir = loadglobalFsDir();
reportfile = strcat(type_of_matrix,'_report-ttest-adj.txt'); 
txt_reportfile = fullfile(globalFsDir, reportfile);
reportfile = fopen(txt_reportfile,'w');   

if strncmpi(type_of_matrix,'pow',3) == 1
    wiring_matrix = wiring_matrix_structure.powerconn_matrix.power_matrix;
    freq_list = wiring_matrix_structure.powerconn_matrix.freqsl;
    pats_list = wiring_matrix_structure.powerconn_matrix.patientsl;
    fprintf(reportfile,' Calculating ttest for Connectivity-based Power \n');
elseif strncmpi(type_of_matrix,'pli',3) == 1
    fprintf(reportfile,' Calculating ttest for Phase-based Power \n');
    wiring_matrix = wiring_matrix_structure.phaseconn_matrix.pli_matrix;
    freq_list = wiring_matrix_structure.phaseconn_matrix.freqsl;
    pats_list = wiring_matrix_structure.phaseconn_matrix.patientsl;
elseif strncmpi(type_of_matrix,'ispc',4) == 1
    fprintf(reportfile,' Calculating ttest for Phase-based Power \n');
    wiring_matrix = wiring_matrix_structure.phaseconn_matrix.ispc_matrix;
    freq_list = wiring_matrix_structure.phaseconn_matrix.freqsl;
    pats_list = wiring_matrix_structure.phaseconn_matrix.patientsl;
end

for f=1:length(freq_list)
    fprintf(reportfile,'\n\n Calculating ttest for Freq:%s\n',num2str(freq_list(f)));
    counter_power=0;counter_power_channels=0;
    for pat=1:length(pats_list)
        if strncmpi(type_of_matrix,'pli',3) + strncmpi(type_of_matrix,'ispc',4) >0 
            c1 = wiring_matrix(pat,1);
            c2 = wiring_matrix(pat,2);
            c1 = c1{1};
            c2 = c2{1};
            c1 = c1(:,:,f);
            c2 = c2(:,:,f);
        else
            c1= wiring_matrix{pat,1,f};
            c2= wiring_matrix{pat,2,f};
        end
        fprintf('connectivity-based ttest for patient:%s\n', num2str(pat));
        [counter_pats, counter_channels] = ttest2andreport(alpha,c1,c2,freq_list(f),pat,reportfile,'conn-based-Power');
        counter_power= counter_power + counter_pats;
        counter_power_channels = counter_power_channels + counter_channels;
    end
    fprintf(reportfile,'\n Total patients (/11) with some chaanel p<alpha at Freq:%s is:  %s\n', num2str(freq_list(f)),num2str(counter_power));
    fprintf(reportfile,'\n Total channels p<alpha at Freq:%s is: = %s\n', num2str(freq_list(f)), num2str(counter_power_channels));
    fprintf(reportfile,'---------------------------------------------------------------------------------------');
end
fprintf('\n Report has been generated at :%s\n',txt_reportfile);
fclose(reportfile);
end
