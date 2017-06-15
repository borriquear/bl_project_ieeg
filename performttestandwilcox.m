function [] = performttestandwilcox(wiring_matrices, type_of_matrix)
%Perform ttest and Man whitney

alpha = 0.05;
globalFsDir = loadglobalFsDir();
freq_list = wiring_matrices.frequencylist;
pats_list = wiring_matrices.patientslist;
reportfile = strcat(type_of_matrix,'_report-ttest-adj.txt'); 
txt_reportfile = fullfile(globalFsDir, reportfile);
reportfile = fopen(txt_reportfile,'w');
fprintf(reportfile,'Calculating ttest and Utest for matrix:%s\n',type_of_matrix);
for f=1:length(freq_list)
    fprintf(reportfile,'\n\n Calculating ttest and Utest for Freq:%s\n',num2str(freq_list(f)));
    %patients with some channel that reject the null hypothesis
    counter_pli=0;counter_pli_channels =0;
    counter_ispc=0;counter_ispc_channels =0;
    counter_power=0;counter_power_channels =0;
    for pat=1:length(pats_list)
        c1= wiring_matrices.wiring_ispc{pat,1,f};
        c2= wiring_matrices.wiring_ispc{pat,2,f};
        fprintf('wiring_matrices.wiring_ispc ttest for patient%s\n', num2str(pat));
        [counter_pats, counter_channels] = ttest2andreport(alpha,c1,c2,freq_list(f),pat,reportfile,'ISPC');
        counter_ispc= counter_ispc + counter_pats;
        counter_ispc_channels = counter_ispc_channels + counter_channels;
  
        c1=wiring_matrices.wiring_pli{pat,1,f};
        c2=wiring_matrices.wiring_pli{pat,2,f};
        fprintf('wiring_matrices.wiring_pli: ttest for patient%s\n', num2str(pat));
        [counter_pats, counter_channels] = ttest2andreport(alpha,c1,c2,freq_list(f),pat,reportfile,'PLI');
        counter_pli= counter_pli + counter_pats;
        counter_pli_channels = counter_pli_channels + counter_channels;
        c1=wiring_matrices.wiring_power{pat,1,f};
        c2=wiring_matrices.wiring_power{pat,2,f};
        fprintf('wiring_matrices.wiring_power ttest for patient%s\n', num2str(pat));
        [counter_pats, counter_channels] = ttest2andreport(alpha,c1,c2,freq_list(f),pat,reportfile,'POWER');
        counter_power= counter_power + counter_pats;
        counter_power_channels = counter_power_channels + counter_channels;
    end
    fprintf(reportfile,'\n Total patients (/11) with some chnanel p<alpha at Freq:%s is: PLI= %s,  ISPC= %s, Power= %s \n', num2str(freq_list(f)), num2str(counter_pli), num2str(counter_ispc), num2str(counter_power));
    fprintf(reportfile,'\n Total channels p<alpha at Freq:%s is: PLI= %s,  ISPC= %s, Power= %s \n', num2str(freq_list(f)), num2str(counter_pli_channels), num2str(counter_ispc_channels), num2str(counter_power_channels));
    fprintf(reportfile,'---------------------------------------------------------------------------------------');
end
fprintf('\n Report has been generated at :%s\n',txt_reportfile);
fclose(reportfile);

return
end

