function [counter_pats,counter_channels] = ttest2andreport(alpha,c1,c2,f,pat,reportfile, type_of_m)
%ttest between two samples and write report in file txt
counter_pats=0;counter_channels=0;
[h,p] = ttest2(c1,c2);
[cha_idx] = find(p<alpha);
[cha_val] = p(cha_idx);
if any(h)
    counter_pats= counter_pats + 1;
    counter_channels = counter_channels + length(cha_idx);
    fprintf(type_of_m, '** %s Reject Null hypothesis for frequency %s subject:%s\n', type_of_m, num2str(f), num2str(pat));
    fprintf(reportfile,'** %s Reject Null hypothesis for frequency %s subject:%s\n', type_of_m, num2str(f), num2str(pat));
    fprintf(reportfile,'\t ** %s Channels idx for %s subject:%s = %s\n', type_of_m, num2str(f), num2str(pat), num2str(cha_idx));
    fprintf(reportfile,'\t ** %s p-values for %s subject:%s = %s\n', type_of_m, num2str(f), num2str(pat), num2str(cha_val));
    %disp(p)
else
    %fprintf(' Fail to reject Null hypothesis for frequency %s subject:%s\n', num2str(freq_list(f)), num2str(pat))

end
end