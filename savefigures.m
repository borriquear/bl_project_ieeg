%>>print -f3 -djpeg 'D:\BIAL PROJECT\patients\figure_results\ECEO-Allroisfigi1'
%>>savefig('D:\BIAL PROJECT\patients\figure_results\ECEO-Allroisfigi1')
folderderst = 'D:\BIALPROJECT\patients\figure_results\';
nbpats = 11;
for i=1:nbpats
    %fileorig = sprintf('pat-%d-WC-coherence-EC-EO', i);
    fileorig = sprintf('pat-%d-Binary-WC-coherence-EC-EO', i);
    filedest =fullfile(folderderst,fileorig );
    figf = sprintf('-f%d',i);
    fprintf('Saving file %s\n',filedest)
    eval(['print ',figf, ' -djpeg ' filedest ])
    savefig(filedest)
end
%fileorig = sprintf('MeanWCDiff-coherence-patient-EC-EO', i);
fileorig = sprintf('MeanWCDiff-coherence-Binary-patient-EC-EO', i);
filedest =fullfile(folderderst,fileorig );
figf = sprintf('-f%d',12);
fprintf('Saving file %s\n',filedest)
eval(['print ',figf, ' -djpeg ' filedest ])
savefig(filedest)

%fileorig = sprintf('MeanWCDiff-coherence-frequency-EC-EO', i);
fileorig = sprintf('MeanWCDiff-coherence-Binary-frequency-EC-EO', i);
filedest =fullfile(folderderst,fileorig );
figf = sprintf('-f%d',13);
fprintf('Saving file %s\n',filedest)
eval(['print ',figf, ' -djpeg ' filedest ])
savefig(filedest)
%%
folderderst = 'D:\BIALPROJECT\patients\figure_results\';
%fileorig = sprintf('pat2-networkmetrics-coherence-EC-EO', i);
fileorig = sprintf('pat1-networplot-EC-fq1-max');
filedest =fullfile(folderderst,fileorig );
figf = sprintf('-f%d',1);
fprintf('Saving file %s\n',filedest)
eval(['print ',figf, ' -djpeg ' filedest ])
savefig(filedest)


