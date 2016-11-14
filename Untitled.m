%>>print -f3 -djpeg 'D:\BIAL PROJECT\patients\figure_results\ECEO-Allroisfigi1'
%>>savefig('D:\BIAL PROJECT\patients\figure_results\ECEO-Allroisfigi1')
folderderst = 'D:\BIALPROJECT\patients\figure_results\';
nbpats = 11;
for i=1:nbpats
    fileorig = sprintf('pat-%d-WC-coherence-EC-EO', i);
    filedest =fullfile(folderderst,fileorig );
    figf = sprintf('-f%d',i);
    fptimtf('Savingfile %s\n',filedest)
    eval(['print ',figf, ' -djpeg ' filedest ])
    savefig(filedest)
end
