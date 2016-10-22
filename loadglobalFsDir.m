function globalFsDir = loadglobalFsDir()
%load global variable globalFsDir
if ~exist('globalFsDir','var')
    %fprintf('globalFsDir not found, loading it...\n');
    eval('global globalFsDir');
    myp = 'D:\BIALPROJECT\patients\';
    eval(['globalFsDir=' 'myp;']);
end
end