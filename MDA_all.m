function [Mfuncthat, distance] = MDA_all(patientid)
%load MPhys for patientid
%load Mfunctvector
[globalFsDir] = loadglobalFsDir();
direlec = fullfile(globalFsDir, patientid, 'elec_recon');
fileeucldis = strcat(patientid, '_euclideandistance.mat');
filedistance = fullfile(direlec, fileeucldis);
load(filedistance); %loaf matrixP
%load vector functmatrices of functional connectivity networks
direlecfunc =  fullfile(globalFsDir, patientid, 'data', 'figures');
%connfile =  'powerconnectivity_freq_40_HYP_TWH028_10202015_s1.mat'
%connfile = fullfile(direlecfunc,connfile);
randmatrix = randn(size(matrixP,2));
functmatrices{1} = matrixP + randmatrix;
%
typeoffunctmatrices = {'TWH038_HYP_power_alpha', 'TWH038_HYP_power_delta', 'TWH038_HYP_phase_alpha'}
for i=1:length(functmatrices)
    Mfunct = functmatrices(i);
    Mfunct = Mfunct{1}
    [Mfuncthat, distance] = MDA(matrixP, Mfunct);
    ditancev{i} = distance;
    Mfuncthatv{i} = Mfuncthat;
    %save Mfuncthat and distance
end
connfile =  strcat('physicotopologialdistance_', patientid,'.mat')
connfile = fullfile(direlecfunc,connfile);
save(connfile, 'typeoffunctmatrices', 'ditancev','Mfuncthatv');

end
function [Mfuncthat, distance] = MDA(Mphys, Mfunct)
%MDA multidimensional analysis, calculates the monotonic regression between
%Mphys, Mfunct 
%IN: Mphys is the euclidean distance if the elctrodes, Mfunctis the
%functional distance that we are taking into account.
% It can be the white matter distance tract between two electrodes
% binarized by the functional connectivity matrix or the product of Mphys
% weighted by the functional connectivity matrix, Mfunct=
% Mphys*1/weightedconnmatrix
%OUT: Mfuncthat regression , distance = D(Mphys,Mfuncthat)

%http://www.mathworks.com/matlabcentral/fileexchange/47196-graph-based-clustering-and-data-visualization-algorithms/content/improve_JP/toolbox_imp_JP/lsqisotonic.m
%
Mfuncthat = lsqisotonic(Mphys,Mfunct);
normaldistance = norm(Mfuncthat-Mphys, 'fro');
normaldistance_f2hat = norm(Mfuncthat-Mfunct, 'fro');
distance =  normaldistance;
distance = normaldistance_f2hat;
%eucliedandistancematrix = L2_distance(Mfunct,Mfuncthat); %is against the physical, i guess
fprintf('The Frobenius norm is:%s\n', num2str(normaldistance));

end

