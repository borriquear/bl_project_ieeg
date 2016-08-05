function [threshold, nstds, corrMatrix] = calculatethresholdmatrix(corrMatrix)
%calculatethresholdmatrix(corrMatrix) returns the threshold and the
%threshold matrix given a correlation matrix
%IN: corrMatrix [0,1] if [-1,1] the function calculates abs(corrMatrix)
%OUT: threshold, new threshold matrix

%meanmatrix = mean2(abs(corrMatrix)); stdmatrix = std2(abs(corrMatrix));
meanmatrix = mean(abs(corrMatrix(:))); stdmatrix = std(abs(corrMatrix(:)));
nstds = 2; %number of standard deviations
threshold = meanmatrix + nstds*stdmatrix;
%threshold = 0;
fprintf('Corr matrix mean=%2.4f +(n)%d*(std)%2.4f = %2.4f\n', meanmatrix, nstds,stdmatrix, threshold);
%if not symmetric make (upper or lower daigonal)
% if issym(corrMatrix) < 1
%     for i=1:size(corrMatrix,1)
%         for j =1:i-1
%             corrMatrix(i,j)= corrMatrix(j,i);
%         end   
%     end
% end
corrMatrix(corrMatrix <= threshold) = 0;
corrMatrix(corrMatrix >  threshold) = 1;

end