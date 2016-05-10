function matrixP = calculatemetricdistance(patientid)
% calculatemetricdistance calculate the metric distance for the electrode
% implant of patientid
%IN: patientid  eg TWH038
%OUT:NxNmtrix containing the eucliden distance between any two pair of
%electrodes
[globalFsDir] = loadglobalFsDir();
direlec = fullfile(globalFsDir, patientid, 'elec_recon');
fileelec = strcat(patientid, 'PostimpLocFormat.txt');
fileelec = fullfile(direlec, fileelec);
fileID = fopen(fileelec, 'r');
if fileID < 0
    errormsg = 'ERROR: PostimpLocFormat file with electrode cordinates did not find in:';
    disp(fileelec);
    error(errormsg);
end
k = 1;
locationvector = [];
electrodesL = [];
while ~feof(fileID)
 tline = fgets(fileID);
 tlinespaced =strread(tline,'%s','delimiter',' ');
 locationvector{k} = tlinespaced(1:4);
 electrodesL{k} = tlinespaced{1};
 k = k + 1;
end
disp('Closing distance file')
fclose(fileID);
disp('Calculating the euclidean matrix from locationvector');
matrixP = zeros(k-1,k-1);
for i=1:k-1
    for j=i+1:k-1
        point1 = locationvector(i); 
        p1aux = point1{1};
        point1 = p1aux(2:4);point1 = str2double(point1);
        point2 = locationvector(j); 
        p2aux = point2{1};
        point2 = p2aux(2:4);point2 = str2double(point2);
        distance = euclideandistance(point1, point2);
        matrixP(i,j)= distance; 
    end
    for j=1:i-1
         matrixP(i,j)= matrixP(j,i); 
    end
end
disp('Saving mat file with for M(P) physical euclidean distance matrix in:')
fileeucldis = strcat(patientid, '_euclideandistance.mat');
filedistance = fullfile(direlec, fileeucldis);
save(filedistance, 'electrodesL','matrixP');
end

function distance = euclideandistance(point1, point2)
% euclideandistance returns the euclidean distance between two points in
% the plane
% IN: point1, point2 eg, point1(1) = x coordinate
% OUT: euclidean distance
distance = sqrt( (point1(1) - point2(1))^2 + (point1(2) - point2(2))^2 +(point1(3) - point2(3))^2 );
end