function [distMatrix] = calculateEuclideandistancematrix(patientid)
%calculateEuclideandistancematrix returns the matrix containing the
%euclidean distance matrix for pair of electrodesbelonging to a label e.g. BiTemp, HD, F... 
%IN: patientid (NOTE patient 31 has not loc file with electrodes coordinates)
%OUT: distance matrix
%patientid  = 'TWH042';
global globalFsDir;
globalFsDir =loadglobalFsDir();
label = 'BiTemp';
%the distMatrix is casted with the vector of electrodes indices [2 3 4 ...66] 
[channelclassindexes] = getindexesfromlabel(patientid, label);
%size of the subset electrodes requested
sizeofReqdistmatrix = length(channelclassindexes);
eegcondition = 'HYP';
[myfullname, EEG, channel_labels, eegdate, eegsession] = initialize_EEG_variables(patientid,eegcondition);
%open file with euclidean distance globalFsDir\patientid\[patienid]PostimpLoc.txt

patpath = strcat(globalFsDir,patientid,'\elec_recon');
fileloc = strcat(patientid,'PostimpLoc.txt');
fileloc = fullfile(patpath,fileloc);
fhloc = fopen(fileloc);
%LAT4 ID X Y Z L S
%channelid= field1=field2
formatfile = '%s %d %f %f %f %c %c';
M = textscan(fhloc, formatfile);
fclose(fhloc);
%celldisp(M);
%M{1} = LAT, M{2} = 4; M{3} =X;M{4} = Y, M{5}= Z,M{6}= L|R ,M{7}= S|D, %
labelsid = M{1}; labelsnb = M{2}; Xpos = M{3};Ypos = M{4}; Zpos = M{5}; 
Side = M{6}; Depth = M{7};
for i = 1:length(labelsid)
    channelsid{i} = strcat(labelsid{i},num2str(labelsnb(i)));
end
% generate euclidean distance matrix

distMatrix = zeros(sizeofReqdistmatrix,sizeofReqdistmatrix);
eventi = 2;
for i = 1:sizeofReqdistmatrix
    labelrow = channel_labels{i+1};
    [isin1 irow] = ismember(labelrow, channelsid);
    %irow = channelclassindexes(i);
    for j = i+1:sizeofReqdistmatrix
        labelcol = channel_labels{j+1};
        [isin2 jcol] = ismember(labelcol, channelsid);
        %jcol = channelclassindexes(j);
        %d(i,j) = sqrt((x-x')^2 + (y-y')^2 + (z-z')^2)
        if isin1 + isin2 > 1
            x = Xpos(irow); y =Ypos(irow); z =Zpos(irow);xx = Xpos(jcol); yy =Ypos(jcol); zz =Zpos(jcol);
            distanceij = sqrt((x-xx)^2 + (y-yy)^2 + (z-zz)^2);
            fprintf('Distance between %s and %s is %f\n',labelrow,labelcol,distanceij );
            distMatrix(i,j) = distanceij;
        end
    end
end
%disp(distMatrix);
patpath = strcat(globalFsDir,patientid,'\elec_recon');
fileloc = strcat(patientid,'euclideandistanceMatrix.mat');
savedmatrix = fullfile(patpath, fileloc);
fprintf('Saving euclidena distance matrix at %s\n',savedmatrix);
save(savedmatrix,'distMatrix','channel_labels','channelclassindexes')
fprintf('Patient = %s. The mean distance across all %s channels is = %.2f\n', patientid, label, mean2(distMatrix));
end