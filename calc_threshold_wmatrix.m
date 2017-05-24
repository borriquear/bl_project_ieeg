function [Badj,thresholdv ] =  calc_threshold_wmatrix(WM, binary)
%% calc_threshold_m calculates a vector of binary networks for each possible threshold
%comprised between the minimum and the maximum value of the connectivity
%netwotk
% INPUT: WM weight matrix of functional connectivity, binary = 1
% OUTPUT: {Binary matrix threshold =min...., Binary matrix threshold =max}
%        : thresholdv. Binary if binary == 1 , Weihghted if binary ==0
[nbrows, nbcols] = size(WM);
%  minWM= min(WM(:));
%  maxWM= max(WM(:));
%  stepWM = (maxWM - minWM)/ (nbrows*nbcols);
%  thresholdv = minWM:stepWM:maxWM;
thresholdv = []; nonzero_thresholdv= [];

%thresholdv = sort(reshape(WM,1,nbrows*nbcols));
thresholdv = linspace(min(WM(:)), max(WM(:)), nbrows*nbcols);
%remove the 0 elements (WM is triangular superior so half of the elements are 0)
% only for PLI and R, for Power the connecitivty can be <0
if min(thresholdv) >= 0
    nonzero_thresholdv = find(thresholdv > 0);
    thresholdv = thresholdv(nonzero_thresholdv:end);
else %power Spearman correlation -1, 1
    thresholdv = sort(abs(thresholdv));
    WM = abs(WM);
end
Badj = cell(2,length(thresholdv));

for i=1:length(thresholdv)
    BM = zeros(nbrows, nbcols);
    if binary == 1
        BM(WM < thresholdv(i) & WM > 0) = 1;
    else
        [row, col, index] = find( WM >= thresholdv(i));
        for irow=1:length(row)
            icol=irow;
            BM(row(irow), col(icol)) = WM(row(irow), col(icol));
        end
        %BM(WM >= thresholdv(i));
    end
    %BM(logical(eye(size(BM)))) = 0;
    %BM(reshape(BM,nbrows,nbcols));
    Badj{1,i}= BM;
    %Calculate the network metric for each matrix Badj{i}
    Badj{2,i} = calculatemetrics(BM, binary);
end
disp(Badj)
end
function netwmets = calculatemetrics(Badj, binary)
%% calculatemetrics calculates the netwotk metrix for one binary matrix
% INPUT: Badj binary matrix of functional connectivity
% OUTPUT: Hash vector network metric for the binary matrix

%keySet =   {'clustering', 'density', 'pathlength', 'richclub', 'transitivity','B0'};
%keySet =   {'clustering', 'density', 'pathlength', 'transitivity','B0'};
if binary == 0
    keySet =   {'wiringcost'};
else
    keySet =   {'B0','clustering', 'wiringcost', 'pathlength'}; 
    %wiring cost for binary matrix is the same as the number of edges
end

valueSet = zeros(1,length(keySet));
for i=1:length(keySet)
    valueSet(i) = calculatemetricshash(Badj, keySet{i});
end

netwmets = containers.Map(keySet,valueSet);

end

function netwmets_el = calculatemetricshash(corrmatrix, label)
%% calculatemetricshash calculates the metric value for a binary matrix and a criterion
%INPUT: corrmatrix, label
%OUTPUT: netwmets (realnumber)
%Set the diagonal to 0

if strcmp(label, 'clustering')==1
    netwmets_el = mean(clustering_coef_bu(corrmatrix));
elseif strcmp(label, 'density')==1
    %Density: Density is the fraction of present connections to possible connections
    netwmets_el = density_und(corrmatrix);
elseif strcmp(label, 'pathlength')==1
    Dmatrix = distance_bin(corrmatrix);
    % The average shortest path length is the  characteristic path length of the network.
    charpathlength_coeff= mean2(Dmatrix(~isinf(Dmatrix)));%avoid NaN for disconnected nodes
    netwmets_el = charpathlength_coeff;
    % elseif strcmp(label, 'richclub')==1
    %     richclub_v = rich_club_bu(corrmatrix);
    %     netwmets = mean(richclub_v);
    % elseif strcmp(label, 'transitivity')==1
    %     netwmets_el = transitivity_bu(corrmatrix);
elseif strcmp(label, 'B0')==1
    % Count number of voxels for B0 and other Betti number use Topological
    % analysis Toolbox
    %s=graph_spectrum(corrmatrix);
    %nc=numel(find(s<10^(-5)));
    %normalize nb of components
    %netwmets_el = nc/size(corrmatrix,1);
    %     [components_v, componentsizes_v] = get_components(corrmatrix);
    %     netwmets_el = mean(componentsizes_v);
    %[S, c] = graphconncomp(sparse(triu(corrmatrix)), 'Directed', 'False')
    nc = numedges(corrmatrix);
    netwmets_el = nc/(size(corrmatrix,1)*size(corrmatrix,1));
    
elseif  strcmp(label, 'wiringcost')==1
    netwmets_el = sum(corrmatrix(:))/ (size(corrmatrix,1)* size(corrmatrix,2));
    %disp(netwmets_el)
else
    fprintf('ERROR label %s do not found!!', label );
end
if isnan(netwmets_el) == 1
    fprintf('WARNING: isnan value in %s , forcing to 0 \n',label );
    netwmets_el(isnan(netwmets_el))=0;
end
end
