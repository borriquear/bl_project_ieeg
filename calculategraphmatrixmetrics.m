function [allmetrics] = calculategraphmatrixmetrics(corrmatrix, legend)
%calculategraphmatrixmetrics calculate network metrics from corrmatrix
%which is created in graphtheoryanalysis for a given threshold
% sizeofmats = size(corrmatrix);
% sizeoflegend = size(legend);
%Degree and Similarity
fprintf('Degree and Similarity metrics')
degree_v = degrees_und(corrmatrix);
numSteps = 1;
%gtom,for each pair of nodes, the fraction of neighbors the two
%nodes share in common, where "neighbors" are one step removed.
%As 'numSteps' is increased, neighbors that are furter out are considered.
%Elements of 'gt' are bounded between 0 and 1.
gtommatrix = gtom(corrmatrix,numSteps);
%Matching index: The matching index computes for any two nodes u and v, 
%the amount of overlap in the connection patterns of u and v. 
%Self-connections and u-v connections are ignored. 
%The matching index is a symmetric quantity, similar to a correlation or a dot product
matchingmatrix = matching_ind_und(corrmatrix);
%Density and Scaling
fprintf('Density and Scaling metrics')
%Density: Density is the fraction of present connections to possible connections
density_coeff = density_und(corrmatrix);

%Rentian scaling
%XYZ  Vector of node placement coordinates
% must be Mx3 matrix, where M is the number of nodes of corrmat
%n Number of partitions to compute. Each partition is a data
%  point. You want a large enough number to adequately estimate
%  the Rent's exponent.
%Efficient physical embedding of topologically complex information processing
%networks in brains and computer circuits
dorentialscaling = 0;
if dorentialscaling ==1
    npartitions = 5000;%n=5000 was used in Bassett et al. 2010 in PLoS CB.
    [Nrs, Ers] = rentian_scaling(corrmatrix,XYZ,npartitions);
    loglog(Ers,Nrs,'*');
end
%Clustering and Community Structure
fprintf('Clustering and Community Structure metrics');
clusteringcoeff_v = clustering_coef_bu(corrmatrix);%returns scalar
transcoeff = transitivity_bu(corrmatrix);%returns scalar
[components_v,componentsizes_v] = get_components(corrmatrix);
% Iterative community finetuning.
% corrmatrix is the input connection matrix.
nnodes  = size(corrmatrix,1);             % number of nodes
M  = 1:nnodes;                   % initial community affiliations
Q0 = -1; Q1 = 0;            % initialize modularity values
while Q1-Q0>1e-5;           % while modularity increases
    Q0 = Q1;                % perform community detection
    [communityLovaffi_v, Q1] = community_louvain(corrmatrix, [], M);
end
%Assortativity
%Assortativity: The assortativity coefficient is a correlation coefficient
%between the degrees of all nodes on two opposite ends of a link. 
%A positive assortativity coefficient indicates that nodes tend to 
%link to other nodes with the same or similar degree.
%Assortativity and Core structure
fprintf('Assortativity and Core structure metrics')
assort_coeff = assortativity_bin(corrmatrix,0);
richclub_v = rich_club_bu(corrmatrix);
klevel = 6;
% rich club coefficient at level k is the fraction of edges that 
%connect nodes of degree k or higher out of the maximum 
%number of edges that such nodes might share.
[R,Nk,Ek] = rich_club_bu(corrmatrix,klevel);
%The optimal core/periphery subdivision is a partition of the network into two non-overlapping groups of nodes
%C, binary vector of optimal core structure
% C = 1 represents nodes in the core
% C = 0 represents nodes in the periphery
try
coreperip_v  = core_periphery_dir(corrmatrix) 
catch ME
    warning('Problem using function.  core_periphery_dir(corrmatrix)');
    coreperip_v = 0;
end
    
%Paths and Distances
fprintf('Paths and Distances metrics')
%Dmatrix = lengths of shortest paths between all
%pairs of nodes. 
%An entry (u,v) represents the length of shortest path 
%from node u to node v.
Dmatrix = distance_bin(corrmatrix);
% The average shortest path length is the  characteristic path length of the network.
charpathlength_coeff= mean2(Dmatrix(~isinf(Dmatrix)));%avoid NaN for disconnected nodes
%Centrality 
%Betweenness centrality: Node betweenness centrality is the fraction of 
%all shortest paths in the network that contain a given node. 
%Nodes with high values of betweenness centrality participate in a 
%large number of shortest paths
BCvector = betweenness_bin(corrmatrix);
%Note: Betweenness centrality may be normalised to the range [0,1] as
%BC/[(N-1)(N-2)], where N is the number of nodes in the network.
BCvector_normalized = BCvector/((nnodes -1)*(nnodes-2));
% The PageRank centrality is a variant of eigenvector centrality. This
%function computes the PageRank centrality of each vertex in a graph.
%PageRank is defined as the stationary distribution achieved
%by instantiating a Markov chain on a graph. The PageRank centrality of
%a given vertex, then, is proportional to the number of steps (or amount
%of time) spent at that vertex as a result of such a process.
%d = prob the damping factor specifies the
%     fraction of the time that a random walker will transition to one of its
%     current state's neighbors. The remaining fraction of the time the
%     walker is restarted at a random vertex. A common value for the damping
% factor is dampingf = 0.85.
dampingf = 0.85;%falff = 0; %si pongo damp = 1 NaN!
pagerank_vector = pagerank_centrality(corrmatrix, dampingf);
allmetrics = {degree_v, gtommatrix,matchingmatrix,density_coeff,clusteringcoeff_v,transcoeff,components_v,componentsizes_v,communityLovaffi_v,assort_coeff,richclub_v,coreperip_v,Dmatrix,charpathlength_coeff,BCvector_normalized,pagerank_vector}
disp('degree_v is'); disp(degree_v);
disp('gtommatrix is'); disp(gtommatrix);
disp('matchingmatrix is'); disp(matchingmatrix);
disp('density_coeff is'); disp(density_coeff);
disp('clusteringcoeff_v is'); disp(clusteringcoeff_v);
disp('transcoeff is'); disp(transcoeff);
disp('components_v is'); disp(components_v);
disp('componentsizes_v is'); disp(componentsizes_v);
disp('communityLovaffi_v is'); disp(communityLovaffi_v);
disp('assort_coeff is'); disp(assort_coeff);
disp('richclub_v is'); disp(richclub_v);
disp('coreperip_v is'); disp(coreperip_v);
disp('charpathlength_coeff is'); disp(charpathlength_coeff);
disp('BCvector_normalized is'); disp(BCvector_normalized);
disp('pagerank_vector is'); disp(pagerank_vector);
end