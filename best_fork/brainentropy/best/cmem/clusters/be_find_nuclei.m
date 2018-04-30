function Seed = be_find_nuclei(APM,Vertices,VertConn,AlreadyClustered,Scale,OPTIONS)
%
% BE_FIND_NUCLEI identifies the possible seeds as local 
% max of the MSP associated with one component
%
% INPUTS: 
%        -APM                   : Contains the MSP score associated to each principal component.
%        -Vertices              : List of vertices in the cortical mesh.
%        -VertConn              : Neighborhood matrix. 
%        -AlreadyClustered      : Vertices that are already clustered.
%        -Scale                 : (Set to default 4) Nb. of neighbors used to clusterize cortical surface.
%        -OPTIONS				: Structure of parameters (described in be_cmem_solver.m)
% OUTPUTS: 
%        -Seed                  : the vector containing the indices of the centers of gravity of each cluster
%
%% ==============================================
% Copyright (C) 2011 - MultiFunkIm & LATIS Team
%
%  Authors: MultiFunkIm team, LATIS team, 2011
%
%% ==============================================
% License
%
% BEst is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BEst is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BEst. If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

% Sort APM coefficients
[tmp,I] = sort(-APM);
clear tmp

% Neighbourhood matrix at a specific scale
% Adds connectivity on the diagonal (there should be zeros on the diagonal of nm
%initially)
VertConn=VertConn+speye(size(VertConn));


% Setting the neighbor matrix to the appropriate degree
neighborhood = VertConn^Scale;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1- Pseudo-clusters aggregation, from the sorted MSP values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if OPTIONS.optional.verbose
	fprintf('%s, Pseudo-clusters aggregation from sorted MSP scores ...', OPTIONS.mandatory.pipeline);
end
stop =0;

ToCluster = I;
ToCluster = (~ismember(ToCluster,AlreadyClustered)).*ToCluster;


clusters = cell(0,0);
k = 0;

while stop == 0
    
    SelSeed = ToCluster(min(find(ToCluster ~= 0)));
    
    k = k+1;
    
    %find corresponding neighbourhood at a specific Scale
    ListNeighbours = find(neighborhood(SelSeed,:) ~=0);
    
    clusters{k} = setdiff(ListNeighbours,AlreadyClustered);
    
    
    AlreadyClustered = [AlreadyClustered,clusters{k}];
    ToCluster = (~ismember(ToCluster,AlreadyClustered)).*ToCluster;
    
    if isempty(find(ToCluster ~= 0)) == 1
        stop = 1;
    end
    
end
clear AlreadyClustered  ToCluster SelSeed ;

if OPTIONS.optional.verbose
	fprintf(' done\n', OPTIONS.mandatory.pipeline);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2- Search for seeds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if OPTIONS.optional.verbose
	fprintf('%s, search for seeds candidates ...', OPTIONS.mandatory.pipeline);
end
Seed = [];

for k = 1:length(clusters)
    
    weight = diag(APM(clusters{k}));
    G = sum(weight * Vertices(clusters{k},:),1)./sum(APM(clusters{k}));
    
    if (size(clusters{k},2) ==1)
        Seed(k) = clusters{k};
    else
        distG = sum((Vertices(clusters{k},:) - repmat(G,size(Vertices(clusters{k},:),1),1)).^2,2);
        Seed(k) = clusters{k}(min(find(distG == min(distG))));
    end
end

if OPTIONS.optional.verbose
	fprintf(' done\n', OPTIONS.mandatory.pipeline);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3- Cleaning of the seeds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if OPTIONS.optional.verbose
	fprintf('%s, cleaning of the seeds ...', OPTIONS.mandatory.pipeline);
end
stop = 0;
k = 0;
while stop == 0
    
    k = k + 1;

    
    ListNeighbours = find(neighborhood(Seed(k),:) ~=0);
    
    % remove the seed point from the list for neighbours
    ListNeighbours = setdiff(ListNeighbours,Seed(k));
    
    voi = intersect(ListNeighbours,Seed);
    if isempty(voi) == 0 
        i = find(ismember(Seed,voi)==1);
        Seed(i) = [];
        clear voi i;
    end
    if k >= length(Seed)
        stop = 1;
    end
    
end

if OPTIONS.optional.verbose
	fprintf(' done\n', OPTIONS.mandatory.pipeline);
end

return 





