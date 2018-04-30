function [OPTIONS, selected_source, cellstruct] = be_create_clusters(nm, scores, OPTIONS)
%CREATE_CLUSTERS Creates clusters with a neighbor matrix.
%   C = CREATE_CLUSTERS(NEIGHBORS, SCORES,OPTIONS)
%   returns a 1xNS  matrix in C with a cluster number for each dipole.
%   Cluster 0 is the null parcel.
%   NEIGHBORS is the neighbor matrix (0s on the diagonal).  
%   SCORES is a 1xNS matrix of score for each dipole.
%   OPTIONS.neighborhood_order is the is the level of neighborhood for 
%   each cluster, it represents the spatial extent of a cluster.
%   OPTIONS.MSP_scores_threshold is the threshold above wich every dipole 
%   will be selected, based on their score.  
%   If THRESHOLD is set to 0 every dipole will be part of a cluster, i.e.
%   the null parcel will be empty.
%
%   NOTE1:
%   Clusters will be created until one of the conditions is reached.
%   Conditions are : No more dipoles have a SCORE above the THRESHOLD
%
%   NOTE2:
%   A cluster must have at least 3 dipoles to be created.
%
%% ==============================================
% Copyright (C) 2011 - LATIS Team
%
%  Authors: LATIS team, 2011
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

%mem = [OPTIONS.InverseMethod  ', (clustering)'];

% Default options settings
Def_OPTIONS.clustering = struct(...
    'neighborhood_order',   4, ...
    'MSP_scores_threshold', 0 );

% Return empty OPTIONS structure
if (nargin == 0)
    OPTIONS = Def_OPTIONS;
%     clusters = [];
    return
end
% Check field names of passed OPTIONS and fill missing ones with default values
OPTIONS = be_struct_copy_fields(OPTIONS, Def_OPTIONS, {'clustering'});
clear Def_OPTIONS

%Adds connectivity on the diagonal (there should be zeros on the diagonal of nm
%initially)
nm=nm+speye(size(nm));

% Initialization
% Number of sources
nb_sources = numel(scores);     

% Setting the neighbor matrix to the appropriate degree
neighborhood = nm^OPTIONS.clustering.neighborhood_order;

% Sort the scores then find which indices are greater than the threshold
[sorted_scores, indices] = sort(scores,'descend');

thresh_index = nb_sources;
if OPTIONS.clustering.MSP_scores_threshold
    thresh_index = find(sorted_scores >= OPTIONS.clustering.MSP_scores_threshold, 1, 'last');
end


% The tresh_index will be empty if no score are greater than the threshold
if isempty(thresh_index)
    thresh_index = -1;
end

% intialization
ii = 1;
selected_source = zeros(nb_sources,1);
cluster_no = 1;

% Cluster creation
% Clusters will be created until one of the condition is reached.
% Conditions are the threshold and the number parameters
while (ii <= thresh_index)
    % The node with the highest score is selected
    node = indices(ii);
    % Verification that the node is not part of a cluster
    if selected_source(node) == 0
        
        % Getting the neighbors of the node
        % neighbors = unique([find(neighborhood(node,:)) node]);
        neighbors = find(neighborhood(node,:));
        
        % If a node is already in a cluster, it stays in the old cluster.
        neighbors(selected_source(neighbors) ~= 0) = [];
        
        
        % Saving the cluster if its big enough (minimum 3 dipoles)
        if numel(neighbors) >= 3
            selected_source(neighbors) = cluster_no;
%             clusters{1,cluster_no} = neighbors;
            cluster_no = cluster_no + 1;
        end
    end
    ii = ii + 1;
end

% After the first pass is completed, we make sure that all the dipoles over
% the threshold were selected.
% The 'free' dipoles are merged to the lowest nearest cluster.
free_nodes = indices(selected_source(indices(1:thresh_index))==0);
for ii = 1:numel(free_nodes)
    % Getting the neighbors of the free_node
    %neighbors = unique([find(neighborhood(free_nodes(ii),:)) free_nodes(ii)]);
    neighbors = find(neighborhood(free_nodes(ii),:));
    
    % Selecting the lowers closest cluster (removing the cluster 0)
    neighbors(selected_source(neighbors) == 0) = [];
    cluster_no = min(selected_source(neighbors));
    
    % Adding the free node to the cluster
    if ~isempty(cluster_no )
%         clusters{1 ,cluster_no} = ...
%             [clusters{1 ,cluster_no} free_nodes(ii)];
        selected_source(free_nodes(ii)) = cluster_no;
    end
end

if nargout > 2
    cellstruct = cell(1, max(selected_source));
    for ii = 1 : max(selected_source)
        cellstruct(ii) = {find(selected_source == ii)'};
    end
end

return
