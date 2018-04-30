function N = get_neighbor_matrix(nb_vertices, faces)
%GET_NEIGHBOR_MATRIX computes the neighborhood matrix from a surface mesh
%
%   INPUTS:
%       -   nb_vertices : total number of vertices in the mesh
%       -   faces       : mesh faces
%
%   OUTPUTS:
%       -   N
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


% Preallocation of the sparse matrix
N = sparse(nb_vertices, nb_vertices);

% Creation of the matrix
x = faces(:,1); 
y = faces(:,2); 
z = faces(:,3); 
Q1 = nb_vertices*(z-1);
N(nb_vertices*(y-1)+x) = 1;
N(Q1+x) = 1;
N(Q1+y) = 1;
N = ceil((N+N')/2);
