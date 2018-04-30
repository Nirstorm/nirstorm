function [OPTIONS, W] = be_spatial_priorw(OPTIONS, neighbors)
%   This function returns the W Green matrix from which the local 
%   covariance matrices will be obtained. 
%
%   INUPTS:
% 		- OPTIONS    	: structure of parameters
%       - neighbors		: neighbor matrix (with 0s on the diag)
%
%   OUTPUTS:
%		- OPTIONS		: Keep track of parameters
%		- W				: Green Matrix
%   
% 	Reference
%       Harrison et al., NIMG 38(4), 677-695 (2007
%       Friston et al., NIMG 39(3), 1104-1120 (2008)
%
%   Authors: LATIS team, 2011.
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
                                              
                                              
% Default options settings
Def_OPTIONS.solver = struct('spatial_smoothing', 0.6);

% Return empty OPTIONS structure
if (nargin == 0)
    OPTIONS = Def_OPTIONS;
    return
end
% Check field names of passed OPTIONS and fill missing ones with default values
OPTIONS = be_struct_copy_fields(OPTIONS, Def_OPTIONS, {'solver'}, 0);
clear Def_OPTIONS
%%
% Parameters:
nb_vertices = size(OPTIONS.automatic.Modality(1).gain,2);
rho = OPTIONS.solver.spatial_smoothing; % scalar that weight the adjacency matrix 
W   = speye(nb_vertices);


% Add comment to result
if rho && ~isempty(neighbors)
    OPTIONS.automatic.Comment = [OPTIONS.automatic.Comment ' | smooth=' num2str(rho)];
else
    return
end

% Preallocation of the sparse matrix
A  = neighbors - spdiags(sum(neighbors,2),0,nb_vertices,nb_vertices);
A0 = rho*A/2;
for i = 1:7
    W = W + A0;
    A0 = rho/2*A0*A / (i+1);
end
W = W.*(W > exp(-8));
W = W'*W; 

return