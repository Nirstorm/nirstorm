function [W,C] = be_whittening_dewhittening_matrices(data)
%BE_DEWHITTENING_MATRICES computes dewhitening matrix for a given dataset
%
%   INPUTS:
%       -   data    : data matrix
%
%   OUTPUTS:
%       -   W       : weights vector
%       -   C       : dewhitening matrix
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

Nt = size(data,2);

data = data-repmat(mean(data,2),1,Nt);
COV  = data*data'/(Nt-1);

[U,D,V] = svd(COV);
W = diag(sqrt(1./diag(D)))*U';
C = U*diag(sqrt(diag(D)));