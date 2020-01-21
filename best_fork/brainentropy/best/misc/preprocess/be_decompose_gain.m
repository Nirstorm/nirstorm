function [Gstruct] = be_decompose_gain(G)
% Normalize the forward operator. 
% Output: a structure with
% - the normalized lead field matrix
% - the orthogonal matrix U (size ncap x ncap)
% - the singular values (size ncap x 1)
% ==============================================
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

Gstruct = struct;

if isempty(G)
    return
end

Gn = bsxfun(@rdivide, G, sqrt(sum(G.^2, 1)));
[temp,lambda,U] = svd(Gn',0);
Gstruct.Gn = Gn;
Gstruct.lambda = diag(lambda);
Gstruct.U = U;

return
