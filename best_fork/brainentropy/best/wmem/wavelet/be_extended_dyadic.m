function [s, info] = be_extended_dyadic( x )
% dyadic extent of the data
%   INPUTS: 
%       -   x       :   initial data
%
%   OUTPUTS: 
%       -   s       :   (extended signal)
%       -   info    :   (structure that keeps the information of the transform)
%
%% ==============================================   
% Copyright (C) 2012 - LATIS Team
%
%  Authors: JM Lina, 2012, jan 1st
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


N = ceil(log2(size(x,2)));
s = zeros(size(x,1),2^N);
offset = floor((2^N-size(x,2))/2);
s(:,1:1+offset)=repmat(x(:,1),1,1+offset);
s(:,1+offset:offset+size(x,2))=x;
s(:,offset+size(x,2)-1:end)=repmat(x(:,end),1,2^N-offset-size(x,2)+2);
info.start = 1+offset;
info.end = offset+size(x,2);
end


