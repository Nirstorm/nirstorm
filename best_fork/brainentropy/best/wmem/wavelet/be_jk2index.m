function [index] = jk2index(N, j, k)
% JK2INDEX obtains the wavelet index at position (j,k) 
%
%   INPUTS:
%       -   N   :   data length (Nsamples)
%       -   j   :   scale
%       -   k   :   time
%
%   OUTPUS:
%       -   index:  coefficient index
%
%% ==============================================   
% Copyright (C) 2012 - LATIS Team
%
%  Authors: LATIS, 2012
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
index = N/2^j + k;
return