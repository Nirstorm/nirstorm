function [time] = jk2time(fs, j, k)
% JK2TIME obtains the time (0 being the start) at position (j,k) 
%
%   INPUTS:
%       -   fs  :   sampling frequency
%       -   j   :   scale
%       -   k   :   time (coefficient)
%
%   OUTPUS:
%       -   time:  time (signal)
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


% 
dt = 1/fs;
time = (2^j-1)*dt/2  + (k-1)*2^j*dt;
return