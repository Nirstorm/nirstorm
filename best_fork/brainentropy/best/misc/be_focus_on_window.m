function [PLAN, params, OPTIONS] = focus_on_window( PLAN, OPTIONS, REMOV, varargin )
%FOCUS_ON_WINDOW extracts data segment corresponding to a ridge line
%
%   INPUTS:
%       -   PLAN    :   ridge plane
%       -   OPTIONS :   OPTIONS structure
%       -   REMOV   :   flag (0=set data outside window to 0, 1=delete)
%
%   OUTPUTS:
%       -   OPTIONS
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

% Optional input
iD      = 1;
if numel(varargin)==1
    iD  = varargin{1};
end
[nbL, nbS] 	=   size( PLAN );

% Get appropriate time segment
if isfield(OPTIONS, 'temporary') && isfield(OPTIONS.temporary, 'timeSMP')
    timeSMP     =   OPTIONS.temporary.timeSMP;
else
    timeSMP     =   be_closest(OPTIONS.optional.TimeSegment([1 end]), OPTIONS.mandatory.DataTime);
    OPTIONS.temporary.timeSMP = timeSMP;
end

% Get appropriate frequency range
if isfield(OPTIONS, 'temporary') && isfield(OPTIONS.temporary, 'frqSMP')
    frqSMP     =   OPTIONS.temporary.frqSMP;
else
    frqSMP      =   sort( be_closest( OPTIONS.ridges.frequency_range, fliplr(OPTIONS.wavelet.freqs_analyzed) ) );
    OPTIONS.temporary.frqSMP = frqSMP;
end

initF2      =   frqSMP(1) - 1;
nbL2        =   frqSMP(end)-frqSMP(1)+1;

if REMOV
    PLAN( :, [1:max(1,timeSMP(1)-1) min(timeSMP(end)+1,nbS):end] ) 	= [];
    PLAN( [1:max(1,frqSMP(1)-1) min(frqSMP(end)+1,nbL):end], : )  	= [];
else
    PLAN( :, [1:max(1,timeSMP(1)-1) min(timeSMP(end)+1,nbS):end] ) 	= 0;
    PLAN( [1:max(1,frqSMP(1)-1) min(frqSMP(end)+1,nbL):end], : )  	= 0;
end
PLAN( PLAN < OPTIONS.ridges.strength_threshold(iD) ) 	= 0;
if ~sum( PLAN(:) )
    PLAN = [];
end

% Output
params.initF2       =   initF2;
params.nbL          =   nbL;
params.nbL2         =   nbL2;
params.smp1         =   timeSMP(1);

return