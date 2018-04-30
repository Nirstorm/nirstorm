function [OPTIONS] = be_rdg_threshold( ridges, OPTIONS, iD, defT )
%BE_RDG_THRESHOLD Returns the ridge selection threshold for the data signal 
% that includes a baseline. Threshold computation is described in accompanying 
% paper.
%
%   Example:
%       thresh = be_rdg_threshod( baseline );
%
%   Example:
%       data = be.c_data(alldata, SF, t0);
%       data.markers{1} = be.c_marker(t, 'after', .5);
%       baseline = data.subsignals(1);
%       thresh = rdg_threshod( baseline );
%
%   INPUTS:
%       -   ridges  :   ridge plane
%       -   OPTIONS :   options structure
%       -   iD      :   number of the current modality
%
%   Ref.:
%       Zerouali, Herry, Jemel, Lina (2011) Localization of cortical 
%       synchronous neural sources, IEEE Trans. Biomed. Eng.
%
%   See also be.c_marker, be.c_data.
%
%   Reference page in Help browser
%       <a href = "matlab: be.matlabdoc.help('+be/@c_data/rdg_thershold')">
%       matlabdoc be.c_data.rdg_thershold</a>
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


% Import LATIS libraries
import be.programming.*

% Statistical threshold 
thresh  = .9;

% Extracted ridge lines
[rdgL]    = OPTIONS.automatic.Modality(iD).ridges_baseline;
    
% Default threshold
T   =   0.2;
if exist('defT', 'var')
    T   =   defT;
end

% Compute strengths
if ~isempty(rdgL)
    STRS    =   cellfun( @(a) mean(ridges(a)), rdgL, 'uni', false );
    STRS    =   sort( [STRS{:}] );
    [y,x]   =   hist( STRS, max(5,numel(rdgL)/10) );
    STRSc   =   cumsum( y ) / sum( y );
    VAL1    =   find(STRSc<thresh);
    VAL2    =   find(STRSc>=thresh);
    if isempty(VAL1)
        VAL1    =   VAL2(1);
        VAL2    =   VAL2(2:end);
    end
    if ~isempty(VAL2)
        slope   =   diff( STRSc([VAL2(1) VAL1(end)]) ) / diff( x([VAL2(1) VAL1(end)]) );
        T       =   ( thresh - STRSc(VAL2(1)) ) / slope + x(VAL2(1));
    end
end

% Apply threshold
idR         =   cellfun( @(a) median(ridges(a)), OPTIONS.automatic.Modality(iD).ridges_data, 'uni', 0 );
idR         =   cellfun( @(a) a<T, idR, 'uni', 0 );
idR         =   cell2mat( idR );
OPTIONS.automatic.Modality(iD).ridges_data(idR) =   [];
OPTIONS.ridges.strength_threshold(iD)           =   T;

% Store baseline
% BS  =   OPTIONS.automatic.Modality(iD).baseline;
% BS  =   be_bpfilter( BS, OPTIONS.automatic.sampling_rate, OPTIONS.ridges.frequency_range );    
% OPTIONS.automatic.Modality(iD).baseline = BS;

return