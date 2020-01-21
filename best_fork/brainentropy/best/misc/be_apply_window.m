function [OPTIONS, obj] = be_apply_window( OPTIONS, obj )
% BE_APPLY_WINDOW sets the appropriate data window
%
%   INPUTS:
%       -   obj
%       -   OPTIONS
%
%   OUTPUTS:
%       -   OPTIONS
%       -   obj
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

nMod = numel(OPTIONS.mandatory.DataTypes);

% Data
DTs = be_closest( OPTIONS.optional.TimeSegment(1), OPTIONS.mandatory.DataTime );
DTn = be_closest( OPTIONS.optional.TimeSegment(end), OPTIONS.mandatory.DataTime );

% Baseline
BSs = be_closest( OPTIONS.optional.BaselineSegment(1), OPTIONS.optional.BaselineTime );
BSn = be_closest( OPTIONS.optional.BaselineSegment(end), OPTIONS.optional.BaselineTime );

if isempty(obj) % - restrict time window, before MEM
    
    % Loop on modalities
    for ii = 1 : nMod
        OPTIONS.automatic.Modality(ii).data     = OPTIONS.automatic.Modality(ii).data(:, DTs:DTn);
        OPTIONS.automatic.Modality(ii).baseline = OPTIONS.automatic.Modality(ii).baseline(:, BSs:BSn);
    end
    OPTIONS.automatic.Comment = [OPTIONS.automatic.Comment ' | timewindow: ' num2str(OPTIONS.optional.TimeSegment(1)) ' to ' num2str(OPTIONS.optional.TimeSegment(end)) 's'];
    
else % - expand result to fit data time window, after MEM
    nr              =   zeros( obj.nb_sources, numel(OPTIONS.mandatory.DataTime) );
    nr(:, DTs:DTn)  =   obj.ImageGridAmp;
    obj.ImageGridAmp=   nr;
    
end

return