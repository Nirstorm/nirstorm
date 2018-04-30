function [OPTIONS, obj] = be_main_sources(obj, OPTIONS)
% BE_MAIN_SOURCES removes sources for which the leadfield contains Nan
% The gain matrix(ces) is in OPTIONS.Modality{}.gain
% There is at least one modality.
%
% Inputs:
% -------
%
%	obj			:	MEM obj structure
%   OPTIONS     :   structure (see bst_sourceimaging.m)
%
%
% Outputs:
% --------
%
%   OPTIONS     :   Updated options fields
%	obj			:	Updated structure
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

if OPTIONS.optional.verbose
    fprintf('%s, be_main_source ...', OPTIONS.mandatory.pipeline); 
end 

obj.nb_dipoles = size(OPTIONS.automatic.Modality(1).gain,2);
obj.iModS = 1:obj.nb_dipoles;

% remove sources with bad leadfields
for ii = 1 : numel(OPTIONS.mandatory.DataTypes)
    [iS] = be_check_gain( OPTIONS.automatic.Modality(ii).gain, OPTIONS.mandatory.DataTypes(ii) );
    obj.iModS = intersect(obj.iModS, iS);
end

% what we keep:
obj.nb_sources = numel(obj.iModS);

if OPTIONS.optional.verbose
    fprintf(' done.\n'); 
end  

return