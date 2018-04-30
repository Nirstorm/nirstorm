function [OPTIONS, obj] = be_main_leadfields(obj, OPTIONS)
% BE_MAIN_LEADFIELDS preprocesses the lead fields for faster MSP computation
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

% We normalize and compute the singular values and eigenvectors of the lead
% fields.

if OPTIONS.optional.verbose
    fprintf('%s, be_main_leadfields ...', OPTIONS.mandatory.pipeline); 
end 

for ii = 1 : numel(OPTIONS.mandatory.DataTypes)
    OPTIONS.automatic.Modality(ii).gain = OPTIONS.automatic.Modality(ii).gain(:, obj.iModS);
    OPTIONS.automatic.Modality(ii).gain_struct = be_decompose_gain(OPTIONS.automatic.Modality(ii).gain);
end

if OPTIONS.optional.verbose
    fprintf(' done.\n'); 
end 


return
