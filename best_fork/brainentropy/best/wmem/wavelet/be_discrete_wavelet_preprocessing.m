function [obj, OPTIONS] = be_discrete_wavelet_preprocessing(obj, OPTIONS)
% BE_DISCRETE_WAVELET_PREPROCESSING prepares data for discrete wavelet
% transform and shrinkage
%
%   INPUTS:
%       -   obj :   data structure
%       -   OPTIONS :   options structure (see be_main.m)
%
%   OUTPUTS:    
%       -   obj
%       -   OPTIONS
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

obj.data_type = 'discrete_wavelet';
% we then wavelet transform the data (that are xtended to the next power 2)
for ii = 1 : numel(OPTIONS.mandatory.DataTypes)
    [obj.data{ii}, obj.info_extension] = ...
        be_extended_dyadic(OPTIONS.automatic.Modality(ii).data);
    [obj.data{ii}, obj.scaling_data{ii}, OPTIONS] = ...
        be_discrete_wavelet_transform(obj.data{ii}, OPTIONS);
end

% we wavelet denoise (if selected) 
% Equation (3) in ref.[1]
if OPTIONS.wavelet.shrinkage
    for ii = 1 : numel(OPTIONS.mandatory.DataTypes);
        if OPTIONS.optional.verbose, fprintf('%s, wavelet processing of the %s\n',...
                OPTIONS.mandatory.pipeline,OPTIONS.mandatory.DataTypes{ii}); 
        end        
        % 2 ways to shrink
        % If there's Emptyroom_data -> AminGhafari method
        % otherwise use the old shrink
        if ~isempty(OPTIONS.automatic.Modality(ii).emptyroom)
        	[obj.data{ii}, OPTIONS] = be_denoise_csoft_multivar(...
                obj.data{ii},obj.info_extension,obj.scaling_data{ii},OPTIONS.automatic.Modality(ii).emptyroom, OPTIONS);
        else
            [W,C] = be_whittening_dewhittening_matrices(OPTIONS.automatic.Modality(ii).data);
            [obj.data{ii}, OPTIONS] = be_wdenoise_csoft(obj.data{ii}, obj.info_extension, W, C, OPTIONS);
        end
    end
end

% selection of boxes with power
[OPTIONS] = be_selected_coeff(obj.data, obj, OPTIONS);
return
