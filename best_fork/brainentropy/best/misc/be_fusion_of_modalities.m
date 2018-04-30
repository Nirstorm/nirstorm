function [obj] = be_fusion_of_modalities(data, obj, OPTIONS)
%BE_FUSION_OF_MODALITIES fuses data and leadfields from EEG and MEG for 
% multimodal sources estimation using MEM
%
%   INPUTS:
%       -   data
%       -   obj
%       -   OPTIONS
%
%   OUTPUTS:
%       -   OPTIONS
%       - obj
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


if isempty(data)
    data        = {OPTIONS.automatic.Modality.data};
    obj.data    = data{1};

end
if isfield(OPTIONS.automatic.Modality(1),'idata')
    idata       = {OPTIONS.automatic.Modality.idata};
    obj.idata   = idata{1};
end

if isfield( OPTIONS.automatic.Modality(1), 'gain' )
    obj.gain        = OPTIONS.automatic.Modality(1).gain;
end

% If the covariance matrix is scale dependent, no multimodality possible
if size(OPTIONS.automatic.Modality(1).covariance,3)>1
    obj.noise_var = OPTIONS.automatic.Modality(1).covariance;
    if length(OPTIONS.mandatory.DataTypes)>1
        error('No multimodality possible with a scale dependent noise_cov');
    end
else
    obj.noise_var   = diag( OPTIONS.automatic.Modality(1).covariance );
end
    
obj.baseline    = OPTIONS.automatic.Modality(1).baseline;
obj.channels    = OPTIONS.automatic.Modality(1).channels;

if length(OPTIONS.mandatory.DataTypes)>1 % fusion of modalities if requested
    if OPTIONS.optional.verbose
        fprintf('%s, MULTIMODAL data ... ',OPTIONS.mandatory.pipeline);
    end
    for ii=2:length(OPTIONS.mandatory.DataTypes)
        obj.data = [obj.data ; data{ii}];
        if exist('idata', 'var'); obj.idata = [obj.idata; idata{ii}]; end
        data{ii}        = []; 
        idata{ii}       = [];

        if isfield( OPTIONS.automatic.Modality(1), 'gain' )
            obj.gain        = [obj.gain ; OPTIONS.automatic.Modality(ii).gain];
        end
        obj.noise_var   = [obj.noise_var; diag(OPTIONS.automatic.Modality(ii).covariance)];
        obj.baseline    = [obj.baseline; OPTIONS.automatic.Modality(ii).baseline];
        obj.channels    = [obj.channels; OPTIONS.automatic.Modality(ii).channels];
        if OPTIONS.optional.verbose
            fprintf('... %s found ', OPTIONS.mandatory.DataTypes{ii})
        end
    end

else
    if OPTIONS.optional.verbose
        fprintf('%s, No multimodalities ...',OPTIONS.mandatory.pipeline);
    end
end

% If the covariance matrix is scale dependent, noise_var is already a matrix
if size(OPTIONS.automatic.Modality(1).covariance,3)==1
    obj.noise_var = diag( obj.noise_var );
end

if OPTIONS.optional.verbose, fprintf('\n'); end