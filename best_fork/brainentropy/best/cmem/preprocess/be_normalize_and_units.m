function [OPTIONS] = be_normalize_and_units(OPTIONS)
%BE_NORMALIZE_UNITS normalizes leadfield, baseline and data units to enhance numerical computations
%
%   INPUTS:
%       -OPTIONS     : Structure of parameters (described in be_main.m)
%          -mandatory          : Structure containing data and channel information
%               |- DataTypes              : 'MEG' or 'EEG' or 'MEEG'
%          -automatic.Modality : Structure containing the data and gain matrix for each modality
%   OUTPUTS:
%       -OPTIONS     : Keep track of parameters
%          -automatic.Modality
%               |-units: structure with new units for data and leadfield
%
% -------------------------------------------------------------------------
%   Author: LATIS 2012
%
% ==============================================
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

%% Normalization by using the mean standard deviation of baseline for each modalities
for ii = 1 : numel(OPTIONS.mandatory.DataTypes) %For every Modality (Data Type)
    
    % Std deviation for every channels on a modality
    SD = std(OPTIONS.automatic.Modality(ii).baseline');
    
    % Define the mean standard deviation (MSD) of the present modality
    MSD = mean(SD);
    
    %Normalize datas, baseline, gain matrix and EmptyRoom_data with the mean std dev
    OPTIONS.automatic.Modality(ii).data         =   OPTIONS.automatic.Modality(ii).data./MSD;
    OPTIONS.automatic.Modality(ii).baseline     =   OPTIONS.automatic.Modality(ii).baseline./MSD;
    OPTIONS.automatic.Modality(ii).gain         =   OPTIONS.automatic.Modality(ii).gain./MSD;
    OPTIONS.automatic.Modality(ii).emptyroom    =   OPTIONS.automatic.Modality(ii).emptyroom/MSD;
    OPTIONS.automatic.Modality(ii).covariance   =   OPTIONS.automatic.Modality(ii).covariance/(MSD^2);
end

%% Normalization on units
switch OPTIONS.optional.normalization
    
    case 'fixed'
        units_dipoles = 1e-9; % nAm
        for ii = 1 : numel(OPTIONS.mandatory.DataTypes)
            ratioG  =   1/max( max(OPTIONS.automatic.Modality(ii).gain) );
            OPTIONS.automatic.Modality(ii).units.Gain_units     =   ratioG;
            OPTIONS.automatic.Modality(ii).units.Data_units     =   ratioG/(units_dipoles);
            OPTIONS.automatic.Modality(ii).units.Cov_units      =   (ratioG/(units_dipoles))^2;
            OPTIONS.automatic.Modality(ii).units_dipoles             =   units_dipoles;

        end
        
    case 'adaptive'
        
        %Local fusion of data and gain matrix to compute the
        %regularisation parameter (J)        
        if numel(OPTIONS.mandatory.DataTypes)>1
            M = [OPTIONS.automatic.Modality(1).data;OPTIONS.automatic.Modality(2).data];
            G = [OPTIONS.automatic.Modality(1).gain;OPTIONS.automatic.Modality(2).gain];
        else
            M   =   OPTIONS.automatic.Modality(1).data;
            G   =   OPTIONS.automatic.Modality(1).gain;
        end
        
        if strcmp(OPTIONS.mandatory.DataTypes{1}, 'NIRS')
        J   =   be_jmne_NIRS(G,M,OPTIONS); % use depth weighting MNE for NIRS MEM 
        else
        J   =   be_jmne(G,M,OPTIONS);
        end
        
        ratioAmp = 1 / max(max(abs(J))); %Same for both modalities
        
        for ii  =   1 : numel(OPTIONS.mandatory.DataTypes)
            
            ratioG  = 1 / max(max(OPTIONS.automatic.Modality(ii).gain));
            OPTIONS.automatic.Modality(ii).units.Data_units = ratioAmp*ratioG;
            OPTIONS.automatic.Modality(ii).units.Cov_units = (ratioAmp*ratioG)^2;
            OPTIONS.automatic.Modality(ii).units.Gain_units = ratioG;
            
            OPTIONS.automatic.Modality(ii).Jmne = J * ratioAmp;
            OPTIONS.automatic.Modality(ii).ratioAmp = ratioAmp;
        end
end

%% Now, we normalize:
% - the leadfields and data
for ii = 1 : numel(OPTIONS.mandatory.DataTypes)
    OPTIONS.automatic.Modality(ii).gain = ...
        OPTIONS.automatic.Modality(ii).gain*OPTIONS.automatic.Modality(ii).units.Gain_units;
    OPTIONS.automatic.Modality(ii).data = ...
        OPTIONS.automatic.Modality(ii).data*OPTIONS.automatic.Modality(ii).units.Data_units;
    OPTIONS.automatic.Modality(ii).baseline = ...
        OPTIONS.automatic.Modality(ii).baseline*OPTIONS.automatic.Modality(ii).units.Data_units;
    
    % we check for the presence of empty room data
    if isfield(OPTIONS.automatic, 'Emptyroom_data')
        OPTIONS.automatic.Modality(ii).emptyroom = ...
            OPTIONS.automatic.Modality(ii).emptyroom*OPTIONS.automatic.Modality(ii).units.Data_units;
    end
    
    % we check the existence of a covariance matrix:
    if ~isempty(OPTIONS.automatic.Modality(ii).covariance) && (OPTIONS.solver.NoiseCov_recompute==0)
        OPTIONS.automatic.Modality(ii).covariance   =   OPTIONS.automatic.Modality(ii).covariance * ...     % BUG FIX (YZ, 29-09-2016)
            OPTIONS.automatic.Modality(ii).units.Cov_units;
%           NbjNoiseCov = size(OPTIONS.automatic.Modality(ii).covariance,3);
%         for i_sc = 1: NbjNoiseCov
%             OPTIONS.automatic.Modality(ii).covariance(:,:,i_sc) = ...
%                 OPTIONS.automatic.Modality(ii).covariance(OPTIONS.automatic.Modality(ii).channels,...
%                 OPTIONS.automatic.Modality(ii).channels,i_sc) * ...
%                 OPTIONS.automatic.Modality(ii).units.Cov_units;
%         end
    end
    
    % we check the existence of imaginary signal part:
    if isfield(OPTIONS.automatic.Modality(ii), 'idata') && ~isempty(OPTIONS.automatic.Modality(ii).idata)
        OPTIONS.automatic.Modality(ii).idata = ...
            OPTIONS.automatic.Modality(ii).idata*OPTIONS.automatic.Modality(ii).units.Data_units;
    end
end

return