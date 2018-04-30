function [OPTIONS, obj] = be_wdata_preprocessing(obj, OPTIONS)
% BE_WDATA_PREPROCESSING wavelet transforms data and sets the appropriate 
% data and covariance units accroding to the chosen MEM pipeline
%
%   INPUTS:
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

obj.t0      = OPTIONS.mandatory.DataTime(1);

if OPTIONS.optional.verbose
    if ~strcmp(OPTIONS.mandatory.pipeline,'wMEM')
        fprintf('%s, wavelet processing not correctly called\n',OPTIONS.mandatory.pipeline);
    else
        fprintf('%s, wavelet pre-processing (new)\n',OPTIONS.mandatory.pipeline);
    end
end

% ====  this is the wavelet-MEM (Lina and co.)
        % we first re-organize the data with respect to the modalities
        % concerned (the order being given by the OPTIONS.DataTypes)
        % normalization/wavelet/denoise
        [obj, OPTIONS] = be_discrete_wavelet_preprocessing(obj, OPTIONS);

        [noise_var] = covariance_processing(OPTIONS);

        % The number of j decomposition in the noise_var
        % noise_var:[NbSensors x NbSensors x Nbj]
        noise_var_jn = size(noise_var,3);

        % we nomalize the cov of each modalities (we regularize the cov matrices)
        for ii = 1 : numel( OPTIONS.automatic.Modality )
            for isc=1:noise_var_jn
                OPTIONS.automatic.Modality(ii).covariance(:,:,isc) = ...
                    noise_var(OPTIONS.automatic.Modality(ii).channels,...
                    OPTIONS.automatic.Modality(ii).channels, isc );
            end
        end
return

% =========================================================================

function [noise_var] = covariance_processing(OPTIONS)
% Noise estimation from the baseline if any. Used in MEM and rMEM
% OPTIONS is locally modified, it is not an output of the function

% If a noise cov is loaded
if ~isempty(OPTIONS.solver.NoiseCov)
    % If the cov mat is scale dependent or it is not identity
    if length(size(OPTIONS.solver.NoiseCov))==3 || ...
            any(any( OPTIONS.solver.NoiseCov - eye( length(OPTIONS.solver.NoiseCov) ) ))
        noise_var   = OPTIONS.solver.NoiseCov;
    end
else
    
    for ii = 1 : numel(OPTIONS.automatic.Modality)
        
        if ~isempty( OPTIONS.automatic.Modality(ii).emptyroom )
            idNb = any(isnan(OPTIONS.automatic.Modality(ii).emptyroom));
            OPTIONS.automatic.Modality(ii).emptyroom(:,idNb)	= 0;
            OPTIONS.automatic.Modality(ii).baseline          	= 0;
            OPTIONS.automatic.Modality(ii).data             	= OPTIONS.automatic.Modality(ii).emptyroom;
            OPTIONS.mandatory.DataTime                       	= OPTIONS.automatic.Emptyroom_time;
            fprintf('%s, New noise covariance is computed using empty room data\n',OPTIONS.mandatory.pipeline);
            
        elseif ~isempty( OPTIONS.automatic.Modality(ii).baseline )
            idNb = any(isnan(OPTIONS.automatic.Modality(ii).baseline));
            OPTIONS.automatic.Modality(ii).baseline(:,idNb)     = 0;
            OPTIONS.automatic.Modality(ii).data                 = OPTIONS.automatic.Modality(ii).baseline;
            if ~isempty( OPTIONS.optional.BaselineTime )
                OPTIONS.mandatory.DataTime                      = OPTIONS.optional.BaselineTime;
            else
                TLM                                             =   be_closest(OPTIONS.optional.BaselineSegment([1 end]), OPTIONS.mandatory.DataTime);
                OPTIONS.mandatory.DataTime                      = OPTIONS.mandatory.DataTime( TLM(1):TLM(end) );      
            end
            fprintf('%s, New noise covariance is computed using baseline\n',OPTIONS.mandatory.pipeline);
            
        else
        	fprintf('%s, No baseline nor covariance matrix provided. Covariance set to identity\n',OPTIONS.mandatory.pipeline);
            OPTIONS.automatic.Modality(ii).baseline             =   eye(OPTIONS.automatic.Modality(ii).channels);
            OPTIONS.solver.NoiseCov_method                      =   0;
        end
        
    end

    OPTIONS.optional.TimeSegment        =   OPTIONS.mandatory.DataTime( [1 end] );
    noise_var                           =   estimate_noise_var(OPTIONS); % NEW. Bug fix TH 30/12/2013   
    
end
    
return

% ====

function [noise_var] = estimate_noise_var(OPTIONS)

    noise_var = [];

    for ii =  1: numel( OPTIONS.automatic.Modality )
    
        iD =   OPTIONS.automatic.Modality(ii).channels;
    
        % Compute covariance matrix

        % Wavelet-based noise estimation
        O = OPTIONS;

        O.automatic.Modality        = O.automatic.Modality(ii);            
        O.mandatory.DataTypes       = O.mandatory.DataTypes(ii);

        % Construct wavelet
        O.wavelet.type                = 'rdw';
        O.wavelet.vanish_moments      = 4;
        O.wavelet.selected_scales     = [];
        O.wavelet.single_box          = 0;
        O.automatic.scales            = [];

        O.optional.verbose            = 0;

        % Make transform
        obj.t0              =   OPTIONS.mandatory.DataTime(1);
        wavelet_obj         =   be_discrete_wavelet_preprocessing(obj, O);
        % If the empty room is available the cov mat is scale dependent
        if ~isempty( OPTIONS.automatic.Modality(ii).emptyroom )
            for isc=1:size(OPTIONS.automatic.scales,2)
                variance = var( wavelet_obj.data{1}(:,end/(2^isc)+1:end/(2^(isc-1)))');
                switch OPTIONS.solver.NoiseCov_method
                    % Diagonale
                    case 4
                        noise_var(iD, iD, isc) = diag(variance);
                    % Diagonale averaged
                    case 5
                        noise_var(iD, iD, isc) = diag(ones(1,length(variance))*mean(variance));
                    % If the NoiseCov_method is not 4 neither 5 we force 5
                    otherwise
                        noise_var(iD, iD, isc) = diag(ones(1,length(variance))*mean(variance));
                end
            end
        % If not a unique cov matrix is computed from the baseline
        else
            variance = var( wavelet_obj.data{1}(:,end/2+1:end)' );
            switch OPTIONS.solver.NoiseCov_method
                % Diagonale
                case 4
                    noise_var(iD, iD) = diag(variance);
                % Diagonale averaged    
                case 5
                    noise_var(iD, iD) = diag(ones(1,length(variance))*mean(variance));
               % If the NoiseCov_method is not 4 neither 5 we force 5
                otherwise
                    noise_var(iD, iD) = diag(ones(1,length(variance))*mean(variance)); 
            end
        end
    end

return
