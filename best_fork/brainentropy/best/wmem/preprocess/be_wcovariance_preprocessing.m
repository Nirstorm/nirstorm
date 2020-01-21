function [OPTIONS, obj] = be_Wcov_preprocessing(obj, OPTIONS)
% BE_MAIN_DATA_PROCESSING sets the covariance units 
% accroding to the wMEM pipeline
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

% ===== call of specific pre-processing depending on the pipeline:
switch OPTIONS.mandatory.pipeline
    
    case 'cMEM'  % === this is the standard MEM (Chowdhurry and Grova) ===== 
        [noise_var] = covariance_processing(OPTIONS);
        % we nomalize the cov of each modalities (we regularize the cov matrices)
        for ii = 1 : numel( OPTIONS.automatic.Modality )
            OPTIONS.automatic.Modality(ii).covariance = noise_var( OPTIONS.automatic.Modality(ii).channels,OPTIONS.automatic.Modality(ii).channels ) * OPTIONS.automatic.Modality(ii).units.Cov_units;
        end

        
    case 'wMEM'  % === this is the wavelet-MEM (Lina and co.) =============
        % we first re-organize the data with respect to the modalities
        % concerned (the order being given by the OPTIONS.DataTypes)
        % normalization/wavelet/denoise 
        [obj, OPTIONS] = be_discrete_wavelet_preprocessing(obj, OPTIONS);
        % we nomalize the cov of each modalities (we regularize the cov matrices)
        for ii = 1 : numel(OPTIONS.mandatory.DataTypes)
            W1 = obj.data{ii}(:,size(obj.data{ii},2)/2+1:end);
            OPTIONS.automatic.Modality(ii).covariance = diag(var(W1')) * OPTIONS.automatic.Modality(ii).units.Cov_units;
            obj.data{ii} = obj.data{ii} * OPTIONS.automatic.Modality(ii).units.Data_units;
            obj.scaling_data{ii} = obj.scaling_data{ii} * OPTIONS.automatic.Modality(ii).units.Data_units;
        end
        
    case 'rMEM' % === this is the ridge based MEM (Zerouali and Herry) ====
        [noise_var] = covariance_processing(OPTIONS);
        % we nomalize the cov of each modalities (we regularize the cov matrices)
        for ii = 1 : numel( OPTIONS.automatic.Modality )
            if OPTIONS.solver.NoiseCov_method == 3
                [CovU CovD CovV] = svd(noise_var( OPTIONS.automatic.Modality(ii).channels,OPTIONS.automatic.Modality(ii).channels ));
                OPTIONS.automatic.Modality(ii).covariance = CovU(:,1:5)*CovD(1:5,1:5)*CovV(:,1:5)'*OPTIONS.automatic.Modality(ii).units.Cov_units;
            else
                OPTIONS.automatic.Modality(ii).covariance = noise_var( OPTIONS.automatic.Modality(ii).channels,OPTIONS.automatic.Modality(ii).channels ) * OPTIONS.automatic.Modality(ii).units.Cov_units;
            end
        end
       
end

% ======= data processing: 
% We first normalize data units consitently 
% with the leadfields normalization 
for ii = 1 : numel( OPTIONS.automatic.Modality )
    %OPTIONS.Data(obj.iMod{ii},:) = OPTIONS.Data(obj.iMod{ii},:)*OPTIONS.Units(ii).Data_units;
    OPTIONS.automatic.Modality(ii).data = OPTIONS.automatic.Modality(ii).data*OPTIONS.automatic.Modality(ii).units.Data_units;
    if isfield(OPTIONS.automatic.Modality(ii), 'idata')&& ~isempty(OPTIONS.automatic.Modality(ii).idata)  % A MODIFIER <--------------------------------------------
        OPTIONS.automatic.Modality(ii).idata = OPTIONS.automatic.Modality(ii).idata*OPTIONS.automatic.Modality(ii).units.Data_units; % <---------------
    end
end
% ========================================

return

% =========================================================================

function [noise_var] = covariance_processing(OPTIONS)
% Noise estimation from the baseline if any. Used in MEM and rMEM

if ~isempty(OPTIONS.solver.NoiseCov) && any(any( OPTIONS.solver.NoiseCov - eye( length(OPTIONS.solver.NoiseCov) ) ))    
    noise_var   = OPTIONS.solver.NoiseCov;
    
elseif ~isempty(OPTIONS.optional.Baseline)    
    
    for ii = 1 : numel(OPTIONS.automatic.Modality)
        idNd = any(isnan(OPTIONS.automatic.Modality(ii).data));
        idNb = any(isnan(OPTIONS.automatic.Modality(ii).baseline));
        OPTIONS.automatic.Modality(ii).baseline(:,idNb) = 0;
        OPTIONS.automatic.Modality(ii).data = OPTIONS.automatic.Modality(ii).baseline;
        OPTIONS.mandatory.Data              = OPTIONS.automatic.Modality(ii).baseline;
        if ~isempty( OPTIONS.optional.BaselineTime )
            OPTIONS.mandatory.DataTime  = OPTIONS.optional.BaselineTime;
        else
            TLM     =   be_closest(OPTIONS.optional.BaselineSegment([1 end]), OPTIONS.mandatory.DataTime);
            OPTIONS.mandatory.DataTime  = OPTIONS.mandatory.DataTime( TLM(1):TLM(end) );      
        end
        OPTIONS.optional.TimeSegment    = OPTIONS.mandatory.DataTime( [1 end] );        
    end
    
    sprintf('MEM : New noise covariance is computed using baseline\n');
    noise_var   = estimate_noise_var(OPTIONS);
    
else
    sprintf('MEM : No baseline nor covariance matrix provided. Covariance set to identity\n');
    OPTIONS.solver.NoiseCov_method  =   0;
    noise_var   = estimate_noise_var(OPTIONS);
    
end
noise_var   =   noise_var / OPTIONS.solver.covariance_scale;
return

% ====

function [noise_var] = estimate_noise_var(OPTIONS)


noise_var = [];

for ii =  1: numel( OPTIONS.automatic.Modality )
    
    iD          =   OPTIONS.automatic.Modality(ii).channels;
    
    % Compute covariance matrix
    switch OPTIONS.solver.NoiseCov_method
        case 0
            noise_var(iD, iD)    =   eye( numel(OPTIONS.automatic.GoodChannel) );
            
        case 1
            noise_var(iD, iD)    =   diag( mean(diag( real( cov(OPTIONS.automatic.Modality(ii).baseline') ) ) ) * ones( numel(OPTIONS.automatic.Modality(ii).channels),1 ) );
            
        case 2
            noise_var(iD, iD)    =   diag( diag( real( cov(OPTIONS.automatic.Modality(ii).baseline') ) ) );
            
        case 3
            noise_var(iD, iD)    =   real( cov(OPTIONS.automatic.Modality(ii).baseline') );
            
        case 4
            % Wavelet-based noise estimation
            O = OPTIONS;
            
            O.automatic.Modality        = O.automatic.Modality(ii);            
            O.mandatory.DataTypes       = O.mandatory.DataTypes(ii);
            
            % Construct wavelet
            O.wavelet.type                = 'RDW';
            O.wavelet.vanish_moments      = 4;
            O.wavelet.shrinkage           = 0;
            O.optional.verbose            = 0;
            O.wavelet.selected_scales     = [];
            O.wavelet.single_box          = 0;
            
            % Make transform
            obj.t0              =   OPTIONS.mandatory.DataTime(1);
            wavelet_obj         =   be_wavelet_preprocessing_data(obj, O);
            noise_var(iD, iD)   =   diag( var( abs( wavelet_obj.data{1}(:,end/2+1:end)' ) ) );
            
    end

end

return
