function [OPTIONS, obj] = be_main_data_preprocessing(obj, OPTIONS)
% BE_MAIN_DATA_PROCESSING sets the appropriate data and covariance units 
% according to the chosen MEM pipeline
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


obj.t0      = OPTIONS.optional.TimeSegment(1);

% ===== call of specific pre-processing depending on the pipeline:
switch OPTIONS.mandatory.pipeline
    
    case 'cMEM'  % === this is the standard MEM (Chowdhurry and Grova) =====
        [noise_var] = covariance_processing(OPTIONS);
        % we nomalize the cov of each modalities (we regularize the cov matrices)
        for ii = 1 : numel( OPTIONS.automatic.Modality )
            if isempty(OPTIONS.automatic.Modality(ii).covariance)
                OPTIONS.automatic.Modality(ii).covariance = noise_var( OPTIONS.automatic.Modality(ii).channels,OPTIONS.automatic.Modality(ii).channels );
            end
        end

    case 'wMEM'  % === this is the wavelet-MEM (Lina and co.) =============
    if OPTIONS.optional.verbose, fprintf('%s, wavelet processing not correctly called\n',OPTIONS.mandatory.pipeline); end
        
    case 'rMEM' % === this is the ridge based MEM (Zerouali and Herry) ====
        if isempty(OPTIONS.automatic.Modality.covariance)
            [noise_var] = covariance_processing(OPTIONS);
            % we normalize the cov of each modalities (we regularize the cov matrices)
            for ii = 1 : numel( OPTIONS.automatic.Modality )
                if OPTIONS.solver.NoiseCov_method == 3
                    [CovU CovD CovV] = svd(noise_var( OPTIONS.automatic.Modality(ii).channels,OPTIONS.automatic.Modality(ii).channels ));
                    OPTIONS.automatic.Modality(ii).covariance = CovU(:,1:5)*CovD(1:5,1:5)*CovV(:,1:5)';
                else
                    OPTIONS.automatic.Modality(ii).covariance = noise_var( OPTIONS.automatic.Modality(ii).channels,OPTIONS.automatic.Modality(ii).channels );
                end
            end
        end
end

return

% =========================================================================

function [noise_var] = covariance_processing(OPTIONS)
% Noise estimation from the baseline if any. Used in MEM and rMEM

if ~isempty(OPTIONS.solver.NoiseCov) && any(any( OPTIONS.solver.NoiseCov - eye( length(OPTIONS.solver.NoiseCov) ) ))    
    noise_var   = OPTIONS.solver.NoiseCov;
    
else    
    
    for ii = 1 : numel(OPTIONS.automatic.Modality)
        
        if ~isempty( OPTIONS.automatic.Modality(ii).emptyroom )
            idNb = any(isnan(OPTIONS.automatic.Modality(ii).emptyroom));
            OPTIONS.automatic.Modality(ii).emptyroom(:,idNb)	= 0;
            OPTIONS.automatic.Modality(ii).baseline          	= OPTIONS.automatic.Modality(ii).emptyroom;
            % OPTIONS.mandatory.Data                            = OPTIONS.automatic.Modality(ii).emptyroom;
            OPTIONS.automatic.Modality(ii).data             	= OPTIONS.automatic.Modality(ii).emptyroom;
            OPTIONS.mandatory.DataTime                       	= OPTIONS.automatic.Emptyroom_time;  
            fprintf('MEM : New noise covariance is computed using empty room data\n');
            
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
            fprintf('MEM : New noise covariance is computed using baseline\n');
            
        else
        	fprintf('MEM : No baseline nor covariance matrix provided. Covariance set to identity\n');
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
    
    iD          =   OPTIONS.automatic.Modality(ii).channels;
    
    % Compute covariance matrix
    switch OPTIONS.solver.NoiseCov_method
        case 0
            noise_var(iD, iD)    =   eye( size(OPTIONS.mandatory.Data,1) );
            
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
            O.wavelet.shrinkage           = 1;
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