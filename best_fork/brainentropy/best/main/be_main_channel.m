function [OPTIONS, obj] = be_main_channel(HeadModel, obj, OPTIONS)
% BE_MAIN_CHANNEL retrieves the indices of channels for each modality contained
%   in OPTIONS.DataTypes. The main objective is to fill OPTIONS.Modality with
%   the appropriate information
%
% Inputs:
% -------
%
%	HeadModel	:	structure of HeadModel used in brainstorm
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
% -------------------------------------------------------------------------
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


% ====== we initialize the Modality basic structure
% Number of modalities:
nMod = numel(OPTIONS.mandatory.DataTypes);
for ii = 1:nMod
    OPTIONS.automatic.Modality(ii).name           = OPTIONS.mandatory.DataTypes{ii};
    OPTIONS.automatic.Modality(ii).data           = [];
    OPTIONS.automatic.Modality(ii).covariance     = [];
    OPTIONS.automatic.Modality(ii).channels       = [];
    OPTIONS.automatic.Modality(ii).baseline       = [];
    OPTIONS.automatic.Modality(ii).paramH0        = [];
    OPTIONS.automatic.Modality(ii).mspDATA        = struct;
    OPTIONS.automatic.Modality(ii).emptyroom      = [];
end



%% ---- MAIN LOOP ---- ALL MODALITIIES ---- %%
cnt  = 0;
obj.nb_channels = 0;
    
for ii = 1 : nMod
    
    
    
    %% ============================ CHANNELS =========================== %%
    % retrieve MEG/EEG channels (i.e get indices of channels for all modalities)
    % the order will be the order the modalities have been listed in
    % OPTIONS.DataTypes
    
    % stand-alone
    OPTIONS.automatic.Modality(ii).channels   = find( strcmpi(OPTIONS.mandatory.ChannelTypes,OPTIONS.mandatory.DataTypes{ii}) );
    OPTIONS.automatic.Modality(ii).nchannels  = cnt + ( 1 : numel( OPTIONS.automatic.Modality(ii).channels ) )';
    cnt     = cnt + numel( OPTIONS.automatic.Modality(ii).channels );
    if isempty( OPTIONS.automatic.Modality(ii).channels )
        error(['MEM > Unable to find appropriate data. No '  OPTIONS.mandatory.DataTypes{ii} ' channels found.']);
    end
    obj.nb_channels = cnt;
    
    if isempty(HeadModel)
        CH = OPTIONS.automatic.Modality(ii).nchannels;
    else
        CH = OPTIONS.automatic.Modality(ii).channels;
    end
   
    
    %% ============================== DATA ============================= %%
    % ====== we pick the data of each modality
    % we then eliminate the OPTIONS.Data (we may introduce the iData field in
    % the OPTIONS.automatic.Modality()
    
    % Case of a ridge-filtered signal
    if isfield(OPTIONS.mandatory, 'Data') && ~isempty(OPTIONS.mandatory.Data)
        OPTIONS.automatic.Modality(ii).data = OPTIONS.mandatory.Data( CH,: );
    end
    % Case of a ridge-filtered signal
    if isfield(OPTIONS.optional, 'iData') && ~isempty(OPTIONS.optional.iData)
        OPTIONS.automatic.Modality(ii).idata          =   OPTIONS.optional.iData( CH,: );
    end
    
    % Case of a wavelet-adaptive clustering
    if isfield( OPTIONS.automatic, 'mspDATA' ) && ~isempty( fieldnames(OPTIONS.automatic.mspDATA) )
        OPTIONS.automatic.Modality(ii).mspDATA.FRQs   =   OPTIONS.automatic.mspDATA.FRQs;
        OPTIONS.automatic.Modality(ii).mspDATA.Time   =   OPTIONS.automatic.mspDATA.Time;
        OPTIONS.automatic.Modality(ii).mspDATA.F      =   OPTIONS.automatic.mspDATA.F( CH,: );
        OPTIONS.automatic   =   rmfield( OPTIONS.automatic,  'mspDATA' );
    end
    
    
    %% ============================ BASELINE =========================== %%
    % ====== we pick the baseline (if any) of each modality
    % we then eliminate the OPTIONS.Baseline
    if ~isempty(OPTIONS.optional.Baseline) 
        OPTIONS.automatic.Modality(ii).baseline = OPTIONS.optional.Baseline( CH,: );
    end


    %% ============================ EMPTY ROOM =========================== %%
    % ====== we pick the empty room data (if any) of each modality    
    if isfield(OPTIONS.automatic, 'Emptyroom_data') 
        ERD     =   OPTIONS.automatic.Emptyroom_data( CH,: );
        if any(ERD(:))
            % emptyroom data available for this modality
            OPTIONS.automatic.Modality(ii).emptyroom = ERD;
        end
    end
    
    
    %% =========================== COVARIANCE ========================== %%
    % ====== we pick the covariance (if any) of each modality
    % we then eliminate the OPTIONS.Baseline
    if ~isempty(OPTIONS.solver.NoiseCov)
        NbjNoiseCov = size(OPTIONS.solver.NoiseCov,3);
        for i_sc = 1: NbjNoiseCov
            OPTIONS.automatic.Modality(ii).covariance(:,:,i_sc) = ...
                OPTIONS.solver.NoiseCov(CH,CH,i_sc );
        end
    end
    
    
    %% ============================== GAIN ============================= %%
    % ====== We pick the gain matrix of each modality 
    % the gain matrix is associated with each modality; we no more keep it in
    % obj
    if isfield(HeadModel, 'Gain') && ~isempty(HeadModel.Gain)  % stand-alone
        % WE PICK THE GAIN MATRIX IN THE ORDER OF THE MODALITIES INDICATED IN
        % OPTIONS.DATATYPES. IN THE HEADMODEL, WE INDICATE THE TYPE OF THE
        % MATIX
        
        Gain = HeadModel.Gain(CH,:);
        
        if isstruct(HeadModel) && isfield(HeadModel, 'GridOrient')
            Gain = bst_gain_orient(Gain, HeadModel.GridOrient);
        end
        OPTIONS.automatic.Modality(ii).gain = Gain;
        
        % average reference the gain matrix
        if strcmpi(OPTIONS.mandatory.DataTypes{ii}, 'eeg')
            avgref  =   mean( OPTIONS.automatic.Modality(ii).gain );
            OPTIONS.automatic.Modality(ii).gain     =   OPTIONS.automatic.Modality(ii).gain - ...
                                                        ones( numel(CH),1 ) * avgref;
        end        
     end      
    
    
    %% ============================ MSP data =========================== %%
    % ====== we pick the temporal data needed for MSP
    % This only applies to rMEM where OPTIONS.mandatory.Data be_main_channel.mis complex
    if ~numel( fieldnames(OPTIONS.automatic.Modality(ii).mspDATA) ) && isfield( OPTIONS, 'temporary' ) && isfield( OPTIONS.temporary, 'mspDATA' )
        OPTIONS.automatic.Modality(ii).mspDATA      = OPTIONS.temporary.mspDATA;
        OPTIONS.automatic.Modality(ii).mspDATA.F    = OPTIONS.automatic.Modality(ii).mspDATA.F( CH, : );
    end    
    
    
end

