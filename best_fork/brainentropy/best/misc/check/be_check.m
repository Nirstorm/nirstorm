function stand_alone = check_caller
% Check is main script is called from BST environment

stand_alone = 1;

% Check caller
caller  =   dbstack(1);
isglob  =   whos('global', 'GlobalData');

if ~isempty(strfind( {caller.name}, 'bst_sourceimaging' )) && ~isempty(isglob)
    stand_alone =  0;
end

return               
        
function [OPTIONS, isRF]    = checkDataPipeline(OPTIONS)

% Check data type
warning('OFF')
isRF    =   0;

% Load whole data - not optimized for standalone call
try
    load( be_fullfile(OPTIONS.automatic.iProtocol.STUDIES, OPTIONS.optional.DataFile), 'History', 'F', 'Time');
    OPTIONS.mandatory.DataTime  =   Time;
    OPTIONS.mandatory.Data      =   F;
    isRF    =   strcmp( History(:,2), 'Ridge-based filtering' );

catch
    if ~strcmpi( OPTIONS.optional.DataFile, {'', 'data'} )
        fprintf('In be_main :\tData file cannot be loaded. Check baseline definition\n');
    end
    
end

% Only useful is signal is already ridgefiltered
if ~any(isRF)
    OPTIONS     =   PIPELINEoptions( OPTIONS );
    
elseif strcmp(OPTIONS.mandatory.pipeline,'wMEM')
    OPTIONS     =   PIPELINEoptions( OPTIONS );
    OPTIONS.mandatory.pipeline    = 'rwMEM';
    
elseif any(isRF)
    OPTIONS.mandatory.pipeline    = 'rMEM';
    OPTIONS.automatic.rMEMfiles   = { strrep(OPTIONS.optional.DataFile, OPTIONS.automatic.iProtocol.STUDIES, '') };
    OPTIONS.mandatory.Data        = [];
    
    % Check if ridges are computed for all modalities 
    iD = find( ~any( ~isnan(F(OPTIONS.automatic.GoodChannel,:)),2 ) );
    if ~isempty(iD)
        NMS = unique( OPTIONS.mandatory.ChannelTypes(iD) );
        NMS = cellfun(@(a) [a ' '], NMS, 'uni', false);
        error(['rMEM : Ridges were not computed for ' NMS{:}]) 
    end
    isRF    =   sum(isRF);
    OPTIONS =   PIPELINEoptions( OPTIONS );
    
end

return

function [OPTIONS]          = check_timedef(OPTIONS, isRF)

% Get baseline data segment
OPTIONS.automatic.sampling_rate     =   round( 1 / diff( OPTIONS.mandatory.DataTime([1 2]) ) );

% Fills MSP_data field
if isempty(OPTIONS.optional.TimeSegment)
    OPTIONS.optional.TimeSegment    =   OPTIONS.mandatory.DataTime([1 end]);
end
OPTIONS.optional.TimeSegment        =   be_closest( OPTIONS.optional.TimeSegment([1 end]), OPTIONS.mandatory.DataTime );
OPTIONS.optional.TimeSegment        =   OPTIONS.mandatory.DataTime(OPTIONS.optional.TimeSegment(1):OPTIONS.optional.TimeSegment(end));

% No baseline
if isempty(OPTIONS.optional.Baseline) && ~any(isRF)
    % Baseline segment definition
    OPTIONS.optional.BaselineTime   =   OPTIONS.mandatory.DataTime;
    OPTIONS.optional.Baseline       =   OPTIONS.mandatory.Data;
    if isempty(OPTIONS.optional.BaselineSegment)
        OPTIONS.optional.BaselineSegment = OPTIONS.mandatory.DataTime(1):1/OPTIONS.automatic.sampling_rate:OPTIONS.mandatory.DataTime(end);
    else
        STb     = be_closest( OPTIONS.optional.BaselineSegment(1), OPTIONS.optional.BaselineTime );
        NDb     = be_closest( OPTIONS.optional.BaselineSegment(end), OPTIONS.optional.BaselineTime );
        OPTIONS.optional.Baseline       =   OPTIONS.optional.Baseline(:, STb:NDb);
        OPTIONS.optional.BaselineTime   =   OPTIONS.optional.BaselineTime(STb:NDb);
    end
    
elseif ~any(isRF)
    % Baseline time definition
    OPTIONS.automatic.BaselineType  =   'independent';
    if isempty(OPTIONS.optional.BaselineTime)
        OPTIONS.optional.BaselineTime       = 1 : size(OPTIONS.optional.Baseline,2);
    end
    
    % Baseline segment definition
    if isempty(OPTIONS.optional.BaselineSegment)
        OPTIONS.optional.BaselineSegment    = OPTIONS.optional.BaselineTime;
        STb     = 1;
        NDb     = size(OPTIONS.optional.Baseline,2);
    else
        STb     = be_closest( OPTIONS.optional.BaselineSegment(1), OPTIONS.optional.BaselineTime );
        NDb     = be_closest( OPTIONS.optional.BaselineSegment(end), OPTIONS.optional.BaselineTime );
    end
    
    % Check
    if STb > NDb || NDb > size(OPTIONS.optional.Baseline,2)
        error('In be_main : bad time definition for the baseline')
    end
    
    % Baseline Channels
    CH = OPTIONS.automatic.GoodChannel;
    if isfield(OPTIONS.optional, 'BaselineChannels') && ~isempty(OPTIONS.optional.BaselineChannels)
        CH = [];
        for ii = 1 : numel( OPTIONS.automatic.GoodChannel)
            iD  = strcmp( {OPTIONS.optional.BaselineChannels.Channel.Name},OPTIONS.optional.Channel(OPTIONS.automatic.GoodChannel(ii)).Name );
            if sum(iD); CH(ii) = find(iD); end
        end
    end
    
    if numel(CH) ~= numel(OPTIONS.automatic.GoodChannel)
        error('MEM : Channels of baseline file do not match channels in data file');
    end
    nBaseline =   zeros( size(OPTIONS.mandatory.Data,1), NDb-STb+1 );
    nBaseline(OPTIONS.automatic.GoodChannel,:) = OPTIONS.optional.Baseline(CH,STb:NDb);
    OPTIONS.optional.Baseline = nBaseline;
    
end
    
return

function OPT                = PIPELINEoptions(OPT, pipeline)
%% ===== method-specific I/O arguments % ====== %%

if nargin==1
    pipeline = OPT.mandatory.pipeline;
end

switch pipeline
    
    case 'cMEM'
        % clustering
        DEF.clustering.clusters_type    = 'blockwise';
        DEF.clustering.MSP_window       = [];
        DEF.clustering.MSP_scores_threshold = 0.7;
        
        % optional
        DEF.optional.Comment            = 'cMEM : MEG';
        
        % automatic 
        DEF.automatic.selected_samples  = [];
        
    case 'wMEM'
        % clustering
        DEF.clustering.clusters_type        = 'wfdr';
        DEF.clustering.MSP_scores_threshold = 'fdr';
        
        % wavelet processing
        DEF.wavelet.type                = 'RDW';
        DEF.wavelet.vanish_moments      = 4;
        DEF.wavelet.shrinkage           = 1;
        DEF.wavelet.selected_scales     = 0;
        DEF.wavelet.verbose             = 0;
        DEF.wavelet.single_box          = 0;
        
        % automatic
        DEF.automatic.selected_samples  = [];
        DEF.automatic.selected_jk       = [];
        DEF.automatic.selected_values   = [];
        DEF.automatic.Mod_in_boxes      = [];
        DEF.automatic.scales            = [];
        
        % optional
        DEF.optional.Comment            = 'wMEM : MEG';
        
    case 'rMEM'
        % clustering
        DEF.clustering.clusters_type    = 'wfdr';
        DEF.clustering.MSP_scores_threshold = 0;
                
        % wavelet processing
        DEF.wavelet.type                = 'CWT';
        DEF.wavelet.vanish_moments      = 4;
        DEF.wavelet.order               = 10;
        DEF.wavelet.nb_levels           = 128;
        DEF.wavelet.verbose             = 0;
        
        % ridge processing
        DEF.ridges.scalo_threshold      = .95;      
        DEF.ridges.energy_threshold     = .95; 
        DEF.ridges.strength_threshold   = [];
        DEF.ridges.skim_map             = 1;        
        DEF.ridges.frequency_range      = [40 100];
        DEF.ridges.min_duration         = 8; %(ms)
        DEF.ridges.cycles_in_window     = 2;
        DEF.ridges.fdrMOD               = {};
        
        % automatic
        DEF.automatic.selected_samples  = [];
        DEF.automatic.selected_jk       = [];    
        DEF.automatic.rMEMfiles         = {};
        DEF.automatic.Mod_in_ridges     = {};
        
        % optional
        DEF.optional.Comment            = 'rMEM : MEG';
        
end
% ========================================== %%

OPT = be_struct_copy_fields(OPT, DEF, []);
OPT.mandatory.pipeline  =   pipeline;

return

function [HeadModel, OPTIONS, FLAG]    = checkIO( HeadModel, OPTIONS, verbose )
    
    FLAG = 0;
    
    % Mandatory
    MFs     =   {'DataTypes', 'ChannelTypes', 'pipeline', 'stand_alone', 'DataTime'};
    if ~any( isfield(OPTIONS.mandatory, MFs) )
        idF     =   find( ~any( isfield(OPTIONS.mandatory, MFs) ) );
        fprintf('In be_main :\tmandatory field(s) %s undefined.\n', MFs{idF} );
        FLAG = 1;
    end
    if size( OPTIONS.mandatory.ChannelTypes, 1) == 1
        OPTIONS.mandatory.ChannelTypes = OPTIONS.mandatory.ChannelTypes';
    end

    % Data, baseline
    if ~isempty(OPTIONS.mandatory.Data)
        % Check data/time match
        if size(OPTIONS.mandatory.DataTime,2) ~= size(OPTIONS.mandatory.Data,2)
            fprintf('In be_main :\tData does not match time definition.\n');
            FLAG = 1;
        end

        % GoodChannel
        if OPTIONS.automatic.stand_alone
            % Check the consistency of # of channels
            GAIN    =   zeros( numel(OPTIONS.mandatory.ChannelTypes), size(HeadModel.Gain(1).matrix,2) );
            nC1     =   numel(OPTIONS.mandatory.ChannelTypes);
            nC2     =   size(OPTIONS.mandatory.Data,1);
            if (nC1~=nC2)
                fprintf('In be_main :\tChannels definition does not match data size.\n');
                FLAG = 1;
            end
            
            % Check the channels of selected modalities
            OPTIONS.automatic.GoodChannel = zeros( nC1,1 );
            for ii = 1 : numel(OPTIONS.mandatory.DataTypes)
                
                % Check if data present for modality ii
                iDc     =   strcmpi(OPTIONS.mandatory.ChannelTypes, OPTIONS.mandatory.DataTypes{ii});
                if ~sum(iDc)
                    fprintf('In be_main :\tNo channels found for data of type %s.\n', OPTIONS.mandatory.DataTypes{ii});
                    FLAG = 1;
                end
                OPTIONS.automatic.GoodChannel = OPTIONS.automatic.GoodChannel + iDc;
                
                % Check if gain present for modality ii
                MODii   =   strcmp(OPTIONS.mandatory.DataTypes{ii}, {HeadModel.Gain.modality});
                if ~any(MODii)
                    fprintf('In be_main :\tNo gain found for data of type %s.\n', OPTIONS.mandatory.DataTypes{ii});
                    FLAG = 1;
                end
                if size(HeadModel.Gain(MODii).matrix,1) ~= sum(iDc)
                    fprintf('In be_main :\tNo match between data and gain channel for modality : %s.\n', OPTIONS.mandatory.DataTypes{ii});
                    FLAG = 1;
                end
                GAIN(iDc,:)   =  HeadModel.Gain(ii).matrix;
                
            end
            OPTIONS.automatic.GoodChannel = logical(OPTIONS.automatic.GoodChannel);           
            OPTIONS.ChannelFlag     = ones( nC1,1 );
            HeadModel.Gain          =   GAIN;
            
        end
        
    elseif isempty(OPTIONS.optional.Baseline) && verbose
        fprintf('In be_main :\tNo data nor baseline were found. Aborted.\n');
        FLAG    =   1;
        
    end

    % Check manual clustering
    if ~isempty(OPTIONS.optional.clustering.clusters) 
        [FLAG, OPTIONS] = be_check_clustering(OPTIONS, HeadModel);    
    end
            
return
