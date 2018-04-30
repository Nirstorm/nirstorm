function [OPTIONS, OutputFiles] = be_ridgefilter(OPTIONS, varargin)
%RIDGEFILTER Returns the ridge-filetered data within time window and
% frequency band specified in options. The ridge strength threshold is obtained 
% from the baseline.
%
%   BE_RIDGEFILTER(OPTIONS) returns new signal
%   that contains a ridge-filtered portion of the original signal. 
%
%   INPUTS:
%       -   OPTIONS :   options structure
%
%   OUTPUTS:
%       -   OPTIONS 
%       -   OutputFiles:names of the new ridge-signal files
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


%% ===== SET RIDGE-FILTERING PRE-OPTIONS ===== %% 
[isstandalone, isprocess]    =   be_check_caller;
HM          =   [];
if isprocess || isstandalone
    HM     =    set_ridgefilter_options( OPTIONS );
end
OPTIONS    =    be_main_call( HM, OPTIONS );
OPTIONS.ridges.min_duration         =   fix( OPTIONS.ridges.min_duration* ...
    OPTIONS.automatic.sampling_rate/1000 );
OPTIONS.automatic.threshold_only    =   isempty(OPTIONS.mandatory.Data) && ...
    ~isempty(OPTIONS.optional.Baseline); 
% -------------------------------------------------------------------------



%% ====== Structure data ====== %%
OutputFiles = {}; sMat = {};

% Time limits
if diff(OPTIONS.optional.TimeSegment([1 end]))*OPTIONS.automatic.sampling_rate+1 < OPTIONS.ridges.min_duration
    fprintf('%s, In be_ridgefilter :\tRidges duration criterion is superior to data length. No ridges found\n',OPTIONS.mandatory.pipeline);
    return
end
% -------------------------------------------------------------------------



%% ===== COMPUTE RIDGES ===== %%

fprintf('\n%s, Ridge-based filtering process starts.\n',OPTIONS.mandatory.pipeline);

% Check if TF already computed
[d, WData, SCALO, F]  = look_for_TF(OPTIONS);
if isstruct(d)
    compute_TF = 0;
    OPTIONS.automatic.TFfile    =   F;
    OPTIONS.automatic.TFtime    =   load( F, 'Time' );
    OPTIONS.automatic.TFtime    =   OPTIONS.automatic.TFtime.Time;
    OPTIONS                     =   be_struct_copy_fields( OPTIONS, d, {'wavelet','ridges'}, 0 );
else
    compute_TF  =   1;
end

% If independent baseline 
thresh  =   zeros( 1, numel(OPTIONS.automatic.Modality) );
if strcmpi(OPTIONS.automatic.BaselineType, 'independent') && isempty(OPTIONS.ridges.strength_threshold)
    for ii = 1 : numel(OPTIONS.automatic.Modality)
        fprintf('%s, Processing %s baseline ...\n', OPTIONS.mandatory.pipeline, OPTIONS.mandatory.DataTypes{ii})
        [d, SCALO{ii}]  =   be_computeridges( ii, 'baseline', OPTIONS );
        % Store new wavelet field
        if ii == 1; OPTIONS =   be_struct_copy_fields(OPTIONS, d, {'wavelet'}); end
        thresh(ii)      =   be_rdg_threshold( SCALO{ii}, OPTIONS, ii );
        fprintf('\tRidge-strength threshold: %.2f', thresh(ii) );
    end      
end 

% In case only threshold is relevant
OPTIONS.automatic.nb_channels       =   numel(OPTIONS.mandatory.ChannelTypes);
if OPTIONS.automatic.threshold_only
    OPTIONS.ridges.strength_threshold   =   thresh;
    return;
end

% Compute ridges
if isstandalone && ~isprocess
    OPTIONS     =   be_main_channel([], [], OPTIONS);
elseif isprocess
    OPTIONS     =   be_main_channel(1, [], OPTIONS);
end
nScalo      =   zeros( OPTIONS.wavelet.nb_levels+1, size(OPTIONS.automatic.Modality(1).data,2) );
for ii  = 1 : numel(OPTIONS.mandatory.DataTypes)
    fprintf('%s, Processing %s data ...\n',OPTIONS.mandatory.pipeline, OPTIONS.mandatory.DataTypes{ii});
    if compute_TF
        [OPTIONS, SCALO{ii}, WData{ii}]   =   be_computeridges( ii, 'data', OPTIONS );
    else
        % Group ridges
        OPTIONS     =   be_ridge_lines( be_avgcell(SCALO{ii}), ii, OPTIONS);
    end

    % Get ridge strength trheshold
    OPTIONS.automatic.Modality(ii).ridgeMap   =   be_avgcell( SCALO{ii} );
    if isnan( OPTIONS.ridges.strength_threshold(ii) ) 
        fprintf('\tBaseline within data,[%3.2f,%3.2f] seconds\n', OPTIONS.optional.BaselineSegment(1), OPTIONS.optional.BaselineSegment(2) ); 
        OPTIONS = be_rdg_threshold( OPTIONS.automatic.Modality(ii).ridgeMap, OPTIONS, ii );        
    end 
    fprintf('\tRidge-strength threshold: %f\n', OPTIONS.ridges.strength_threshold(ii) );
    
    % Make ridges for modality
    rdg                 =   OPTIONS.automatic.Modality(ii).ridges_data;
    nScalo( [rdg{:}] )  =   1;
    fprintf('\t%i ridge lines found.\n', numel(rdg) )
    
end   
% -------------------------------------------------------------------------



%% ===== SAVE VARIABLES TO DATABASE ===== %%

if ~sum( nScalo(:) )
    fprintf('%s, Ridge-based filtering process : ABORTED\n',OPTIONS.mandatory.pipeline)
    return
end

% Save wavelet transform
if compute_TF
    SCALO   =   cellfun(@(a) cellfun(@(b) sparse(b),a,'uni', 0), SCALO, 'uni', 0);
    [OPTIONS.automatic.TFfile] = keepCWT(WData, SCALO, OPTIONS); 
end

% Get multimodal ridges lines
%[R, C] =   be_make_groups2( nScalo, OPTIONS, 1 );
[R, C] =   be_group_ridges(nScalo, OPTIONS.ridges.min_duration, 1);
fprintf( '%s, In RidgeFilter: %i ridge line(s) extracted.\n', OPTIONS.mandatory.pipeline, numel(R));

% Save ridge-filtered signals
for ii = 1:numel(R)
    [OutputFiles{ii}, sMat{ii}] = keepRDG(R, C, ii, WData, OPTIONS);
end
fprintf('%s, In RidgeFilter: all ridge-based signals stored.\n',OPTIONS.mandatory.pipeline);

% Add save signals to DB
if ~OPTIONS.automatic.stand_alone
    for ii = 1:numel(OutputFiles)
        db_add_data(OPTIONS.automatic.iStudy, OutputFiles{ii}, sMat{ii});
    end
    % Refresh display
%    db_reload_studies(OPTIONS.automatic.iStudy, 1);
    fprintf('%s, In RidgeFilter: new signals added to database.\n',OPTIONS.mandatory.pipeline);
end

% -------------------------------------------------------------------------
fprintf('%s, Ridge-based filtering process : COMPLETE \n',OPTIONS.mandatory.pipeline);

return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------- HELPER FUNCTIONS --------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [FileName]         =   keepCWT(WData, Scalo, OPTIONS)
% ===== KEEP CONTINUOUS WAVELET TRANSFORM =====

    fprintf('\nSaving continous wavelet transform ... ')
    % Create file structure
    if ~OPTIONS.automatic.stand_alone
        FileMat      	= db_template('timefreqmat');
    end 
    FileMat.Comment     = OPTIONS.automatic.TFcomment;
    FileMat.DataType    = OPTIONS.optional.FileType;
    FileMat.DataTypes   = OPTIONS.mandatory.DataTypes;
    FileMat.TF          = WData;
    FileMat.Time        = OPTIONS.automatic.TFtime;
    FileMat.TimeBands   = [];
    FileMat.Freqs       = OPTIONS.wavelet.freqs_analyzed;
    FileMat.RowNames    = OPTIONS.optional.ChannelNames;
    FileMat.Measure     = 'Power';
    FileMat.ChannelFlag = OPTIONS.optional.ChannelFlags;
    FileMat.nAvg        = 1;
    FileMat.ridgeMAP    = Scalo;
    
    % Options
    FileMat.Options.Measure         = 'Power';
    FileMat.Options.WaveletType     = 'Morse';
    FileMat.Options.VanishMoments   = OPTIONS.wavelet.vanish_moments;
    FileMat.Options.Order           = OPTIONS.wavelet.order;
    FileMat.Options.nLevels         = OPTIONS.wavelet.nb_levels;
    FileMat.Options.MaximaThresh    = OPTIONS.ridges.scalo_threshold;
    FileMat.Options.EnergyThresh    = OPTIONS.ridges.energy_threshold;
    FileMat.Options.RidgesThresh    = OPTIONS.ridges.strength_threshold;
    FileMat.Options.RidgesLength    = OPTIONS.ridges.min_duration;  
    FileMat.Options.RidgesFRange    = OPTIONS.ridges.frequency_range;
    FileMat.Options.RidgesWindow    = OPTIONS.optional.TimeSegment;
    FileMat.Options.RidgesBaseline  = OPTIONS.optional.BaselineSegment;
    
    % History: Computation
    if ~isfield(FileMat, 'History') | isempty(FileMat.History)
		FileMat.History = {datestr(date), 'compute', 'Time-frequency decomposition'};
    else
       FileMat.History(size(FileMat.History,1)+1,1:3) = {datestr(date), 'compute', 'Time-frequency decomposition'}; 
    end    

    if OPTIONS.automatic.stand_alone && ~OPTIONS.automatic.process
        FileName = FileMat;
    else
        iS  =   bst_get('Study', OPTIONS.automatic.iStudy);
        FileMat.DataFile  = file_win2unix(strrep(OPTIONS.optional.DataFile, OPTIONS.automatic.iProtocol.STUDIES, ''));
        % Output filename
        FileName = bst_process('GetNewFilename', bst_fileparts(iS.FileName), 'timefreq');
        % Save file
        save(FileName, '-struct', 'FileMat');
        % Add file to database structure
        db_add_data(OPTIONS.automatic.iStudy, FileName, FileMat);
    end   
    fprintf('done.\n');
    
return

function [OutputF, sMat]    =   keepRDG(R, C, ii, WData, OPTIONS)
    
    fprintf('\tStoring ridge #%i ...\n', ii)
        
    % Fuse signals
    idCH        =   1 : size( OPTIONS.mandatory.Data,1);
    if ~OPTIONS.automatic.stand_alone || OPTIONS.automatic.process
        FUS     =   be_fusion_of_modalities( [], [], OPTIONS );
        OPTIONS.automatic.nb_channels
        % Get study channel info
        nChan   =   bst_get('ChannelForStudy', OPTIONS.automatic.iStudy);        
        idCH    =   OPTIONS.automatic.GoodChannel;
    else
        Blims   =   be_closest( OPTIONS.optional.BaselineSegment, OPTIONS.mandatory.DataTime ); 
        FUS.baseline    =  OPTIONS.mandatory.Data(idCH, Blims(1):Blims(2) ); 
    end
    % Adjust truncated TFtime
    dt          =   ( be_closest( OPTIONS.automatic.TFtime(1),OPTIONS.mandatory.DataTime ) - 1 ) * (OPTIONS.wavelet.nb_levels+1);
    
    % Create a new signal from the extracted ridge lines
    sMat.F                  = NaN( nChan.nbChannels, numel(R{ii}) );
    sMat.Fi                 = NaN( nChan.nbChannels, numel(R{ii}) );
    sMat.Baseline           = NaN( nChan.nbChannels, size(FUS.baseline,2) );
    cpxCoeffs               = be_get_coeffs( WData,  OPTIONS.wavelet.freqs_analyzed, R{ii}-dt );
    sMat.F(idCH,:)          = real(cpxCoeffs);
    sMat.Fi(idCH,:)         = imag(cpxCoeffs);
    sMat.Baseline(idCH,:)   = FUS.baseline;
    sMat.BaselineSegment    = OPTIONS.optional.BaselineSegment;
    sMat.Comment            = ['Stand-alone call | ridgefiltered - (# ' num2str(ii) ')'];
    sMat.CentralFrq         = round( OPTIONS.wavelet.freqs_analyzed( round(C{ii}) ) );
    sMat.TFpath             = R{ii};
    sMat.rdgMOD             = {OPTIONS.automatic.Modality.ridgeMap};
    sMat.Time               = OPTIONS.mandatory.DataTime( fix(R{ii}/(OPTIONS.wavelet.nb_levels+1)) + 1 );
    sMat.ridgeID            = ii;
    sMat.TFfile             = OPTIONS.automatic.TFfile;
       
    
    if ~OPTIONS.automatic.stand_alone || OPTIONS.automatic.process
        
        % History
        TM                      = fix(clock);
        sMat.History            = OPTIONS.automatic.DataInfo.History;
        sMat.History{end+1,1}   = [date ' ' num2str(TM(4)) ':' num2str(TM(4)) ':' num2str(TM(4))];
        sMat.History{end,2}     = 'Ridge-based filtering';

        [iS,idS]            = bst_get('Study', OPTIONS.automatic.iStudy);
        sMat.Comment        = [OPTIONS.automatic.DataInfo.Comment ' | ridgefiltered (#' num2str(sMat.ridgeID) ')'];
        sMat.DataType       = OPTIONS.automatic.DataInfo.DataType;
        sMat.Device         = OPTIONS.automatic.DataInfo.Device;
        sMat.nAvg           = OPTIONS.automatic.DataInfo.nAvg;
        if isfield(OPTIONS.automatic.DataInfo, 'BadTrial')
            sMat.BadTrial       = OPTIONS.automatic.DataInfo.BadTrial;
        end
        sMat.ChannelFlag    = OPTIONS.automatic.DataInfo.ChannelFlag;
        sMat.History{end,3} = iS.Data(OPTIONS.automatic.iItem).FileName;
        sMat.TFfile         = file_win2unix(strrep(OPTIONS.automatic.TFfile, OPTIONS.automatic.iProtocol.STUDIES, ''));
        
        % save new file
        [d1, d2, d3]        = bst_fileparts( iS.Data(OPTIONS.automatic.iItem).FileName );
        OutputF             = be_fullfile(d1, [d2, '_RidgeFilter_', num2str(ii), d3]);
        
        % Save file
        FileName            = be_fullfile(OPTIONS.automatic.iProtocol.STUDIES, OutputF);
        bst_save(FileName, sMat, 'v6');
        
        % Add file to database structure
        nData               =   struct('FileName', OutputF, 'Comment', sMat.Comment, 'DataType', sMat.DataType, 'BadTrial', 0); 
        iData               =   length(iS.Data) + 1;
        iS.Data(iData)      =   nData;
        bst_set('Study', idS, iS);    
        panel_protocols('UpdateNode', 'Study', idS, 0);
        db_save();
       
    else
        sMat.TFfile.TF  =   [];       
        OutputF         =   sMat;
    end
    
return
           
function [OPTIONS]          =   get_selected_channels(OPTIONS)


    %% ====== Retrieve channels definition ============================= %%
    CH  = bst_get('ChannelForStudy', OPTIONS.automatic.iStudy);
    if ~(CH.nbChannels)
        error('Ridge-filter : No channel definition available.');
    end
    ChannelMat  = load( be_fullfile(OPTIONS.automatic.iProtocol.STUDIES, CH.FileName) , 'Channel');
    
    
    %% ====== Sort modalities according to BS ========================== %%
    TypeList    = {ChannelMat.Channel.Type};
    [dum, pos]  = unique(TypeList);
    TypeList    = TypeList( sort(pos) );
    [dum, pos]  = ismember( OPTIONS.mandatory.DataTypes, TypeList );
    [dum, pos]  = sort(pos);
    OPTIONS.mandatory.DataTypes = OPTIONS.mandatory.DataTypes(pos);
    
    
    %% ====== Get available channels =================================== %%
    iCHAN    = repmat( {'no good'}, 1, numel(ChannelMat.Channel) );
    if ~isempty(OPTIONS.mandatory.DataTypes)
        % Get sensor names
        for i = 1:numel(OPTIONS.mandatory.DataTypes)
            iChan   =   find(strcmpi({ChannelMat.Channel.Type}, OPTIONS.mandatory.DataTypes{i}));
            if isempty(iChan)
                warning(['Ridge-filter: No sensors are available for ' OPTIONS.mandatory.DataTypes{i}]) ;
            else
                iCHAN(iChan) = repmat( OPTIONS.mandatory.DataTypes(i), 1, numel(iChan) );
            end
        end
        % No sensors: error
        if isempty(iCHAN)
            error('No sensors are selected.')
        end
    end
    OPTIONS.mandatory.ChannelTypes = iCHAN;
    
    
return

function [d, W, S, F]       =   look_for_TF(OPTIONS)

W = {};
S = {};
d = [];
F = '';

if ~OPTIONS.automatic.stand_alone || OPTIONS.automatic.process
    
    fprintf('Looking for precomputed TF planes :\n')
    global GlobalData
    iP      =   bst_get('ProtocolInfo');
    
    % get T-F files datafiles
    DTfs    =   {GlobalData.DataBase.ProtocolStudies(GlobalData.DataBase.iProtocol).Study(OPTIONS.automatic.iStudy).Timefreq.FileName};
    
    % check if they match defined ridge properties
    TF      =   [];
    ii      =   1;
    stop    =   0;
    while ii <= numel(DTfs) && ~stop
        tf      =   load( be_fullfile(iP.STUDIES, DTfs{ii}), 'Options', 'DataTypes', 'Time', 'DataFile' );
        sameS   =   strcmpi( tf.DataFile, GlobalData.DataBase.ProtocolStudies(GlobalData.DataBase.iProtocol).Study(OPTIONS.automatic.iStudy).Data(OPTIONS.automatic.iItem).FileName );          C       =   isfield(tf.Options, 'WaveletType') && strcmpi(tf.Options.WaveletType, 'morse');
        if C && isempty(F) && sameS
            C(2)    =   tf.Options.nLevels >= OPTIONS.wavelet.nb_levels;
            C(3)    =   tf.Options.RidgesFRange(1) <= OPTIONS.ridges.frequency_range(1);
            C(4)    =   tf.Options.RidgesFRange(2) >= OPTIONS.ridges.frequency_range(2);
            C(5)    =   tf.Time(1)      <= OPTIONS.optional.TimeSegment(1);
            C(6)    =   tf.Time(end)    >= OPTIONS.optional.TimeSegment(end);
            
            RNG     =   [];
            for jj = 1 : numel(OPTIONS.mandatory.DataTypes)
                RNG = [RNG find( strcmp( OPTIONS.mandatory.DataTypes{jj}, tf.DataTypes ) )];
            end
            
            C(7)    =   numel(OPTIONS.mandatory.DataTypes)==numel(RNG); 
        end
            
        if sum(C)==7
            TF  =   load( be_fullfile(iP.STUDIES, DTfs{ii}), 'TF', 'ridgeMAP', 'DataFile', 'Freqs' );
            W   =   TF.TF(RNG);
            S   =   TF.ridgeMAP(RNG);
            
            % Make options structure
            d   =   be_main;
            d.wavelet.freqs_analyzed    =   TF.Freqs; 
            d.wavelet.vanish_moments    =   tf.Options.VanishMoments;
            d.wavelet.order             =   tf.Options.Order;
            d.wavelet.nb_levels         =   tf.Options.nLevels;
            d.ridges.scalo_threshold	=   tf.Options.MaximaThresh;
            d.ridges.energy_threshold	=   tf.Options.EnergyThresh;
            d.ridges.strength_threshold	=   tf.Options.RidgesThresh;
            d.ridges.frequency_range	=   tf.Options.RidgesFRange;
            d.ridges.min_duration       =   tf.Options.RidgesLength;
                
            F   =   be_fullfile( iP.STUDIES, DTfs{ii} );
            fprintf('\tTF plane found. Using : '' %s ''\n', F);
            stop = 1;
        else
            fprintf('\tTF plane %d : no match\n', ii)
            
        end
        ii = ii + 1;
    end
    
    if isempty(TF)
        fprintf('\tFound none. Recomputing ...\n')
    end
 
else
    fprintf('Stand-alone call: computing new TF plane\n')
    
end

return
    
function [HM, OPTIONS]      =   set_ridgefilter_options( OPTIONS )
    
    % Fake headmodel
    nMod    =   numel(OPTIONS.mandatory.DataTypes);
    for ii  =   1 : nMod
        HM.Gain(ii).modality    =   OPTIONS.mandatory.DataTypes{1};
        HM.Gain(ii).matrix      =   zeros( numel(find(strcmp(OPTIONS.mandatory.ChannelTypes,HM.Gain(ii).modality))), 1 );
    end
    HM.vertex_connectivity      =   0;
    
return