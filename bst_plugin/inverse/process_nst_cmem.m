function varargout = process_nst_cmem( varargin )

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Edouard Delaire (2022)
% Thomas Vincent, ZhengChen Cai (2015-2016)

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>

    % Description the process
    sProcess.Comment     = 'Compute sources: BEst';
    sProcess.FileTag     = '';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = {'NIRS', 'Sources'};
    sProcess.isSeparator = 0;
    sProcess.Index       = 1502;
    sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/NIRSTORM#Inverse_problem_using_cMEM';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw', 'data'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'data', 'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
       % Definition of the options
    sProcess.options.thresh_dis2cortex.Comment = 'Reconstruction Field of view (distance to montage border)';
    sProcess.options.thresh_dis2cortex.Type    = 'value';
    sProcess.options.thresh_dis2cortex.Value   = {3, 'cm', 2};


    sProcess.options.mem.Comment = {'panel_brainentropy', 'Source estimation options: '};
    sProcess.options.mem.Type    = 'editpref';
    sProcess.options.mem.Value   = be_main;
    

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

OutputFiles = {};
MethodOptions.MEMpaneloptions = sProcess.options.mem.Value.MEMpaneloptions;

% Canceled by user
if isempty(MethodOptions)
    return
end

% Install/load brainentropy plugin
[isInstalled, errMessage] = bst_plugin('Load', 'brainentropy');
if ~isInstalled
    bst_error('The Brainentropy toolbox is required to use MEM');
    return;
end

% Get plugin information
PluginDescription  = bst_plugin('GetInstalled', 'brainentropy'); 
if isempty(PluginDescription.GetVersionFcn) || bst_plugin('CompareVersions', PluginDescription.GetVersionFcn(), '3.0.0') < 0
   bst_error('Please update the BrainEntropy toolbox to the verson 3.0.0 or higher');
   return;
end

%% Backward compatibility
if isfield(sProcess.options, 'depth_weightingMNE') && isfield(sProcess.options, 'depth_weightingMEM')
    bst_report('Warning', sProcess, sInputs, 'Options for depth-weighting was moved to MEM panel. Please update your script')
    sProcess.options.mem.Value.MEMpaneloptions.model.depth_weigth_MNE      = sProcess.options.depth_weightingMNE.Value{1};
    sProcess.options.mem.Value.MEMpaneloptions.model.depth_weigth_MEM      = sProcess.options.depth_weightingMEM.Value{1};
end

%% Load head model
sStudy = bst_get('Study', sInputs.iStudy);
if isempty(sStudy.iHeadModel)
    bst_error('No head model found. Consider running "NIRS -> Compute head model"');
    return;
end
HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;


%% Load data & baseline
% Load recordings
if strcmp(sInputs.FileType, 'data')     % Imported data structure
    sDataIn = in_bst_data(sInputs(1).FileName);
elseif strcmp(sInputs.FileType, 'raw')  % Continuous data file
    sDataIn = in_bst(sInputs(1).FileName, [], 1, 1, 'no');
end

ChannelMat = in_bst_channel(sInputs(1).ChannelFile);
if ~isfield(ChannelMat.Nirs, 'Wavelengths')
    bst_error(['cMEM source reconstruction works only for dOD data ' ... 
               ' (eg do not use MBLL prior to this process)']);
    return;
end

OPTIONS         = getOptions(sProcess, HeadModelFile, sInputs(1).FileName);
pipeline        = OPTIONS.MEMpaneloptions.mandatory.pipeline;

if  strcmp(pipeline,'rMEM')
    bst_report('Warning', sProcess, sInputs, 'rMEM was not tested for fNIRS data, proceed with caution')
end

%% Run MEM
bst_progress('start', ['Reconstruction by ' pipeline], sprintf('Launching %s...', pipeline));
sResults = Compute(OPTIONS,ChannelMat, sDataIn );

%% Save results
bst_progress('text', 'Saving Results...');

for iMap = 1:length(sResults)

    ResultFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName),  ['results_NIRS_' nst_protect_fn_str(sResults(iMap).Comment)]);

    ResultsMat = sResults(iMap);
    ResultsMat.DataFile   = sInputs.FileName;

    bst_save(ResultFile, ResultsMat, 'v6');
    db_add_data( sInputs.iStudy, ResultFile, ResultsMat);

    OutputFiles{end+1} = ResultFile;

end

bst_progress('stop', 'Reconstruction by MNE', 'Finishing...');
end

function sResults = Compute(OPTIONS, ChannelMat, sDataIn )


    nirs_head_model = in_bst_headmodel(OPTIONS.HeadModelFile, 1);
    if ~isfield(nirs_head_model, 'NIRSMethod') && ndims(nirs_head_model.Gain) == 3
        nirs_head_model = process_nst_import_head_model('convert_head_model', ChannelMat, nirs_head_model, 0);
    end

    sCortex = in_tess_bst(nirs_head_model.SurfaceFile);
    
    nb_nodes        = size(sCortex.Vertices, 1);
    nb_wavelengths  = length(ChannelMat.Nirs.Wavelengths);
    HM.SurfaceFile = nirs_head_model.SurfaceFile;

    %% define the reconstruction FOV
    thresh_dis2cortex       = OPTIONS.thresh_dis2cortex;
    valid_nodes             = nst_headmodel_get_FOV(ChannelMat, sCortex, thresh_dis2cortex, sDataIn.ChannelFlag );

    OPTIONS.MEMpaneloptions.optional.cortex_vertices = sCortex.Vertices(valid_nodes, :); 
    HM.vertex_connectivity = sCortex.VertConn(valid_nodes, valid_nodes);

    isReconstructed = true(1, nb_wavelengths); 
    for iwl=1:nb_wavelengths
        swl = ['WL' num2str(ChannelMat.Nirs.Wavelengths(iwl))];
        selected_chans = strcmpi({ChannelMat.Channel.Group}, swl) & (sDataIn.ChannelFlag>0)';
        
        if ~any(selected_chans)
            isReconstructed(iwl) = false;
            continue
        end

        OPTIONS.GoodChannel     = ones(sum(selected_chans), 1);
        OPTIONS.ChannelFlag     = ones(sum(selected_chans), 1);
        OPTIONS.Channel         = ChannelMat.Channel(selected_chans);
        OPTIONS.DataTime        = round(sDataIn.Time,6);
        OPTIONS.Data            = sDataIn.F(selected_chans,:);
    
        
        HM.Gain = nirs_head_model.Gain(selected_chans, valid_nodes); 
        
        %% launch MEM (cMEM only in current version)
        bst_progress('text', ['Running cMEM for wavelength #' num2str(iwl) '...']);
        [result, sOptions(iwl)] = be_main_call(HM, OPTIONS);

        if strcmp(OPTIONS.MEMpaneloptions.mandatory.pipeline ,'wMEM')
            selected_samples = sOptions(iwl).automatic.selected_samples;

            % sort the sample by time instead of energy
            [~,ia] = sort(selected_samples(1,:));
            result.ImageGridAmp{1} = result.ImageGridAmp{1}(:,ia);
            result.ImageGridAmp{2} = result.ImageGridAmp{2}(ia,:);

            sOptions(iwl).automatic.selected_samples = selected_samples(:,ia);
        end

        result.Time    = OPTIONS.DataTime;
        result.Options =  sOptions(iwl);
        result.Function = result.Options.FunctionName;
        result.Comment =  [sOptions(iwl).Comment ' | ' swl 'nm'];
        result.History  = OPTIONS.History;
        result.HeadModelFile = OPTIONS.HeadModelFile;
        result.HeadModelType = nirs_head_model.HeadModelType;
        result.SurfaceFile   = file_short(nirs_head_model.SurfaceFile);

        result = bst_history('add', result, 'compute', sOptions(iwl).Comment );
        result.DisplayUnits   =  'OD';
        
        sResults(iwl) = result;
    end
    
    % Filter reconstructed wavelengh
    sResults = sResults(isReconstructed);
    sOptions = sOptions(isReconstructed);

    % Compute MBLL 
    if length(sResults) > 1

        sResults_hb = nst_mbll_source(sResults, ChannelMat.Nirs.Wavelengths(isReconstructed));
        sResults    = [ sResults, sResults_hb];

    end

    isSaveFactor = isfield(sOptions(1), 'output') && sOptions(1).output.save_factor;
    sResults = nst_misc_FOV_to_cortex(sResults, nb_nodes, valid_nodes, isSaveFactor);
end

function OPTIONS = getOptions(sProcess, HeadModelFileName, DataFile)

    MethodOptions.MEMpaneloptions =   sProcess.options.mem.Value.MEMpaneloptions;
    MethodOptions.SourceOrient{1} = 'fixed';

    OPTIONS = process_inverse_2016('Compute');
    OPTIONS.InverseMethod = 'mem';
    OPTIONS = struct_copy_fields(OPTIONS, MethodOptions, 1);
    OPTIONS.DataTypes = {'NIRS'};
    OPTIONS.NoiseCov = [];
    OPTIONS.MEMpaneloptions.solver.NoiseCov_recompute   = 1;

    OPTIONS.thresh_dis2cortex = sProcess.options.thresh_dis2cortex.Value{1}.*0.01;
    
    OPTIONS.Comment         = 'MEM';
    OPTIONS.FunctionName    = 'mem';
    OPTIONS.DataFile        = DataFile;
    OPTIONS.ResultFile      = [];
    OPTIONS.HeadModelFile   =  HeadModelFileName;
    
    sDataIn = in_bst_data(DataFile, 'History');
    OPTIONS.History       = sDataIn.History;

end
