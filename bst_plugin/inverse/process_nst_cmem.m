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
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/NIRSFingerTapping';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw', 'data'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'data', 'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
     
    % Options: MNE options
    sProcess.options.mem.Comment = {'panel_brainentropy', 'Source estimation options: '};
    sProcess.options.mem.Type    = 'editpref';
    sProcess.options.mem.Value   = be_main;
    
    % Definition of the options
    sProcess.options.thresh_dis2cortex.Comment = 'Reconstruction Field of view (distance to montage border)';
    sProcess.options.thresh_dis2cortex.Type    = 'value';
    sProcess.options.thresh_dis2cortex.Value   = {3, 'cm',2};
    
    sProcess.options.NoiseCov_recompute.Comment = 'Compute noise convariance for MNE';
    sProcess.options.NoiseCov_recompute.Type    = 'checkbox';
    sProcess.options.NoiseCov_recompute.Value   = 1;
        

    %TODO: tell neighborhood order will be ignored after
    sProcess.options.auto_neighborhood_order.Comment = 'Set neighborhood order automatically (default)';
    sProcess.options.auto_neighborhood_order.Type    = 'checkbox';
    sProcess.options.auto_neighborhood_order.Value   = 1;
    
    sProcess.options.store_sparse_results.Comment = 'Store sparse results';
    sProcess.options.store_sparse_results.Type    = 'checkbox';
    sProcess.options.store_sparse_results.Value   = 0;
    sProcess.options.store_sparse_results.Group   = 'output';

   
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
[isInstalled, errMessage] = bst_plugin('Install', 'brainentropy', 1);
if ~isInstalled
    bst_error('The Brainentropy toolbox is required to use MEM');
    return;
end

% Get plugin information
PluginDescription  = bst_plugin('GetInstalled', 'brainentropy'); 
if isempty(PluginDescription.GetVersionFcn) || bst_plugin('CompareVersions', PluginDescription.GetVersionFcn(), '2.7.4') < 0
   bst_error('Please update the BrainEntropy toolbox to the verson 2.7.4 or higher');
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

nirs_head_model = in_bst_headmodel(sStudy.HeadModel(sStudy.iHeadModel).FileName);
nirs_head_model.FileName = sStudy.HeadModel(sStudy.iHeadModel).FileName;

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

OPTIONS         = getOptions(sProcess,nirs_head_model, sInputs(1).FileName);
pipeline        = OPTIONS.MEMpaneloptions.mandatory.pipeline;

if strcmp(pipeline,'wMEM') || strcmp(pipeline,'rMEM')
    bst_report('Warning', sProcess, sInputs, sprintf('%s was not tested for fNIRS data, proceed with caution',pipeline ))
end

%% Run MEM
bst_progress('start', ['Reconstruction by ' pipeline], sprintf('Launching %s...', pipeline));
sResults = Compute(OPTIONS,ChannelMat, sDataIn );

%% Save results

bst_progress('text', 'Saving Results...');

for iMap = 1:length(sResults)
     [sStudy, ResultFile] = add_surf_data(sResults(iMap).ImageGridAmp , sDataIn.Time, nirs_head_model, ...
                                          sResults(iMap).Comment, sInputs, sStudy, ...
                                          sResults(iMap).History, sResults(iMap).Units , ...
                                          sResults(iMap).MEMoptions);

    OutputFiles{end+1} = ResultFile;

end


bst_progress('stop', ['Reconstruction by ' pipeline], 'Finishing...');
% Update Brainstorm database
bst_set('Study', sInputs.iStudy, sStudy);
end

function sResults = Compute(OPTIONS,ChannelMat, sDataIn )


    nirs_head_model = in_bst_headmodel(OPTIONS.HeadModelFile);
    cortex = in_tess_bst(nirs_head_model.SurfaceFile);
    
    nb_nodes        = size(cortex.Vertices, 1);
    nb_wavelengths  = length(ChannelMat.Nirs.Wavelengths);

    HM.SurfaceFile = nirs_head_model.SurfaceFile;

    %% define the reconstruction FOV
    thresh_dis2cortex       = OPTIONS.thresh_dis2cortex;
    valid_nodes             = nst_headmodel_get_FOV(ChannelMat, cortex, thresh_dis2cortex,sDataIn.ChannelFlag );

    OPTIONS.MEMpaneloptions.optional.cortex_vertices = cortex.Vertices(valid_nodes, :); 
    HM.vertex_connectivity = cortex.VertConn(valid_nodes, valid_nodes);

    %% estimate the neighborhood order for cMEM  (goal: # of clusters ~= # of good channels) 
    if OPTIONS.flag_auto_nbo
        swl = ['WL' num2str(ChannelMat.Nirs.Wavelengths(1))];
        n_channel = sum(strcmpi({ChannelMat.Channel.Group}, swl) & (sDataIn.ChannelFlag>0)');
    
        nbo = estimate_nbo(cortex, valid_nodes, n_channel, 1 );
        OPTIONS.MEMpaneloptions.clustering.neighborhood_order = nbo;
    end

       
    for iwl=1:nb_wavelengths
        swl = ['WL' num2str(ChannelMat.Nirs.Wavelengths(iwl))];
        selected_chans = strcmpi({ChannelMat.Channel.Group}, swl) & (sDataIn.ChannelFlag>0)';
        
        OPTIONS.GoodChannel     = ones(sum(selected_chans), 1);
        OPTIONS.ChannelFlag     = ones(sum(selected_chans), 1);
        OPTIONS.Channel         = ChannelMat.Channel(selected_chans);
        OPTIONS.DataTime        = round(sDataIn.Time,6);
        OPTIONS.Data            = sDataIn.F(selected_chans,:);
    
        gain = nst_headmodel_get_gains(nirs_head_model,iwl, ChannelMat.Channel, find(selected_chans));
        
        % Remove 0 from the gain matrix
        HM.Gain = gain(:,valid_nodes);
        HM.Gain(HM.Gain==0) = min(HM.Gain(HM.Gain>0));
    
        bst_progress('text', ['WL' num2str(iwl) ', kept ' num2str(length(valid_nodes)) ...
                 ' nodes that were in VOI and have non-zero sensitivity']);
    
        %% launch MEM (cMEM only in current version)
        bst_progress('text', ['Running cMEM for wavelength #' num2str(iwl) '...']);
        [result, sOptions(iwl)] = be_main_call(HM, OPTIONS);

        if strcmp(OPTIONS.MEMpaneloptions.mandatory.pipeline ,'wMEM')
            selected_samples = result.MEMoptions.automatic.selected_samples;

            % sort the sample by time instead of energy
            [~,ia] = sort(selected_samples(1,:));
            result.ImageGridAmp{1} = result.ImageGridAmp{1}(:,ia);
            result.ImageGridAmp{2} = result.ImageGridAmp{2}(ia,:);

            result.MEMoptions.automatic.selected_samples = result.MEMoptions.automatic.selected_samples(:,ia);
        end

        result.Comment =  [result.MEMoptions.automatic.Comment ' | ' swl 'nm'];
        result.History =  [result.MEMoptions.automatic.Comment];
        result.Units   =  'OD';
        result.MEMoptions.automatic.neighborhood_order = sOptions(iwl).MEMpaneloptions.clustering.neighborhood_order;
        result.MEMoptions.automatic.valid_nodes = valid_nodes;
        
        sResults(iwl) = result;
    end
    
    bst_progress('text', 'Calculating HbO/HbR/HbT in source space...');
    % Compute dHb
    hb_extinctions = nst_get_hb_extinctions(ChannelMat.Nirs.Wavelengths);
    hb_extinctions = hb_extinctions ./10;% mm-1.mole-1.L

    if ~iscell(sResults(1).ImageGridAmp)
        dOD_sources =  permute( cat(3, sResults.ImageGridAmp), [1 3 2]);
    else
        isConsistent = 1;
        for iResult = 2:length(sResults)
            isConsistent =  isConsistent && isequal( sResults(1).ImageGridAmp{2}, sResults(iResult).ImageGridAmp{2});
        end

        if isConsistent
            dOD_sources = zeros(size(sResults(1).ImageGridAmp{1}, 1) ,length(sResults),size(sResults(1).ImageGridAmp{1}, 2));
            for iResult = 1:length(sResults)
                dOD_sources(:, iResult, :)  = sResults(iResult).ImageGridAmp{1};
            end
        else
            % If not consistent, go back to full time-course
            dOD_sources = zeros(size(sResults(1).ImageGridAmp{1}, 1) ,length(sResults),size(sResults(1).ImageGridAmp{2}, 2));
            for iResult = 1:length(sResults)
                sResults(iResult).ImageGridAmp = sResults(iResult).ImageGridAmp{1} * sResults(iResult).ImageGridAmp{2};
                dOD_sources(:, iResult, :)  = sResults(iResult).ImageGridAmp;
            end
        end
    end

    Hb_sources = zeros(length(valid_nodes), 3, size(dOD_sources,3));
    for inode=1:length(valid_nodes)
        Hb_sources(inode, 1:2, :) = pinv(hb_extinctions) * ...
                                    squeeze(dOD_sources(inode, :, :));
    
    end
    Hb_sources(:,3,:) = squeeze(sum(Hb_sources, 2));

    hb_unit_factor = 1e6;
    hb_unit = '\mumol.l-1';
    hb_types = {'HbO', 'HbR','HbT'};

    sResults_hb = repmat(sResults(1), 1, 3);
    for iHb = 1:3
        sResults_hb(iHb).Comment = [ sResults(end).History     ' | ' hb_types{iHb}];
        sResults_hb(iHb).History = sResults(end).History;
        sResults_hb(iHb).Units = hb_unit;
        if iscell(sResults_hb(iHb).ImageGridAmp )
            sResults_hb(iHb).ImageGridAmp{1} = squeeze(Hb_sources(:,iHb,:)) .* hb_unit_factor;
        else
            sResults_hb(iHb).ImageGridAmp = squeeze(Hb_sources(:,iHb,:)) .* hb_unit_factor;
        end
    end

    sResults = [ sResults, sResults_hb];

    mapping = zeros(nb_nodes, length(valid_nodes)); 
    for iNode = 1:length(valid_nodes)
        mapping(valid_nodes(iNode), iNode) = 1;
    end
    mapping = sparse(mapping);

    for iMap = 1:length(sResults)
        if iscell(sResults(iMap).ImageGridAmp)
            sResults(iMap).ImageGridAmp = [ {mapping} sResults(iMap).ImageGridAmp ];
        else
            sResults(iMap).ImageGridAmp  = {mapping ,  sResults(iMap).ImageGridAmp};
        end
    end
end


function OPTIONS = getOptions(sProcess,HeadModel, DataFile)
    MethodOptions.MEMpaneloptions =   sProcess.options.mem.Value.MEMpaneloptions;
    % Add fields that are not defined by the options of the MEM interface
    if ~isempty(MethodOptions)
        switch (HeadModel.HeadModelType)
            case {'surface', 'ImageGrid'}
                MethodOptions.SourceOrient{1} = 'fixed';
            case 'volume'
                MethodOptions.SourceOrient{1} = 'free';
                MethodOptions.flagSourceOrient = [0 0 2 0];
        end
    end
    % Canceled by user
    if isempty(MethodOptions)
        return
    end
    % Add options to list
    OPTIONS = process_inverse_2016('Compute');
    OPTIONS.InverseMethod = 'mem';
    OPTIONS = struct_copy_fields(OPTIONS, MethodOptions, 1);
    OPTIONS.DataTypes = {'NIRS'};
    OPTIONS.NoiseCov = [];
    OPTIONS.MEMpaneloptions.solver.NoiseCov_recompute   = 1;
    OPTIONS.MEMpaneloptions.model.NoiseCov_recompute    = sProcess.options.NoiseCov_recompute.Value;

    OPTIONS.thresh_dis2cortex = sProcess.options.thresh_dis2cortex.Value{1}.*0.01;
    OPTIONS.flag_auto_nbo = sProcess.options.auto_neighborhood_order.Value;

    %% Run cMEM
    
    OPTIONS.Comment         = 'MEM';
    OPTIONS.DataFile        = DataFile;
    OPTIONS.ResultFile      = [];
    OPTIONS.HeadModelFile   =  HeadModel.FileName;
    OPTIONS.FunctionName    = 'mem';

end


function nbo = estimate_nbo (cortex, valid_nodes, nClust,isRandom)

    if nargin < 4
        isRandom = 1; 
    end

    bst_progress('text', 'Setting neighborhood order for clustering automatically...');
    
    % see panel_scout.m line 2019
    % Split hemispheres

    [rH, lH] = tess_hemisplit(cortex);
    rH       = intersect(rH,valid_nodes);
    lH       = intersect(lH,valid_nodes);

    if ~isempty(rH) && ~isempty(lH)
        isConnected = 0;
    else
        isConnected = 1;
    end

    % Clustering
    if isConnected
        Labels = tess_cluster(cortex.VertConn(valid_nodes,valid_nodes), nClust, isRandom);
    else
        ratiolH = length(lH)./(length(lH)+length(rH));
        ratiorH = length(rH)./(length(lH)+length(rH));
        Labels = ones(length(valid_nodes), 1); 
        Labels(ismember(valid_nodes,lH)) = tess_cluster(cortex.VertConn(lH,lH), max(1,round(nClust.*ratiolH)), isRandom);
        Labels(ismember(valid_nodes,rH)) = tess_cluster(cortex.VertConn(rH,rH), max(1,round(nClust.*ratiorH)), isRandom) + max(Labels(ismember(valid_nodes,lH)));
    end    
    uniqueLabels = unique(Labels);
    
    % neighborhood order list
    VoisinsOA = adj2Voisins(cortex.VertConn(valid_nodes,valid_nodes));
    neighborhood_order = zeros(length(uniqueLabels),1);
    for iScout = 1:length(uniqueLabels)
        vertices = find(Labels == uniqueLabels(iScout));
        seed = vertices(1);
        [~,nbo_idx] = min(abs(VoisinsOA(:,seed) - length(vertices)));
        neighborhood_order(iScout) = nbo_idx;     
    end

    % MEM clustering is different to the one above +2 is an emperical value
    % best way should be tuning nbo in be_create_clusters.m but it takes a
    % lot of time to compute. 
    nbo = round(mean(neighborhood_order))+2; 
end

function [sStudy, ResultFile] = add_surf_data(data, time, head_model, name, ...
                                              sInputs, sStudy, history_comment, ...
                                              data_unit, diagnosis)
                                          
    if nargin < 8
        data_unit = '';
    end
    
    ResultFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ...
                            ['results_NIRS_' protect_fn_str(name)]);
                        
    % ===== CREATE FILE STRUCTURE =====
    ResultsMat = db_template('resultsmat');
    ResultsMat.Comment       = name;
    ResultsMat.Function      = '';
    ResultsMat.ImageGridAmp  = data;

    if nargin >= 9 && ~isempty(diagnosis) 
        ResultsMat.diagnosis = diagnosis;
    end

    ResultsMat.DisplayUnits  = data_unit;
    ResultsMat.Time          = time;
    ResultsMat.DataFile      = sInputs.FileName;
    ResultsMat.HeadModelFile = head_model.FileName;
    ResultsMat.HeadModelType = head_model.HeadModelType;
    ResultsMat.ChannelFlag   = [];
    ResultsMat.GoodChannel   = [];
    ResultsMat.SurfaceFile   = file_short(head_model.SurfaceFile);
    ResultsMat.GridLoc       = [];
    ResultsMat.GridOrient    = [];
    ResultsMat.nAvg          = 1;

    % History
    ResultsMat = bst_history('add', ResultsMat, 'compute', history_comment);
    % Save new file structure
    bst_save(ResultFile, ResultsMat, 'v6');
    % ===== REGISTER NEW FILE =====
    % Create new results structure
    newResult = db_template('results');
    newResult.Comment       = name;
    newResult.FileName      = file_short(ResultFile);
    newResult.DataFile      = sInputs.FileName;
    newResult.isLink        = 0;
    newResult.HeadModelType = ResultsMat.HeadModelType;
    % Add new entry to the database
    iResult = length(sStudy.Result) + 1;
    sStudy.Result(iResult) = newResult;
    % Update Brainstorm database
    bst_set('Study', sInputs.iStudy, sStudy);
end

function sfn = protect_fn_str(s)
sfn = strrep(s, ' ', '_');
end

function VoisinsOA = adj2Voisins(adj)
% Convert the adjacency matrix 'adj' to 'VoisinsOA' neighbor cell vector
len = length(adj);

VoisinsOA = zeros(12,len);

%h = waitbar(0,'Please wait...');
for iScale = 1:12
    adj_i = logical(adj^iScale);
    adj_i(logical(eye(size(adj_i)))) = 0;
    for iSource = 1:len
       VoisinsOA(iScale,iSource) = length(find(adj_i(iSource,:)));
    end
    %waitbar(iScale / 12)
end
%close(h) 

end
