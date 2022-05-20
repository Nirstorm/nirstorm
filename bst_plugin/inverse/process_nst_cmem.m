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
    
    sProcess.options.depth_weightingMNE.Comment = 'Depth weighting factor for <B>MNE</B>';
    sProcess.options.depth_weightingMNE.Type    = 'value';
    sProcess.options.depth_weightingMNE.Value   = {0.5, '', 1};
    
    sProcess.options.depth_weightingMEM.Comment = 'Depth weighting factor for <B>MEM</B>';
    sProcess.options.depth_weightingMEM.Type    = 'value';
    sProcess.options.depth_weightingMEM.Value   = {0.3, '', 1};

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
   
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
OutputFiles = {};
MethodOptions.MEMpaneloptions = sProcess.options.mem.Value.MEMpaneloptions;
store_sparse_results          = sProcess.options.store_sparse_results.Value;

% Canceled by user
if isempty(MethodOptions)
    return
end

% Install/load brainentropy plugin
[isInstalled, errMessage] = bst_plugin('Install', 'brainentropy', 1);
if ~isInstalled
    return;
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
    sDataRaw = in_bst_data(sInputs(1).FileName, 'F');
end

ChannelMat = in_bst_channel(sInputs(1).ChannelFile);
if ~isfield(ChannelMat.Nirs, 'Wavelengths')
    bst_error(['cMEM source reconstruction works only for dOD data ' ... 
               ' (eg do not use MBLL prior to this process)']);
    return;
end

nb_wavelengths  = length(ChannelMat.Nirs.Wavelengths);
measure_tag     = 'WL';
OPTIONS = getOptions(sProcess,nirs_head_model, sInputs(1).FileName);

%% Run cMEM
bst_progress('start', 'Reconstruction by cMEM', 'Launching cMEM...');
[dOD_sources_cMEM,Hb_sources, O_updated] = Compute(OPTIONS,ChannelMat, sDataIn );

bst_progress('text', 'Saving Results...');

for iwl=1:nb_wavelengths
     swl = [measure_tag num2str(ChannelMat.Nirs.Wavelengths(iwl))];
     [sStudy, ResultFile] = add_surf_data(dOD_sources_cMEM(:,iwl,:), sDataIn.Time, nirs_head_model, ...
         [OPTIONS.Comment     ' sources - ' swl 'nm'], ...
         sInputs, sStudy, O_updated.Comment, ...
         'OD', store_sparse_results);
     OutputFiles{end+1} = ResultFile;
end

hb_unit_factor = 1e6;
hb_unit = '\mumol.l-1';
hb_types = {'HbO', 'HbR','HbT'};
for ihb=1:3
    [sStudy, ResultFile] = add_surf_data(squeeze(Hb_sources(:,ihb,:)) .* hb_unit_factor,...
                                         sDataIn.Time, nirs_head_model, ...
                                         [OPTIONS.Comment  ' sources - ' hb_types{ihb}], ...
                                         sInputs, sStudy, O_updated.Comment, ...
                                         hb_unit, store_sparse_results);    
    OutputFiles{end+1} = ResultFile;
end

bst_progress('stop', 'Reconstruction by cMEM', 'Finishing...');
% Update Brainstorm database
bst_set('Study', sInputs.iStudy, sStudy);
end

function [dOD_sources,Hb_sources, O_updated] = Compute(OPTIONS,ChannelMat, sDataIn )

    nirs_head_model = in_bst_headmodel(OPTIONS.HeadModelFile);
    cortex = in_tess_bst(nirs_head_model.SurfaceFile);
    
    nb_nodes = size(cortex.Vertices, 1);
    nb_samples = length(sDataIn.Time);
    nb_wavelengths  = length(ChannelMat.Nirs.Wavelengths);

    HM.SurfaceFile = nirs_head_model.SurfaceFile;

    %% define the reconstruction FOV
    thresh_dis2cortex       = OPTIONS.thresh_dis2cortex;
    valid_nodes             = nst_headmodel_get_FOV(ChannelMat, cortex, thresh_dis2cortex);

    OPTIONS.MEMpaneloptions.optional.cortex_vertices = cortex.Vertices(valid_nodes, :); 
    HM.vertex_connectivity = cortex.VertConn(valid_nodes, valid_nodes);

    %% estimate the neighborhood order for cMEM  (goal: # of clusters ~= # of good channels) 
    if OPTIONS.flag_auto_nbo
        swl = ['WL' num2str(ChannelMat.Nirs.Wavelengths(1))];
        n_channel = sum(strcmpi({ChannelMat.Channel.Group}, swl) & (sDataIn.ChannelFlag>0)');
    
        nbo = estimate_nbo(cortex, valid_nodes, n_channel, 1 );
        OPTIONS.MEMpaneloptions.clustering.neighborhood_order = nbo;
    end


    dOD_sources = zeros(nb_nodes, nb_wavelengths, nb_samples);

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
        [Results, O_updated] = be_main_call_NIRS(HM, OPTIONS);
        
        %cMEM results
        grid_amp = zeros(nb_nodes, nb_samples); 
        grid_amp(valid_nodes,:) = Results.ImageGridAmp;
        
        dOD_sources(:, iwl, :) = grid_amp; 
    end
    
    bst_progress('text', 'Calculating HbO/HbR/HbT in source space...');
    % Compute dHb
    hb_extinctions = nst_get_hb_extinctions(ChannelMat.Nirs.Wavelengths);
    hb_extinctions = hb_extinctions ./10;% mm-1.mole-1.L

    Hb_sources = zeros(nb_nodes, 3, nb_samples);
    for idx=1:length(valid_nodes)
        inode = valid_nodes(idx);
        Hb_sources(inode, 1:2, :) = pinv(hb_extinctions) * ...
                                    squeeze(dOD_sources(inode, :, :));
    
    end
    Hb_sources(:,3,:) = squeeze(sum(Hb_sources, 2));
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
    OPTIONS.MEMpaneloptions.model.depth_weigth_MNE      = sProcess.options.depth_weightingMNE.Value{1};
    OPTIONS.MEMpaneloptions.model.depth_weigth_MEM      = sProcess.options.depth_weightingMEM.Value{1};
    OPTIONS.MEMpaneloptions.model.NoiseCov_recompute    = sProcess.options.NoiseCov_recompute.Value;

    OPTIONS.thresh_dis2cortex = sProcess.options.thresh_dis2cortex.Value{1}.*0.01;
    OPTIONS.flag_auto_nbo = sProcess.options.auto_neighborhood_order.Value;

    %% Run cMEM
    
    OPTIONS.Comment         = 'cMEM';
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
                                              data_unit, store_sparse)
                                          
    if nargin < 8
        data_unit = '';
    end
    
    ResultFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ...
                            ['results_NIRS_' protect_fn_str(name)]);
                        
    % ===== CREATE FILE STRUCTURE =====
    ResultsMat = db_template('resultsmat');
    ResultsMat.Comment       = name;
    ResultsMat.Function      = '';
    if store_sparse
        ResultsMat.ImageGridAmp = sparse(data); %TODO TOCHECK with FT: sparse data seem not well handled. Eg while viewing (could not reproduce)
    else
        ResultsMat.ImageGridAmp = data;
    end
    ResultsMat.DisplayUnits = data_unit;
    ResultsMat.Time          = time;
    ResultsMat.DataFile      = sInputs.FileName;
    ResultsMat.HeadModelFile = head_model.FileName;
    ResultsMat.HeadModelType = head_model.HeadModelType;
    ResultsMat.ChannelFlag   = [];
    ResultsMat.GoodChannel   = [];
    ResultsMat.SurfaceFile   = file_short(head_model.SurfaceFile);
    ResultsMat.GridLoc    = [];
    ResultsMat.GridOrient = [];
    ResultsMat.nAvg      = 1;
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
