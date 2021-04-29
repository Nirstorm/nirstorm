function varargout = process_nst_cmem_fusion( varargin )

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
% Authors: Thomas Vincent, ZhengChen Cai (2015-2016)
%
% TODO: fix baseline as external data input

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Compute sources: BEst'; % in panel_brainentropy.m comment name has to be "Compute sources: BEst" to call BEST panel
sProcess.FileTag     = '';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = {'NIRS', 'Work in progress'};
sProcess.Index       = 1801; %0: not shown, >0: defines place in the list of processes
sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/NIRSFingerTapping';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'raw', 'data'};
% Definition of the outputs of this process
sProcess.OutputTypes = {'data', 'data'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
sProcess.isSeparator = 1;

% Options: MNE options
sProcess.options.mem.Comment = {'panel_brainentropy', 'Estimation options: '};
sProcess.options.mem.Type    = 'editpref';
sProcess.options.mem.Value   = be_main();

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


% Install/load brainentropy plugin
[isInstalled, errMessage] = bst_plugin('Install', 'brainentropy', 1);
if ~isInstalled
    return;
end
            
%% compute fusion MSP using all input trials


for iFile = 1:size(sInputs,2)
    sStudy = bst_get('Study', sInputs(iFile).iStudy);
    
    % Load head model
    if isempty(sStudy.iHeadModel)
        bst_error('No head model found. Consider running "NIRS -> Compute head model"');
        return;
    end
    nirs_head_model = in_bst_headmodel(sStudy.HeadModel(sStudy.iHeadModel).FileName);
    nirs_head_model.FileName = sStudy.HeadModel(sStudy.iHeadModel).FileName;
    cortex = in_tess_bst(nirs_head_model.SurfaceFile);
    HM.SurfaceFile = nirs_head_model.SurfaceFile;
    nb_nodes = size(cortex.Vertices, 1);
    
    % Load data & baseline
    % Load recordings
    
    if strcmp(sInputs(iFile).FileType, 'data')     % Imported data structure
        sDataIn = in_bst_data(sInputs(iFile).FileName);
    elseif strcmp(sInputs(iFile).FileType, 'raw')  % Continuous data file
        sDataIn = in_bst(sInputs(iFile).FileName, [], 1, 1, 'no');
        sDataRaw = in_bst_data(sInputs(iFile).FileName, 'F');
    end
    
    % TODO handle raw/data
    nb_samples = length(sDataIn.Time);
    ChannelMat = in_bst_channel(sInputs(iFile).ChannelFile);
    
    if ~isfield(ChannelMat.Nirs, 'Wavelengths')
        bst_error(['cMEM source reconstruction works only for dOD data ' ...
            ' (eg do not use MBLL prior to this process)']);
        return;
    end
    nb_wavelengths = length(ChannelMat.Nirs.Wavelengths);
    
    measure_tag = 'WL';
    
    MethodOptions.MEMpaneloptions =   sProcess.options.mem.Value.MEMpaneloptions;
    % Add fields that are not defined by the options of the MEM interface
    if ~isempty(MethodOptions)
        switch (nirs_head_model.HeadModelType)
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
    OPTIONS.MEMpaneloptions.solver.NoiseCov_recompute = 1;
    OPTIONS.MEMpaneloptions.model.depth_weigth_MNE = sProcess.options.depth_weightingMNE.Value{1};
    OPTIONS.MEMpaneloptions.model.depth_weigth_MEM = sProcess.options.depth_weightingMEM.Value{1};
    OPTIONS.Comment = 'cMEM';
    OPTIONS.DataFile      = sInputs(iFile).FileName;
    OPTIONS.DataTime      = round(sDataIn.Time,6);
    OPTIONS.ResultFile    = [];
    OPTIONS.HeadModelFile =  sStudy.HeadModel(sStudy.iHeadModel).FileName;
    OPTIONS.FunctionName  = 'mem';
    
    dOD_sources_cMEM = zeros(nb_nodes, nb_wavelengths, nb_samples);
    
    %sensitivity_thresh = str2double(sProcess.options.sensitivity_thresh.Value); %TODO: safe numerical conversion
    thresh_dis2cortex = sProcess.options.thresh_dis2cortex.Value{1}.*0.01;
    store_sparse_results = sProcess.options.store_sparse_results.Value;
    pbar = bst_progress('start', 'Reconstruction by cMEM', 'Launching cMEM...');
    % HM.SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
    
    for iwl=1:nb_wavelengths
        swl = [measure_tag num2str(ChannelMat.Nirs.Wavelengths(iwl))];
        selected_chans = strcmpi({ChannelMat.Channel.Group}, swl) & (sDataIn.ChannelFlag>0)';
        
        OPTIONS.GoodChannel = ones(sum(selected_chans), 1);
        OPTIONS.ChannelFlag   = ones(sum(selected_chans), 1);
        OPTIONS.Channel = ChannelMat.Channel(selected_chans);
        
        %TODO: fix gain matrix to avoid zeros
        gain = get_gains(nirs_head_model.Gain, nirs_head_model.pair_names, iwl, ...
            ChannelMat.Channel, find(selected_chans));
        %% define the reconstruction FOV
        [~, ~, ~, ~, src_locs, ~, ~, det_locs, ~, ~] = explode_channels(ChannelMat);
        optodes_pos = [src_locs;det_locs];
        % inflate surface 100% to calculate distances to optodes (see BST folder figure_3d.m line 2595)
        iVertices = 1:length(cortex.Vertices);
        % Smoothing factor
        SurfSmoothIterations = ceil(300 * 1 * length(iVertices) / 100000);
        % Calculate smoothed vertices locations
        Vertices_sm = cortex.Vertices;
        Vertices_sm(iVertices,:) = tess_smooth(cortex.Vertices(iVertices,:), 1, SurfSmoothIterations, cortex.VertConn(iVertices,iVertices), 1);
        dis2cortex = pdist2(Vertices_sm,optodes_pos);
        valid_nodes = find(min(dis2cortex,[],2)<thresh_dis2cortex);
        
        
        OPTIONS.MEMpaneloptions.optional.cortex_vertices = cortex.Vertices(valid_nodes, :);
        HM.vertex_connectivity = cortex.VertConn(valid_nodes, valid_nodes);
        
        HM.Gain = gain(:,valid_nodes);
        HM.Gain(HM.Gain==0) = min(HM.Gain(HM.Gain>0));
        pbar = bst_progress('text', ['WL' num2str(iwl) ', kept ' num2str(length(valid_nodes)) ...
            ' nodes that were in VOI and have non-zero sensitivity']);
        OPTIONS.Data = sDataIn.F(selected_chans,:);
        
        % estimate the neighborhood order for cMEM  (goal: # of clusters ~= # of good channels)
        flag_auto_nbo = sProcess.options.auto_neighborhood_order.Value;
        if flag_auto_nbo
            pbar = bst_progress('text', 'Setting neighborhood order for clustering automatically...');
            % see panel_scout.m line 2019
            % Progress bar
            % Split hemispheres
            [rH, lH] = tess_hemisplit(cortex);
            rH = intersect(rH,valid_nodes);
            lH = intersect(lH,valid_nodes);
            if ~isempty(rH) || ~isempty(lH)
                isConnected = 0;
            else
                isConnected = 1;
            end
            % Clustering
            nClust = length(OPTIONS.GoodChannel);
            isRandom = 1;
            if isConnected
                Labels = tess_cluster(cortex.VertConn(valid_nodes,valid_nodes), nClust, isRandom);
            else
                ratiolH = length(lH)./(length(lH)+length(rH));
                ratiorH = length(rH)./(length(lH)+length(rH));
                Labels = ones(length(valid_nodes), 1);
                Labels(ismember(valid_nodes,lH)) = tess_cluster(cortex.VertConn(lH,lH), round(nClust.*ratiolH), isRandom);
                Labels(ismember(valid_nodes,rH)) = tess_cluster(cortex.VertConn(rH,rH), round(nClust.*ratiorH), isRandom) + max(Labels(ismember(valid_nodes,lH)));
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
            OPTIONS.MEMpaneloptions.clustering.neighborhood_order = round(mean(neighborhood_order))+2;
        end
        
        % launch MEM (cMEM only in current version)
        pbar = bst_progress('text', ['Fusing MSP - wavelength #' num2str(iwl) '...']);
        [scr_wav, O_updated] = be_main_call_NIRS_MSP_fusion(HM, OPTIONS);
        SCR{iFile,iwl} = scr_wav;
    end
    
end
%% reconstruction using MEM and fusion MSP
for iFile = 1:size(sInputs,2)
    sStudy = bst_get('Study', sInputs(iFile).iStudy);
    
    %% Load head model
    if isempty(sStudy.iHeadModel)
        bst_error('No head model found. Consider running "NIRS -> Compute head model"');
        return;
    end
    nirs_head_model = in_bst_headmodel(sStudy.HeadModel(sStudy.iHeadModel).FileName);
    nirs_head_model.FileName = sStudy.HeadModel(sStudy.iHeadModel).FileName;
    cortex = in_tess_bst(nirs_head_model.SurfaceFile);
    HM.SurfaceFile = nirs_head_model.SurfaceFile;
    nb_nodes = size(cortex.Vertices, 1);
    
    %% Load data & baseline
    % Load recordings
    
    if strcmp(sInputs(iFile).FileType, 'data')     % Imported data structure
        sDataIn = in_bst_data(sInputs(iFile).FileName);
    elseif strcmp(sInputs(iFile).FileType, 'raw')  % Continuous data file
        sDataIn = in_bst(sInputs(iFile).FileName, [], 1, 1, 'no');
        sDataRaw = in_bst_data(sInputs(iFile).FileName, 'F');
    end
    
    % TODO handle raw/data
    nb_samples = length(sDataIn.Time);
    ChannelMat = in_bst_channel(sInputs(iFile).ChannelFile);
    
    if ~isfield(ChannelMat.Nirs, 'Wavelengths')
        bst_error(['cMEM source reconstruction works only for dOD data ' ...
            ' (eg do not use MBLL prior to this process)']);
        return;
    end
    nb_wavelengths = length(ChannelMat.Nirs.Wavelengths);
    
    measure_tag = 'WL';
    
    % baseline_id = find(cell2mat(cellfun(@(a) ~isempty(a), ...
    %                    strfind({sStudy.Data.Comment}, 'baseline'),...
    %                    'uni', false )));
    % if isempty(baseline_id)
    %     error('no baseline found');
    % end
    % baseline_data = in_bst_data(sStudy.Data(baseline_id).FileName);
    
    % iItem = [];
    % iItem = [iItem sInputs.iItem];
    % if length(iItem) == 1
    % % Default options
    % MethodOptions = be_main();
    % % Interface to edit options
    % MethodOptions = gui_show_dialog('MEM options', @panel_brainentropy, [], [], MethodOptions);
    % end
    
    MethodOptions.MEMpaneloptions =   sProcess.options.mem.Value.MEMpaneloptions;
    % Add fields that are not defined by the options of the MEM interface
    if ~isempty(MethodOptions)
        switch (nirs_head_model.HeadModelType)
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
    OPTIONS.MEMpaneloptions.solver.NoiseCov_recompute = 1;
    OPTIONS.MEMpaneloptions.model.depth_weigth_MNE = sProcess.options.depth_weightingMNE.Value{1};
    OPTIONS.MEMpaneloptions.model.depth_weigth_MEM = sProcess.options.depth_weightingMEM.Value{1};
    
    %% Old baseline definition code -> see how it can be defined from here and injected through new panel option
    % O.optional.TimeSegment                = [sDataIn.Time(1) sDataIn.Time(end)];
    % O.optional.BaselineSegment            = [baseline_data.Time(1) baseline_data.Time(end)];
    % O.optional.BaselineTime = baseline_data.Time;
    
    % Verify how to consider 2 differents baseline length (useless ???
    % Maybe we can check and take the smallest baseline ...
    % Define Baseline length
    % s = size(EEG_baseline,2);   % Number of sample on the baseline
    % ds = Time(2)-Time(1);       % Interval of time between these samples
    % base_time = [0:s-1]*ds;     % Define a fictive window of time of the baseline
    % O.optional.BaselineTime               = base_time;
    %*************************************************************************
    
    %% Run cMEM
    
    % % Separate NIRS channels from others (NIRS_AUX etc.)
    % [fnirs, fchannel_def, nirs_other, channel_def_other] = ...
    %     filter_data_by_channel_type(sDataIn.F', ChannelMat, 'NIRS');
    % HM.Gain.modality = 'MEG';
    
    OPTIONS.Comment = 'cMEM';
    OPTIONS.DataFile      = sInputs(iFile).FileName;
    OPTIONS.DataTime      = round(sDataIn.Time,6);
    OPTIONS.ResultFile    = [];
    OPTIONS.HeadModelFile =  sStudy.HeadModel(sStudy.iHeadModel).FileName;
    OPTIONS.FunctionName  = 'mem';
    
    dOD_sources_cMEM = zeros(nb_nodes, nb_wavelengths, nb_samples);
    
    %sensitivity_thresh = str2double(sProcess.options.sensitivity_thresh.Value); %TODO: safe numerical conversion
    thresh_dis2cortex = sProcess.options.thresh_dis2cortex.Value{1}.*0.01;
    store_sparse_results = sProcess.options.store_sparse_results.Value;
    pbar = bst_progress('start', 'Reconstruction by cMEM', 'Launching cMEM...');
    % HM.SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
    for iwl=1:nb_wavelengths
        swl = [measure_tag num2str(ChannelMat.Nirs.Wavelengths(iwl))];
        selected_chans = strcmpi({ChannelMat.Channel.Group}, swl) & (sDataIn.ChannelFlag>0)';
        
        OPTIONS.GoodChannel = ones(sum(selected_chans), 1);
        OPTIONS.ChannelFlag   = ones(sum(selected_chans), 1);
        OPTIONS.Channel = ChannelMat.Channel(selected_chans);
        
        %TODO: fix gain matrix to avoid zeros
        gain = get_gains(nirs_head_model.Gain, nirs_head_model.pair_names, iwl, ...
            ChannelMat.Channel, find(selected_chans));
        %% define the reconstruction FOV
        [~, ~, ~, ~, src_locs, ~, ~, det_locs, ~, ~] = explode_channels(ChannelMat);
        optodes_pos = [src_locs;det_locs];
        % inflate surface 100% to calculate distances to optodes (see BST folder figure_3d.m line 2595)
        iVertices = 1:length(cortex.Vertices);
        % Smoothing factor
        SurfSmoothIterations = ceil(300 * 1 * length(iVertices) / 100000);
        % Calculate smoothed vertices locations
        Vertices_sm = cortex.Vertices;
        Vertices_sm(iVertices,:) = tess_smooth(cortex.Vertices(iVertices,:), 1, SurfSmoothIterations, cortex.VertConn(iVertices,iVertices), 1);
        dis2cortex = pdist2(Vertices_sm,optodes_pos);
        valid_nodes = find(min(dis2cortex,[],2)<thresh_dis2cortex);
        
        
        OPTIONS.MEMpaneloptions.optional.cortex_vertices = cortex.Vertices(valid_nodes, :);
        HM.vertex_connectivity = cortex.VertConn(valid_nodes, valid_nodes);
        
        HM.Gain = gain(:,valid_nodes);
        HM.Gain(HM.Gain==0) = min(HM.Gain(HM.Gain>0));
        pbar = bst_progress('text', ['WL' num2str(iwl) ', kept ' num2str(length(valid_nodes)) ...
            ' nodes that were in VOI and have non-zero sensitivity']);
        %     disp(['WL' num2str(iwl) ', kept ' num2str(length(valid_nodes)) ...
        %              ' nodes that were in VOI and have non-zero sensitivity']);
        
        OPTIONS.Data = sDataIn.F(selected_chans,:);
        
        %% estimate the neighborhood order for cMEM  (goal: # of clusters ~= # of good channels)
        flag_auto_nbo = sProcess.options.auto_neighborhood_order.Value;
        if flag_auto_nbo
            pbar = bst_progress('text', 'Setting neighborhood order for clustering automatically...');
            % see panel_scout.m line 2019
            % Progress bar
            % Split hemispheres
            [rH, lH] = tess_hemisplit(cortex);
            rH = intersect(rH,valid_nodes);
            lH = intersect(lH,valid_nodes);
            if ~isempty(rH) || ~isempty(lH)
                isConnected = 0;
            else
                isConnected = 1;
            end
            % Clustering
            nClust = length(OPTIONS.GoodChannel);
            isRandom = 1;
            if isConnected
                Labels = tess_cluster(cortex.VertConn(valid_nodes,valid_nodes), nClust, isRandom);
            else
                ratiolH = length(lH)./(length(lH)+length(rH));
                ratiorH = length(rH)./(length(lH)+length(rH));
                Labels = ones(length(valid_nodes), 1);
                Labels(ismember(valid_nodes,lH)) = tess_cluster(cortex.VertConn(lH,lH), round(nClust.*ratiolH), isRandom);
                Labels(ismember(valid_nodes,rH)) = tess_cluster(cortex.VertConn(rH,rH), round(nClust.*ratiorH), isRandom) + max(Labels(ismember(valid_nodes,lH)));
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
            OPTIONS.MEMpaneloptions.clustering.neighborhood_order = round(mean(neighborhood_order))+2;
        end
        
        %% launch MEM (cMEM only in current version)
        pbar = bst_progress('text', ['Running cMEM for wavelength #' num2str(iwl) '...']);
        
        %%%%% fusion option %%%%%
        OPTIONS.MEMpaneloptions.SCR = (SCR{1,1}+SCR{1,2}-SCR{1,1}.*SCR{1,2}).*(SCR{2,1}+SCR{2,2}-SCR{2,1}.*SCR{2,2});
        OPTIONS.MEMpaneloptions.SCR = OPTIONS.MEMpaneloptions.SCR./max(OPTIONS.MEMpaneloptions.SCR(:));
%         if iwl == 1
%             OPTIONS.MEMpaneloptions.SCR = SCR{1,1}.*SCR{2,1};
%             OPTIONS.MEMpaneloptions.SCR = OPTIONS.MEMpaneloptions.SCR./max(OPTIONS.MEMpaneloptions.SCR(:));
%         else
%             OPTIONS.MEMpaneloptions.SCR = SCR{1,2}.*SCR{2,2};
%             OPTIONS.MEMpaneloptions.SCR = OPTIONS.MEMpaneloptions.SCR./max(OPTIONS.MEMpaneloptions.SCR(:));
%         end
        
        [Results, O_updated] = be_main_call_NIRS_fusion(HM, OPTIONS);
        
        %cMEM results
        grid_amp = zeros(nb_nodes, nb_samples); %zeros(nb_nodes, nb_samples);
        grid_amp(valid_nodes,:) = Results.ImageGridAmp;
        
        dOD_sources_cMEM(:, iwl, :) = grid_amp;
        
        [sStudy, ResultFile] = add_surf_data(grid_amp, sDataIn.Time, nirs_head_model, ...
            ['cMEM sources - ' swl 'nm'], ...
            sInputs(iFile), sStudy, 'cMEM sources reconstruction', ...
            'OD', store_sparse_results,Results.MEMoptions);
        OutputFiles{end+1} = ResultFile;
        
        %MNE results
        %     grid_amp = zeros(nb_nodes, nb_samples);
        %     grid_amp(valid_nodes,:) = Results.MEMoptions.automatic.minimum_norm;
        %
        %     dOD_sources_MNE(:, iwl, :) = grid_amp;
        %
        %     [sStudy, ResultFile] = add_surf_data(grid_amp, sDataIn.Time, nirs_head_model, ...
        %                                         ['MNE sources - ' swl 'nm'], ...
        %                                         sInputs, sStudy, 'MNE sources reconstruction', ...
        %                                         'OD', store_sparse_results);
        %
        %     OutputFiles{end+1} = ResultFile;
    end
    pbar = bst_progress('text', 'Calculating HbO/HbR/HbT in source space...');
    % Compute dHb
    hb_extinctions = nst_get_hb_extinctions(ChannelMat.Nirs.Wavelengths);
    % hb_extinctions = [2321.4 1791.7;... % 830HbO 830HbR cm-1.mole-1.L
    %     956.89 4932.00];   % 690HbO 690HbR cm-1.mole-1.L
    hb_extinctions = hb_extinctions ./10;% mm-1.mole-1.L
    nirs_hbo_hbr = zeros(nb_nodes, 2, nb_samples);
    %nirs_hbo_hbr_mne = zeros(nb_nodes, 2, nb_samples);
    for idx=1:length(valid_nodes)
        inode = valid_nodes(idx);
        nirs_hbo_hbr(inode, :, :) = pinv(hb_extinctions) * ...
            squeeze(dOD_sources_cMEM(inode, :, :));
        %     %nirs_hbo_hbr_mne(inode, :, :) = pinv(hb_extinctions) * ...
        %                                 squeeze(dOD_sources_MNE(inode, :, :));
    end
    nirs_hbt = squeeze(sum(nirs_hbo_hbr, 2));
    
    hb_unit_factor = 1e6;
    hb_unit = '\mumol.l-1';
    
    hb_types = {'HbO', 'HbR'};
    for ihb=1:2
        [sStudy, ResultFile] = add_surf_data(squeeze(nirs_hbo_hbr(:,ihb,:)) .* hb_unit_factor,...
            sDataIn.Time, nirs_head_model, ...
            ['cMEM sources - ' hb_types{ihb}], ...
            sInputs(iFile), sStudy, 'cMEM sources reconstruction - dHb', ...
            hb_unit, store_sparse_results);
        OutputFiles{end+1} = ResultFile;
    end
    
    [sStudy, ResultFile] = add_surf_data(nirs_hbt .* hb_unit_factor, sDataIn.Time, nirs_head_model, ...
        'cMEM sources - HbT', ...
        sInputs(iFile), sStudy, 'cMEM sources reconstruction - dHb', ...
        hb_unit, store_sparse_results);
    OutputFiles{end+1} = ResultFile;
    pbar = bst_progress('stop', 'Reconstruction by cMEM', 'Finishing...');
    % Update Brainstorm database
    bst_set('Study', sInputs(iFile).iStudy, sStudy);
end
end

function gains = get_gains(gain_in, gain_pair_names, iwl, channel_def, selected_chans)

% hash gain pair names
gpair_idx = containers.Map();
for gpic=1:length(gain_pair_names)
    gpair_idx(gain_pair_names{gpic}) = gpic;
end

if ~isempty(strfind(channel_def(1).Name, 'WL'))
    measure_tag = 'WL';
else
    measure_tag = 'Hb';
end

gains = zeros(length(selected_chans), size(gain_in, 3));
for ic=1:length(selected_chans)
    ichan = selected_chans(ic);
    chan_name = channel_def(ichan).Name;
    pair_name = chan_name(1:strfind(chan_name, 'WL')-1);
    if ~isempty(pair_name) && gpair_idx.isKey(pair_name)
        gains(ic, :) = squeeze(gain_in(gpair_idx(pair_name), iwl, :));
    end
end

end

function [sStudy, ResultFile] = add_surf_data(data, time, head_model, name, ...
    sInputs, sStudy, history_comment, ...
    data_unit, store_sparse, MEMoptions)

if nargin < 8
    data_unit = '';
end

if nargin < 10
    MEMoptions = [];
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
ResultsMat.MEMoptions = MEMoptions;
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

function [pair_names, pair_loc, pair_ichans, pair_sd_indexes, ...
    src_coords, src_ids, src_ichans, ...
    det_coords, det_ids, det_ichans] = explode_channels(channel_def)
%% Explode channel data according to pairs, sources and detectors
% Args
%    - channel_def: struct
%        Definition of channels as given by brainstorm
%        Used fields: Channel
%
% TOCHECK WARNING: uses containers.Map which is available with matlab > v2008
%
%  Outputs:
%     - pair_names: cell array of str, size: nb_pairs
%         Pair names, format: SXDX
%     - pair_loc: array of double, size: nb_pairs x 3 x 2
%         Pair localization (coordinates of source and detector)
%     - pair_ichans: matrix of double, size: nb_pairs x nb_wavelengths
%         Input channel indexes grouped by pairs
%     - pair_sd_indexes: matrix of double, size: nb_pairs x 2
%         1-based continuours indexes of sources and detectors for each
%         sources.
%     - src_coords:   nb_sources x 3
%         Source coordinates, indexed by 1-based continuous index
%         To access via source ID, as read from pair name:
%             src_coords(src_id2idx(src_ID),:)
%     - src_ids: 1d array of double, size: nb_sources
%         vector of source ids (as used in pair name)
%     - src_chans: cellarray of 1d array of double, size: nb_sources
%         Channel indexes to which the source belongs (indexed by 1-based
%         continuous index).
%     - det_coords:   nb_detectors x 3
%         Detector coordinates, indexed by 1-based continuous index
%         To access via detector ID, as used in pair name:
%             det_coords(det_id2idx(det_ID),:)
%     - det_ids: 1d array of double, size: max_detector_id (hashing vector)
%         vector of detector ids (as used in pair name)
%     - det_chans: cellarray of 1d array of double, size: nb_sources
%         Channel indexes to which the detector belongs (indexed by 1-based
%         continuous index).

MT_OD = 1;
MT_HB = 2;

if isfield(channel_def.Nirs, 'Wavelengths')
    nb_measures = length(channel_def.Nirs.Wavelengths);
    measure_type = MT_OD;
else
    nb_measures = length(channel_def.Nirs.Hb);
    measure_type = MT_HB;
end

pair_to_chans = containers.Map();
pair_to_sd = containers.Map();
src_to_chans = containers.Map('KeyType', 'double', 'ValueType', 'any');
src_coords_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
det_to_chans = containers.Map('KeyType', 'double', 'ValueType', 'any');
det_coords_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
for ichan=1:length(channel_def.Channel)
    if strcmp(channel_def.Channel(ichan).Type, 'NIRS')
        chan_name = channel_def.Channel(ichan).Name;
        if measure_type == MT_OD
            iwl = strfind(chan_name, 'WL');
            pair_name = chan_name(1:iwl-1);
            wl = str2double(chan_name(iwl+2:end));
            imeasure = channel_def.Nirs.Wavelengths==wl;
        else
            ihb = strfind(chan_name, 'Hb');
            pair_name = chan_name(1:ihb-1);
            imeasure = strcmp(chan_name(ihb:end), channel_def.Nirs.Hb);
        end
        
        if pair_to_chans.isKey(pair_name)
            measures = pair_to_chans(pair_name);
        else
            measures = zeros(1, nb_measures);
        end
        measures(imeasure) = ichan;
        pair_to_chans(pair_name) = measures;
        
        
        [src_id, det_id] = split_pair_name(pair_name);
        pair_to_sd(pair_name) = [src_id, det_id];
        if src_to_chans.isKey(src_id)
            src_to_chans(src_id) = [src_to_chans(src_id) ichan];
        else
            src_to_chans(src_id) = ichan;
            src_coords_map(src_id) = channel_def.Channel(ichan).Loc(:, 1);
        end
        if det_to_chans.isKey(det_id)
            det_to_chans(det_id) = [det_to_chans(det_id) ichan];
        else
            det_to_chans(det_id) = ichan;
            det_coords_map(det_id) = channel_def.Channel(ichan).Loc(:, 2);
        end
        
    end
end

src_coords = cell2mat(src_coords_map.values)';
src_ichans = src_to_chans.values;
src_ids = cell2mat(src_coords_map.keys);

det_coords = cell2mat(det_coords_map.values)';
det_ichans = det_to_chans.values;
det_ids = cell2mat(det_coords_map.keys);

nb_pairs = pair_to_chans.size(1);
pair_names = pair_to_chans.keys;
pair_ichans = zeros(nb_pairs, nb_measures);
pair_loc = zeros(nb_pairs, 3, 2);
pair_sd_indexes = zeros(nb_pairs, 2);
for ipair=1:nb_pairs
    p_indexes = pair_to_chans(pair_names{ipair});
    pair_ichans(ipair, :) = p_indexes;
    pair_loc(ipair, : , :) = channel_def.Channel(pair_ichans(ipair, 1)).Loc;
    sdi = pair_to_sd(pair_names{ipair});
    pair_sd_indexes(ipair, 1) = find(src_ids==sdi(1));
    pair_sd_indexes(ipair, 2) = find(det_ids==sdi(2));
end

end

function [isrc, idet] = split_pair_name(pair_name)
pair_re = 'S([0-9]{1,2})D([0-9]{1,2})';
toks = regexp(pair_name, pair_re , 'tokens');
isrc = str2double(toks{1}{1});
idet = str2double(toks{1}{2});
end
