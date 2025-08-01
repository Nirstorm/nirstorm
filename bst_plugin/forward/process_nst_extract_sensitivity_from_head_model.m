function varargout = process_nst_extract_sensitivity_from_head_model( varargin )

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
% Authors: Edouard Delaire (2023) 
%          Thomas Vincent, ZhengChen Cai (2017)

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Extract sensitivity surfaces from head model';
    sProcess.FileTag     = '';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = {'NIRS', 'Sources'};
    sProcess.Index       = 1405;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'results', 'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;


    % === Process description
    sProcess.options.label1.Comment = ['<b>Light sensitivity scale.</b> <BR>' ...
                                        'Export the nirs global and channel specific sensitivity based on the head model<BR>', ...
                                       'Note: if the sensitivity is exported in the log scale (db), the map is thresholded at -2db.'];
    sProcess.options.label1.Type = 'label';


    sProcess.options.method.Comment = {['Linear'], ...
                                       ['Scale with the global max (dB): <FONT color=#7F7F7F>&nbsp;&nbsp;&nbsp;' ...
                                            'sensitivity = log10(x / max(sensitivity))'],... 
                                       ['Scale with the local max (dB) (only for channel specific sensitivity map): <FONT color=#7F7F7F>&nbsp;&nbsp;&nbsp;' ...
                                            'For each channel, sensitivity = log10(x / max(sensitivity))</FONT>'] ;...
                                       'linear', 'db_global','db_local'};
    sProcess.options.method.Type    = 'radio_label';
    sProcess.options.method.Value   = 'linear';


    sProcess.options.label2.Comment = ['<b>Measure of overlap</b> <BR>'];
    sProcess.options.label2.Type = 'label';

    sProcess.options.export_overlap.Comment = 'Export measure of overlap';
    sProcess.options.export_overlap.Type    = 'checkbox';
    sProcess.options.export_overlap.Value   = 1;    


    
    sProcess.options.label3.Comment = ['<b>NIRS field of view (FOV)</b> <BR>' ...
                                       'Export the FOV used for source reconstruction as new atlas'];
    sProcess.options.label3.Type = 'label';

    sProcess.options.export_FOV.Comment = 'Export NIRS FOV';
    sProcess.options.export_FOV.Type    = 'checkbox';
    sProcess.options.export_FOV.Value   = 1;    
    sProcess.options.export_FOV.Controller = 'FOV';

    % Definition of the options
    sProcess.options.thresh_dis2cortex.Comment = 'Reconstruction Field of view (distance to montage border)';
    sProcess.options.thresh_dis2cortex.Type    = 'value';
    sProcess.options.thresh_dis2cortex.Value   = {3, 'cm',2};
    sProcess.options.thresh_dis2cortex.Class   = 'FOV';


end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

OutputFiles = {};

% Load Head Model
sStudy = bst_get('Study', sInputs.iStudy);
if isempty(sStudy.iHeadModel)
    bst_error('No head model found. Consider process "Compute head model from fluence"');
    return;
end

head_model = in_bst_headmodel(sStudy.HeadModel(sStudy.iHeadModel).FileName, 1);
if ~isfield(head_model, 'NIRSMethod') && ndims(head_model.Gain) == 3
    head_model = process_nst_import_head_model('convert_head_model', ChannelMat, head_model, 0);
end

if ~strcmp(head_model.HeadModelType, 'surface')
    bst_error('Extraction only works for surface head model');
    return;
end

% Initialize results
sResults        = [];

% Load ChannelFlag
bst_chan_data   = load(file_fullpath(sInputs.FileName), 'ChannelFlag');
ChannelFlag     = bst_chan_data.ChannelFlag;
ChannelMat      = in_bst_channel(sInputs(1).ChannelFile);


% Load Cortex 
sSubject    = bst_get('Subject', sInputs.SubjectName);
sCortex     = in_tess_bst(head_model.SurfaceFile);

% Compute sensitivity map
sSensitivity    = get_sensitivity_map(head_model, ChannelMat, ChannelFlag, sProcess.options.method.Value);
sResults        = [sResults,  sSensitivity];


% Estimate Coverage
voronoi_fn  = process_nst_compute_voronoi('get_voronoi_fn', sSubject);
if ~exist(voronoi_fn, 'file')
    error('Could not find the required Voronoi file.');
end

%threshold for coverage
p_thresh    = 1;
act_vol     = 1000; % A definir comme un parametre donne par l'utilisateur
sVoronoi    = in_mri_bst(voronoi_fn);

median_voronoi_volume = process_nst_compute_voronoi('get_median_voronoi_volume', sVoronoi);  
delta_mu_a = 0.1;
threshold = compute_threshold(p_thresh, act_vol, median_voronoi_volume, delta_mu_a);

sCoverage = getCoverage(head_model, ChannelMat, ChannelFlag, threshold);
sResults  = [sResults,  sCoverage];



% Save results
for iMap = 1:length(sResults)

    ResultFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName),  ['results_NIRS_' nst_protect_fn_str(sResults(iMap).Comment)]);
    ResultsMat          = sResults(iMap);
    %ResultsMat.Options  = OPTIONS;

    bst_save(ResultFile, ResultsMat, 'v6');
    db_add_data( sInputs.iStudy, ResultFile, ResultsMat);

    OutputFiles{end+1} = ResultFile;
end


% Save the NIRS FOV
thresh_dis2cortex           = sProcess.options.thresh_dis2cortex.Value{1}*0.01;
[valid_nodes,dis2cortex]    = nst_headmodel_get_FOV(ChannelMat, sCortex, thresh_dis2cortex, ChannelFlag);

if any(strcmp({sCortex.Atlas.Name},'NIRS-FOV'))
    iAtlas = find(strcmp({sCortex.Atlas.Name},'NIRS-FOV'));
else
    sCortex.Atlas(end+1).Name = 'NIRS-FOV';
    iAtlas = length( sCortex.Atlas);
end


sCortex.Atlas(iAtlas).Scouts(end+1)              = db_template('Scout'); 
sCortex.Atlas(iAtlas).Scouts(end).Vertices       = valid_nodes;
sCortex.Atlas(iAtlas).Scouts(end).Seed           = valid_nodes(1);
sCortex.Atlas(iAtlas).Scouts(end).Label          = sprintf('NIRS FOV (%d cm)',sProcess.options.thresh_dis2cortex.Value{1} );
sCortex.Atlas(iAtlas).Scouts(end)                = panel_scout('SetColorAuto',sCortex.Atlas(iAtlas).Scouts(end), length(sCortex.Atlas(iAtlas).Scouts));

bst_save(file_fullpath(head_model.SurfaceFile), sCortex)
end

function threshold = compute_threshold(p_thresh, act_vol, V_hat, delta_mu_a)
% @========================================================================
% compute_threshold computes the threshold to determine the coverage matrix 
% Formula used : threshold = log((100 + p_thresh) / 100) / ((act_vol / V_hat) * delta_mu_a)
% ========================================================================@
        
        numerator = log10((100 + p_thresh) / 100);
        denomimator = (act_vol / V_hat) * delta_mu_a;
        threshold = numerator / denomimator;

end

function  sResults = get_sensitivity_map(head_model, ChannelMat, ChannelFlag, normalization)
    
    % load cortex
    sCortex     = in_tess_bst(head_model.SurfaceFile);


    % Get Montage information
    montage_info    = nst_montage_info_from_bst_channels(ChannelMat.Channel, ChannelFlag);
    max_sources     = max(montage_info.src_ids);
    max_dets        = max(montage_info.det_ids);
    nb_nodes        = size(sCortex.Vertices   , 1);
    nb_Wavelengths  = length(ChannelMat.Nirs.Wavelengths);

    time        = 1:(max_sources*100 + max_dets);
    isUsedTime  = zeros(1, length(time));
    
    sensitivity_surf        = zeros(nb_nodes, nb_Wavelengths, length(time));
    sensitivity_surf_sum    = zeros(nb_nodes, nb_Wavelengths);

    %% Compute sensitivity
    for iwl = 1:nb_Wavelengths
    
        swl = ['WL' num2str(ChannelMat.Nirs.Wavelengths(iwl))];
        selected_chans = strcmpi({ChannelMat.Channel.Group}, swl) & (ChannelFlag>0)';
        idx_chan       = find(selected_chans);
    
        sensitivity         = head_model.Gain(idx_chan, :);
    
    
        for iChan = 1:length(idx_chan)
            chan = ChannelMat.Channel(idx_chan(iChan));
            [src_id, det_id] = nst_unformat_channel(chan.Name );
    
            sensitivity_surf(:, iwl, det_id + src_id*100) = squeeze(sensitivity(iChan,:));
            isUsedTime(det_id + src_id*100)              = 1;
        end
    
        sensitivity_surf_sum(:,iwl) = sum(sensitivity,  1) ;
    end

    %% Normalize values and threshold 
    threshold_value = -2; % in db
    
    if contains(normalization,'db')
        for iwl=1:nb_Wavelengths
            sensitivity_surf_sum(:,iwl) = log10(sensitivity_surf_sum(:,iwl) ./ ( eps + max(sensitivity_surf_sum(:,iwl))));
            mask = zeros(size(sensitivity_surf_sum));
            mask(:,iwl) = sensitivity_surf_sum(:,iwl) < threshold_value;
            sensitivity_surf_sum(mask == 1) = 0;
    
            k = zeros(1,  size(sensitivity_surf,3));
    
            if strcmp(normalization,'db_local') % channel wise normalization
                k(1,:) =  squeeze(max(sensitivity_surf(:, iwl, :))) + eps;     
            else % Global normalization 
                k(1,:) = max(max(sensitivity_surf(:, iwl, :))) + eps;   
            end    
            
            sensitivity_surf(:,iwl,isUsedTime == 1) = log10( squeeze(sensitivity_surf(:, iwl, isUsedTime == 1)) ./ repmat(k(isUsedTime == 1),nb_nodes,1));
            
            mask = zeros(size(sensitivity_surf));
            mask(:,iwl,:) = sensitivity_surf(:,iwl,:) < threshold_value;
            sensitivity_surf(mask == 1 ) = 0;
        end

        displayUnit = 'dB';
    else % linear

        for iwl=1:nb_Wavelengths
            mask = zeros(size(sensitivity_surf_sum));
            mask(:,iwl) = sensitivity_surf_sum(:,iwl) < 10^(threshold_value)*max(sensitivity_surf_sum(:,iwl));
            sensitivity_surf_sum(mask == 1) = 0;
    
            mask = zeros(size(sensitivity_surf));
            mask(:,iwl,:) = sensitivity_surf(:,iwl,:) <  10^(threshold_value)*max(sensitivity_surf(:,iwl,:),[],'all');
            sensitivity_surf(mask == 1 ) = 0;
        end

        displayUnit = 'mm';

    end   

    sResults               = repmat(db_template('resultsmat'), 1, 2*nb_Wavelengths);

    % Save sensitivity maps
    for iwl =  1:nb_Wavelengths
        % Store restult as brainstrom structure 
        sResults(iwl).Comment       = ['Sensitivities - WL' num2str(iwl)];
        sResults(iwl).ImageGridAmp  = squeeze(sensitivity_surf(:,iwl,:));
        sResults(iwl).Time          = time;

        sResults(iwl).DisplayUnits  = displayUnit;
        sResults(iwl).ChannelFlag   = ChannelFlag;
        %sResults(iwl).HeadModelFile = OPTIONS.HeadModelFile;
        sResults(iwl).HeadModelType = head_model.HeadModelType;
        sResults(iwl).SurfaceFile   = file_short(head_model.SurfaceFile);
        %sResults(1).History       = OPTIONS.History;
        sResults(iwl) = bst_history('add', sResults(iwl), 'compute', 'sensitivity imported from MCXlab');
    end

    % Save summed sensitivity maps
    for iwl =  1:nb_Wavelengths

        % Store restult as brainstrom structure 
        sResults(nb_Wavelengths + iwl).Comment       =  ['Summed sensitivities - WL' num2str(iwl)];
        sResults(nb_Wavelengths + iwl).ImageGridAmp  = squeeze(sensitivity_surf_sum(:,iwl));
        sResults(nb_Wavelengths + iwl).Time          = [1];

        sResults(nb_Wavelengths + iwl).DisplayUnits  = displayUnit;
        sResults(nb_Wavelengths + iwl).ChannelFlag   = ChannelFlag;
        %sResults(iwl).HeadModelFile = OPTIONS.HeadModelFile;
        sResults(nb_Wavelengths + iwl).HeadModelType = head_model.HeadModelType;
        sResults(nb_Wavelengths + iwl).SurfaceFile   = file_short(head_model.SurfaceFile);
        %sResults(1).History       = OPTIONS.History;
        sResults(nb_Wavelengths + iwl) = bst_history('add', sResults(nb_Wavelengths + iwl), 'compute', 'Summed sensitivity imported from MCXlab');
    end


end

function sResults = getCoverage(head_model, ChannelMat, ChannelFlag, threshold)
    
    sResults = [];

    % load cortex
    sCortex     = in_tess_bst(head_model.SurfaceFile);


    % Get Montage information
    montage_info    = nst_montage_info_from_bst_channels(ChannelMat.Channel, ChannelFlag);
    max_sources     = max(montage_info.src_ids);
    max_dets        = max(montage_info.det_ids);
    nb_nodes        = size(sCortex.Vertices   , 1);
    nb_Wavelengths  = length(ChannelMat.Nirs.Wavelengths);

    time        = 1:(max_sources*100 + max_dets);
    isUsedTime  = zeros(1, length(time));
    
    coverage_channel        = zeros(nb_nodes, nb_Wavelengths, length(time));
    overlap                 = zeros(nb_nodes, nb_Wavelengths);

    %% Compute sensitivity
    for iwl = 1:nb_Wavelengths
    
        swl = ['WL' num2str(ChannelMat.Nirs.Wavelengths(iwl))];
        selected_chans = strcmpi({ChannelMat.Channel.Group}, swl) & (ChannelFlag>0)';
        idx_chan       = find(selected_chans);
    
        sensitivity         = head_model.Gain(idx_chan, :);
    
    
        for iChan = 1:length(idx_chan)
            chan = ChannelMat.Channel(idx_chan(iChan));
            [src_id, det_id] = nst_unformat_channel(chan.Name );
    
            coverage_channel(:, iwl, det_id + src_id*100) = squeeze(sensitivity(iChan,:) > threshold);
            isUsedTime(det_id + src_id*100)               = any(sensitivity(iChan,:) > threshold);
        end
    
        overlap(:,iwl) = sum(coverage_channel,  3) ;
    end
    

    % Save coverage maps
    sResults               = repmat(db_template('resultsmat'), 1, 2*nb_Wavelengths);

    for iwl =  1:nb_Wavelengths
        % Store restult as brainstrom structure 
        sResults(iwl).Comment       = ['Coverage - WL' num2str(iwl)];
        sResults(iwl).ImageGridAmp  = squeeze(coverage_channel(:,iwl,:));
        sResults(iwl).Time          = time;

        sResults(iwl).DisplayUnits  = '';
        sResults(iwl).ChannelFlag   = ChannelFlag;
        %sResults(iwl).HeadModelFile = OPTIONS.HeadModelFile;
        sResults(iwl).HeadModelType = head_model.HeadModelType;
        sResults(iwl).SurfaceFile   = file_short(head_model.SurfaceFile);
        %sResults(1).History       = OPTIONS.History;
        sResults(iwl) = bst_history('add', sResults(iwl), 'compute', 'Computed channel coverage');
    end

    % Save summed sensitivity maps
    for iwl =  1:nb_Wavelengths

        % Store restult as brainstrom structure 
        sResults(nb_Wavelengths + iwl).Comment       =  ['Overlap measure - WL' num2str(iwl)];
        sResults(nb_Wavelengths + iwl).ImageGridAmp  = squeeze(overlap(:,iwl));
        sResults(nb_Wavelengths + iwl).Time          = [1];

        sResults(nb_Wavelengths + iwl).DisplayUnits  = '# channel';
        sResults(nb_Wavelengths + iwl).ChannelFlag   = ChannelFlag;
        %sResults(iwl).HeadModelFile = OPTIONS.HeadModelFile;
        sResults(nb_Wavelengths + iwl).HeadModelType = head_model.HeadModelType;
        sResults(nb_Wavelengths + iwl).SurfaceFile   = file_short(head_model.SurfaceFile);
        %sResults(1).History       = OPTIONS.History;
        sResults(nb_Wavelengths + iwl) = bst_history('add', sResults(nb_Wavelengths + iwl), 'compute', 'Computed overlap measure');
    end


end