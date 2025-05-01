function varargout = process_nst_import_head_model( varargin )

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
% Authors: Thomas Vincent, ZhengChen Cai (2017-2018)

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Compute head model from fluence';
    sProcess.FileTag     = '';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = {'NIRS', 'Sources'};
    sProcess.Index       = 1404;
    sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/Compute-head-model-from-fluence';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;


    sProcess.options.label.Comment = '<B>Fluences:</B>';
    sProcess.options.label.Type    = 'label';

    % Definition of the options    
    sProcess.options.data_source.Comment = 'Fluence Data Source (URL or path)';
    sProcess.options.data_source.Type    = 'text';
    sProcess.options.data_source.Value = [nst_get_repository_url(), '/fluence/'];
    

    % === FWHM (kernel size)

    sProcess.options.label1.Comment = '<B>Smoothing Method:</B>';
    sProcess.options.label1.Type    = 'label';

    sProcess.options.method.Comment = {'<FONT color="#777777">Before 2023 (not recommended)</FONT>', ...
                                       '<B> Geodesic</B> <FONT color="#777777">(recommended)</FONT>'; ...
                                       'surfstat_before_2023', 'geodesic_dist'};
    sProcess.options.method.Type    = 'radio_label';
    sProcess.options.method.Value   = 'surfstat_before_2023';

    sProcess.options.smoothing_fwhm.Comment = 'Spatial smoothing FWHM: ';
    sProcess.options.smoothing_fwhm.Type    = 'value';
    sProcess.options.smoothing_fwhm.Value   = {0, 'mm', 2};
  
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

OutputFiles = {sInputs.FileName};


ChannelMat = in_bst_channel(sInputs(1).ChannelFile);
if ~isfield(ChannelMat.Nirs, 'Wavelengths')
    bst_error('Head model importation works only for dOD data (eg do not use MBLL prior to this process)');
    return;
end


%% Retrieve list of vertex indexes corresponding to current montage
% Load channel file

% Load subject
[sSubject, iSubject] = bst_get('Subject', sInputs.SubjectName);

% Load MRI, Head and Cortex. 
if isempty(sSubject.iAnatomy)
    bst_error(['No anatomical data found for ' sInputs.SubjectName]);
end

voronoi_fn = process_nst_compute_voronoi('get_voronoi_fn', sSubject);
    
if ~exist(voronoi_fn, 'file')
    sProcess.options.do_grey_mask.Value   = 1; 
    process_nst_compute_voronoi('Run', sProcess, sInputs);
end


OPTIONS = struct();

% Subject Informations 
OPTIONS.SubjectName     = sInputs.SubjectName;
OPTIONS.MriFile         = sSubject.Anatomy(sSubject.iAnatomy).FileName;
OPTIONS.VoronoiFile     = voronoi_fn;
OPTIONS.HeadFile        = sSubject.Surface(sSubject.iScalp  ).FileName;
OPTIONS.CortexFile      = sSubject.Surface(sSubject.iCortex ).FileName;
OPTIONS.ChannelFile     = sInputs(1).ChannelFile;

% Use defined options : 
OPTIONS.FluenceFolder       = sProcess.options.data_source.Value;
OPTIONS.smoothing_method    = sProcess.options.method.Value;
OPTIONS.smoothing_fwhm      = sProcess.options.smoothing_fwhm.Value{1};

% Compute Head model
HeadModelMat = Compute(OPTIONS);


% Save Head Model
sStudy = bst_get('Study', sInputs.iStudy);

HeadModelFile = bst_fullfile(bst_fileparts(file_fullpath(sStudy.FileName)), 'headmodel_nirs_mcx_fluence.mat');
HeadModelFile = file_unique(HeadModelFile);
bst_save(HeadModelFile, HeadModelMat, 'v7');

% Update Study structure
newHeadModel = db_template('HeadModel');
newHeadModel.FileName = file_short(HeadModelFile);
newHeadModel.Comment = 'NIRS head model from fluence';
newHeadModel.HeadModelType  = 'surface';    

iHeadModel = length(sStudy.HeadModel) + 1;
if ~isempty(sStudy.HeadModel)
    sStudy.HeadModel(end+1) = newHeadModel(1);
else
    sStudy.HeadModel = newHeadModel;
end
sStudy.iHeadModel = iHeadModel;

% Update DataBase
bst_set('Study', sInputs.iStudy, sStudy);
panel_protocols('UpdateNode', 'Study', sInputs.iStudy);
db_save();

end

function [HeadModelMat, err] = Compute(OPTIONS)
    
    HeadModelMat = [];
    err = '';

    % Retrieve optode coordinates
    % Load channel file
    sMri    = in_mri_bst (OPTIONS.MriFile);
    sHead   = in_tess_bst(OPTIONS.HeadFile);
    sCortex = in_tess_bst(OPTIONS.CortexFile);

    ChannelMat = in_bst_channel(OPTIONS.ChannelFile);
    sChannels       = ChannelMat.Channel;
    
    % Get Voronoi
    voronoi_bst = in_mri_bst(OPTIONS.VoronoiFile);
    voronoi = voronoi_bst.Cube;
    voronoi_mask = (voronoi > -1) & ~isnan(voronoi);

    % Load montage informations
    montage_info    = nst_montage_info_from_bst_channels(sChannels);
    pair_names      = montage_info.pair_names;
    src_locs        = montage_info.src_pos;
    src_ids         = montage_info.src_ids;
    det_locs        = montage_info.det_pos;
    det_ids         = montage_info.det_ids;
    pair_sd_idx     = montage_info.pair_sd_indexes;
    
    nb_wavelengths = length(ChannelMat.Nirs.Wavelengths);
    nb_sources     = size(src_locs, 1);
    nb_dets        = size(det_locs, 1);
    nb_pairs       = length(pair_names);
    
    
    % Find closest head vertices (for which we have fluence data)    
    [src_hvidx, det_hvidx] = get_head_vertices_closest_to_optodes(sMri, sHead, src_locs, det_locs);
    
    %% Load fluence data from local .brainstorm folder (download if not available)
    use_closest_wl = 0;
    [all_fluences_flat_sparse, all_reference_voxels_index] = request_fluences([src_hvidx ; det_hvidx], sMri.Comment, ...
                                                                             ChannelMat.Nirs.Wavelengths, OPTIONS.FluenceFolder, nan, nan, [], '',...
                                                                             use_closest_wl);
    if isempty(all_fluences_flat_sparse)
        return;
    end
    
    
    mri_zeros       = zeros(size(sMri.Cube));
    src_fluences    = cell(1, nb_sources);
    det_fluences    = cell(1, nb_dets);

    bst_progress('start', 'Export fluences','Unpacking volumic fluences...', 1, (nb_sources+nb_dets)*nb_wavelengths);
    for iwl = 1:nb_wavelengths
        src_iv = 1;
        det_iv = 1;
        for ivertex=1:length([src_hvidx ; det_hvidx])
            if ivertex <= nb_sources
                src_fluences{src_iv}{iwl} = mri_zeros;
                src_fluences{src_iv}{iwl}(:) = all_fluences_flat_sparse{ivertex}{iwl};
                src_iv = src_iv + 1;
            else
                det_fluences{det_iv}{iwl} = mri_zeros;
                det_fluences{det_iv}{iwl}(:) = all_fluences_flat_sparse{ivertex}{iwl};
                det_iv = det_iv + 1;
            end
            bst_progress('inc',1);
        end
    end
    
    det_reference_voxels_index = all_reference_voxels_index((nb_sources+1):(nb_sources+nb_dets));
    
    % Compute sensitivity matrix (volume) and interpolate on cortical surface
    % using Voronoi partitionning
    
    bst_progress('start', 'Sensitivity computation','Interpolate sensitivities on cortex...', 1, nb_pairs);
    
    nb_nodes            = size(sCortex.Vertices, 1);
    sensitivity_surf    = zeros(nb_nodes, nb_pairs, nb_wavelengths);
    
    if length(src_ids) == 1
        pair_ids  = [ src_ids(pair_sd_idx(:, 1))]; 
    else 
        pair_ids  = [ src_ids(pair_sd_idx(:, 1))]'; 
    end    
    if length(det_ids) == 1
        pair_ids  = [pair_ids, det_ids(pair_sd_idx(:, 2))];
    else 
        pair_ids  = [pair_ids, det_ids(pair_sd_idx(:, 2))'];
    end    
    separations_by_pairs = process_nst_separations('Compute', ChannelMat.Channel, pair_ids);
    
    
    for ipair  = 1:nb_pairs

        isrc = pair_sd_idx(ipair, 1);
        idet = pair_sd_idx(ipair, 2);
        separation = separations_by_pairs(ipair);

        for iwl=1:nb_wavelengths
    
            ref_fluence = src_fluences{isrc}{iwl}(det_reference_voxels_index{idet}{iwl}(1),...
                                                  det_reference_voxels_index{idet}{iwl}(2),...
                                                  det_reference_voxels_index{idet}{iwl}(3));
    
            
            separation_threshold = 0.055; % Below which fluence normalization fixing is allowed
            if ref_fluence == 0 && separation > separation_threshold
                sensitivity_vol = mri_zeros;
            else 
                if ref_fluence==0
                    normalization_factor = min(src_fluences{isrc}{iwl}(src_fluences{isrc}{iwl}>0));  
                    msg = sprintf('Fluence of S%02d is null at position of D%02d (wavelength=%dnm, separation=%1.2fcm).\n Using default normalization.', ...
                                   src_ids(isrc), det_ids(idet), ChannelMat.Nirs.Wavelengths(iwl), separation*100);
                    bst_report('Warning', 'process_nst_import_head_model', sInputs, msg);
                else
                    normalization_factor = ref_fluence;
                end
                sensitivity_vol = src_fluences{isrc}{iwl} .* ...
                                  det_fluences{idet}{iwl}./normalization_factor ;
            end

            % modified by zhengchen to normalize the sensitivity
            %sensitivity_vol = sensitivity_vol./max(sensitivity_vol(:)); 
            
            sens_tmp = accumarray(voronoi(voronoi_mask), sensitivity_vol(voronoi_mask), [nb_nodes+1 1],@(x)sum(x)/numel(x)); 
            sens_tmp(end)=[]; % trash last column
    
    
            sensitivity_surf(:,ipair,iwl) = sens_tmp;
        end
        bst_progress('inc',1);
    end
    
    [sensitivity_surf, warmInfo] = smooth_sensitivity_map(OPTIONS.CortexFile, sensitivity_surf, OPTIONS.smoothing_method, OPTIONS.smoothing_fwhm);
    if ~isempty(warmInfo)
        bst_report('Warning', 'process_nst_import_head_model', sInputs, warmInfo);
    end
    
    sensitivity_surf = permute(sensitivity_surf,[2,3,1]);
    
    % we finally format as brainstorm format 
    sensitivity_surf_bst = zeros(length(sChannels), nb_nodes*3);
    
    for iWavelenth = 1:nb_wavelengths
    
        gain_in         = sensitivity_surf;
        gain_pair_names = pair_names;
    
        % hash gain pair names
        gpair_idx = containers.Map();
        for gpic=1:length(gain_pair_names)
            gpair_idx(gain_pair_names{gpic}) = gpic;
        end
    
        selected_chans = channel_find(sChannels, sprintf('WL%d',ChannelMat.Nirs.Wavelengths(iWavelenth)));
    
        gains = zeros(length(selected_chans), size(gain_in, 3));
        for ic=1:length(selected_chans)
            ichan = selected_chans(ic);
            chan_name = sChannels(ichan).Name;
            pair_name = chan_name(1:strfind(chan_name, 'WL')-1);
            if ~isempty(pair_name) && gpair_idx.isKey(pair_name)
                gains(ic, :) = squeeze(gain_in(gpair_idx(pair_name), iwl, :));
            end
        end
    
        sensitivity_surf_bst(selected_chans, 1:3:end) = gains;
        sensitivity_surf_bst(selected_chans, 2:3:end) = gains;
        sensitivity_surf_bst(selected_chans, 3:3:end) = gains;
    
    end
    
    bst_progress('stop');
    
    
    %% Outputs
    % Save the new head model
    
    % Create structure
    HeadModelMat                = db_template('headmodelmat');
    HeadModelMat.NIRSMethod     = 'MCXlab';
    HeadModelMat.Gain           = sensitivity_surf_bst;
    HeadModelMat.HeadModelType  = 'surface';
    HeadModelMat.SurfaceFile    = OPTIONS.CortexFile;
    HeadModelMat.Comment        = 'NIRS head model';
    HeadModelMat.Param          = struct('FluenceFolder',    OPTIONS.FluenceFolder , ...
                                         'smoothing_method', OPTIONS.smoothing_method, ...
                                         'smoothing_fwhm',   OPTIONS.smoothing_fwhm);
    
    HeadModelMat                = bst_history('add', HeadModelMat, 'compute', 'Compute NIRS head model from MCX fluence results');

end

function sensitivity = get_sensitivity_from_chans(head_model, pair_names)

    head_model_pair_ids = containers.Map();
    for ipair=1:length(head_model.pair_names)
        head_model_pair_ids(head_model.pair_names{ipair}) = ipair;
    end

    pairs_not_found = {};
    for ipair=1:length(pair_names)
        pair_name = pair_names{ipair};
        if ~head_model_pair_ids.isKey(pair_name)
            pairs_not_found{end+1} = pair_name;
        end
    end
    if ~isempty(pairs_not_found)
        throw(MException('NIRSTORM:HeadmodelMismatch', ....
              ['Sensitivity not found for pairs: ', ...
               strjoin(pairs_not_found, ', ')]));
    end

    head_model_size = size(head_model.Gain);
    sensitivity = zeros(length(pair_names), head_model_size(2), head_model_size(3));
    for ipair=1:length(pair_names)
        ipair_head_model = head_model_pair_ids(pair_names{ipair});
        sensitivity(ipair, :, :) = head_model.Gain(ipair_head_model, :, :);
    end

end


function [src_head_vertex_ids, det_head_vertex_ids] = get_head_vertices_closest_to_optodes(sMri, sHead, src_locs, det_locs)

    head_vertices_mri = cs_convert(sMri, 'scs', 'mri', sHead.Vertices) * 1000;
    src_locs_mri = cs_convert(sMri, 'scs', 'mri', src_locs) * 1000;
    det_locs_mri = cs_convert(sMri, 'scs', 'mri', det_locs) * 1000;
    src_head_vertex_ids = nst_knnsearch(head_vertices_mri, src_locs_mri);
    det_head_vertex_ids = nst_knnsearch(head_vertices_mri, det_locs_mri);

end

function [fluences, reference] = request_fluences(head_vertices, anat_name, wavelengths, data_source, sparse_threshold, voi_mask, cube_size, local_cache_dir, use_closest_wl, sInput)

if nargin < 5
    sparse_threshold = nan;
end
 
if nargin < 6
    voi_mask = nan;
end

if nargin < 7
    cube_size = [];
end

if nargin < 8 || isempty(local_cache_dir)
    local_cache_dir = bst_fullfile(nst_get_local_user_dir(), ...
                                   'fluence', nst_protect_fn_str(anat_name));
end
%TODO: assert local cache directory exists

if nargin < 9 
   use_closest_wl = 0; 
end

if nargin < 10
    sInput = [];
end

fluence_fns = {};
fluences = {};
reference = {};
if ~isempty(strfind(data_source, 'http'))
    
    % Checking if URL repository can be found
    jurl = java.net.URL(data_source);
    conn = openConnection(jurl);
    status = getResponseCode(conn);
    if status == 404
        default = [nst_get_repository_url(), '/fluence/'];
        msg = sprintf('Fluence repository not found at %s.\n Switching to default: %s', data_source, default);
        bst_report('Warning', 'process_nst_import_head_model', sInput, msg);
    end
    
    if ~fluence_is_available(anat_name)
        bst_error(['Precomputed fluence data not available for anatomy "' anat_name '"']);
        return;
    end
    if ~strcmp(data_source(end), '/')
        data_source = [data_source '/'];
    end
    if ~exist(local_cache_dir, 'dir')
        mkdir(local_cache_dir);
    end
    
    if use_closest_wl
        wavelengths_for_dl = get_template_closest_wl(wavelengths);
        diff = (wavelengths_for_dl ~= wavelengths);
        if any(diff)
            ndiffs = sum(diff);
            diff_idx = find(diff);
            repl = '';
            for ii=1:ndiffs
                iwl = diff_idx(ii);
                repl = [repl sprintf('  - %dnm -> %dnm\n', wavelengths(iwl), wavelengths_for_dl(iwl))];
            end
            msg = sprintf('Some wavelengths do not have precomputed fluence data.\nReplacing them with the closest available:\n%s', repl);
            bst_report('Warning', 'process_nst_import_head_model', sInput, msg);
        end
    else
        wavelengths_for_dl = wavelengths;
    end
    
    to_download_urls = {};
    to_download_spec = {};
    dest_fns = {};
    default_fluence_file_size = 1000000; %bytes
    idownload = 1;
    total_download_size = 0;
    nb_files_to_check = length(head_vertices)*length(wavelengths);
    bst_progress('start', 'Retrieving server info', sprintf('Checking required downloads from server (%d files)...', nb_files_to_check), 1, nb_files_to_check);
    check_url_existence = 1;
    for ivertex=1:length(head_vertices)
        vertex_id = head_vertices(ivertex);
        for iwl=1:length(wavelengths)
            wl = wavelengths(iwl);
            fluence_bfn =  get_fluence_fn(vertex_id, wavelengths_for_dl(iwl));
            fluence_fn = fullfile(local_cache_dir, fluence_bfn);
            fluence_fns{ivertex}{iwl} = fluence_fn;
            
            if ~file_exist(fluence_fn)
                url = [data_source nst_protect_fn_str(anat_name) '/' fluence_bfn];
                to_download_urls{idownload} = url;
                if check_url_existence
                    tstart = tic();
                    jurl = java.net.URL(url);
                    conn = openConnection(jurl);
                    status = getResponseCode(conn);
                    if status == 404
                        bst_error(['Fluence file not available at ' url]);
                        return;
                    end
                    query_duration = toc(tstart);
                    if query_duration > 0.15
                        fprintf('Quit checking fluence file existence (too much time: %1.2f s / file).\n', query_duration);
                        check_url_existence = 0;
                    end
                end
                to_download_spec{idownload} = [vertex_id, wl];
                dest_fns{idownload} = fluence_fn;
                idownload = idownload + 1;
                
                download_size = default_fluence_file_size;
                total_download_size = total_download_size + download_size;
            end
            bst_progress('inc',1);
        end
    end
    bst_progress('stop');
    if total_download_size > 0
        % Ask user confirmation
        if ~java_dialog('confirm', ['Warning: ' format_file_size(total_download_size) ...
                ' of fluence data will be downloaded to ' local_cache_dir '.' 10 10 ...
                'Confirm download?' 10 10], 'Download warning')
            return;
        end
        % Process downloads
        msg = sprintf('Downloading %d fluence files (%s)...', ...
                      length(to_download_urls), ...
                      format_file_size(total_download_size));
        bst_progress('start', 'Downloading fluences', msg, 1, length(to_download_urls));
        for idownload=1:length(to_download_urls)
            vertex_id = to_download_spec{idownload}(1);
            wl = to_download_spec{idownload}(2);
            download_msg = ['Download fluence for v' num2str(vertex_id) ...
                            ', ' num2str(wl) 'nm'];
            if ~nst_download(to_download_urls{idownload}, dest_fns{idownload}, download_msg)
                return;
            end
            bst_progress('inc',1);
        end
        bst_progress('stop');
    end
else
    misssing_fluences = {}; 
    scout_missing = db_template('Scout'); 
    scout_missing.Label='Missing';
    
    for ivertex=1:length(head_vertices)
        vertex_id = head_vertices(ivertex);
        for iwl=1:length(wavelengths)
            wl = wavelengths(iwl);
            fluence_bfn =  get_fluence_fn(vertex_id, wl);
            fluence_fn = fullfile(data_source, fluence_bfn);
            
            if ~exist(fluence_fn, 'file')
                misssing_fluences{end+1}= [vertex_id wl];
                scout_missing.Vertices(end+1)= vertex_id;
            end
            fluence_fns{ivertex}{iwl} = fluence_fn;
        end
    end
    
    if ~isempty(misssing_fluences)
        scout_missing.Seed = scout_missing.Vertices(1);
        for k = 1:length(misssing_fluences)
            
            tmp = misssing_fluences{k};
            vertex_id = tmp(1); 
            wl = tmp(2);
            fluence_bfn =  get_fluence_fn(vertex_id, wl);
            fluence_fn = fullfile(data_source, fluence_bfn);
            disp(['Fluence file not found for v' num2str(vertex_id) ...
                           ', ' num2str(wl) 'nm (' fluence_fn ')']);
                
        end 
        bst_error('Missing fluences, see comand windows');
        return;
    end           
end

fluences = cell(length(head_vertices), 1);
if any(isnan(voi_mask))
    bst_progress('start', 'Get fluences', sprintf('Load volumic fluences (%d files)...', ...
                 length(head_vertices)*length(wavelengths)), 1, length(head_vertices)*length(wavelengths));
    reference_voxels_index = cell(length(head_vertices), 1);
    for ivertex=1:length(head_vertices)
        for iwl=1:length(wavelengths)
            fluence = load(fluence_fns{ivertex}{iwl});
            reference_voxel_index = fluence.reference_voxel_index;
            fluence = fluence.fluence_flat_sparse_vol;
            %         if ~issparse(fluence) && ~isnan(sparse_threshold)
            %             fluence(fluence < sparse_threshold) = 0;
            %             %fluence = sparse(fluence(:));
            %         end
            fluences{ivertex}{iwl} = fluence;
            reference_voxels_index{ivertex}{iwl} = reference_voxel_index;
            bst_progress('inc',1);
        end
    end
    reference = reference_voxels_index;
    bst_progress('stop');
else
    bst_progress('start', 'Get fluences', 'Load fluences head mask...',...
                 1, length(head_vertices));
    ref_voxel_indexes = zeros(size(head_vertices));
    for ivertex=1:length(head_vertices)
        data = load(fluence_fns{ivertex}{1}, 'reference_voxel_index');
        reference_voxel_index = data.reference_voxel_index;
        ref_voxel_indexes(ivertex) = sub2ind(cube_size, reference_voxel_index(1), ...
                                             reference_voxel_index(2), ...
                                             reference_voxel_index(3));
        bst_progress('inc',1);
    end
    bst_progress('stop');
    bst_progress('start', 'Get fluences', ...
                 sprintf('Load & mask volumic fluences (%d files)...', length(head_vertices)*length(wavelengths)), ...
                 1, length(head_vertices)*length(wavelengths));
    ref_fluences = cell(length(head_vertices), 1);
    for ivertex=1:length(head_vertices)
        for iwl=1:length(wavelengths)
            data = load(fluence_fns{ivertex}{iwl}, 'fluence_flat_sparse_vol');
            fluences{ivertex}{iwl} = data.fluence_flat_sparse_vol(voi_mask);
            ref_fluences{ivertex}{iwl} = data.fluence_flat_sparse_vol(ref_voxel_indexes);
            bst_progress('inc',1);
        end
    end
    bst_progress('stop');
    reference = ref_fluences;
end

end

function fluence_fn = get_fluence_fn(vertex_id, wl)
   fluence_fn = sprintf('fluence_%dnm_v%06d.mat', wl, vertex_id);
end


function closest_wavelengths = get_template_closest_wl(wavelengths)

    template_wls = [685];
    closest_wavelengths = template_wls(arrayfun(@(wl) iclosest(template_wls, wl), wavelengths));

end

function ic = iclosest(catalog, value)
    [v,ic] = min(abs(catalog-value));
end

function str_size = format_file_size(size)
    if size < 1e3
        str_size = sprintf('%do', size);
    elseif size < 1e6
        str_size = sprintf('%1.1fko', size/1e3);
    elseif size < 1e9
        str_size = sprintf('%1.1fMo', size/1e6);
    else
        str_size = sprintf('%1.2fGo', size/1e9);
    end
end

function flag = fluence_is_available(anat_name)
    flag = any(strcmp(strtrim(anat_name), {'MRI: Colin27 4NIRS'}));
end


%% ===== surface smoothing  =====
function [sensitivity_surf, warmInfo] = smooth_sensitivity_map(surface, sensitivity_surf, method, smoothing_fwhm)
    warmInfo = [];

    if smoothing_fwhm > 0
        FWHM = smoothing_fwhm / 1000;
        
        if  strcmp(method, 'surfstat_before_2023')
            [sensitivity_surf, msgInfo, warmInfo] = process_ssmooth_surfstat('compute',  surface,  sensitivity_surf, FWHM, 'before_2023');
        elseif strcmp(method, 'geodesic_dist')
            [sensitivity_surf, msgInfo, warmInfo] = process_ssmooth('compute',  surface, sensitivity_surf, FWHM, 'geodesic_dist');
        else
            [sensitivity_surf, msgInfo, warmInfo] = process_ssmooth_surfstat('compute',  surface,  sensitivity_surf, FWHM, 'before_2023');
        end

    end

end