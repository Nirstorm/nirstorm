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
    
    sProcess.options.use_closest_wl.Comment = 'Use closest available wavelength';
    sProcess.options.use_closest_wl.Type    = 'checkbox';
    sProcess.options.use_closest_wl.Value   = 0;


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
    
    sProcess.options.label2.Comment = '<B>Extra options:</B>';
    sProcess.options.label2.Type    = 'label';

    sProcess.options.use_all_pairs.Comment = 'Use all possible pairs: ';
    sProcess.options.use_all_pairs.Type    = 'checkbox';
    sProcess.options.use_all_pairs.Value   = 0;

    sProcess.options.normalize_fluence.Comment = 'Normalize by source fluence at detector position';
    sProcess.options.normalize_fluence.Type    = 'checkbox';
    sProcess.options.normalize_fluence.Value   = 1;

    sProcess.options.force_median_spread.Comment = 'Force median spread<FONT color="#777777">(not recommended)</FONT>';
    sProcess.options.force_median_spread.Type    = 'checkbox';
    sProcess.options.force_median_spread.Value   = 0; 
    
    sProcess.options.sensitivity_threshold_pct.Comment = 'Threshold (% of max-min): ';
    sProcess.options.sensitivity_threshold_pct.Type    = 'value';
    sProcess.options.sensitivity_threshold_pct.Value   = {0, '%', 2};
    
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

OutputFiles = {sInputs.FileName};

use_closest_wl      = sProcess.options.use_closest_wl.Value;
use_all_pairs       = sProcess.options.use_all_pairs.Value;
sens_thresh_pct     = sProcess.options.sensitivity_threshold_pct.Value{1};
normalize_fluence   = sProcess.options.normalize_fluence.Value;
data_source         = sProcess.options.data_source.Value;

ChannelMat = in_bst_channel(sInputs(1).ChannelFile);
if ~isfield(ChannelMat.Nirs, 'Wavelengths')
    bst_error(['Head model importation works only for dOD data ' ... 
               ' (eg do not use MBLL prior to this process)']);
    return;
end


%% Retrieve list of vertex indexes corresponding to current montage
% Load channel file
ChannelMat = in_bst_channel(sInputs(1).ChannelFile);

% Load subject
[sSubject, iSubject] = bst_get('Subject', sInputs.SubjectName);

% Load MRI, Head and Cortex. 
if isempty(sSubject.iAnatomy)
    bst_error(['No anatomical data found for ' sInputs.SubjectName]);
end

sMri    = in_mri_bst (sSubject.Anatomy(sSubject.iAnatomy).FileName);
sHead   = in_tess_bst(sSubject.Surface(sSubject.iScalp  ).FileName);
sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex ).FileName);

% Retrieve optode coordinates
% Load channel file

montage_info    = nst_montage_info_from_bst_channels(ChannelMat.Channel);
pair_names      = montage_info.pair_names;
src_locs        = montage_info.src_pos;
src_ids         = montage_info.src_ids;
det_locs        = montage_info.det_pos;
det_ids         = montage_info.det_ids;
pair_sd_idx     = montage_info.pair_sd_indexes;

nb_wavelengths = length(ChannelMat.Nirs.Wavelengths);
nb_sources = size(src_locs, 1);
nb_dets = size(det_locs, 1);
if use_all_pairs
    [gs, gd] = meshgrid(1:nb_sources, 1:nb_dets);
    pair_sd_idx = [gs(:) gd(:)];
    pair_names = {};
    for ipair = 1:size(pair_sd_idx, 1)
        pair_names{ipair} = nst_format_channel(src_ids(pair_sd_idx(ipair, 1)),...
                                               det_ids(pair_sd_idx(ipair, 2)));
    end
end

nb_pairs = length(pair_names);
% Find closest head vertices (for which we have fluence data)
% Put everything in mri referential

[src_hvidx, det_hvidx] = get_head_vertices_closest_to_optodes(sMri, sHead, src_locs, det_locs);

%% Load fluence data from local .brainstorm folder (download if not available)
        
[all_fluences_flat_sparse, all_reference_voxels_index]= request_fluences([src_hvidx ; det_hvidx], sMri.Comment, ...
                                                                         ChannelMat.Nirs.Wavelengths, data_source, nan, nan, [], '',...
                                                                         use_closest_wl);
if isempty(all_fluences_flat_sparse)
    return;
end

bst_progress('start', 'Export fluences','Unpacking volumic fluences...', 1, (nb_sources+nb_dets)*nb_wavelengths);

mri_zeros       = zeros(size(sMri.Cube));
src_fluences    = cell(1, nb_sources);
det_fluences    = cell(1, nb_dets);

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
sensitivity_surf    = zeros(nb_nodes,nb_pairs, nb_wavelengths);

voronoi = get_voronoi(sProcess, sInputs);
voronoi_mask = (voronoi > -1) & ~isnan(voronoi);



% separations_chans = process_nst_separations('Compute', ChannelMat.Channel);

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
separations_by_pairs = process_nst_separations('Compute', ChannelMat.Channel,pair_ids);


for ipair=1:nb_pairs
    isrc = pair_sd_idx(ipair, 1);
    idet = pair_sd_idx(ipair, 2);
    separation = separations_by_pairs(ipair);
    for iwl=1:nb_wavelengths
        if normalize_fluence
            ref_fluence = src_fluences{isrc}{iwl}(det_reference_voxels_index{idet}{iwl}(1),...
                                                  det_reference_voxels_index{idet}{iwl}(2),...
                                                  det_reference_voxels_index{idet}{iwl}(3));
        else
            ref_fluence = 1;
        end
        
        separation_threshold = 0.055; % Below which fluence normalization fixing is allowed
        if ref_fluence==0 && separation > separation_threshold
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

if sProcess.options.smoothing_fwhm.Value{1} > 0
    FWHM = sProcess.options.smoothing_fwhm.Value{1} / 1000;

    if ~isfield(sProcess.options, 'method') || strcmp(sProcess.options.method.Value, 'surfstat_before_2023')
        [sensitivity_surf, msgInfo, warmInfo] = process_ssmooth_surfstat('compute', ... 
                                        sSubject.Surface(sSubject.iCortex).FileName, ...
                                        sensitivity_surf, FWHM, 'before_2023');
    elseif strcmp(sProcess.options.method.Value, 'geodesic_dist')
        [sensitivity_surf, msgInfo, warmInfo] = process_ssmooth('compute', ... 
                                        sSubject.Surface(sSubject.iCortex).FileName, ...
                                        sensitivity_surf, FWHM, 'geodesic_dist');
         bst_report('Warning', 'process_nst_import_head_model', sInputs, warmInfo);
    end
end

sensitivity_surf = permute(sensitivity_surf,[2,3,1]);

bst_progress('stop');

%% Sensitivity thresholding
if sens_thresh_pct > 0

    nz_sensitivity = sensitivity_surf(sensitivity_surf>0);
    min_sens =  min(nz_sensitivity);
    max_sens =  max(nz_sensitivity);
    sens_thresh = min_sens + sens_thresh_pct/100 * (max_sens - min_sens);

    sensitivity_surf(sensitivity_surf<sens_thresh) = 0;
end

if sProcess.options.force_median_spread.Value
    if ~isfield(sCortex, 'VertArea') || isempty(sCortex.VertArea)
        [tmp, sCortex.VertArea] = tess_area(sCortex.Vertices, sCortex.Faces);
    end
    
    areas = zeros(nb_pairs, nb_wavelengths);
    for ipair=1:nb_pairs
        for iwl=1:nb_wavelengths
            areas(ipair, iwl) = sum(sCortex.VertArea(squeeze(sensitivity_surf(ipair, iwl, :) > 0)));
        end
    end
    median_area = median(areas(areas>0));
    
    for ipair=1:nb_pairs
        for iwl=1:nb_wavelengths
            area = areas(ipair, iwl);
            if area > 0
                vi_orig = find(sensitivity_surf(ipair, iwl, :) > 0);
                if area > median_area
                    [sorted_ss, vertex_stack] = sort(squeeze(sensitivity_surf(ipair, iwl, :)));
                    vertex_stack = intersect(vertex_stack, vi_orig, 'stable');
                    iv = 1;
                    while sum(sCortex.VertArea(vertex_stack(iv:end))) > median_area
                        iv = iv + 1;
                    end
                    sensitivity_surf(ipair, iwl, vertex_stack(1:iv)) = 0;
                    
                else
                    seedXYZ = mean(sCortex.Vertices(vi_orig, :));
                    to_add = [];
                    vi_growth = vi_orig;
                    while area < median_area
                        % Grow one vertex at a time
                        % Get closest neighbours
                        viNew = setdiff(tess_scout_swell(vi_growth, sCortex.VertConn), vi_growth);
                        if ~isempty(viNew)
                            % Compute the distance from each point to the seed
                            distFromSeed = sqrt(sum(bst_bsxfun(@minus, sCortex.Vertices(viNew,:), seedXYZ) .^ 2, 2));
                            % Get the minimum distance
                            [minVal, iMin] = min(distFromSeed);
                            iMin = iMin(1);
                            % Add the closest vertex to scout vertices
                            vi_growth = union(vi_growth, viNew(iMin));
                            area = sum(sCortex.VertArea(vi_growth));
                            to_add(end+1) = viNew(iMin); %#ok<AGROW>
                        else
                            warning('Sensitivity cluster growth stopped before reaching target area');
                            break;
                        end
                    end
                    sensitivity_surf(ipair, iwl, to_add) = min(sensitivity_surf(ipair, iwl, vi_orig));
                end
            end
        end
    end
end


%% Outputs
% Save the new head model
sStudy = bst_get('Study', sInputs.iStudy);

% Create structure
HeadModelMat = db_template('headmodelmat');
HeadModelMat.Gain           = sensitivity_surf;
HeadModelMat.HeadModelType  = 'surface';
HeadModelMat.SurfaceFile    = sSubject.Surface(sSubject.iCortex).FileName;
HeadModelMat.Comment        = 'NIRS head model';
if use_all_pairs
    HeadModelMat.Comment = [HeadModelMat.Comment ' [all pairs]'];
end
% newHeadModelMat.VoiNodes = voi_nodes;
HeadModelMat.pair_names = pair_names;
HeadModelMat  = bst_history('add', HeadModelMat, 'compute', 'Compute NIRS head model from MCX fluence results');
% Output file name
HeadModelFile = bst_fullfile(bst_fileparts(file_fullpath(sStudy.FileName)), 'headmodel_nirs_mcx_fluence.mat');
HeadModelFile = file_unique(HeadModelFile);
% Save file
bst_save(HeadModelFile, HeadModelMat, 'v7');

newHeadModel = db_template('HeadModel');
newHeadModel.FileName = file_short(HeadModelFile);
newHeadModel.Comment = 'NIRS head model from fluence';
if use_all_pairs
    newHeadModel.Comment = [newHeadModel.Comment ' [all pairs]'];
end
newHeadModel.HeadModelType  = 'surface';    
% Update Study structure
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

% Save database
db_save();
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
                                   'fluence', protect_fn_str(anat_name));
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
                url = [data_source protect_fn_str(anat_name) '/' fluence_bfn];
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
% TODO: check online
flag = any(strcmp(strtrim(anat_name), {'MRI: Colin27 4NIRS'}));
end



function voronoi = get_voronoi(sProcess, sInputs)
[sSubject, iSubject] = bst_get('Subject', sInputs.SubjectName);
voronoi_fn = process_nst_compute_voronoi('get_voronoi_fn', sSubject);

if ~exist(voronoi_fn, 'file')
    sProcess.options.do_grey_mask.Value   = 1; 
    process_nst_compute_voronoi('Run', sProcess, sInputs);
end
voronoi_bst = in_mri_bst(voronoi_fn);
voronoi = voronoi_bst.Cube;
end

function fn = protect_fn_str(sfn)
fn = strrep(sfn, ' ', '_');
fn = strrep(fn, '"', '');
fn = strrep(fn, ':', '_');
fn = strrep(fn, '(', '_');
fn = strrep(fn, ')', '_');
end

%% ===== surface smoothing  =====
function sens_smoothed = surface_smooth(FWHM, SurfaceMat, sens_temp, dispInfo) 
    % ===== PROCESS =====
    % Convert surface to SurfStat format
    cortS.tri = SurfaceMat.Faces;
    cortS.coord = SurfaceMat.Vertices';

    % Get the average edge length
    [vi,vj] = find(SurfaceMat.VertConn);
    Vertices = SurfaceMat.VertConn;
    meanDist = mean(sqrt((Vertices(vi,1) - Vertices(vj,1)).^2 + (Vertices(vi,2) - Vertices(vj,2)).^2 + (Vertices(vi,3) - Vertices(vj,3)).^2));
    % FWHM in surfstat is in mesh units: Convert from millimeters to "edges"
    FWHMedge = FWHM ./ meanDist;
    
    % Display the result of this conversion
    msgInfo = ['Average distance between two vertices: ' num2str(round(meanDist*10000)/10) ' mm' 10 ...
               'SurfStatSmooth called with FWHM=' num2str(round(FWHMedge * 1000)/1000) ' edges'];
    
    %bst_report('Info', sProcess, sInput, msgInfo);
    if dispInfo
        disp(msgInfo);
        disp(['SMOOTH> ' strrep(msgInfo, char(10), [10 'SMOOTH> '])]); 
    end
      
    if round(FWHMedge * 1000)/1000 ==0
        disp('WARNING: FWHM too small, smoothing will not be performed.')
    end
    % Smooth surface
    
    sens_smoothed = SurfStatSmooth(sens_temp', cortS, FWHMedge)';
    
    
    
end
