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
%          Edouard Delaire (2025)

eval(macro_method);
end

function sProcess = GetDescription() 
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
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput) 

    OutputFiles = {sInput.FileName};
    
    
    ChannelMat = in_bst_channel(sInput.ChannelFile);
    if ~isfield(ChannelMat.Nirs, 'Wavelengths')
        bst_error('Head model importation works only for dOD data (eg do not use MBLL prior to this process)');
        return;
    end

    % Load subject
    sSubject = bst_get('Subject', sInput.SubjectName);
    voronoi_fn = process_nst_compute_voronoi('get_voronoi_fn', sSubject);

    %  Check that all the information are present
    if isempty(sSubject.iAnatomy) ||  isempty(sSubject.iCortex) || isempty(sSubject.iScalp) || ~exist(voronoi_fn,"file")
        bst_error(sprintf('Anatomical data for subject %s are incomplete ', sInput.SubjectName));
    end

    % Load cortical surface.
    sCortex  = in_tess_bst(sSubject.Surface(sSubject.iCortex ).FileName);

    % Load channels
    ChannelMat  = in_bst_channel(sInput.ChannelFile);

    OPTIONS = struct();
    
    % Subject Informations 
    OPTIONS.MriFile             = sSubject.Anatomy(sSubject.iAnatomy).FileName;
    OPTIONS.VoronoiFile         = voronoi_fn;
    OPTIONS.HeadFile            = sSubject.Surface(sSubject.iScalp  ).FileName;
    OPTIONS.CortexFile          = sSubject.Surface(sSubject.iCortex ).FileName;
    OPTIONS.Channel             = ChannelMat.Channel;

    % Use defined options : 
    OPTIONS.FluenceFolder       = sProcess.options.data_source.Value;
    OPTIONS.smoothing_method    = sProcess.options.method.Value;
    OPTIONS.smoothing_fwhm      = sProcess.options.smoothing_fwhm.Value{1};
    
    % Compute Head model
    [Gain, error_message, warning_message]  = Compute(OPTIONS);

    if  ~isempty(error_message)
        bst_report('Error',   sProcess, sInput, error_message);
        return;
    end
    
    if  ~isempty(warning_message)
        for iMessage = 1:length(warning_message)
            bst_report('Warning',   sProcess, sInput, warning_message{iMessage});
        end
    end
    

    % Outputs
    HeadModelMat                = db_template('headmodelmat');
    HeadModelMat.NIRSMethod     = 'MCXlab';
    HeadModelMat.Gain           = Gain;
    HeadModelMat.HeadModelType  = 'surface';
    HeadModelMat.SurfaceFile    = OPTIONS.CortexFile;
    HeadModelMat.GridOrient     = sCortex.VertNormals;
    HeadModelMat.Comment        = 'NIRS head model';
    HeadModelMat.Param          = struct('FluenceFolder',    OPTIONS.FluenceFolder , ...
                                         'smoothing_method', OPTIONS.smoothing_method, ...
                                         'smoothing_fwhm',   OPTIONS.smoothing_fwhm);
    
    HeadModelMat                = bst_history('add', HeadModelMat, 'compute', 'Compute NIRS head model from MCX fluence results');


    % Save Head Model
    sStudy = bst_get('Study', sInput.iStudy);
    
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
    bst_set('Study', sInput.iStudy, sStudy);
    panel_protocols('UpdateNode', 'Study', sInput.iStudy);
    db_save();

end

function [Gain, error_message, warning_message] = Compute(OPTIONS)
% Compute the nirs head model.
%
% Input: 
%
%  OPTIONS: Struct with the following fields
%  |--  OPTIONS.MriFile             : path to MRI
%  |--  OPTIONS.VoronoiFile         : path to voronoi volume to surface interpolator
%  |--  OPTIONS.HeadFile            : path to scalp surface
%  |--  OPTIONS.CortexFile          : path to cortex surface
%  |--  OPTIONS.Channel             : full brainstorm channel structure
%
%  |--  OPTIONS.FluenceFolder       ; path to the fluence folder. Can be a server
%  |--  OPTIONS.smoothing_method    : method used to smooth the fluence {'surfstat_before_2023' or  'geodesic_dist' (recomended) }
%  |--  OPTIONS.smoothing_fwhm      : parameter of the smoothing kernel
%
% Output: 
% -- HeadModelMat: Brainstorm structure for the head model
% -- error_message: Description of the error
% -- warning_message: Description of the warning

    Gain            = nan;
    error_message   = '';
    warning_message = {};

    % Retrieve optode coordinates
    % Load channel file
    sMri    = in_mri_bst (OPTIONS.MriFile);
    sHead   = in_tess_bst(OPTIONS.HeadFile);
    sCortex = in_tess_bst(OPTIONS.CortexFile);

    sChannels   = OPTIONS.Channel;
    
    iNIRS         = channel_find(sChannels, 'NIRS');
    sChannelsNIRS = sChannels(iNIRS);

    % find used wavelentgh 
    wavelengths = unique({sChannelsNIRS.Group});
    wavelengths = cellfun( @(x) str2double(x(3:end)), wavelengths) ;

    % Get Voronoi
    voronoi_bst = in_mri_bst(OPTIONS.VoronoiFile);
    voronoi     = voronoi_bst.Cube;
    voronoi_mask = (voronoi > -1) & ~isnan(voronoi);

    % Load montage informations
    montage_info    = nst_montage_info_from_bst_channels(sChannelsNIRS);
    src_locs        = montage_info.src_pos;
    src_ids         = montage_info.src_ids;
    det_locs        = montage_info.det_pos;
    det_ids         = montage_info.det_ids;
    nb_sources      = size(src_locs, 1);
    nb_dets         = size(det_locs, 1);
    
    
    % Find closest head vertices (for which we have fluence data)    
    [src_hvidx, det_hvidx] = get_head_vertices_closest_to_optodes(sMri, sHead, src_locs, det_locs);
    
    % Load fluence data from local .brainstorm folder (download if not available)
    local_cache_dir = bst_fullfile(nst_get_local_user_dir(),  'fluence', nst_protect_fn_str(sMri.Comment));

    % Request Fluences files 
    [fluence_fns, missing_fluences] = request_fluences(     OPTIONS.FluenceFolder  , ...
                                        [src_hvidx ; det_hvidx] , ...
                                        wavelengths, ...
                                        local_cache_dir);

    % If missing fluences, list them, and return.
    if ~isempty(missing_fluences)
        error_message = list_missing_fluences(missing_fluences);
        return;
    end

    % Load Fluences
    [all_fluences, all_reference_voxels_index] = load_fluences(fluence_fns, size(sMri.Cube));

    if isempty(all_fluences)
        error_message = ' Some fluence could not be loaded';
        return;
    end

    % Unpack the fluences, and separate sources and detectors.
    src_fluences    = all_fluences(1:nb_sources);
    det_fluences    = all_fluences((nb_sources+1):(nb_sources+nb_dets));
    det_reference_voxels_index = all_reference_voxels_index((nb_sources+1):(nb_sources+nb_dets));

    separations = process_nst_separations('Compute', sChannelsNIRS);

    % Compute sensitivity matrix (volume) and interpolate on cortical surface
    % using Voronoi partitionning
    bst_progress('start', 'Sensitivity computation','Interpolate sensitivities on cortex...', 1, length(sChannelsNIRS));
    nb_nodes             = size(sCortex.Vertices, 1);
    sensitivity_surf = zeros(nb_nodes, length(sChannels));

    for iChannel = 1:length(sChannelsNIRS)

        [source_id, det_id, measure_id] = nst_unformat_channel(sChannelsNIRS(iChannel).Name);
        
        isrc    = find( src_ids ==  source_id);
        idet    = find( det_ids ==  det_id);
        iwl     = find( wavelengths == measure_id);

        ref_fluence = src_fluences{isrc}{iwl}(det_reference_voxels_index{idet}{iwl}(1),...
                                              det_reference_voxels_index{idet}{iwl}(2),...
                                              det_reference_voxels_index{idet}{iwl}(3));

        separation_threshold = 0.055; % Below which fluence normalization fixing is allowed
        if ref_fluence == 0 && separations(iChannel) > separation_threshold
            sensitivity_vol = nan(size(sMri.Cube));
            warning_message{end+1} = sprintf('Fluence for channel cannot be computed S%02dD%02d (wavelength=%dnm, separation=%1.2fcm).\n', ...
                                                               source_id, det_id, ChannelMat.Nirs.Wavelengths(iwl), separations(iChannel)*100);
            
        else 
            if ref_fluence==0
                normalization_factor = min(src_fluences{isrc}{iwl}(src_fluences{isrc}{iwl}>0));  

                warning_message{end+1} = sprintf('Fluence of S%02d is null at position of D%02d (wavelength=%dnm, separation=%1.2fcm).\n Using default normalization.', ...
                                                                                source_id, det_id, ChannelMat.Nirs.Wavelengths(iwl), separations(iChannel)*100);
            else
                normalization_factor = ref_fluence;
            end
    
            sensitivity_vol = src_fluences{isrc}{iwl} .* ...
                              det_fluences{idet}{iwl} ./ normalization_factor ;
        end
        
        sens_tmp = accumarray(voronoi(voronoi_mask), sensitivity_vol(voronoi_mask), [nb_nodes+1 1],@(x)sum(x)/numel(x)); 
        sens_tmp(end)=[]; % trash last column
    
        sensitivity_surf(:, iNIRS(iChannel)) = sens_tmp;
        
        bst_progress('inc',1);
    end
        
    [sensitivity_surf, warmInfo] = smooth_sensitivity_map(OPTIONS.CortexFile, sensitivity_surf, OPTIONS.smoothing_method, OPTIONS.smoothing_fwhm);
    if ~isempty(warmInfo)
        warning_message{end+1} =  warmInfo;
    end
    
    % sensitivity_surf is now nChannel x nNodes
    sensitivity_surf = sensitivity_surf';

    % We convert to 3 orientation for saving in brainstorm database  
    Gain = normal_to_xyz(sensitivity_surf, sCortex.VertNormals);
  
end

function [src_head_vertex_ids, det_head_vertex_ids] = get_head_vertices_closest_to_optodes(sMri, sHead, src_locs, det_locs)

    head_vertices_mri = cs_convert(sMri, 'scs', 'mri', sHead.Vertices) * 1000;
    src_locs_mri = cs_convert(sMri, 'scs', 'mri', src_locs) * 1000;
    det_locs_mri = cs_convert(sMri, 'scs', 'mri', det_locs) * 1000;
    src_head_vertex_ids = nst_knnsearch(head_vertices_mri, src_locs_mri);
    det_head_vertex_ids = nst_knnsearch(head_vertices_mri, det_locs_mri);

end

function [fluence_fns, missing_fluences] = request_fluences(data_source, head_vertices, wavelengths, local_cache_dir)
    
    if nargin < 4
        local_cache_dir = '';
    end    

    if contains(data_source, 'http')

        if ~exist(local_cache_dir, 'dir')
            mkdir(local_cache_dir);
        end

        fluence_folder = local_cache_dir;
    else

        fluence_folder = data_source;
    end

    % List fluences Files
    assert(isfolder(fluence_folder));
    fluence_fns = list_fluences(fluence_folder, head_vertices, wavelengths);
    

    % Check for missing fluences
    flat_fluence_fns = [fluence_fns{:}]';
    missing_fluences = flat_fluence_fns(cellfun(@(x) exist(x, 'file') ~=2 , flat_fluence_fns));
    
    % Download missing fluences 
    if ~isempty(missing_fluences) && contains(data_source, 'http')

        download_fluences(data_source, local_cache_dir,  missing_fluences);
        missing_fluences = missing_fluences(cellfun(@(x) exist(x, 'file') ~=2 , missing_fluences));

    end

end

function download_fluences(data_source, local_cache_dir,  missing_fluences)

   % Checking if URL repository can be found
    jurl    = java.net.URL(data_source);
    conn    = openConnection(jurl);
    status  = getResponseCode(conn);

    if status == 404
        default = [nst_get_repository_url(), '/fluence/'];
        msg = sprintf('Fluence repository not found at %s.\n Switching to default: %s', data_source, default);
        bst_report('Warning', 'process_nst_import_head_model', sInput, msg);
    end
    
    % Check files existance on the server
    can_be_downloaded = true(1,length(missing_fluences));

    for iFluence = 1:length(missing_fluences)
        [~, fluence_name, ext] = fileparts(missing_fluences{iFluence});
        url = fullfile(data_source, [fluence_name,ext] );
        url = strrep(url, '\', '/');

        tstart = tic();


        jurl    = java.net.URL(url);
        conn    = openConnection(jurl);
        status  = getResponseCode(conn);

    
        if status == 404
            can_be_downloaded(iFluence) = false;
        end

        query_duration = toc(tstart);

        if query_duration > 0.15
            fprintf('Quit checking fluence file existence (too much time: %1.2f s / file).\n', query_duration);
            break;
        end
    end
    
    if any(~can_be_downloaded)
        disp('Some fluences are missing on the server : ')
        disp(list_missing_fluences(missing_fluences(~can_be_downloaded)))
    end

    missing_fluences = missing_fluences(can_be_downloaded);
    if isempty(missing_fluences)
        return;
    end

    default_fluence_file_size = 1000000; %bytes
    total_download_size = sum(can_be_downloaded) * default_fluence_file_size;

    % Ask user confirmation
    if ~java_dialog('confirm', ['Warning: ' format_file_size(total_download_size) ...
                                ' of fluence data will be downloaded to ' local_cache_dir '.' 10 10 ...
                                'Confirm download?' 10 10], 'Download warning')
        return;
    end

    % Process downloads
    msg = sprintf('Downloading %d fluence files (%s)...', ...
                  length(missing_fluences), ...
                  format_file_size(total_download_size));


    bst_progress('start', 'Downloading fluences', msg, 1, length(missing_fluences));
    for idownload=1:length(missing_fluences)
        
        [~, fluence_name, ext] = fileparts(missing_fluences{idownload});
        url = fullfile(data_source, [fluence_name,ext] );
        url = strrep(url, '\', '/');

        if ~nst_download(url, missing_fluences{idownload}, ['Downloading fluence '  missing_fluences{idownload}])
            return;
        end

        bst_progress('inc',1);
    end

    bst_progress('stop');
end

function fluence_fn = get_fluence_fn(vertex_id, wl)
    fluence_fn = sprintf('fluence_%dnm_v%06d.mat', wl, vertex_id);
end

function fluence_fns = list_fluences(data_source, head_vertices, wavelengths)
% Return the list of fluences files for the list of head vertices and
% wavelentgh. 


    fluence_fns = cell(1, length(head_vertices));

    for ivertex=1:length(head_vertices)
        vertex_id = head_vertices(ivertex);

        for iwl=1:length(wavelengths)

            fluence_bfn =  get_fluence_fn(vertex_id, wavelengths(iwl));
            fluence_fn  = fullfile(data_source, fluence_bfn);
            
            fluence_fns{ivertex}{iwl} = fluence_fn;
        end
    end
end

function warning_msg = list_missing_fluences(missing_fluences)
    
    warning_msg = ''; 
    tokens = cellfun( @(x) regexp(x, 'fluence_(\d+)nm_v(\d+)\.mat', 'tokens'), missing_fluences);
    
    % Flatten each inner token (each is a 1x1 cell containing a 1x2 cell)
    misssing_wavelength = cellfun(@(x) x{1}, tokens, 'UniformOutput', false);
    misssing_vertex = cellfun(@(x) x{2}, tokens, 'UniformOutput', false);

    [counts, list_missing_wavelength] = groupcounts(misssing_wavelength);

    for iWavelength = 1:length(list_missing_wavelength)

        warning_msg = [ warning_msg, ...
                        sprintf('%d fluences missing for wavelength %s : \n',counts(iWavelength), list_missing_wavelength{iWavelength})];
        idx_fluences = find(cellfun(@(x) strcmp(x,  list_missing_wavelength{iWavelength}), misssing_wavelength));

        for k = 1:length(idx_fluences)
            warning_msg = [ warning_msg, ...
                            sprintf('Vertex %s : %s \n', misssing_vertex{idx_fluences(k)},missing_fluences{idx_fluences(k)})];
        end
        warning_msg  = [warning_msg, newline];
    end

end

function [fluences, reference] = load_fluences(fluence_fns, cube_size)

    flat_fluence_fns = [fluence_fns{:}]';

    nFluences = length(flat_fluence_fns);
    nVertex    = length(fluence_fns);
    nWavelength = nFluences / nVertex; 

    assert(nWavelength == round(nWavelength))


    fluences    = repmat( {repmat({zeros(cube_size)}, 1, nWavelength)} ,  nVertex, 1);
    reference   = cell(nVertex, 1);


    bst_progress('start', 'Load fluences', sprintf('Loading %d fluences...',  nFluences), 1, nFluences);
    for ivertex=1:nVertex
        for iwl=1:nWavelength

            fluence = load(fluence_fns{ivertex}{iwl});
            reference_voxel_index = fluence.reference_voxel_index;
            fluence = fluence.fluence_flat_sparse_vol;


            fluences{ivertex}{iwl}(:)  = fluence;
            reference{ivertex}{iwl} = reference_voxel_index;
            bst_progress('inc',1);
        end
    end
    bst_progress('stop');
end


function [fluences, reference] = load_fluence_with_mask(fluence_fns, cube_size, mask)

    flat_fluence_fns = [fluence_fns{:}]';

    nFluences = length(flat_fluence_fns);
    nVertex    = length(fluence_fns);
    nWavelength = nFluences / nVertex; 

    assert(nWavelength == round(nWavelength))

    fluences    = cell(nVertex, 1);
    reference   = cell(nVertex, 1);

    bst_progress('start', 'Get fluences', 'Load fluences head mask...', 1, nVertex);

    ref_voxel_indexes = zeros(1, nVertex);

    for ivertex=1:nVertex

        data = load(fluence_fns{ivertex}{1}, 'reference_voxel_index');

        reference_voxel_index      = data.reference_voxel_index;

        ref_voxel_indexes(ivertex) = sub2ind(   cube_size, ...
                                                reference_voxel_index(1), ...
                                                reference_voxel_index(2), ...
                                                reference_voxel_index(3));
        bst_progress('inc',1);
    end
    bst_progress('stop');



    bst_progress('start', 'Get fluences', sprintf('Load & mask volumic fluences (%d files)...', nFluences), 1, nFluences);
    
    for ivertex=1:nVertex
        for iwl=1:nWavelength

            data = load(fluence_fns{ivertex}{iwl}, 'fluence_flat_sparse_vol');

            fluences{ivertex}{iwl}  = data.fluence_flat_sparse_vol(mask);
            reference{ivertex}{iwl} = data.fluence_flat_sparse_vol(ref_voxel_indexes);
            
            bst_progress('inc',1);
        end
    end
    
    bst_progress('stop');

end


function headmodel = convert_head_model( sChannel, headmodel, apply_orientation)   
% Convert the head model from old nirstorm structure to the brainstorm
% structure 

    if nargin < 3
        apply_orientation = 1;
    end

    gain = zeros(length(sChannel.Channel), size(headmodel.Gain, 3));

    iNIRS = good_channel(sChannel.Channel, [], 'NIRS');
    sChannelNirs = sChannel.Channel(iNIRS);

    for iChannel = 1:length(iNIRS)
        [source_id, det_id, measure_id] = nst_unformat_channel(sChannelNirs(iChannel).Name);
        

        iPair =  find(strcmp(headmodel.pair_names, sprintf('S%dD%d',source_id, det_id )));
        iWl   =  find(sChannel.Nirs.Wavelengths == measure_id );

        gain(iNIRS(iChannel), : ) = headmodel.Gain(iPair, iWl, :);
    end

    if apply_orientation
        
        % COnvert the nornmal vector to x,y,z 
        sCortex = in_tess_bst(headmodel.SurfaceFile);

        gain = normal_to_xyz(gain, sCortex.VertNormals);
        headmodel.GridOrient =  sCortex.VertNormals;

    else
        % Ensute to empty the Grid Orient so the orientation can't be
        % applied accidentally later. 

        headmodel.GridOrient = [];
    end
    
    % Update and complete the head model
    headmodel.Gain = gain;
    headmodel = rmfield(headmodel, 'pair_names');
    headmodel.NIRSMethod     = 'MCXlab';
    headmodel                = bst_history('add', headmodel, 'compute', 'Convert head model structure');

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
            [sensitivity_surf, ~, warmInfo] = process_ssmooth_surfstat('compute',  surface,  sensitivity_surf, FWHM, 'before_2023');
        elseif strcmp(method, 'geodesic_dist')
            [sensitivity_surf, ~, warmInfo] = process_ssmooth('compute',  surface, sensitivity_surf, FWHM, 'geodesic_dist');
        else
            [sensitivity_surf, ~, warmInfo] = process_ssmooth_surfstat('compute',  surface,  sensitivity_surf, FWHM, 'before_2023');
        end

    end

end

function A = normal_to_xyz(sensitivity_surf, VertNormals)
% Expend the sensitivity map to x,y,z coordinates 
% from the normal magnitude

    A  = zeros( size(sensitivity_surf, 1), 3 * size(sensitivity_surf, 2) );
    for k = 1:3
        A(:, k:3:end) = sensitivity_surf .* VertNormals(:, k)';
    end

end