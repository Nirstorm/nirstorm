function varargout = process_nst_OM( varargin )

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
% Authors: 
% Edouard Delaire (2021-2025), Jean-Eudes Bornert (2025)
% Thomas Vincent, Alexis Machado, ZhengChen Cai (2017)

eval(macro_method);
end

function sProcess = GetDescription() 

% Description of the process
sProcess.Comment     = 'Compute optimal montage';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = {'NIRS', 'Sources'};
sProcess.Index       = 1406;
sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/NIRS_Optimal_montage';
sProcess.InputTypes  = {'import'};
sProcess.OutputTypes = {'import'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 0;

sProcess.options.subjectname.Comment = 'Subject name:';
sProcess.options.subjectname.Type    = 'subjectname';
sProcess.options.subjectname.Value   = '';

sProcess.options.fluencesCond.Comment = {'panel_nst_OM', 'Optimal montage parameters'};
sProcess.options.fluencesCond.Type    = 'editpref';
sProcess.options.fluencesCond.Value   = [];

end

function s = str_pad(s,padsize)
    if (length(s) < padsize)
        s = [repmat('&nbsp;', 1, padsize - length(s)), s];
    end
    s = ['<FONT FACE="monospace">' s '</FONT>'];
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFile = Run(sProcess, sInput) 

    OutputFile = {};
    
    if bst_iscompiled()
        bst_error('Optimum montage is not available in the compiled version of brainstorm');
        return;
    end        
    cplex_url = 'https://www.ibm.com/us-en/marketplace/ibm-ilog-cplex/resources';
    
    if ~check_cplex(cplex_url)
        bst_error(['CPLEX >12.3 required. See ' cplex_url]);
    end

    options     = sProcess.options.fluencesCond.Value;
    if ~isfield(options, 'condition_name') || isempty(options.condition_name)
        options.condition_name = 'planning_optimal_montage';
    end
    
    SubjectName = options.SubjectName;
    sProcess.options.subjectname.Value = SubjectName;
    
    sSubject    = bst_get('Subject',SubjectName);
    if isempty(sSubject.iCortex) || isempty(sSubject.iScalp)
        bst_error('No available Cortex and Head surface for this subject.');
        return;
    end
    
    [status, error, options] = check_user_inputs(options);
    if ~status
        err_msg = sprintf("%d errors occured : \n%s", length(error), strjoin(" - " + error, '\n'));
        bst_error(err_msg);
        return;
    end
    
    [ROI_cortex, options.ROI_head] = get_regions_of_interest(sSubject, options);   

    % Obtain the anatomical MRI
    sMri     = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
    options.cubeSize        = size(sMri.Cube);
    
    options = get_weight_tables(sSubject, sProcess, sInput, options, ROI_cortex);
    if isempty(options.sensitivity_mat) || nnz(options.sensitivity_mat) == 0
        bst_error(sprintf('Weight table is null for ROI: %s', ROI_cortex.Label));
        return
    end

    if options.flag_display
        options = display_weight_table(options);
    end
    
    % Denoise of the weight table
    [options, voxels_changed, msg] = denoise_weight_table(options);
    
    if ~isempty(voxels_changed)
        bst_report('Warning', sProcess, sInput, msg);
    end
    
    % Compute Optimal Montage
    [ChannelMats, montageSufix, infos] = compute_optimal_montage(options);
    OutputFile = cell(1, length(ChannelMats));
    for iChannel = 1 :length(ChannelMats)
        ChannelMat = ChannelMats(iChannel);
        
        bst_report('Info', sProcess, sInput, infos{iChannel});
        
        %Deal with multiple versions of the folders
        sSubjStudies      = bst_get('StudyWithSubject', sSubject.FileName, 'intra_subject', 'default_study');
        condition_name    = [options.condition_name, '_', montageSufix{iChannel}];

        if any(strcmp(condition_name, {sSubjStudies.Name}))
            list_cond = union({sSubjStudies.Name}, [condition_name, '_']);

            [~, tag]    = file_unique([condition_name, '_'], list_cond, 1);
            condition_name = [condition_name, tag];
        end
        
        iStudy = db_add_condition(sSubject.Name, condition_name);
        sStudy = bst_get('Study', iStudy);
        
        % Save channel definition

        [~, iChannelStudy] = bst_get('ChannelForStudy', iStudy);
        db_set_channel(iChannelStudy, ChannelMat, 1, 0);
            
        % Save time-series data
        sDataOut              = db_template('data');
        sDataOut.F            = process_nst_separations('Compute',ChannelMat.Channel) * 100;
        sDataOut.Comment      = 'Separations';
        sDataOut.ChannelFlag  = ones(length(ChannelMat.Channel),1);
        sDataOut.Time         = (1);
        sDataOut.DataType     = 'recordings';
        sDataOut.nAvg         = 1;
        sDataOut.DisplayUnits = 'cm';
        sDataOut = bst_history('Add', sDataOut, 'Compute', 'Computed optimal montage'  );

        % Generate a new file name in the same folder
        OutputFile{iChannel} = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_chan_dist');
        bst_save(OutputFile{iChannel} , sDataOut, 'v7');
        
        % Register in database
        db_add_data(iStudy, OutputFile{iChannel} , sDataOut);
    end
    bst_report('Open', 'current'); 
end

function succeeded = check_cplex(cplex_url)
% @========================================================================
% check_cplex Function that checks if the cplex software is accessible by 
% matlab 
% ========================================================================@
    num_try     = 3;
    succeeded   = false;
    while num_try > 0 & ~succeeded
    
        num_try = num_try -1;
        try
            cplx = Cplex();
    
            if bst_plugin('CompareVersions', cplx.getVersion(),'12.3')  < 0 
                throw(['CPLEX >12.3 required. See ' cplex_url]);
            end
    
            succeeded = true;
        catch e
    
            selpath = uigetdir([], sprintf('%s. Please select the cplex directory', e.message));
            if selpath
                addpath(genpath(selpath))
            end
        end
    end
end

function [status, error, options] = check_user_inputs(options)
    status = 1;
    error = {};
    if isfield(options, 'ROI_head')
        if ~isempty(options.ROI_head)
            mandatory_fields = {'surface', 'ROI_cortex', 'Atlas_cortex', 'ROI_head', 'Atlas_head', 'SubjectName',  'outputdir', 'nb_sources', 'nb_detectors', 'nAdjacentDet', 'sep_optode', 'sepmin_SD', 'wavelengths', 'condition_name', 'data_source', 'exist_weight'};
        else
            mandatory_fields = {'surface', 'ROI_cortex', 'Atlas_cortex', 'ROI_head', 'Atlas_head', 'Extent', 'SubjectName',  'outputdir', 'nb_sources', 'nb_detectors', 'nAdjacentDet', 'sep_optode', 'sepmin_SD', 'wavelengths', 'condition_name', 'data_source', 'exist_weight'};
        end
    end
    

    C = setdiff(mandatory_fields, fieldnames(options));
    
    if ~isempty(C)
        status = 0;
        error{end+1} = sprintf('Missing option fields: %s', strjoin(C, ', '));
        
        return;
    end
    
    options.sep_optode_min  = options.sep_optode(1);
    options.sep_optode_max  = options.sep_optode(2);
    options.sep_SD_min      = options.sepmin_SD;
    
    if ~isfield(options, 'include_coverage') || ~isfield(options, 'lambda_coverage') 
        options.include_coverage = 0;
        options.lambda_coverage  = [];
    end

    % Check wavelength input
    try
        options.wavelengths = str2double(strsplit(options.wavelengths,','));
    catch
        options.wavelengths = [];
    end
    
    if(length(options.wavelengths) ~= 1) 
        status = 0;
        error{end+1} = 'User must specify a wavelength. ';
        
    end
    
    % Check coverage input
    if options.include_coverage
        if options.lambda_coverage(1) > options.lambda_coverage(3) && options.lambda_coverage(3) > 0
            status = 0;
            error{end+1} = 'Min Lambda value cannot be greater than final value. ';
        end

        if options.lambda_coverage(2) == 0
            status = 0;
            error{end+1} = 'Step value of lambda cannot be 0 : Infinite loop. ';
        end
    end
    
    is_positive_integer = @(x) ~isnan(x) && x > 0 && round(x) == x;
    
    if ~is_positive_integer(options.nb_sources)
        status = 0;
        error{end+1} = 'Number of sources must be a positive integer. ';
    end
    
    if ~is_positive_integer(options.nb_detectors)
        status = 0;
        error{end+1} = 'Number of detectors must be a positive integer. ';
    end
    
    if ~is_positive_integer(options.nAdjacentDet)
        status = 0;
        error{end+1} = 'Number of adjacence must be a positive integer. ';
    end
    
    if ~isempty(options.outputdir) && ~isfolder(options.outputdir)
        status = 0;
        error{end+1} = 'Folder for weight table must exist. ';
    end
    
    if options.sep_optode_min > options.sep_SD_min
        status = 0;
        error{end+1} = 'The minimum distance between source and detector has to be larger than the minimum optodes distance';
    end
end

function options = get_weight_tables(sSubject, sProcess, sInput, options, ROI_cortex)

    sensitivity_mat = [];
    coverage_mat    = [];
    listVertexSeen  = [];
    maxVertexSeen   = 0;
    ROI_head = options.ROI_head;

    sMri        = in_mri_bst (sSubject.Anatomy(sSubject.iAnatomy).FileName);
    voronoi_fn  = process_nst_compute_voronoi('get_voronoi_fn', sSubject);

    if ~exist(voronoi_fn, 'file')
        bst_error('ERROR: Missing voronoi volume to surface interpolator.');
    end
        
    sVoronoi = in_mri_bst(voronoi_fn);
    voronoi     = sVoronoi.Cube;
        
    % Load segmentation
    segmentation_name = 'segmentation_5tissues';
    iseg = find(strcmp({sSubject.Anatomy.Comment}, segmentation_name));
    if isempty(iseg)
        bst_error(sprintf('ERROR: Please import segmentation file as MRI and rename it as "%s"', segmentation_name));
        return;
    end
    
    GM_mask = process_nst_compute_voronoi('get_grey_matter_mask',sSubject.Anatomy(iseg).FileName);
    voronoi_mask = (voronoi > -1) &  ~isnan(voronoi) & GM_mask & ismember(voronoi, ROI_cortex.Vertices);
    
    voi               = zeros(options.cubeSize(1), options.cubeSize(2), options.cubeSize(3));
    voi(voronoi_mask) = 1;
    
    voi_flat    = voi(:);
    vois        = sparse(voi_flat > 0);
    
    if nnz(vois) == 0
        bst_error('VOI mask is empty after intersecting with Voronoi and GM masks');
        return;
    end

    weight_cache = struct('name', {}, 'sensitivity_mat', {}, 'coverage_mat', {}, 'listVertexSeen', {}, 'maxVertexSeen', {}, 'options', {});
    
    if exist(fullfile(options.outputdir , 'weight_tables.mat'), 'file')
        tmp = load (fullfile(options.outputdir, 'weight_tables.mat'));
        
        if isfield(tmp,'weight_cache') && ~isempty(tmp.weight_cache)
            weight_cache = tmp.weight_cache;
        end
        
        if options.exist_weight && isfield(weight_cache, 'sensitivity_mat') && any(strcmp( {weight_cache.name}, ROI_cortex.Label))
            
            tmp = weight_cache(strcmp( {weight_cache.name}, ROI_cortex.Label));

            if isequal(tmp.options.head_vertex_ids, ROI_head.head_vertex_ids) && tmp.options.sep_SD_min == options.sep_SD_min && tmp.options.sep_optode_max == options.sep_optode_max
                sensitivity_mat = tmp.sensitivity_mat;
                coverage_mat    = tmp.coverage_mat;
                listVertexSeen  = tmp.listVertexSeen;
                maxVertexSeen   = tmp.maxVertexSeen;
            end
            
        elseif options.exist_weight && ~isfield(weight_cache, 'sensitivity_mat')
            file_delete(fullfile(options.outputdir, 'weight_tables.mat'), 1, 1);
            bst_report('Warning', sProcess, sInput, 'Weight table format updated. Old file has been deleted and WT has been recomputed.');
        end    
    end
    
    if isempty(sensitivity_mat) || isempty(coverage_mat)
        
        % Request fluences fliles
        local_cache_dir = bst_fullfile(nst_get_local_user_dir(),  'fluence', nst_protect_fn_str(sMri.Comment));
        if contains(options.data_source, 'http')
            if ~process_nst_import_head_model('fluence_is_available', sMri.Comment)
                bst_error(['Precomputed fluence data not available for anatomy "' sMri.Comment '"']);
                return;
            end

            options.data_source = fullfile(options.data_source, nst_protect_fn_str(sMri.Comment));
        end


        [fluence_fns, missing_fluences] = process_nst_import_head_model('request_fluences', ...
                                                                        options.data_source, ...
                                                                        ROI_head.head_vertex_ids, ...
                                                                        options.wavelengths, ...
                                                                        local_cache_dir);

        % If missing fluences, list them, and return.
        if ~isempty(missing_fluences)   
            bst_error(process_nst_import_head_model('list_missing_fluences', missing_fluences));
            return;
        end

        % Load and mask fluences
        [fluence_volumes, reference] = process_nst_import_head_model('load_fluence_with_mask', ...
                                                                     fluence_fns, ...
                                                                     options.cubeSize, ...
                                                                     vois);

        % Compute weight table
        [sensitivity_mat, coverage_mat, listVertexSeen, maxVertexSeen] = compute_weights(fluence_volumes, ...
                                                                                         ROI_head.head_vertices_coords, ...
                                                                                         reference, ...
                                                                                         options);
        
        % Save the weight table in cache
        if ~isempty(options.outputdir)
            options_out = options;
            options_out.head_vertex_ids = ROI_head.head_vertex_ids;

            tmp = struct();
            tmp.name = ROI_cortex.Label; tmp.sensitivity_mat = sensitivity_mat; tmp.coverage_mat = coverage_mat; 
            tmp.listVertexSeen = listVertexSeen; tmp.maxVertexSeen = maxVertexSeen; tmp.options = options_out;
            
            if isempty(weight_cache)
                weight_cache = tmp;
            else
                
                % Force the new structure
                try
                    weight_cache(end+1) = tmp;
                catch
                    weight_cache = tmp;
                end

            end
            
            save(fullfile(options.outputdir, 'weight_tables.mat'), 'weight_cache');
        end
    end
    options.sensitivity_mat = sensitivity_mat;
    options.coverage_mat    = coverage_mat;
    options.listVertexSeen  = listVertexSeen;
    options.maxVertexSeen   = maxVertexSeen;
end

function [sensitivity_mat, coverage_mat, listVertexSeen, maxVertexSeen] = compute_weights(fluence_volumes, head_vertices_coords, reference, options)

    holder_distances = nst_pdist(head_vertices_coords, head_vertices_coords).*1000; % mm
    nHolders = size(head_vertices_coords, 1);
    iwl = 1;

    mat_sensitivity_idx = zeros(2,nHolders);
    mat_sensitivity_val = zeros(1,nHolders);
    n_sensitivity_val   = 1;

    mat_coverage_idx = zeros(2, nHolders);
    mat_coverage_val = zeros(1, nHolders);
    n_coverage_val   = 1;
    
    bst_progress('start', 'Compute weights','Preparation of fluences...', 1, 2);

        
    fluences =  zeros(size(full(fluence_volumes{1}{iwl}),1),nHolders);
    for isrc=1:nHolders
          fluences(:,isrc) =  full(fluence_volumes{isrc}{iwl});
    end
    bst_progress('inc', 1);
    
    ref = zeros(nHolders,nHolders);
    for isrc=1:nHolders
        for idet=1:nHolders
            if holder_distances(isrc, idet) > options.sep_SD_min && holder_distances(isrc, idet) < options.sep_optode_max                   
                ref(isrc,idet) = full(reference{isrc}{iwl}(idet));  
            end
        end    
    end
    
    bst_progress('inc', 1);

    %threshold for coverage
    p_thresh = 1;
    act_vol = 1000; % A definir comme un parametre donne par l'utilisateur
    delta_mu_a = 0.1;

    sSubject    = bst_get('Subject', options.SubjectName);
    voronoi_fn  = process_nst_compute_voronoi('get_voronoi_fn', sSubject);

    if ~exist(voronoi_fn, 'file')
        error('Could not find the required Voronoi file.');
    end
    sVoronoi = in_mri_bst(voronoi_fn);

    median_volume = process_nst_compute_voronoi('get_median_voronoi_volume', sVoronoi);
    
    %TODO : tune parameters
    threshold = process_nst_extract_sensitivity_from_head_model('compute_threshold', p_thresh, act_vol, median_volume, delta_mu_a);
    
    bst_progress('start', 'Compute weight tables','Computing summed sensitivities of holder pairs...', 1, nHolders^2);    
    
    isVertexSeen = false( size(fluences,1), 1);
    listVertexSeen = cell(nHolders, nHolders);
    
    for isrc=1:nHolders
        fluenceSrc = fluences(:,isrc);
        for idet=1:nHolders
            if holder_distances(isrc, idet) > options.sep_SD_min && holder_distances(isrc, idet) < options.sep_optode_max
                if ref(isrc,idet) ~=0 
                    fluenceDet = fluences(:,idet);
                    sensitivity     = (fluenceSrc .* fluenceDet) /  ref(isrc,idet); 
                    isVertexSeen    = isVertexSeen | (sensitivity > 0);
                    
                    mat_sensitivity_idx(1,n_sensitivity_val) = isrc; mat_sensitivity_idx(2,n_sensitivity_val) = idet; mat_sensitivity_val(n_sensitivity_val) = sum(sensitivity);
                    n_sensitivity_val = n_sensitivity_val + 1;
                    
                    coverage = sum(sensitivity > threshold);
                    listVertexSeen{isrc, idet} = find(sensitivity > threshold);

                    mat_coverage_idx(1, n_coverage_val) = isrc; mat_coverage_idx(2, n_coverage_val) = idet; mat_coverage_val(n_coverage_val) = coverage;
                    n_coverage_val = n_coverage_val + 1;
                end
            end    
        end
        bst_progress('inc', nHolders);
    end
    
    maxVertexSeen = sum(isVertexSeen);
    % Convert mat_coverage_val to % of the ROI that can be seen
    mat_coverage_val = mat_coverage_val ./ maxVertexSeen;
    
    % Generate the matrices
    sensitivity_mat = sparse(mat_sensitivity_idx(1,1:n_sensitivity_val-1),mat_sensitivity_idx(2,1:n_sensitivity_val-1), mat_sensitivity_val(1:n_sensitivity_val-1),nHolders, nHolders); 
    coverage_mat    = sparse(mat_coverage_idx(1,1:n_coverage_val-1),      mat_coverage_idx(2,1:n_coverage_val-1),       mat_coverage_val(1:n_coverage_val-1), nHolders, nHolders);
    
    bst_progress('stop');  
end

function [ChannelMat, montageSufix, infos] = compute_optimal_montage(options)
    infos = {};
    head_vertices_coords = options.ROI_head.head_vertices_coords;
    %======================================================================
    % 1) Compute OM by maximizing sensitivity only
    %======================================================================
    
    % Define the cplex problem
    [cplex, options] = define_prob(options.sensitivity_mat, head_vertices_coords, options);
    
    % Progress bar
    bst_progress('start', 'Optimization','Running optimization with Cplex. May take several minutes (see matlab console) ...');
  
    results = cplex.solve();
    if ~isfield(results, 'x')
        bst_error(['OM computation failed  at Cplex step:', results.statusstring]);
        return;
    end

    bst_progress('stop');

    % Calculation of montage_pairs matrix, montage_sensitivity and montage_coverage vector
    [montage_pairs_simple, montage_sensitivity_simple, montage_coverage_simple, channels_coverage_simple] = montage_pairs_and_weight(results,options);


    ChannelMat      = create_channelMat_from_montage(montage_pairs_simple, head_vertices_coords, options.wavelengths);
    montageSufix{1} = 'simple';
    
    % Premature ending in case coverage constraint is not asked
    if ~options.include_coverage
        str = [sprintf('Only sensitivity : \n'), ...
              display_channel_info(montage_pairs_simple, montage_sensitivity_simple, montage_coverage_simple, channels_coverage_simple, head_vertices_coords)];
        
        infos{end+1} = str;
        return;
    end
    
    %======================================================================
    % 2) Compute OM by maximizing sensitivity and coverage
    %======================================================================
    
    cov_min  = options.lambda_coverage(1);
    cov_step = options.lambda_coverage(2);
    cov_max  = options.lambda_coverage(3);
    lambda2  = (cov_min:cov_step:cov_max);
        
    % Display
    str = [sprintf('Only sensitivity : \n'), ...
          display_channel_info(montage_pairs_simple, montage_sensitivity_simple, montage_coverage_simple, channels_coverage_simple, head_vertices_coords)];

    infos{end+1} = str;

    % Define the cplex problem
    for iLambda = 1:length(lambda2)
        options.lambda1 = 1/sum(montage_sensitivity_simple);
        options.lambda2 = lambda2(iLambda);

        wt =  options.lambda1  * options.sensitivity_mat  + options.lambda2 * options.coverage_mat;
        
        if options.flag_display
            options = display_weight_table(options);
        end
        
        [cplex, options] = define_prob(wt, head_vertices_coords, options);
    
        % Progress bar
        bst_progress('start', 'Optimization','Running optimization with Cplex. May take several minutes (see matlab console) ...');
        
        results = cplex.solve();
        if ~isfield(results, 'x')
            bst_error(['OM computation failed  at Cplex step:', results.statusstring]);
            return;
        end
    
        bst_progress('stop');
                                                        
        % Calculation of montage_pairs matrix, montage_sensitivity and montage_coverage vector
        [montage_pairs, montage_sensitivity, montage_coverage, channels_coverage] = montage_pairs_and_weight(results, options);

        str = [sprintf('Sensitivity and Coverage (lambda = %d):\n', lambda2(iLambda)), ...
              display_channel_info(montage_pairs, montage_sensitivity, montage_coverage, channels_coverage, head_vertices_coords)];

        infos{end+1} = str;

        % Convert Montage to Brainstorm structure
        ChannelMat = [ ChannelMat , create_channelMat_from_montage(montage_pairs, head_vertices_coords, options.wavelengths)];
        montageSufix{end+1} = sprintf('complex_lambda_%d', lambda2(iLambda));
    end
end

function [head_vertices, sHead, sSubject] = proj_cortex_scout_to_scalp(cortex_scout, extent_m, save_in_db)

    if nargin < 3
        save_in_db = 0;
    end

    sSubject        = cortex_scout.sSubject;
    sHead           = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);
    sCortex         = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);

    dis2head        = nst_pdist(sHead.Vertices, sCortex.Vertices(cortex_scout.sScout.Vertices,:));
    head_vertices   = find(min(dis2head,[],2) < extent_m); 
    
    iAtlasHead = find(strcmp('User scouts', {sHead.Atlas.Name}));
    exclude_scout = sHead.Atlas(iAtlasHead).Scouts(strcmp('FluenceExclude', {sHead.Atlas(iAtlasHead).Scouts.Label}));
    if ~isempty(exclude_scout)
        head_vertices = setdiff(head_vertices, exclude_scout.Vertices);
    end
    
    limiting_scout = sHead.Atlas(iAtlasHead).Scouts(strcmp('FluenceRegion', {sHead.Atlas(iAtlasHead).Scouts.Label}));
    if ~isempty(limiting_scout)
        head_vertices = intersect(head_vertices, limiting_scout.Vertices);
    end
    
    if save_in_db 

        new_scout_name = ['From cortical ' cortex_scout.sScout.Label '(' num2str(extent_m*100) ' cm)'];
        scout_idx = find(strcmp(new_scout_name, {sHead.Atlas(iAtlasHead).Scouts.Label}));

        if isempty(scout_idx)
            scout_idx = size(sHead.Atlas(iAtlasHead).Scouts,2) + 1;
        end

        sHead.Atlas(iAtlasHead).Scouts(scout_idx) = db_template('Scout');
        sHead.Atlas(iAtlasHead).Scouts(scout_idx).Vertices = head_vertices';
        sHead.Atlas(iAtlasHead).Scouts(scout_idx).Seed = head_vertices(1);
        sHead.Atlas(iAtlasHead).Scouts(scout_idx).Color = [0,0,0];
        sHead.Atlas(iAtlasHead).Scouts(scout_idx).Label = new_scout_name;
        bst_save(file_fullpath(sSubject.Surface(sSubject.iScalp).FileName), sHead, 'v7');

        db_save();

    end

end

function [cplex, options] = define_prob(weight_table, head_vertices_coords, options)
% @========================================================================
% define_prob_simple Initializes the problem
% Added in options by function : holder_distances, nH, thresh_sep_optode_optode 
% ========================================================================@
    
    holder_distances = nst_pdist(head_vertices_coords, head_vertices_coords).*1000; % mm

    nHolders = size(head_vertices_coords, 1);
    
    nS = options.nb_sources; % number of sources
    nD = options.nb_detectors; % number of detectors
    
    flag_sep_optode_optode = 1;
    thresh_sep_optode_optode = [options.sep_optode_min options.sep_optode_max]; 
    % mm (1) no optodes sep below this thresh (2) Optodes above this thresholds do not form pairs
    
    flag_sep_src_det = 1;
    thresh_min_sep_src_det = options.sep_SD_min; 
    % mm (1) min sep between a source and a det %TODO: why not using thresh_sep_optode_optode(1)?
    
    flag_adjacency = 1;
    nAdjacentDet = options.nAdjacentDet; 
    % minimal number of adjacent detector in the given range
    
    nH   = nHolders; % number of holders
    xp   = zeros(nH,1) ; % binary
    yq   = zeros(nH,1) ; % binary
    wq_V = zeros(nH,1) ; % float
    v    = [xp;yq;wq_V]; nVar = size(v,1);
    
    % Calls function for equation 1 : Sum(xp)=nSrc (11)
    [Aeq_1, E_1] = add_constraint_nSrc(nVar, nH, nS);
    
    % Calls function for equation 2 : Sum(yq)=nDet (11)
    [Aeq_2, E_2] = add_constraint_nDet(nVar, nH, nD);
    
    %......................................................................
    % DEBUG FUNCTION : Display Equality matrix
    % display_eq_matrix(Aeq_1, Aeq_2);
    %......................................................................
    
    %======================================================================
    % Create inequality matrixes
    %======================================================================
    % Calls function for inequation 1 : wq_V-Myq<=0 (12)
    [Aineq_1, I_1] = create_ineq_matrix(nH, nVar);
    
    % Calls function for inequation 2 : wq_V-Sum(Vpq*xp) <= 0 (13)
    [Aineq_2, I_2] = detect_activ_constr(nH, nVar, weight_table);
    
    % Calls function for inequation 3 : Optode/Optode minimal separation constraints
    [Aineq_3, I_3] = opt_opt_min_sep_constr(flag_sep_optode_optode, holder_distances, nH, thresh_sep_optode_optode, nS, nD, nVar);
    
    % Calls function for inequation 4 : Source/det min separation constraints
    [Aineq_4, I_4] = src_det_min_sep_constr(flag_sep_src_det, holder_distances, nH, thresh_min_sep_src_det, nD, nVar);
    
    % Calls function for inequation 5 : Adjacency constraints
    [Aineq_5, I_5] = adj_constr(flag_adjacency, holder_distances, nH, thresh_sep_optode_optode, nVar, nAdjacentDet);
    
    %......................................................................
    % DEBUG FUNCTION : Display InEquality matrix Aineq_1 to 5
    % display_ineq_matrix(Aineq_1, Aineq_2, Aineq_3, Aineq_4, Aineq_5);
    %......................................................................
    
    Aeq = [Aeq_1 ; Aeq_2];
    E   = [E_1 ; E_2];
    
    Aineq   = [Aineq_1 ; Aineq_2 ; Aineq_3 ; Aineq_4 ; Aineq_5];
    I       = [I_1 ; I_2 ; I_3 ; I_4 ; I_5];
    
    f       = [zeros(1, nH) , zeros(1, nH) ones(1, nH)]';
    lb      = []; 
    ub      = []; 
    ctype   = [repmat('B', 1, 2*nH) repmat('S', 1, nH)];
    
    %Cplex optimisation
    prob        = cplexcreateprob('cplexmilp');
    prob.f      = f;
    prob.lb     = lb;       prob.ub = ub;
    prob.ctype  = ctype;
    prob.Aineq  = Aineq;    prob.bineq = I;
    prob.Aeq    = Aeq;      prob.beq = E;
    prob.x0         = [];
    prob.options    = [];

    cplex=Cplex(prob);
    cplex.Model.sense = 'maximize';
    cplex.Param.timelimit.Cur=300;

    % Delete clone[number].log files created by Cplex
    cplex.Param.output.clonelog.Cur = 0;

    options.holder_distances = holder_distances;
    options.thresh_sep_optode_optode = thresh_sep_optode_optode; 
    options.nH = nH;
    
end

%==========================================================================
% FUNCTIONS USED TO DEFINE THE OPTIMISATION PROBLEM
%==========================================================================

function [A, E] = add_constraint_nSrc(nVar, nH, nS)
% @========================================================================
% add_constraint_nSrc Adds the constraint for the total number of sources.
%   See equation 11 : Sum(xp) = nSrc
% ========================================================================@

    A = zeros(1, nVar);
    A(1, 1:nH) = 1;

    E = nS;
end

function [A, E] = add_constraint_nDet(nVar, nH, nD)
% @========================================================================
% add_constraint_nDet Adds the constraint for the total number of detectors.
%   See equation 11 : Sum(yq)=nDet
% ========================================================================@

    A = zeros(1, nVar);
    A(1, nH+1:2*nH) = 1;

    E = nD;
end

function [A, I] = create_ineq_matrix(nH, nVar)
% @========================================================================
% create_ineq_matrix Creates the inequality matrix. 
%   See equation 12 : wq_V-Myq <= 0
% ========================================================================@

    M = realmax;
    A = zeros(nH, nVar);
    for iH = 1:nH
        A(iH, nH+iH) = -M;
        A(iH, 2*nH+iH) = 1;
    end
    I = zeros(nH, 1);
end

function [A, I] = detect_activ_constr(nH, nVar, weight_table)
% @========================================================================
% detect_activ_constr Generates inequality constraints for detector activation.
%   See equation 13 : wq_V-Sum(Vpq*xp) <= 0
% ========================================================================@

    A=zeros(nH, nVar);
    for iH = 1:nH
        A(iH, 2*nH+iH) = 1;
        A(iH, 1:nH) = -weight_table(:, iH); %Note: this takes into account weight assymetry
    end
    I = zeros(nH, 1);
end

function [A, I] = opt_opt_min_sep_constr(flag_sep_optode_optode, holder_distances, nH, thresh_sep_optode_optode, nS, nD, nVar)
% @========================================================================
% opt_opt_min_sep_constr Optode/Optode minimal separation constraints
% ========================================================================@

    if flag_sep_optode_optode
        separations = holder_distances;
        forbid_sets = zeros(nH);
        forbid_sets(separations < thresh_sep_optode_optode(1)) = 1; %TODO expose as separate parameter
    
        M2 = nS+nD;
        A = sparse(2*nH, nVar);
        A(1:nH, 1:2*nH) = [(M2-1)*eye(nH) + forbid_sets  forbid_sets] ;
        A(nH +1:2*nH, 1:2*nH) = [forbid_sets (M2-1)*eye(nH) + forbid_sets] ;
        I = M2*ones(2*nH, 1);
    else
        A = [];
        I = [];
    end
end

function [A, I] = src_det_min_sep_constr(flag_sep_src_det, holder_distances, nH, thresh_min_sep_src_det, nD, nVar)
% @========================================================================
% src_det_min_sep_constr Source/det minimal separation constraints
% ========================================================================@

    if flag_sep_src_det
        separations = holder_distances;
        forbid_sets = zeros(nH);
        forbid_sets(separations < thresh_min_sep_src_det) = 1;
    
        M2 = nD;
        A = sparse(nH, nVar);
        A(1:nH, 1:2*nH) = [M2*eye(nH) forbid_sets] ;
        I = M2*ones(nH, 1);
    else
        A = [];
        I = [];
    end
end

function [A, I] = adj_constr(flag_adjacency, holder_distances, nH, thresh_sep_optode_optode, nVar, nAdjacentDet)
% @========================================================================
% adj_constr Adjacency constraints
% ========================================================================@

    if flag_adjacency
        separations = holder_distances;
        Allow_sets = zeros(nH);
        Allow_sets(separations >= thresh_sep_optode_optode(1) & separations <= thresh_sep_optode_optode(2)) = 1;
    
        A = sparse(nH,nVar);
        A(1:nH,1:2*nH) = [eye(nH)*nAdjacentDet, -Allow_sets];
        I = zeros(nH, 1);
    else
        A = [];
        I = [];
    end
end

function [montage_pairs, montage_sensitivity, montage_coverage, channels_coverage] = montage_pairs_and_weight(results,options)
% @========================================================================
% montage_pairs_and_weight Calculation of montage pairs matrix and montage
% weight vector
% ========================================================================@
    
    x=results.x;
    x=round(x);
    isources = find(x(1:options.nH)==1);
    idetectors = find(x(options.nH+1:2*options.nH)==1);
    
    % Memory management
    max_pairs = length(isources) * length(idetectors);
    montage_pairs = zeros(max_pairs, 2);
    montage_sensitivity = zeros(max_pairs, 1);
    channels_coverage = zeros(max_pairs, 1);
    
    ipair = 0;

    for isrc = 1:length(isources)
        for idet = 1:length(idetectors)
            if options.holder_distances(isources(isrc), idetectors(idet)) > options.thresh_sep_optode_optode(1) && ...
                    options.holder_distances(isources(isrc), idetectors(idet)) < options.thresh_sep_optode_optode(2) && ...
                    full(options.sensitivity_mat(isources(isrc), idetectors(idet)))
                
                ipair = ipair + 1;

                montage_pairs(ipair,:)          = [isources(isrc) idetectors(idet)];
                montage_sensitivity(ipair,:)    = full(options.sensitivity_mat(isources(isrc), idetectors(idet)));
                channels_coverage(ipair,:)      = full(options.coverage_mat(isources(isrc), idetectors(idet)));
            end
        end
    end
    
     % Make sure the matrix is the right size
    montage_pairs       = montage_pairs(1:ipair, :);
    montage_sensitivity = montage_sensitivity(1:ipair, :);
    coverage_mat        = options.listVertexSeen;
    
    % Compute montage coverage
    list_vertex_seen = {};
    for iPair = 1:size(montage_pairs, 1)
        pair = montage_pairs(iPair, :);

        if isempty(list_vertex_seen)
            list_vertex_seen = coverage_mat{pair(1), pair(2)}';
        else
            list_vertex_seen = union(list_vertex_seen, coverage_mat{pair(1), pair(2)}' );
        end
    end

    montage_coverage    = length(list_vertex_seen) / options.maxVertexSeen;
    channels_coverage   = channels_coverage(1:ipair, :);
end

function info = display_channel_info(montage_pairs, montage_sensitivity,  montage_coverage, channels_coverage, head_vertices_coords)
% @========================================================================
% display_channel_info Used to create the string containing the channels
% informations
% ========================================================================@

    info = '';
    src_indexes = zeros(max(montage_pairs(:, 1)), 1);
    det_indexes = zeros(max(montage_pairs(:, 2)), 1);
    det_next_idx = 1;
    src_next_idx = 1;
    tab_dist_mm = zeros(size(montage_pairs, 1), 1);

    for ipair = 1:size(montage_pairs, 1)
        ihead_vertex_src = montage_pairs(ipair, 1);
        ihead_vertex_det = montage_pairs(ipair, 2);
        
        if src_indexes(ihead_vertex_src) == 0
            idx_src = src_next_idx;
            src_indexes(ihead_vertex_src) = src_next_idx;
            src_next_idx = src_next_idx + 1;
        else
            idx_src = src_indexes(ihead_vertex_src);
        end
        
        if det_indexes(ihead_vertex_det) == 0
            idx_det = det_next_idx;
            det_indexes(ihead_vertex_det) = det_next_idx;
            det_next_idx = det_next_idx + 1;
        else
            idx_det = det_indexes(ihead_vertex_det);
        end
        
        distance_mm = nst_pdist(head_vertices_coords(ihead_vertex_src, :), head_vertices_coords(ihead_vertex_det, :)) * 1000;
        tab_dist_mm(ipair) = distance_mm;
        sensitivity = montage_sensitivity(ipair, :);
        coverage = channels_coverage(ipair, :);

        info = [info, sprintf('Channel S%02dD%02d >>> Distance: %4.1f mm    Sensitivity: %6.3f    Coverage: %5.2f%%  \n', idx_src, idx_det, distance_mm, sensitivity, coverage * 100)];
    end
    
    mean_distance = mean(tab_dist_mm);
    distance_range = [min(tab_dist_mm), max(tab_dist_mm)];
    total_sensitivity = sum(montage_sensitivity);
    
    percentage_overlap = 1 - montage_coverage / sum(channels_coverage);
    if ~isempty(distance_range)
        info = [info, sprintf('TOTAL          >>> Distance (mean/range): %.1f mm [%.1f-%.1f]    Total sensitivity: %.3f Total coverage: %.3f%%  Overlap measure: %.3f%% \n', mean_distance, distance_range(1), distance_range(2), total_sensitivity, 100*montage_coverage, 100*percentage_overlap)];
        info = strrep(info, '  ', '&nbsp;&nbsp;');
    end
end

%==========================================================================
% DEBUG FUNCTIONS
%==========================================================================

function options = display_weight_table(options)
% @========================================================================
% display_weight_table Displays multiple graphs useful to see a
% representation of the sensitivity & coverage matrices
% ========================================================================@
    if ~isfield(options, 'hFig')
        options.hFig = figure;
        options.hFigTab = uitabgroup; drawnow;
    end

    sensitivity_mat = options.sensitivity_mat;
    coverage_mat    = options.coverage_mat;
    ROI_head        = options.ROI_head;
    hFigTab         = options.hFigTab;
    isDefinedLambda = 0;
    
    if isfield(options, 'lambda1')
        sensitivity_mat = options.lambda1 * sensitivity_mat;
        isDefinedLambda = isDefinedLambda + 1;
    end
    if isfield(options, 'lambda2') && options.lambda2 > 0
        coverage_mat = options.lambda2 * coverage_mat;
        isDefinedLambda = isDefinedLambda + 1;
    end
    
    if isDefinedLambda == 0
        onglet = uitab(hFigTab,'title','Weight tables');
    elseif isDefinedLambda == 1
        onglet = uitab(hFigTab,'title','Lambda = 0');
    else
        onglet = uitab(hFigTab,'title',['Lambda = ', num2str(options.lambda2)]);
    end

    distances = squareform(pdist(ROI_head.head_vertices_coords));
    [~, I] = min(median(distances));
    [~, order] = sort( abs(distances(I, :) - median(distances(I, :))));
    sensitivity_mat = sensitivity_mat(order,order);
    coverage_mat  = coverage_mat(order,order);
    
    hpc = uipanel('Parent', onglet, ...
              'Units', 'Normalized', ...
              'Position', [0.01 0.01 0.98 0.98], ...
              'FontWeight','demi');

    ax1 = subplot(1, 3, 1, 'parent', hpc);
    imagesc(ax1, sensitivity_mat);
    title(ax1, 'Sensitivity matrix');
    colorbar(ax1);
    
    ax2 = subplot(1, 3, 2, 'parent', hpc);
    imagesc(ax2, coverage_mat);
    title(ax2, 'Coverage matrix');
    colorbar(ax2);
    
    ax3 = subplot(1, 3, 3, 'parent', hpc);
    if isDefinedLambda == 0
        plot(ax3, sensitivity_mat(:), coverage_mat(:), '+')
        xlabel(ax3, 'Sensitivity');
        ylabel(ax3, 'Coverage');
        title(ax3, 'Sensitivity VS. Coverage');
        set(hpc,'Title',' Sensitivity & Coverage Matrices ','FontSize',8);
    elseif isDefinedLambda == 1
        plot(ax3, sensitivity_mat(:), coverage_mat(:), '+')
        xlabel(ax3, 'Sensitivity');
        ylabel(ax3, 'Coverage');
        title(ax3, 'Sensitivity VS. Coverage');
        set(hpc,'Title',' Lambda = 0 ','FontSize',8);
    else
        ratio = zeros(size(coverage_mat));
        idx_ratio = sensitivity_mat > 0;
        ratio(idx_ratio) =  log10( 1 +  coverage_mat(idx_ratio) ./ sensitivity_mat(idx_ratio));
        
        imagesc(ax3, ratio);
        title(ax3, 'Ratio matrix');
        colorbar(ax3);
    end
    
    %......................................................................
    % Save figures in the wt folder
    % save_figure(options, onglet, hpc);
    %......................................................................

end

function save_figure(options, onglet, hpc)
% @========================================================================
% save_figure Function used to save any figure created during the process
% ========================================================================@
    tabTitle = onglet.Title;
    safeFilename = matlab.lang.makeValidName(tabTitle);
    fullFilePath = fullfile(options.outputdir, [safeFilename, '.png']);
    exportgraphics(hpc, fullFilePath, 'Resolution', 300);
end

function display_eq_matrix(Aeq_1, Aeq_2)
    matrixes2show={Aeq_1, Aeq_2, [Aeq_1 ; Aeq_2]};
    for ii=1:2
        figure ('name', 'equality matrix')
        mat2show=matrixes2show{ii} ;
        colormap(gray(3))
        imagesc((mat2show));
        %set(gca, 'PlotBoxAspectRatio', [size(mat2show, 1) size(mat2show, 2) 1])
        ylim([0.5 size((mat2show), 1)+0.5])
    end
end

function display_ineq_matrix(Aineq_1, Aineq_2, Aineq_3, Aineq_4, Aineq_5)
    matrixes2show={Aineq_1, Aineq_2, Aineq_3, Aineq_4, [Aineq_1 ; Aineq_2 ; Aineq_3 ; Aineq_4 ; Aineq_5]};
    for ii=1:5
        figure('name','inequality matrix')
        %colormap(flipud(autumn(3)))
        mat2show=matrixes2show{ii};
        imagesc(mat2show);
        %set(gca, 'PlotBoxAspectRatio', [size(mat2show, 1) size(mat2show, 2) 1])
        ylim([0.5 size(mat2show, 1)+0.5])
    end
end

function ChannelMat = create_channelMat_from_montage(montage_pairs, head_vertices_coords, wavelengths)

    % montage_pairs contains head vertex indexes -> remap to 1-based consecutive indexes
    src_indexes = zeros(max(montage_pairs(:, 1)), 1);
    det_indexes = zeros(max(montage_pairs(:, 2)), 1);
    
    
    % Forge channels from head points pairs & wavelengths
    nChannels       = size(montage_pairs, 1) * length(wavelengths);
    iChan           = 1;
    det_next_idx    = 1;
    src_next_idx    = 1;

    ChannelMat = db_template('channelmat');
    ChannelMat.Comment = 'NIRS-BRS channels';
    ChannelMat.Channel = repmat(db_template('channeldesc'), [1, nChannels]);
    ChannelMat.Nirs.Wavelengths = wavelengths;
    ChannelMat = bst_history('Add', ChannelMat, 'Compute', 'Computed optimal montage'  );

    for ipair=1:size(montage_pairs, 1)
        ihead_vertex_src = montage_pairs(ipair, 1);
        ihead_vertex_det = montage_pairs(ipair, 2);
        
        if src_indexes(ihead_vertex_src) == 0
            idx_src = src_next_idx;
            src_indexes(ihead_vertex_src) = src_next_idx;
            src_next_idx = src_next_idx + 1;
        else
            idx_src = src_indexes(ihead_vertex_src);
        end
        if det_indexes(ihead_vertex_det) == 0
            idx_det = det_next_idx;
            det_indexes(ihead_vertex_det) = det_next_idx;
            det_next_idx = det_next_idx + 1;
        else
            idx_det = det_indexes(ihead_vertex_det);
        end

        for iwl=1:length(wavelengths)
            
            ChannelMat.Channel(iChan).Name      = sprintf('S%dD%dWL%d', idx_src, idx_det, ChannelMat.Nirs.Wavelengths(iwl));
            ChannelMat.Channel(iChan).Type      = 'NIRS';
            ChannelMat.Channel(iChan).Loc(:,1)  = head_vertices_coords(ihead_vertex_src, :);
            ChannelMat.Channel(iChan).Loc(:,2)  = head_vertices_coords(ihead_vertex_det, :);
            ChannelMat.Channel(iChan).Orient    = [];
            ChannelMat.Channel(iChan).Weight    = 1;
            ChannelMat.Channel(iChan).Comment   = [];
            ChannelMat.Channel(iChan).Group     = sprintf('WL%d', round(ChannelMat.Nirs.Wavelengths(iwl)));
            
            iChan = iChan + 1;
        end
    end
end

function [options, voxels_changed, msg] = denoise_weight_table(options)
% @========================================================================
% denoise_weight_table Denoises the sensitivity matrix and displays the
% comparaison
% ========================================================================@
    ROI_head        = options.ROI_head;
    sensitivity_mat = options.sensitivity_mat;
    coverage_mat    = options.coverage_mat;
    
    distances = squareform(pdist(ROI_head.head_vertices_coords));
    [~,I] = min(median(distances));
    [~, order] = sort( abs(distances(I, :) - median(distances(I, :))));
    sensitivity_mat = sensitivity_mat(order,order);
    coverage_mat  = coverage_mat(order,order);

    % Denoise weight table based on 4 neighbors (cross patern)
    stvty_mat_full = full(sensitivity_mat);
    cvge_mat_full  = full(coverage_mat);

    [max_r, max_c] = size(stvty_mat_full);
    
    sensitivity_mat_denoised = zeros(max_r, max_c);
    coverage_mat_denoised = cvge_mat_full;

    for r = 1 : max_r
        for c = 1 : max_c
            neighbors = list_neighbors(stvty_mat_full, r, c, max_r, max_c);
            sensitivity_mat_denoised(r, c) = median(neighbors);
            
            if sensitivity_mat_denoised(r, c) ~= stvty_mat_full(r, c)
                neighbors = list_neighbors(cvge_mat_full, r, c, max_r, max_c);
                coverage_mat_denoised(r, c) = median(neighbors);
            end
        end
    end

    max_original = max(sensitivity_mat(:));
    max_filtered = max(sensitivity_mat_denoised(:));
    nnz_original = nnz(sensitivity_mat);
    nnz_filtered = nnz(sensitivity_mat_denoised);
    tresh = 10;

    % Early return if conditions are not respected
    if max_original < tresh * max_filtered || nnz_original > tresh * nnz_filtered
        voxels_changed = [];
        msg = '';
        return;
    end

    ratio       = zeros(size(sensitivity_mat));
    idx_ratio   = sensitivity_mat > 0;
    ratio(idx_ratio)    = sensitivity_mat_denoised(idx_ratio) ./ sensitivity_mat(idx_ratio);
    voxels_changed      = any(ratio > tresh);

    sensitivity_mat = sparse(sensitivity_mat_denoised);
    coverage_mat    = sparse(coverage_mat_denoised);
    
    msg = ['The sensitivity matrix has been denoised to compensate for abnormally high values. ' ...
           'For greater accuracy when computing the head model, please recalculate the fluences ' ...
           'for the generated montage with higher number of photons.'];

    [~, inverse_order] = sort(order);
    options.sensitivity_mat = sensitivity_mat(inverse_order, inverse_order);
    options.coverage_mat    = coverage_mat(inverse_order, inverse_order);
    voxels_changed          = find(voxels_changed(inverse_order));

    if options.flag_display
        options = display_denoise_weight(options);
    end

end

function neighbors = list_neighbors(mat, r, c, max_r, max_c)
% @========================================================================
% list_neighbors Extracts the list of values arround a pixel in a + shape
% mat : Matrix studied ; (r, c) coordinates of the pixel studied ; 
% max_r : max rows ; max_c : max columns
% ========================================================================@
    neighbors = [mat(r, c)];
    if (r-1) >=1
        neighbors(end+1) = mat(r-1, c);
    end
    if (r+1) <= max_r
        neighbors(end+1) = mat(r+1, c);
    end
    if (c-1) >=1
        neighbors(end+1) = mat(r, c-1);
    end
    if (c+1) <= max_c
        neighbors(end+1) = mat(r, c+1);
    end
end

function options = display_denoise_weight(options)
% @========================================================================
% disp_denoise_weight Function used to display the denoised version of
% sentitivity and coverage matrices
% ========================================================================@
    % For display
    if ~isfield(options, 'hFig')
        options.hFig = figure;
        options.hFigTab = uitabgroup; drawnow;
    end
    
    onglet = uitab(options.hFigTab,'title','Denoise');
    
    hpc = uipanel('Parent', onglet, ...
        'Units', 'Normalized', ...
        'Position', [0.01 0.01 0.98 0.98], ...
        'FontWeight','demi');
    set(hpc,'Title',' Sensitivity matrix ','FontSize',8);
    
    ax1 = subplot(1, 3, 1, 'parent', hpc);
    imagesc(ax1, options.sensitivity_mat);
    title(ax1, 'Denoised sensitivity matrix');
    colorbar(ax1);
    
    ax2 = subplot(1, 3, 2, 'parent', hpc);
    imagesc(ax2, options.coverage_mat);
    title(ax2, 'Denoised coverage matrix');
    colorbar(ax2);
    
    ax3 = subplot(1, 3, 3, 'parent', hpc);
    plot(ax3, options.sensitivity_mat(:), options.coverage_mat(:), '+')
    xlabel(ax3, 'Sensitivity');
    ylabel(ax3, 'Coverage');
    title(ax3, 'Sensitivity VS. Coverage');
    set(hpc,'Title',' Sensitivity & Coverage Matrices ','FontSize',8);

    %......................................................................
    % Save figures in the wt folder
    % save_figure(options,onglet,hpc);
    %......................................................................
end

function [ROI_cortex, ROI_head] = get_regions_of_interest(sSubject, options)

    sHead   = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);
    sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);

    i_atlas_cortex = strcmp({sCortex.Atlas.Name},options.Atlas_cortex);
    i_scout_cortex = strcmp({sCortex.Atlas(i_atlas_cortex).Scouts.Label}, options.ROI_cortex);
    ROI_cortex     = sCortex.Atlas(i_atlas_cortex).Scouts(i_scout_cortex);   
    
    if isempty(options.Atlas_head) && isempty(options.ROI_head)
        
        cortex_to_scalp_extent = options.Extent;
        cortex_scout.sSubject = sSubject;
        cortex_scout.sScout   = ROI_cortex;
        
        [head_vertex_ids, ~, ~] = proj_cortex_scout_to_scalp(cortex_scout, cortex_to_scalp_extent.*0.01, 1); 

    else    

        i_atlas_head = strcmp({sHead.Atlas.Name},options.Atlas_head);
        i_scout_head = strcmp({sHead.Atlas(i_atlas_head).Scouts.Label}, options.ROI_head);
        head_vertex_ids = sHead.Atlas(i_atlas_head).Scouts(i_scout_head).Vertices;  
        
        exclude_scout = sHead.Atlas(i_atlas_head).Scouts(strcmp('FluenceExclude', {sHead.Atlas(i_atlas_head).Scouts.Label}));
        if ~isempty(exclude_scout)
            head_vertex_ids = setdiff(head_vertex_ids, exclude_scout.Vertices);
        end
        
    end
    
    % Extract head point coordinates
    head_vertices_coords = sHead.Vertices(head_vertex_ids, :);
    
    ROI_head = struct('head_vertex_ids',head_vertex_ids, 'head_vertices_coords', head_vertices_coords);
end