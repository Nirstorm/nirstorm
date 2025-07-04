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

% Description the process
sProcess.Comment     = 'Compute optimal montage';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = {'NIRS', 'Sources'};
sProcess.Index       = 1406;
sProcess.Description = '';
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
    try
        cplx = Cplex();
        if bst_plugin('CompareVersions', cplx.getVersion(),'12.3')  < 0 
            bst_error(['CPLEX >12.3 required. See ' cplex_url]);
            return
        end
    catch
        bst_error(['CPLEX >12.3 required. See ' cplex_url]);
        return
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
    
    [ROI_cortex, ROI_head] = get_regions_of_interest(sSubject, options);
        
    
    try
        wavelengths = str2double(strsplit(options.wavelengths,','));
    catch
        bst_report('Error', sProcess, [], 'List of wavelengths must be integers separated by comas');
        return
    end
    
    if(length(wavelengths) < 1) % TODO: at least 2 wavelengths ?
        bst_report('Error', sProcess, [], 'List of wavelengths must not be empty');
        return
    end
    
    % Obtain the anatomical MRI
    sMri     = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
    cubeSize = size(sMri.Cube);
    
    options.sep_optode_min  = options.sep_optode(1);
    options.sep_optode_max  = options.sep_optode(2);
    options.sep_SD_min      = options.sepmin_SD;
    options.cubeSize        = cubeSize;
    options.wavelengths     = wavelengths;
    
    if options.sep_optode_min > options.sep_SD_min
        bst_error(sprintf('ERROR: The minimum distance between source and detector has to be larger than the minimum optodes distance'));
        return;
    end
    
    [options.sensitivity_mat, options.coverage_mat] = get_weight_tables(sSubject, options, ROI_cortex, ROI_head);
    if isempty(options.sensitivity_mat) || nnz(options.sensitivity_mat) == 0
        bst_error(sprintf('Weight table is null for ROI: %s', ROI_cortex.Label));
        return
    end
    
    % Experimental : Denoise the weight table. (be carefull)
    % options.sensitivity_mat = denoise_weight_table(options.sensitivity_mat, threshold);
    
    % Compute Optimal Montage
    montage_pairs_simple, montage_sensitivity_simple, montage_coverage_simple, montage_pairs, montage_sensitivity, montage_coverage = compute_optimal_montage(ROI_head.head_vertices_coords,options);
    
    % Convert Montage to Brainstorm structure
    %A MODIFIER
    ChannelMatSimple                = create_channelMat_from_montage(montage_pairs_simple, montage_sensitivity_simple, montage_coverage_simple, ROI_head.head_vertices_coords, options.wavelengths);
    ChannelMat                      = create_channelMat_from_montage(montage_pairs, montage_sensitivity, montage_coverage, ROI_head.head_vertices_coords, options.wavelengths);
    
    
    sSubjStudies      = bst_get('StudyWithSubject', sSubject.FileName, 'intra_subject', 'default_study');
    options.condition_name    = file_unique(options.condition_name, {sSubjStudies.Name}, 1);
    
    iStudy = db_add_condition(sSubject.Name, options.condition_name);
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
    
    % Generate a new file name in the same folder
    OutputFile{1} = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_chan_dist');
    bst_save(OutputFile{1} , sDataOut, 'v7');
    
    % Register in database
    db_add_data(iStudy, OutputFile{1} , sDataOut);
    
end

function [sensitivity_mat, coverage_mat] = get_weight_tables(sSubject, options, ROI_cortex, ROI_head)

    sensitivity_mat = [];
    coverage_mat    = [];

    voronoi_fn = process_nst_compute_voronoi('get_voronoi_fn', sSubject);
    
    if ~exist(voronoi_fn, 'file')
        voronoi_fn = bst_process('CallProcess', 'process_nst_compute_voronoi', [], [], ...
                    'subjectname',        sSubject.Name, ...
                    'segmentation_label', options.segmentation_label);
    end
        
    voronoi_bst = in_mri_bst(voronoi_fn);
    voronoi     = voronoi_bst.Cube;
        
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

    weight_cache = struct();
    
    if exist(fullfile(options.outputdir , 'weight_tables.mat'), 'file')
        tmp = load (fullfile(options.outputdir, 'weight_tables.mat'));
        
        if isfield(tmp,'weight_cache') && ~isempty(tmp.weight_cache)
            weight_cache = tmp.weight_cache;
        end
        
        if options.exist_weight && isfield(weight_cache,  matlab.lang.makeValidName(ROI_cortex.Label))
            tmp = weight_cache.(matlab.lang.makeValidName(ROI_cortex.Label));    
            if isequal(tmp.options.head_vertex_ids, ROI_head.head_vertex_ids) && tmp.options.sep_SD_min == options.sep_SD_min &&  tmp.options.sep_optode_max == options.sep_optode_max   
                sensitivity_mat = tmp.sensitivity_mat;
                coverage_mat = tmp.coverage_mat;

            end    
        end    
    end
    
    if isempty(sensitivity_mat)
        
        % Compute weight table
        [fluence_volumes, reference] = process_nst_import_head_model('request_fluences', ...
            ROI_head.head_vertex_ids, ...
            sSubject.Anatomy(sSubject.iAnatomy).Comment, ...
            options.wavelengths, ...
            options.data_source, ...
            [], ...
            vois, ...
            options.cubeSize);
        
        [sensitivity_mat, coverage_mat] = compute_weights(fluence_volumes, ROI_head.head_vertices_coords, reference, options);
        
        if ~isempty(options.outputdir)
            options_out = options;
            options_out.head_vertex_ids = ROI_head.head_vertex_ids;
            tmp = struct('sensitivity_mat', sensitivity_mat,'coverage_mat', coverage_mat,'options',options_out);
            weight_cache.(matlab.lang.makeValidName(ROI_cortex.Label)) = tmp;
            save(fullfile(options.outputdir, 'weight_tables.mat'), 'weight_cache');
        end
        
    end

end

function [sensitivity_mat, coverage_mat] = compute_weights(fluence_volumes, head_vertices_coords, reference, options)

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
    V_hat = 1; % changer par une fonction permettant de le calculer
    delta_mu_a = 0.1;
    threshold = process_nst_extract_sensitivity_from_head_model('compute_threshold', p_thresh, act_vol, V_hat, delta_mu_a);
    
    bst_progress('start', 'Compute weight tables','Computing summed sensitivities of holder pairs...', 1, nHolders^2);
    for isrc=1:nHolders
        fluenceSrc = fluences(:,isrc);
        for idet=1:nHolders
            if holder_distances(isrc, idet) > options.sep_SD_min && holder_distances(isrc, idet) < options.sep_optode_max
                if ref(isrc,idet) ~=0 
                    fluenceDet = fluences(:,idet);
                    sensitivity = (fluenceSrc .* fluenceDet) /  ref(isrc,idet); 
                    
                    mat_sensitivity_idx(1,n_sensitivity_val) = isrc; mat_sensitivity_idx(2,n_sensitivity_val) = idet; mat_sensitivity_val(n_sensitivity_val) = sum(sensitivity) ;
                    n_sensitivity_val = n_sensitivity_val + 1;

                    coverage = sum(sensitivity > threshold) / length(sensitivity);

                    mat_coverage_idx(1, n_coverage_val) = isrc; mat_coverage_idx(2, n_coverage_val) = idet; mat_coverage_val(n_coverage_val) = coverage;
                    n_coverage_val = n_coverage_val + 1;
                end
            end    
        end
        bst_progress('inc', nHolders);
    end

    sensitivity_mat     = sparse(mat_sensitivity_idx(1,1:n_sensitivity_val-1),mat_sensitivity_idx(2,1:n_sensitivity_val-1), mat_sensitivity_val(1:n_sensitivity_val-1),nHolders, nHolders); 
    coverage_mat        = sparse(mat_coverage_idx(1,1:n_coverage_val-1), mat_coverage_idx(2,1:n_coverage_val-1), mat_coverage_val(1:n_coverage_val-1), nHolders, nHolders);
    
    bst_progress('stop');  
end

function [montage_pairs_simple, montage_sensitivity_simple, montage_coverage_simple, montage_pairs, montage_sensitivity, montage_coverage] = compute_optimal_montage(head_vertices_coords, options)
    
    %======================================================================
    % 1) Compute OM by maximizing sensitivity only
    %======================================================================
    
    % Define the cplex problem
    [prob, options] = define_prob(options.sensitivity_mat, head_vertices_coords, options);
   
    cplex=Cplex(prob);
    cplex.Model.sense = 'maximize';
    cplex.Param.timelimit.Cur=300;

    % Delete clone[number].log files created by Cplex
    cplex.Param.output.clonelog.Cur = 0;

    % Progress bar
    bst_progress('start', 'Optimization','Running optimization with Cplex. May take several minutes (see matlab console) ...');
    
    %======================================================================
    % initial condition : find sources whose pairs have maximum energy then
    % complete with detectors
    % Uncomment the line to use initial condition
    % cplex = init_solution(cplex, options.sensitivity_mat, options.nH, options.nb_sources, options.nb_detectors);
    %======================================================================
  
    results = cplex.solve();
    if ~isfield(results, 'x')
        bst_error(['OM computation failed  at Cplex step:', results.statusstring]);
        return;
    end

    bst_progress('stop');

    % Calculation of montage_pairs matrix, montage_sensitivity and montage_coverage vector
    [montage_pairs_simple, montage_sensitivity_simple, montage_coverage_simple] = montage_pairs_and_weight(results,options);

    %======================================================================
    % 2) Compute OM by maximizing sensitivity and coverage
    %======================================================================
    % Define the cplex problem

    lambda1 = 1/sum(montage_sensitivity_simple);
    lambda2 = 2; % options.lambda 
    wt = lambda1 * options.sensitivity_mat  + lambda2 * options.coverage_mat;
    [prob, options] = define_prob(wt, head_vertices_coords, options);

    cplex=Cplex(prob);
    cplex.Model.sense = 'maximize';
    
    % Faire varier le temps pour verifier
    cplex.Param.timelimit.Cur=300;

    % Delete clone[number].log files created by Cplex
    cplex.Param.output.clonelog.Cur = 0;

    % Progress bar
    bst_progress('start', 'Optimization','Running optimization with Cplex. May take several minutes (see matlab console) ...');
    
    results = cplex.solve();
    if ~isfield(results, 'x')
        bst_error(['OM computation failed  at Cplex step:', results.statusstring]);
        return;
    end

    bst_progress('stop');

    % Calculation of montage_pairs matrix, montage_sensitivity and montage_coverage vector
    [montage_pairs, montage_sensitivity, montage_coverage] = montage_pairs_and_weight(results, options);

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

function [prob, options] = define_prob(weight_table, head_vertices_coords, options)
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
        
    %......................................................................
    % DEBUG FUNCTION : normalize weigth table
    % norm_w_table(weight_table);
    %......................................................................
    
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

function cplex = init_solution(cplex, weight_table, nH, nS, nD)
%@=========================================================================
% INITIAL SOLUTION FUNCTIONS (cplex)
%=========================================================================@
    
    xp_0 = zeros(nH, 1);
    yq_0 = zeros(nH, 1);
    [~,idx] = sort(weight_table(:), 'descend');
    % Seems to only needed if holder list is not the same size as weight_table
    
    [src_pos_idx, ~] = ind2sub(size(weight_table), idx(1:nS));
    
    src_weight_table = weight_table(src_pos_idx,:);
    [~, idx] = sort(src_weight_table(:), 'descend');
    
    [~, det_pos_idx] = ind2sub(size(src_weight_table), idx(1:nD));
    
    xp_0(src_pos_idx) = 1;
    yq_0(det_pos_idx) = 1;
    
    wq_V_O = zeros(nH, 1);
    for ii = 1:numel(yq_0)
        if yq_0(ii) == 1
            wq_V_O(ii) = sum(weight_table(src_pos_idx, ii));
        end
    end
    
    incumbent_x0 = full(sum(reshape(weight_table(src_pos_idx, det_pos_idx), [], 1)));
    fprintf('incumbent=%g\n', incumbent_x0)
    x0 = [xp_0 ; yq_0 ; wq_V_O];

    cplex.MipStart(1).name='optim';
    cplex.MipStart(1).effortlevel=1;
    cplex.MipStart(1).x=x0;
    cplex.MipStart(1).xindices=int32((1:numel(x0))');
end

function [montage_pairs, montage_sensitivity, montage_coverage] = montage_pairs_and_weight(results,options)
% @========================================================================
% montage_pairs_and_weight Calculation of montage pairs matrix and montage
% weight vector
% ========================================================================@
    
    x=results.x;
    x=round(x);
    isources = find(x(1:options.nH)==1);
    idetectors = find(x(options.nH+1:2*options.nH)==1);
    
    ipair = 1;
    
    % Memory management
    max_pairs = length(isources) * length(idetectors);
    montage_pairs = zeros(max_pairs, 2);
    montage_sensitivity = zeros(max_pairs, 1);
    montage_coverage = zeros(max_pairs, 1);
    
    for isrc = 1:length(isources)
        for idet = 1:length(idetectors)
            if options.holder_distances(isources(isrc), idetectors(idet)) > options.thresh_sep_optode_optode(1) && ...
                    options.holder_distances(isources(isrc), idetectors(idet)) < options.thresh_sep_optode_optode(2) && ...
                    full(options.sensitivity_mat(isources(isrc), idetectors(idet)))
                
                if ipair <= max_pairs
                    montage_pairs(ipair,:) = [isources(isrc) idetectors(idet)];
                    montage_sensitivity(ipair,:) = full(options.sensitivity_mat(isources(isrc), idetectors(idet)));
                    montage_coverage(ipair,:) = full(options.coverage_mat(isources(isrc), idetectors(idet)));

                    ipair = ipair + 1;
                else
                    warning('Memory management error : The variables are not correctly sized. ');   
                end
            end
        end
    end
    
    % Make sure the matrix is the right size
    montage_pairs = montage_pairs(1:ipair-1, :);
    montage_sensitivity = montage_sensitivity(1:ipair-1, :);
    montage_coverage = montage_coverage(1:ipair-1, :);
end

%==========================================================================
% DEBUG FUNCTIONS
%==========================================================================

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

function weight_table = norm_w_table(weight_table)
    weight_table = weight_table/max(weight_table(:));
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

function ChannelMat = create_channelMat_from_montage(montage_pairs, montage_sensitivity, montage_coverage, head_vertices_coords, wavelengths)

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
        
        disp(['Channel S', num2str(idx_src), 'D' num2str(idx_det) ...
            ' >>> Distance: ', num2str(round(nst_pdist(head_vertices_coords(ihead_vertex_src, :),head_vertices_coords(ihead_vertex_det, :)).*1000,1)), 'mm    ', ...
            'Sensitivity: ', num2str(round(montage_sensitivity(ipair,:),3)), '    ' ...
            'Coverage : ', num2str(round(montage_coverage(ipair,:),3))]);

        formatSpec = 'Channel S%dD%d >>> Distance: %4.1fmm    Sensitivity: %6.3f    Coverage: %.3f';

        % 2. Calculer les valeurs
        distance_mm = nst_pdist(head_vertices_coords(ihead_vertex_src, :), head_vertices_coords(ihead_vertex_det, :)) * 1000;
        sensitivity = montage_sensitivity(ipair, :);
        coverage = montage_coverage(ipair, :);
        
        % 3. Créer et afficher la chaîne formatée
        fprintf(formatSpec, idx_src, idx_det, distance_mm, sensitivity, coverage);

        
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

function weight_table = denoise_weight_table(weight_table, threshold)
    
    if nargin < 2 || isempty(threshold)
        nrqr = 5;
        threshold = median(weight_table(weight_table>0)) + nrqr* iqr(weight_table(weight_table>0));
    end

    weight_table_new = weight_table;
    weight_table_new(weight_table > threshold) = 0;

    figure;
    subplot(221)
    imagesc(full(weight_table))
    title('Before Removing supirous node')

    subplot(222)
    imagesc(full(weight_table_new))
    title('After Removing supirous node')

    subplot(223); histogram(weight_table(weight_table <= threshold)); 
    subplot(224); histogram(weight_table(weight_table > threshold));
    title(sprintf('Threshold %.f ', threshold))


    weight_table  = weight_table_new;
    
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