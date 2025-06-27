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
% Edouard Delaire (2021), Jean-Eudes Bornert (2025)
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
% Definition of the input accepted by this process
sProcess.InputTypes  = {'import'};
% Definition of the outputs of this process
sProcess.OutputTypes = {'import'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 0;

% Subject name
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
function OutputFile = Run(sProcess, ~) 
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

options = sProcess.options.fluencesCond.Value;
    
SubjectName =  options.SubjectName;
sProcess.options.subjectname.Value = SubjectName;
sSubject = bst_get('Subject',SubjectName);

if isempty(sSubject.iCortex) || isempty(sSubject.iScalp)
    bst_error('No available Cortex and Head surface for this subject.');
    return;
end    

sHead = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);
sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);

i_atlas_cortex = strcmp({sCortex.Atlas.Name},options.Atlas_cortex);
i_scout_cortex = strcmp({sCortex.Atlas(i_atlas_cortex).Scouts.Label}, options.ROI_cortex);
ROI_cortex = sCortex.Atlas(i_atlas_cortex).Scouts(i_scout_cortex);   

if isempty(options.Atlas_head) && isempty(options.ROI_head)
    
    cortex_to_scalp_extent = options.Extent;
    cortex_scout.sSubject = sSubject;
    cortex_scout.sScout   = ROI_cortex;

    [head_vertex_ids, ~, ~] = proj_cortex_scout_to_scalp(cortex_scout,cortex_to_scalp_extent.*0.01, 1); 
else    
    i_atlas_head = strcmp({sHead.Atlas.Name},options.Atlas_head);
    i_scout_head = strcmp({sHead.Atlas(i_atlas_head).Scouts.Label}, options.ROI_head);
    head_vertex_ids = sHead.Atlas(i_atlas_head).Scouts(i_scout_head).Vertices;  
    
exclude_scout = sHead.Atlas(i_atlas_head).Scouts(strcmp('FluenceExclude', {sHead.Atlas(i_atlas_head).Scouts.Label}));
    if ~isempty(exclude_scout)
        head_vertex_ids = setdiff(head_vertex_ids, exclude_scout.Vertices);
    end
    
end


condition_name = options.condition_name;
if isempty(condition_name)
    condition_name = 'planning_optimal_montage';
end

% Extract head point coordinates
head_vertices_coords = sHead.Vertices(head_vertex_ids, :);

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
sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
cubeSize = size(sMri.Cube);

options.sep_optode_min = options.sep_optode(1);
options.sep_optode_max  = options.sep_optode(2);
options.sep_SD_min = options.sepmin_SD;
options.cubeSize = cubeSize;
if options.sep_optode_min > options.sep_SD_min
    bst_error(sprintf('ERROR: The minimum distance between source and detector has to be larger than the minimum optodes distance'));
    return;
end


sROI = extract_scout_surface(sCortex, ROI_cortex);
sROI.Vertices = cs_convert(sMri, 'scs', 'voxel', sROI.Vertices');
    
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
iseg = find(strcmp({sSubject.Anatomy.Comment},segmentation_name));
if isempty(iseg)
    bst_error(sprintf('ERROR: Please import segmentation file as MRI and rename it as "%s"', segmentation_name));
    return;
end

GM_mask = process_nst_compute_voronoi('get_grey_matter_mask',sSubject.Anatomy(iseg).FileName);
voronoi_mask = (voronoi > -1) &  ~isnan(voronoi) & GM_mask & ismember(voronoi,ROI_cortex.Vertices);

voi               = zeros(cubeSize(1), cubeSize(2), cubeSize(3));
voi(voronoi_mask) = 1;

voi_flat    = voi(:);
vois        = sparse(voi_flat > 0);

if nnz(vois)==0
    bst_error('VOI mask is empty after intersecting with Voronoi and GM masks');
    return;
end

weight_table = [];
weight_cache = struct();

if exist(fullfile(options.outputdir , 'weight_tables.mat'), 'file')
    tmp = load (fullfile(options.outputdir, 'weight_tables.mat'));
    
    if isfield(tmp,'weight_cache') && ~isempty(tmp.weight_cache)
        weight_cache = tmp.weight_cache;
    end

    if options.exist_weight && isfield(weight_cache,  strrep(ROI_cortex.Label,' ','_'))
        tmp = weight_cache.(strrep(ROI_cortex.Label,' ','_'));    
        if isequal(tmp.options.head_vertex_ids,head_vertex_ids) && tmp.options.sep_SD_min == options.sep_SD_min &&  tmp.options.sep_optode_max == options.sep_optode_max   
            weight_table = tmp.weight_table;
        end    
    end    
end

if isempty(weight_table)
    % Compute weight table
    sparse_threshold = 1e-6;
    fluence_data_source = options.data_source;

    [fluence_volumes,reference] = process_nst_import_head_model('request_fluences', head_vertex_ids, ...
        sSubject.Anatomy(sSubject.iAnatomy).Comment, ...
        wavelengths, fluence_data_source, sparse_threshold,vois,cubeSize);
    if isempty(fluence_volumes)
        return;
   end
    
    weight_table = compute_weights(fluence_volumes,head_vertices_coords,reference, options);
    
    if(nnz(weight_table) == 0)
        bst_error(sprintf('Weight table is null for ROI: %s', ROI_cortex.Label));
        return
    end
    if ~isempty(options.outputdir)
        options_out = options;
        options_out.head_vertex_ids = head_vertex_ids;
        tmp = struct('weight_table',weight_table,'options',options_out);
        weight_cache.(strrep(ROI_cortex.Label,' ','_')) = tmp;
        save(fullfile(options.outputdir, 'weight_tables.mat'), 'weight_cache');
    end
end

% weight_table_new = weight_table;
% nrqr = 5;
% spurious = median(weight_table(weight_table>0)) + nrqr* iqr(weight_table(weight_table>0));
% weight_table_new(weight_table > spurious) = 0;
%     
% figure;
% subplot(221)
% imagesc(full(weight_table))
% title('Before Removing supirous node')
% subplot(222)
% imagesc(full(weight_table_new))
% title('After Removing supirous node')
% subplot(223); hist(weight_table(weight_table < spurious)); subplot(224); hist(weight_table(weight_table > spurious));
% title(sprintf('Threshold %.f (M interquartile %d)',spurious,nrqr))
% 
% weight_table = weight_table_new;

options.weight_tables = weight_table;
[montage_pairs,montage_weight] = compute_optimal_montage(head_vertices_coords,options);

% montage_pairs contains head vertex indexes -> remap to 1-based consecutive indexes
src_indexes = zeros(max(montage_pairs(:, 1)), 1);
det_indexes = zeros(max(montage_pairs(:, 2)), 1);


% Forge channels from head points pairs & wavelengths
nChannels = size(montage_pairs, 1) * length(wavelengths);
ChannelMat = db_template('channelmat');
ChannelMat.Comment = 'NIRS-BRS channels';
ChannelMat.Channel = repmat(db_template('channeldesc'), [1, nChannels]);
ChannelMat.Nirs.Wavelengths = wavelengths;
iChan = 1;
det_next_idx = 1;
src_next_idx = 1;

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
    
    disp(['Channel S', num2str(idx_src), 'D' num2str(idx_det) ' >>> Distance: ',...
        num2str(round(nst_pdist(head_vertices_coords(ihead_vertex_src, :),head_vertices_coords(ihead_vertex_det, :)).*1000,1)), 'mm    ', 'Weight: ',...
        num2str(round(montage_weight(ipair,:),3))]);
    
    for iwl=1:length(wavelengths)
        
        ChannelMat.Channel(iChan).Name    = sprintf('S%dD%dWL%d', idx_src, idx_det, ...
        ChannelMat.Nirs.Wavelengths(iwl));
        ChannelMat.Channel(iChan).Type    = 'NIRS';
        
        ChannelMat.Channel(iChan).Loc(:,1)  = head_vertices_coords(ihead_vertex_src, :);
        ChannelMat.Channel(iChan).Loc(:,2)  = head_vertices_coords(ihead_vertex_det, :);
        ChannelMat.Channel(iChan).Orient  = [];
        ChannelMat.Channel(iChan).Weight  = 1;
        ChannelMat.Channel(iChan).Comment = [];
        ChannelMat.Channel(iChan).Group = sprintf('WL%d', round(ChannelMat.Nirs.Wavelengths(iwl)));
        
        iChan = iChan + 1;
    end
end

[sSubjStudies, ~] = bst_get('StudyWithSubject', sSubject.FileName, 'intra_subject', 'default_study');
condition_name    = file_unique(condition_name, {sSubjStudies.Name}, 1);

iStudy = db_add_condition(sSubject.Name, condition_name);
sStudy = bst_get('Study', iStudy);

% Save channel definition
[~, iChannelStudy] = bst_get('ChannelForStudy', iStudy);
db_set_channel(iChannelStudy, ChannelMat, 1, 0);
    

separations = process_nst_separations('Compute',ChannelMat.Channel) * 100; %convert to cm
% Save time-series data
sDataOut              = db_template('data');
sDataOut.F            = separations;
sDataOut.Comment      = 'Separations';
sDataOut.ChannelFlag  = ones(length(separations),1);
sDataOut.Time         = (1);
sDataOut.DataType     = 'recordings';
sDataOut.nAvg         = 1;
sDataOut.DisplayUnits = 'cm';

% Generate a new file name in the same folder
OutputFile{1} = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_chan_dist');
sDataOut.FileName = OutputFile{1} ;
bst_save(OutputFile{1} , sDataOut, 'v7');
% Register in database
db_add_data(iStudy, OutputFile{1} , sDataOut);
    

end

function weight_tables = compute_weights(fluence_volumes,head_vertices_coords,reference, options)
    holder_distances = nst_pdist(head_vertices_coords, head_vertices_coords).*1000; % mm
    nHolders = size(head_vertices_coords, 1);
    iwl = 1;

    mat_idx = zeros(2,nHolders);
    mat_val = zeros(1,nHolders);
    n_val = 1;
    
    bst_progress('start', 'Compute weights','Preparation of fluences...', 1, 2);

        
    fluences =  zeros(size(full(fluence_volumes{1}{iwl}),1),nHolders);
    for isrc=1:nHolders
          fluences(:,isrc) =  full(fluence_volumes{isrc}{iwl});
    end
    bst_progress('inc', 1);
    
    ref = zeros(nHolders,nHolders);
    for isrc=1:nHolders
        for idet=1:nHolders
            if holder_distances(isrc, idet) > options.sep_SD_min && holder_distances(isrc, idet)< options.sep_optode_max                   
                ref(isrc,idet) = full(reference{isrc}{iwl}(idet));  
            end
        end    
    end
    
    bst_progress('inc', 1);

    bst_progress('start', 'Compute weights','Computing summed sensitivities of holder pairs...', 1, nHolders^2);
    for isrc=1:nHolders
        fluenceSrc = fluences(:,isrc);
        for idet=1:nHolders
            if holder_distances(isrc, idet) > options.sep_SD_min && holder_distances(isrc, idet)< options.sep_optode_max
                %A=normFactor*fluenceSrc.*fluenceDet./diff_mask(idx_vox);
                fluenceDet = fluences(:,idet);

                if ref(isrc,idet) ~=0 
                    % sensitivity = fluenceSrc .* fluenceDet ./ fluence_volumes{isrc}{iwl}(ref_det_pos); % Asymmetry
                    % ED: fluenceSrc' * fluenceDet is slighty faster than sum ( fluenceSrc .* fluenceDet )
                    sensitivity = fluenceSrc' * fluenceDet; % ./ fluence_volumes{isrc}{iwl}(ref_det_pos); % Asymmetry
                    
                    mat_idx(1,n_val) = isrc;mat_idx(2,n_val) = idet; mat_val(n_val) = sensitivity / ref(isrc,idet);
                    n_val = n_val +1;
                end
            end    
        end
        bst_progress('inc', nHolders);
    end
    weight_tables = sparse(mat_idx(1,1:n_val-1),mat_idx(2,1:n_val-1), mat_val(1:n_val-1),nHolders,nHolders); 
    bst_progress('stop');
end


function [montage_pairs,montage_weight] = compute_optimal_montage(head_vertices_coords,options)
    
    %Define the cplex problem
    [prob, options] = define_prob_simple(head_vertices_coords, options);
   
    cplex=Cplex(prob);
    cplex.Model.sense = 'maximize';
    cplex.Param.timelimit.Cur=300;

    % Delete clone[number].log files created by Cplex
    cplex.Param.output.clonelog.Cur = 0;

    %Progress bar
    bst_progress('start', 'Optimization','Running optimization with Cplex. May take several minutes (see matlab console) ...', 1, 1000);
    
    %======================================================================
    % initial condition : find sources whose pairs have maximum energy then
    % complete with detectors
    % Uncomment the line to use initial condition
    % cplex = init_solution(cplex, options.weight_tables, options.nH, options.nb_sources, options.nb_detectors);
    %======================================================================
  
    results = cplex.solve;
    if ~isfield(results, 'x')
        bst_error(['OM computation failed  at Cplex step:', results.statusstring]);
        return;
    end
   
    addConstraint = 0;
    if addConstraint
       
    end

    bst_progress('stop');

    %Calculation of montage_pairs matrix and montage_weight vector
    [montage_pairs, montage_weight] = montage_pairs_and_weight(results, options.nH, options.holder_distances, ...
        options.thresh_sep_optode_optode, options.weight_tables);
end


function sSurfNew = extract_scout_surface(sSurf, sScouts)

    iRemoveVert = setdiff(1:size(sSurf.Vertices,1), [sScouts.Vertices]);
    tag = 'OM_scout_extract';
    
    % Unload everything
    bst_memory('UnloadAll', 'Forced');
    
    % === REMOVE VERTICES FROM SURFACE FILE ===
    % Remove vertices
    [Vertices, Faces, Atlas] = tess_remove_vert(sSurf.Vertices, sSurf.Faces, iRemoveVert, sSurf.Atlas);
    % Remove the handles of the scouts
    for iAtlas = 1:length(Atlas)
        for is = 1:length(Atlas(iAtlas).Scouts)
            Atlas(iAtlas).Scouts(is).Handles = [];
        end
        %         if isfield(Atlas(iAtlas).Scouts, 'Handles');
        %             Atlas(iAtlas).Scouts = rmfield(Atlas(iAtlas).Scouts, 'Handles');
        %         end
    end
    % Build new surface
    sSurfNew = db_template('surfacemat');
    sSurfNew.Comment  = [sSurf.Comment '_' tag];
    sSurfNew.Vertices = Vertices;
    sSurfNew.Faces    = Faces;
    sSurfNew.Atlas    = Atlas;
    sSurfNew.iAtlas   = sSurf.iAtlas;
    
    % % === SAVE NEW FILE ===
    % % Output filename
    % NewTessFile = strrep(file_fullpath(sSurf.FileName), '.mat', ['_' tag '.mat']);
    % NewTessFile = file_unique(NewTessFile);
    % % Save file back
    % bst_save(NewTessFile, sSurfNew, 'v7');
    % % Get subject
    % [sSubject, iSubject] = bst_get('SurfaceFile', sSurf.FileName);
    % % Register this file in Brainstorm database
    % db_add_surface(iSubject, NewTessFile, sSurfNew.Comment);

end


function [head_vertices, sHead, sSubject] = proj_cortex_scout_to_scalp(cortex_scout, extent_m, save_in_db)

    if nargin < 3
        save_in_db = 0;
    end
    sSubject = cortex_scout.sSubject;
    sHead = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);
    sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);
    dis2head = nst_pdist(sHead.Vertices, sCortex.Vertices(cortex_scout.sScout.Vertices,:));
    head_vertices = find(min(dis2head,[],2) < extent_m); 
    
    % TODO: properly select atlas
    exclude_scout = sHead.Atlas.Scouts(strcmp('FluenceExclude', {sHead.Atlas.Scouts.Label}));
    if ~isempty(exclude_scout)
        head_vertices = setdiff(head_vertices, exclude_scout.Vertices);
    end
    
    limiting_scout = sHead.Atlas.Scouts(strcmp('FluenceRegion', {sHead.Atlas.Scouts.Label}));
    if ~isempty(limiting_scout)
        head_vertices = intersect(head_vertices, limiting_scout.Vertices);
    end
    
    if save_in_db && ...,
       ~any(strcmp(['From cortical ' cortex_scout.sScout.Label '(' num2str(extent_m*100) ' cm)']...,
        ,{sHead.Atlas.Scouts.Label}))
        scout_idx = size(sHead.Atlas.Scouts,2) + 1;
        sHead.Atlas.Scouts(scout_idx) = db_template('Scout');
        sHead.Atlas.Scouts(scout_idx).Vertices = head_vertices';
        sHead.Atlas.Scouts(scout_idx).Seed = head_vertices(1);
        sHead.Atlas.Scouts(scout_idx).Color = [0,0,0];
        sHead.Atlas.Scouts(scout_idx).Label = ['From cortical ' cortex_scout.sScout.Label ...
                                               '(' num2str(extent_m*100) ' cm)'];
        bst_save(file_fullpath(sSubject.Surface(sSubject.iScalp).FileName), sHead, 'v7');
        db_save();
    end

end

function [prob, options] = define_prob_simple(head_vertices_coords, options)
% @========================================================================
% define_prob_simple Initializes the problem
% Added in options by function : holder_distances, nH, thresh_sep_optode_optode 
% ========================================================================@
    
    holder_distances = nst_pdist(head_vertices_coords, head_vertices_coords).*1000; % mm

    nHolders = size(head_vertices_coords, 1);
    weight_table = options.weight_tables;
    
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
    
    
    bst_progress('start', 'Optimization setup', 'Setting up optimization...', 1, 4);
    
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
    
    bst_progress('inc', 1);
    
    Aeq = [Aeq_1 ; Aeq_2];
    E = [E_1 ; E_2];
    
    Aineq = [Aineq_1 ; Aineq_2 ; Aineq_3 ; Aineq_4,; Aineq_5];
    I = [I_1 ; I_2 ; I_3 ; I_4 ; I_5];
    
    f  = [zeros(1, 2*nH) ones(1, nH)]';
    lb = [];%[zeros(1,2*nH) zeros(1,nH)]';
    ub = [];%[ones(1,2*nH)  realmax.*ones(1,nH)]';
    ctype = [repmat('B', 1, 2*nH) repmat('S', 1, nH)];
    
    %Cplex optimisation
    prob = cplexcreateprob('cplexmilp');
    prob.f = f;
    prob.lb = lb; prob.ub = ub;
    prob.ctype = ctype;
    prob.Aineq = Aineq; prob.bineq = I;
    prob.Aeq = Aeq; prob.beq = E;
    prob.x0 = [];
    prob.options = [];
    bst_progress('inc', 1);
    bst_progress('stop');
    
    options.holder_distances = holder_distances;
    options.thresh_sep_optode_optode = thresh_sep_optode_optode; 
    options.nH = nH;
end

function [A, E] = add_constraint_nSrc(nVar, nH, nS)
% @========================================================================
% add_constraint_nSrc Adds the constraint for the total number of sources.
%   See equation 11 : Sum(xp) = nSrc
% ========================================================================@

    A = zeros(1, nVar);
    A(1, 1:nH) = 1;

    E = nS;
end

%==========================================================================
% FUNCTIONS USED BY "compute_optimal_montage" FUNCTION
%==========================================================================

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
    bst_progress('inc', 1);
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
    
    bst_progress('inc', 1);
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

function [montage_pairs, montage_weight] = montage_pairs_and_weight(results, nH, holder_distances, thresh_sep_optode_optode, weight_table)
% @========================================================================
% montage_pairs_and_weight Calculation of montage pairs matrix and montage
% weight vector
% ========================================================================@

    x=results.x;
    x=round(x);
    isources = find(x(1:nH)==1);
    idetectors = find(x(nH+1:2*nH)==1);
    
    ipair = 1;
    
    % Memory management
    max_pairs = length(isources) * length(idetectors) -1;
    montage_pairs = zeros(max_pairs, 2);
    montage_weight = zeros(max_pairs, 1);
    
    for isrc = 1:length(isources)
        for idet = 1:length(idetectors)
            if holder_distances(isources(isrc), idetectors(idet)) > thresh_sep_optode_optode(1) && ...
                    holder_distances(isources(isrc), idetectors(idet)) < thresh_sep_optode_optode(2) && ...
                    full(weight_table(isources(isrc), idetectors(idet)))
                
                if ipair <= max_pairs
                    montage_pairs(ipair,:) = [isources(isrc) idetectors(idet)];
                    montage_weight(ipair,:) = full(weight_table(isources(isrc), idetectors(idet)));
                    ipair = ipair + 1;
                else
                    warning('Memory management error : The variables are not correctly sized. ');   
                end
            end
        end
    end
    
    % Make sure the matrix is the right size
    montage_pairs = montage_pairs(1:ipair-1, :);
    montage_weight = montage_weight(1:ipair-1, :);
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