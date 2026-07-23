function varargout = process_nst_compare_montage( varargin )

% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2016 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
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
% Authors: Edouard Delaire, Jean-Eudes Bornert 2025

eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
    % Description the process
    sProcess.Comment     = 'Compare montages';
    sProcess.FileTag     = '';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = {'NIRS', 'Sources'};
    sProcess.Index       = 1408;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'data', 'raw'}; 
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    
    % === CLUSTERS
    sProcess.options.scouts.Comment = '';
    sProcess.options.scouts.Type    = 'scout';
    sProcess.options.scouts.Value   = {};

    % Definition of the options
    sProcess.options.thresh_dis2cortex.Comment = 'Reconstruction Field of view (distance to montage border)';
    sProcess.options.thresh_dis2cortex.Type    = 'value';
    sProcess.options.thresh_dis2cortex.Value   = {10, 'cm',2};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end

%% ===== RUN ======

function OutputFile = Run(sProcess, sInput)
    OutputFile = '';
    OutputFiles = {};

    % Load subject
    sSubject = bst_get('Subject', sInput.SubjectName);

    % Load ROI
    sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);
    sHead   = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);
    ChannelMat = in_bst_channel(sInput(1).ChannelFile);

    % Load Head model
    sStudy = bst_get('Study', sInput.iStudy);
    HeadModelFile   = sStudy.HeadModel(sStudy.iHeadModel).FileName;
    nirs_head_model = in_bst_headmodel(HeadModelFile, 1);

    % Compute ground truth
    ROI         = sProcess.options.scouts.Value;
    iAtlas      = find(strcmp( {sCortex.Atlas.Name},ROI{1}));
    iRois       = cellfun(@(x)find(strcmp( {sCortex.Atlas(iAtlas).Scouts.Label},x)),   ROI{2});
    
    assert(length(iRois) == 1, 'Please select only one ROI');
    ROI_select = sCortex.Atlas(iAtlas).Scouts(iRois);
    sGroundTruth = false(size(sCortex.Vertices,1), 1);
    sGroundTruth(ROI_select.Vertices) = 1;

    % Get only vertex close to the surface in the ROI
    distance = process_nst_cpt_cortex_to_head_distance('Compute', sCortex, sHead)*1000;
    threshold_distance = 20;
    idx_exclude = find(distance > threshold_distance);

    sGroundTruth(idx_exclude) = 0;
    fprintf('Removing %d vertex from the ground truth (too deep) \n', length(idx_exclude));

    % Load Head Model
    sStudy = bst_get('Study', sInput.iStudy);
    if isempty(sStudy.iHeadModel)
        bst_error('No head model found. Consider process "Compute head model from fluence"');
        return;
    end
    
    head_model = in_bst_headmodel(sStudy.HeadModel(sStudy.iHeadModel).FileName, 1);
    G = head_model.Gain;

    coverage = compute_coverage(sSubject, G);
    metrics_montage = compute_metrics_montage(G, coverage, sGroundTruth);
    
    [data_simul, groundTruth, valid_nodes, ChannelFlag, event, PSF, metrics_simul] = compute_metrics_simulation(sCortex, ROI_select, ChannelMat, nirs_head_model, sProcess.options);

    % === SAVE DATA ===
    % 1. Save Simulated Data
    sDataOut = db_template('data');
    sDataOut.F            = data_simul; 
    sDataOut.Comment      = 'Simulated data';
    sDataOut.ChannelFlag  = ChannelFlag;
    sDataOut.Time         = [1];
    sDataOut.DataType     = 'recordings'; 
    sDataOut.nAvg         = 1;
    sDataOut.Events       = [];
    sDataOut = bst_history('add', sDataOut, 'process', sProcess.Comment);
    sDataOut.DisplayUnits = 'delta OD';
    
    OutputFile_data = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_sim');
    sDataOut.FileName = file_short(OutputFile_data);
    bst_save(OutputFile_data, sDataOut, 'v7');
    db_add_data(sInput.iStudy, OutputFile_data, sDataOut);
    OutputFiles{end+1} = OutputFile_data;

    ResultsMat = db_template('resultsmat');
    ResultsMat.Comment       = 'Ground Truth';
    ResultsMat.DataFile      = file_short(OutputFile_data);
    ResultsMat.Function      = '';
    ResultsMat.Time          = [1];
    ResultsMat.ImageGridAmp  = groundTruth;
    ResultsMat.ChannelFlag   = [];
    ResultsMat.GoodChannel   = [];
    ResultsMat.DisplayUnits  = 'delta OD';
    ResultsMat.SurfaceFile   = nirs_head_model.SurfaceFile;
    % % Save new file structure
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'results_ground_truth_simul');
    bst_save(OutputFile, ResultsMat, 'v6');
    % Update database
    db_add_data(sInput.iStudy, OutputFile, ResultsMat);

    ResultsMat = db_template('resultsmat');
    ResultsMat.Comment       = 'PSF';
    ResultsMat.DataFile      = file_short(OutputFile_data);
    ResultsMat.Function      = '';
    ResultsMat.Time          = [1];
    ResultsMat.ImageGridAmp  = PSF;
    ResultsMat.ChannelFlag   = [];
    ResultsMat.GoodChannel   = [];
    ResultsMat.DisplayUnits  = 'delta OD';
    ResultsMat.SurfaceFile   = nirs_head_model.SurfaceFile;
    % % Save new file structure
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'results_PSF_simul');
    bst_save(OutputFile, ResultsMat, 'v6');
    % Update database
    db_add_data(sInput.iStudy, OutputFile, ResultsMat);
    
    metric_table = [struct2table(metrics_montage), struct2table(metrics_simul)];
    % name of each column
    metric_table.Properties.VariableNames = {'Total sensitivity (mm)', 'Coverage (%)', 'Overlap (# channels)', 'DLE (mm)', 'SD (mm)', 'AUC (%)'};
    metric_table.Properties.RowNames = {sInput.Condition};
    OutputFile = save_table(metric_table, sInput.SubjectName, sInput.Condition, 'Metrics results');

end

%% Coverage, overlap, sensitivity

function [coverage_channel]  = compute_coverage(sSubject, G)

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
    threshold = process_nst_extract_sensitivity_from_head_model('compute_threshold', p_thresh, act_vol, median_voronoi_volume, delta_mu_a);

    coverage_channel  = G  > threshold;
end

function metrics = compute_metrics_montage(sensitivity, coverage_channel, ROI)

    metrics = struct();

    % 1. Total sensitivity
    metrics.total_sensitivity = sum(sum(sensitivity(:, ROI)));

    %2. Overlap measure - % coverage of the ROI
    overlap = sum(coverage_channel,  1)' ;
    nVertex = sum(ROI);
    nVertexCovered = sum(ROI == 1 &  overlap >= 1);
    metrics.coverage = 100 * nVertexCovered / nVertex;
    
    %3. Overlap. Median number of overlap in the ROI
    metrics.nOverlap = median(overlap(ROI));

end

%% AUC, SD, DLE

function [data_simul, groundTruth, valid_nodes, ChannelFlag, event, PSF, metrics] = compute_metrics_simulation(sCortex, ROI_select, ChannelMat, nirs_head_model, OPTIONS)

    % Set up activation struct
    activation = struct(); 
    activation.Label    = ROI_select(1).Label;
    activation.Vertices = ROI_select(1).Vertices;
    activation.ampmode  = 'unif';
    activation.options  = get_default_options([1], 'oscilation');
    activation.options.SNR = 1;

    % Run simulation
    [data_simul, groundTruth, valid_nodes, ChannelFlag, event, PSF] =  simulNirs(sCortex, nirs_head_model, activation, ChannelMat, OPTIONS);
    [metrics] = compute_metric(sCortex, groundTruth, PSF);
end

function sFile = save_table(t, subject_name, condition, comment, extra, displayUnits)
    % Save a table as a Brainstorm matrix item using the 2D Description format
    if nargin < 5 
        extra = struct();
    end
    if nargin < 6
        displayUnits = '';
    end
    
    MatNew = db_template('matrixmat');
    MatNew.Value = table2array(t);
    [nRows, nCols] = size(MatNew.Value);
    
    MatNew.Std = [];
    MatNew.Time = [];
    MatNew.DisplayUnits = displayUnits;
    MatNew.Comment = comment;
    
    if ~isempty(t.Properties.RowNames)
        rowNames = t.Properties.RowNames(:);
    else
        rowNames = arrayfun(@(n) sprintf('row_%d', n), 1:nRows, 'UniformOutput', 0)';
    end
    colNames = t.Properties.VariableNames(:)';
    
    MatNew.Description = cell(nRows, nCols);
    MatNew.Description{1, 1} = {rowNames{1}, colNames{1}};
    
    if nCols > 1
        MatNew.Description(1, 2:nCols) = colNames(2:nCols);
    end
    if nRows > 1
        MatNew.Description(2:nRows, 1) = rowNames(2:nRows);
    end
    
    sSubject = bst_get('Subject', subject_name, 1);
    if isempty(sSubject)
        db_add_subject(subject_name, []);
    end
    
    [sStudy, iStudy] = bst_get('StudyWithCondition', [subject_name '/' condition]);
    if isempty(sStudy) 
        iStudy = db_add_condition(subject_name, condition);
        sStudy = bst_get('Study', iStudy);
    end
    
    extra_fields = fieldnames(extra);
    for ifield = 1:length(extra_fields)
        MatNew.(extra_fields{ifield}) = extra.(extra_fields{ifield});
    end
    
    sFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'matrix_table');
    bst_save(sFile, MatNew, 'v6');
    db_add_data(iStudy, sFile, MatNew);
end

function [data_simul, groundTruth, valid_nodes, ChannelFlag, event, PSF] = simulNirs(sCortex, head_model, activation, ChannelMat, OPTIONS)
    % Load head model
    iwl = 1;
    swl = ['WL' num2str(ChannelMat.Nirs.Wavelengths(iwl))];
    selected_chans = strcmpi({ChannelMat.Channel.Group}, swl);

    thresh_dis2cortex   = OPTIONS.thresh_dis2cortex.Value{1} .* 0.01;
    valid_nodes         = nst_headmodel_get_FOV(ChannelMat, sCortex, thresh_dis2cortex);
    gain                = head_model.Gain(selected_chans , valid_nodes);

    Time       = [1];  

    [y, event] = simulate_oscilations(Time, activation.options);

    nodes = zeros(1, size(sCortex.Vertices, 1));
    nodes(activation.Vertices) = 1; 

    nodes = nodes(valid_nodes);
    activated_vertices = find(nodes);

    data_cortex = zeros(length(valid_nodes), length(Time));
    data_cortex(activated_vertices, :) = repmat(y, length(activated_vertices), 1); 
    
    data_head = gain * data_cortex;

    % Final simulated data (Noise + Scaled Signal)
    data_simul =  data_head;

    groundTruth = zeros(size(sCortex.Vertices, 1), length(Time));
    groundTruth(valid_nodes, :) = data_cortex; 

    ChannelFlag = -1 * ones(size(data_simul, 1), 1);
    ChannelFlag(selected_chans, :) = 1;
    
    % Invert Operator (G) using standard Minimum Norm
    L = gain;
    lambda = 0.1 * trace(L*L') / size(L,1); % Simple regularization
    G = L' / (L * L' + lambda * eye(size(L, 1))); 
    
    % PSF of the ROI (R = G * L)
    PSF_ROI = G * L * nodes'; 
    
    PSF = zeros(size(groundTruth));
    PSF(valid_nodes) = PSF_ROI;
end

function [results] = compute_metric(sCortex, groundTruth, PSF)
    
    results = struct();
    idx_truth = find(groundTruth);
    
    % Estimate distance from all the vertex to the ROI, in milimeter
    distances = 1000 * min(nst_pdist(sCortex.Vertices,sCortex.Vertices(idx_truth,:)),[],2);

    % 1. DLE --  Distance from the max vertex to the ground truth
    [~, maxVert] = max(abs(PSF));
    results.DLE = distances(maxVert);

    % 2. Spatial Dispersion
    results.SD = sqrt(sum(distances.^2 .* PSF.^2) / sum(PSF.^2));
    
    % prepare ROC
    ROC_Struct = prepare_ROC(sCortex);

    % 3. AUC 
    [~, Res_close_summary, ~] =  Compute_ALL_AUC_global(0, ...
                                                         groundTruth, PSF, 1, ...
                                                         ROC_Struct.VoisinsOA, ...
                                                         ROC_Struct.mycluster, ...
                                                         ROC_Struct.nb_resampling,...
                                                         ROC_Struct.ordreVoisinage,...
                                                         ROC_Struct.thresholds, ...
                                                         []);
    

    results.auc_close       = 100*Res_close_summary.AUC_mean;
end

function [output] = prepare_ROC(sCortex)
    
    output = struct();
    output.VoisinsOA   = adj2Voisins(sCortex.VertConn);
    output.nClusters   = 30;
    output.isRandom    = 1;
    output.VERBOSE     = 0;
    output.Labels      = tess_cluster(sCortex.VertConn, output.nClusters, output.isRandom, output.VERBOSE);
    
    mycluster = cell(1,output.nClusters);
    for i = 1:output.nClusters
        mycluster{i} = find(output.Labels == i);
    end
    output.mycluster = mycluster;

    output.nb_resampling    = 100;
    output.ordreVoisinage   = 5;
    output.thresholds = linspace(0, 1, 100);
end

function VoisinsOA = adj2Voisins(adj)
% Convert the adjacency matrix 'adj' to 'VoisinsOA' neighbor cell vector
len = length(adj);

VoisinsOA = cell(12,len);
bst_progress('start', 'AUC computation', 'Please wait...', 1, 12);
for iScale = 1:12
    adj_i = logical(adj^iScale);
    adj_i(logical(eye(size(adj_i)))) = 0;
    for iSource = 1:len
       VoisinsOA{iScale,iSource} = find(adj_i(iSource,:));
    end
    bst_progress('set', iScale);
end
bst_progress('stop');

end

%% ===== AUC COMPUTATION =====
function Jfictif = AUC_create_JFictif_cluster(Jmeas, clusters, support, nbvertices, nTarget)
    %    Generation des sources fictives dans repartie dans le clustering !
    %    (pour le calcul d AUC_spurious)


    Isources_fictives_selected = [];
    IndiceCluster = randperm(length(clusters));
    i = 0;
    indice_decroissant = 1; % on ajoutre d abor le plus actif dans chaque cluster
    % puis le second plus actif ... et ainsi de suite
    
    while (length(Isources_fictives_selected) < nTarget)
        i = i+1;
        if (i > length(clusters))
            % on est deja passe une fois par chaque cluster
            % on prendra donc le 2nd point le plus actif
            i=1;
            indice_decroissant = indice_decroissant+1;
        end
        SelectCluster = IndiceCluster(i);
        
        [Jsorted, Indices_sorted ] = sort(Jmeas(clusters{SelectCluster}), 'descend');
        
        %Indice_source_in_cluster = randperm(length(clusters{SelectCluster}));
        %Candidat =  clusters{SelectCluster}(Indice_source_in_cluster(1));
        
        if (indice_decroissant  <= length(Indices_sorted))
            Candidat = clusters{SelectCluster}(Indices_sorted(indice_decroissant));
            
            if (~isempty(intersect(Candidat, support)) & ...
                    isempty(intersect(Candidat, Isources_fictives_selected))) 
                % Candidat bien dans le support fictif et pas deja pris
                Isources_fictives_selected = [Isources_fictives_selected Candidat];
            end
        end % on verifie que le point d indice indice_decroissant est toujours dans le clusters 
        
    end % end iwhile


    Jfictif = zeros(nbvertices,1);
    Jfictif(Isources_fictives_selected) = 1;
end

function Jfictif = AUC_create_JFictif_random(support, nbvertices, nTarget)
    IndicesSourcesF     = randperm(length(support));
    
    
    if size(IndicesSourcesF,2) > nTarget
        Isources_fictives_selected = support(IndicesSourcesF(1:nTarget));
    else
        Isources_fictives_selected = support(IndicesSourcesF(1:size(IndicesSourcesF,2)));
    end
    
    
    Jfictif = zeros(nbvertices,1);
    Jfictif(Isources_fictives_selected) = 1;
end

function [AUC] = ComputeAUC(Vec1,Vec2)
% ***********************************************************************
% [AUC] = ComputeAUC(Vec1,Vec2)
%  compute area under the curve using the trapezoids method
% ***********************************************************************
% Inputs : 
%   Vec1 : list of abscissae
%   Vec2 : list of ordinates
%
% ***********************************************************************
% Outputs : 
%    AUC : area under the curve
% ***********************************************************************
% C. Grova - Montreal Neurological Institute - 09 02 2004 
% ***********************************************************************

    if (length(Vec1) ~= length(Vec2))
        disp('Both vectors should have the same length');
        exit;
    end
    
    [X, I]  = sort(Vec1);
    Y       =  Vec2(I); 
    
    AUC = trapz(X,Y);
end

function [Res_summary,Res_close_summary,Res_far_summary  ] = Compute_ALL_AUC_global(Mesh, Jtheo, J, peak,  VoisinsOA, clusters, nb_resampling, ordreVoisinage, thresholds, area)
% ***********************************************************************
% [AUC,AUC_close,AUC_far] = Compute_ALL_AUC_global(Mesh, Jtheo, Jmeas, peak,  VoisinsOA, clusters, ordreVois, visu)
% ***********************************************************************
% Inputs : 
% ***********************************************************************
% Outputs : 
% ***********************************************************************
% C. Grova - Montreal Neurological Institute - 21 03 2005 
% ***********************************************************************


if nargin < 9
    thresholds = [];
end
if nargin < 10
    area = [];
end

Jmeas = abs(J(:,peak)) / max(abs(J(:,peak))); % normalisation entre 0 et 1 du vecteur resultat a tester


for i = 1:nb_resampling
    [Res_close(i),Res_far(i)]   = Compute_AUC_global(Jtheo, Jmeas,  VoisinsOA, clusters, ordreVoisinage, thresholds, area);
end

Res_close_summary   = struct();
Res_far_summary     = struct();
Res_summary         = struct();

fields = fieldnames(Res_close);

for i = 1:length(fields)
    data_close =  cell2mat({Res_close.(fields{i})});
    if size(data_close,2) == nb_resampling
        Res_close_summary.(sprintf('%s_mean', fields{i}))   = mean(data_close,2);
        Res_close_summary.(sprintf('%s_std', fields{i}))    = std(data_close,[],2);
    else
        Res_close_summary.(sprintf('%s', fields{i})) = Res_close(1).(fields{i});
    end
    
    data_far =  cell2mat({Res_far.(fields{i})});
    if size(data_close,2) == nb_resampling
        Res_far_summary.(sprintf('%s_mean', fields{i}))   = mean(data_far,2);
        Res_far_summary.(sprintf('%s_std', fields{i}))    = std(data_far,[],2);
    else
        Res_far_summary.(sprintf('%s', fields{i})) = Res_close(1).(fields{i});
    end
    if size(data_close,2) == nb_resampling
        Res_summary.(sprintf('%s_mean', fields{i}))   = (mean(data_close,2) + mean(data_far,2)) / 2;
        Res_summary.(sprintf('%s_std', fields{i}))    = ( std(data_close,[],2) + std(data_far,[],2)) / 2;

    else
        Res_summary.(sprintf('%s', fields{i})) = Res_close(1).(fields{i});
    end
end
end

function [Res_close,Res_far] = Compute_AUC_global(Jtheo, Jmeas,  VoisinsOA, clusters, ordreVois, thresholds, area)
% ***********************************************************************
% [ROC_Result, AUC_Result] = Compute_AUC_global(Jtheo, Jmeas,  VOISINS, ListMethods, clusters, ordreVois)
% Draw a ROC curve given a a Patch structure used ofr simulation ans a SourceLoc structure for 
% ***********************************************************************
% Inputs : 
%   - Jtheo: Ground truth. Binnary map
%   - Jmeas: Measured map. Scaled between 0 and 1.
%   - VoisinsOA
%   - clusters
%   - ordreVois
%   - thresholds. Threshold used for the computqtion of the ROC curve. By
%   default: linspace(0, 1, 100)
%   - area: 
% ***********************************************************************
% Outputs : 
%   - Res_close, Res_far. structure containing the results for the close
%   and far sources. 
% ***********************************************************************
% C. Grova - Montreal Neurological Institute - 21 03 2005 
% ***********************************************************************

if nargin < 8 || isempty(thresholds)
    thresholds = linspace(0, 1, 100);  
end
if nargin < 9 || isempty(area)
    area = [];  
end

nbvertices = length(Jtheo);

Itheo       = find(Jtheo ~=0);
Itheo       = unique(Itheo);

% Target number for the number of vertex in JFictif
nTarget = length(Itheo);

% --------------------------------------------------------------------
% Computation of AUC_close
% --------------------------------------------------------------------

% definition du Support proche a l aide d un voinnage a l ordre 10 du patch
SupportClose = [];
for i = 1:length(Itheo)
    SupportClose = unique([ SupportClose VoisinsOA{ordreVois,Itheo(i)}]);
end

% Tirage des sources fictives
Isources_fictives   = intersect(Itheo, SupportClose);
Isources_fictives   = setxor(SupportClose, Isources_fictives);

if 0 && ~isempty(clusters)
    Jfictif = AUC_create_JFictif_cluster(Jmeas, clusters, Isources_fictives, nbvertices, nTarget);
else
    Jfictif = AUC_create_JFictif_random(Isources_fictives, nbvertices, nTarget);
end

Res_close       = Compute_RocParam(Jmeas, Jtheo, Jfictif,  thresholds, area);
Res_close.AUC   = ComputeAUC(1-Res_close.specificity,Res_close.sensitivity);
if ~isempty(area)
    Res_close.AUC_area = ComputeAUC(1-Res_close.specificity_area,Res_close.sensitivity_area);
end

% --------------------------------------------------------------------
% Computation of AUC_far
% --------------------------------------------------------------------

SupportFar = setdiff(1:length(Jmeas), SupportClose);

Isources_fictives = intersect(Itheo, SupportFar);
Isources_fictives = setxor(SupportFar, Isources_fictives);

% Analyse par courbe ROC : 
if ~isempty(clusters)
    Jfictif = AUC_create_JFictif_cluster(Jmeas, clusters, Isources_fictives, nbvertices, nTarget);
else
    Jfictif = AUC_create_JFictif_random(Isources_fictives, nbvertices, nTarget);
end

Res_far     = Compute_RocParam(Jmeas, Jtheo, Jfictif,  thresholds, area);
Res_far.AUC = ComputeAUC(1-Res_far.specificity,Res_far.sensitivity);
if ~isempty(area)
    Res_far.AUC_area = ComputeAUC(1-Res_far.specificity_area,Res_far.sensitivity_area);
end
end

function [Roc_struct] = Compute_RocParam(J, Jtheo, Jfictif, thresholds, area)

% ***********************************************************************
% [Roc_struct] = Compute_RocParam(J, Jth, thresholds)
% Compute ROC parameters form the computed distribution J and the theoretical distribution Jth
% ***********************************************************************
% Inputs : 
%    J : vector (nbvertices) containing the probability of each source to be active (to be thresholded 
%    between 0 and 1
%    Jtheo : vector (nbvertices) containing 1 if the dipole is active 0 otherwise
%    thresholds : list of thresholds to apply to J
%
% ***********************************************************************
% Outputs :
%   the sutructure Roc_struct containing following fields
%    .specificity : tn / (tn + fp) nb de ce qui n a pas ete detecte, sur tout ce qui ne devrait pas etre detecte
%    .sensitivity : tp / (tp + fn) nb de ce qui a ete detecte, sur tout ce qui devrait etre detecte
%    .ppv  : positive predictive value : tp  / (tp + fp)
%    .npv  : negative predictive value : tn / (fn + tn)
%    .dice : Dice: 2*TP/(2*TP + FP + FN)
%    .tp : true positive rate
%    .fp : false positive rate
%    .tn : true negative rate
%    .fn : false negative rate
%
% ***********************************************************************
% C. Grova - Montreal Neurological Institute - 05 01 2004 
% ***********************************************************************

if nargin < 5
    area = [];
end

nbtest = length(thresholds);


Roc_struct.tp = zeros(nbtest,1);
Roc_struct.tn = zeros(nbtest,1);
Roc_struct.fp = zeros(nbtest,1);
Roc_struct.fn = zeros(nbtest,1);
Roc_struct.specificity = zeros(nbtest,1);
Roc_struct.sensitivity = zeros(nbtest,1);
Roc_struct.ppv = zeros(nbtest,1);
Roc_struct.npv = zeros(nbtest,1);
Roc_struct.dice = zeros(nbtest,1);

if ~isempty(area)
    Roc_struct.tp_area = zeros(nbtest,1);
    Roc_struct.tn_area = zeros(nbtest,1);
    Roc_struct.fp_area = zeros(nbtest,1);
    Roc_struct.fn_area = zeros(nbtest,1);
    Roc_struct.specificity_area = zeros(nbtest,1);
    Roc_struct.sensitivity_area = zeros(nbtest,1);
    Roc_struct.ppv_area = zeros(nbtest,1);
    Roc_struct.npv_area = zeros(nbtest,1);
    Roc_struct.dice_area = zeros(nbtest,1);
end

for i = 1:nbtest
    Jthresh = (J >= thresholds(i));

    Roc_struct.tp(i) = sum(Jthresh .* Jtheo);  
    Roc_struct.tn(i) = sum(~(Jthresh) .* Jfictif);  
    Roc_struct.fp(i) = sum(Jthresh .* Jfictif);  
    Roc_struct.fn(i) = sum(~(Jthresh) .* (Jtheo)); 
    
    
    Roc_struct.specificity(i) =   Roc_struct.tn(i) /  (Roc_struct.tn(i) + Roc_struct.fp(i));
    Roc_struct.sensitivity(i) =   Roc_struct.tp(i) /  (Roc_struct.tp(i) + Roc_struct.fn(i));
    Roc_struct.ppv(i)         =   Roc_struct.tp(i) /  (Roc_struct.tp(i) + Roc_struct.fp(i));
    Roc_struct.npv(i)         =   Roc_struct.tn(i) /  (Roc_struct.tn(i) + Roc_struct.fn(i));
    Roc_struct.dice(i)        = 2*Roc_struct.tp(i) /  (2*Roc_struct.tp(i) + Roc_struct.fp(i) + Roc_struct.fn(i));
    
    
    if ~isempty(area)
        Roc_struct.tp_area(i) = sum(area (Jthresh .* Jtheo)); 
        Roc_struct.tn_area(i) = sum(area (~(Jthresh) .* Jfictif));
        Roc_struct.fp_area(i) = sum(area (Jthresh .* Jfictif)); 
        Roc_struct.fn_area(i) = sum(area (~(Jthresh) .* (Jtheo))); 
    
    
        Roc_struct.specificity_area(i) =   Roc_struct.tn_area(i) /  (Roc_struct.tn_area(i) + Roc_struct.fp_area(i));
        Roc_struct.sensitivity_area(i) =   Roc_struct.tp_area(i) /  (Roc_struct.tp_area(i) + Roc_struct.fn_area(i));
        Roc_struct.ppv_area(i)         =   Roc_struct.tp_area(i) /  (Roc_struct.tp_area(i) + Roc_struct.fp_area(i));
        Roc_struct.npv_area(i)         =   Roc_struct.tn_area(i) /  (Roc_struct.tn_area(i) + Roc_struct.fn_area(i));
        Roc_struct.dice_area(i)        = 2*Roc_struct.tp_area(i) /  (2*Roc_struct.tp_area(i) + Roc_struct.fp_area(i) + Roc_struct.fn_area(i));
    end

end
Roc_struct.thresholds = thresholds;
end

%% ===== HELPER FUNCTIONS =====
function options = get_default_options(Time, type)
    options = struct(); 
    options.type = type;
    switch(type)
        case 'oscilation'
            options.freq        = 0.1;
            options.peak_time   = Time(round(length(Time)/2));
            options.duration    = 40;
        case 'task' 
            hrf_types   = process_nst_glm_fit('get_hrf_types');
            options.hrf = hrf_types.GAMMA;
            options.task_duration = 10;
            options.rest_duration = [30, 40];
    end
end

function [y, event] = simulate_oscilations(Time, options)
    event = db_template('Event');
    event.label = 'oscilations';
    event.times = options.peak_time;
    event.color = [ .4    .4    1];
    event.epochs = 1;
    sigma   = options.duration / 2.354;
    Tc      = (Time - options.peak_time);
    y = cos(2 * pi * options.freq * Tc) .* exp(-Tc.^2 ./ (2*sigma^2));
end