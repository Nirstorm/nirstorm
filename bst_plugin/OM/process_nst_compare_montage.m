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

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end

function OutputFile = Run(sProcess, sInput)
    OutputFile = '';

    % Load subject
    sSubject = bst_get('Subject', sInput.SubjectName);

    % Load ROI
    sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);
    sHead   = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);

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
    metrics = compute_metrics_montage(G, coverage, sGroundTruth);
    
    metric_table = struct2table(metrics);
    metric_table.Properties.RowNames{1} = sInput.Condition;

    OutputFile = nst_save_table_in_bst(metric_table, sInput.SubjectName, sInput.Condition, 'Metrics results');

end



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