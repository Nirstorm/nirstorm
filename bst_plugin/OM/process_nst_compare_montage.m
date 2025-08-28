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
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS', 'Sources'};
    sProcess.Index       = 1408;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes  = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    
    % === CLUSTERS
    sProcess.options.scouts.Comment = '';
    sProcess.options.scouts.Type    = 'scout';
    sProcess.options.scouts.Value   = {};


    % Definition of the options
    % === TARGET
    % File selection options
    SelectOptions = {...
        '', ...                            % Filename
        '', ...                            % FileFormat
        'save', ...                        % Dialog type: {open,save}
        'Save text file...', ... % Window title
        'ExportData', ...                  % LastUsedDir: {ImportData,ImportChannel,ImportAnat,ExportChannel,ExportData,ExportAnat,ExportProtocol,ExportImage,ExportScript}
        'single', ...                      % Selection mode: {single,multiple}
        'files', ...                        % Selection mode: {files,dirs,files_and_dirs}
        {{'.mat'}, 'text file', 'mat'}, ... % Available file formats
        []};                          % DefaultFormats: {ChannelIn,DataIn,DipolesIn,EventsIn,MriIn,NoiseCovIn,ResultsIn,SspIn,SurfaceIn,TimefreqIn}
    % Option: MRI file
    sProcess.options.textFile.Comment = 'Output folder:';
    sProcess.options.textFile.Type    = 'filename';
    sProcess.options.textFile.Value   = SelectOptions;
       
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end

function OutputFile = Run(sProcess, sInput)
    OutputFile = '';

    % Load subject
    sSubject = bst_get('Subject', sInput.SubjectName);

    % identify inputs
    iMaps       = find(contains(lower({sInput.Comment}),{'sensitivities'}) & ...
                      ~contains(lower({sInput.Comment}),{'summed'}));

    % Compute distance from head to cortex
    sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);
    sHead   = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);

    distance = process_nst_cpt_cortex_to_head_distance('Compute', sCortex, sHead)*1000;
    threshold_distance = 20;
    idx_FOV = find(distance < threshold_distance);

    % Compute ground truth
    ROI         = sProcess.options.scouts.Value;
    iAtlas      = find(strcmp( {sCortex.Atlas.Name},ROI{1}));
    iRois       = cellfun(@(x)find(strcmp( {sCortex.Atlas(iAtlas).Scouts.Label},x)),   ROI{2});
    
    assert(length(iRois) == 1, 'Please select only one ROI');

    ROI_select = sCortex.Atlas(iAtlas).Scouts(iRois);
    sGroundTruth = zeros(size(sCortex.Vertices,1), 1);
    sGroundTruth(ROI_select.Vertices) = 1;

    % Prepare new surface - only FOV
    iRemoveVert = setdiff(1:size(sCortex.Vertices,1), idx_FOV);
    if ~isempty(iRemoveVert)
        [sCortex.Vertices, sCortex.Faces, sCortex.Atlas] = tess_remove_vert(sCortex.Vertices, sCortex.Faces, iRemoveVert, sCortex.Atlas);
        sCortex.VertConn = sCortex.VertConn(idx_FOV,idx_FOV);
    end

    % Compute area
    [~, VertArea] = tess_area(sCortex.Vertices, sCortex.Faces);

    % prepare new ground truth map based on the FOV
    Jtheo           = sGroundTruth;
    Jtheo           = Jtheo(idx_FOV);
    idx_truth       = find(Jtheo); 
    

    % Depth of the ROI
    sScalp = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);
    depth = 1000 * min(nst_pdist(sCortex.Vertices(idx_truth,:),sScalp.Vertices),[],2);

    % Excentricity of the ROI
    sMri        = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);

    NZ = sMri.SCS.NAS;
    OD = sMri.SCS.RPA;
    OG = sMri.SCS.LPA;

    Mesh.vertices = sCortex.Vertices;
    Mesh.faces = sCortex.Faces;
    mriCoord = cs_convert(sMri, 'scs', 'mri', Mesh.vertices)' * 1000;
    Mesh.vertices = mriCoord';
    Mesh.faces = sCortex.Faces;

    Mesh.vertices(:,1)=size(sMri.Cube,1)*sMri.Voxsize(1)-Mesh.vertices(:,1) ;%*1;%sMri.Voxsize(1);
    Mesh.vertices(:,2)=size(sMri.Cube,2)*sMri.Voxsize(2)-Mesh.vertices(:,2) ;%*1;%sMri.Voxsize(2);
    Mesh.vertices(:,3)=size(sMri.Cube,3)*sMri.Voxsize(3)-Mesh.vertices(:,3) ;%*1;%sMri.Voxsize(3);

    %eccentricity
    [~, all_eccentricity] = eccentricity(Mesh,NZ,OD,OG,eye(4,3), idx_truth, 1);
    
    % save information about the simulation     
    for iFile = 1:length(iMaps)

        results = struct();
        results.comment         = string( sInput(iMaps(iFile)).Condition);

        sData           = in_bst_results(sInput(iMaps(iFile)).FileName);


        % Metric on overlap
        overlap         = getCoverage(sSubject, sData.ImageGridAmp);
        sOverlap        = overlap(idx_FOV, :);

        % Information about the ground truth
        results.depth             = mean(depth);
        summary_func_eccentricity = {@mean, @median, @min, @max};
        for iFun = 1:length(summary_func_eccentricity)
            label = sprintf('eccentricity_%s', func2str(summary_func_eccentricity{iFun}));
            value = summary_func_eccentricity{iFun}(all_eccentricity(idx_truth));
            results.(label) = value;
        end

        results.NVertex = length(idx_truth);

        % Compute spatial metrics (at the time of the peak)
        % idx_truth  -- sMap_max
        results.coverage_vertex = sum(sOverlap(idx_truth) >= 1) / length(idx_truth);
        results.coverage_area = sum(VertArea(sOverlap(idx_truth) >= 1)) / sum(VertArea(idx_truth));

        results.median_overlap = median(sOverlap(idx_truth));
        results.max_overlap = max(sOverlap(idx_truth));

        threshopld = 1:(1+max(sOverlap(idx_truth)));
        cov_tmp = zeros(1, length(threshopld));
        for iThreshold = 1:length(threshopld)
            coverage_area = sum(VertArea(sOverlap(idx_truth) >= threshopld(iThreshold))) / sum(VertArea(idx_truth));
            cov_tmp(iThreshold) = coverage_area;
        end
    
        results.threshold = threshopld;
        results.cov_tmp = cov_tmp;

        % Metric on summed sensitivity
        sumSensitivity          = getSumSensitivity(sSubject, sData.ImageGridAmp);
        sSumSensitivity         = sumSensitivity(idx_FOV, :);

        results.sum_sensitivity = sum(sSumSensitivity(idx_truth));
        results.std_sensitivity = std(sSumSensitivity(idx_truth));
        results.median_sensitivity = median(sSumSensitivity(idx_truth));

        all_results(iFile,:)    = results;
    end


    figure;
    c = jet(length(all_results) + 10);
    
    for iFile = 1:length(all_results)
        subplot(4,3,iFile)
        plot(all_results(iFile).threshold, all_results(iFile).cov_tmp ,   'DisplayName',all_results(iFile).comment, 'Color',c(iFile, :) )
    
       
        area        = trapz(all_results(iFile).threshold, all_results(iFile).cov_tmp);
        max_val     = max(all_results(iFile).threshold);
        sumSensitivity = all_results(iFile).sum_sensitivity;
    
        comment = strrep(strrep(all_results(iFile).comment, '_', ' '), 'OM ROI ', '');
    
        ylim([0, 0.8])
        xlim([0,16])
        title(sprintf( '%s (M: %d,  A: %.2f, S: %.2f)', comment , max_val, area,sumSensitivity))
    end

    %3. Save the results 
    save(sProcess.options.textFile.Value{1}, 'all_results');
end


function overlap = getCoverage(sSubject, sensitivity)

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
    

    coverage_channel  = squeeze(sensitivity  > threshold);
    overlap = sum(coverage_channel,  2) ;

   
end

function sumSensitivity         = getSumSensitivity(sSubject, sensitivity)

    sumSensitivity = sum(sensitivity, 2);

end
