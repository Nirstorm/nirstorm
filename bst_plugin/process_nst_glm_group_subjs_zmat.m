function varargout = process_nst_glm_group_subjs_zmat( varargin )
% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2017 University of Southern California & McGill University
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
% Authors: Thomas Vincent, 2019
%  
eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'GLM - group roi-based subjects z-matrix';
    sProcess.Category    = 'Custom2';
    sProcess.SubGroup    = 'NIRS - wip';
    sProcess.Index       = 1403;
    sProcess.isSeparator = 0;
    sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/%5BWIP%5D-GLM-implementation';
    % todo add a new tutorial
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'matrix'};
    
    sProcess.nInputs     = 2;
    sProcess.nMinFiles   = 1;
    sProcess.isPaired    = 0;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)  %#ok<DEFNU>
    Comment = sProcess.Comment; 
end

function [OutputFiles1, OutputFiles2]  = Run(sProcess, sInputs1, sInput2) %#ok<DEFNU>
    OutputFiles1 = {};
    OutputFiles2 = {};
    
    % Retrieve all subject con maps from sInputs1
    con_data = in_bst_results(sInputs1(1).FileName);
    all_con = zeros(length(sInputs1), size(con_data.ImageGridAmp, 1));
    all_con_std = zeros(size(all_con));
    for isubj=1:length(sInputs1)
        con_data = in_bst_results(sInputs1(isubj).FileName);
        all_con(isubj, :) = con_data.ImageGridAmp;
        all_con_std(isubj, :) = con_data.contrast_std;
    end
    
    if length(sInput2) > 1
        error('Expected only 1 item for 2nd input (roi mask)');
    end
    
    % Retrieve roi mask from sInput2
    mask_data = in_bst_results(sInput2.FileName);
    mask = mask_data.ImageGridAmp;
  
    % Define set of ROIs
    if isfield(mask_data, 'fov_roi_indexes') && ~isempty(mask_data.fov_roi_indexes) && ...
            ~isempty(mask_data.from_atlas)
        roi_indexes = mask_data.fov_roi_indexes;
    else
        roi_indexes = sort(unique(mask));
        if roi_indexes(1) == 0
            roi_indexes = roi_indexes(2:end);
        end
    end
    if  ~isempty(mask_data.from_atlas)
        roi_names = {mask_data.from_atlas.Scouts(roi_indexes).Label};
    else
        roi_names = arrayfun(@(n) num2str(n), roi_indexes, 'UniformOutput', 0);
    end
    
    % Loop over full list of ROIs within FOV and fill matrix with mean of 
    % z-values (beta / beta_std) within each region
    all_roi_z = zeros(length(sInputs1), length(roi_names));
    for iroi=1:length(roi_indexes)
        roi_index = roi_indexes(iroi);
        roi_mask = mask==roi_index;
        if sum(roi_mask) > 0
            mean_std = mean(all_con_std(:, roi_mask), 2);
            if mean_std > 0
                all_roi_z(:, iroi) = mean(all_con(:, roi_mask), 2) ./ mean_std;
            end
        end
    end

    
    
    % Save
    DataMat = db_template('matrixmat');
    DataMat.Value = all_roi_z;
    DataMat.Comment = [sInputs1(1).Comment ' | subjects roi z-values'];
    DataMat.DisplayUnits = 'z-score';
    DataMat.Time = 1:length(sInputs1);
    DataMat.RowNames = {sInputs1.SubjectName};
    DataMat.ColNames = roi_names;
    DataMat.Description = roi_names;
    
    OutputFiles1{1} = bst_process('GetNewFilename', fileparts(sInput2.FileName), 'matrix');
    % Save file
    bst_save(OutputFiles1{1}, DataMat, 'v6');
    % Add file to database structure
    db_add_data(sInput2.iStudy, OutputFiles1{1}, DataMat);
end
