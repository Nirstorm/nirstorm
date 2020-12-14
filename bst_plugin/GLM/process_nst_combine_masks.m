function varargout = process_nst_combine_masks( varargin )
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
    sProcess.Comment     = 'Combine masks';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS', 'Work in progress'};
    sProcess.Index       = 1813;
    sProcess.isSeparator = 0;
    sProcess.Description = '';
    % todo add a new tutorial
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'results'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    sProcess.options.operation.Comment = 'Operation';
    sProcess.options.operation.Type    = 'combobox';
    operations = get_mask_combinations();
    sProcess.options.operation.Value   = {operations.intersection,...
                                          fieldnames(operations)};

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)  %#ok<DEFNU>
    Comment = sProcess.Comment; 
end

function operations = get_mask_combinations()
   operations.intersection = 1;
   operations.union = 2;
end

function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles={};
    
    operations = get_mask_combinations();
    operation_idx = sProcess.options.operation.Value{1};
    operation_str = sProcess.options.operation.Value{2}{sProcess.options.operation.Value{1}};

    for iInput=1:length(sInputs)
        ResultsMat = in_bst_results(sInputs(iInput).FileName, 0);
        % TODO: make some checks on input mask (binary values)

        if iInput==1
            combined_mask = ResultsMat.ImageGridAmp;
            roi_ids = unique(combined_mask);
        else
            roi_ids = union(roi_ids, ResultsMat.ImageGridAmp);
            for iroi=1:length(roi_ids)
                roi_id = roi_ids(iroi);
                if roi_id ~=0 % ignore background
                    roi_mask = combined_mask == roi_id;
                    switch operation_idx
                        case operations.intersection
                            cm = roi_mask & ResultsMat.ImageGridAmp == roi_id;
                        case operations.union
                            cm = roi_mask | ResultsMat.ImageGridAmp == roi_id;
                        otherwise
                            error('Invalid mask operation');
                    end
                    combined_mask(roi_mask) = 0;
                    combined_mask(cm) = roi_id;
                end
            end
        end
    end
       
    DataMat = db_template('resultsmat');
    
    
    % TODO check consistency along inputs
    if isfield(ResultsMat, 'from_atlas')
        DataMat.from_atlas = ResultsMat.from_atlas;
        DataMat.fov_mask = ResultsMat.fov_mask; %TODO
        DataMat.fov_roi_indexes = ResultsMat.fov_roi_indexes;
    end
    if isfield(ResultsMat, 'contrast_name')
        DataMat.contrast_name = ResultsMat.contrast_name;
    end
    
    DataMat.ImageGridAmp  = combined_mask;
    DataMat.ImagingKernel = [];
    
    [varying, common_pref, common_suf] = str_remove_common({sInputs.Comment});
    varying(cellfun(@(c) isempty(c), varying)) = [];
    [common_pref, common_suf] = remove_mask_str(common_pref, common_suf);
    
    DataMat.Comment       = [common_pref strjoin(varying, '_') common_suf ' | ' operation_str ' mask'];
    DataMat.Function      = operation_str;
    DataMat.Time          = ResultsMat.Time; %TODO insure to properly handle time
    DataMat.DataFile      = [];
    DataMat.HeadModelFile = [];
    DataMat.HeadModelType = ResultsMat.HeadModelType;
    DataMat.nComponents   = ResultsMat.nComponents;
    DataMat.GridLoc       = ResultsMat.GridLoc;
    DataMat.SurfaceFile   = ResultsMat.SurfaceFile;
    DataMat.GoodChannel   = ResultsMat.GoodChannel;
    DataMat.ChannelFlag   = ResultsMat.ChannelFlag;
    DataMat.History       = ResultsMat.History;
    DataMat.ColormapType  = ResultsMat.ColormapType;
   
    % Output file tag
    % fileTag = bst_process('GetFileTag', sInputs(1).FileName);
    fileTag = ['results_' operation_str '_mask'];
    % Output filename
    DataFile = bst_process('GetNewFilename', bst_fileparts(sInputs(1).FileName), fileTag);
    % Save on disk
    bst_save(DataFile, DataMat, 'v6');
    % Register in database
    db_add_data(sInputs(1).iStudy, DataFile, DataMat);
    % Return data file
    OutputFiles{1} = DataFile;
end


function [s1, s2] = remove_mask_str(s1,  s2)

seps = {' | ', ' - ', '_'};
for isep=1:length(seps)
    s1 = replace(s1, [seps{isep} 'mask'], '');
    s2 = replace(s2, [seps{isep} 'mask'], '');
end
s1 = replace(s1, 'mask', '');
s2 = replace(s2, 'mask', '');
end