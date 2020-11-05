function varargout = process_nst_glm_contrast_mask( varargin )
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
    sProcess.Comment     = 'GLM - group t-contrast mask';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = {'NIRS', 'GLM'};
    sProcess.Index       = 1603;
    sProcess.isSeparator = 0;
    sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/%5BWIP%5D-GLM-implementation';
    % todo add a new tutorial
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'presults'};
    sProcess.OutputTypes = {'results'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    sProcess = process_extract_pthresh('DefineOptions', sProcess);
    
    sProcess.options.do_atlas_inter.Comment = 'Produce ROIs from atlas';
    sProcess.options.do_atlas_inter.Type    = 'checkbox';
    sProcess.options.do_atlas_inter.Value   =  0;

    sProcess.options.atlas.Comment  = 'Atlas';
    sProcess.options.atlas.Type     = 'atlas';
    sProcess.options.atlas.Value    = '';
    
    sProcess.options.min_atlas_roi_size.Comment = 'Min size of atlas ROI:';
    sProcess.options.min_atlas_roi_size.Type    = 'value';
    sProcess.options.min_atlas_roi_size.Value   = {2, '', 0};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)  %#ok<DEFNU>
    Comment = sProcess.Comment; 
end

function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>
    OutputFiles={};
    
    [StatThreshOptions, strCorrect] = process_extract_pthresh('GetOptions', sProcess);
    ResultsMat = in_bst_results(sInput.FileName, 0);
  
    if sProcess.options.do_atlas_inter.Value
        % from process_source_atlas.m
        SurfaceMat = in_tess_bst(ResultsMat.SurfaceFile);
        iAtlas = find(strcmpi({SurfaceMat.Atlas.Name}, sProcess.options.atlas.Value));
        if isempty(iAtlas)
            error(['Atlas not found: "' sProcess.options.atlas.Value '"']);
        end
        atlas = SurfaceMat.Atlas(iAtlas);
        atlas_mask = zeros(size(ResultsMat.tmap));
        for iscout=1:length(atlas.Scouts)
            atlas_mask(atlas.Scouts(iscout).Vertices) = iscout;
        end
        %TODO: precompute atlas label surface map. Maybe take a look at
        %ResultsMat.GridAtlas?
    else
        atlas = [];
    end
     
    fov_mask = ResultsMat.df > 0;
    
    func_mask = process_extract_pthresh('Compute', ResultsMat, StatThreshOptions) ~= 0;
    
    if ~isempty(atlas)
        comment = [ResultsMat.Comment ' | atlas-based activity mask'];
        func_mask = func_mask .* atlas_mask;
        roi_indexes = sort(unique(func_mask))';
        if roi_indexes(1) == 0
            roi_indexes = roi_indexes(2:end);
        end
        for iroi=roi_indexes
            mm = func_mask==iroi;
            if sum(mm) < sProcess.options.min_atlas_roi_size.Value{1}
                func_mask(mm) = 0;
            end
        end
        
        fov_roi_mask = fov_mask .* atlas_mask;
        fov_roi_indexes = sort(unique(fov_roi_mask))';
        if fov_roi_indexes(1) == 0
            fov_roi_indexes = fov_roi_indexes(2:end);
        end
    else
        fov_roi_indexes = [];
        comment = [ResultsMat.Comment ' | activity mask'];
    end
    
    DataMat = db_template('resultsmat');

    DataMat.from_atlas = atlas;
    DataMat.fov_mask = fov_mask;
    DataMat.fov_roi_indexes = fov_roi_indexes;
    
    DataMat.contrast_name = ResultsMat.contrast_name;
    
    DataMat.ImageGridAmp  = func_mask;
    DataMat.ImagingKernel = [];
    DataMat.Comment       = comment;
    DataMat.Function      = 'pthresh';
    DataMat.Time          = ResultsMat.Time;
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
    
    % Add history entry
    DataMat = bst_history('add', DataMat, 'pthresh', ['Setting the stat threshold: ' strCorrect]);
    DataMat = bst_history('add', DataMat, 'pthresh', ['Original file: ' sInput.FileName]);
    % Output file tag
    fileTag = bst_process('GetFileTag', sInput.FileName);
    fileTag = [fileTag(2:end) '_pthresh_mask'];
    % Output filename
    DataFile = bst_process('GetNewFilename', bst_fileparts(sInput.FileName), fileTag);
    % Save on disk
    bst_save(DataFile, DataMat, 'v6');
    % Register in database
    db_add_data(sInput.iStudy, DataFile, DataMat);
    % Return data file
    OutputFiles{1} = DataFile;
end
