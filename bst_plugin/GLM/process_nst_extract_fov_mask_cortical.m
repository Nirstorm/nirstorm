function varargout = process_nst_extract_fov_mask_cortical( varargin ) %#ok<STOUT,>

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
% Authors: Thomas Vincent (2019)

eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Extract cortical FOV mask';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS', 'Work in progress'};
    sProcess.Index       = 1815;
    sProcess.isSeparator = 1;
    sProcess.Description = '';
    % todo add a new tutorial
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'results'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    sProcess.options.save_surf_mask.Comment = 'Save surface mask';
    sProcess.options.save_surf_mask.Type    = 'checkbox';
    sProcess.options.save_surf_mask.Value   =  0;
 
    sProcess.options.save_atlas.Comment = 'Save to scout atlas (disabled if empty): ';
    sProcess.options.save_atlas.Type    = 'text';
    sProcess.options.save_atlas.Value   =  '';
    
    sProcess.options.do_atlas_inter.Comment = 'Use atlas';
    sProcess.options.do_atlas_inter.Type    = 'checkbox';
    sProcess.options.do_atlas_inter.Value   =  0;

    sProcess.options.atlas.Comment  = 'Atlas';
    sProcess.options.atlas.Type     = 'atlas';
    sProcess.options.atlas.Value    = '';
    
    sProcess.options.keep_full_atlas_rois.Comment = 'Keep full atlas regions';
    sProcess.options.keep_full_atlas_rois.Type    = 'checkbox';
    sProcess.options.keep_full_atlas_rois.Value   =  0;    
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

OutputFiles = {};

ResultsMat = in_bst_results(sInputs(1).FileName);
if ~isempty(ResultsMat.HeadModelFile) && strcmp(ResultsMat.HeadModelType, 'surface')
    sHeadModel = in_bst_headmodel(ResultsMat.HeadModelFile);
    assert(strcmp(sHeadModel.SurfaceFile, ResultsMat.SurfaceFile));
    fov_mask = any(squeeze(any(sHeadModel.Gain > 0)))';
else
    fov_mask = any(full(ResultsMat.ImageGridAmp)>0, 2);
end

if sProcess.options.do_atlas_inter.Value
     % from process_source_atlas.m
    SurfaceMat = in_tess_bst(ResultsMat.SurfaceFile);
    iAtlas = find(strcmpi({SurfaceMat.Atlas.Name}, sProcess.options.atlas.Value));
    if isempty(iAtlas)
        error(['Atlas not found: "' sProcess.options.atlas.Value '"']);
    end
    atlas = SurfaceMat.Atlas(iAtlas);
    atlas_mask = zeros(size(fov_mask));
    for iscout=1:length(atlas.Scouts)
        atlas_mask(atlas.Scouts(iscout).Vertices) = iscout;
    end
    fov_mask = fov_mask .* atlas_mask;
else
    atlas = [];
end

if sProcess.options.save_surf_mask.Value
    sStudy = bst_get('Study', sInputs(1).iStudy);
   [sStudy, ResultFile] = nst_bst_add_surf_data(fov_mask, 1, sHeadModel, 'mask', ...
                                                [ResultsMat.Comment ' | FOV mask'], ...
                                                sInputs, sStudy, 'FOV mask extraction');
    OutputFiles{end+1} = ResultFile; 
end

if ~isempty(sProcess.options.save_atlas.Value)
    sCortex = in_tess_bst(ResultsMat.SurfaceFile);
 
    out_atlas_name = sProcess.options.save_atlas.Value;
    atlas_idx = find(strcmp({sCortex.Atlas.Name}, out_atlas_name), 1);
    if isempty(atlas_idx) % new one
        atlas_idx = size(sCortex.Atlas, 2) + 1;
    end
    sCortex.Atlas(atlas_idx).Name = out_atlas_name;
    
    fov_roi_indexes = sort(unique(fov_mask))';
    if fov_roi_indexes(1) == 0
        fov_roi_indexes = fov_roi_indexes(2:end);
    end
    
    sCortex.Atlas(atlas_idx).Scouts = repmat(db_template('Scout'), 0);
    for iroi=1:length(fov_roi_indexes)
        roi_index = fov_roi_indexes(iroi);
        roi_scout = db_template('Scout');
        if sProcess.options.keep_full_atlas_rois.Value && ~isempty(atlas)
            roi_scout =  atlas.Scouts(roi_index);
        else
            roi_mask = fov_mask==roi_index;
            roi_scout.Vertices = find(roi_mask);
            roi_scout.Seed = roi_scout.Vertices(1);
            if ~isempty(atlas)
                roi_scout.Color = atlas.Scouts(roi_index).Color;
                roi_scout.Label = atlas.Scouts(roi_index).Label;
            else
                roi_scout.Color = rand(1, 3);
                roi_scout.Label = sptrinf('roi_%d');
            end
        end
        
        sCortex.Atlas(atlas_idx).Scouts(iroi) = roi_scout;
    end
    
    bst_save(file_fullpath(ResultsMat.SurfaceFile), sCortex, 'v7');
    db_save();
    
end

end
