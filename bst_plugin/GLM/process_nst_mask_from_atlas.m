function varargout = process_nst_mask_from_atlas( varargin )

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

function sProcess = GetDescription() 
% Description the process
sProcess.Comment     = 'Cortical mask from atlas scouts';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'NIRS';
sProcess.Index       = 1901;
sProcess.Description = '';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data', 'raw', 'results'};
% Definition of the outputs of this process
sProcess.OutputTypes = {'results', 'results', 'results'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;


sProcess.options.surface_name.Comment  = 'Surface: ';
sProcess.options.surface_name.Type     = 'text';
sProcess.options.surface_name.Value    = '';

sProcess.options.atlas.Comment  = 'Atlas';
sProcess.options.atlas.Type     = 'atlas';
sProcess.options.atlas.Value    = '';

sProcess.options.scout_names.Comment  = 'Regions (comma-separated list of scout labels): ';
sProcess.options.scout_names.Type     = 'text';
sProcess.options.scout_names.Value    = '';

sProcess.options.output_comment.Comment  = 'Output label: ';
sProcess.options.output_comment.Type     = 'text';
sProcess.options.output_comment.Value    = 'Atlas-based mask';

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) 

% Get scout vertices & load head mesh
atlas_name = sProcess.options.atlas.Value;
surface_name = sProcess.options.surface_name.Value;
scout_names = cellfun(@(s) strtrim(s), strsplit(sProcess.options.scout_names.Value, ','), ...
                      'UniformOutput', false);
                  
[sSubject, iSubject] = bst_get('Subject', sInputs(1).SubjectName);
if ~isempty(sSubject.Surface)
    i_surface = find(strcmp({sSubject.Surface.Comment}, surface_name));
    if ~isempty(i_surface)
        if length(i_surface) > 1
            warning('Multiple surfaces with name %s, taking 1st one.', surface_name);
        end
        surface_file = sSubject.Surface(i_surface(1)).FileName;
        sSurf = in_tess_bst(surface_file);
    else
        error('Surface %s not found.', surface_name);
    end
else
    error('No surface available');
end

[mask, atlas] = Compute(sSubject, sSurf, atlas_name, scout_names);

DataMat.ImageGridAmp  = mask;
DataMat = db_template('resultsmat');

DataMat.from_atlas = atlas;

% DataMat.fov_mask = fov_mask;
% DataMat.fov_roi_indexes = fov_roi_indexes;

DataMat.ImageGridAmp  = mask;
DataMat.ImagingKernel = [];
DataMat.Comment       = sProcess.options.output_comment.Value;
DataMat.Function      = '';
DataMat.Time          = 1;
DataMat.DataFile      = [];
DataMat.HeadModelFile = [];
DataMat.SurfaceFile   = surface_file;

% Output filename
DataFile = bst_process('GetNewFilename', bst_fileparts(sInputs(1).FileName), 'results_mask_');
% Save on disk
bst_save(DataFile, DataMat, 'v6');
% Register in database
db_add_data(sInputs(1).iStudy, DataFile, DataMat);
% Return data file
OutputFiles{1} = DataFile;
    
end

function [mask, atlas] = Compute(sSubject, sSurf, atlas_name, scout_names)

iatlas_selected = find(strcmp({sSurf.Atlas.Name}, atlas_name));

if length(iatlas_selected) > 1
    warning('Multiple atlas with name %s, taking 1st one.', atlas_name);
end

if isempty(iatlas_selected)
    error('Atlas %s not found', atlas_name);
else
    atlas = sSurf.Atlas(iatlas_selected(1));
    scouts = atlas.Scouts;
    
    scouts_found = ismember(scout_names, {scouts.Label});
    if ~all(scouts_found)
        error('Scouts not found: %s', strjoin(scout_names(~scouts_found), ', '));
    end
    i_selected_scouts = find(ismember({scouts.Label}, scout_names));
    mask = zeros(size(sSurf.Vertices,1), 1);
    for ii=1:length(i_selected_scouts)
        i_selected_scout = i_selected_scouts(ii);  
        mask(scouts(i_selected_scout).Vertices) = i_selected_scout;
    end
end
end