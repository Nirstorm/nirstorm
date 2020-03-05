function varargout = process_nst_sub_headmodel( varargin )

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
% Authors: Thomas Vincent (2018)

%TODO: put this in extra installation scenario

eval(macro_method);
end


function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Head model from precomputed one';
sProcess.FileTag     = '';
sProcess.Category    = 'File2';
sProcess.SubGroup    = 'NIRS';
sProcess.Index       = 1409;
sProcess.Description = '';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data', 'raw'};
% Definition of the outputs of this process
sProcess.OutputTypes = {'data', 'raw'};
sProcess.nInputs     = 2;
sProcess.nMinFiles   = 1;
sProcess.isPaired    = 0;

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputA, sInputB) %#ok<DEFNU>
OutputFiles = {};

% sInputA: for which to compute sub head models based on actual montage
% sInputB: associated with larger head model where to take sub head models from

% Load head model from sInputB
sStudyB = bst_get('Study', sInputB.iStudy);
if isempty(sStudyB.iHeadModel)
    bst_error('No head model found in File B. Consider process "Compute head model from fluence"');
    return;
end
parent_head_model_fn = sStudyB.HeadModel(sStudyB.iHeadModel).FileName;
parent_head_model = in_bst_headmodel(parent_head_model_fn);

% Get pair names from sInputA
ChannelMatA = in_bst_channel(sInputA(1).ChannelFile);

montage_info_A = nst_montage_info_from_bst_channels(ChannelMatA.Channel);
pair_names = montage_info_A.pair_names;

sub_sensitivity_mat = process_nst_import_head_model('get_sensitivity_from_chans', ...
                                                    parent_head_model, pair_names);

% Save sub head model
sStudyA = bst_get('Study', sInputA.iStudy);
[sSubjectA, iSubjectA] = bst_get('Subject', sInputA.SubjectName);

% Create structure
HeadModelMat = db_template('headmodelmat');
HeadModelMat.Gain           = sub_sensitivity_mat;
HeadModelMat.HeadModelType  = 'surface';
HeadModelMat.SurfaceFile    = parent_head_model.SurfaceFile;
HeadModelMat.Comment       = 'NIRS head model';

HeadModelMat.pair_names = pair_names;
HeadModelMat = bst_history('add', HeadModelMat, 'compute', ...
                          ['Sub head model from ' file_short(parent_head_model_fn)]);
% Output file name
from_fluence = ~isempty(strfind(parent_head_model_fn,'mcx_fluence'));
if from_fluence
    suffix = '_mcx_fluence';
else
    suffix = '';
end
HeadModelFile = bst_fullfile(bst_fileparts(file_fullpath(sStudyA.FileName)), ...
                                           sprintf('headmodel_nirs%s.mat', suffix));
HeadModelFile = file_unique(HeadModelFile);
% Save file
bst_save(HeadModelFile, HeadModelMat, 'v7');

newHeadModel = db_template('HeadModel');
newHeadModel.FileName = file_short(HeadModelFile);
if from_fluence
    newHeadModel.Comment = 'NIRS head model from fluence';
else
    newHeadModel.Comment = 'NIRS head model';
end

newHeadModel.HeadModelType  = 'surface';    
% Update Study structure
iHeadModel = length(sStudyA.HeadModel) + 1;
if ~isempty(sStudyA.HeadModel)
    sStudyA.HeadModel(end+1) = newHeadModel;
else
    sStudyA.HeadModel = newHeadModel;
end
sStudyA.iHeadModel = iHeadModel;
% Update DataBase
bst_set('Study', sInputA.iStudy, sStudyA);
panel_protocols('UpdateNode', 'Study', sInputA.iStudy);

% Save database
db_save();
end
